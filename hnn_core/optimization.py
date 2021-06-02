"""Parameter optimization functions."""

# Authors: Blake Caldwell <blake_caldwell@brown.edu>
#          Mainak Jas <mjas@mgh.harvard.edu>

from math import ceil, floor
from collections import OrderedDict

import numpy as np
import scipy.stats as stats
from scipy.optimize import fmin_cobyla

from .network import Network
from .dipole import rmse


def _get_range(val, multiplier):
    range_min = max(0, val - val * multiplier / 100.)
    range_max = val + val * multiplier / 100.
    ranges = {'initial': val, 'minval': range_min, 'maxval': range_max}
    return ranges


def split_by_evinput(params, sigma_range_multiplier, timing_range_multiplier,
                     synweight_range_multiplier):
    """ Sorts parameter ranges by evoked inputs into a dictionary

    Parameters
    ----------
    params: an instance of Params
        Full set of simulation parameters

    Returns
    -------
    sorted_evinput_params: Dict
        Dictionary with parameters grouped by evoked inputs.
        Keys for each evoked input are
        'mean', 'sigma', 'ranges', 'start',
        'end', 'cdf', 'weights', 'opt_start', 'opt_end'.
        sorted_evinput_params['evprox1']['ranges'] is a dict
        with keys as the parameters to be optimized and values
        indicating the ranges over which they should be
        optimized. E.g.,
        sorted_evinput['evprox1']['ranges'] =
        {
            'gbar_evprox_1_L2Pyr_ampa':
                {
                    'initial': 0.05125,
                    'minval': 0,
                    'maxval': 0.0915
                }
        }
        Elements are sorted by their start time.
    """

    evinput_params = {}
    for evinput_t in params['t_*']:
        id_str = evinput_t.lstrip('t_')
        timing_mean = float(params[evinput_t])

        if f'sigma_{evinput_t}' in params:
            timing_sigma = float(params['sigma_' + evinput_t])
        else:
            timing_sigma = 3.0
            print("Couldn't fing timing_sigma. Using default %f" %
                  timing_sigma)

        if timing_sigma == 0.0:
            # sigma of 0 will not produce a CDF
            timing_sigma = 0.01

        evinput_params[id_str] = {'mean': timing_mean, 'sigma': timing_sigma,
                                  'ranges': {}}

        evinput_params[id_str]['ranges']['sigma_' + evinput_t] = \
            _get_range(timing_sigma, sigma_range_multiplier)

        # calculate range for time
        timing_bound = timing_sigma * timing_range_multiplier
        range_min = max(0, timing_mean - timing_bound)
        range_max = min(float(params['tstop']), timing_mean + timing_bound)

        evinput_params[id_str]['start'] = range_min
        evinput_params[id_str]['end'] = range_max
        evinput_params[id_str]['ranges'][evinput_t] =  \
            {'initial': timing_mean, 'minval': range_min, 'maxval': range_max}

        # calculate ranges for syn. weights
        for label in params[f'gbar_{id_str}*']:
            value = params[label]
            ranges = _get_range(value, synweight_range_multiplier)    
            if value == 0.0:
                ranges['minval'] = value
                ranges['maxval'] = 1.0
            evinput_params[id_str]['ranges'][label] = ranges

    sorted_evinput_params = OrderedDict(sorted(evinput_params.items(),
                                               key=lambda x: x[1]['start']))
    return sorted_evinput_params


def _generate_weights(evinput_params, params, decay_multiplier):
    """Calculation of weight function for wRMSE calcuation

    Returns
    -------
    evinput_params : dict
        Adds the keys 'weights' and 'cdf' to evinput_params[input_name]
        and find evinput_params['opt_start'] and evinput_params['opt_end']
        based on weights.
    """
    num_step = ceil(params['tstop'] / params['dt']) + 1
    times = np.linspace(0, params['tstop'], num_step)

    for evinput_this in evinput_params.values():
        # calculate cdf using start time (minival of optimization range)
        evinput_this['cdf'] = stats.norm.cdf(
            times, evinput_this['start'], evinput_this['sigma'])

    for input_name, evinput_this in evinput_params.items():
        evinput_this['weights'] = evinput_this['cdf'].copy()

        for other_input, evinput_other in evinput_params.items():
            # check ordering to only use inputs after us
            # and don't subtract our own cdf(s)
            if (evinput_other['mean'] < evinput_this['mean'] or
                    input_name == other_input):
                continue

            decay_factor = decay_multiplier * \
                (evinput_other['mean'] - evinput_this['mean']) / params['tstop']
            evinput_this['weights'] -= evinput_other['cdf'] * decay_factor

        # weights should not drop below 0
        np.clip(evinput_this['weights'], a_min=0, a_max=None,
                out=evinput_this['weights'])

        # start and stop optimization where the weights are insignificant
        indices = np.where(evinput_this['weights'] > 0.01)
        evinput_this['opt_start'] = min(evinput_this['start'],
                                        times[indices][0])
        evinput_this['opt_end'] = max(evinput_this['end'],
                                      times[indices][-1])

        # convert to multiples of dt
        evinput_this['opt_start'] = floor(
            evinput_this['opt_start'] / params['dt']) * params['dt']
        evinput_params[input_name]['opt_end'] = ceil(
            evinput_this['opt_end'] / params['dt']) * params['dt']

    return evinput_params


def create_last_chunk(input_chunks):
    """ This creates a chunk that combines parameters for
    all chunks in input_chunks (final step)

    Parameters
    ----------
    input_chunks: List
        List of ordered chunks for optimization

    Returns
    -------
    chunk: dict
        Dictionary of with parameters for combined
        chunk (final step)
    """
    chunk = {'inputs': [], 'ranges': {}, 'opt_start': 0.0,
             'opt_end': 0.0}

    for evinput in input_chunks:
        chunk['inputs'].extend(evinput['inputs'])
        chunk['ranges'].update(evinput['ranges'])
        if evinput['opt_end'] > chunk['opt_end']:
            chunk['opt_end'] = evinput['opt_end']

    # wRMSE with weights of 1's is the same as regular RMSE.
    chunk['weights'] = np.ones(len(input_chunks[-1]['weights']))

    return chunk


def consolidate_chunks(inputs):
    """Consolidates inputs into optimization "chunks" defined by
    opt_start and opt_end

    Parameters
    ----------
    inputs: Dict
        Sorted dictionary of inputs with their parameters
        and weight functions

    Returns
    -------
    chunks: list
        Combine the evinput_params whenever the end is overlapping
        with the next.
    """
    chunks = list()
    for input_name in inputs:
        input_dict = inputs[input_name].copy()
        input_dict['inputs'] = [input_name]

        if (len(chunks) > 0 and
                input_dict['start'] <= chunks[-1]['end']):
            # update previous chunk
            chunks[-1]['inputs'].extend(input_dict['inputs'])
            chunks[-1]['end'] = input_dict['end']
            chunks[-1]['ranges'].update(input_dict['ranges'])
            chunks[-1]['opt_end'] = max(chunks[-1]['opt_end'],
                                        input_dict['opt_end'])
            # average the weights
            chunks[-1]['weights'] = (chunks[-1]['weights'] +
                                     input_dict['weights']) / 2
        else:
            # new chunk
            chunks.append(input_dict)

    # add one last chunk to the end
    if len(chunks) > 1:
        last_chunk = create_last_chunk(chunks)
        chunks.append(last_chunk)

    return chunks


def optrun(new_params, opt_params, params, opt_dpls):
    """ This is the function to run a simulation

    Parameters
    ----------
    new_params: Array
        List or numpy array with the parameters chosen by
        optimization engine. Order is consistent with
        opt_params['ranges'].
    opt_params : dict
        The optimization parameters
    params : dict
        The params dictionary.
    opt_dpls : dict
        Dictionary with keys exp_dpl and best for
        the experimental dipole and best dipole.

    Returns
    -------
    avg_rmse: float
        Weighted RMSE between data in dpl and exp_dpl
    """
    print("Optimization step %d, iteration %d" % (opt_params['cur_step'] + 1,
                                                  opt_params['optiter'] + 1))

    from .parallel_backends import _BACKEND, JoblibBackend
    if _BACKEND is None:
        _BACKEND = JoblibBackend(n_jobs=1)

    # set parameters
    for param_name, test_value in zip(opt_params['ranges'].keys(), new_params):
        params[param_name] = test_value

    # run the simulation, but stop early if possible
    params['tstop'] = opt_params['opt_end']
    net = Network(params, add_drives_from_params=True)
    dpls = _BACKEND.simulate(net, n_trials=1)
    # avg_dpl = average_dipoles(dpls)
    avg_dpl = dpls[0].copy()
    avg_rmse = rmse(avg_dpl, opt_dpls['exp_dpl'],
                    tstart=opt_params['opt_start'],
                    tstop=opt_params['opt_end'],
                    weights=opt_params['weights'])

    if avg_rmse < opt_params['stepminopterr']:
        best = "[best] "
        opt_params['stepminopterr'] = avg_rmse
        opt_dpls['best_dpl'] = avg_dpl
    else:
        best = ""

    print("%sweighted RMSE: %.2f over range [%3.3f-%3.3f] ms" %
          (best, avg_rmse, opt_params['opt_start'], opt_params['opt_end']))

    opt_params['optiter'] += 1

    return avg_rmse  # nlopt expects error


def run_optimization(maxiter, param_ranges, optrun):

    cons = list()
    x0 = list()
    for idx, param_name in enumerate(param_ranges):
        x0.append(param_ranges[param_name]['initial'])
        cons.append(
            lambda x: x[idx] <= param_ranges[param_name]['maxval'])
        cons.append(
            lambda x: x[idx] >= param_ranges[param_name]['minval'])

    result = fmin_cobyla(optrun, cons=cons, rhoend=1e-4,
                         x0=x0, maxfun=maxiter)
    return result


def optimize_evoked(params, exp_dpl, maxiter,
                    timing_range_multiplier=3.0, sigma_range_multiplier=50.0,
                    synweight_range_multiplier=500.0, decay_multiplier=1.6):
    """Optimize drives to generate evoked response.

    Parameters
    ----------
    params : dict
        The initial params
    exp_dpl : instance of Dipole
        The target experimental dipole.
    maxiter : int
        The maximum number of of simulations to run for optimizing
        one "chunk".
    timing_range_multiplier : float
        The scale of timing values to sweep over.
    sigma_range_multiplier : float
        The scale of sigma values to sweep over.
    synweight_range_multiplier : float
        The scale of input synaptic weights to sweep over.
    decay_multiplier : float
        The decay multiplier.

    Returns
    -------
    params : dict
        The optimized params dictionary.
    """
    from .parallel_backends import _BACKEND, JoblibBackend

    if _BACKEND is None:
        _BACKEND = JoblibBackend(n_jobs=1)

    print("Running simulation with initial parameters")
    net = Network(params, add_drives_from_params=True)
    initial_dpl = _BACKEND.simulate(net, n_trials=1)

    # Create a sorted dictionary with the inputs and parameters
    # belonging to each.
    # Then, calculate the appropriate weight function to be used
    # in RMSE using the CDF defined by the input timing's mean and
    # std. deviation parameters.
    # Lastly, use the weights to define non-overlapping "chunks" of
    # the simulation timeframe to optimize. Chunks are consolidated if
    # more than one input should
    # be optimized at a time.
    evinput_params = split_by_evinput(params, sigma_range_multiplier,
                                      timing_range_multiplier,
                                      synweight_range_multiplier)
    evinput_params = _generate_weights(evinput_params, params,
                                       decay_multiplier)
    param_chunks = consolidate_chunks(evinput_params)

    avg_rmse = rmse(initial_dpl[0], exp_dpl, tstop=params['tstop'])
    print("Initial RMSE: %.2f" % avg_rmse)

    opt_params = dict()
    for step in range(len(param_chunks)):
        opt_params['cur_step'] = step
        total_steps = len(param_chunks)

        # param_chunks is the optimization information for all steps.
        # opt_params is a pointer to the params for each step
        opt_params.update(param_chunks[step])

        if maxiter == 0:
            print("Skipping optimization step %d (0 simulations)" % (step + 1))
            continue

        if (opt_params['cur_step'] > 0 and
                opt_params['cur_step'] == total_steps - 1):
            # update currently used params
            for var_name in opt_params['ranges']:
                params[var_name] = opt_params['ranges'][var_name]['initial']

            # For the last step (all inputs), recalculate ranges and update
            # param_chunks. If previous optimization steps increased
            # std. dev. this could result in fewer optimization steps as
            # inputs may be deemed too close together and be grouped in a
            # single optimization step.
            #
            # The purpose of the last step (with regular RMSE) is to clean up
            # overfitting introduced by local weighted RMSE optimization.

            evinput_params = split_by_evinput(params, sigma_range_multiplier,
                                              timing_range_multiplier,
                                              synweight_range_multiplier)
            evinput_params = _generate_weights(evinput_params, params,
                                               decay_multiplier)
            param_chunks = consolidate_chunks(evinput_params)

            # reload opt_params for the last step in case the number of
            # steps was changed by updateoptparams()
            opt_params.update(param_chunks[total_steps - 1])

        print("Starting optimization step %d/%d" % (step + 1, total_steps))

        opt_params['optiter'] = 0
        opt_params['stepminopterr'] = 1e9  # min optimization error so far
        opt_dpls = dict(best_dpl=None, exp_dpl=exp_dpl)

        def _optrun(new_params):
            return optrun(new_params, opt_params,
                          params, opt_dpls=opt_dpls)

        print('Optimizing from [%3.3f-%3.3f] ms' % (opt_params['opt_start'],
                                                    opt_params['opt_end']))
        opt_results = run_optimization(maxiter=maxiter,
                                       param_ranges=opt_params['ranges'],
                                       optrun=_optrun)

        # update opt_params for the next round
        for var_name, value in zip(opt_params['ranges'], opt_results):
            opt_params['ranges'][var_name]['initial'] = value

    # save the optimized params
    for var_name in opt_params['ranges']:
        params[var_name] = opt_params['ranges'][var_name]['initial']

    avg_rmse = rmse(opt_dpls['best_dpl'], exp_dpl, tstop=params['tstop'])
    print("Final RMSE: %.2f" % avg_rmse)
    return params
