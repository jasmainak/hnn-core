"""Objective functions for parameter optimization."""

# Authors: Carolina Fernandez <cxf418@miami.edu>
#          Nick Tolley <nicholas_tolley@brown.edu>
#          Ryan Thorpe <ryan_thorpe@brown.edu>
#          Mainak Jas <mjas@mgh.harvard.edu>

from hnn_core import simulate_dipole
from ..dipole import _rmse


def _rmse_evoked(initial_net, initial_params, set_params, predicted_params,
                 update_params, obj_values, scale_factor, smooth_window_len,
                 tstop, **kwargs):
    """The objective function for evoked responses.

    Parameters
    ----------
    initial_net : instance of Network
        The network object.
    initial_params : dict
        Keys are parameter names, values are initial parameters.
    set_params : func
        User-defined function that sets network drives and parameters.
    predicted_params : list
        Parameters selected by the optimizer.
    update_params : func
        Function to update params.
    scale_factor : float
        The dipole scale factor.
    smooth_window_len : float
        The smooth window length.
    tstop : float
        The simulated dipole's duration.
    target : instance of Dipole
        A dipole object with experimental data.

    Returns
    -------
    obj : float
        Normalized RMSE between recorded and simulated dipole.
    """

    params = update_params(initial_params, predicted_params)

    # simulate dpl with predicted params
    new_net = initial_net.copy()
    set_params(new_net, params)
    dpl = simulate_dipole(new_net, tstop=tstop, n_trials=1)[0]

    # smooth & scale
    dpl.scale(scale_factor)
    if smooth_window_len is not None:
        dpl.smooth(smooth_window_len)

    obj = _rmse(dpl, kwargs['target'], tstop=tstop)

    obj_values.append(obj)

    return obj


def _maximize_psd(initial_net, initial_params, set_params, predicted_params,
                  update_params, obj_values, scale_factor, smooth_window_len,
                  tstop, **kwargs):
    """The objective function for evoked responses.

    Parameters
    ----------
    initial_net : instance of Network
        The network object.
    initial_params : dict
        Keys are parameter names, values are initial parameters.
    set_params : func
        User-defined function that sets network drives and parameters.
    predicted_params : list
        Parameters selected by the optimizer.
    update_params : func
        Function to update params.
    scale_factor : float
        The dipole scale factor.
    smooth_window_len : float
        The smooth window length.
    tstop : float
        The simulated dipole's duration.
    f_bands : list of tuples
        Lower and higher limit for each frequency band.
    relative_bandpower : tuple
        Weight for each frequency band.

    Returns
    -------
    obj : float
        Normalized RMSE between recorded and simulated dipole.
    """

    import numpy as np

    from scipy.signal import periodogram

    params = update_params(initial_params, predicted_params)

    # simulate dpl with predicted params
    new_net = initial_net.copy()
    set_params(new_net, params)
    dpl = simulate_dipole(new_net, tstop=tstop, n_trials=1)[0]

    # smooth & scale
    dpl.scale(scale_factor)
    if smooth_window_len is not None:
        dpl.smooth(smooth_window_len)

    # resample?

    # get psd of simulated dpl
    freqs_simulated, psd_simulated = periodogram(dpl.data['agg'], dpl.sfreq,
                                                 window='hamming')

    # for each f band
    f_bands_psds = list()
    for idx, f_band in enumerate(kwargs['f_bands']):
        f_band_idx = np.where(np.logical_and(freqs_simulated >= f_band[0],
                                             freqs_simulated <= f_band[1]))[0]
        f_bands_psds.append((-kwargs['relative_bandpower'][idx] * sum(
            psd_simulated[f_band_idx])) / sum(psd_simulated))

    # grand sum
    obj = sum(f_bands_psds)

    obj_values.append(obj)

    return obj
