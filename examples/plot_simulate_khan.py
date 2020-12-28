"""
==========================
Simulate steady state beta
==========================
"""

# Authors: Mainak Jas <mjas@mgh.harvard.edu>
#          Sheraz Khan

import os.path as op
import tempfile

###############################################################################
# Let us import hnn_core

import hnn_core
from hnn_core import simulate_dipole, read_params, Network, read_spikes
from hnn_core.viz import plot_dipole

hnn_core_root = op.dirname(hnn_core.__file__)

###############################################################################
# Then we read the parameters file
params_fname = op.join(hnn_core_root, 'param', 'default.json')
params = read_params(params_fname)
params['tstop'] = 1000

###############################################################################
# Let us first create our network from the params file and visualize the cells
# inside it. The default behaviour of Network is to add and instantiate six
# 'default' drives, but we will suppress that by setting the
# `initialise_hnn_drives`-argument to `False`.
net = Network(params, initialise_hnn_drives=False)

###############################################################################
# The network of cells is now defined, to which we add external drives as
# required. Weights are prescribed separately for AMPA and NMDA receptors
# (receptors that are not used can be omitted or set to zero)

weights_ampa = {'L2_basket': 0.006562, 'L2_pyramidal': .000007,
                'L5_pyramidal': 0.142300}
weights_nmda = {'L2_basket': 0.019482, 'L2_pyramidal': 0.004317,
                'L5_pyramidal': 0.080074}

t0 = 0
numspikes = 2
location = 'proximal'
net.add_bursty_drive(
        'bursty1', distribution='normal', t0=t0, sigma_t0=2.,
        T=params['tstop'], burst_f=25, burst_sigma_f=0.,
        numspikes=numspikes, repeats=1,
        weights_ampa=weights_ampa, weights_nmda=weights_nmda,
        location=location, seedcore=2)

t0 = 14
# XXX: original paper says 3, also 7 ms instead of 10 ms ISI
numspikes = 2
location = 'distal'
net.add_bursty_drive(
        'bursty2', distribution='normal', t0=t0, sigma_t0=2.,
        T=params['tstop'], burst_f=25, burst_sigma_f=0., numspikes=numspikes,
        repeats=1, weights_ampa=weights_ampa, weights_nmda=weights_nmda,
        location=location, seedcore=2)

dpls = simulate_dipole(net, n_trials=1, postproc=True)

###############################################################################
# and then plot it
import matplotlib.pyplot as plt
fig, axes = plt.subplots(2, 1, sharex=True, figsize=(6, 6))
plot_dipole(dpls, ax=axes[0], layer='agg', show=False)
net.cell_response.plot_spikes_hist(ax=axes[1],
                                   spike_types=['evprox', 'evdist'])
net.cell_response.plot_spikes_raster()
