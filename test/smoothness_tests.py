import numpy as np
import h5py
import argparse
import os
import time
import matplotlib.pyplot as P

import surfinBH


# Number of data points to generate for each 1d plot
NUM_PARAMS = 50

# Number of 1d plots to generate per fit
NUM_TESTS = 10

# plot settings
marker_size=100
marker_size_star=9
label_fontsize = 14
title_fontsize = 16
ticks_fontsize = 16
line_width = 1.5
legend_size = 12


# --------------------------------------------------------------------------
def test_smoothness_wrapper(ax_pair, x_list, function, ylabel, \
        params_list, label, y_index=None):
    """ Wrapper to make plots of fit errors along different 1d directions.
        x_list is an array of param values along the 1d direction.
    """

    # Get list of q, chiA, and chiB
    if aligned_spin_only:
        q_list, chiAz_list, chiBz_list = params_list
        chiA_list = [[0,0,tmp] for tmp in chiAz_list]
        chiB_list = [[0,0,tmp] for tmp in chiBz_list]
    else:
        q_list, chiA_mag_list, chiA_th_list, chiA_ph_list, \
            chiB_mag_list, chiB_th_list, chiB_ph_list = params_list

        chiA_list = np.array(\
                [chiA_mag_list*np.sin(chiA_th_list)*np.cos(chiA_ph_list),
                 chiA_mag_list*np.sin(chiA_th_list)*np.sin(chiA_ph_list),
                 chiA_mag_list*np.cos(chiA_th_list)]\
                        ).T

        chiB_list = np.array(\
                [chiB_mag_list*np.sin(chiB_th_list)*np.cos(chiB_ph_list),
                 chiB_mag_list*np.sin(chiB_th_list)*np.sin(chiB_ph_list),
                 chiB_mag_list*np.cos(chiB_th_list)]\
                        ).T

    # Plot function and error estimates in a row of subplots
    y_list = []
    y_err_list = []
    for i in range(len(q_list)):
        y, y_err = function(q_list[i], chiA_list[i], chiB_list[i])
        if y_index is None:
            y_list.append(y)
            y_err_list.append(y_err)
        elif type(y_index) == int:
            y_list.append(y[y_index])
            y_err_list.append(y_err[y_index])
        elif y_index == 'magnitude':
            y_list.append(np.sqrt(np.sum(y**2)))
            y_err_list.append(np.sqrt(np.sum(y_err**2)))

    ax_pair[0].plot(x_list, y_list, label=label)
    ax_pair[1].semilogy(x_list, y_err_list, label=label)

    ax_pair[0].set_ylabel('%s'%ylabel, fontsize=label_fontsize)
    ax_pair[1].set_ylabel('$\Delta$%s'%ylabel, fontsize=label_fontsize)



# --------------------------------------------------------------------------
def generate_params_along_line_aligned(x_param):
    """ Generates params along a 1d line for aligned-spin params.
        The given x_param is varied but all other params are kept const at
        some randomly chosen values.
        For each of the aligned-spin params:
            if x_param is this particular param, returns a uniform list of
                values in the range of this param. Also saves this list as
                x_list.
            else, returns a list, each element of which is the same, a
                randomly chosen value in the range of this param.
        Also returns a label with all the fixed params.
    """

    label = ''

    if x_param == 'q':
        q_list = np.linspace(1, param_lims['q'], NUM_PARAMS)
        x_list = q_list
    else:
        q = np.random.uniform(1, param_lims['q'])
        q_list = [q for i in range(NUM_PARAMS)]
        label += '$q=%.2f$ '%q

    if x_param == 'chi1':
        chiAz_list = np.linspace(-param_lims['chiAmag'], \
            param_lims['chiAmag'], NUM_PARAMS)
        x_list = chiAz_list
    else:
        chiAz = np.random.uniform(-param_lims['chiAmag'], param_lims['chiAmag'])
        chiAz_list = [chiAz for i in range(NUM_PARAMS)]
        label += '$\chi_{1z}=%.2f$ '%chiAz

    if x_param == 'chi2':
        chiBz_list = np.linspace(-param_lims['chiBmag'], \
            param_lims['chiBmag'], NUM_PARAMS)
        x_list = chiBz_list
    else:
        chiBz = np.random.uniform(-param_lims['chiBmag'], param_lims['chiBmag'])
        chiBz_list = [chiBz for i in range(NUM_PARAMS)]
        label += '$\chi_{2z}=%.2f$ '%chiBz

    params_list = [q_list, chiAz_list, chiBz_list]

    return x_list, label, params_list


# --------------------------------------------------------------------------
def generate_params_along_line(x_param):
    """ Generates params along a 1d line for 7d params.
        The given x_param is varied but all other params are kept const at
        some randomly chosen values.
        For each of the 7d params:
            if x_param is this particular param, returns a uniform list of
                values in the range of this param. Also saves this list as
                x_list.
            else, returns a list, each element of which is the same, a
                randomly chosen value in the range of this param.
        Also returns a label with all the fixed params.
    """

    label = ''

    if x_param == 'q':
        q_list = np.linspace(1, param_lims['q'], NUM_PARAMS)
        x_list = q_list
    else:
        q = np.random.uniform(1, param_lims['q'])
        q_list = [q for i in range(NUM_PARAMS)]
        label += '$q=%.2f$ '%q

    if x_param == 'chi1':
        chiA_mag_list = np.linspace(0, param_lims['chiAmag'], NUM_PARAMS)
        x_list = chiA_mag_list
    else:
        chiA_mag = np.random.uniform(0, param_lims['chiAmag'])
        chiA_mag_list = [chiA_mag for i in range(NUM_PARAMS)]
        label += '$\chi_1=%.2f$ '%chiA_mag

    if x_param == 'chi1_th':
        chiA_th_list = np.linspace(0, np.pi, NUM_PARAMS)
        x_list = chiA_th_list
    else:
        chiA_th = np.random.uniform(0, np.pi)
        chiA_th_list = [chiA_th for i in range(NUM_PARAMS)]
        label += '$\chi_{1\\theta}=%.2f$ '%chiA_th

    if x_param == 'chi1_ph':
        chiA_ph_list = np.linspace(0, 2*np.pi, NUM_PARAMS)
        x_list = chiA_ph_list
    else:
        chiA_ph = np.random.uniform(0, 2*np.pi)
        chiA_ph_list = [chiA_ph for i in range(NUM_PARAMS)]
        label += '$\chi_{1\\phi}=%.2f$ '%chiA_ph

    if x_param == 'chi2':
        chiB_mag_list = np.linspace(0, param_lims['chiBmag'], NUM_PARAMS)
        x_list = chiB_mag_list
    else:
        chiB_mag = np.random.uniform(0, param_lims['chiBmag'])
        chiB_mag_list = [chiB_mag for i in range(NUM_PARAMS)]
        label += '$\chi_2=%.2f$ '%chiB_mag

    if x_param == 'chi2_th':
        chiB_th_list = np.linspace(0, np.pi, NUM_PARAMS)
        x_list = chiB_th_list
    else:
        chiB_th = np.random.uniform(0, np.pi)
        chiB_th_list = [chiB_th for i in range(NUM_PARAMS)]
        label += '$\chi_{2\\theta}=%.2f$ '%chiB_th

    if x_param == 'chi2_ph':
        chiB_ph_list = np.linspace(0, 2*np.pi, NUM_PARAMS)
        x_list = chiB_ph_list
    else:
        chiB_ph = np.random.uniform(0, 2*np.pi)
        chiB_ph_list = [chiB_ph for i in range(NUM_PARAMS)]
        label += '$\chi_{2\\phi}=%.2f$ '%chiB_ph

    params_list = [q_list, chiA_mag_list, chiA_th_list, chiA_ph_list, \
        chiB_mag_list, chiB_th_list, chiB_ph_list]

    return x_list, label, params_list


# --------------------------------------------------------------------------
def test_smoothness(x_param, x_param_label):
    """ Tests smoothness in the direction of x_param.
    Does NUM_TESTS number of tests, for each tests the rest of the
    params are fixed at some randomly chosen values.
    """

    if aligned_spin_only:
        # Don't need spin direction plots for aligned-spin fits
        if x_param[-3:] in ['_th', '_ph']:
            return
        fig, axarr = P.subplots(4,2,figsize=(10,10))
    else:
        fig, axarr = P.subplots(7,2,figsize=(10,15))

    P.subplots_adjust(hspace=0.25, wspace=0.35)
    axarr = axarr.reshape(-1, order='C')

    comp_labels = ['x', 'y', 'z']
    for i in range(NUM_TESTS):

        # Get parameters along a 1d line, where x_param is varied but
        # all other params are kept const at some randomly chosen values
        if aligned_spin_only:
            x_list, label, params_list \
                = generate_params_along_line_aligned(x_param)
        else:
            x_list, label, params_list = generate_params_along_line(x_param)

        # final mass plots along the 1d line
        test_smoothness_wrapper(axarr[:2], x_list, fit.mC, '$m$', params_list, \
            label)

        # final spin, but plot only z values for aligned-spins
        chi_idx = 0
        for idx in range(3):
            if (idx == 2) or (not aligned_spin_only):
                test_smoothness_wrapper(\
                    axarr[2+2*chi_idx:4+2*chi_idx], x_list, fit.chiC, \
                    '$\chi_{%s}$'%comp_labels[idx], params_list, label, \
                    y_index=idx)
                chi_idx += 1

        # final kick, but plot only x,y values for aligned-spins
        vel_idx = 0
        for idx in range(3):
            if (idx in [0,1]) or (not aligned_spin_only):
                test_smoothness_wrapper(\
                    axarr[2+2*chi_idx+2*vel_idx:4+2*chi_idx+2*vel_idx], \
                    x_list, fit.velC, '$v_{%s}$'%comp_labels[idx], \
                    params_list, label, y_index=idx)
                vel_idx += 1


    axarr[-1].set_xlabel(x_param_label, fontsize=label_fontsize)
    axarr[-2].set_xlabel(x_param_label, fontsize=label_fontsize)

    axarr[1].legend(loc=(1.1,-1))
    P.savefig('%s/%s_smoothness_1d_%s.png'%(outdir, args.name, x_param), \
        bbox_inches='tight')
    P.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Tests smoothness of fit ' \
        'along different 1d directions', \
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--name", "-n", type=str, required=True, \
        help="Fit name without the surfinBH prefix. Eg. 7dq2.")
    args = parser.parse_args()

    # Load fit
    fit_name = 'surfinBH%s'%args.name
    fit = surfinBH.LoadFits(fit_name)

    # get allowed range of params
    param_lims = fit.hard_param_lims

    # get bool for whether aligned_spin_only
    aligned_spin_only = fit.aligned_spin_only

    outdir = 'smoothness_tests'
    os.system('mkdir -p %s'%outdir)

    # test smoothness along different 1d directions
    test_smoothness('q', '$q$')
    test_smoothness('chi1', '$\chi_{1}$')
    test_smoothness('chi1_th', '$\chi_{1\\theta}$')
    test_smoothness('chi1_ph', '$\chi_{1\\phi}$')
    test_smoothness('chi2', '$\chi_{2}$')
    test_smoothness('chi2_th', '$\chi_{2\\theta}$')
    test_smoothness('chi2_ph', '$\chi_{2\\phi}$')
