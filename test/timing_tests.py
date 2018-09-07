import numpy as np
import h5py
import argparse
import os
import time

import surfinBH


def eval_fits(fit, param_lims, function):

    start_time = time.time()
    for i in range(num_tests):
        # Generate params randomly within allowed values
        q = np.random.uniform(1, param_lims['q'])
        chiAmag= np.random.uniform(0, param_lims['chiAmag'])
        chiBmag= np.random.uniform(0, param_lims['chiBmag'])
        if aligned_spin_only:
            chiAth, chiBth, chiAph, chiBph = 0,0,0,0
        else:
            chiAth = np.random.uniform(0, np.pi)
            chiBth = np.random.uniform(0, np.pi)
            chiAph = np.random.uniform(0, 2*np.pi)
            chiBph = np.random.uniform(0, 2*np.pi)

        chiA = [chiAmag*np.sin(chiAth)*np.cos(chiAph),
                chiAmag*np.sin(chiAth)*np.sin(chiAph),
                chiAmag*np.cos(chiAth)]

        chiB = [chiBmag*np.sin(chiBth)*np.cos(chiBph),
                chiBmag*np.sin(chiBth)*np.sin(chiBph),
                chiBmag*np.cos(chiBth)]

        # evaluate fit
        y = function(q, chiA, chiB)

    end_time = time.time()
    return (end_time - start_time)/num_tests


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Tests evaluation times for'
            ' fits.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--name", "-n", type=str, required=True, \
            help="Fit name without the surfinBH prefix. Eg. 7dq2.")
    args = parser.parse_args()

    # Load fit
    fit_name = 'surfinBH%s'%args.name
    fit = surfinBH.LoadFits(fit_name)

    # number of test evaluations per fit
    num_tests = 1000

    # get allowed range of params
    param_lims = fit.hard_param_lims

    # get bool for whether aligned_spin_only
    aligned_spin_only = fit.aligned_spin_only

    # remnant mass
    eval_time = eval_fits(fit, param_lims, fit.mC)
    print('Time per evaluation of mC: %.3e s'%eval_time)

    # remnant spin
    eval_time = eval_fits(fit, param_lims, fit.chiC)
    print('Time per evaluation of chiC: %.3e s'%eval_time)

    # remnant kick
    eval_time = eval_fits(fit, param_lims, fit.velC)
    print('Time per evaluation of velC: %.3e s'%eval_time)
