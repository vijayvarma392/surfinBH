import numpy as np
import h5py
import argparse
import os

import surfinBH

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate regression data for'
            ' fits.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--name", "-n", type=str, required=True, \
            help="Fit name without the surfinBH prefix. Eg. 7dq2.")
    args = parser.parse_args()

    # If it get's overwritten, we always have git
    h5file = h5py.File('test/regression_data/fit_%s.h5'%args.name, 'w')

    # Load fit, assuming fit data is already present in surfinBH_data
    fit_keys = ['mC', 'chiC', 'velC']
    fit_name = 'surfinBH%s'%args.name
    fit = surfinBH.LoadFits(fit_name)

    # valid range for params
    param_lims = fit.hard_param_lims

    # number of test evaluations per fit
    num_tests = 10
    for i in range(num_tests):
        # save each test evaluation as group
        test_h5grp = h5file.create_group('test_%d'%i)
        x = []
        for lims in param_lims:
            x.append(np.random.uniform(lims[0], lims[1]))

        # save param
        test_h5grp.create_dataset('x', data=x)

        # save evaluations as a group
        y_h5grp = test_h5grp.create_group('y')

        for key in fit_keys:
            y = fit(key, x)
            # save fit evaluation
            y_h5grp.create_dataset(key, data=y)

    h5file.close()
