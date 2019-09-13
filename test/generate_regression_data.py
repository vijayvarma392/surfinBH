import numpy as np
import h5py
import argparse
import os

import surfinBH

def save_data(h5grp, fit, num_tests, kwargs={}):

    for i in range(num_tests):
        # save each test evaluation as group
        test_h5grp = h5grp.create_group('test_%d'%i)

        # Generate params randomly within allowed values
        q, chiA, chiB = fit._generate_random_params_for_tests()

        # save params
        test_h5grp.create_dataset('q', data=q)
        test_h5grp.create_dataset('chiA', data=chiA)
        test_h5grp.create_dataset('chiB', data=chiB)

        # save evaluations as a group
        y_h5grp = test_h5grp.create_group('y')

        # remnant mass
        y = fit.mf(q, chiA, chiB, **kwargs)
        y_h5grp.create_dataset('mf', data=y)

        # remnant spin
        y = fit.chif(q, chiA, chiB, **kwargs)
        y_h5grp.create_dataset('chif', data=y)

        # remnant kick
        y = fit.vf(q, chiA, chiB, **kwargs)
        y_h5grp.create_dataset('vf', data=y)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate regression data for'
            ' fits.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--name", "-n", type=str, required=True, \
            help="Fit name. Eg. NRSur7dq4Remnant or surfin7dq2.")
    args = parser.parse_args()

    # allow for both naming formats surfinBH7dq2 and NRSur7dq4Remnant
    if 'surfinBH' in args.name:
        name_tag = args.name.split('surfinBH')[-1]
    else:
        name_tag = args.name.split('NRSur')[-1].split('Remnant')[0]

    # If it gets overwritten, we always have git
    h5file = h5py.File('test/regression_data/fit_%s.h5'%name_tag, 'w')

    # Load fit, assuming fit data is already present in surfinBH_data
    fit = surfinBH.LoadFits(args.name)

    # number of test evaluations per fit
    num_tests = 10

    # save regression without optional kwargs
    h5grp = h5file.create_group('No_kwargs')
    save_data(h5grp, fit, num_tests)

    # save regression with optional kwargs.
    # NOTE: this would be empty unless _extra_regression_kwargs() is
    # overridden
    extra_kwargs_list = fit._extra_regression_kwargs()
    for i in range(len(extra_kwargs_list)):
        kwargs = extra_kwargs_list[i]
        h5grp = h5file.create_group('kwargs_set_%d'%i)
        save_data(h5grp, fit, num_tests, kwargs=kwargs)

    h5file.close()
