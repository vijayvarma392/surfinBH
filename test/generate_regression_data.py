import numpy as np
import h5py
import argparse
import os

import surfinBH

def save_data(h5grp, fit, num_tests, kwargs={}):

    # get allowed range of params
    param_lims = fit.hard_param_lims

    # get bool for whether aligned_spin_only
    aligned_spin_only = fit.aligned_spin_only

    for i in range(num_tests):
        # save each test evaluation as group
        test_h5grp = h5grp.create_group('test_%d'%i)

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


        # save params
        test_h5grp.create_dataset('q', data=q)
        test_h5grp.create_dataset('chiA', data=chiA)
        test_h5grp.create_dataset('chiB', data=chiB)

        # save evaluations as a group
        y_h5grp = test_h5grp.create_group('y')

        # remnant mass
        y = fit.mC(q, chiA, chiB, **kwargs)
        y_h5grp.create_dataset('mC', data=y)

        # remnant spin
        y = fit.chiC(q, chiA, chiB, **kwargs)
        y_h5grp.create_dataset('chiC', data=y)

        # remnant kick
        y = fit.velC(q, chiA, chiB, **kwargs)
        y_h5grp.create_dataset('velC', data=y)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate regression data for'
            ' fits.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--name", "-n", type=str, required=True, \
            help="Fit name without the surfinBH prefix. Eg. 7dq2.")
    args = parser.parse_args()

    # If it gets overwritten, we always have git
    h5file = h5py.File('test/regression_data/fit_%s.h5'%args.name, 'w')

    # Load fit, assuming fit data is already present in surfinBH_data
    fit_name = 'surfinBH%s'%args.name
    fit = surfinBH.LoadFits(fit_name)

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
