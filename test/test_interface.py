import numpy as np
import h5py
import os
import warnings

import surfinBH

# Ignore these specific warnings.
warnings.filterwarnings("ignore", message="Mass ratio outside training range.")
warnings.filterwarnings("ignore", message="Spin magnitude of BhA outside"
        " training range.")
warnings.filterwarnings("ignore", message="Spin magnitude of BhB outside"
        " training range.")
warnings.filterwarnings("ignore", message="Extrapolating dynamics to")

def single_kwargs_test(fit, num_tests, kwargs={}):

    # get range of params
    param_lims = fit.soft_param_lims

    # get bool for whether aligned_spin_only
    aligned_spin_only = fit.aligned_spin_only

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

        # Check evaluation of different call modules

        # remnant mass
        mC_and_err = fit.mC(q, chiA, chiB, **kwargs)

        # remnant spin
        chiC_and_err = fit.chiC(q, chiA, chiB, **kwargs)

        # remnant kick
        velC_and_err = fit.velC(q, chiA, chiB, **kwargs)

        # all together
        #mC, chiC, velC, mC_err_est, chiC_err_est, velC_err_est \
        all_data = fit.all(q, chiA, chiB, **kwargs)

        # Check that all_data has the right ordering
        rtol = 1e-11

        if hasattr(mC_and_err, '__len__'):      # has both fit val and err_est
            np.testing.assert_allclose(mC_and_err[0], all_data[0], rtol=rtol)
            np.testing.assert_allclose(mC_and_err[1], all_data[3], rtol=rtol)
        else:
            np.testing.assert_allclose(mC_and_err, all_data[0], rtol=rtol)

        if hasattr(chiC_and_err, '__len__'):      # has both fit val and err_est
            np.testing.assert_allclose(chiC_and_err[0], all_data[1], rtol=rtol)
            np.testing.assert_allclose(chiC_and_err[1], all_data[4], rtol=rtol)
        else:
            np.testing.assert_allclose(chiC_and_err, all_data[1], rtol=rtol)

        if hasattr(velC_and_err, '__len__'):      # has both fit val and err_est
            np.testing.assert_allclose(velC_and_err[0], all_data[2], rtol=rtol)
            np.testing.assert_allclose(velC_and_err[1], all_data[5], rtol=rtol)
        else:
            np.testing.assert_allclose(velC_and_err, all_data[2], rtol=rtol)


def test_interface():
    """ Tests that the call modules are evaluating without breaking.
        Also tests that fit.all is doing the right thing.
    """
    # List of all available fits
    fit_names = surfinBH.fits_collection.keys()
    for name in fit_names:
        short_name = name.split('surfinBH')[-1]

        # Load fit
        fit = surfinBH.LoadFits(name)

        num_tests = 3

        # Test without any kwargs
        single_kwargs_test(fit, num_tests)

        # Test with different optional kwargs
        # NOTE: this would be empty unless _extra_regression_kwargs() is
        # overridden
        extra_kwargs_list = fit._extra_regression_kwargs()
        for i in range(len(extra_kwargs_list)):
            kwargs = extra_kwargs_list[i]
            single_kwargs_test(fit, num_tests, kwargs=kwargs)
