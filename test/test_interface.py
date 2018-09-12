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

    for i in range(num_tests):
        # Generate params randomly within allowed values
        q, chiA, chiB = fit._generate_random_params_for_tests() 

        # Check evaluation of different call modules

        # remnant mass
        mf_and_err = fit.mf(q, chiA, chiB, **kwargs)

        # remnant spin
        chif_and_err = fit.chif(q, chiA, chiB, **kwargs)

        # remnant kick
        vf_and_err = fit.vf(q, chiA, chiB, **kwargs)

        # all together
        #mf, chif, vf, mf_err_est, chif_err_est, vf_err_est \
        all_data = fit.all(q, chiA, chiB, **kwargs)

        # Check that all_data has the right ordering
        rtol = 1e-11

        if hasattr(mf_and_err, '__len__'):      # has both fit val and err_est
            np.testing.assert_allclose(mf_and_err[0], all_data[0], rtol=rtol)
            np.testing.assert_allclose(mf_and_err[1], all_data[3], rtol=rtol)
        else:
            np.testing.assert_allclose(mf_and_err, all_data[0], rtol=rtol)

        if hasattr(chif_and_err, '__len__'):      # has both fit val and err_est
            np.testing.assert_allclose(chif_and_err[0], all_data[1], rtol=rtol)
            np.testing.assert_allclose(chif_and_err[1], all_data[4], rtol=rtol)
        else:
            np.testing.assert_allclose(chif_and_err, all_data[1], rtol=rtol)

        if hasattr(vf_and_err, '__len__'):      # has both fit val and err_est
            np.testing.assert_allclose(vf_and_err[0], all_data[2], rtol=rtol)
            np.testing.assert_allclose(vf_and_err[1], all_data[5], rtol=rtol)
        else:
            np.testing.assert_allclose(vf_and_err, all_data[2], rtol=rtol)


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
