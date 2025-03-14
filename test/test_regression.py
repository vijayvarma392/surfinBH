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

# lower rtol seems to make the test fail for the GPR error
# estimates when testing on a different machine. This
# seems to happen even with the same version of sklearn.
rtol = 1e-6

def test_fit_regression():
    """ Compares all existing fits against saved regression data.
        Regression data should already have been generated by
        generate_regression_data.py
    """
    # List of all available fits
    fit_names = surfinBH.fits_collection.keys()
    for name in fit_names:

        if name == 'NRSur7dq4EmriRemnant':
            # FIXME Somehow, the error estimate for this fit seems particularly
            # finnicky accross different machines. For now, essentially
            # hacking it.
            atol = 1e-3
        else:
            atol = 0

        # allow for both naming formats surfinBH7dq2 and NRSur7dq4Remnant
        if 'surfinBH' in name:
            name_tag = name.split('surfinBH')[-1]
        else:
            name_tag = name.split('NRSur')[-1].split('Remnant')[0]

        # Load fit
        fit = surfinBH.LoadFits(name)

        # Load regression data
        regression_h5file = h5py.File('test/regression_data/fit_%s.h5'%(
                name_tag), 'r')

        extra_kwargs_list = fit._extra_regression_kwargs()

        # Each group in the h5file can have different kwargs to test
        kwargs_grp_keys = regression_h5file.keys()
        for kw_grp_key in kwargs_grp_keys:
            kw_h5grp = regression_h5file[kw_grp_key]

            print('\nrunning %s'%kw_grp_key)

            if kw_grp_key == 'No_kwargs':
                kwargs = {}
            else:
                kwargs = extra_kwargs_list[int(kw_grp_key.split('_set_')[-1])]

            # Compare fit and regression data for each regression test
            test_keys = kw_h5grp.keys()
            for test in test_keys:
                print('... running %s'%test)
                test_h5grp = kw_h5grp[test]
                q = test_h5grp['q'][()]
                chiA = test_h5grp['chiA'][()]
                chiB = test_h5grp['chiB'][()]
                
                if name == 'NRSur3dq8BMSRemnant':
                    #compute fits
                    alpha, boost,  alpha_err, boost_err = fit.all(q, chiA, chiB, **kwargs)
                    
                    # supertranslation
                    y_reg = test_h5grp['y/alpha'][()]
                    y_fit = alpha, alpha_err
                    np.testing.assert_allclose(y_fit, y_reg, rtol=rtol, atol=atol)

                    # boost velocity
                    y_reg = test_h5grp['y/boost'][()]
                    y_fit = boost, boost_err
                    np.testing.assert_allclose(y_fit, y_reg, rtol=rtol, atol=atol)
                    
                else:
                    # remnant mass
                    y_reg = test_h5grp['y/mf'][()]
                    y_fit = fit.mf(q, chiA, chiB, **kwargs)
                    np.testing.assert_allclose(y_fit, y_reg, rtol=rtol, atol=atol)
    
                    # remnant spin
                    y_reg = test_h5grp['y/chif'][()]
                    y_fit = fit.chif(q, chiA, chiB, **kwargs)
                    np.testing.assert_allclose(y_fit, y_reg, rtol=rtol, atol=atol)
                    
                    # remnant kick
                    # Needed for NRSur7dq4EmriRemnant
                    if 'vf' in test_h5grp['y'].keys():
                        y_reg = test_h5grp['y/vf'][()]
                        y_fit = fit.vf(q, chiA, chiB, **kwargs)
                        np.testing.assert_allclose(y_fit, y_reg, rtol=rtol, atol=atol)
   
                     
        regression_h5file.close()
