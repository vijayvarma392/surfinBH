import numpy as np
from surfinBH import surfinBH

#=============================================================================
class Fit3dq8(surfinBH.SurFinBH):
    """
A class for the 3dq8 surrogare final BH fits presented in Varma et al. 2018,
in prep. This model predicts the final mass mC, final spin chiC and final
kick velocity velC, for the remnants of nonprecessing binary black hole systems.
The fits are done using Gaussian Process Regression (GPR) and also provide an
error estimate along with the fit value.

Usage:

import surfinBH

# First load the model
fit_3dq8 = surfinBH.LoadFits('3dq8')

# Mass ratio and component spins along orbital angular momentum direction
q = 6.7
chi1z = 0.74
chi2z = -0.6
x = [q, chi1z, chi2z]

## Evaluate the fits and GPR error estimate.

# Final mass
mC, mC_err_est = fit_3dq8('mC', x)

# Final spin along orbital angular momentum direction
chiCz, chiCz_err_est = fit_3dq8('chiC', x)

# Final kick velocity (velCz is zero for nonprecessing systems)
velC, velC_err_est = fit_3dq8('velC', x)
velCx, velCy = velC
velCx_err_est, velCy_err_est = velC_err_est
    """

    #-------------------------------------------------------------------------
    def __init__(self, name, **kwargs):
        surfinBH.SurFinBH.__init__(self, name)
        self.fit_keys = ['mC', 'chiC', 'velC']


    #-------------------------------------------------------------------------
    def load_fits(self, h5file):
        """ Loads fits from h5file and returns a dictionary of fits. """
        fits = {}
        for key in ['mC', 'chiCz', 'velCx', 'velCy']:
            fits[key] = self.load_scalar_fit(fit_key=key, h5file=h5file)
        return fits


    #-------------------------------------------------------------------------
    def __call__(self, fit_key, x, **kwargs):
        """
Evaluates fits for the 3dq8 model.

Transforms the input parameter to fit parameters before evaluating the fit.
That is, maps from [q, chi1z, chi2z] to [np.log(q), chiHat, chi_a]
chiHat is defined in Eq.(3) of 1508.07253.
chi_a = (chi1z - chi2z)/2.
        """
        q, chi1z, chi2z = x

        eta = q/(1.+q)**2
        chi_wtAvg = (q*chi1z+chi2z)/(1.+q)
        chiHat = (chi_wtAvg - 38.*eta/113.*(chi1z + chi2z))/(1. - 76.*eta/113.)
        chi_a = (chi1z - chi2z)/2.

        mapped_x = [np.log(q), chiHat, chi_a]

        if fit_key == 'mC':
            mC, mC_err = self.evaluate_fits(mapped_x, 'mC')
            return mC, mC_err
        elif fit_key == 'chiC':
            chiCz, chiCz_err = self.evaluate_fits(mapped_x, 'chiCz')
            return chiCz, chiCz_err
        elif fit_key == 'velC':
            velCx, velCx_err = self.evaluate_fits(mapped_x, 'velCx')
            velCy, velCy_err = self.evaluate_fits(mapped_x, 'velCy')
            return np.array([velCx, velCy]), np.array([velCx_err, velCy_err])
        else:
            raise ValueError('Invalid fit_key')
