import numpy as np
from surfinBH import surfinBH

#=============================================================================
class Fit3dq8(surfinBH.SurFinBH):
    """
A class for the surfinBH3dq8 model presented in Varma et al., 2018, in prep.
This model predicts the final mass mC, final spin chiC and final kick velocity
velC, for the remnants of nonprecessing binary black hole systems.  The fits
are done using Gaussian Process Regression (GPR) and also provide an error
estimate along with the fit value.

This model has been trained in the parameter space:
    q <= 8, |chi1z| <= 0.8, |chi2z| <= 0.8

However, it extrapolates reasonably to:
    q <= 10, |chi1z| <= 1, |chi2z| <= 1

IMPORTANT NOTE: The kick vector is defined in the followin frame. The z-axis is
along the orbital angular momentum. The x-axis is along the line of separation
from the smaller BH to the larger BH at t=-100M, when t=0 occurs at the peak of
the waveform. The y-axis completes the triad.

Usage:

import surfinBH

fit_name = 'surfinBH3dq8'

# Get data for the fit. This only needs to done **once, ever**.
surfinBH.DownloadData(fit_name)

# Load the fit. This only needs to be done **once** at the start of your script.
fit = surfinBH.LoadFits(fit_name)

# Mass ratio and component spins along orbital angular momentum direction
q = 6.7
chi1z = 0.74
chi2z = -0.6
x = [q, chi1z, chi2z]

## Evaluate the fits and GPR error estimate.

# Final mass and its 1-sigma error etimate
mC, mC_err_est = fit('mC', x)

# Final spin and its 1-sigma error estimate
chiCz, chiCz_err_est = fit('chiC', x)

# Final kick vector and its 1-sigma error estimate
# NOTE: velCz is zero for nonprecessing systems
velC, velC_err_est = fit('velC', x)
velCx, velCy = velC
velCx_err_est, velCy_err_est = velC_err_est
    """

    #-------------------------------------------------------------------------
    def __init__(self, name):

        # soft_lims -> raise warning when outside lims
        # hard_lim -> raise error when outside lims
        # Same order as x in the call function. Each element is
        # a [minVal, maxVal] pair.
        soft_param_lims = [[0.99, 8.01], [-0.801, 0.801], [-0.801, 0.801]]
        hard_param_lims = [[0.99, 10.01], [-1, 1], [-1, 1]]
        super(Fit3dq8, self).__init__(name, soft_param_lims, hard_param_lims)


    #-------------------------------------------------------------------------
    def load_fits(self, h5file):
        """ Loads fits from h5file and returns a dictionary of fits. """
        fits = {}
        for key in ['mC', 'chiCz', 'velCx', 'velCy']:
            fits[key] = self.load_scalar_fit(fit_key=key, h5file=h5file)
        return fits

    #-------------------------------------------------------------------------
    def get_fit_params(self, x, fit_key):
        """
Transforms the input parameter to fit parameters for the 3dq8 model.
That is, maps from [q, chi1z, chi2z] to [np.log(q), chiHat, chi_a]
chiHat is defined in Eq.(3) of 1508.07253.
chi_a = (chi1z - chi2z)/2.
        """
        q, chi1z, chi2z = x
        eta = q/(1.+q)**2
        chi_wtAvg = (q*chi1z+chi2z)/(1.+q)
        chiHat = (chi_wtAvg - 38.*eta/113.*(chi1z + chi2z))/(1. - 76.*eta/113.)
        chi_a = (chi1z - chi2z)/2.
        fit_params = [np.log(q), chiHat, chi_a]
        return fit_params


    #-------------------------------------------------------------------------
    def __call__(self, fit_key, x, **kwargs):
        """
Evaluates fits for the 3dq8 model.

        """
        # Warn/Exit if extrapolating
        self.check_param_limits(x)

        if fit_key == 'mC':
            mC, mC_err = self.evaluate_fits(x, 'mC')
            return mC, mC_err
        elif fit_key == 'chiC':
            chiCz, chiCz_err = self.evaluate_fits(x, 'chiCz')
            return chiCz, chiCz_err
        elif fit_key == 'velC':
            velCx, velCx_err = self.evaluate_fits(x, 'velCx')
            velCy, velCy_err = self.evaluate_fits(x, 'velCy')
            return np.array([velCx, velCy]), np.array([velCx_err, velCy_err])
        else:
            raise ValueError('Invalid fit_key')
