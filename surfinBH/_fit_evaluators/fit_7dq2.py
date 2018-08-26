import numpy as np
from surfinBH import surfinBH
import warnings

#=============================================================================
class Fit7dq2(surfinBH.SurFinBH):
    """
A class for the surfinBH7dq2 model presented in Varma et al., 2018, in prep.
This model predicts the final mass mC, final spin chiC and final kick velocity
velC, for the remnants of precessing binary black hole systems.  The fits are
done using Gaussian Process Regression (GPR) and also provide an error estimate
along with the fit value.

This model has been trained in the parameter space:
    q <= 2, |chiA| <= 0.8, |chiB| <= 0.8

However, it extrapolates reasonably to:
    q <= 3, |chiA| <= 1, |chiB| <= 1

IMPORTANT NOTE: The component spins, remnant spin and kick vectors are defined
in the followin frame. The z-axis is along the orbital angular momentum at
t=-100M, when t=0 occurs at the peak of the waveform. The x-axis is along the
line of separation from the smaller BH to the larger BH at this time. The
y-axis completes the triad.

Usage:

import surfinBH

fit_name = 'surfinBH7dq2'

# Get data for the fit. This only needs to done **once, ever**.
surfinBH.DownloadData(fit_name)

# Load the fit. This only needs to be done **once** at the start of your script.
fit = surfinBH.LoadFits(fit_name)

# Mass ratio and component spins
q = 1.2
chiA = [0.1, 0.2, 0.3]
chiB = [0.2, -0.5, 0.3]
x = [q] + chiA + chiB

## Evaluate the fits and GPR error estimate.

# Final mass and its 1-sigma error etimate
mC, mC_err_est = fit('mC', x)

# Final spin vector and its 1-sigma error estimate
chiC, chiC_err_est = fit('chiC', x)

# Final kick vector and its 1-sigma error estimate
velC, velC_err_est = fit('velC', x)
    """

    #-------------------------------------------------------------------------
    def __init__(self, name):

        #NOTE: These are not the actual limits.
        # We override check_param_limits() to set limits, the
        # only purpose these serve is to generate regression data in
        # surfinBH/test/generate_regression_data.py
        soft_param_lims = None
        hard_param_lims = [[0.99, 3.01],
                [-0.6, 0.6],
                [-0.6, 0.6],
                [-0.6, 0.6],
                [-0.6, 0.6],
                [-0.6, 0.6],
                [-0.6, 0.6]]

        super(Fit7dq2, self).__init__(name, soft_param_lims, hard_param_lims)

    #-------------------------------------------------------------------------
    def load_fits(self, h5file):
        """ Loads fits from h5file and returns a dictionary of fits. """
        fits = {}
        for key in ['mC']:
            fits[key] = self.load_scalar_fit(fit_key=key, h5file=h5file)
        for key in ['chiC', 'velC']:
            fits[key] = self.load_vector_fit(key, h5file)
        return fits

    #-------------------------------------------------------------------------
    def get_fit_params(self, x, fit_key):
        """
Transforms the input parameter to fit parameters for the 7dq2 model.
That is, maps from
x = [q, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z]
fit_params = [np.log(q), chi1x, chi1y, chiHat, chi2x, chi2y, chi_a]

chiHat is defined in Eq.(3) of 1508.07253, but with chi1z and chi2z instead
of chi1 and chi2.
chi_a = (chi1z - chi2z)/2.
        """
        q, chi1z, chi2z = x[0], x[3], x[6]
        eta = q/(1.+q)**2
        chi_wtAvg = (q*chi1z+chi2z)/(1.+q)
        chiHat = (chi_wtAvg - 38.*eta/113.*(chi1z + chi2z))/(1. - 76.*eta/113.)
        chi_a = (chi1z - chi2z)/2.

        fit_params = x
        fit_params[0] = np.log(q)
        fit_params[3] = chiHat
        fit_params[6] = chi_a

        return fit_params

    #-------------------------------------------------------------------------
    def check_param_limits(self, x):
        """ Checks that x is within allowed range of paramters.
            Raises a warning if outside self.soft_param_lims and
            raises an error if outside self.hard_param_lims.
            If these are None, skips the checks.
        """
        q = x[0]
        chiAmag = np.sqrt(np.sum(x[1:4]**2))
        chiBmag = np.sqrt(np.sum(x[4:7]**2))

        if q < 1:
            raise ValueError('Mass ratio should be >= 1.')
        elif q > 3.01:
            raise Exception('Mass ratio outside allowed range.')
        elif q > 2.01:
            warnings.warn('Mass ratio outside training range.')

        if chiAmag > 1.:
            raise Exception('Spin magnitude of BhA outside allowed range.')
        elif chiAmag > 0.81:
            warnings.warn('Spin magnitude of BhA outside training range.')

        if chiBmag > 1.:
            raise Exception('Spin magnitude of BhB outside allowed range.')
        elif chiBmag > 0.81:
            warnings.warn('Spin magnitude of BhB outside training range.')



    #-------------------------------------------------------------------------
    def __call__(self, fit_key, x, **kwargs):

        x = np.array(x)

        # Warn/Exit if extrapolating
        self.check_param_limits(x)

        if fit_key == 'mC':
            mC, mC_err = self.evaluate_fits(x, fit_key)
            return mC, mC_err
        elif fit_key == 'chiC' or fit_key == 'velC':
            res = self.evaluate_fits(x, fit_key)
            return res.T[0], res.T[1]
        else:
            raise ValueError('Invalid fit_key')
