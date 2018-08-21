import numpy as np
from surfinBH import surfinBH

#=============================================================================
class Fit7dq2(surfinBH.SurFinBH):
    """
A class for the surfinBH7dq2 model presented in Varma et al., 2018, in prep.
This model predicts the final mass mC, final spin chiC and final kick velocity
velC, for the remnants of precessing binary black hole systems.  The fits are
done using Gaussian Process Regression (GPR) and also provide an error estimate
along with the fit value.

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
    def __init__(self, name, **kwargs):
        surfinBH.SurFinBH.__init__(self, name)
        self.fit_keys = ['mC', 'chiC', 'velC']

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
    def __call__(self, fit_key, x, **kwargs):
        if fit_key == 'mC':
            mC, mC_err = self.evaluate_fits(x, fit_key)
            return mC, mC_err
        elif fit_key == 'chiC' or fit_key == 'velC':
            res = self.evaluate_fits(x, fit_key)
            return res.T[0], res.T[1]
        else:
            raise ValueError('Invalid fit_key')
