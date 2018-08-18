import numpy as np
from surfinBH import surfinBH

#=============================================================================
class Fit3dq8(surfinBH.SurFinBH):
    """ Evaluates Surrogate fits for Final BH properties. """

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
