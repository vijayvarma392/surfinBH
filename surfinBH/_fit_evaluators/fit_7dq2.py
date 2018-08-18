import numpy as np
from surfinBH import surfinBH

#=============================================================================
class Fit7dq2(surfinBH.SurFinBH):
    """ Evaluates Surrogate fits for Final BH properties. """

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
