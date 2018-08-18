import numpy as np
import os, sys
import h5py

import _eval_pysur

#=============================================================================
class SurFinBH(object):
    """ Evaluates Surrogate fits for Final BH properties. """

    #-------------------------------------------------------------------------
    def __init__(self, name):
        self.name = name
        h5file = h5py.File('%s/fit_%s.h5'%(self.data_path(), name), 'r')
        self.fits = self.load_fits(h5file)
        h5file.close()

    #-------------------------------------------------------------------------
    def data_path(self):
        """ Return the default path for fit data h5 files"""
        return os.path.abspath('%s/../data'%(os.path.dirname( \
            os.path.realpath(__file__))))

    #-------------------------------------------------------------------------
    def _read_dict(self, f):
        """ Converts h5 groups to dictionaries
        """
        d = {}
        for k, item in f.iteritems():
            if type(item) == h5py._hl.dataset.Dataset:
                v = item.value
                if type(v) == np.string_:
                    v = str(v)
                if type(v) == str and v == "NONE":
                    d[k] = None
                elif type(v) == str and v == "EMPTYARR":
                    d[k] = np.array([])
                else:
                    d[k] = v
            elif k[:5] == "DICT_":
                d[k[5:]] = self._read_dict(item)
            elif k[:5] == "LIST_":
                tmpD = self._read_dict(item)
                d[k[5:]] = [tmpD[str(i)] for i in range(len(tmpD))]
        return d

    #-------------------------------------------------------------------------
    def load_scalar_fit(self, fit_key=None, h5file=None, fit_data=None):
        """ Loads a single fit
        """
        if (fit_key is None) ^ (h5file is None):
            raise ValueError("Either specify both fit_key and h5file, or"
                " neither")

        if not ((fit_key is None) ^ (fit_data is None)):
            raise ValueError("Specify exactly one of fit_key and fit_data.")

        if fit_data is None:
            fit_data = self._read_dict(h5file[fit_key])

        if 'fitType' in fit_data.keys() and fit_data['fitType'] == 'GPR':
            fit = _eval_pysur.evaluate_fit.getGPRFitAndErrorEvaluator(fit_data)
        else:
            fit = _eval_pysur.evaluate_fit.getFitEvaluator(fit_data)

        return fit

    #-------------------------------------------------------------------------
    def load_vector_fit(self, fit_key, h5file):
        """ Loads a vector of fits
        """
        vector_fit = []
        for i in range(len(h5file[fit_key].keys())):
            fit_data = self._read_dict(h5file[fit_key]['comp_%d'%i])
            vector_fit.append(self.load_scalar_fit(fit_data=fit_data))
        return vector_fit

    #-------------------------------------------------------------------------
    def evaluate_fits(self, x, fit_key):
        """ Evaluates a particular fit by passing fit_key to self.fits """
        fit = self.fits[fit_key]
        if type(fit) == list:
            res = []
            for i in range(len(fit)):
                res.append(fit[i](x))
            return np.array(res)
        else:
            return fit(x)

    #-------------------------------------------------------------------------
    def load_fits(self, h5file):
        """ Loads fits from h5file and returns a dictionary of fits. """
        raise NotImplementedError("Please override me.")
        return fits

    #-------------------------------------------------------------------------
    def __call__(self, fit_key, x, **kwargs):
        """ Evaluates a particular fit. This varies for each surrogate.
            Allowed values for fit_key are 'mC', 'chiC' and 'velC'.
            If fit_key is a scalar, returns the fit value and error
            estimate (if available).
            If fit_key is a scalar, returns two arrays. The first array
            is a vector of fit values and the second is the corresponding
            vector of errors estimates (if available).

            See fit_evaluators.fit_7dq2.py for an example.
        """
        raise NotImplementedError("Please override me.")
