"""surfinBH
========

Surrogate final black hole properties for mergers of binary black holes.
See https://pypi.org/project/surfinBH/ for more details.
"""
__copyright__ = "Copyright (C) 2018 Vijay Varma"
__email__ = "vvarma@caltech.edu"
__status__ = "testing"
__author__ = "Vijay Varma"
__version__ = "0.0.8.dev5"
__license__ = """
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import numpy as np
import os, sys
import h5py
import warnings

import _eval_pysur
from _dataPath import DataPath


#=============================================================================
class SurFinBH(object):
    """
Class to load and evaluate surrogate fits for final BH properties.
Each derived class should do the following:

    1. define _load_fits(self, h5file)
    2. define _get_fit_params(self, x, fit_key)
    3. define _eval_wrapper(self, fit_key, x, **kwargs)
    4. define _check_param_limits(self, q, chiA, chiB)

See _fit_evaluators.fit_7dq2.py for an example.
    """

    #-------------------------------------------------------------------------
    def __init__(self, name):
        """
        name:           Name of the fit excluding the surfinBH prefix. Ex: 7dq2.
        """
        self.name = name
        h5file = h5py.File('%s/fit_%s.h5'%(DataPath(), name), 'r')
        self.fits = self._load_fits(h5file)
        h5file.close()

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
    def _load_scalar_fit(self, fit_key=None, h5file=None, fit_data=None):
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
    def _load_vector_fit(self, fit_key, h5file):
        """ Loads a vector of fits
        """
        vector_fit = []
        for i in range(len(h5file[fit_key].keys())):
            fit_data = self._read_dict(h5file[fit_key]['comp_%d'%i])
            vector_fit.append(self._load_scalar_fit(fit_data=fit_data))
        return vector_fit

    #-------------------------------------------------------------------------
    def _evaluate_fits(self, x, fit_key):
        """ Evaluates a particular fit by passing fit_key to self.fits.
            Assumes self._get_fit_params() has been overriden.
        """
        fit = self.fits[fit_key]
        fit_params = self._get_fit_params(np.copy(x), fit_key)
        if type(fit) == list:
            res = []
            for i in range(len(fit)):
                res.append(fit[i](fit_params))
            return np.array(res)
        else:
            return fit(fit_params)

    #-------------------------------------------------------------------------
    def _check_unused_kwargs(self, kwargs):
        """ Call this at the end of call module to check if all the kwargs have
        been used. Assumes kwargs were extracted using pop.
        """
        if len(kwargs.keys()) != 0:
            unused = ""
            for k in kwargs.keys():
                unused += "'%s', "%k
            if unused[-2:] == ", ":     # get rid of trailing comma
                unused = unused[:-2]
            raise Exception('Unused keys in kwargs: %s'%unused)

    #-------------------------------------------------------------------------
    #----------------------  Override these  ---------------------------------
    #-------------------------------------------------------------------------

    #-------------------------------------------------------------------------
    def _load_fits(self, h5file):
        """ Loads fits from h5file and returns a dictionary of fits. """
        raise NotImplementedError("Please override me.")
        return fits

    #-------------------------------------------------------------------------
    def _get_fit_params(self, x, fit_key):
        """ Maps from input params x to the fit_params used to evaluate the
            fit.
        """
        raise NotImplementedError("Please override me.")
        return fit_params

    def _check_param_limits(self, q, chiA, chiB, **kwargs):
        """ Checks that the params are withing allowed range.
        """
        raise NotImplementedError("Please override me.")

    #-------------------------------------------------------------------------
    def _eval_wrapper(self, fit_key, q, chiA, chiB, **kwargs):
        """ Evaluates a particular fit. This varies for each surrogate.
            Allowed values for fit_key are 'mC', 'chiC' and 'velC' and 'all'.
            chiA and chiB should have size 3.

            Each derived class should have its own _eval_wrapper function but
            call self._check_param_limits() first to do some sanity checks.
            See _fit_evaluators.fit_7dq2.py for an example.
        """
        raise NotImplementedError("Please override me.")


    #-------------------------------------------------------------------------
    #----------------------   Call methods   ---------------------------------
    #-------------------------------------------------------------------------

    def mC(self, *args, **kwargs):
        """ Evaluates fit and 1-sigma error estimate for remnant mass.
        Returns:
            mC, mC_err_est
        """
        return self._eval_wrapper('mC', *args, **kwargs)

    def chiC(self, *args, **kwargs):
        """ Evaluates fit and 1-sigma error estimate for remnant spin.
        Returns:
            chiC, chiC_err_est

        chiC and chiC_err_est are arrays of size 3.
        """
        return self._eval_wrapper('chiC', *args, **kwargs)

    def velC(self, *args, **kwargs):
        """ Evaluates fit and 1-sigma error estimate for remnant kick velocity.
        Returns:
            velC, velC_err_est

        velC and velC_err_est are arrays of size 3.
        """
        return self._eval_wrapper('velC', *args, **kwargs)

    def all(self, *args, **kwargs):
        """ Evaluates fit and 1-sigma error estimate for remnant mass, spin
        and kick velocity.
        Returns:
            mC, chiC, velC, mC_err_est, chiC_err_est, velC_err_est

        chiC, velC, chiC_err_est and velC_err_est are arrays of size 3.
        """
        return self._eval_wrapper('all', *args, **kwargs)

