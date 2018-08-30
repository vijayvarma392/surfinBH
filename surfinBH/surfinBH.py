__copyright__ = "Copyright (C) 2018 Vijay Varma"
__email__ = "vvarma@caltech.edu"
__status__ = "testing"
__author__ = "Vijay Varma"
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
    3. define __call__(self, fit_key, x, **kwargs)
    4. define soft_param_lims and hard_param_lims, the limits for parameters
       beyond which warnings/errors are raised.

See _fit_evaluators.fit_7dq2.py for an example.
    """

    #-------------------------------------------------------------------------
    def __init__(self, name, soft_param_lims, hard_param_lims):
        """
        name:           Name of the fit excluding the surfinBH prefix. Ex: 7dq2.
        soft_param_lims: Parameter bounds beyond which a warning is raised.
        hard_param_lims: Parameter bounds beyond which an error is raised.
                         Same order and len as x in the __call__ function.
                         Each element is a [minVal, maxVal] pair.
                         Setting soft_param_lims/hard_param_lims to None will
                         skip that particular check.
        """
        self.name = name
        h5file = h5py.File('%s/fit_%s.h5'%(DataPath(), name), 'r')
        self.fits = self._load_fits(h5file)
        h5file.close()

        self.soft_param_lims = soft_param_lims
        self.hard_param_lims = hard_param_lims

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
        fit_params = self._get_fit_params(x, fit_key)
        if type(fit) == list:
            res = []
            for i in range(len(fit)):
                res.append(fit[i](fit_params))
            return np.array(res)
        else:
            return fit(fit_params)

    #-------------------------------------------------------------------------
    def _check_param_limits(self, x):
        """ Checks that x is within allowed range of paramters.
            Raises a warning if outside self.soft_param_lims and
            raises an error if outside self.hard_param_lims.
            If these are None, skips the checks.
        """
        if self.soft_param_lims is not None:
            if len(x) != len(self.soft_param_lims):
                raise Exception("Expected x to be of len=%d"\
                    %len(self.soft_param_lims))

        if self.hard_param_lims is not None:
            for i in range(len(x)):
                if x[i] < self.hard_param_lims[i][0] \
                        or x[i] > self.hard_param_lims[i][1]:

                    raise Exception('Parameter x[%d] outside allowed range.'%i)

        if self.soft_param_lims is not None:
            for i in range(len(x)):
                if x[i] < self.soft_param_lims[i][0] \
                        or x[i] > self.soft_param_lims[i][1]:

                    warnings.warn('Parameter x[%d] outside training range.'%i)

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

    #-------------------------------------------------------------------------
    def __call__(self, fit_key, x, **kwargs):
        """ Evaluates a particular fit. This varies for each surrogate.
            Allowed values for fit_key are 'mC', 'chiC' and 'velC'.
            If fit_key is a scalar, returns the fit value and error
            estimate (if available).
            If fit_key is a scalar, returns two arrays. The first array
            is a vector of fit values and the second is the corresponding
            vector of errors estimates (if available).

            Each derived class should have its own __call__ function but
            call self._check_param_limits() first to do some sanity checks.
            See _fit_evaluators.fit_7dq2.py for an example.
        """
        raise NotImplementedError("Please override me.")
