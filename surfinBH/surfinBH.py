"""surfinBH
========

Surrogate final black hole properties for mergers of binary black holes.
See https://pypi.org/project/surfinBH/ for more details.
"""
__copyright__ = "Copyright (C) 2018 Vijay Varma"
__email__ = "vvarma@caltech.edu"
__status__ = "testing"
__author__ = "Vijay Varma"
__version__ = "0.1.3.dev0"
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

from . import _eval_pysur
from ._dataPath import DataPath


#=============================================================================
class SurFinBH(object):
    """
Class to load and evaluate surrogate fits for final BH properties.
Each derived class should do the following:

    1. define _load_fits(self, h5file)
    2. define _get_fit_params(self, x, fit_key)
    3. define _eval_wrapper(self, fit_key, x, **kwargs)
    4. define soft_param_lims and hard_param_lims.
    5. define _extra_regression_kwargs, to test any additional kwargs used in
          the _eval_wrapper method.

See _fit_evaluators.fit_7dq2.py for an example.
    """

    #-------------------------------------------------------------------------
    def __init__(self, name, soft_param_lims, hard_param_lims,
            aligned_spin_only=False):
        """
        name: Name of the fit excluding the surfinBH prefix. Ex: 7dq2.
        soft_param_lims: param limits beyond which to raise a warning.
        hard_param_lims: param limits beyond which to raise an error.
        aligned_spin_only: raise an error if given precessing spins.
        See _fit_evaluators.fit_7dq2.py for an example.
        """
        self.name = name
        self.soft_param_lims = soft_param_lims
        self.hard_param_lims = hard_param_lims
        self.aligned_spin_only = aligned_spin_only

        h5file = h5py.File('%s/fit_%s.h5'%(DataPath(), name), 'r')
        self.fits = self._load_fits(h5file)
        h5file.close()

    #-------------------------------------------------------------------------
    def _read_dict(self, f):
        """ Converts h5 groups to dictionaries
        """
        d = {}
        for k, item in f.items():
            if type(item) == h5py._hl.dataset.Dataset:
                v = item.value
                if type(v) == np.string_:
                    v = str(v)
                if type(v) == str and v == "NONE":
                    d[k] = None
                elif type(v) == str and v == "EMPTYARR":
                    d[k] = np.array([])
                elif isinstance(v, bytes):
                    d[k] = v.decode('utf-8')
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
    def _check_param_limits(self, q, chiA, chiB, allow_extrap):
        """ Checks that params are within allowed range of paramters.
        Raises a warning if outside self.soft_param_lims limits and
        raises an error if outside self.hard_param_lims.
        If allow_extrap=True, skips these checks.
        """

        if q < 1:
            raise ValueError('Mass ratio should be >= 1.')

        chiAmag = np.sqrt(np.sum(chiA**2))
        chiBmag = np.sqrt(np.sum(chiB**2))
        if chiAmag > 1:
            raise ValueError('Spin magnitude of BhA > 1.')
        if chiBmag > 1:
            raise ValueError('Spin magnitude of BhB > 1.')

        if self.aligned_spin_only:
            if np.sqrt(np.sum(chiA[:2]**2)) > 1e-10:
                raise ValueError('The x & y components of chiA should be zero.')
            if np.sqrt(np.sum(chiB[:2]**2)) > 1e-10:
                raise ValueError('The x & y components of chiB should be zero.')

        # Do not check param limits if allow_extrap=True
        if allow_extrap:
            return

        if q > self.hard_param_lims['q']:
            raise ValueError('Mass ratio outside allowed range.')
        elif q > self.soft_param_lims['q']:
            warnings.warn('Mass ratio outside training range.')

        if chiAmag > self.hard_param_lims['chiAmag']:
            raise ValueError('Spin magnitude of BhA outside allowed range.')
        elif chiAmag > self.soft_param_lims['chiAmag']:
            warnings.warn('Spin magnitude of BhA outside training range.')

        if chiBmag > self.hard_param_lims['chiBmag']:
            raise ValueError('Spin magnitude of BhB outside allowed range.')
        elif chiBmag > self.soft_param_lims['chiBmag']:
            warnings.warn('Spin magnitude of BhB outside training range.')

    #-------------------------------------------------------------------------
    def _generate_random_params_for_tests(self):
        """ Generate random parameters to use in tests.
        """
        # Generate params randomly within allowed values
        q = np.random.uniform(1, self.hard_param_lims['q'])
        chiAmag= np.random.uniform(0, self.hard_param_lims['chiAmag'])
        chiBmag= np.random.uniform(0, self.hard_param_lims['chiBmag'])
        if self.aligned_spin_only:
            chiAph, chiBph = 0, 0
            chiAth, chiBth = np.random.choice([0, np.pi]), \
                np.random.choice([0, np.pi])
        else:
            chiAth = np.random.uniform(0, np.pi)
            chiBth = np.random.uniform(0, np.pi)
            chiAph = np.random.uniform(0, 2*np.pi)
            chiBph = np.random.uniform(0, 2*np.pi)

        chiA = [chiAmag*np.sin(chiAth)*np.cos(chiAph),
                chiAmag*np.sin(chiAth)*np.sin(chiAph),
                chiAmag*np.cos(chiAth)]

        chiB = [chiBmag*np.sin(chiBth)*np.cos(chiBph),
                chiBmag*np.sin(chiBth)*np.sin(chiBph),
                chiBmag*np.cos(chiBth)]

        return q, chiA, chiB


    #-------------------------------------------------------------------------
    #----------------------  Override these  ---------------------------------
    #-------------------------------------------------------------------------

    #-------------------------------------------------------------------------
    def _load_fits(self, h5file):
        """ Loads fits from h5file and returns a dictionary of fits. """
        raise NotImplementedError("Please override me.")
        return fits


    #-------------------------------------------------------------------------
    def _extra_regression_kwargs(self):
        """ Add additional kwargs for regression tests. If not overriden,
            this will be empty. See _fit_evaluators.fit_7dq2.py for an example.
        """
        return []

    #-------------------------------------------------------------------------
    def _get_fit_params(self, x, fit_key):
        """ Maps from input params x to the fit_params used to evaluate the
            fit.
        """
        raise NotImplementedError("Please override me.")
        return fit_params

    #-------------------------------------------------------------------------
    def _eval_wrapper(self, fit_key, q, chiA, chiB, **kwargs):
        """ Evaluates a particular fit. This varies for each surrogate.
            Allowed values for fit_key are 'mf', 'chif' and 'vf' and 'all'.
            chiA and chiB should have size 3.

            Each derived class should have its own _eval_wrapper function but
            call self._check_param_limits() first to do some sanity checks.
            See _fit_evaluators.fit_7dq2.py for an example.
        """
        raise NotImplementedError("Please override me.")


    #-------------------------------------------------------------------------
    #----------------------   Call methods   ---------------------------------
    #-------------------------------------------------------------------------

    def mf(self, *args, **kwargs):
        """ Evaluates fit and 1-sigma error estimate for remnant mass.
        Returns:
            mf, mf_err_est
        """
        return self._eval_wrapper('mf', *args, **kwargs)

    def chif(self, *args, **kwargs):
        """ Evaluates fit and 1-sigma error estimate for remnant spin.
        Returns:
            chif, chif_err_est

        chif and chif_err_est are arrays of size 3.
        """
        return self._eval_wrapper('chif', *args, **kwargs)

    def vf(self, *args, **kwargs):
        """ Evaluates fit and 1-sigma error estimate for remnant kick velocity.
        Returns:
            vf, vf_err_est

        vf and vf_err_est are arrays of size 3.
        """
        return self._eval_wrapper('vf', *args, **kwargs)

    def all(self, *args, **kwargs):
        """ Evaluates fit and 1-sigma error estimate for remnant mass, spin
        and kick velocity.
        Returns:
            mf, chif, vf, mf_err_est, chif_err_est, vf_err_est

        chif, vf, chif_err_est and vf_err_est are arrays of size 3.
        """
        return self._eval_wrapper('all', *args, **kwargs)

