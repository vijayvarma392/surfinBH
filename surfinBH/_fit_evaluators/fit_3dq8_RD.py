import numpy as np
from surfinBH import surfinBH
import warnings

#=============================================================================
class Fit3dq8_RD(surfinBH.SurFinBH):
    """A class for the NRSur3dq8_RD model presented in Magana Zertuche et al.,
    arxiv:TODO.

    This model predicts the final mass mf, final spin chif, and
    complex QNM amplitudes A_(l,m,n,p), for the remnants of
    nonprecessing binary black hole systems. The fits are done using
    Gaussian Process Regression (GPR) and also provide an error
    estimate along with the fit value.

    This model has been trained in the parameter space:
        q <= 8, |chiAz| <= 0.8, |chiBz| <= 0.8

    However, it extrapolates reasonably to:
        q <= 10, |chiAz| <= 1, |chiBz| <= 1

    =========================================================================
    Usage:

    import surfinBH

    # Load the fit
    fit = surfinBH.LoadFits('NRSur3dq8_RD')

    We provide the following call methods:
        # remnant mass and 1-sigma error estimate
        mf, mf_err = fit.mf(q, chiA, chiB, **kwargs)

        # remnant spin and 1-sigma error estimate
        chif, chif_err = fit.chif(q, chiA, chiB, **kwargs)

        # remnant mass, spin, QNM amplitudes, and error estimates
        mf, chif, QNM_dict, mf_err, chif_err, QNM_err_est_dict
            = fit.all(q, chiA, chiB, modes=mode_list,
                      **kwargs)

    The arguments for each of these call methods are as follows:
    Arguments:
        q:      Mass ratio (q>=1)

        chiA:   Dimensionless spin of the larger BH (array of size 3).

        chiB:   Dimensionless spin of the smaller BH (array of size 3).
                This model allows only nonprecessing spins, so only the
                z-components of these arrays should be non-zero.

    Optional arguments:
        modes:
            A list of mode labels of the form (l,m,n,p) which are a
            subset of the modeled modes. The default is all the
            modeled modes, namely
            modes = [ (2,2,0,1),(2,-2,0,-1),(2,2,1,1),(2,-2,1,-1),
                      (2,0,0,1),(2,0,0,-1),(4,4,0,1),(4,-4,0,-1),
                      (3,2,0,1),(3,-2,0,-1) ]
        allow_extrap:
            If False, raises a warning when q > 8.1 or |chiA|,|chiB| > 0.81,
                and raises an error when q > 10.1 or |chiA|,|chiB| > 1.
            If True, allows extrapolation to any q and |chiA|,|chiB| <= 1.
                Use at your own risk.
            Default: False.

    The remnant spin is only the magnitude. The complex QNM amplitudes
    are given in the superrest frame of the remnant at 20M after the
    peak.

    """

    #-------------------------------------------------------------------------
    def __init__(self, name):

        # Set of modes available in this model
        self.modeled_modes = [ (2,2,0,1),(2,-2,0,-1),(2,2,1,1),(2,-2,1,-1),
                               (2,0,0,1),(2,0,0,-1),(4,4,0,1),(4,-4,0,-1),
                               (3,2,0,1),(3,-2,0,-1) ]

        # Param limits beyond which to raise a warning
        soft_param_lims = {
            'q': 8.1,
            'chiAmag': 0.81,
            'chiBmag': 0.81,
                }

        # Param limits beyond which to raise an error
        hard_param_lims = {
            'q': 10.1,
            'chiAmag': 1,
            'chiBmag': 1,
                }
        super(Fit3dq8_RD, self).__init__(name, soft_param_lims, hard_param_lims,
                aligned_spin_only=True)

    #-------------------------------------------------------------------------
    def _mode2keys(self, mode):
        """Utility function to convert (l,m,n,p)
        into the pair e.g. ("A_l,m,n,p_r", "A_l,m,n,p_i")"""
        mode_str = ','.join([str(i) for i in mode])
        return ("A_" + mode_str + "_r", "A_" + mode_str + "_i")

    #-------------------------------------------------------------------------
    def _load_fits(self, h5file):
        """ Loads fits from h5file and returns a dictionary of fits. """
        fits = {}

        mode_keys = [s for mode in self.modeled_modes
                       for s in self._mode2keys(mode)]
        keys = ["chi_f", "M_f"] + mode_keys

        for key in keys:
            fits[key] = self._load_scalar_fit(fit_key=key, h5file=h5file)
        return fits

    #-------------------------------------------------------------------------
    def _get_fit_params(self, x, fit_key):
        """ Transforms the input parameter to fit parameters for the 3dq8 model.
    That is, maps from [q, chiAz, chiBz] to [np.log(q), chiHat, chi_a]
    chiHat is defined in Eq.(3) of 1508.07253.
    chi_a = (chiAz - chiBz)/2.
        """
        q, chiAz, chiBz = x
        eta = q/(1.+q)**2
        chi_wtAvg = (q*chiAz+chiBz)/(1.+q)
        chiHat = (chi_wtAvg - 38.*eta/113.*(chiAz + chiBz))/(1. - 76.*eta/113.)
        chi_a = (chiAz - chiBz)/2.
        fit_params = [np.log(q), chiHat, chi_a]
        return fit_params




    #-------------------------------------------------------------------------
    def vf(self, *args, **kwargs):
        """vf is not implemented in this model."""
        raise NotImplementedError("vf is not implemented in this model.")

    def all(self, *args, **kwargs):
        """ Evaluates fit and 1-sigma error estimate for remnant mass, spin
        and QNM amplitudes.
        Returns:
            mf, chif, QNM_dict, mf_err_est, chif_err_est, QNM_err_est_dict

        chif and chif_err_est are only the spin magnitude.
        """
        return self._eval_wrapper('all', *args, **kwargs)


    #-------------------------------------------------------------------------
    def _eval_wrapper(self, fit_key, q, chiA, chiB, **kwargs):
        """ Evaluates the NRSur3dq8_RD model.
        """
        chiA = np.array(chiA)
        chiB = np.array(chiB)

        # Warn/Exit if extrapolating
        allow_extrap = kwargs.pop('allow_extrap', False)
        self._check_param_limits(q, chiA, chiB, allow_extrap)

        modes = kwargs.pop('modes', self.modeled_modes)

        self._check_unused_kwargs(kwargs)

        x = [q, chiA[2], chiB[2]]
        if fit_key == 'mf' or fit_key == 'all':
            mf, mf_err = self._evaluate_fits(x, 'M_f')
            if fit_key == 'mf':
                return mf, mf_err
        if fit_key == 'chif' or fit_key == 'all':
            chif, chif_err = self._evaluate_fits(x, 'chi_f')
            if fit_key == 'chif':
                return chif, chif_err

        QNM_dict = {}
        QNM_err_est_dict = {}

        for mode in modes:
            mode_str_r, mode_str_i = self._mode2keys(mode)
            A_r, A_r_err = self._evaluate_fits(x, mode_str_r)
            A_i, A_i_err = self._evaluate_fits(x, mode_str_i)
            A     = A_r + 1.j * A_i
            A_err = A_r_err + 1.j * A_i_err
            QNM_dict[mode] = A
            QNM_err_est_dict[mode] = A_err

        if fit_key == 'all':
            return mf, chif, QNM_dict, mf_err, chif_err, QNM_err_est_dict
