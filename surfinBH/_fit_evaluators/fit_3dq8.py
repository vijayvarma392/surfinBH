import numpy as np
from surfinBH import surfinBH
import warnings

#=============================================================================
class Fit3dq8(surfinBH.SurFinBH):
    """ A class for the surfinBH3dq8 model presented in Varma et al., 2018,
    in prep. This model predicts the final mass mf, final spin chif and final
    kick velocity vf, for the remnants of nonprecessing binary black hole
    systems. The fits are done using Gaussian Process Regression (GPR) and
    also provide an error estimate along with the fit value.

    This model has been trained in the parameter space:
        q <= 8, |chiAz| <= 0.8, |chiBz| <= 0.8

    However, it extrapolates reasonably to:
        q <= 10, |chiAz| <= 1, |chiBz| <= 1

    =========================================================================
    Usage:

    import surfinBH

    # Load the fit
    fit = surfinBH.LoadFits('surfinBH3dq8')

    We provide the following call methods:
        # remnant mass and 1-sigma error estimate
        mf, mf_err = fit.mf(q, chiA, chiB, **kwargs)

        # remnant spin and 1-sigma error estimate
        chif, chif_err = fit.chif(q, chiA, chiB, **kwargs)

        # remnant recoil kick and 1-sigma error estimate
        vf, vf_err = fit.vf(q, chiA, chiB, **kwargs)

        # All of these together
        mf, chif, vf, mf_err, chif_err, vf_err
            = fit.all(q, chiA, chiB, **kwargs)

    The arguments for each of these call methods are as follows:
    Arguments:
        q:      Mass ratio (q>=1)

        chiA:   Dimensionless spin of the larger BH (array of size 3).

        chiB:   Dimensionless spin of the smaller BH (array of size 3).
                This model allows only nonprecessing spins, so only the
                z-components of these arrays should be non-zero.

    Optional arguments:
        allow_extrap:
            If False, raises a warning when q > 8.1 or |chiA|,|chiB| > 0.81,
                and raises an error when q > 10.1 or |chiA|,|chiB| > 1.
            If True, allows extrapolation to any q and |chiA|,|chiB| <= 1.
                Use at your own risk.
            Default: False.

    The spin and kick vectors are defined in the coorbital frame at t=-100 M
    from the peak of the waveform. This frame is defined as:
    The z-axis is along the orbital angular momentum direction of the binary.
    The x-axis is along the line of separation from the smaller BH to
        the larger BH at this time.
    The y-axis completes the triad.
    We obtain this frame from the waveform as defined in arxiv:1705.07089.
    """

    #-------------------------------------------------------------------------
    def __init__(self, name):

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
        super(Fit3dq8, self).__init__(name, soft_param_lims, hard_param_lims,
                aligned_spin_only=True)

    #-------------------------------------------------------------------------
    def _load_fits(self, h5file):
        """ Loads fits from h5file and returns a dictionary of fits. """
        fits = {}
        for key in ['mf', 'chifz', 'vfx', 'vfy']:
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
    def _eval_wrapper(self, fit_key, q, chiA, chiB, **kwargs):
        """ Evaluates the surfinBH3dq8 model.
        """
        chiA = np.array(chiA)
        chiB = np.array(chiB)

        # Warn/Exit if extrapolating
        allow_extrap = kwargs.pop('allow_extrap', False)
        self._check_param_limits(q, chiA, chiB, allow_extrap)

        self._check_unused_kwargs(kwargs)

        x = [q, chiA[2], chiB[2]]
        if fit_key == 'mf' or fit_key == 'all':
            mf, mf_err = self._evaluate_fits(x, 'mf')
            if fit_key == 'mf':
                return mf, mf_err
        if fit_key == 'chif' or fit_key == 'all':
            chifz, chifz_err = self._evaluate_fits(x, 'chifz')
            chif = np.array([0,0,chifz])
            chif_err = np.array([0,0,chifz_err])
            if fit_key == 'chif':
                return chif, chif_err
        if fit_key == 'vf' or fit_key == 'all':
            vfx, vfx_err = self._evaluate_fits(x, 'vfx')
            vfy, vfy_err = self._evaluate_fits(x, 'vfy')
            vf = np.array([vfx, vfy, 0])
            vf_err = np.array([vfx_err, vfy_err, 0])
            if fit_key == 'vf':
                return vf, vf_err
        if fit_key == 'all':
            return mf, chif, vf, mf_err, chif_err, vf_err
