import numpy as np
from surfinBH import surfinBH

#=============================================================================
class Fit3dq8(surfinBH.SurFinBH):
    """ A class for the surfinBH3dq8 model presented in Varma et al., 2018,
    in prep. This model predicts the final mass mC, final spin chiC and final
    kick velocity velC, for the remnants of nonprecessing binary black hole
    systems. The fits are done using Gaussian Process Regression (GPR) and
    also provide an error estimate along with the fit value.

    This model has been trained in the parameter space:
        q <= 8, |chiAz| <= 0.8, |chiBz| <= 0.8

    However, it extrapolates reasonably to:
        q <= 10, |chiAz| <= 1, |chiBz| <= 1

    See __call__ method for evaluation.
    """

    #-------------------------------------------------------------------------
    def __init__(self, name):

        # soft_lims -> raise warning when outside lims
        # hard_lim -> raise error when outside lims
        # Same order as x in the call function. Each element is
        # a [minVal, maxVal] pair.
        soft_param_lims = [[0.99, 8.01], [-0.801, 0.801], [-0.801, 0.801]]
        hard_param_lims = [[0.99, 10.01], [-1, 1], [-1, 1]]
        super(Fit3dq8, self).__init__(name, soft_param_lims, hard_param_lims)


    #-------------------------------------------------------------------------
    def _load_fits(self, h5file):
        """ Loads fits from h5file and returns a dictionary of fits. """
        fits = {}
        for key in ['mC', 'chiCz', 'velCx', 'velCy']:
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
    def __call__(self, fit_key, x, **kwargs):
        """ Evaluates the surfinBH3dq8 model.

    Arguments:
        fit_key:
            'mC', 'chiC' or 'velC', for remnant mass, spin and kick velocity,
            respectively.

        x:
            Array of [q, chiAz, chiBz], where q >= 1 is the mass ratio and
            chiAz (chiBz) is the dimensionless spin of the larger (smaller) BH
            along the z-direction.

    Returns:
        If fit_key='mC':
            returns mC, mC_err_est
            The value and 1-sigma error estimate in remnant mass.
        If fit_key='chiC':
            returns chiCz, chiCz_err_est:
            The value and 1-sigma error estimate in remnant spin along the
            z-direction.
        If fit_key='velC':
            returns [velCx, velCy], [velCx_err_est, velCy_err_est]:
            The value and 1-sigma error estimate in the in-plane remnant kick
            vector.

    The spins and kick are defined in the coorbital frame at t=-100 M from the
    peak of the waveform. This frame is defined as:
    The z-axis is along the orbital angular momentum direction of the binary.
    The x-axis is along the line of separation from the smaller BH to
        the larger BH at this time.
    The y-axis completes the triad.
    We obtain this frame from the waveform as defined in arxiv:1705.07089.
        """
        # Warn/Exit if extrapolating
        self._check_param_limits(x)

        if fit_key == 'mC':
            mC, mC_err = self._evaluate_fits(x, 'mC')
            return mC, mC_err
        elif fit_key == 'chiC':
            chiCz, chiCz_err = self._evaluate_fits(x, 'chiCz')
            return chiCz, chiCz_err
        elif fit_key == 'velC':
            velCx, velCx_err = self._evaluate_fits(x, 'velCx')
            velCy, velCy_err = self._evaluate_fits(x, 'velCy')
            return np.array([velCx, velCy]), np.array([velCx_err, velCy_err])
        else:
            raise ValueError('Invalid fit_key')
