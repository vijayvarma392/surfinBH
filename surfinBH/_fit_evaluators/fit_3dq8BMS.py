import numpy as np
from surfinBH import surfinBH
import warnings


#=============================================================================
class Fit3dq8BMS(surfinBH.SurFinBH):
    """ A class for the NRSur3dq8BMSRemnant model presented in Da Re et al.,
    arxiv:????.?????.
    
    This model predicts the supertransation modes up to ell = 8 and the 
    3 components of the boost velocity of the BMS transformation from the 
    inspiral (PN) BMS frame to the remnant black hole BMS frame. The boost 
    velocity coincides with the remnant black hole kick velocity as observed 
    from the inspiral (PN) BMS frame.
    
    The model was trained on nonprecessing binary black hole systems. The fits 
    are done using Gaussian Process Regression (GPR) and also provide an error 
    estimate along with the fit value.

    This model has been trained in the parameter space:
        q <= 8, |chiAz| <= 0.8, |chiBz| <= 0.8

    =========================================================================
    Usage:

    import surfinBH

    # Load the fit
    fit = surfinBH.LoadFits('NRSur3dq8BMSRemnant')

    #To evaluate the supertranslation parameter "alpha" and the boost velocity 
     "boost" together with their respective 1-sigma error estimates call:
    alpha, boost,  alpha_err, boost_err = fit.all(q, chiA, chiB, **kwargs)

    # NOTE: The ell=0 mode, corresponding to a time translation, 
      is set to be identically zero.

    # alpha is expressed as a complex array of spherical harmonics modes in the 
      order 
      (0,0),(1,-1),(1,0),(1,1),(2,-2),(2,-1),(2,0),...
      (same convention as the spherical_functions and scri packages)

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

    The inspiral frame can be defined up to an arbitrary U(1) rotation about 
    the z-axis. We fix this gauge freedom using the same frame alignment choice
    used in NrSur3dq8Remnant. The inspiral frame at -100M from the peak of the 
    waveform is defined as:
    The z-axis is along the orbital angular momentum direction of the binary.
    The x-axis is along the line of separation from the smaller BH to
        the larger BH at this time.
    The y-axis completes the triad.
    The rotation leading to this frame choice is included in the BMS transformation
    as outlined in arxiv:????.?????.
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
        super(Fit3dq8BMS, self).__init__(name, soft_param_lims, hard_param_lims,
                aligned_spin_only=True)

    #-------------------------------------------------------------------------
    def _load_fits(self, h5file):
        """ Loads fits from h5file and returns a dictionary of fits. """
        fits = {}
        ell_max = 8
        keys_sup_real = [f'S_{ell}{m}_real' for ell in range(1,ell_max + 1)
                                       for m in range(0,ell + 1)]
        keys_sup_imag = [f'S_{ell}{m}_imag' for ell in range(1,ell_max + 1)
                                       for m in range(1,ell + 1)]
        keys_vel = ['boost_x','boost_y','boost_z']
        for key in keys_vel + keys_sup_real + keys_sup_imag:
            fits[key] = self._load_scalar_fit(fit_key=key, h5file=h5file)
        return fits

    #-------------------------------------------------------------------------
    def _get_fit_params(self, x, fit_key):
        """ Transforms the input parameter to fit parameters for the 3dq8BMS model.
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
    def mf(self, *args, **kwargs):
        """mf is not implemented in this model. Will return (None, None),"""
        return self._eval_wrapper('mf', *args, **kwargs)
    
    def chif(self, *args, **kwargs):
        """chif is not implemented in this model. Will return (None, None),"""
        return self._eval_wrapper('chif', *args, **kwargs)
    
    def vf(self, *args, **kwargs):
        """vf is not implemented in this model. Will return (None, None),"""
        return self._eval_wrapper('vf', *args, **kwargs)

    def all(self, *args, **kwargs):
        """ Evaluates fit and 1-sigma error estimate for supertranslation parameter 
        alpha, and boost velocity of the BMS transformation from the inspiral BMS 
        frame to the remnant BMS frame.
        Returns:
            alpha, boost,  alpha_err, boost_err
        """
        return self._eval_wrapper('all', *args, **kwargs)

    #-------------------------------------------------------------------------
    def _eval_wrapper(self, fit_key, q, chiA, chiB, **kwargs):
        """ Evaluates the NRSur3dq8BMSRemnant model.
        """

        # mf,chif,vf not implemented for this model
        if fit_key == 'mf':
            mf, mf_err = None, None
            return mf, mf_err
        
        if fit_key == 'chif':
            chif, chif_err = None, None
            return chif, chif_err
        
        if fit_key == 'vf':
            vf, vf_err = None, None
            return vf, vf_err
        
        chiA = np.array(chiA)
        chiB = np.array(chiB)

        # Warn/Exit if extrapolating
        allow_extrap = kwargs.pop('allow_extrap', False)
        self._check_param_limits(q, chiA, chiB, allow_extrap)

        self._check_unused_kwargs(kwargs)

        x = [q, chiA[2], chiB[2]]

        ell_max = 8
        LM_total_size = ell_max * (ell_max + 2) + 1
        alpha = np.zeros((LM_total_size,), dtype=complex)
        alpha_err = np.zeros((LM_total_size,), dtype=complex)
        LMs = [(ell,m) for ell in range(1,ell_max + 1) for m in range(-ell,ell + 1)]
        for i, (ell, m) in enumerate(LMs):
            k = i + 1 #first slot is reserved for ell = 0 mode
            sgn = int(np.sign(m))
            # the fits are trained with the variable S_{ell |m|} = 1/2(alpha_{ell |m|} + (-1)^m alpha_{ell -|m|}) 
            # the variable S_{ell |m|} coincides with alpha_{ell |m|} because alpha is real
            S_real, S_real_err = self._evaluate_fits(x, f'S_{ell}{abs(m)}_real')
            if m == 0:
                alpha[k] = S_real
                alpha_err[k] = S_real_err
                # avoid the else to be faster
                continue
            S_imag, S_imag_err = self._evaluate_fits(x, f'S_{ell}{abs(m)}_imag')
            # general formula alpha_{ell m} = (Re{alpha_{ell |m|}} + i sgn(m) Im{alpha_{ell |m|}}) sgn(m)^m for all ell,m
            alpha[k] = (S_real + sgn * 1j * S_imag) * sgn**int(m)
            alpha_err[k] = (S_real_err + sgn * 1j * S_imag_err) * sgn**int(m)

        boost_x, boost_x_err = self._evaluate_fits(x, 'boost_x')
        boost_y, boost_y_err = self._evaluate_fits(x, 'boost_y')
        boost_z, boost_z_err = self._evaluate_fits(x, 'boost_z')
        boost = np.array([boost_x, boost_y, boost_z])
        boost_err = np.array([boost_x_err, boost_y_err, boost_z_err])
        
        if fit_key == 'all':
            return alpha, boost, alpha_err, boost_err 



