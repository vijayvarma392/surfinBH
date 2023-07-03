import numpy as np
from surfinBH import surfinBH
import warnings

#=============================================================================
class Fit7dq4Emri(surfinBH.SurFinBH):
    """ A class for the NRSur7dq4EmriRemnant model presented in Boschini et al.,
    arxiv:[] (hereafter THE PAPER). 
    
    This model predicts the final mass mf and final spin chif, for the remnants 
    of precessing binary black hole systems. 
    The fits are done using Gaussian Process Regression (GPR) and also 
    provide an error estimate along with the fit value.

    This model has been trained in the parameter space:
        q <= 4 and 100 <= q <= 1000, |chiAz| <= 0.8, |chiBz| <= 0.8

    It interpolates remnant properties at:
        4 < q < 100

    It extrapolates reasonably to:
        q > 1000 , |chiAz| <= 1, |chiBz| <= 1

    =========================================================================
    Usage:

    import surfinBH

    # Load the fit
    fit = surfinBH.LoadFits('NRSur7dq4EmriRemnant')

    We provide the following call methods:
        # remnant mass and 1-sigma error estimate
        mf, mf_err = fit.mf(q, chiA, chiB, **kwargs)

        # remnant spin and 1-sigma error estimate
        chif, chif_err = fit.chif(q, chiA, chiB, **kwargs)

    The arguments for each of these call methods are as follows:
    Arguments:
        q:     Mass ratio (q = mA/mB >=1)

        chiA:  Dimensionless spin vector of the heavier black hole at
                reference epoch.

        chiB:  Dimensionless spin vector of the lighter black hole at
                reference epoch.

                This follows the same convention as LAL, where the spin
                components are defined as:
                \chi_z = \chi \cdot \hat{L}, where L is the orbital angular
                    momentum vector at the epoch.
                \chi_x = \chi \cdot \hat{n}, where n = body2 -> body1 is the
                    separation vector at the epoch. body1 is the heavier body.
                \chi_y = \chi \cdot \hat{L \cross n}.
                These spin components are frame-independent as they are defined
                using vector inner products. This is equivalent to specifying
                the spins in the coorbital frame at the reference epoch. 
                The reference epoch is assumed to be at t=-100 M from the peak 
                of the waveform, see THE PAPER for additional details.

    Optional arguments:
        allow_extrap:
            If False, raises a warning when q > 1000.1 or|chiA|,|chiB| > 0.81,
                and raises an error when |chiA|,|chiB| > 1.
            If True, allows extrapolation to any q and |chiA|,|chiB| <= 1.
                Use at your own risk.
            Default: False.

    Inertial frame for returned values:

        The returned chif are in the LAL inertial frame defined as follows:
            The +ve z-axis is along the orbital angular momentum at the
            reference epoch. The separation vector from the lighter BH to the
            heavier BH at the reference epoch is along the +ve x-axis. The
            y-axis completes the right-handed triad.

            Note that the default reference epoch corresponds to t=-100M.
            See LIGO DCC document T1800226 for the LAL frame diagram.
    """

    #-------------------------------------------------------------------------
    def __init__(self, name):

        # Param limits beyond which to raise a warning
        soft_param_lims = {
            'q': 1000.1,
            'chiAmag': 0.81,
            'chiBmag': 0.81,
                }

        # Param limits beyond which to raise an error
        hard_param_lims = {
            'chiAmag': 1,
            'chiBmag': 1,
                }
        super(Fit7dq4Emri, self).__init__(name, soft_param_lims, \
            hard_param_lims)
        self.qEmriMax = 1000

    #-------------------------------------------------------------------------
    def _load_fits(self, h5file):
        """ Loads fits from h5file and returns a dictionary of fits. """
        fits = {}
        for key in ['mf']:
            fits[key] = self._load_scalar_fit(fit_key=key, h5file=h5file)
        for key in ['chif']:
            fits[key] = self._load_vector_fit(key, h5file)
        return fits

    #-------------------------------------------------------------------------
    def _get_fit_params(self, x, fit_key):
        """ Transforms the input parameter to fit parameters for the 7dq4Emri 
    model.
    That is, maps from
    x = [q, chiAx, chiAy, chiAz, chiBx, chiBy, chiBz]
    fit_params = [np.log(q), chiAx, chiAy, chiHat, chiBx, chiBy, chi_a]

    chiHat is defined in Eq.(3) of 1508.07253, but with chiAz and chiBz instead
    of chiA and chiB.
    chi_a = (chiAz - chiBz)/2.
        """
        q, chiAz, chiBz = x[0], x[3], x[6]
        eta = q/(1.+q)**2
        chi_wtAvg = (q*chiAz+chiBz)/(1.+q)
        chiHat = (chi_wtAvg - 38.*eta/113.*(chiAz + chiBz))/(1. - 76.*eta/113.)
        chi_a = (chiAz - chiBz)/2.

        fit_params = [np.log(q), x[1], x[2], chiHat, x[4], x[5], chi_a]
        return fit_params

    #-------------------------------------------------------------------------
    def _eval_wrapper(self, fit_key, q, chiA, chiB, **kwargs):
        """ Evaluates the NRSur7dq4EmriRemnant model.
        """
        chiA = np.array(chiA)
        chiB = np.array(chiB)

        # Warn/Exit if extrapolating
        allow_extrap = kwargs.pop('allow_extrap', False)
        self._check_param_limits(q, chiA, chiB, allow_extrap)

        x = np.concatenate(([q], chiA, chiB))

        def eval_r_isco(chiz):
            """
            Estimate ISCO radius following Bardeen et al. 
            (https://doi.org/10.1086/151796)
            """
 
            Z1 = 1+(1-chiz**2)**(1/3)*((1+chiz)**(1/3)+(1-chiz)**(1/3))
            Z2 = (3*chiz**2+Z1**2)**(1/2)
            r_ISCO = 3+Z2-np.sign(chiz)*((3-Z1)*(3+Z1+2*Z2))**(1/2)
    
            return r_ISCO

        def eval_emri(q, chiA):
            """
            Estimate remnant mass and spin for EMRI binary following 
            Boschini & al. (). Additional details in THE PAPER.    
            """
            
            E_ISCO = np.sqrt(1-2/(3*eval_r_isco(chi[:,2])))
            L_ISCO = np.array([0,0,2/(3*3**(1/2))*(1+2*np.sqrt(3*r_ISCO-2))])
                      
            return chiA + 1/q*(L_ISCO-2*chiA*E_ISCO)  

        def get_extrap(x, chiA):
            """Extrapolate remnant spin at q > 1000"""
            x_i = x.copy()
            x_i[0] = np.log(self.qEmriMax)

            if x[0] < np.log(2*self.qEmriMax):
                tmp = self._evaluate_fits(x, fit_key)
                y_i, y_i_err = tmp.T[0], tmp.T[1]
                y_f = eval_emri(np.exp(x[0]), chiA)
                y_f_err = 2*np.abs(chiA - y_f)
                y = [(y_f[i]-y_i[i])*np.sin(np.pi/(2*self.qEmriMax)* \
                    (np.exp(x[0])-self.qEmriMax))**2+y_i[i] for i in range(3)]
                y_err = [(y_f_err[i]-y_i_err[i])* \
                    np.sin(np.pi/(2*self.qEmriMax)* \
                    (np.exp(x[0])-self.qEmriMax))**2+ \
                    y_i_err[i] for i in range(3)]
            else:
                y = eval_emri(np.exp(x[0]), chiA)
                y_err = 2*np.abs(chiA - y_f)

            return np.array(y), np.array(y_err)

        def eval_vector_fit(x, fit_key, chiA):
            if fit_key == 'chif':
                if np.exp(x[0]) > self.qEmriMax:
                    fit_val, fit_err = get_extrap(x, chiA)
                else:
                    res = self._evaluate_fits(x, fit_key)
                    fit_val = res.T[0]
                    fit_err = res.T[1]
            
            return fit_val, fit_err

        if fit_key == 'mf' or fit_key == 'all':
            mf, mf_err = self._evaluate_fits(x, 'mf')
            mf, mf_err = np.tanh(mf), 1/np.cosh(np.arctanh(mf))**2*mf_err
            if fit_key == 'mf':
                return mf, mf_err

        if fit_key == 'chif' or fit_key == 'all':
            chif, chif_err = eval_vector_fit(x, 'chif', chiA)
            if fit_key == 'chif':
                return chif, chif_err

        if fit_key == 'all':
            return mf, chif, mf_err, chif_err
