import numpy as np
from surfinBH import surfinBH
import warnings

#=============================================================================
class Fit7dq4Emri(surfinBH.SurFinBH):
    """ A class for the NRSur7dq4EmriRemnant model presented in Boschini et al.,
    arxiv:2307.03435 (hereafter THE PAPER). 
    
    This model predicts the final mass mf and final spin chif, for the remnants 
    of precessing binary black hole systems,  extending to arbitrarily large 
    mass ratios.
    The fits are done using Gaussian Process Regression (GPR) and also
    provide an error estimate along with the fit value.

    This model has been trained in the parameter space:
        NR: q <= 4, |chiA| <= 0.8, |chiB| <= 0.8
        EMRI: 100 <= q <= 1000 , |chiA| <= 1, |chiB| <= 1

    But, the model can be evaluated at arbitrary mass ratios and spins.
    For the remnant mass, the model uses a GPR fit at all mass ratios, which
    includes an error estimate.
    For the remnant spin, the GPR fit only covers mass ratios up to q=1000.
    Therefore, at q<=1000, the error estimate is provided by GPR. At q >= 2000
    the model simply returns the EMRI limit and the error is estimated as the
    absolute difference between this value and chiA (which is the limit of chif
    at q -> inf). There is a transition region (1000 < q < 2000) to smoothly
    connect these two regimes, and the error estimate also has a similar
    transition. See Sec. III in THE PAPER for details.

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

        # All of these together (this fit returns vf = vf_err = None)
        # This output is needed for compatibility with the other fits
        mf, chif, vf, mf_err, chif_err, vf_err =
            fit.all(q, chiA, chiB, **kwargs)

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
                For comparable-mass binaries (q<=6) the reference epoch is
                assumed to be at t=-100 M from the peak of the waveform
                amplitude. In the EMRI limit (q>=100) the reference epoch is the
                ISCO (innermost stable circular orbit). In the intermediate
                regime the reference epoch is a transition between these two
                frames. See Sec. III in THE PAPER for additional details.

     Inertial frame for returned values:

        The returned chif are in the LAL inertial frame defined as follows:
            The +ve z-axis is along the orbital angular momentum at the
            reference epoch. The separation vector from the lighter BH to the
            heavier BH at the reference epoch is along the +ve x-axis. The
            y-axis completes the right-handed triad.

            Note that the default reference epoch corresponds to t=-100M for
            comparable-mass binaries, the ISCO for EMRI binaries, and an
            intemediate epoch for binaries in between the two limits.
            See LIGO DCC document T1800226 for the LAL frame diagram.
    """

    #-------------------------------------------------------------------------
    def __init__(self, name):

        # Param limits beyond which to raise a warning
        # With current settings a warning is never raised
        soft_param_lims = {
            'q': np.Inf,
            'chiAmag': 1,
            'chiBmag': 1,
                }

        # Param limits beyond which to raise an error
        # With current settings an error is never raised
        hard_param_lims = {
            'q' : np.Inf,
            'chiAmag': 1,
            'chiBmag': 1,
                }
        super(Fit7dq4Emri, self).__init__(name, soft_param_lims, \
            hard_param_lims)

        # The upper bound of the EMRI training region. It defines when the model
        # switches from the GPR fit to the EMRI limit (with a smooth transition)
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
    def _generate_random_params_for_tests(self):
        """ Generate random parameters to use in tests."""

        # Generate params randomly within allowed values. NRSur7dq4EmriRemnant
        # is valid at arbitrary mass ratios, but here an upper value is needed.
        q = 10**np.random.uniform(0, 10)
        chiAmag= np.random.uniform(0, self.hard_param_lims['chiAmag'])
        chiBmag= np.random.uniform(0, self.hard_param_lims['chiBmag'])
        if self.aligned_spin_only:
            chiAph, chiBph = 0, 0
            chiAth, chiBth = np.random.choice([0, np.pi]), \
                np.random.choice([0, np.pi])
        else:
            chiAth = np.arccos(np.random.uniform(-1., 1.))
            chiBth = np.arccos(np.random.uniform(-1., 1.))
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
    def _eval_wrapper(self, fit_key, q, chiA, chiB, **kwargs):
        """ Evaluates the NRSur7dq4EmriRemnant model. """

        chiA = np.array(chiA)
        chiB = np.array(chiB)

        # Warn if trying to pass allow_extrap
        allow_extrap = kwargs.pop('allow_extrap', False)
        if allow_extrap:
            warnings.warn('Optional argument allow_extrap is unused for this'
                          'model. It works at arbitrary mass ratios and spins.')
        self._check_param_limits(q, chiA, chiB, True)

        x = np.concatenate(([q], chiA, chiB))

        def eval_r_isco(chiz):
            """
            Estimate ISCO radius following Eq. (2.21) by Bardeen et al.
            (https://doi.org/10.1086/151796)
            """

            Z1 = 1+(1-chiz**2)**(1/3)*((1+chiz)**(1/3)+(1-chiz)**(1/3))
            Z2 = (3*chiz**2+Z1**2)**(1/2)
            r_ISCO = 3+Z2-np.sign(chiz)*((3-Z1)*(3+Z1+2*Z2))**(1/2)

            return r_ISCO

        def eval_chif_emri(q, chiA):
            """
            Compute remnant spin for EMRI binary following Eq. (7)
            presented in THE PAPER.
            """

            r_ISCO = eval_r_isco(chiA[2])
            E_ISCO = np.sqrt(1-2/(3*r_ISCO))
            L_ISCO = np.array([0,0,2/(3*3**(1/2))*(1+2*np.sqrt(3*r_ISCO-2))])

            return chiA + 1/q*(L_ISCO-2*chiA*E_ISCO)

        def eval_chif(x, fit_key):
            """
            The model evaluates remnant spin using GPR fit up to the upper bound
            of the EMRI training region (q = 1000). At q > 1000, it returns the
            expected EMRI limit (Eq. (7) in THE PAPER). There is a transition
            function (Eq. (9) in THE PAPER) connecting smoothly these
            two regions. The transition takes place from q=1000 to q=2000.
            """

            if x[0] <= self.qEmriMax:
                # From the GPR fit
                res = self._evaluate_fits(x, fit_key)
                y, y_err = res.T[0], res.T[1]

            #This if condition defines the width of the transition region
            elif x[0] < 2*self.qEmriMax:
                x_i = x.copy()
                x_i[0] = self.qEmriMax
                # GPR fit is used to compute chif at the lower end of the
                # transition region (q=1000)
                res = self._evaluate_fits(x_i, fit_key)
                y_i, y_i_err = res.T[0], res.T[1]
                # EMRI limit is used to compute chif at the upper end of the
                # transition region (q=2000)
                y_f = eval_chif_emri(2*self.qEmriMax, x[1:4])
                # Error estimate when evaluating the EMRI limit.
                # The model returns the difference between remnant spin at
                # q = 2000 and the limit value (chif->chiA when q->inf)
                # as an estimate of the error.
                y_f_err = np.abs(x[1:4] - y_f)
                # transition function to smoothly connect GPR fit and EMRI limit
                g = np.sin(np.pi/(2*self.qEmriMax) * (x[0] - self.qEmriMax))**2
                y = (1 - g) * y_i + g * y_f
                # Adding errors in quadrature
                y_err = np.sqrt((1 - g)**2 * y_i_err**2 + g**2 * y_f_err**2)

            else:
                y = eval_chif_emri(x[0], x[1:4])
                # Error estimate when evaluating the EMRI limit.
                # The model returns the difference between remnant spin at a
                # given q (>2000) and the limit value (chif->chiA when q->inf)
                # as an estimate of the error.
                y_err = np.abs(x[1:4] - y)

            return y, y_err

        def eval_vector_fit(x, fit_key):
            if fit_key == 'chif':
                fit_val, fit_err = eval_chif(x, fit_key)

            return fit_val, fit_err

        # The remnant mass fit is enforced through a mapping procedure in
        # order to include physical constraint (mf/M < 1, M = m1+m2). The fit is
        # trained passing arctanh(mf) and it returs the tanh, so the values are
        # constrained in the interval [0,1). See Sec. III in THE PAPER for
        # additional details.
        if fit_key == 'mf' or fit_key == 'all':
            mf, mf_err = self._evaluate_fits(x, 'mf')
            mf, mf_err = np.tanh(mf), 1/np.cosh(mf)**2*mf_err
            if fit_key == 'mf':
                return mf, mf_err

        if fit_key == 'chif' or fit_key == 'all':
            chif, chif_err = eval_vector_fit(x, 'chif')
            if fit_key == 'chif':
                return chif, chif_err

        if fit_key == 'vf' or fit_key == 'all':
            vf, vf_err = None, None
            if fit_key == 'vf':
                return vf, vf_err

        if fit_key == 'all':
            return mf, chif, vf, mf_err, chif_err, vf_err
