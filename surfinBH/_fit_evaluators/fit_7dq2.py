import numpy as np
from surfinBH import surfinBH
from surfinBH._lal_spin_evolution import evolve_pn_spins
import surfinBH._utils as utils
import warnings

#=============================================================================
class Fit7dq2(surfinBH.SurFinBH):
    """ A class for the surfinBH7dq2 model presented in Varma et al., 2018,
    in prep. This model predicts the final mass mC, final spin vector chiC and
    final kick velocity vector velC, for the remnants of precessing binary
    black hole systems.  The fits are done using Gaussian Process Regression
    (GPR) and also provide an error estimate along with the fit value.

    This model has been trained in the parameter space:
        q <= 2, |chiA| <= 0.8, |chiB| <= 0.8

    However, it extrapolates reasonably to:
        q <= 3, |chiA| <= 1, |chiB| <= 1

    See __call__ method for evaluation.
    """

    #-------------------------------------------------------------------------
    def __init__(self, name, load_nrsur=False):

        #NOTE: These are not the actual limits.
        # We override _check_param_limits() to set limits, the
        # only purpose these serve is to generate regression data in
        # surfinBH/test/generate_regression_data.py
        soft_param_lims = None
        hard_param_lims = [[0.99, 3.01],
                [-0.6, 0.6],
                [-0.6, 0.6],
                [-0.6, 0.6],
                [-0.6, 0.6],
                [-0.6, 0.6],
                [-0.6, 0.6]]

        super(Fit7dq2, self).__init__(name, soft_param_lims, hard_param_lims)
        self.nrsur = None

    def _load_NRSur7dq2(self):
        import NRSur7dq2
        self.nrsur = NRSur7dq2.NRSurrogate7dq2()
        print('Loaded NRSur7dq2 waveform model')

    #-------------------------------------------------------------------------
    def _load_fits(self, h5file):
        """ Loads fits from h5file and returns a dictionary of fits. """
        fits = {}
        for key in ['mC']:
            fits[key] = self._load_scalar_fit(fit_key=key, h5file=h5file)
        for key in ['chiC', 'velC']:
            fits[key] = self._load_vector_fit(key, h5file)
        return fits

    #-------------------------------------------------------------------------
    def _get_fit_params(self, x, fit_key):
        """ Transforms the input parameter to fit parameters for the 7dq2 model.
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

        fit_params = x
        fit_params[0] = np.log(q)
        fit_params[3] = chiHat
        fit_params[6] = chi_a

        return fit_params

    #-------------------------------------------------------------------------
    def _check_param_limits(self, x):
        """ Checks that x is within allowed range of paramters.
        Raises a warning if outside training limits and
        raises an error if outside allowed limits.
        Training limits: q <= 2.01, chiAmag <= 0.81, chiBmag <= 0.81.
        Allowed limits: q <= 3.01, chiAmag <= 1, chiBmag <= 1.
        """
        q = x[0]
        chiAmag = np.sqrt(np.sum(x[1:4]**2))
        chiBmag = np.sqrt(np.sum(x[4:7]**2))

        if q < 1:
            raise ValueError('Mass ratio should be >= 1.')
        elif q > 3.01:
            raise Exception('Mass ratio outside allowed range.')
        elif q > 2.01:
            warnings.warn('Mass ratio outside training range.')

        if chiAmag > 1.:
            raise Exception('Spin magnitude of BhA outside allowed range.')
        elif chiAmag > 0.81:
            warnings.warn('Spin magnitude of BhA outside training range.')

        if chiBmag > 1.:
            raise Exception('Spin magnitude of BhB outside allowed range.')
        elif chiBmag > 0.81:
            warnings.warn('Spin magnitude of BhB outside training range.')

    def _evolve_spins(self, q, chiA0, chiB0, omega0, PN_approximant,
            PN_dt, PN_spin0, PN_phase0):
        """ Evolves spins of the component BHs from an initial orbital
        frequency = omega0 until t=-100 M from the peak of the waveform.
        If omega0 < 0.018, use PN to evolve the spins until
        orbital frequency = omega0. Then evolves further with the NRSur7dq2
        waveform model until t=-100M from the peak.

        Assumes chiA0 and chiB0 are defined in the inertial frame defined
        at orbital frequency = omega0 as:
            The z-axis is along the Newtonian orbital angular momentum when the
                PN orbital frequency = omega0.
            The x-axis is along the line of separation from the smaller BH to
                the larger BH at this frequency.
            The y-axis completes the triad.

        Returns spins in the coorbital frame at t=-100M, as well as the
        coprecessing frame quaternion and orbital phase in the coprecessing
        frame at this time.
        """

        # obrbital frequency beyond which we use the NRSur7dq2 model.
        omega0_nrsur = 0.018

        # If omega0 is below the NRSur7dq2 start frequency, we use PN
        # to evolve the spins until orbital frequency = omega0_nrsur.
        if omega0 < omega0_nrsur:
            # Note that we update omega0_nrsur here with the PN
            # frequency that was closest to the input omega0_nrsur.
            chiA0_nrsur_copr, chiB0_nrsur_copr, quat0_nrsur_copr, phi0_nrsur, \
                omega0_nrsur \
                = evolve_pn_spins(q, chiA0, chiB0, omega0,
                    omega0_nrsur, approximant=PN_approximant,
                    dt=PN_dt, spinO=PN_spin0,
                    phaseO=PN_phase0)

        # Load NRSur7dq2 if needed
        if self.nrsur is None:
            self._load_NRSur7dq2()

        # evaluate NRSur7dq2 dynamics
        quat, orbphase, chiA_copr, chiB_copr = self.nrsur.get_dynamics(q,
            chiA0_nrsur_copr, chiB0_nrsur_copr, init_quat=quat0_nrsur_copr,
            init_phase=phi0_nrsur, omega_ref=omega0_nrsur)

        # get data at time node where remnant fits are done
        fitnode_time = -100
        nodeIdx = np.argmin(np.abs(self.nrsur.tds - fitnode_time))
        quat_fitnode = quat.T[nodeIdx]
        orbphase_fitnode = orbphase[nodeIdx]

        # get coorbital frame spins at the time node
        chiA_coorb_fitnode = utils.rotate_in_plane(chiA_copr[nodeIdx],
                orbphase_fitnode)
        chiB_coorb_fitnode = utils.rotate_in_plane(chiB_copr[nodeIdx],
                orbphase_fitnode)

        return chiA_coorb_fitnode, chiB_coorb_fitnode, quat_fitnode, \
            orbphase_fitnode


    #-------------------------------------------------------------------------
    def __call__(self, fit_key, x, **kwargs):
        """Evaluates the surfinBH7dq2 model.

    Arguments:
        fit_key:
            'mC', 'chiC' or 'velC', for remnant mass, spin vector and kick
            vector respectively.

        x:
            Array of [q, chiAx, chiAy, chiAz, chiBx, chiBy, chiBz], where q >=
            1 is the mass ratio and chiA (chiB) is the dimensionless spin of
            the larger (smaller) BH.

            By default, the spins are assumed to be the component spins at
            t=-100 M from the peak of the waveform, and in the coorbital frame,
            defined as:
            The z-axis is along the orbital angular momentum at t=-100M.
            The x-axis is along the line of separation from the smaller BH to
                the larger BH at this time.
            The y-axis completes the triad.
            We obtain this frame from the waveform as defined in
            arxiv:1705.07089.

            If 'omega0' is given, instead the spins are assumed to be the
            component spins when the PN orbital frequency = omega0. The
            spins are assumed to be in the inertial frame, defined as:
            The z-axis is along the Newtonian orbital angular momentum when the
                PN orbital frequency = omega0.
            The x-axis is along the line of separation from the smaller BH to
                the larger BH at this frequency.
            The y-axis completes the triad.
            We obtain this frame from PN.
            Given the spins at omega0, we evolve the spins using PN until
            the orbital frequency = 0.018 M, then we further evolve the spins
            using the NRSur7dq2 model (arxiv:1705.07089) until t=-100M from the
            peak of the waveform.  Then we evaluate the fits using the spins at
            t=-100 M.  Finally, we transform the remnant spin and kick vectors
            back to the inertial frame defined above.

    Optional PN arguments:

        omega0:
            Initial dimensionless orbital frequency in units of 1/M, where M is
            the total mass. If this is given, the spins in x are assumed to be
            the component spins at this orbital frequency, and in the inertial
            frame as defined above. The returned remnant spin and kick vectors
            are also in the same frame.
            Default: None.

        PN_approximant:
            Approximant used to do the PN spin evolution. Choose from
            'SpinTaylorT1', 'SpinTaylorT4' or 'SpinTaylorT5'.
            Default: 'SpinTaylorT4'.

        PN_dt:
            Dimensionless time step size in units of M (total mass), used for
            the PN evolution. You may need to increase this if omega0 is very
            low.
            Default: 0.1

        PN_spin_order:
            Twice the PN order of spin effects. E.g., use 7 for 3.5PN.
            Default: 7

        PN_phase_order
            Twice the PN order in phase. E.g., use 7 for 3.5PN.
            Default: 7

    Returns:
        If fit_key='mC':
            returns mC, mC_err_est
            The value and 1-sigma error estimate in remnant mass.
        If fit_key='chiC':
            returns chiC, chiC_err_est
            The value and 1-sigma error estimate in remnant spin vector.
        If fit_key='velC':
            returns velC, velC_err_est
            The value and 1-sigma error estimate in remnant kick vector.

        By default, these vectors are defined in the coorbital frame at
            t=-100 M, as described above.
        If 'omega0' is not None, these are defined in the inertial frame at
            orbital frequency = omega0, as described above.
        """

        x = np.array(x)

        # Warn/Exit if extrapolating
        self._check_param_limits(x)

        omega0 = kwargs.pop('omega0', None)
        PN_approximant = kwargs.pop('PN_approximant', 'SpinTaylorT4')
        PN_dt = kwargs.pop('PN_dt', 0.1)
        PN_spin_order = kwargs.pop('PN_spin_order', 7)
        PN_phase_order = kwargs.pop('PN_phase_order', 7)

        # If omega0 is given, evolve the spins from omega0
        # to t = -100 M from the peak. Replace x with these spins.
        if omega0 is not None:
            q = x[0]
            chiA0 = x[1:4]
            chiB0 = x[4:7]
            chiA_coorb_fitnode, chiB_coorb_fitnode, quat_fitnode, \
                orbphase_fitnode \
                = self._evolve_spins(q, chiA0, chiB0, omega0,
                    PN_approximant, PN_dt, PN_spin_order,
                    PN_phase_order)

            # Update x to contain coorbital frame spins at t=-100M
            x = np.concatenate(([q], chiA_coorb_fitnode, chiB_coorb_fitnode))

        if fit_key == 'mC':
            fit_val, fit_err = self._evaluate_fits(x, fit_key)
        elif fit_key == 'chiC' or fit_key == 'velC':
            res = self._evaluate_fits(x, fit_key)
            fit_val = res.T[0]
            fit_err = res.T[1]
            if omega0 is not None:
                # If spins were given in inertial frame at omega0,
                # transform vectors and errors back to the same frame.
                fit_val = utils.transform_vector_coorb_to_inertial(fit_val,
                    orbphase_fitnode, quat_fitnode)
                fit_err = utils.transform_error_coorb_to_inertial(fit_val,
                    fit_err, orbphase_fitnode, quat_fitnode)
        else:
            raise ValueError('Invalid fit_key')

        if len(kwargs.keys()) != 0:
            unused = ""
            for k in kwargs.keys():
                unused += "'%s', "%k
            if unused[-2:] == ", ":     # get rid of trailing comma
                unused = unused[:-2]
            raise Exception('Unused keys in kwargs: %s'%unused)

        return fit_val, fit_err
