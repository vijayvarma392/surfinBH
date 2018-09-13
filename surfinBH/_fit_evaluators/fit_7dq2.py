import numpy as np
from surfinBH import surfinBH
from surfinBH._lal_spin_evolution import evolve_pn_spins
import surfinBH._utils as utils
import warnings

#=============================================================================
class Fit7dq2(surfinBH.SurFinBH):
    """ A class for the surfinBH7dq2 model presented in Varma et al., 2018,
    in prep. This model predicts the final mass mf, final spin vector chif and
    final kick velocity vector vf, for the remnants of precessing binary
    black hole systems.  The fits are done using Gaussian Process Regression
    (GPR) and also provide an error estimate along with the fit value.

    This model has been trained in the parameter space:
        q <= 2, |chiA| <= 0.8, |chiB| <= 0.8

    However, it extrapolates reasonably to:
        q <= 4, |chiA| <= 1, |chiB| <= 1

    =========================================================================
    Usage:

    import surfinBH

    # Load the fit
    fit = surfinBH.LoadFits('surfinBH7dq2')

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

        chiA:   Dimensionless spin of the larger BH (array of size 3)

        chiB:   Dimensionless spin of the smaller BH (array of size 3)

            By default, the spins are assumed to be the component spins at
            t=-100 M from the peak of the waveform, and in the coorbital frame,
            defined as:
            The z-axis is along the orbital angular momentum at t=-100M.
            The x-axis is along the line of separation from the smaller BH to
                the larger BH at this time.
            The y-axis completes the triad.
            We obtain this frame from the waveform as defined in
            arxiv:1705.07089.
            The returned spin and kick vectors are also in the same frame.

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
            the orbital frequency = omega_switch, then we further evolve the
            spins using the NRSur7dq2 model (arxiv:1705.07089) until t=-100M
            from the peak of the waveform. Then we evaluate the fits using the
            spins at t=-100 M. Finally, we transform the remnant spin and kick
            vectors back to the inertial frame defined above.

    Optional arguments:
        allow_extrap:
            If False, raises a warning when q > 2.1 or |chiA|,|chiB| > 0.81,
                and raises an error when q > 4.1 or |chiA|,|chiB| > 1.
            If True, allows extrapolation to any q and |chiA|,|chiB| <= 1.
                Use at your own risk.
            Default: False.

        omega0:
            Initial dimensionless orbital frequency in units of 1/M, where M is
            the total mass. If this is given, the spins in x are assumed to be
            the component spins at this orbital frequency, and in the inertial
            frame as defined above. The returned remnant spin and kick vectors
            are also in the same frame.
            Default: None.

        PN_approximant:
            Approximant used to do the PN spin evolution. Choose from
            'SpinTaylorT4', 'SpinTaylorT1' or 'SpinTaylorT2'.
            Default: 'SpinTaylorT4'.

        PN_dt:
            Dimensionless time step size in units of M (total mass), used for
            the PN evolution. You may need to increase this if omega0 is very
            low.
            Default: 0.1

        PN_spin_order:
            Twice the PN order of spin effects. E.g., use 7 for 3.5PN.
            Default: 7

        PN_phase_order:
            Twice the PN order in phase. E.g., use 7 for 3.5PN.
            Default: 7

        omega_switch:
            Dimensionless orbital frequency at which to switch from PN to
            NRSur7dq2 model. You may need to increase this if the NRSur7dq2
            model raises an exception like:
            "Got omega_ref = 0.0180 < 0.0184 = omega_0, too small!"
            Default: 0.018
    """

    #-------------------------------------------------------------------------
    def __init__(self, name, load_nrsur=False):

        # Param limits beyond which to raise a warning
        soft_param_lims = {
            'q': 2.1,
            'chiAmag': 0.81,
            'chiBmag': 0.81,
                }

        # Param limits beyond which to raise an error
        hard_param_lims = {
            'q': 4.1,
            'chiAmag': 1,
            'chiBmag': 1,
                }

        super(Fit7dq2, self).__init__(name, soft_param_lims, hard_param_lims)
        self.nrsur = None

    #-------------------------------------------------------------------------
    def _load_NRSur7dq2(self):
        import NRSur7dq2
        self.nrsur = NRSur7dq2.NRSurrogate7dq2()
        print('Loaded NRSur7dq2 waveform model')

    #-------------------------------------------------------------------------
    def _load_fits(self, h5file):
        """ Loads fits from h5file and returns a dictionary of fits. """
        fits = {}
        for key in ['mf']:
            fits[key] = self._load_scalar_fit(fit_key=key, h5file=h5file)
        for key in ['chif', 'vf']:
            fits[key] = self._load_vector_fit(key, h5file)
        return fits

    #-------------------------------------------------------------------------
    def _extra_regression_kwargs(self):
        """ List of additional kwargs to use in regression tests.
        """

        # larger than default sometimes needed when extrapolating
        omega_switch_test = 0.019

        extra_args = []
        extra_args.append({
            'omega0': 5e-3,
            'PN_approximant': 'SpinTaylorT4',
            'PN_dt': 0.1,
            'PN_spin_order': 7,
            'PN_phase_order': 7,
            'omega_switch': omega_switch_test,
            })


        extra_args.append({
            'omega0': 6e-3,
            'PN_approximant': 'SpinTaylorT1',
            'PN_dt': 0.5,
            'PN_spin_order': 5,
            'PN_phase_order': 7,
            'omega_switch': omega_switch_test,
            })

        extra_args.append({
            'omega0': 7e-3,
            'PN_approximant': 'SpinTaylorT2',
            'PN_dt': 1,
            'PN_spin_order': 7,
            'PN_phase_order': 5,
            'omega_switch': omega_switch_test,
            })

        # These should be pure NRSur7dq2
        extra_args.append({'omega0': 3e-2})
        extra_args.append({'omega0': 5e-2})

        return extra_args

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
    def _evolve_spins(self, q, chiA0, chiB0, omega0, PN_approximant,
            PN_dt, PN_spin0, PN_phase0, omega0_nrsur):
        """ Evolves spins of the component BHs from an initial orbital
        frequency = omega0 until t=-100 M from the peak of the waveform.
        If omega0 < omega0_nrsur, use PN to evolve the spins until
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

        if omega0 < omega0_nrsur:
            # If omega0 is below the NRSur7dq2 start frequency, we use PN
            # to evolve the spins until orbital frequency = omega0_nrsur.

            # Note that we update omega0_nrsur here with the PN
            # frequency that was closest to the input omega0_nrsur.
            chiA0_nrsur_copr, chiB0_nrsur_copr, quat0_nrsur_copr, \
                phi0_nrsur, omega0_nrsur \
                = evolve_pn_spins(q, chiA0, chiB0, omega0,
                    omega0_nrsur, approximant=PN_approximant,
                    dt=PN_dt, spinO=PN_spin0,
                    phaseO=PN_phase0)
        else:
            # If omega0>= omega0_nrsur, we evolve spins directly with NRSur7dq2
            # waveform model. We set the coprecessing frame quaternion to
            # identity and orbital phase to 0 at omega=omega0, hence the
            # coprecessing frame is the same as the inertial frame here.

            # Note that we update omega0_nrsur here and set it to omega0
            chiA0_nrsur_copr, chiB0_nrsur_copr, quat0_nrsur_copr, \
                phi0_nrsur, omega0_nrsur \
                = chiA0, chiB0, [1,0,0,0], 0, omega0

        # Load NRSur7dq2 if needed
        if self.nrsur is None:
            self._load_NRSur7dq2()

        # evaluate NRSur7dq2 dynamics
        # We set allow_extrapolation=True always since we test param limits
        # independently
        quat, orbphase, chiA_copr, chiB_copr = self.nrsur.get_dynamics(q,
            chiA0_nrsur_copr, chiB0_nrsur_copr, init_quat=quat0_nrsur_copr,
            init_phase=phi0_nrsur, omega_ref=omega0_nrsur,
            allow_extrapolation=True)

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
    def _eval_wrapper(self, fit_key, q, chiA, chiB, **kwargs):
        """Evaluates the surfinBH7dq2 model.
        """
        chiA = np.array(chiA)
        chiB = np.array(chiB)

        # Warn/Exit if extrapolating
        allow_extrap = kwargs.pop('allow_extrap', False)
        self._check_param_limits(q, chiA, chiB, allow_extrap)

        omega0 = kwargs.pop('omega0', None)
        PN_approximant = kwargs.pop('PN_approximant', 'SpinTaylorT4')
        PN_dt = kwargs.pop('PN_dt', 0.1)
        PN_spin_order = kwargs.pop('PN_spin_order', 7)
        PN_phase_order = kwargs.pop('PN_phase_order', 7)
        omega_switch = kwargs.pop('omega_switch', 0.018)

        self._check_unused_kwargs(kwargs)

        if omega0 is None:
            # If omega0 is given, assume chiA, chiB are the coorbital frame
            # spins at t=-100 M.
            x = np.concatenate(([q], chiA, chiB))
        else:
            # If omega0 is given, evolve the spins from omega0
            # to t = -100 M from the peak.
            chiA_coorb_fitnode, chiB_coorb_fitnode, quat_fitnode, \
                orbphase_fitnode \
                = self._evolve_spins(q, chiA, chiB, omega0,
                    PN_approximant, PN_dt, PN_spin_order,
                    PN_phase_order, omega_switch)
            # x should contain coorbital frame spins at t=-100M
            x = np.concatenate(([q], chiA_coorb_fitnode, chiB_coorb_fitnode))


        def eval_vector_fit(x, fit_key):
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
            return fit_val, fit_err


        if fit_key == 'mf' or fit_key == 'all':
            mf, mf_err = self._evaluate_fits(x, 'mf')
            if fit_key == 'mf':
                return mf, mf_err

        if fit_key == 'chif' or fit_key == 'all':
            chif, chif_err = eval_vector_fit(x, 'chif')
            if fit_key == 'chif':
                return chif, chif_err

        if fit_key == 'vf' or fit_key == 'all':
            vf, vf_err = eval_vector_fit(x, 'vf')
            if fit_key == 'vf':
                return vf, vf_err

        if fit_key == 'all':
            return mf, chif, vf, mf_err, chif_err, vf_err
