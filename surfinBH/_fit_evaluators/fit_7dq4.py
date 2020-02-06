import numpy as np
import sys
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from surfinBH import surfinBH
from surfinBH._lal_spin_evolution import lal_spin_evloution_wrapper
import surfinBH._utils as utils
import warnings

#=============================================================================
class Fit7dq4(surfinBH.SurFinBH):
    """ A class for the NRSur7dq4Remnant model presented in Varma et al.,
    arxiv:1905.09300, hereafter referred to as THE PAPER.

    This model predicts the final mass mf, final spin vector
    chif and final kick velocity vector vf, for the remnants of precessing
    binary black hole systems.  The fits are done using Gaussian Process
    Regression (GPR) and also provide an error estimate along with the fit
    value.

    This model has been trained in the parameter space:
        q <= 4, |chiA| <= 0.8, |chiB| <= 0.8

    However, it extrapolates reasonably to:
        q <= 6, |chiA| <= 1, |chiB| <= 1

    =========================================================================
    Usage:

    import surfinBH

    # Load the fit
    fit = surfinBH.LoadFits('NRSur7dq4Remnant')

    We provide the following call methods:
        # remnant mass and 1-sigma error estimate
        mf, mf_err = fit.mf(q, chiA, chiB, **kwargs)

        # remnant spin and 1-sigma error estimate
        chif, chif_err = fit.chif(q, chiA, chiB, **kwargs)

        # remnant recoil kick and 1-sigma error estimate (units of c)
        vf, vf_err = fit.vf(q, chiA, chiB, **kwargs)

        # All of these together
        mf, chif, vf, mf_err, chif_err, vf_err
            = fit.all(q, chiA, chiB, **kwargs)

    The arguments for each of these call methods are as follows:
    Arguments:
        q:      Mass ratio (q = mA/mB >= 1)

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
                the spins in the coorbital frame at the reference epoch. See
                THE PAPER for a definition of the coorbital frame.


    Optional arguments:

        omega0: Orbital frequency used to set the reference epoch.
                Default: None.

                If omega0 is None, the reference epoch is assumed to be at
                t=-100 M from the peak of the waveform, see THE PAPER for
                definition of the peak.

                If 'omega0' is given, the reference epoch is take to be the
                time at which the orbital frequency in the coprecessing frame
                equals omega0. omega0 should be in dimensionless units of
                rad/M, where M is the total mass.

                See THE PAPER for how the orbital frequency is
                computed as well as the definition of the coprecessing frame.

        allow_extrap:
            If False, raises a warning when q > 4.1 or |chiA|,|chiB| > 0.81,
                and raises an error when q > 6.1 or |chiA|,|chiB| > 1.
            If True, allows extrapolation to any q and |chiA|,|chiB| <= 1.
                Use at your own risk.
            Default: False.

    Optional PN evolution arguments:

        If the omega0 option is used, the spins need to be evolved from omega0
        until t=-100M, where the fits will be evaluated. For the late inspiral
        part, we use the internal spin evolution of NRSur7dq4 (also described
        in THE PAPER), which is very accurate. However, this surrogate is not
        long enough for small values of omega0 as it only has data starting at
        t=-4300M. Therefore, whenever the input omega0 is smaller than
        omega_switch_IG (defined below), we use PN evolution to go from omega0
        to about t=-4300M, beyond which we use NRSur7dq4 for spin evolution.

        PN_approximant:
            Approximant used to do the PN spin evolution. Choose from
            'SpinTaylorT4', 'SpinTaylorT1' or 'SpinTaylorT5'.
            Default: 'SpinTaylorT4'.

        PN_dt:
            Dimensionless time step size in units of M, used for the PN
            evolution. You may need to increase this if omega0 is very low.
            Default: 0.1

        PN_spin_order:
            Twice the PN order of spin effects. E.g., use 7 for 3.5PN.
            Default: 7

        PN_phase_order:
            Twice the PN order in phase. E.g., use 7 for 3.5PN.
            Default: 7

        t_sur_switch:
            The dimensionless time (from the peak) at which we switch from PN
            to the surrogate. Should be something larger than -4300.
            Default: -4000.

        omega_switch_IG:
            Initial guess for dimensionless orbital frequency, using which the
            switch will be made from PN to NRSur7dq4. This should be large
            enough to work for generic parts of the surrogate parameter space.
            You may need to increase this if the NRSur7dq4 model raises an
            exception like: "Got omega_ref=0.03 < 0.031=omega_0, too small!"
            Default: 0.03

            How t_sur_switch and omega_switch_IG work: The PN data is first
            generated starting at omega0, then the PN spins at omega_switch_IG
            are used to generate the NRSur7dq4 dynamics. NRSur7dq4 integrate
            the dynamics both forwards and backwards, so it will have omega and
            spins as a time series starting from -4300M. This is used to pick
            the omega0_sur and spins at t_sur_switch. Then the surrogate
            is reevaluated using omega0_sur and spins at t_sur_switch, thus
            ensuring that the switch always happens at t_sur_switch, even if
            omega_switch_IG corresponds to a later time.

    Inertial frame for returned values:

        The returned chif/vf are in the LAL inertial frame defined as follows:
            The +ve z-axis is along the orbital angular momentum at the
            reference epoch. The separation vector from the lighter BH to the
            heavier BH at the reference epoch is along the +ve x-axis. The
            y-axis completes the right-handed triad.

            Note that the default reference epoch corresponds to t=-100M, but
            if omega0 is given the reference epoch is taken to be the time at
            which the orbital frequency in the coprecessing frame is equal to
            omega0. This agrees with the LAL convention. See LIGO DCC document
            T1800226 for the LAL frame diagram.
    """

    #-------------------------------------------------------------------------
    def __init__(self, name, load_nrsur=False):

        # Param limits beyond which to raise a warning
        soft_param_lims = {
            'q': 4.1,
            'chiAmag': 0.81,
            'chiBmag': 0.81,
                }

        # Param limits beyond which to raise an error
        hard_param_lims = {
            'q': 6.1,
            'chiAmag': 1,
            'chiBmag': 1,
                }

        super(Fit7dq4, self).__init__(name, soft_param_lims, hard_param_lims)
        self.nrsur = None
        self.fitnode_time = -100      # Time at which the fits are constructed

    #-------------------------------------------------------------------------
    def _load_NRSur7dq4(self):
        import gwsurrogate
        from gwsurrogate.new.precessing_surrogate import splinterp_many
        self.nrsur = gwsurrogate.LoadSurrogate('NRSur7dq4')
        self.splinterp_many = splinterp_many

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

        extra_args = []
        extra_args.append({
            'omega0': 5e-3,
            'PN_approximant': 'SpinTaylorT4',
            'PN_dt': 0.1,
            'PN_spin_order': 7,
            'PN_phase_order': 7,
            })


        extra_args.append({
            'omega0': 6e-3,
            'PN_approximant': 'SpinTaylorT1',
            'PN_dt': 0.5,
            'PN_spin_order': 5,
            'PN_phase_order': 7,
            })

        extra_args.append({
            'omega0': 7e-3,
            'PN_approximant': 'SpinTaylorT5',
            'PN_dt': 1,
            'PN_spin_order': 7,
            'PN_phase_order': 5,
            })

        # These should be pure NRSur7dq4
        extra_args.append({'omega0': 3e-2})
        extra_args.append({'omega0': 5e-2})

        return extra_args

    #-------------------------------------------------------------------------
    def _get_fit_params(self, x, fit_key):
        """ Transforms the input parameter to fit parameters for the 7dq4 model.
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

    def _get_coorbital_frame_spins_at_idx(self, chiA, chiB, omega, lNhat, phi, \
            idx):
        """ Computes PN spins and dynamics at a given idx.

        Inputs:
            chiA:       Dimless spin evolution of BhA in inertial frame.
            chiB:       Dimless spin evolution of BhB in inertial frame.
            omega:      Orbital frequency evolution in dimless units.
            lNhat:      Orbital angular momentum direction evolution.
            phi:        Orbital phase evolution.
            idx:        Index for output.

        Outputs (all are time series):
            chiA_at_idx_coorb:   Spin of BhA at idx, in coorbital frame.
            chiB_at_idx_coorb:   Spin of BhB at idx, in coorbital frame.
            quat_copr_at_idx:   Coprecessing frame quaternion at idx.
            phi_at_idx:         Orbital phase in the coprecessing frame at idx.
            omega_at_idx        Dimensionless orbital frequency at idx.

        The inertial frame is assumed to be aligned to the coorbital frame at
        the first index.
        """

        # Compute omega, inertial spins, angular momentum direction and orbital
        # phase at idx
        omega_at_idx = omega[idx]
        chiA_at_idx = chiA[idx]
        chiB_at_idx = chiB[idx]
        lNhat_at_idx = lNhat[idx]
        phi_at_idx = phi[idx]

        # Align the z-direction along orbital angular momentum direction
        # at idx. This moves us in to the coprecessing frame.
        quat_copr_at_idx = utils.alignVec_quat(lNhat_at_idx)
        chiA_at_idx_copr = utils.transformTimeDependentVector(
                np.array([quat_copr_at_idx]).T,
                np.array([chiA_at_idx]).T, inverse=1).T[0]
        chiB_at_idx_copr = utils.transformTimeDependentVector(
                np.array([quat_copr_at_idx]).T,
                np.array([chiB_at_idx]).T, inverse=1).T[0]

        # get coorbital frame spins at idx
        chiA_at_idx_coorb = utils.rotate_in_plane(chiA_at_idx_copr, phi_at_idx)
        chiB_at_idx_coorb = utils.rotate_in_plane(chiB_at_idx_copr, phi_at_idx)

        return chiA_at_idx_coorb, chiB_at_idx_coorb, quat_copr_at_idx, \
            phi_at_idx, omega_at_idx

    def _get_PN_spins_at_surrogate_start(self, PN_approximant, q, omega0, \
            chiA0, chiB0, PN_dt, PN_spin_order, PN_phase_order, \
            omega_switch_IG, t_sur_switch):
        """ Computes PN spins and frame dynamics at a time close to the start
            of the surrogate waveform model.

            Generates PN spins and frame quantities using spins at omega0.
            Then uses the PN spins at omega_switch_IG to generate the surrogate
            dynamics.
            Then gets the surrogate orbital frequency at t_sur_switch, let's
            call this omega_init_sur.
            Then use the PN spins at omega_init_sur to regenerate the surrogate
            dynamics.
        """

        # Get PN spin evolution starting at omega0
        omega_PN, phi_PN, chiA_PN, chiB_PN, lNhat_PN, e1_PN \
                = lal_spin_evloution_wrapper(PN_approximant, q, omega0, \
                chiA0, chiB0, PN_dt, PN_spin_order, PN_phase_order)

        # Get PN coorbital frame spins and frame dynamics at
        # omega_PN=omega_switch_IG
        idx = np.argmin(np.abs(omega_PN - omega_switch_IG))
        chiA_PN_at_idx_coorb, chiB_PN_at_idx_coorb, quat_PN_copr_at_idx, \
            phi_PN_at_idx, omega_PN_at_idx \
            = self._get_coorbital_frame_spins_at_idx(chiA_PN, chiB_PN, \
            omega_PN, lNhat_PN, phi_PN, idx)

        # Now evaluate the surrogate dynamics (both forwards and backwards)
        # using PN spins at omega_switch_IG
        quat_sur, orbphase_sur, chiA_copr_sur, chiB_copr_sur \
            = self.nrsur._sur_dimless.get_dynamics(q, chiA_PN_at_idx_coorb, \
            chiB_PN_at_idx_coorb, init_quat=quat_PN_copr_at_idx,
            init_orbphase=phi_PN_at_idx, omega_ref=omega_switch_IG)
        dyn_times = self.nrsur._sur_dimless.tds
        omega_sur = np.gradient(orbphase_sur, dyn_times)

        # Get surrogate orbital frequency at t_sur_switch, which is
        # close to the start of the surrogate data
        omega_init_sur = omega_sur[np.argmin(np.abs( \
                dyn_times - t_sur_switch))]

        # Get PN coorbital frame spins and frame dynamics at omega_init_sur
        idx = np.argmin(np.abs(omega_PN - omega_init_sur))
        chiA_PN_at_idx_coorb, chiB_PN_at_idx_coorb, quat_PN_copr_at_idx, \
            phi_PN_at_idx, omega_PN_at_idx \
            = self._get_coorbital_frame_spins_at_idx(chiA_PN, chiB_PN, \
            omega_PN, lNhat_PN, phi_PN, idx)

        return chiA_PN_at_idx_coorb, chiB_PN_at_idx_coorb, \
            quat_PN_copr_at_idx, phi_PN_at_idx, omega_PN_at_idx, \
            chiA_PN, chiB_PN, omega_PN


    #-------------------------------------------------------------------------
    def _evolve_spins(self, q, chiA0, chiB0, omega0, \
            return_spin_evolution=False, **kwargs):
        """ Evolves spins of the component BHs from an initial orbital
        frequency = omega0 until t=-100 M from the peak of the waveform.  If
        omega0 < omega_switch_IG, use PN to evolve the spins until
        t=t_sur_switch. Then evolves further with the NRSur7dq4 waveform model
        until t=-100M from the peak.

        Returns spins in the coorbital frame at t=-100M, as well as the
        coprecessing frame quaternion and orbital phase in the coprecessing
        frame at this time.

        If return_spin_evolution is given, also returns the PN and surrogate
        spin times series.
        """

        PN_approximant = kwargs.pop('PN_approximant', 'SpinTaylorT4')
        PN_dt = kwargs.pop('PN_dt', 0.1)
        PN_spin_order = kwargs.pop('PN_spin_order', 7)
        PN_phase_order = kwargs.pop('PN_phase_order', 7)
        # Initial guess for surrogate omega0, this should be large enough for
        # all q=6 cases
        omega_switch_IG = kwargs.pop('omega_switch_IG', 0.03)
        # The surrogate begins at -4300, use -4000 to be safe
        t_sur_switch = kwargs.pop('t_sur_switch', -4000)
        self._check_unused_kwargs(kwargs)

        # Load NRSur7dq4 if not previously loaded
        if self.nrsur is None:
            self._load_NRSur7dq4()

        # If omega0 is below the NRSur7dq4 initial guess frequency, we use PN
        # to evolve the spins. We get the initial spins and omega_init_sur such
        # that should go into the surrogate such that the inital time is
        # t_sur_switch.
        if omega0 < omega_switch_IG:
            chiA0_nrsur_coorb, chiB0_nrsur_coorb, quat0_nrsur_copr, \
                phi0_nrsur, omega_init_sur, chiA_PN, chiB_PN, omega_PN \
                = self._get_PN_spins_at_surrogate_start(PN_approximant, q, \
                omega0, chiA0, chiB0, PN_dt, PN_spin_order, PN_phase_order, \
                omega_switch_IG, t_sur_switch)

        # If omega0 >= omega_switch_IG, we evolve spins directly with NRSur7dq4
        # waveform model. We set the coprecessing frame quaternion to identity
        # and orbital phase to 0 at omega=omega0, hence the coprecessing frame
        # is the same as the inertial frame here.
        else:
            # Note that here we set omega_init_sur to omega0
            chiA0_nrsur_coorb, chiB0_nrsur_coorb, quat0_nrsur_copr, \
                phi0_nrsur, omega_init_sur, chiA_PN, chiB_PN, omega_PN \
                = chiA0, chiB0, [1,0,0,0], 0, omega0, None, None, None

        # Now evaluate the surrogate dynamics using PN spins at omega_init_sur
        quat_sur, orbphase_sur, chiA_copr_sur, chiB_copr_sur \
            = self.nrsur._sur_dimless.get_dynamics(q, chiA0_nrsur_coorb, \
            chiB0_nrsur_coorb, init_quat=quat0_nrsur_copr,
            init_orbphase=phi0_nrsur, omega_ref=omega_init_sur)

        # get data at time node where remnant fits are done
        dyn_times = self.nrsur._sur_dimless.tds
        nodeIdx = np.argmin(np.abs(dyn_times - self.fitnode_time))
        quat_fitnode = quat_sur.T[nodeIdx]
        orbphase_fitnode = orbphase_sur[nodeIdx]

        # get coorbital frame spins at the time node
        chiA_coorb_fitnode = utils.rotate_in_plane(chiA_copr_sur[nodeIdx],
                orbphase_fitnode)
        chiB_coorb_fitnode = utils.rotate_in_plane(chiB_copr_sur[nodeIdx],
                orbphase_fitnode)

        if return_spin_evolution:
            # Transform spins to the reference inertial frame
            chiA_inertial_sur = utils.transformTimeDependentVector(quat_sur, \
                    chiA_copr_sur.T).T
            chiB_inertial_sur = utils.transformTimeDependentVector(quat_sur, \
                    chiB_copr_sur.T).T
            spin_evolution = {
                    't_sur': dyn_times,
                    'chiA_sur': chiA_inertial_sur,
                    'chiB_sur': chiB_inertial_sur,
                    'orbphase_sur': orbphase_sur,
                    'quat_sur': quat_sur,
                    'omega_PN': omega_PN,
                    'chiA_PN': chiA_PN,
                    'chiB_PN': chiB_PN,
                    'omega_init_sur': omega_init_sur,
                }
        else:
            spin_evolution = None

        return chiA_coorb_fitnode, chiB_coorb_fitnode, quat_fitnode, \
            orbphase_fitnode, spin_evolution

    #-------------------------------------------------------------------------
    def _eval_wrapper(self, fit_key, q, chiA, chiB, **kwargs):
        """Evaluates the NRSur7dq4Remnant model.
        """
        chiA = np.array(chiA)
        chiB = np.array(chiB)

        # Warn/Exit if extrapolating
        allow_extrap = kwargs.pop('allow_extrap', False)
        omega0 = kwargs.pop('omega0', None)
        self._check_param_limits(q, chiA, chiB, allow_extrap)

        if omega0 is None:
            # If omega0 is given, assume chiA, chiB are the coorbital frame
            # spins at t=-100 M.
            x = np.concatenate(([q], chiA, chiB))
        else:
            # If omega0 is given, evolve the spins from omega0
            # to t = -100 M from the peak.
            chiA_coorb_fitnode, chiB_coorb_fitnode, quat_fitnode, \
                orbphase_fitnode, _ \
                = self._evolve_spins(q, chiA, chiB, omega0, **kwargs)

            # x should contain coorbital frame spins at t=-100M
            x = np.concatenate(([q], chiA_coorb_fitnode, chiB_coorb_fitnode))

        def eval_vector_fit(x, fit_key):
            res = self._evaluate_fits(x, fit_key)
            fit_val = res.T[0]
            fit_err = res.T[1]
            if omega0 is not None:
                # The fit are constructed in the coorbital frame at t=-100M,
                # now we transform the remnant vectors into the LAL inertial
                # frame, which is the same as the coorbital frame at omega0.
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
