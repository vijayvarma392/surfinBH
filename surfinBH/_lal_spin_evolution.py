import numpy as np
import lalsimulation as lalsim
from lal import MSUN_SI, MTSUN_SI, PC_SI, C_SI
from . import _utils

def lal_spin_evloution_wrapper(approximant, q, omega0, chiA0, chiB0,
        dt, spinO, phaseO):
    """
    Inputs:
        approximant:    'SpinTaylorT1/T2/T4'
        q:              Mass ratio (q>=1)
        omega0:         Initial orbital frequency in dimless units.
        chiA0:          Dimless spin of BhA at initial freq.
        chiB0:          Dimless spin of BhB at initial freq.
        dt:             Dimless step time for evolution.
        spinO:          Twice PN order of spin effects.
        phaseO:         Twice PN order in phase.

    Outputs (all are time series):
        Omega:    Dimensionless orbital frequency.
        Phi:            Orbital phase (radians)
        ChiA:           Dimensionless spin of BhA
        ChiB:           Dimensionless spin of BhB
        LNhat:          Orbital angular momentum direction
        E1:             Orbital plane basis vector

    The frame is defined at the initial frequency, as follows:
        z-axis is set by the orbital angular momentum direction.
        x-axis is the separation vector from BhB to BhA.
        y-axis completes the triad by right-hand rule.
        All quantities are defined in this fixed frame, including initial spins,
        returned spins, other vectors like LNhat, etc.
    """

    approxTag = lalsim.SimInspiralGetApproximantFromString(approximant)

    # Total mass in solar masses
    M = 100     # This does not affect the returned values as they are
                # dimension less

    # time step and initial GW freq in SI units
    MT = M*MTSUN_SI
    deltaT = dt*MT
    fStart = omega0/np.pi/MT

    # component masses of the binary
    m1_SI =  M*MSUN_SI*q/(1.+q)
    m2_SI =  M*MSUN_SI/(1.+q)

    # spins at fStart
    s1x, s1y, s1z = chiA0
    s2x, s2y, s2z = chiB0

    # integrate as far forward as possible
    fEnd = 0

    # initial value of orbital angular momentum unit vector, i.e at fStart
    lnhatx, lnhaty, lnhatz = 0,0,1

    # initial value of orbital plane basis vector, i.e at fStart
    e1x, e1y, e1z = 1, 0, 0

    # tidal deformability parameters
    lambda1, lambda2 = 0, 0
    quadparam1, quadparam2 = 1, 1

    # twice PN order of tidal effects
    tideO = 0

    # include some L-S terms
    lscorr = 1

    ### This function evolves the orbital equations for a precessing binary
    ### using the "TaylorT1/T2/T4" approximant for solving the orbital dynamics
    ### (see arXiv:0907.0700 for a review of the various PN approximants).
    ###
    ### It returns time series of the "orbital velocity", orbital phase,
    ### and components for both individual spin vectors, the "Newtonian"
    ### orbital angular momentum (which defines the instantaneous plane)
    ### and "E1", a basis vector in the instantaneous orbital plane. Note that
    ### LNhat and E1 completely specify the instantaneous orbital plane.
    ### It also returns the time and phase of the final time step
    ###
    ### For input, the function takes the two masses, the initial orbital phase,
    ### Values of S1, S2, LNhat, E1 vectors at starting time,
    ### the desired time step size, the starting GW frequency,
    ### and PN order at which to evolve the phase,
    ###
    ### NOTE: All vectors are given in the frame
    ### where the z-axis is set by the angular momentum at reference frequency,
    ### the x-axis is chosen orthogonal to it, and the y-axis is given by the
    ### RH rule. Initial values must be passed in this frame, and the time
    ### series of the vector components will also be returned in this frame.
    ###
    ###
    ### V,            post-Newtonian parameter [returned]
    ### Phi,          orbital phase            [returned]
    ### S1x,          Spin1 vector x component [returned]
    ### S1y,          "    "    "  y component [returned]
    ### S1z,          "    "    "  z component [returned]
    ### S2x,          Spin2 vector x component [returned]
    ### S2y,          "    "    "  y component [returned]
    ### S2z,          "    "    "  z component [returned]
    ### LNhatx,       unit orbital ang. mom. x [returned]
    ### LNhaty,       "    "    "  y component [returned]
    ### LNhatz,       "    "    "  z component [returned]
    ### E1x,          orb. plane basis vector x[returned]
    ### E1y,          "    "    "  y component [returned]
    ### E1z,          "    "    "  z component [returned]
    ###
    ###
    ### deltaT,       sampling interval (s)
    ### m1_SI,        mass of companion 1 (kg)
    ### m2_SI,        mass of companion 2 (kg)
    ### fStart,       starting GW frequency
    ### fEnd,         ending GW frequency, fEnd=0 means integrate as far
    ###                 forward as possible
    ### s1x,          initial value of S1x
    ### s1y,          initial value of S1y
    ### s1z,          initial value of S1z
    ### s2x,          initial value of S2x
    ### s2y,          initial value of S2y
    ### s2z,          initial value of S2z
    ### lnhatx,       initial value of LNhatx
    ### lnhaty,       initial value of LNhaty
    ### lnhatz,       initial value of LNhatz
    ### e1x,          initial value of E1x
    ### e1y,          initial value of E1y
    ### e1z,          initial value of E1z
    ### lambda1,      tidal deformability of mass 1
    ### lambda2,      tidal deformability of mass 2
    ### quadparam1,   phenom. parameter describing induced quad.
    ###                 moment of body 1 (=1 for BHs, ~2-12 for NSs)
    ### quadparam2,   phenom. parameter describing induced quad.
    ###                 moment of body 2 (=1 for BHs, ~2-12 for NSs)
    ### spinO,        twice PN order of spin effects
    ### tideO,        twice PN order of tidal effects
    ### phaseO,       twice post-Newtonian order
    ### lscorr,       Flag to include L-S correction terms
    ### approx        PN approximant (SpinTaylorT1/T2/T4)
    V, Phi, S1x, S1y, S1z, S2x, S2y, S2z, LNhatx, LNhaty, LNhatz, \
        E1x, E1y, E1z = lalsim.SimInspiralSpinTaylorPNEvolveOrbit(deltaT, \
        m1_SI, m2_SI, fStart, fEnd, s1x, s1y, s1z, s2x, s2y, s2z, \
        lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2, \
        quadparam1, quadparam2, spinO, tideO, phaseO, lscorr, approxTag)

    V = np.array(V.data.data)
    Phi = np.array(Phi.data.data)
    ChiA = np.array([S1x.data.data, S1y.data.data, S1z.data.data]).T
    ChiB = np.array([S2x.data.data, S2y.data.data, S2z.data.data]).T
    LNhat = np.array([LNhatx.data.data, LNhaty.data.data, LNhatz.data.data]).T
    E1 = np.array([E1x.data.data, E1y.data.data, E1z.data.data]).T

    Omega = V**3
    return Omega, Phi, ChiA, ChiB, LNhat, E1

# -------------------------------------------------------------------------
def evolve_pn_spins(q, chiA0, chiB0, omega0, omegaTimesM_final,
        approximant='SpinTaylorT4', dt=0.1, spinO=7, phaseO=7):
    """ Evolves PN spins from a starting orbital frequency and spins to a final
    frequency.

    Inputs:
        q:                  Mass ratio (q>=1)
        chiA0:              Dimless spin of BhA at initial freq.
        chiB0:              Dimless spin of BhB at initial freq.
        omega0:       Initial orbital frequency in dimless units.
        omegaTimesM_final:  Final orbital frequency in dimless units.
        approximant:        'SpinTaylorT1/T2/T4'. Default: 'SpinTaylorT4'.
        dt:        Dimless step time for evolution. Default: 0.1 .
        spinO:              Twice PN order of spin effects. Default: 5 .
        phaseO:             Twice PN order in phase. Default: 8 .

    Outputs (all are time series):
        chiA_end_copr:      Spin of BhA at final frequency, in coprecessing
                                frame.
        chiB_end_copr:      Spin of BhB at final frequency, in coprecessing
                                frame.
        q_copr_end:         Coprecessing frame quaternion at final frequency.
        phi_end:            Orbital phase in the coprecessing frame at final
                                frequency.
        omegaTimesM_end     Dimensionless final frequency. Should agree with
                                omegaTimesM_final.

    The inertial frame is assumed to be aligned
    to the coorbital frame at orbital frequency = omega0. chiA0 and chiB0
    are the inertial/coorbital frame spins at omega0.
    """
    omega, phi, chiA, chiB, lNhat, e1 = lal_spin_evloution_wrapper(approximant,
            q, omega0, chiA0, chiB0, dt, spinO, phaseO)

    # Compute omega, inertial spins, angular momentum direction and orbital
    # phase when omega = omegaTimesM_final
    end_idx = np.argmin(np.abs(omega - omegaTimesM_final))
    omegaTimesM_end = omega[end_idx]
    chiA_end = chiA[end_idx]
    chiB_end = chiB[end_idx]
    lNhat_end = lNhat[end_idx]
    phi_end = phi[end_idx]

    # Align the z-direction along orbital angular momentum direction
    # at end_idx. This moves us in to the coprecessing frame.
    q_copr_end = _utils.alignVec_quat(lNhat_end)
    chiA_end_copr = _utils.transformTimeDependentVector(
            np.array([q_copr_end]).T,
            np.array([chiA_end]).T, inverse=1).T[0]
    chiB_end_copr = _utils.transformTimeDependentVector(
            np.array([q_copr_end]).T,
            np.array([chiB_end]).T, inverse=1).T[0]

    return chiA_end_copr, chiB_end_copr, q_copr_end, phi_end, omegaTimesM_end
