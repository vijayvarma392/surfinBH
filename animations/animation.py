desc = """Animations of binary black hole scattering.
Generates an animation of a binary black hole merger and the final remnant.

Example usage:
python animation.py --q 2 --omega0 1e-2 --chiA0 0.6 0.5 -0.1 --chiB0 0.7 0.3 0.1

Note: Time values displayed in the plot are non-uniform and non-linear:
During the inspiral there are 20 frames per orbit.
After the merger each frame corresponds to a time step of 100M.
"""

import numpy as np
import matplotlib.pyplot as P
import argparse

from scipy.interpolate import UnivariateSpline
from scipy.interpolate import InterpolatedUnivariateSpline

import surfinBH
import NRSur7dq2
#from gwtools import rotations

from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d import proj3d
import matplotlib.animation as animation
from matplotlib.patches import FancyArrowPatch


#----------------------------------------------------------------------------
def spline_interp(newX, oldX, oldY, allowExtrapolation=False):
    """ Interpolates using splnes.
        If allowExtrapolation=True, extrapolates to zero.
    """
    if len(oldY) != len(oldX):
        raise Exception('Lengths dont match.')

    if not allowExtrapolation:
        if np.min(newX) < np.min(oldX) or np.max(newX) > np.max(oldX):
            print np.min(newX), np.min(oldX), np.max(newX), np.max(oldX)
            print np.min(newX) < np.min(oldX)
            print np.max(newX) > np.max(oldX)
            raise Exception('Trying to extrapolate, but '\
                'allowExtrapolation=False')

    if not np.all(np.diff(oldX) > 0):
        raise Exception('oldX must have increasing values')

    # returns 0 when extrapolating
    newY = InterpolatedUnivariateSpline(oldX, oldY, ext=1)(newX)
    return newY

#----------------------------------------------------------------------------
def get_trajectory(separation, quat_nrsur, orbphase_nrsur, bh_label):
    """ Gets trajectory of a component BH in a binary given the separation,
    the coprecessing frame quaternion and orbital phase in the coprecessing
    frame.
    """
    #FIXME COM should not be equidistant from both BHs

    if bh_label == 'A':
        offset = 0
    else:
        offset = np.pi

    x_copr = separation * np.cos(orbphase_nrsur+offset)
    y_copr = separation * np.sin(orbphase_nrsur+offset)
    z_copr = np.zeros(len(x_copr))

    Bh_traj_copr = np.array([x_copr, y_copr, z_copr])
    Bh_traj = surfinBH._utils.transformTimeDependentVector(quat_nrsur, \
        Bh_traj_copr, inverse=0)
    return Bh_traj


#-----------------------------------------------------------------------------
def get_uniform_in_orbits_times(t, phi_orb, pts_per_orbit):
    """
    returns sparse time array such that there are pts_per_orbit points
    in each orbit.
    """
    # get numer of orbits
    n_orbits = int(abs((phi_orb[-1] - phi_orb[0])/(2*np.pi)))

    # get sparse times such that there are pts_per_orbit points in each orbit
    n_pts = int(n_orbits*pts_per_orbit)
    phi_orb_sparse = np.linspace(phi_orb[0], phi_orb[-1], n_pts)
    t_sparse = np.interp(phi_orb_sparse, phi_orb, t)

    return t_sparse

#----------------------------------------------------------------------------
def get_omegaOrb_from_sparse_data(t_sparse, phiOrb_sparse):
    """ Computes orbital frequency from sparse data using splines.
    """
    # spline interpolant for phase
    phiOrb_spl = UnivariateSpline(t_sparse, phiOrb_sparse, s=0)

    # spline for phase derivative
    omegaOrb_spl = phiOrb_spl.derivative()

    return omegaOrb_spl(t_sparse)

#----------------------------------------------------------------------------
def get_separation_from_omega(omega):
    """ Let's do zeroth order approx: Keppler's law """
    separation = omega**(-2./3)
    return separation


#----------------------------------------------------------------------------
def get_quivers(Bh_loc, chi_vec, scale_factor=10):
    """ Gets quivers for spins on a BH
    """
    X, Y, Z =  Bh_loc
    u, v, w =  chi_vec
    segments = (X, Y, Z, X+v*scale_factor, Y+u*scale_factor, Z+w*scale_factor)
    segments = np.array(segments).reshape(6,-1)
    return [[[x, y, z], [u, v, w]] for x, y, z, u, v, w in zip(*list(segments))]

#----------------------------------------------------------------------------
def update_lines(num, lines, hist_frames):
    """ The function that goes into animation
    """
    current_time = t[num]
    if current_time < 0:

        if num == 0:
            # Clear remnant stuff
            line = lines[6]
            line.set_data([], [])
            line.set_3d_properties([])
            line = lines[7]
            line.set_segments([])

        for idx in range(len(dataLines_binary)):
            time_text.set_text('time = %.1f M'%current_time)

            line = lines[idx]
            data = dataLines_binary[idx]


            if idx < 4:
                if idx < 2:
                    start = max(0, num-hist_frames)
                else:
                    start = max(0, num-1)

                # NOTE: there is no .set_data() for 3 dim data...
                line.set_data(data[0:2, start:num])
                line.set_3d_properties(data[2, start:num])
            else:
                if idx == 4:
                    Bh_loc = BhA_traj[:,num-1]
                    chi_vec = chiA_nrsur[num-1]
                elif idx == 5:
                    Bh_loc = BhB_traj[:,num-1]
                    chi_vec = chiB_nrsur[num-1]

                line.set_segments(get_quivers(Bh_loc, chi_vec))
    else:
        num = num - zero_idx
        if num == 0:
            # Clear binary stuff
            for idx in range(4):
                line = lines[idx]
                line.set_data([], [])
                line.set_3d_properties([])
            for idx in range(4,6):
                line = lines[idx]
                line.set_segments([])

        for idx in range(len(dataLines_remnant)):
            time_text.set_text('time = %.1f M'%current_time)
            line = lines[6+idx]
            data = dataLines_remnant[idx]

            if idx == 0:
                # NOTE: there is no .set_data() for 3 dim data...
                line.set_data(data[0:2, num-1:num])
                line.set_3d_properties(data[2, num-1:num])
            else:
                Bh_loc = BhC_traj[:,num-1]
                chi_vec = chif
                line.set_segments(get_quivers(Bh_loc, chi_vec))

    return lines

#############################    main    ##################################

parser = argparse.ArgumentParser(description=desc,
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--omega0', type=float, required=True,
    help='Starting orbital frequency. Currently, > 0.018.')
parser.add_argument('--q', type=float, required=True,
    help='Mass ratio.')
parser.add_argument('--chiA0', type=float, required=True, nargs=3,
    help='Spin of BhA at omega0. Array of size 3.')
parser.add_argument('--chiB0', type=float, required=True, nargs=3,
    help='Spin of BhB at omega0. Array of size 3.')
parser.add_argument('--fit_name', type=str, default='surfinBH7dq2', \
    help='Fit to use.')

args = parser.parse_args()

q = args.q
chiA0 = np.array(args.chiA0)
chiB0 = np.array(args.chiB0)
omega0 =  args.omega0

# evaluate remnant fit
fit_name = args.fit_name
fit = surfinBH.LoadFits(fit_name)
mf, chif, vf, mf_err, chif_err, vf_err = fit.all(q, chiA0, chiB0, omega0=omega0)

mA = q/(1.+q)
mB = 1./(1.+q)

nr_sur = NRSur7dq2.NRSurrogate7dq2()

# get NRSur dynamics
quat_nrsur, orbphase_nrsur, _, _ \
    = nr_sur.get_dynamics(q, chiA0, chiB0, omega_ref=omega0, \
    allow_extrapolation=True)

pts_per_orbit = 20
t_binary = get_uniform_in_orbits_times(nr_sur.tds, orbphase_nrsur, \
    pts_per_orbit)

# interpolate dynamics on to t_binary
quat_nrsur = np.array([spline_interp(t_binary, nr_sur.tds, tmp) \
    for tmp in quat_nrsur])
orbphase_nrsur = spline_interp(t_binary, nr_sur.tds, orbphase_nrsur)

omega_nrsur = get_omegaOrb_from_sparse_data(t_binary, orbphase_nrsur)

h_nrsur, chiA_nrsur, chiB_nrsur = nr_sur(q, chiA0, chiB0, \
    f_ref=omega0/np.pi, return_spins=True, allow_extrapolation=True, LMax=2,
    t=t_binary)

#LHat = rotations.lHat_from_quat(quat_nrsur)

separation = get_separation_from_omega(omega_nrsur)
BhA_traj = get_trajectory(separation, quat_nrsur, orbphase_nrsur, 'A')
BhB_traj = get_trajectory(separation, quat_nrsur, orbphase_nrsur, 'B')


# time array for remnant
t_remnant = np.arange(0, 10000, 100)

# assume merger is at origin
BhC_traj = np.array([tmp*t_remnant for tmp in vf])

# Attaching 3D axis to the figure
fig = P.figure()
ax = axes3d.Axes3D(fig)

# FIXME check that this makes sense
markersize_BhA = mA*50
markersize_BhB = mB*50
markersize_BhC = mf*50

time_text = ax.text2D(0.05, 0.9, '', transform=ax.transAxes, fontsize=14)

# NOTE: Can't pass empty arrays into 3d version of plot()
dataLines_binary = [BhA_traj, BhB_traj, BhA_traj, BhB_traj, 1, 1]

lines = [\
    # These two are for plotting component tracjectories
    ax.plot(BhA_traj[0,0:1], BhA_traj[1,0:1], BhA_traj[2,0:1])[0], \
    ax.plot(BhB_traj[0,0:1], BhB_traj[1,0:1], BhB_traj[2,0:1])[0], \

    # These two are for plotting component BHs
    ax.plot(BhA_traj[0,0:1], BhA_traj[1,0:1], BhA_traj[2,0:1], \
        marker='o', markersize=markersize_BhA, markerfacecolor='k', \
        markeredgewidth=0)[0], \
    ax.plot(BhB_traj[0,0:1], BhB_traj[1,0:1], BhB_traj[2,0:1], \
        marker='o', markersize=markersize_BhB, markerfacecolor='k',
        markeredgewidth=0)[0], \

    # These two are plotting component BH spins
    ax.quiver(0,0,0,1,1,1),
    ax.quiver(0,0,0,1,1,1),

    # This is for plotting remnant BH
    ax.plot(BhC_traj[0,0:1]-1e10, BhC_traj[1,0:1], BhC_traj[2,0:1], \
        marker='o', markersize=markersize_BhC, markerfacecolor='k', \
        markeredgewidth=0)[0], \
    # This is for plotting remnant spin
    ax.quiver(-1e10,0,0,-1e10,-1e10,-1e10),
    ]

dataLines_remnant = [BhC_traj, 1]

max_range = np.nanmax(separation)

# Setting the axes properties
ax.set_xlim3d([-max_range, max_range])
ax.set_xlabel('X')

ax.set_ylim3d([-max_range, max_range])
ax.set_ylabel('Y')

ax.set_zlim3d([-max_range, max_range])
ax.set_zlabel('Z')

ax.set_title(fit_name, fontsize=16)

# Creating the Animation object
hist_frames = 15


zero_idx = np.argmin(np.abs(t_binary))

# common time array
t = np.append(t_binary[:zero_idx], t_remnant)

line_ani = animation.FuncAnimation(fig, update_lines, len(t), \
        fargs=(lines, hist_frames), interval=50, blit=False, repeat=True,
        repeat_delay=5e3)

# Pause settings
pause = False
def onClick(event):
    global pause
    if pause:
        line_ani.event_source.start()
        pause = False
    else:
        line_ani.event_source.stop()
        pause = True
fig.canvas.mpl_connect('button_press_event', onClick)
P.show()
