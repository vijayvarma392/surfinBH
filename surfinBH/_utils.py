import numpy as np

#-----------------------------------------------------------------------------
def multiplyQuats(q1, q2):
    """q1, q2 must be [scalar, x, y, z] but those may be arrays or scalars"""
    return np.array([
            q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3],
            q1[2]*q2[3] - q2[2]*q1[3] + q1[0]*q2[1] + q2[0]*q1[1],
            q1[3]*q2[1] - q2[3]*q1[1] + q1[0]*q2[2] + q2[0]*q1[2],
            q1[1]*q2[2] - q2[1]*q1[2] + q1[0]*q2[3] + q2[0]*q1[3]])

#-----------------------------------------------------------------------------
def quatInv(q):
    """Returns QBar such that Q*QBar = 1"""
    qConj = -q
    qConj[0] = -qConj[0]
    normSqr = multiplyQuats(q, qConj)[0]
    return qConj/normSqr

#-----------------------------------------------------------------------------
def alignVec_quat(vec):
    """Returns a unit quaternion that will align vec with the z-axis"""
    alpha = np.arctan2(vec[1], vec[0])
    beta = np.arccos(vec[2])
    gamma = -alpha*vec[2]
    cb = np.cos(0.5*beta)
    sb = np.sin(0.5*beta)
    return np.array([cb*np.cos(0.5*(alpha + gamma)),
                     sb*np.sin(0.5*(gamma - alpha)),
                     sb*np.cos(0.5*(gamma - alpha)),
                     cb*np.sin(0.5*(alpha + gamma))])

#-----------------------------------------------------------------------------
def lHat_from_quat(quat):
    qInv = quatInv(quat)
    return multiplyQuats(quat, multiplyQuats(
                    np.array([0., 0., 0., 1.]), qInv))[1:]

#-----------------------------------------------------------------------------
def transformTimeDependentVector(quat, vec, inverse=0):
    """Given (for example) a minimal rotation frame quat, transforms
    vec from the minimal rotation frame to the inertial frame.
    With inverse=1, transforms from the inertial frame to the minimal
    rotation frame."""
    qInv = quatInv(quat)
    if inverse:
        return transformTimeDependentVector(qInv, vec, inverse=0)

    return multiplyQuats(quat, multiplyQuats(np.append(np.array([
            np.zeros(len(vec[0]))]), vec, 0), qInv))[1:]


#-----------------------------------------------------------------------------
def rotate_in_plane(chi, phase):
    """For transforming spins between the coprecessing and coorbital frames"""
    v = chi.T
    sp = np.sin(phase)
    cp = np.cos(phase)
    res = 1.*v
    res[0] = v[0]*cp + v[1]*sp
    res[1] = v[1]*cp - v[0]*sp
    return res.T

#-----------------------------------------------------------------------------
def transform_vector_coorb_to_inertial(vec_coorb, orbPhase, quat_copr):
    """Given a vector (of size 3) in coorbital frame, orbital phase in
    coprecessing frame and a minimal rotation frame quat, transforms
    the vector from the coorbital to the inertial frame.
    """

    # Transform to coprecessing frame
    vec_copr = rotate_in_plane(vec_coorb, -orbPhase)

    # Transform to inertial frame
    vec = transformTimeDependentVector(np.array([quat_copr]).T,
        np.array([vec_copr]).T).T[0]

    return np.array(vec)


def transform_error_coorb_to_inertial(vec_coorb, vec_err_coorb, orbPhase,
        quat_copr):
    """ Transform error in a vector from the coorbital frame to the inertial
    frame. Generates distributions in the coorbital frame, transforms them
    to inertial frame and returns 1-simga widths in the inertial frame.
    """

    # for reproducibility
    np.random.seed(0)

    # Get distribution in coorbital frame
    dist_coorb = np.array([np.random.normal(m, s, 1000)
        for m,s in zip(vec_coorb, vec_err_coorb)]).T

    # Transform distribution to coprecessing frame
    dist_copr = rotate_in_plane(dist_coorb, -orbPhase)

    # Transform distribution to inertial frame
    dist_inertial = transformTimeDependentVector(
        np.array([quat_copr for _ in dist_copr]).T, dist_copr.T).T

    # Get 1sigma width in inertial frame
    vec_err_inertial = np.std(dist_inertial, axis=0)

    return vec_err_inertial
