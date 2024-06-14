'''
module for calculating RMSD using the Kabsch method, QCP method,
or calculating without alignment 
'''

import numpy as np
from scipy import optimize
from typing import Tuple

from CanonizedRMSD.data import RMSDResult

def center_of_geometry(coordinates: np.ndarray) -> np.ndarray:
    # Center of geometry for the given coordinates.
    assert coordinates.shape[1] == 3
    return np.mean(coordinates, axis=0)

def center(coordinates: np.ndarray) -> np.ndarray:
    # Coordinates after recentered.
    return coordinates - center_of_geometry(coordinates)

def M_mtx(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    # Compute inner product between coordinate matrices.
    return B.T @ A

def K_mtx(M):
    # Compute symmetric key matrix.
    assert M.shape == (3, 3)
    S_xx = M[0, 0]
    S_xy = M[0, 1]
    S_xz = M[0, 2]
    S_yx = M[1, 0]
    S_yy = M[1, 1]
    S_yz = M[1, 2]
    S_zx = M[2, 0]
    S_zy = M[2, 1]
    S_zz = M[2, 2]
    # p = plus, m = minus
    S_xx_yy_zz_ppp = S_xx + S_yy + S_zz
    S_yz_zy_pm = S_yz - S_zy
    S_zx_xz_pm = S_zx - S_xz
    S_xy_yx_pm = S_xy - S_yx
    S_xx_yy_zz_pmm = S_xx - S_yy - S_zz
    S_xy_yx_pp = S_xy + S_yx
    S_zx_xz_pp = S_zx + S_xz
    S_xx_yy_zz_mpm = -S_xx + S_yy - S_zz
    S_yz_zy_pp = S_yz + S_zy
    S_xx_yy_zz_mmp = -S_xx - S_yy + S_zz
    return np.array(
        [
            [S_xx_yy_zz_ppp, S_yz_zy_pm, S_zx_xz_pm, S_xy_yx_pm],
            [S_yz_zy_pm, S_xx_yy_zz_pmm, S_xy_yx_pp, S_zx_xz_pp],
            [S_zx_xz_pm, S_xy_yx_pp, S_xx_yy_zz_mpm, S_yz_zy_pp],
            [S_xy_yx_pm, S_zx_xz_pp, S_yz_zy_pp, S_xx_yy_zz_mmp],
        ]
    )


def coefficients(M: np.ndarray, K: np.ndarray) -> Tuple[float, float, float]:
    # Compute quaternion polynomial coefficients.
    c2 = -2 * np.trace(M.T @ M)
    c1 = -8 * np.linalg.det(M)  # TODO: Slow?
    c0 = np.linalg.det(K)  # TODO: Slow?
    return c2, c1, c0


def lambda_max(Ga: float, Gb: float, c2: float, c1: float, c0: float) -> float:
    # Find largest root of the quaternion polynomial.
    def P(x):
        # Quaternion polynomial
        return x ** 4 + c2 * x ** 2 + c1 * x + c0
    def dP(x):
        # Fist derivative of the quaternion polynomial
        return 4 * x ** 3 + 2 * c2 * x + c1
    x0 = (Ga + Gb) * 0.5
    lmax = optimize.newton(P, x0, fprime=dP)
    return lmax


def _lambda_max_eig(K: np.ndarray) -> float:
    # Find largest eigenvalue of :math:`K`.
    e, _ = np.linalg.eig(K)
    return max(e)

def standard_svd(matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    # Perform singular value decomposition on the input matrix. 
    # The singular values are provided as a diagonal matrix.
    u, s, v = np.linalg.svd(matrix)
    S = np.zeros((u.shape[1],v.shape[0]))
    S[:len(s), :len(s)] = np.diag(s)
    return u, S, v

def qcp_rmsd(a: np.ndarray, b: np.ndarray, atol: float = 1e-9) -> RMSDResult:
    # Compute RMSD using the quaternion polynomial method.
    
    A = center(a)
    B = center(b)

    assert A.shape == B.shape

    N = A.shape[0]

    Ga = np.trace(A.T @ A)
    Gb = np.trace(B.T @ B)

    M = M_mtx(A, B)
    K = K_mtx(M)

    c2, c1, c0 = coefficients(M, K)

    try:
        # Fast calculation of the largest eigenvalue of K as root of the characteristic
        # polynomial.
        l_max = lambda_max(Ga, Gb, c2, c1, c0)
    except RuntimeError:  # Newton method fails to converge; see GitHub Issue #35
        # Fallback to (slower) explicit calculation of the largest eigenvalue of K
        l_max = _lambda_max_eig(K)

    s = Ga + Gb - 2 * l_max

    if abs(s) < atol:  # Avoid numerical errors when Ga + Gb = 2 * l_max
        rmsd = 0.0
    else:
        rmsd = np.sqrt(s / N)

    return RMSDResult(rmsd=rmsd)

def kabsch_rmsd(a: np.ndarray, b: np.ndarray):
    # Compute RMSD using the Kabsch algorithm
    center_a = center_of_geometry(a)
    center_b = center_of_geometry(b)
    p = a - center_a
    q = b - center_b
    dot_prod = np.dot(p.T, q)
    v, s, wt = standard_svd(dot_prod)
    d = np.sign(np.linalg.det(np.dot(v, wt).T))
    e = np.eye(3)
    e[2,2] = d
    u = np.dot(np.dot(wt.T,e),v.T)     # rotation matrix
    q = np.dot(q, u)      # transformed q coordinates
    rmsd = np.linalg.norm(p - q) / np.sqrt(p.shape[0])
    return RMSDResult(rmsd=rmsd, 
                      transition=center_b - center_a, 
                      rotation=u, 
                      transformed=q + center_a)
    
def no_alignment_rmsd(a: np.ndarray, b: np.ndarray):
    # Compute RMSD without alignment
    rmsd = np.linalg.norm(a - b) / np.sqrt(a.shape[0])
    return RMSDResult(rmsd=rmsd)

RMSD_CALCULATORS = {
    'kabsch': kabsch_rmsd,
    'qcp': qcp_rmsd,
    'no_alignment': no_alignment_rmsd
}