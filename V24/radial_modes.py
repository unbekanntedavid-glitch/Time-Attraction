"""Fundamental radial mode of a static spherical star: Rayleigh-Ritz / linear finite elements for
Chandrasekhar's (1964) Sturm-Liouville problem, written for xi = dr/r (regular at the centre):
  zeta = r^3 e^{-nuC/2} xi, nuC = 2 nu (g_tt=-e^{nuC}), lambda=-ln(1-2m/r)
  omega^2 int W zeta^2 = int [Pi zeta'^2 - Q zeta^2]
  Pi = e^{(lam+3nuC)/2} r^-2 Gamma P
  Q  = e^{(lam+3nuC)/2} r^-2 (P+eps) [nuC'^2/4 + 2 nuC'/r - 8 pi e^lam P]
  W  = e^{(3lam+nuC)/2} r^-2 (P+eps)
Natural (free) boundary at the surface (Lagrangian Delta P=0), regularity at the centre.
Gamma = (eps+P)/P dP/deps (equilibrium index; then omega^2=0 must coincide with dM/drho_c=0 in GR).
Units: km (G=c=1). omega^2 in km^-2 (observer at infinity); f_kHz = c*sqrt(omega^2)/(2 pi)/1e3.
Variant 'matter' uses matter eps,P and the actual lapse; variant 'metric' replaces (eps+P) and P in
Q, W by the effective Einstein-tensor sources (eps_eff, p_r_eff) of the modified background."""
import numpy as np
from scipy.linalg import eigh, eig, cholesky, LinAlgError, solve_banded
C_KMS = 299792.458

class IllPosed(ValueError):
    pass

def coefficients(p, r, variant='matter'):
    eos = p['eos']; h = p['h'](r); mu = p['mu'](r) if 'mu' in p else p['m'](r)/r**3; m = r**3*mu
    e, P, _, cs2 = eos.at(h)
    f = p['f'](r)
    A = r*(mu+4*np.pi*P)/(1-2*mu*r*r)          # regular form of (m+4pi r^3 P)/(r(r-2m))
    nuCp = 2*A*(f if p['mode'] == 'gradient' else 1.)
    nuC = 2*p['nu'](r)
    lam = -np.log(1-2*mu*r*r)
    Gam = cs2*(e+P)/P
    ee, pp = e, P
    if variant == 'metric':
        if p['mode'] == 'mass':
            ee = e*f
        elif p['mode'] == 'gradient':
            pp = f*P-(1-f)*mu/(4*np.pi)
    a = np.exp((lam+nuC)/2)*r**2*Gam*P          # Pi e^{-nuC} r^4
    b = 3 - r*nuCp/2
    q = np.exp((lam+nuC)/2)*r**4*(pp+ee)*(nuCp**2/4 + 2*nuCp/r - 8*np.pi*np.exp(lam)*pp)
    w = np.exp((3*lam-nuC)/2)*r**4*(pp+ee)
    return a, b, q, w

def nodes(R, N):
    s = np.linspace(0, 1, N+1)
    return R*(1-(1-s)**2)*0.5 + R*s*0.5          # mild refinement towards the surface

def fundamental(p, N=1200, variant='matter', nmodes=2, ng=4, trial=None, on_ill='raise', force_indefinite_solver=False):
    R = p['R']; rn = nodes(R, N)
    xg, wg = np.polynomial.legendre.leggauss(ng); t = (xg+1)/2; wg = wg/2
    dr = np.diff(rn)
    rg = (rn[:-1, None] + dr[:, None]*t[None, :])
    a, b, q, w = coefficients(p, rg.ravel(), variant)
    a, b, q, w = [z.reshape(rg.shape) for z in (a, b, q, w)]
    wneg = w <= 0
    if np.any(wneg) and on_ill == 'raise':
        raise IllPosed(f'non-positive inertia weight W (eps+p_eff<=0) on {np.mean(w <= 0):.3%} of quadrature points; outermost r/R with W<=0 from {rg[w <= 0].min()/R:.4f}')
    phi0 = 1-t; phi1 = t; d0 = -1/dr[:, None]; d1 = 1/dr[:, None]
    def g(i):  # operator (b xi + r xi') applied to basis i
        return (b*(phi0 if i == 0 else phi1) + rg*(d0 if i == 0 else d1))
    K = np.zeros((N+1, N+1)); Mm = np.zeros((N+1, N+1))
    for i in (0, 1):
        for j in (0, 1):
            phi_i = phi0 if i == 0 else phi1; phi_j = phi0 if j == 0 else phi1
            kij = np.sum(wg*(a*g(i)*g(j) - q*phi_i*phi_j), axis=1)*dr
            mij = np.sum(wg*(w*phi_i*phi_j), axis=1)*dr
            idx = np.arange(N)
            K[idx+i, idx+j] += kij; Mm[idx+i, idx+j] += mij
    if trial is not None:   # Rayleigh quotient of a prescribed nodal trial vector xi(r_nodes)
        v = trial(rn); return dict(rayleigh=float(v @ K @ v / (v @ Mm @ v)), r=rn)
    if not np.any(wneg) and not force_indefinite_solver:
        vals, vecs = eigh(K, Mm, subset_by_index=[0, nmodes-1])
        v = vecs[:, 0]
        return dict(omega2=vals, xi=v, r=rn, nodes_fundamental=_nodes(v), pencil_definite=True, w_nonpos_fraction=0.,
                    solver='eigh(K,M)', n_negative_type=0, n_complex=0, omega2_negative_type_max=np.nan)
    if not np.any(wneg):   # validation path only: definite pencil pushed through the indefinite solver
        return _indefinite(K, Mm, rn, nmodes, 0., np.nan, definite=True)
    return _indefinite(K, Mm, rn, nmodes, float(np.mean(wneg)), float(rg[wneg].min()/R))

def _band(A):   # tridiagonal FE matrix -> (1,1) banded storage
    n = A.shape[0]; ab = np.zeros((3, n)); ab[0, 1:] = np.diag(A, 1); ab[1] = np.diag(A); ab[2, :-1] = np.diag(A, -1)
    return ab

def _inv_iter(K, Mm, lam, its=3):
    """Eigenvector of the real eigenvalue lam of K v = lam M v by inverse iteration (tridiagonal solves)."""
    sh = lam + 1e-10*max(1e-6, abs(lam)); ab = _band(K - sh*Mm); n = K.shape[0]
    v = np.cos(np.arange(n)*0.7311) + 0.1
    for _ in range(its):
        v = solve_banded((1, 1), ab, Mm @ v); v /= np.linalg.norm(v)
    return v

def _nodes(v):
    return int(np.sum(np.diff(np.sign(v[np.abs(v) > 1e-12*np.abs(v).max()])) != 0))

def _indefinite(K, Mm, rn, nmodes, frac, rmin, definite=False):
    """Run 3: the inertia matrix M is indefinite (W<=0 somewhere), so K v = omega^2 M v is NOT a
    symmetric-definite problem and has no variational (Rayleigh-Ritz) lower-mode meaning. It is still solved
    in full. Each real eigenpair is classified by the sign of its kinetic norm v^T M v ('positive type' =
    positive kinetic energy, 'negative type' = negative kinetic energy). Reported omega2 = the smallest real
    positive-type eigenvalues (this reduces to the ordinary fundamental mode when M is positive definite).
    If K is positive definite, the equivalent symmetric problem M v = mu K v (mu = 1/omega^2) is used and all
    eigenvalues are real; otherwise the general QZ solver (scipy.linalg.eig) is used and complex pairs are counted."""
    try:
        cholesky(K, lower=True); Kpd = True
    except LinAlgError:
        Kpd = False
    if Kpd:
        mu, V = eigh(Mm, K)
        good = np.abs(mu) > 1e-14*np.abs(mu).max()
        om = 1/mu[good]; V = V[:, good]; typ = np.sign(mu[good]); ncx = 0; solver = 'eigh(M,K), K pos.-definite'
    else:
        lam = eig(K, Mm, right=False)              # QZ eigenvalues of the full pencil
        lam = lam[np.isfinite(lam)]
        isreal = np.abs(lam.imag) <= 1e-9*np.maximum(1., np.abs(lam.real))
        ncx = int(np.sum(~isreal)); om = np.sort(lam.real[isreal])
        V = np.column_stack([_inv_iter(K, Mm, x) for x in om]) if om.size else np.zeros((K.shape[0], 0))
        typ = np.sign(np.sum(V*(Mm@V), axis=0)); solver = 'eig(K,M) QZ eigenvalues + banded inverse iteration, K indefinite'
    pos = typ > 0; neg = typ < 0
    o = np.argsort(om[pos]); omp = om[pos][o]; Vp = V[:, pos][:, o]
    vals = np.full(nmodes, np.nan); vals[:min(nmodes, omp.size)] = omp[:nmodes]
    v = Vp[:, 0] if omp.size else np.full(K.shape[0], np.nan)
    return dict(omega2=vals, xi=v, r=rn, nodes_fundamental=_nodes(v) if omp.size else -1, pencil_definite=False,
                w_nonpos_fraction=frac, w_nonpos_outermost_from_rR=rmin, solver=solver, K_positive_definite=Kpd,
                n_negative_type=int(np.sum(neg)), n_complex=ncx,
                omega2_negative_type_max=float(om[neg].max()) if np.any(neg) else np.nan)

def f_kHz(omega2):
    return np.sign(omega2)*C_KMS*np.sqrt(np.abs(omega2))/(2*np.pi)/1e3
