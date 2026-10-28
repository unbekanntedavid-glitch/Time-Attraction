"""Converged TAL equilibrium profiles, re-using the unchanged V21 code (stellar_tal_v21.py,
byte-identical copy of V21 stellar_tal.py, SHA-256 cb93c90f...1906).
converged() repeats the fixed-point loop of stellar_tal.star() exactly, but also returns the
converged f-profile; dense_profile() re-integrates the same ODE system with that f-profile and
returns callables r->(h,m,eps,P,Gamma,f,nu) in geometrised km units (G=c=1).
Lapse convention as in the package note: g_tt=-e^{2 nu}."""
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.interpolate import PchipInterpolator
import stellar_tal_v21 as st

EOS = st.EOS; average = st.average; GM = st.GM

def converged(eos, logrho, mode, L, rtol=1e-10, ngrid=801, nq=80, tol=2e-7):
    base = st.integrate(eos, logrho, rtol=rtol, ngrid=ngrid)
    if mode == 'GR':
        return base, None, None
    prev = base
    for it in range(1, 31):
        av = average(prev['rg'], prev['er'], L, nq=nq)
        xp = L*L*av
        xinterp = PchipInterpolator(prev['rg'], xp)
        oldR = prev['R']
        fp = (lambda xi, oR: (lambda r: np.exp(-xi(min(r, oR)))))(xinterp, oldR)
        cur = st.integrate(eos, logrho, mode, fp, rtol, ngrid)
        res = max(abs(cur['R']/prev['R']-1), abs(cur['M']/prev['M']-1))
        if res < tol:
            xc = L*L*average(cur['rg'], cur['er'], L, nq=nq)
            xold = np.interp(cur['rg'], prev['rg'], xp)
            if float(np.max(abs(xc-xold))) < tol:
                return cur, fp, (xinterp, oldR)
        prev = cur
    raise RuntimeError('fixed point failed')

def dense_profile(eos, logrho, mode='GR', L=0., rtol=1e-10, nh=12000):
    """Return dict with star summary and interpolants of the converged profile."""
    cur, fp, xi = converged(eos, logrho, mode, L, rtol=rtol)
    hc = eos.h_of_rho(10**logrho); ec, pc, _, _ = eos.at(hc)
    f0 = 1. if fp is None else float(fp(0.))
    fm = f0 if mode == 'mass' else 1.; fg = f0 if mode == 'gradient' else 1.
    r0 = 1e-4; m0 = 4*np.pi/3*ec*fm*r0**3
    dh = 2*np.pi/3*(ec*fm+3*pc)*fg*r0*r0
    def rhs(h, y):   # identical to stellar_tal.integrate
        r, m = y; e, p, _, _ = eos.at(h)
        f = 1. if fp is None else float(fp(r))
        grad = f if mode == 'gradient' else 1.
        dr = -r*(r-2*m)/(m+4*np.pi*r**3*p)/grad
        return [dr, 4*np.pi*r*r*e*(f if mode == 'mass' else 1.)*dr]
    sol = solve_ivp(rhs, (hc-dh, 1e-10), [r0, m0], rtol=rtol, atol=[1e-9, 1e-10],
                    dense_output=True, max_step=hc/100)
    R, M = sol.y[:, -1]
    assert abs(R/cur['R']-1) < 1e-9 and abs(M/cur['M']/GM-1) < 1e-9, (R, cur['R'], M/GM, cur['M'])
    # Profile mapping, regular at the centre (audit fix O1, run 2): h and mu = m/r^3 are even, analytic
    # functions of r, so they are interpolated in s = r^2 with their exact centre values h(0)=hc,
    # mu(0)=(4pi/3) eps_c f_m. Then h'(r)=2r dh/ds -> 0, m=r^3 mu, A=r(mu+4pi P)/(1-2 mu r^2) are regular.
    # (Run 1 interpolated h(r), m(r) in r with a single node at r=0 -> wrong slopes near the centre.)
    u = np.linspace(0, 1, nh)
    hs = np.unique(np.r_[hc-dh, (hc-dh)*(1-u**1.0), (hc-dh)*np.exp(-12*u)*(1-u)+1e-10, 1e-10,
                         (hc-dh)-hc*np.geomspace(1e-12, 1e-3, 400), sol.t])[::-1]
    hs = hs[(hs <= hc-dh) & (hs >= 1e-10)]
    rr, mm = sol.sol(hs)
    o = np.argsort(rr); rr, mm, hh = rr[o], mm[o], hs[o]
    keep = np.r_[True, np.diff(rr) > 0]; rr, mm, hh = rr[keep], mm[keep], hh[keep]
    muc = 4*np.pi/3*ec*fm
    ss = np.r_[0., rr**2]; h_s = PchipInterpolator(ss, np.r_[hc, hh]); mu_s = PchipInterpolator(ss, np.r_[muc, mm/rr**3])
    h_of_r = lambda r: h_s(np.asarray(r, float)**2)
    mu_of_r = lambda r: mu_s(np.asarray(r, float)**2)
    m_of_r = lambda r: np.asarray(r, float)**3*mu_s(np.asarray(r, float)**2)
    dmu_dr = (lambda d: (lambda r: 2*np.asarray(r, float)*d(np.asarray(r, float)**2)))(mu_s.derivative())
    # f: the equilibrium ODE keeps the package f-profile fp (Pchip of x in r, so M, R are reproduced
    # exactly); for the mapped profile x is re-interpolated in s = r^2 on the same nodes (x is even in r).
    fr_pkg = (lambda r: np.ones_like(np.asarray(r, float))) if fp is None else \
             (lambda r: np.exp(-xi[0](np.minimum(np.asarray(r, float), xi[1]))))
    if fp is None:
        fr = fr_pkg; dlnf_dr = lambda r: np.zeros_like(np.asarray(r, float))
    else:
        xn = xi[0].x; xs_ = PchipInterpolator(xn**2, xi[0](xn)); dxs = xs_.derivative(); oR = xi[1]
        fr = lambda r: np.exp(-xs_(np.minimum(np.asarray(r, float), oR)**2))
        dlnf_dr = lambda r: np.where(np.asarray(r, float) < oR, -2*np.asarray(r, float)*dxs(np.minimum(np.asarray(r, float), oR)**2), 0.)
    # surface lapse: e^{2 nu(R)}; GR & mass mode exterior is Schwarzschild.
    nuR = 0.5*np.log(1-2*M/R)
    ext = None
    if mode == 'gradient':
        # exterior f from the top-hat around exterior points (rho=0 outside, ball reaches into star)
        nuR, ext = exterior_lapse(cur, L, M, R)
    def nu_of_r(r):  # nu' = A (GR, mass) or f A (gradient) = -dh/dr  ->  nu = nu(R) - h + h(R)
        return nuR - (h_of_r(r) - 1e-10)
    return dict(R=R, M=M, Msun=M/GM, hc=hc, mode=mode, L=L, logrho=logrho, sol=sol, ec=ec, pc=pc, muc=muc,
                h=h_of_r, m=m_of_r, mu=mu_of_r, dmu_dr=dmu_dr, f=fr, f_pkg=fr_pkg, dlnf_dr=dlnf_dr, nu=nu_of_r, nuR=nuR, eos=eos, cur=cur, xi=xi, ext=ext)

def exterior_lapse(cur, L, M, R, nq=160):
    """nu(R) = -int_R^inf f(r) M/(r(r-2M)) dr with f=exp(-L^2 rhobar_L) from the converged profile;
    f=1 beyond R+L. Returns nu(R) and diagnostic dict."""
    rq = R + L*np.linspace(0, 1, 401)[1:]
    xq = L*L*average(cur['rg'], cur['er'], L, queries=rq, nq=nq)
    xR = L*L*average(cur['rg'], cur['er'], L, queries=[R], nq=nq)[0]
    rq = np.r_[R, rq]; xq = np.r_[xR, xq]
    xs = PchipInterpolator(rq, xq)
    A = lambda r: M/(r*(r-2*M))
    xg, wg = np.polynomial.legendre.leggauss(200)
    rm = 0.5*(R+L - R)*xg + 0.5*(R+L + R)
    wm = 0.5*L*wg
    d = float(np.sum(wm * (-np.expm1(-xs(rm))) * A(rm)))
    nuGR = 0.5*np.log(1-2*M/R)
    return nuGR + d, dict(x_surface=float(xR), delta_nu_exterior=float(d), nu_GR_R=float(nuGR))

def surface_x(cur, L, nq=160):
    return float(L*L*average(cur['rg'], cur['er'], L, queries=[cur['R']], nq=nq)[0])
