#!/usr/bin/env python3
"""Conditional TAL stellar structure audit. Python>=3.10; numpy scipy pandas matplotlib.
Package: prepared for Zenodo V21 (GitHub V19). Licence CC BY-NC 4.0, rights held by David Szabolcsi.
V21: NICER overlay read from observational_inputs.csv; c_s^2=dP/deps causal-onset flags; dM/drho_c heuristic.
Run: python3 user/workspace/co-scientist-outputs/tal-compact-stars/stellar_tal.py
No network required. Results written alongside this file. See note_HU_EN.md.
EOS: Read+2009 PRD79 124032, arXiv:0812.2163 Table II/III.
rho in TAL and dm/dr = total energy density / c^2, not rest density.
Euclidean top-hat in areal coordinates is a declared convention, not covariant.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize_scalar,brentq
from numpy.polynomial.legendre import leggauss
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

OUT=Path(__file__).resolve().parent
G=6.67430e-11; C=299792458.; MS=1.988409870698051e30; RS=6.957e8
GM=G*MS/C**2/1000 # solar mass in km
# cgs density -> km^-2; p/c^2 has same units as cgs density
CONV=G/C**2*1e9
CCGS=(C*100)**2
PARAM={'SLy':(34.384,3.005,2.988,2.851),'APR4':(34.269,2.830,3.445,3.348),'H4':(34.669,2.909,2.246,2.144)}
LS=[.001,.367,1.,10.,100.] # km

class EOS:
    def __init__(self,name):
        lp,g1,g2,g3=PARAM[name]
        # Read Table II coefficients are K/c^2 in cgs density units.
        gc=np.array([1.58425,1.28733,.62223,1.35692])
        bc=np.array([2.44034e7,3.78358e11,2.62780e12])
        kc=[6.80110e-9]
        for i,b in enumerate(bc): kc.append(kc[-1]*b**(gc[i]-gc[i+1]))
        k1=10**lp/CCGS/(10**14.7)**g1
        match=(kc[-1]/k1)**(1/(g1-gc[-1]))
        assert bc[-1]<match<10**14.7
        self.b=np.r_[bc,match,10**14.7,1e15]
        self.g=np.r_[gc,g1,g2,g3]
        self.k=np.r_[kc,k1,k1*(10**14.7)**(g1-g2),k1*(10**14.7)**(g1-g2)*(1e15)**(g2-g3)]
        self.a=np.zeros(7)
        for i,b in enumerate(self.b):
            self.a[i+1]=self.a[i]+self.k[i]*b**(self.g[i]-1)/(self.g[i]-1)-self.k[i+1]*b**(self.g[i+1]-1)/(self.g[i+1]-1)
        self.hb=np.log(1+self.a[:-1]+self.g[:-1]/(self.g[:-1]-1)*self.k[:-1]*self.b**(self.g[:-1]-1))
    def h_causal(self,lo=5.,hi=16.5,n=20001):
        """Enthalpy at which c_s^2/c^2 = dP/d(epsilon) first reaches 1 (bracketed brentq in h); nan if none."""
        if not hasattr(self,'_hc1'):
            lr=np.linspace(lo,hi,n); hh=self.h_of_rho(10**lr); cs=self.at(hh)[3]; k=np.flatnonzero(cs>1.)
            if k.size==0: self._hc1=np.nan
            else:
                i=k[0]; f=lambda h: self.at(h)[3]-1.
                self._hc1=float(brentq(f,hh[i-1],hh[i],xtol=1e-15,rtol=1e-15)) if f(hh[i-1])<0 else float(hh[i])
        return self._hc1
    def h_of_rho(self,rho):
        i=np.searchsorted(self.b,rho)
        return np.log(1+self.a[i]+self.g[i]/(self.g[i]-1)*self.k[i]*rho**(self.g[i]-1))
    def at(self,h):
        h=np.maximum(h,1e-14)
        i=np.searchsorted(self.hb,h); g=self.g[i]; k=self.k[i]; a=self.a[i]
        rho=((np.expm1(h)-a)*(g-1)/(g*k))**(1/(g-1))
        p=k*rho**g; e=(1+a)*rho+p/(g-1)
        return e*CONV,p*CONV,rho,g*p/(e+p)


def core_boundary(eos,sol,hc,dh,r0):
    """Areal radius of the c_s^2/c^2=1 surface: dense ODE solution evaluated at the causality-threshold
    enthalpy h* (single acausal interval for all three fits). 0 if h_c<=h*. Inside the central
    expansion (h*>hc-dh) uses h ~ hc-(dh/r0^2) r^2."""
    hs=eos.h_causal()
    if not (hs==hs) or hc<=hs: return 0.
    if hs>hc-dh: return float(r0*np.sqrt((hc-hs)/dh))
    return float(sol.sol(hs)[0])


def average(rgrid,rho,L,queries=None,nq=80):
    """Exact shell-intersection geometry, numerical 1D Gauss quadrature.
    rho=0 outside star; r,L same length units. No central-density proxy.
    Integrates in normalized t=(s-r)/L to avoid small-L subtraction.
    """
    x,w=leggauss(nq); q=rgrid if queries is None else np.asarray(queries)
    R=rgrid[-1]; ans=np.zeros_like(q,dtype=float)
    for j,r in enumerate(q):
        if r<1e-10:
            b=min(L,R); s=b*(x+1)/2
            ans[j]=3/L**3*b/2*np.sum(w*s*s*np.interp(s,rgrid,rho))
            continue
        # Full shells s < L-r.
        b=min(max(L-r,0.),R)
        if b>0:
            s=b*(x+1)/2
            ans[j]+=3/L**3*b/2*np.sum(w*s*s*np.interp(s,rgrid,rho))
        # Partial shells: cap solid angle pi*(L^2-(s-r)^2)/(r*s).
        lo=max(-1.,1-2*r/L); hi=min(1.,(R-r)/L)
        if hi>lo:
            t=(hi-lo)/2*x+(hi+lo)/2; s=r+L*t
            ans[j]+=3/(4*r)*(hi-lo)/2*np.sum(w*s*np.interp(s,rgrid,rho)*(1-t*t))
    return ans


def integrate(eos,logrho,mode='GR',fprofile=None,rtol=1e-10,ngrid=801):
    hc=eos.h_of_rho(10**logrho); ec,pc,_,_=eos.at(hc)
    f0=1. if fprofile is None else float(fprofile(0.))
    fm=f0 if mode=='mass' else 1.; fg=f0 if mode=='gradient' else 1.
    r0=1e-4; m0=4*np.pi/3*ec*fm*r0**3
    dh=2*np.pi/3*(ec*fm+3*pc)*fg*r0*r0
    def rhs(h,y):
        r,m=y; e,p,_,_=eos.at(h)
        f=1. if fprofile is None else float(fprofile(r))
        grad=f if mode=='gradient' else 1.
        dr=-r*(r-2*m)/(m+4*np.pi*r**3*p)/grad
        return [dr,4*np.pi*r*r*e*(f if mode=='mass' else 1.)*dr]
    sol=solve_ivp(rhs,(hc-dh,1e-10),[r0,m0],rtol=rtol,atol=[1e-9,1e-10],dense_output=True,max_step=hc/100)
    if not sol.success: raise RuntimeError(sol.message)
    R,M=sol.y[:,-1]
    hs=np.r_[hc,np.linspace(hc-dh,1e-10,2401)]
    rr=np.r_[0.,sol.sol(hs[1:])[0]]
    ee=eos.at(hs)[0]
    rg=np.linspace(0,R,ngrid)
    er=PchipInterpolator(rr,ee)(rg); er[-1]=0
    # V21: c_s^2/c^2 = dP/d(epsilon) (=Gamma P/(epsilon+P) for each polytrope piece) at every
    # accepted ODE step (sol.t), every output node (hs) and both sides of each EOS segment boundary.
    hb=eos.hb[eos.hb<hc]
    hh=np.r_[sol.t,hs,hb*(1-1e-8),hb*(1+1e-8)]
    rh=np.r_[sol.y[0],rr,sol.sol(np.minimum(hb*(1-1e-8),hc-dh))[0] if hb.size else [],sol.sol(np.minimum(hb*(1+1e-8),hc-dh))[0] if hb.size else []]
    cs=eos.at(hh)[3]; bad=cs>1.
    return dict(R=R,M=M/GM,rg=rg,er=er,hc=hc,cs2c=float(eos.at(hc)[3]),
                cs2max=float(cs.max()),n_acausal_steps=int(bad.sum()),n_cs2_steps=int(cs.size),
                r_acausal_core_sampled=float(rh[bad].max()) if bad.any() else 0.,
                r_acausal_core=core_boundary(eos,sol,hc,dh,r0))


def star(eos,logrho,mode,L=0,rtol=1e-10,ngrid=801,nq=80,tol=2e-7):
    base=integrate(eos,logrho,rtol=rtol,ngrid=ngrid)
    if mode=='GR': return base|dict(iterations=0,residual=0.,xmax=0.)
    prev=base
    for it in range(1,31):
        av=average(prev['rg'],prev['er'],L,nq=nq)
        xp=L*L*av
        # interpolate exponent; outside grid irrelevant except tiny surface movement.
        xinterp=PchipInterpolator(prev['rg'],xp)
        oldR=prev['R']
        fp=lambda r: np.exp(-xinterp(min(r,oldR)))
        cur=integrate(eos,logrho,mode,fp,rtol,ngrid)
        res=max(abs(cur['R']/prev['R']-1),abs(cur['M']/prev['M']-1))
        if res<tol:
            # true fixed-point residual of the exponent on the new profile
            xc=L*L*average(cur['rg'],cur['er'],L,nq=nq)
            xold=np.interp(cur['rg'],prev['rg'],xp)
            fres=float(np.max(abs(xc-xold)))
            if fres<tol:
                return cur|dict(iterations=it,residual=fres,xmax=float(xc.max()))
        prev=cur
    raise RuntimeError(f'fixed point failed {logrho} {mode} {L}: {res}')


def record(s,eos,mode,L,lr):
    return dict(EOS=eos,mode=mode,L_m=L*1000,log10_rho0c_cgs=lr,M_Msun=s['M'],R_km=s['R'],xmax=s['xmax'],iterations=s['iterations'],fp_residual=s['residual'],cs2c_over_c2=s['cs2c'],
                cs2_max_over_c2_all_steps=s['cs2max'],acausal_point=bool(s['cs2max']>1.),r_acausal_core_km=s['r_acausal_core'],
                r_acausal_core_sampled_lower_km=s['r_acausal_core_sampled'])


def sequences(names=None):
    rows=[]; maxima=[]; tests=[]
    suffix='' if names is None else '_'+names[0]
    for name in (PARAM if names is None else names):
        eos=EOS(name)
        for mode,L in [('GR',0.)]+[(m,L) for m in ['mass','gradient'] for L in LS]:
            t=time.time(); cache={}
            def calc(lr):
                key=round(float(lr),10)
                if key not in cache: cache[key]=star(eos,float(lr),mode,L)
                return cache[key]
            grid=np.linspace(14.55,15.65,53)
            ss=[calc(lr) for lr in grid]; mm=np.array([s['M'] for s in ss]); im=mm.argmax()
            assert 0<im<len(grid)-1
            opt=minimize_scalar(lambda z:-calc(z)['M'],bounds=(grid[im-1],grid[im+1]),method='bounded',options={'xatol':2e-5})
            lrmax=float(opt.x); sm=calc(lrmax)
            # refine radii by root solve, not coarse M-R interpolation
            rec=record(sm,name,mode,L,lrmax)
            for target in [1.4,1.418,1.43,2.,2.073]:
                label=str(target).replace('.','p')
                if mm[0]<target<sm['M']:
                    lr=brentq(lambda z:calc(z)['M']-target,grid[0],lrmax,xtol=2e-6)
                    rec['R_at_'+label+'_km']=calc(lr)['R']
                    rec['log10_rho0c_at_'+label]=lr
                else: rec['R_at_'+label+'_km']=np.nan; rec['log10_rho0c_at_'+label]=np.nan
            maxima.append(rec)
            for lr,s in sorted(cache.items()): rows.append(record(s,name,mode,L,lr))
            pd.DataFrame(rows).to_csv(OUT/('mr_sequences'+suffix+'.csv'),index=False)
            pd.DataFrame(maxima).to_csv(OUT/('stellar_summary'+suffix+'.csv'),index=False)
            print(name,mode,L,'Mmax',round(sm['M'],6),'R1.4',rec['R_at_1p4_km'],'secs',round(time.time()-t,1),flush=True)
            # Convergence checks at peak and 1.4 density for largest correction + small-L.
            if mode=='GR' or (mode!='GR' and L in [.367,10.]):
                fine=star(eos,lrmax,mode,L,rtol=1e-11,ngrid=1601,nq=160,tol=2e-7)
                tests.append(dict(EOS=name,mode=mode,L_m=L*1000,logrho=lrmax,dM_Msun=fine['M']-sm['M'],dR_km=fine['R']-sm['R'],fp_residual=fine['residual']))
    pd.DataFrame(tests).to_csv(OUT/('convergence'+suffix+'.csv'),index=False)


def causal_onset(eos,lo=5.,hi=16.5,n=20001):
    """Lowest rest-mass density (log10 g/cm^3) at which c_s^2/c^2 = dP/d(epsilon) exceeds 1.
    Scans a fine grid, then root-solves inside the bracketing EOS piece. nan if causal up to 10**hi."""
    lr=np.linspace(lo,hi,n); cs=eos.at(eos.h_of_rho(10**lr))[3]
    k=np.flatnonzero(cs>1.)
    if k.size==0: return np.nan,np.nan,int(0)
    i=k[0]; f=lambda z: eos.at(eos.h_of_rho(10**z))[3]-1.
    z=brentq(f,lr[i-1],lr[i],xtol=1e-12) if f(lr[i-1])<0 else lr[i]
    e=eos.at(eos.h_of_rho(10**z))[0]/CONV
    # number of separate acausal intervals on the scan (re-entry check)
    nint=int(np.sum(np.diff((cs>1.).astype(int))==1)+(cs[0]>1.))
    return float(z),float(np.log10(e)),nint


def causality_stability():
    """V21 addendum: causal-onset flags and turning-point (dM/drho_c>0) heuristic for all sequences.
    Adds columns to stellar_summary.csv and writes causality_stability_log.txt."""
    su=pd.read_csv(OUT/'stellar_summary.csv'); df=pd.read_csv(OUT/'mr_sequences.csv'); log=[]
    log.append('V21 causality / turning-point log. c_s^2/c^2 = dP/d(epsilon) evaluated at every accepted ODE step, '
               'every output node and both sides of every EOS piece boundary of each solved star. '
               'Core boundary = areal radius where c_s^2/c^2=1, from the dense ODE solution at the threshold enthalpy.')
    onset={}
    for name in PARAM:
        z,ze,nint=causal_onset(EOS(name)); onset[name]=(z,ze,nint)
        log.append(f'{name}: first c_s^2>c^2 at log10 rho0={z:.6f} (rho0={10**z if z==z else np.nan:.6e} g/cm3), '
                   f'log10 eps/c^2={ze:.6f}; separate acausal intervals on 1e5-10^16.5 g/cm3 scan: {nint}')
    for j,row in su.iterrows():
        eos=EOS(row.EOS); mode=row['mode']; L=row.L_m/1000; z,ze,_=onset[row.EOS]
        pk=star(eos,float(row.log10_rho0c_cgs),mode,L)
        flag=bool(pk['cs2max']>1.)
        su.loc[j,'causal_violation_flag']=flag
        su.loc[j,'cs2_max_over_c2_in_Mmax_star']=pk['cs2max']
        su.loc[j,'cs2_onset_log10_rho0_cgs']=z
        su.loc[j,'cs2_onset_log10_eps_over_c2_cgs']=ze
        su.loc[j,'cs2_onset_percent_of_rho0c_at_Mmax']=100*10**(z-row.log10_rho0c_cgs) if flag else np.nan
        su.loc[j,'r_acausal_core_at_Mmax_km']=pk['r_acausal_core']
        su.loc[j,'r_acausal_core_sampled_lower_at_Mmax_km']=pk['r_acausal_core_sampled']
        if flag:
            cz=star(eos,z,mode,L); mc,rc,lc=cz['M'],cz['R'],z
        else: mc,rc,lc=row.M_Msun,row.R_km,row.log10_rho0c_cgs
        su.loc[j,'M_max_causal_Msun']=mc; su.loc[j,'R_at_M_max_causal_km']=rc; su.loc[j,'log10_rho0c_at_M_max_causal']=lc
        # turning-point heuristic dM/drho_c>0 on the 53-point scan grid; the refined peak gives the edge.
        s=df[(df.EOS==row.EOS)&(df['mode']==mode)&(df.L_m==row.L_m)].sort_values('log10_rho0c_cgs')
        grid=np.round(np.linspace(14.55,15.65,53),10)
        g=s[np.isin(np.round(s.log10_rho0c_cgs.to_numpy(),10),grid)]
        lr=g.log10_rho0c_cgs.to_numpy(); m=g.M_Msun.to_numpy(); dm=np.diff(m)/np.diff(lr)
        changes=int(np.sum(np.sign(dm[1:])!=np.sign(dm[:-1])))
        # all computed points (incl. optimizer/root-finder points): flips within noise near the peak
        la=s.log10_rho0c_cgs.to_numpy(); ma=s.M_Msun.to_numpy(); da=np.diff(ma)
        fl=np.flatnonzero(np.sign(da[1:])!=np.sign(da[:-1]))
        far=[i for i in fl if abs(la[i+1]-row.log10_rho0c_cgs)>2e-4]
        su.loc[j,'stable_branch_log10_rho0c_min']=lr[0]; su.loc[j,'stable_branch_log10_rho0c_max']=row.log10_rho0c_cgs
        su.loc[j,'stable_branch_M_min_Msun']=m[0]; su.loc[j,'stable_branch_M_max_Msun']=row.M_Msun
        su.loc[j,'dMdrho_sign_changes_scan_grid']=changes
        su.loc[j,'dMdrho_sign_changes_all_points']=len(fl)
        su.loc[j,'dMdrho_flips_outside_peak_2e-4dex']=len(far)
        # size of the reversal at each flip: smaller of the two adjacent |Delta M| steps (compare 6.4e-7 Msun refinement floor)
        su.loc[j,'dMdrho_flip_max_reversal_dM_Msun']=float(max([min(abs(da[i]),abs(da[i+1])) for i in fl]+[0.]))
        for lab in ['1p4','2p0']:
            lt=row['log10_rho0c_at_'+lab]
            if lt==lt:
                k=np.searchsorted(lr,lt); k=min(max(k,1),len(lr)-1)
                st=bool(lr[0]<=lt<row.log10_rho0c_cgs and dm[k-1]>0)
                su.loc[j,'stable_branch_at_'+lab]=st
                su.loc[j,'causal_at_'+lab]=bool(not (lt>=z)) if z==z else True
            else:
                su.loc[j,'stable_branch_at_'+lab]=np.nan; su.loc[j,'causal_at_'+lab]=np.nan
        r=su.loc[j]
        log.append(f"{row.EOS:4s} {mode:8s} L={row.L_m:8g} m: Mmax={row.M_Msun:.6f} at log10 rho0c={row.log10_rho0c_cgs:.6f}; "
                   f"cs2max={pk['cs2max']:.6f}; acausal={flag}"
                   +(f"; onset at {r.cs2_onset_percent_of_rho0c_at_Mmax:.3f}% of rho0c, acausal core boundary (c_s^2=c^2) r={pk['r_acausal_core']:.6f} km (sampled lower estimate {pk['r_acausal_core_sampled']:.6f}); "
                     f"M_max_causal={mc:.6f} R={rc:.4f} km" if flag else '')
                   +f"; R1.4={row.R_at_1p4_km:.4f} (causal {r.causal_at_1p4}, stable {r.stable_branch_at_1p4})"
                   +f"; R2.0={row.R_at_2p0_km:.4f} (causal {r.causal_at_2p0}, stable {r.stable_branch_at_2p0})"
                   +f"; dM/drho_c sign changes on scan grid={changes}, on all points={len(fl)} (outside +-2e-4 dex of peak: {len(far)})")
    su.to_csv(OUT/'stellar_summary.csv',index=False)
    log.append('Turning-point criterion dM/drho_c>0 is a GR-motivated heuristic; for the modified (mass/gradient) equations '
               'it is NOT a proven radial-stability result. A radial-mode (Chandrasekhar/Sturm-Liouville) analysis remains open.')
    (OUT/'causality_stability_log.txt').write_text('\n'.join(log)+'\n')
    print('\n'.join(log))


def lane_emden(n,N=2401):
    x0=1e-6
    def rhs(x,y): return [y[1],-2*y[1]/x-max(y[0],0.)**n]
    def end(x,y): return y[0]
    end.terminal=True;end.direction=-1
    sol=solve_ivp(rhs,(x0,20),[1-x0*x0/6,-x0/3],events=end,rtol=2e-10,atol=1e-12,dense_output=True,max_step=.02)
    x1=sol.t_events[0][0]; x=x1*(1-(1-np.linspace(0,1,N))**3)
    d=np.maximum(sol.sol(np.maximum(x,x0))[0],0.)**n;d[-1]=0
    d/=4*np.pi*np.trapezoid(d*x*x,x)
    return x/x1,d*x1**3


def wd():
    M=1.018;sigM=.011;R=.00803*RS;sigR=.00011*RS;vobs=80.65;sigv=.77
    comp=G*M*MS/(R*C*C);zgr=(1-2*comp)**-.5-1;vgr=C*zgr/1000
    siggr=vgr*np.hypot(sigM/M,sigR/R);sigma=np.hypot(sigv,siggr);res=vobs-vgr
    # Second, conservative allowance using paper's quoted redshift mass error.
    mass_allowance=vgr*((1.017-M)+1.96*np.hypot(.025,sigM))/M
    # Additive local clock multiplication: 1+z=(1+z_GR)exp(x_em-x_obs).
    # Differential A/B measurement: compare extra B shift minus A shift.
    rrA,ddA=lane_emden(3.);RA=1.711*RS;MA=2.042
    lvals=np.unique(np.r_[np.geomspace(1,1e10,361),[1,367,1e3,1e4,1e5]])
    records=[];summary=[]
    for n in [1.5,3.]:
        rr,dd=lane_emden(n)
        for L in lvals:
            rb=average(rr*R,dd*M*MS/R**3,L,queries=[R],nq=160)[0]
            ra=average(rrA*RA,ddA*MA*MS/RA**3,L,queries=[RA],nq=160)[0]
            xb=G*L*L/C**2*rb;xa=G*L*L/C**2*ra
            vb=C/1000*(1+zgr)*np.expm1(xb); va=C/1000*np.expm1(xa)
            records.append(dict(n_polytrope=n,L_m=L,xB_surface=xb,xA_surface=xa,delta_vB_kms=vb,delta_vA_kms=va,delta_v_differential_kms=vb-va,excluded_approx95=(vb-va>res+1.96*sigma),excluded_quoted_mass95=(vb-va>mass_allowance)))
        sub=pd.DataFrame(records);sub=sub[sub.n_polytrope==n];ex=sub[sub.excluded_approx95]
        ex2=sub[sub.excluded_quoted_mass95]
        summary.append(dict(n=n,max_v_kms=float(sub.delta_v_differential_kms.max()),peak_L_m=float(sub.loc[sub.delta_v_differential_kms.idxmax(),'L_m']),excluded_sampled_band_m=[] if ex.empty else [float(ex.L_m.min()),float(ex.L_m.max())],excluded_quoted_mass_sampled_band_m=[] if ex2.empty else [float(ex2.L_m.min()),float(ex2.L_m.max())]))
    pd.DataFrame(records).to_csv(OUT/'sirius_redshift.csv',index=False)
    (OUT/'wd_summary.json').write_text(json.dumps(dict(M=M,R_m=R,compactness=comp,vgr_kms=vgr,sigma_gr_kms=siggr,residual_kms=res,sigma_combined_kms=sigma,upper_extra_approx95_kms=res+1.96*sigma,quoted_mass_extra_approx95_kms=mass_allowance,models=summary),indent=2))


def check_observational_inputs(ob,tol=1e-9):
    """Each symmetrized error must equal the arithmetic mean of the published lower and upper errors.
    Raises ValueError naming every failing row; returns per-row results."""
    res=[];bad=[]
    for r in ob.itertuples():
        for q,lo,hi,sym in [('M',r.M_minus_Msun,r.M_plus_Msun,r.M_sym_err_Msun),('R',r.R_minus_km,r.R_plus_km,r.R_sym_err_km)]:
            mean=(float(lo)+float(hi))/2; d=float(sym)-mean
            res.append(dict(source=r.source,analysis=r.analysis,quantity=q,minus=float(lo),plus=float(hi),sym=float(sym),mean=mean,diff=d,pass_=abs(d)<=tol))
            if abs(d)>tol: bad.append(f"{r.source} [{r.analysis}] {q}: sym_err={float(sym)!r} but (minus+plus)/2={mean!r}")
    if bad: raise ValueError('observational_inputs.csv: symmetrized error is not the mean of lower/upper errors:\n  '+'\n  '.join(bad))
    return res


def _causal_split(df,su,name,mode,Lm):
    """Return (causal part, acausal continuation, endpoint or None) of one sequence up to its maximum."""
    row=su[(su.EOS==name)&(su['mode']==mode)&(su.L_m==Lm)].iloc[0]
    s=df[(df.EOS==name)&(df['mode']==mode)&(df.L_m==Lm)].sort_values('log10_rho0c_cgs')
    s=s[s.log10_rho0c_cgs<=row.log10_rho0c_cgs]
    if not bool(row.causal_violation_flag): return s[['R_km','M_Msun']].to_numpy(),None,None
    z=row.cs2_onset_log10_rho0_cgs; ep=np.array([[row.R_at_M_max_causal_km,row.M_max_causal_Msun]])
    c=np.vstack([s[s.log10_rho0c_cgs<z][['R_km','M_Msun']].to_numpy(),ep])
    a=np.vstack([ep,s[s.log10_rho0c_cgs>z][['R_km','M_Msun']].to_numpy()])
    return c,a,ep[0]


def plots():
    df=pd.read_csv(OUT/'mr_sequences.csv');su=pd.read_csv(OUT/'stellar_summary.csv')
    cols={'SLy':'#0072B2','APR4':'#D55E00','H4':'#009E73'}
    # V21: all NICER reference values come from observational_inputs.csv (no hard-coded numbers).
    ob=pd.read_csv(OUT/'observational_inputs.csv',keep_default_na=False)
    chk=check_observational_inputs(ob)
    for c in chk: print(f"sym-error check {c['source']} [{c['analysis']}] {c['quantity']}: sym={c['sym']:g} mean={c['mean']:g} diff={c['diff']:.1e} {'PASS' if c['pass_'] else 'FAIL'}")
    ob=ob[ob.used_in_plot.astype(str).str.lower()=='true']
    obs=[(r.plot_label,float(r.R_km),float(r.M_Msun),float(r.R_sym_err_km),float(r.M_sym_err_Msun),r.plot_color,
          [(float(r.ellipse_joint_prob_2),float(r.ellipse_alpha_2)),(float(r.ellipse_joint_prob_1),float(r.ellipse_alpha_1))]) for r in ob.itertuples()]
    colour_key='; '.join(f"{r.plot_color_name}: {r.plot_label}" for r in ob.itertuples())+'.'
    colour_key=colour_key[:1].upper()+colour_key[1:]
    fig,axes=plt.subplots(2,5,figsize=(19,8),sharex=True,sharey=True)
    for i,mode in enumerate(['mass','gradient']):
        for j,L in enumerate(LS):
            ax=axes[i,j]
            for label,r,m,sr,sm,col,levels in obs:
                for prob,alpha in levels:
                    fac=np.sqrt(-2*np.log(1-prob))
                    ax.add_patch(Ellipse((r,m),2*fac*sr,2*fac*sm,facecolor=col,edgecolor=col,alpha=alpha,lw=.8))
                ax.plot(r,m,'+',color=col,ms=5)
            for name in PARAM:
                c,a,ep=_causal_split(df,su,name,'GR',0.)
                ax.plot(c[:,0],c[:,1],color=cols[name],ls='--',lw=1.1)
                if a is not None:
                    ax.plot(a[:,0],a[:,1],color=cols[name],ls=':',lw=1.0,alpha=.45)
                    ax.plot(*ep,marker='s',ms=5,mfc='none',mec=cols[name],mew=1.0,ls='none')
                c,a,ep=_causal_split(df,su,name,mode,L*1000)
                ax.plot(c[:,0],c[:,1],color=cols[name],lw=1.5,label=name)
                if a is not None:
                    ax.plot(a[:,0],a[:,1],color=cols[name],ls=':',lw=1.6,alpha=.55)
                    ax.plot(*ep,marker='o',ms=5.5,mfc='white',mec=cols[name],mew=1.3,ls='none')
            ax.axhline(2.,color='black',ls=':',lw=.7);ax.set_title(f'L = {L*1000:g} m');ax.set_xlim(9,16);ax.set_ylim(.85,2.45);ax.grid(alpha=.18)
            if j==0:ax.set_ylabel(('Published mass-source ansatz' if i==0 else 'Alternative gradient/lapse ansatz')+'\nM / solar mass')
            if i==1:ax.set_xlabel('R [km]')
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    h=[Line2D([],[],color=cols[n],lw=1.5,label=n) for n in PARAM]+[
       Line2D([],[],color='0.3',ls=':',lw=1.6,alpha=.6,label='acausal continuation'),
       Line2D([],[],marker='o',mfc='white',mec='0.3',ls='none',label='TAL causal endpoint'),
       Line2D([],[],marker='s',mfc='none',mec='0.3',ls='none',label='GR causal endpoint')]
    axes[0,0].legend(handles=h,fontsize=7,loc='upper left')
    axes[1,0].legend(handles=[Patch(facecolor=col,edgecolor=col,alpha=.35,label=label) for label,r,m,sr,sm,col,levels in obs],fontsize=7,loc='upper left',title='NICER (GR-derived)',title_fontsize=7)
    fl=su[su.causal_violation_flag.astype(bool)]
    acausal_txt='; '.join(f"{n}: {len(g)} of {len(su[su.EOS==n])} sequences" for n,g in fl.groupby('EOS',sort=False)) or 'none'
    fig.suptitle('Conditional TAL M–R sequences: solid TAL, dashed GR within the causal range; dotted = acausal continuation | SLy / APR4 / H4 Read fits',fontsize=13,y=.995)
    fig.text(.5,.012,'Shading: independent, symmetrized Gaussian 68% / 95% JOINT approximations (not published NICER contours). '+colour_key+'\nNonrotating TOV; GR-derived NICER posteriors are illustrative references, not a TAL likelihood. Modified-branch stability is not established.\nCurves are solid/dashed up to the causal endpoint (c_s^2=c^2 first reached at the centre; open circle TAL, open square GR); the dotted, faded continuation up to the sequence maximum\nis a potentially acausal benchmark of the published fit ('+acausal_txt+'). No causal completion was imposed.',ha='center',fontsize=8.5)
    fig.tight_layout(rect=[0,.08,1,.965]);fig.savefig(OUT/'mr_nicer.png',dpi=180);fig.savefig(OUT/'mr_nicer.pdf');plt.close(fig)
    # Fixed-mass shifts versus selected EOS envelope.
    rows=[]
    for tar in ['1p4','1p418','1p43','2p0','2p073']:
        key='R_at_'+tar+'_km'
        bg=su[su['mode']=='GR'];spread=bg[key].max()-bg[key].min() if bg[key].notna().sum()>=2 else np.nan
        for _,r in su[su['mode']!='GR'].iterrows():
            base=bg[bg.EOS==r.EOS].iloc[0]
            rows.append(dict(EOS=r.EOS,mode=r['mode'],L_m=r.L_m,M_Msun=float(tar.replace('p','.')),delta_R_km=r[key]-base[key],GR_EOS_spread_km=spread,n_GR_EOS_available=int(bg[key].notna().sum())))
    pd.DataFrame(rows).to_csv(OUT/'radius_shifts.csv',index=False)
    wddata=pd.read_csv(OUT/'sirius_redshift.csv');w=json.loads((OUT/'wd_summary.json').read_text())
    fig,ax=plt.subplots(figsize=(8,4.5))
    for n,g in wddata.groupby('n_polytrope'):ax.semilogx(g.L_m,g.delta_v_differential_kms,label=f'WD n={n:g}; Sirius A n=3')
    ax.axhline(w['upper_extra_approx95_kms'],ls='--',color='k',label='Approx. 95% extra-shift upper allowance')
    ax.set(xlabel='L [m]',ylabel='Additional differential redshift [km/s]',title='Sirius B: conditional additive-clock test, GR exterior retained')
    ax.grid(alpha=.2);ax.legend(fontsize=8);fig.tight_layout();fig.savefig(OUT/'sirius_redshift.png',dpi=180);plt.close(fig)


def tests():
    # Uniform sphere overlap: center L<R; all enclosed; at surface L=R cap fraction 5/16.
    r=np.linspace(0,10,3001);d=np.ones_like(r)
    a=average(r,d,2,queries=[0,3,10],nq=160)
    assert np.max(abs(a[:2]-1))<1e-10
    surf=average(r,d,10,queries=[10],nq=160)[0]
    assert abs(surf-5/16)<1e-10
    large=average(r,d,100,queries=[0,10],nq=160)
    assert np.max(abs(large-.001))<1e-10
    eos_tests={}
    for name in PARAM:
        e=EOS(name);rs=np.geomspace(1e5,5e15,1000);hh=e.h_of_rho(rs)
        inv=e.at(hh)[2]
        err=float(np.max(abs(inv/rs-1)))
        assert err<1e-6
        eos_tests[name]=dict(max_enthalpy_inverse_relative_error=err,core_crust_rho_cgs=float(e.b[3]),parameters=list(PARAM[name]))
    (OUT/'geometry_tests.json').write_text(json.dumps(dict(uniform_small=a.tolist(),uniform_surface_L_R=float(surf),uniform_large=large.tolist(),EOS=eos_tests),indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--plots-only',action='store_true');p.add_argument('--wd-only',action='store_true');p.add_argument('--eos',choices=list(PARAM));p.add_argument('--combine',action='store_true');p.add_argument('--causal-only',action='store_true');a=p.parse_args()
    tests()
    if a.eos:sequences([a.eos])
    elif a.combine:
        for stem in ['mr_sequences','stellar_summary','convergence']:
            pd.concat([pd.read_csv(OUT/(stem+'_'+n+'.csv')) for n in PARAM],ignore_index=True).to_csv(OUT/(stem+'.csv'),index=False)
        causality_stability();wd();plots()
    elif a.wd_only:wd()
    elif a.plots_only:plots()
    elif a.causal_only:causality_stability()
    else:sequences();causality_stability();wd();plots()
