"""Extra diagnostics for stellar_tal.py; run after the main sequence calculation.
V21 (Zenodo V21 / GitHub V19): adds causal-endpoint and core-boundary refinement tests and the NICER sym-error check.
Licence CC BY-NC 4.0, rights held by David Szabolcsi."""
import json,sys,hashlib,platform
from pathlib import Path
import numpy as np,pandas as pd,scipy,matplotlib
from scipy.integrate import simpson
from scipy.optimize import brentq
from scipy.interpolate import PchipInterpolator
import stellar_tal as s
O=Path(__file__).resolve().parent
D=pd.read_csv(O/'mr_sequences.csv');S=pd.read_csv(O/'stellar_summary.csv')
profiles=[];diagnostics=[]
for name in s.PARAM:
    e=s.EOS(name)
    for mode in ['mass','gradient']:
        sub=D[(D.EOS==name)&(D['mode']==mode)&(D.L_m==10000)]
        row=sub.iloc[np.argmin(abs(sub.M_Msun-1.4))]
        st=s.star(e,float(row.log10_rho0c_cgs),mode,10)
        x=100*s.average(st['rg'],st['er'],10,nq=160)
        bare=4*np.pi*simpson(st['er']*st['rg']**2,x=st['rg'])/s.GM
        for r,er,xx in zip(st['rg'][::4],st['er'][::4],x[::4]):
            profiles.append(dict(EOS=name,mode=mode,L_m=10000,r_km=r,rho_energy_kg_m3=er/s.CONV*1000,x=xx,f=np.exp(-xx)))
        diagnostics.append(dict(EOS=name,mode=mode,M_Msun=st['M'],R_km=st['R'],bare_coordinate_mass_Msun=bare,x_center=float(x[0]),x_surface=float(x[-1]),central_proxy_x=float(100*st['er'][0]),fp_residual=st['residual']))
# Independent spherical-cap quadrature split at every profile knot; piecewise-linear density.
# Simpson/Gauss production integral tested against analytic 4-point polynomial-exact shell integrals.
def exact_shell_average(rg,er,r,L):
    from numpy.polynomial.legendre import leggauss
    xx,ww=leggauss(4)
    cuts=np.unique(np.r_[rg,abs(r-L),r+L,max(L-r,0)])
    cuts=cuts[(cuts>=0)&(cuts<=rg[-1])]
    total=0
    for a,b in zip(cuts[:-1],cuts[1:]):
        z=(a+b)/2+(b-a)*xx/2
        if r==0:weight=4*np.pi*z*z*(z<L)
        else:
            cost=(z*z+r*r-L*L)/(2*z*r)
            weight=2*np.pi*z*z*np.clip(1-cost,0,2)
        total+=(b-a)/2*np.sum(ww*np.interp(z,rg,er)*weight)
    return total/(4*np.pi*L**3/3)
geometry_cross=[]
e=s.EOS('SLy');st=s.star(e,15.,'mass',10)
for L in [.001,.367,1,10,100]:
    for frac in [0.,.5,1.]:
        r=st['R']*frac
        v=s.average(st['rg'],st['er'],L,queries=[r],nq=80)[0]
        ref=exact_shell_average(st['rg'],st['er'],r,L)
        geometry_cross.append(dict(L_m=L*1000,r_over_R=frac,relative_error_80_vs_exact_linear=float(v/ref-1)))
# Causal flags: maximum sound speed on each maximum-mass star, including segment jumps.
for j,row in S.iterrows():
    e=s.EOS(row.EOS);hc=e.h_of_rho(10**row.log10_rho0c_cgs)
    hs=np.r_[np.linspace(1e-9,hc,2001),e.hb[e.hb<hc]*(1-1e-8),e.hb[e.hb<hc]*(1+1e-8)]
    S.loc[j,'max_cs2_over_c2']=np.max(e.at(hs)[3])
    S.loc[j,'reaches_2_Msun']=bool(row.M_Msun>=2.)
    S.loc[j,'reaches_J0740_central_mass']=bool(row.M_Msun>=2.073)
S.to_csv(O/'stellar_summary.csv',index=False)
pd.DataFrame(profiles).to_csv(O/'profiles_L10km.csv',index=False)
conv=[]
for _,row in S.iterrows():
    e=s.EOS(row.EOS);mode=row['mode'];L=row.L_m/1000
    seq=D[(D.EOS==row.EOS)&(D['mode']==mode)&(D.L_m==row.L_m)]
    rr=seq.iloc[np.argmin(abs(seq.M_Msun-1.4))]
    for target,lr in [('peak',row.log10_rho0c_cgs),('M1.4',rr.log10_rho0c_cgs)]:
        coarse=s.star(e,float(lr),mode,L)
        fine=s.star(e,float(lr),mode,L,rtol=1e-11,ngrid=1601,nq=160,tol=2e-7)
        dm=fine['M']-coarse['M'];dr=fine['R']-coarse['R']
        if target=='M1.4':
            branch=seq[seq.log10_rho0c_cgs<row.log10_rho0c_cgs].sort_values('M_Msun').drop_duplicates('M_Msun')
            slope=float(PchipInterpolator(branch.M_Msun,branch.R_km).derivative()(coarse['M']))
            fixed=dr-slope*dm
        else:fixed=np.nan
        conv.append(dict(EOS=row.EOS,mode=mode,L_m=row.L_m,target=target,dM_Msun=dm,dR_km=dr,estimated_dR_fixed_M_km=fixed,fp_residual=fine['residual']))
    print('Convergence',row.EOS,mode,L,flush=True)
pd.DataFrame(conv).to_csv(O/'convergence.csv',index=False)
# V21 audit fix: targeted refinement of every causal endpoint and every c_s^2=c^2 core boundary.
cep=[]
for _,row in S.iterrows():
    e=s.EOS(row.EOS);mode=row['mode'];L=row.L_m/1000
    cases=[('causal_endpoint',float(row.log10_rho0c_at_M_max_causal))]
    if bool(row.causal_violation_flag): cases.append(('Mmax_core_boundary',float(row.log10_rho0c_cgs)))
    for kind,lr in cases:
        coarse=s.star(e,lr,mode,L)
        fine=s.star(e,lr,mode,L,rtol=1e-11,ngrid=1601,nq=160,tol=2e-7)
        hs=e.h_causal()
        cep.append(dict(EOS=row.EOS,mode=mode,L_m=row.L_m,target=kind,endpoint_is_Mmax=not bool(row.causal_violation_flag),
            log10_rho0c=lr,M_Msun=coarse['M'],R_km=coarse['R'],dM_Msun=fine['M']-coarse['M'],dR_km=fine['R']-coarse['R'],
            r_core_km=coarse['r_acausal_core'],dr_core_km=fine['r_acausal_core']-coarse['r_acausal_core'],
            r_core_sampled_lower_km=coarse['r_acausal_core_sampled'],
            cs2_minus_1_at_threshold_h=float(e.at(hs)[3]-1) if hs==hs else np.nan,
            central_cs2_over_c2=coarse['cs2c'],fp_residual=fine['residual']))
    print('Causal endpoint/core refinement',row.EOS,mode,L,flush=True)
CE=pd.DataFrame(cep);CE.to_csv(O/'causal_endpoint_convergence.csv',index=False)
ce_summary=dict(n_endpoint_tests=int((CE.target=='causal_endpoint').sum()),n_truncated_endpoints=int(((CE.target=='causal_endpoint')&~CE.endpoint_is_Mmax).sum()),
    n_core_boundary_tests=int((CE.target=='Mmax_core_boundary').sum()),
    max_abs_dM_endpoint_Msun=float(CE[CE.target=='causal_endpoint'].dM_Msun.abs().max()),
    max_abs_dR_endpoint_km=float(CE[CE.target=='causal_endpoint'].dR_km.abs().max()),
    max_abs_dr_core_km=float(CE[CE.target=='Mmax_core_boundary'].dr_core_km.abs().max()),
    max_abs_cs2_minus_1_at_threshold=float(CE.cs2_minus_1_at_threshold_h.abs().max()))
obs_check=s.check_observational_inputs(pd.read_csv(O/'observational_inputs.csv',keep_default_na=False))
# tabulated results and quoted fit errors, Read Table III (visually checked p20)
bench={'SLy':dict(Mtab=2.049,Mfit_percent=.02,Rtab=11.736,Rfit_percent=-.21),'APR4':dict(Mtab=2.213,Mfit_percent=-1.08,Rtab=11.428,Rfit_percent=-.90),'H4':dict(Mtab=2.032,Mfit_percent=-.85,Rtab=13.774,Rfit_percent=1.34)}
for n,b in bench.items():
    row=S[(S.EOS==n)&(S['mode']=='GR')].iloc[0]
    b['Mfit_Read']=b['Mtab']*(1+b['Mfit_percent']/100)
    b['Rfit_Read']=b['Rtab']*(1+b['Rfit_percent']/100)
    b['M_this']=row.M_Msun;b['R_this']=row.R_at_1p4_km
    b['M_relative_difference']=b['M_this']/b['Mfit_Read']-1
    b['R_relative_difference']=b['R_this']/b['Rfit_Read']-1
w=json.loads((O/'wd_summary.json').read_text())
profile_free=[]
for key in ['upper_extra_approx95_kms','quoted_mass_extra_approx95_kms']:
    exponent=np.log1p(w[key]*1000/(s.C*(1+w['vgr_kms']*1000/s.C)))
    hi=3*s.G*s.MS*1.018/(4*np.pi*s.C**2*exponent)
    profile_free.append(dict(allowance_key=key,lower_m=2*w['R_m'],upper_m=hi))
meta=dict(versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,pandas=pd.__version__,matplotlib=matplotlib.__version__),Read_GR_benchmarks=bench,top_hat_profile_diagnostics=diagnostics,exact_shell_crosschecks=geometry_cross,profile_independent_WD_band=profile_free,causal_endpoint_and_core_refinement=ce_summary,observational_inputs_sym_error_check=obs_check)
(O/'numerical_diagnostics.json').write_text(json.dumps(meta,indent=2))
print(json.dumps(meta,indent=2))
