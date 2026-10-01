import numpy as np

def _vl(dm, dp):
    out = np.zeros_like(dm, dtype=float)
    m = dm * dp > 0.0
    out[m] = 2.0 * dm[m] * dp[m] / (dm[m] + dp[m])
    return out

def _slopes(phi, phi_in):
    n = len(phi)
    g = np.empty(n + 4, dtype=float)
    g[2:-2] = phi
    g[0:2] = phi_in
    g[-2:] = phi[-1]
    return _vl(g[2:n+2] - g[1:n+1], g[3:n+3] - g[2:n+2])

def _uface(u):
    uf = np.empty(len(u)+1, dtype=float)
    uf[0], uf[-1] = u[0], u[-1]
    if len(u) > 1:
        uf[1:-1] = 0.5*(u[:-1] + u[1:])
    return uf

def _cell_T(a, b, dx):
    if a <= 0 or b <= 0:
        raise ValueError("SLFV reference operator requires u > 0.")
    if abs(b-a) <= 1e-14*max(1., abs(a), abs(b)):
        return dx/a
    return dx*np.log(b/a)/(b-a)

def _travel_faces(uf, dx):
    return np.r_[0., np.cumsum([_cell_T(uf[i], uf[i+1], dx)
                                 for i in range(len(uf)-1)])]

def _invert(T, tf, uf, dx):
    if T <= 0: return 0.
    if T >= tf[-1]: return dx*(len(uf)-1)
    i = max(0, min(len(uf)-2, int(np.searchsorted(tf, T, side="right")-1)))
    tau = T-tf[i]
    a,b = uf[i],uf[i+1]
    if abs(b-a) <= 1e-14*max(1.,abs(a),abs(b)):
        x=a*tau
    else:
        s=(b-a)/dx
        x=a*np.expm1(s*tau)/s
    return i*dx+np.clip(x,0.,dx)

def _backtrace(u, dx, dt):
    uf=_uface(u); tf=_travel_faces(uf,dx)
    dep=np.empty_like(tf)
    for j,T in enumerate(tf):
        td=T-dt
        dep[j]=_invert(td,tf,uf,dx) if td >= 0 else td*uf[0]
    return dep

def _primitive(phi,slope,dx,x):
    n=len(phi); L=n*dx
    x=float(np.clip(x,0.,L))
    if x<=0: return 0.
    if x>=L: return float(np.sum(phi)*dx)
    i=min(int(x//dx),n-1); z=x-i*dx
    return float(np.sum(phi[:i])*dx + phi[i]*z
                 + 0.5*slope[i]/dx*((z-0.5*dx)**2-(0.5*dx)**2))

def _mass(phi,slope,dx,a,b,phi_in):
    if b<=a: return 0.
    if b<=0: return phi_in*(b-a)
    if a>=0:
        return _primitive(phi,slope,dx,b)-_primitive(phi,slope,dx,a)
    return phi_in*(-a)+_primitive(phi,slope,dx,b)

def slfv_advection_phi(phi,u,dx,dt,phi_in=0.):
    """Conservative finite-volume semi-Lagrangian transport for u>0."""
    phi=np.asarray(phi,dtype=float); u=np.asarray(u,dtype=float)
    if phi.ndim != 1 or u.ndim != 1 or len(phi)!=len(u):
        raise ValueError("phi and u must be 1-D arrays of equal length")
    if dx<=0 or dt<0: raise ValueError("Require dx>0 and dt>=0")
    if len(phi)==0 or dt==0: return phi.copy()
    dep=_backtrace(u,dx,dt)
    slope=_slopes(phi,float(phi_in))
    out=np.empty_like(phi)
    for i in range(len(phi)):
        out[i]=_mass(phi,slope,dx,dep[i],dep[i+1],float(phi_in))/dx
    if not np.all(np.isfinite(out)):
        raise FloatingPointError("Non-finite SLFV result")
    return out
