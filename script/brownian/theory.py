"""
analytical results
"""
## ----------------------------------
## ----------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import warnings

kT = 300*1.38E-23
PI = 3.1415926

    
######final version of analytical expression ###########

######### constant force ##################
def extRatio_cusp_harmonic(e10, e20, f0, xb1, xb2,Eb_min=1, m=1):
    kT = 300*1.38E-23
    k1 = 2*e10*kT*1.0E18/xb1**2
    k2 = 2*e20*kT*1.0E18/xb2**2
    f1 = np.sqrt(2.0*e10*k1*kT)
    f2 = np.sqrt(2.0*e20*k2*kT)
    f = f0*1.0E-12
    Eb1 = e10*(1-f/f1)**2
    Eb2 = e20*(1-f/f2)**2
    flag=True
    if 1-f/f1<np.sqrt(Eb_min/e10) or 1-f/f2<np.sqrt(Eb_min/e20):
        #print("warning")
        flag=False

    dE = Eb1-Eb2
    alpha = (1+m)*np.sqrt(k2/k1)*np.exp(dE)*(f2-f)/(f1-f)
    return 1.0/(1.0+alpha), flag

def extRatio_linear_cubic(e10, e20, f0, xb1, xb2, output=True, Eb_min=1, m=1):
    """
    xb1: APC-Ag bond length, nm
    xb2: BCR-Ag bond length, nm
    e10: APC-Ag affinity, kT
    e20: BCR-Ag affinity, kT
    f: force, pN
    """
    kT = 4.012
    ## Eq.30, Eq.31, eq.32 in si
    fa = 3*e10*kT/(2*xb1) ## pN
    fb = 3*e20*kT/(2*xb2) ## pN
    if f0>fa or f0>fb:
        return 0, False
    
    ka = 4*(fa/xb1)*np.sqrt(1-f0/fa)*0.001 ## nN/nm
    kb = 4*(fb/xb2)*np.sqrt(1-f0/fb)*0.001 ## nN/nm
    
    ga=gb=1
    
    ## barrier height
    Eb1 = e10*(1-f0/fa)**(3/2)
    Eb2 = e20*(1-f0/fb)**(3/2)
    
    ## prefactor
    tau0_a = (ga*gb/(2*ka*kb))*((kb-ka)/ga + kb/gb + np.sqrt(((kb-ka)/ga + kb/gb)**2+4*ka*kb/(ga*gb)))
    tau0_b = (ga*gb/(2*ka*kb))*((-kb+ka)/ga - kb/gb + np.sqrt(((-kb+ka)/ga - kb/gb)**2+4*ka*kb/(ga*gb)))
    if Eb1<Eb_min or Eb2<Eb_min:
        return 1/(1+tau0_a*np.exp(Eb1-Eb2)/tau0_b), False
    return 1/(1+tau0_a*np.exp(Eb1-Eb2)/tau0_b), True


def bell(e10, e20, f0, xb1, xb2, m=1, func= extRatio_cusp_harmonic):
    eta0, _ = func(e10, e20, 0, xb1, xb2, m)
    tmp = 1/eta0-1
    return 1/(1+tmp*np.exp(f0*(xb2-xb1)/4.012)), True


def tau_cusp_harmonic(e10, e20, f0, xb1, xb2):
    kT = 300*1.38E-23
    k1 = 2*e10*kT*1.0E18/xb1**2
    k2 = 2*e20*kT*1.0E18/xb2**2
    f1 = np.sqrt(2.0*e10*k1*kT)
    f2 = np.sqrt(2.0*e20*k2*kT)
    f = f0*1.0E-12
    Eb1 = e10*(1-f/f1)**2
    Eb2 = e20*(1-f/f2)**2
    flag=True
    ga=gb=1
    gab = ga * gb / (ga + gb)
    if 1-f/f1<np.sqrt(Eb_min/e10) or 1-f/f2<np.sqrt(Eb_min/e20):
        #print("warning")
        flag=False
    
    tau_a = np.sqrt(PI * kT / k1) * ga * np.exp(Eb1) / (k1 * xb1 - f)
    
    tau_b = np.sqrt(PI * kT / k2) * gab * np.exp(Eb2) / (k2 * xb2 - f)
    
    return tau_a, tau_b, tau_a * tau_b / (tau_a + tau_b)


def tau_linear_cubic(xb1, xb2, e10, e20, f0, output=True, Eb_min=3, m=None):
    """
    xb1: APC-Ag bond length, nm
    xb2: BCR-Ag bond length, nm
    e10: APC-Ag affinity, kT
    e20: BCR-Ag affinity, kT
    f: force, pN
    """
    kT = 4.012
    ## Eq.30, Eq.31, eq.32 in si
    fa = 3*e10*kT/(2*xb1) ## pN
    fb = 3*e20*kT/(2*xb2) ## pN
    if f0>fa or f0>fb:
        return np.nan, np.nan, np.nan
        
    
    ka = 4*(fa/xb1)*np.sqrt(1-f0/fa)*0.001 ## nN/nm
    kb = 4*(fb/xb2)*np.sqrt(1-f0/fb)*0.001 ## nN/nm
    
    ga=gb=1
    
    ## barrier height
    Eb1 = e10*(1-f0/fa)**(3/2)
    Eb2 = e20*(1-f0/fb)**(3/2)
    
    ## prefactor
    tau0_a = 2*PI*(ga*gb/(2*ka*kb))*((kb-ka)/ga + kb/gb + np.sqrt(((kb-ka)/ga + kb/gb)**2+4*ka*kb/(ga*gb)))
    tau0_b = 2*PI*(ga*gb/(2*ka*kb))*((-kb+ka)/ga - kb/gb + np.sqrt(((-kb+ka)/ga - kb/gb)**2+4*ka*kb/(ga*gb)))
    if Eb1<Eb_min:
        tau_a=np.nan
    elif f0<=fa:
        tau_a = tau0_a*np.exp(Eb1)
    if Eb2<Eb_min:
        tau_b = np.nan
    elif f0<=fb:
        tau_b = tau0_b*np.exp(Eb2)
    return tau_a, tau_b, tau_a * tau_b / (tau_a + tau_b)

def tau_linear_cubic_single(e0, xb, f0):
    kT = 4.012
    ## Eq.30, Eq.31, eq.32 in si
    fr = 3*e0*kT/(2*xb) ## pN

    if f0>fr:
        return np.nan
        
    
    ka = 4*(fr/xb)*np.sqrt(1-f0/fr)*0.001 ## nN/nm
    g=1
    
    ## barrier height
    Eb = e0*(1-f0/fr)**(3/2)

    return g*2*3.1416*np.exp(Eb)/ka



###################################
#### for time dependend force calculation
def eta0(k1,k2,E1,E2, f0=0, m=1):
    ## extraction without external force
    f1 = np.sqrt(2.0*E1*k1*kT)
    f2 = np.sqrt(2.0*E2*k2*kT)
    f = f0*1.0E-12
    if f>f1*(1.0-1.0/np.sqrt(E1)) or f>f2*(1.0-1.0/np.sqrt(E2)):
        approxWarn()
    dE = E1*(1-f/f1)**2-E2*(1-f/f2)**2
    tmp = (1+m)*np.sqrt(k2/k1)*np.exp(dE)*(f2-f)/(f1-f)
    return 1.0/(1.0+tmp)


def tauCA(f0, E1, k1, m=1, cutoff=False):
    ## f0 : pN
    ## E1 : in kT
    Mg = m
    f1 = np.sqrt(2*k1*E1*kT)
    f = f0*1.0E-12
    
    
    part1 = 2*Mg*np.sqrt(PI)/k1
    part2 = np.exp(E1*(1-f/f1)**2)/(np.sqrt(E1)*(1-f/f1))
        
    if cutoff and (f>f1 or E1*(1-f/f1)**2<0.5):
        return None
    return part1*part2

def tau_cusp(f0, E, xb, m=1):
    '''
    tau0 = 2\pi\gamma/\sqrt{k^2*\pi Eb} *exp(Eb)
    '''
    k = 2*E*kT*1.0E18/xb**2
    tau0 = 2*PI*m/np.sqrt(k**2*PI*E)
    #tau0 = 2*PI/np.sqrt(k)
    fm = np.sqrt(2*k*E*kT)
    f = f0*1.0E-12
    if f<fm:
        part1 = 1-f/fm
        part2 = E*(1-f/fm)**2
        return tau0*np.exp(part2)/part1
    
    return 1


def tauAB(f0, E2, k2, m=0.5, cutoff=False):
    f2 = np.sqrt(2*k2*E2*kT)
    f = f0*1.0E-12
    

    part1 = 2*m*np.sqrt(PI)/(k2*np.sqrt(E2))
    part2 = np.exp(E2*(1-f/f2)**2)/(1-f/f2)
    
    if cutoff and (f>f2 or E2*(1-f/f2)**2<0.5):
        return None
    return part1*part2
    
    
    
def Sab2(k2, E2, fi0, ft0, r0, m=1.0, pot="cusp"):
    ## survival probability for Ag-BCR bond
    ## force in pN
    Mg = 1.0   ## M*gma 
    ## bond parameters
    f2 = np.sqrt(2.0*E2*k2*kT)
    
    fi = fi0*1.0E-12
    ft = ft0*1.0E-12
    r = r0*1.0E-12
    if fi0==ft0 or r0==0:
        return 1.0
    #part1 = -(1+m)*np.sqrt(2*k2**3*kT/PI)/(4*r*Mg)
    part1 = -(1+m)*np.sqrt(2*(k2**3)*kT/PI)/(4*r*Mg)
    part2 = np.exp(-E2*(1-ft/f2)**2)-np.exp(-E2*(1-fi/f2)**2)
    return np.exp(part1*part2)

def Sca2(k1, E1, fi0, ft0, r0, m=1.0, pot="cusp"):
    ## survical probability for CR-Ag bond
    ## force in pN
    Mg = 1.0   ## M*gma 
    ## bond parameters
    f1 = np.sqrt(2.0*E1*k1*kT)
    
    fi = fi0*1.0E-12
    ft = ft0*1.0E-12
    r = r0*1.0E-12
    if fi0==ft0 or r0==0:
        return 1.0
    #part1 = -np.sqrt(2*k1**3*kT/PI)/(4*r*Mg)
    part1 = -np.sqrt(2*(k1**3)*kT/(PI))/(4*r*Mg)
    part2 = np.exp(-E1*(1-ft/f1)**2)-np.exp(-E1*(1-fi/f1)**2)
    return np.exp(part1*part2)


def extRatioFab(E1, E2, xb, zb, r=0, f0=0, fm=30, m=1):
    k1 = 2*E1*kT*1.0E18/xb**2
    k2 = 2*E2*kT*1.0E18/zb**2
    return extRatioF2(E1, E2, k1, k2, r, f0, fm, m), True
    
def extRatioF2(E1, E2, k1, k2, r=0, f0=0, fm=30, m=1.0):
    ## return extraction ratio under increasing force
    ### force in pN
    combine = lambda f: Sab2(k2, E2, f0, f, r, m)*Sca2(k1, E1, f0, f, r, m)/(r*tauCA(f, E1, k1, m))
    if f0 == fm or r==0:
        ret1 = 0
    else:
        ret1 = integrate.quad(combine, f0, fm)[0]
    ret2 = eta0(k1,k2,E1,E2,fm, m)*Sab2(k2, E2, f0, fm, r, m)*Sca2(k1, E1, f0, fm, r, m)
    #print(ret1, ret2)
    return ret1+ret2

def extRatioP(E1, E2, k1, k2, tH, tL, fH, fL, m=1.0):
    print("using periodic force: tH={0:.2f}, tL={1:.2f}, fH={2:.4f}, fL={3:.4f}".format(tH, tL, fH, fL))
    tau1H = tauCA(fH, E1, k1)
    tau1L = tauCA(fL, E1, k1)
    tau2H = tauAB(fH, E2, k2)
    tau2L = tauAB(fL, E2, k2)
    tmp1 = tH*(1.0/tau1H+1.0/tau2H)
    tmp2 = tL*(1.0/tau1L+1.0/tau2L)
    print("t1H={0:.2f}, t2H={1:.2f}".format(tau1H, tau2H))
    print("t1L={0:.2f}, t2L={1:.2f}".format(tau1L, tau2L))
    alpha = (1-np.exp(-tmp1))/(1-np.exp(-tmp1-tmp2))
    etaH = 1.0/(1.0+tau1H/tau2H)
    etaL = 1.0/(1.0+tau1L/tau2L)
    print("etaH={0:.4f}, etaL={1:.4f}, alpha={2:.4f}".format(etaH, etaL, alpha))
    return etaH*alpha+etaL*(1-alpha)

def ruptureForce(E1, E2, k1, k2, r=0, f0=0, fm=80, m=1.0):
    ## force in pN
    df = 0.001
    combine = lambda f: f*(-Sca2(k1, E1, f0, f+df,r,m)*Sab2(k2, E2, f0, f+df, r, m)+Sca2(k1, E1, f0, f,r,m)*Sab2(k2, E2, f0, f, r, m))/df
    if f0 == fm or r ==0:
        ret1 = 0
    else:
        ret1 = integrate.quad(combine, f0, fm)[0]
    ret2 = fm*Sab2(k2, E2, f0, fm, r, m)*Sca2(k1, E1, f0, fm, r, m)
    return ret1 + ret2

def ruptureTime(E1, E2, k1, k2, r=0.01, f0=0, fm=80, m=1.0):
    dt = 0.001
    combine = lambda t: t*(-Sca2(k1, E1, f0, r*(t+dt),r,m)*Sab2(k2, E2, f0, r*(t+dt), r, m)+Sca2(k1, E1, f0, r*t,r,m)*Sab2(k2, E2, f0, r*t, r, m))/dt
    if f0 == fm or r ==0:
        ret1 = 0
    else:
        ret1 = integrate.quad(combine, f0/r, fm/r)[0]
    ret2 = fm*Sab2(k2, E2, f0, fm, r, m)*Sca2(k1, E1, f0, fm, r, m)
    return ret1 + ret2



def find_E50(E1, xb1, xb2, f0, eta0=0.5, tol=0.01, func=extRatio_cusp_harmonic):
    Emin = 10
    Emax = 50
    flag=True
    while Emax-Emin>tol:
        Emid = (Emax+Emin)/2
        eta, flag = func(E1, Emid, f0, xb1, xb2, m=1)
        #print(eta)
        if abs(eta-eta0)<0.001:
            return Emid, flag
        elif eta>eta0:
            Emax = Emid
        else:
            Emin = Emid
        #print(Emax, ", ", Emin, ", ", Emid, ", ", eta)
    #print(eta)
    return Emid, flag

