"""
theoretical calculation of antigen extraction chance

"""

import numpy as np
from scipy import integrate
import warnings

kT = 300*1.38E-23
PI = 3.1415926

"""
extraction chance under a ramping force for the cusp-harmonic potential
"""


    
def approxWarn():
    warnings.warn("--force too large, approximation fails--", DeprecationWarning)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    approxWarn()

def getStiffness(xb, Eb):
    return 2*Eb*kT*1.0E18 / xb**2

def extRatio_rampongF_cuspHarmonic(E1, E2, xa, xb, r=0, f0=0, fm=30, m=1):
    k1 = 2*E1*kT*1.0E18/xa**2
    k2 = 2*E2*kT*1.0E18/xb**2

    return extRatio_rampingF_cuspHarmonic_withStiffness(E1, E2, k1, k2, r, f0, fm, m)

def meanRuptureForce(E1, E2, xa, xb, r=0, f0=0, fm=30, m=1):
    k1 = 2*E1*kT*1.0E18/xa**2
    k2 = 2*E2*kT*1.0E18/xb**2
    return meanRuptureForce_withStiffness(E1, E2, k1, k2, r, f0, fm, m)
    
def meanRuptureForce_withStiffness(E1, E2, k1, k2, r=0, f0=0, fm=80, m=1.0):
    ## force in pN
    df = 0.001
    combine = lambda f: f*(-Sca2(k1, E1, f0, f+df,r,m)*Sab2(k2, E2, f0, f+df, r, m)+Sca2(k1, E1, f0, f,r,m)*Sab2(k2, E2, f0, f, r, m))/df
    if f0 == fm or r ==0:
        ret1 = 0
    else:
        ret1 = integrate.quad(combine, f0, fm)[0]
    ret2 = fm*Sab2(k2, E2, f0, fm, r, m)*Sca2(k1, E1, f0, fm, r, m)
    return ret1 + ret2
    

def extRatio_rampingF_cuspHarmonic_withStiffness(E1, E2, k1, k2, r=0, f0=0, fm=30, m=1.0):
    '''
    calculate the extraction chance under a ramping force, in the cusp-harmonic potential
    @param
        E1, APC-Ag affinity
        E2, BCR-Ag affinity
        k1: APC-Ag stiffness, 2*E1*kT*1.0E18/xa**2, xa in nm, E1 in kT
        k2: BCR-Ag stiffness, 2*E2*kT*1.0E18/xb**2, xb in nm, E2 in kT
        r: ramping rate, pN/s
        f0: initial force, pN
        fm: cut-off force, pN
        m: ratio of damping constant
    @return:
        extraction chance
    '''
    
    
    combine = lambda f: Sab2(k2, E2, f0, f, r, m)*Sca2(k1, E1, f0, f, r, m)/(r*tauCA(f, E1, k1, m))
    if f0 == fm or r==0:
        ret1 = 0
    else:
        ret1 = integrate.quad(combine, f0, fm)[0]
    ret2 = eta0(k1,k2,E1,E2,fm, m)*Sab2(k2, E2, f0, fm, r, m)*Sca2(k1, E1, f0, fm, r, m)
    return ret1+ret2



def Sab2(k2, E2, fi0, ft0, r0, m=1.0, pot="cusp"):
    '''
    return the survival probabily for the BCR-Ag bond when force reaches ft
    @param:
        k2: BCR-Ag stiffness
        E2: BCR-Ag affinity
        fi0: initial force, in pN
        ft0: force at time t, in pN
        m: ratio of damping constant
        pot: potential function
    @return:
        survival probability at t
    '''
    
    Mg = 1.0   ## M*gma 
    ## bond parameters
    f2 = np.sqrt(2.0*E2*k2*kT)
    #print("force f2=", f2)
    fi = fi0*1.0E-12
    ft = ft0*1.0E-12
    r = r0*1.0E-12
    
    if fi0==ft0 or r0==0:
        return 1.0
    part1 = -(1+m)*np.sqrt(2*(k2**3)*kT/PI)/(4*r*Mg)
    part2 = np.exp(-E2*(1-ft/f2)**2)-np.exp(-E2*(1-fi/f2)**2)
    return np.exp(part1*part2)


def Sca2(k1, E1, fi0, ft0, r0, m=1.0, pot="cusp"):
    '''
    return the survival probabily for the APC-Ag bond when force reaches ft
    @param:
        k1: APC-Ag stiffness
        E1: APC-Ag affinity
        fi0: initial force
        ft0: force at time t
        m: ratio of damping constant
        pot: potential function
    @return:
        survival probability at t
    '''
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


def tauCA(f0, E1, k1, m=1, cutoff=False):
    '''
    return the lifetime of APC-Ag bond under a constant force f0
    @param: 
        f0: force in pN
        E1: APC-Ag bond affinity
        k1: APC-Ag bond stiffness
    @return
        the lifetime of APC-Ag bond
    '''
    Mg = m
    f1 = np.sqrt(2*k1*E1*kT)
    f = f0*1.0E-12
    
    
    part1 = 2*Mg*np.sqrt(PI)/k1
    part2 = np.exp(E1*(1-f/f1)**2)/(np.sqrt(E1)*(1-f/f1))
        
    if cutoff and (f>f1 or E1*(1-f/f1)**2<0.5):
        return None
    return part1*part2

def eta0(k1,k2,E1,E2, f0=0, m=1):
    '''
    return the extraction chance under a constant force
    @param:
        k1: APC-Ag bond stiffness
        k2: BCR-Ag bond stiffness
        E1: APC-Ag bond affinity
        E2: BCR-Ag bond affinity
        f0: force
        m: ratio of damping constant
    @return:
        extraction chance
    '''
    f1 = np.sqrt(2.0*E1*k1*kT)
    f2 = np.sqrt(2.0*E2*k2*kT)
    f = f0*1.0E-12
    if f>f1*(1.0-1.0/np.sqrt(E1)) or f>f2*(1.0-1.0/np.sqrt(E2)):
        approxWarn()
    dE = E1*(1-f/f1)**2-E2*(1-f/f2)**2
    tmp = (1+m)*np.sqrt(k2/k1)*np.exp(dE)*(f2-f)/(f1-f)
    return 1.0/(1.0+tmp)



def lifetime_linear_cubic(eb, xb, f0=0, output=False, Eb_min=None, pre_factor=True):
    '''
    lifetime of a single bond with affinity Eb and bond length xb
    based on Kramers' theory
    '''
    kT = 4.012
    const = 2*3.1415
    if f0==0:
        fb = 3*eb*kT/(2*xb)
        k = 4*(fb/xb)*0.001 ## nN/nm
        if not pre_factor:
            const = 1
            k=1
        return const*np.exp(eb)/k
    else:
        fb = 3*eb*kT/(2*xb)
        if Eb_min is not None and f0>fb:
            return np.nan
        Ef = eb*(1-f0/fb)**(3/2)
        
        if output:
            print("     fb=", fb, ", f=", f0)
            print("     Ef=", Ef, "kT")
        if Eb_min is not None and Ef<Eb_min:
            return np.nan
        k = 4*(fb/xb)*0.001 * np.sqrt(1 - f0/fb)
        if not pre_factor:
            const = 1
            k=1
        return const*np.exp(Ef) / k
    
    

def extRatio_linear_cubic(xb1, xb2, e10, e20, f0, output=True, Eb_min=3, pre_factor=True):
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
        return np.nan, False, -1
    
    ka = 4*(fa/xb1)*np.sqrt(1-f0/fa)*0.001 ## nN/nm
    kb = 4*(fb/xb2)*np.sqrt(1-f0/fb)*0.001 ## nN/nm
    
    ga=gb=1
    
    ## barrier height
    Eb1 = e10*(1-f0/fa)**(3/2)
    Eb2 = e20*(1-f0/fb)**(3/2)
    
    ## prefactor
    tau0_a = (ga*gb/(2*ka*kb))*((kb-ka)/ga + kb/gb + np.sqrt(((kb-ka)/ga + kb/gb)**2+4*ka*kb/(ga*gb)))
    tau0_b = (ga*gb/(2*ka*kb))*((-kb+ka)/ga - kb/gb + np.sqrt(((-kb+ka)/ga - kb/gb)**2+4*ka*kb/(ga*gb)))
    
    if not pre_factor:
        tau0_a = tau0_b = 1
    if Eb1<Eb_min or Eb2<Eb_min:
        return 1/(1+tau0_a*np.exp(Eb1-Eb2)/tau0_b), False, min(Eb1, Eb2)
    return 1/(1+tau0_a*np.exp(Eb1-Eb2)/tau0_b), True, min(Eb1, Eb2)


def extRatio(xb1, xb2, e10, e20, f0, m0=1.0, output=True):
    ## parameters
    m = m0
    Mg = 1.0   ## M*gma
    T = 300   # temperature
    kB = 1.38E-23   # Boltzmann constant
    
    ## bond parameters
    E1 = e10*kB*T  ## CR-Ag bond binding energy, in kT
    
    k10 = 2*E1/(xb1*1.0E-9)**2
    k1 = k10       ## CR-Ag bond stiffness, nN/nm
    xb = np.sqrt(2*E1/k1)*1.0E9   ## CR-Ag bond maximum deformation
    
    
    E2 = e20*kB*T  ## BCR-Ag bond binding energy, in kT
    k20 = 2*E2/(xb2*1.0E-9)**2
    k2 = k20       ## BCR-Ag bond stiffness, nN/nm
    zb = np.sqrt(2*E2/k2)*1.0E9  ## BCR-Ag bond range, nN/nm
    
    
    f = f0*1.0E-12   ## convert F from pN to N
    
    
    f1 = np.sqrt(2.0*E1*k1)
    f2 = np.sqrt(2.0*E2*k2)
    
    Eb1 = E1*(1-f/f1)**2   ## energy barrier
    Eb2 = E2*(1-f/f2)**2   ## 

    p = np.exp((Eb1-Eb2)/(kB*T))
    q = (1+m)*np.sqrt(k2/k1)*(f2-f)/(f1-f)

    if f<f1*(1.0-2.0/np.sqrt(E1/(kB*T))) and f<f2*(1.0-2.0/np.sqrt(E2/(kB*T))):
        return 1/(1+q*p), True
    elif f>1*f1 or f>1*f2:
        #write("---- Error: energy barrier vanishes ----")
        if output:
            print("*", end="")
        return 1/(1+q*p), False
    else:
        #write("---- Warrining: energy barrier too small ----")
        if output:
            print("!", end="")
        return 1/(1+q*p), False
    return None, False



