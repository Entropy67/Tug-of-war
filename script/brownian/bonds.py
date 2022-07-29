"""
class to model different Bonds
"""

## ==================================
## ----------------------------------
### import modules
import numpy as np
import matplotlib.pyplot as plt
import warnings

## ----------------------------------
### constants
PI = 3.1415926
kT = 300*1.38E-23


def approxWarn():
    warnings.warn("--force too large, approximation fails--", DeprecationWarning)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    approxWarn()

    
def getBonds(e1=10, e2=10, x1=3.0, x2=2.0, pot="cusp", output=True):

    apc = Bond()
    apc.x1 =x1 ## APC-Ag bond length, nm
    apc.e1= e1 ## APC-Ag affinity, kT
    apc.potential=pot ### potential type
    apc.setup() ### setup
    if output:
        apc.info()


    bcr = Bond()
    bcr.x1 = x2 ## BCR-Ag bond length, nm
    bcr.e1 = e2 ## BCR-Ag affinity, kT
    bcr.potential=pot ### potential type
    bcr.setup() ### setup
    if output:
        bcr.info()
    return apc, bcr


    

class Bond:
    '''
    representation of a bond
    includint its potential landscape information
    '''
    
    def __init__(self):
        
        self.e1 = 6 ### unbinding potential well depth in kT
        self.e2 = 2  ### binding potential barrier in kT
        
        ### we consider the bond length is fixed when changing the affinity
        self.x1 = 2.0 ### unbinding potential length in nm
        self.x2 = 2.0 ### binding potential length in nm
        
        self.f0 = 0  ## external force applied, in pN
        
        self.m = 1.0 ## mass
        self.gma = 1.0 ## viscocity constant

        self.x = 0 ### inital position
        
        self.potential = "cusp" ### potential: can be "cusp", "linear-cubic"
        self.setup()
        return
    
    def setup(self):
        self.mg = self.m*self.gma ### mass * gamma
        self.std = np.sqrt(2*self.gma*kT*1.0E18) ### noise
        
        self.f = self.f0*1.0E-12 ### in N
        self.t_off, self.t_on, self.koff, self.kon = 0, 0, 0, 0
        if self.potential == "cusp":
            self.de = self.e1-self.e2
            self.E1 = self.e1*kT  ## unbinding energy
            self.E2 = self.e2*kT  ## binding energy
            self.dE = self.E1-self.E2
            
            ## compute stiffness k1 and k2 using x1 and x2
            self.k1 = 2*self.E1*1.0E18/self.x1**2  ## unit: nN/nm
            self.k2 = 2*self.E2*1.0E18/self.x2**2  ## unit: nN/nm

            self.xb = self.x1 + self.x2/2
            self.xd = self.x1

            self.x_unbound = self.x1+self.x2

            self.t_off = kramer(self.e1, self.k1, self.f, self.gma)
            self.t_on = kramer(self.e2, self.k2, 0, self.gma)
            self.koff = 1.0/self.t_off
            self.kon = 1.0/self.t_on
            
            self.x_break = self.x1
            
        elif self.potential == "linear-cubic":
            self.xb = self.x1
            self.Eb = self.e1*kT
            
            self.k1 = 6*self.Eb*1.0E18/self.xb**2
            
            self.f_drift = 3*self.Eb/(2*self.xb*1.0E-9)
            
            self.x_break = self.x1
            #self.x_break = 0.5*self.xb*(1+np.sqrt(1-self.f/self.f_drift))
            
        else:
            raise Exception("No such potential!")
            
        return
    
    def noiseOff(self):
        '''
        turn noise off
        '''
        self.std = 0
        return
    
    
    def broken(self, x):
        '''
        @return true if the bond breaks otherwise false
        '''
        if x<self.x_break:
            return False
        else:
            return True
        
    def reflect(self, x):
        '''
        implement reflective boundary condition
        @return: new coordinate after reflection
        '''
        if x<=self.x1:
            return x
        else:
            return 2*self.x1-x

    def info(self):
        print("bond length: x1={0:.4f}, x2={1:.4f}".format(self.x1, self.x2))
        print("bond stiffness: k1={0:.4f}, k2={1:.4f}".format(self.k1, self.k2))
        print("energy barrier: e1={0:.1f}, e2={1:.1f}, de={2:.1f}".format(self.e1, self.e2, self.de))
        print("wait time : t_on={0:.3f}, t_off={1:.3f}".format(self.t_on, self.t_off))
        print("reaction rate: k_on={0:.4e}, k_off={1:.3e}".format(self.kon, self.koff))
        pass
    
    
        
def kramer(e, k, f, m=1.0):
    """
    bond lifetime, used for cusp-harmonic potential
    """
    ## f in N
    tau0 = 2*m*np.sqrt(PI)/(k*np.sqrt(e))
    fs = np.sqrt(2.0*e*k*kT)
    dE=e*(1-f/fs)**2
    if 1-f/fs<0:
        approxWarn()
        #raise Exception("====== Error: energy barrier vanishes! ======")
    return tau0*np.exp(dE)/(1-f/fs)

        