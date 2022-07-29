"""
class to represent different force schemes
"""

class Force_prm:
    """
    A class used to represent force schedule. It includes various force schemes
    ...

    Attributes
    ----------
    scheme: str
        the force scheme, can be one of the following:
            c: constant force
            r: linear ramping force
            p: periodic pulse
            nr: nonlinear ramping force
            sr: sigmoid force
        
    f0 : float
        constant force magnitude F(t) = f0
    r : float
        loading rate of a ramping force, F(t) = r * t
    beta : float
        nonlinearity of a non-linear dynamical force,  F(t) = r * t^beta
    
    tL, tH : float
        low-force duration, high-force duration
    fL, fH : float
        low-force magnitude, high-force magnitude
        
    tS: float
        switching time for sigmodal force, F(t) = f0 * t / ( t + tS)
        
    f: float
        instantenous force
        
    Methods
    -------
    loadprm(prmdict)
        load parameters from json file
        
    get_f(t)
        get instantenous force at time t
    """
    
    def __init__(self, scheme="c"):
        """
        initilize attributes
        """
        self.scheme=scheme
        
        self.f0 = 0
        
        ### ramping, r
        self.r = 0
        
        ### nonlinear ramping, nr
        self.beta = 1.0
        
        
        ## pulse, p
        self.tL = 100
        self.tH = 100
        self.fL = 0
        self.fH = 0
        
        ### sigmoid ramping, sr
        self.tS = 100
        
        ##
        self.f = 0
        
    def loadprm(self, prmdict):
        """
        load parameters from json file
        @arguments:
            prmdict: dict, contain parameters
        """
        self.scheme=prmdict["scheme"]
        self.r = prmdict["r"]
        self.f0 = prmdict["f0"]
        self.beta = prmdict["beta"]
        self.tL = prmdict["tL"]
        self.tH = prmdict["tH"]
        self.fL = prmdict["fL"]
        self.fH = prmdict["fH"]
        self.tS = prmdict["tS"]
        pass
        
    def get_f(self, t):
        """
        get instantenous force at time t
        @arguments:
            t: float, time
        @return:
            force in pN
        """
        if self.scheme=="r":
            ### ramping force
            self.f = self.r*t
            
        elif self.scheme=="c":
            ## const force
            self.f = self.f0
            
        elif self.scheme=="p":
            ### periodic pulse
            t_eff = int(t)%int(self.tL + self.tH)
            if t_eff > self.tH:
                self.f = self.fL
            else:
                self.f = self.fH
                
        elif self.scheme=="nr":
            ## nonlinear ramping
            self.f = self.r*(t**self.beta)
            
        elif self.scheme=="sr":
            ## sigmoid ramping
            self.f = self.f0*t/(t+self.tS)
                
        return self.f