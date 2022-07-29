"""
Plasma cell
"""
from . import theory as the

class PlasmaCell:
    def __init__(self, GC, E0, xb, t0=0):
        """
        a plasma B cell can generate antibodies which enter the feedback antibody pool and act as part of the tether bond
        the lifetime is T
        
        args:
            GC: germinal center class
            E0: BCR-Ag binding energy, barrier height
            xb: BCR-Ag bond length
            t0: the time point when a plasma cell is born
        """
        self.gc = GC
        self.alive = True
        self.E = E0
        self.xb = xb
        self.t0 = t0
        
        self.tau = the.lifetime_linear_cubic(self.E, self.xb) ### force-free bond lifetime
        self.T = GC.Tp   #40
        
    
    def goDie(self):
        """
        the plasma cell go die
        """
        self.alive = False

