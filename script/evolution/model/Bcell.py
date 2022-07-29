"""
B cell class

A B cell can divide, mutate, differentiate, apoptosis and extract antigens
"""
import numpy as np
from . import theory as the

TWO_MUTATION = False ## TWO mutation happen simultenously
FEEDBACK_XB = False ### whether the feedback antibody bond length determines the tether bond length


class Bcell:  
    def __init__(self, gc, E0, xb):
        self.gc = gc
        self.alive = True
        self.mature = False ## the B cell can extract antigen and compete for survival if True
        self.E = E0  ## binding energy, in kBT
        self.xb = xb ### BCR-Ag bond length, in nm
        self.eta = 0 ### average antigen extraction chance
        self.w = 0 ### fitness
        self.setup()
        return
        
    def setup(self):
        self.dE = self.gc.dE  ## binding energy change per mutation, in kBT
        self.dxb = self.gc.dxb ## bong length change per mutation, in nm
        self.pd = self.gc.pd # prob of differentiation to plasma cell, 0.02
        self.pa = self.gc.pa  ## apoptosis probability
        self.pm = self.gc.pm ### mutation probability
        self.pm_xb = self.gc.pm_xb ## xb mutation probability
        self.tau = the.lifetime_linear_cubic(self.E, self.xb) ### force-free bond lifetime
        pass
        
    def divide(self, num_children=2):
        '''
        proliferate num_children daughter cells and go die afterwards
        '''
        sons = []
        for i in range(num_children):
            soni = Bcell(self.gc, self.E, self.xb)
            soni.mutate()
            sons.append(soni)
        self.goDie()
        return sons
    
    def comp_divide(self):
        '''
        divide according to fitness self.w
        num of suns = self.w
        '''
        if self.w == 0:
            return []
        
        num_division =int( min(np.random.poisson(self.w), self.gc.w0)) ### w0 is the maximal birth rate
        sons = self.divide(num_division)
        return sons
    
    def differentiation(self, eta_mean=0):
        '''
        differentiate into plasma B cell with probability pd
        we set eta_mean = 0 by default so all B cells can differentiate
        '''
        if self.eta>eta_mean and np.random.uniform()<self.pd:
            self.goDie()
            self.gc.addPlasma(self.E, self.xb)
        
    def mutate(self):
        ''' 
        a B cell can mutate to a different affinity and/or bond length
        '''
        
        if TWO_MUTATION: #### affinity and xb change simulatenously
            if np.random.uniform()<self.pm:
                self.E += np.random.normal(0, self.dE)
                xb_tmp = self.xb + np.random.normal(0, self.dxb)

                while xb_tmp < 0.1 or xb_tmp > 10:
                    ### restrict the range of xb
                    xb_tmp = self.xb + np.random.normal(0, self.dxb)
                self.xb = xb_tmp
        else: ### affinity and xb change independently
            if np.random.uniform() < self.pm:
                 self.E += np.random.normal(0, self.dE)
                    
            if np.random.uniform() < self.pm_xb:
                ### change in bond length
                xb_tmp = self.xb + np.random.normal(0, self.dxb)

                while xb_tmp < 0.1 or xb_tmp > 10:
                    ### restrict the range of xb
                    xb_tmp = self.xb + np.random.normal(0, self.dxb)
                self.xb = xb_tmp
        return
    
        
    def goDie(self):
        self.alive = False
    
    def extract(self, Ea, xa):
        ## get extraction chance or extracted antigen number
        if FEEDBACK_XB:
            ## the tether bond length is determined by xb of the feedback antibody
            self.eta = self.gc.Ag_extract(Ea, self.E, xa,  self.xb)
        else:
            ## the tether bond length is fixed at gc.xb1
            ### fixing the tether stiffness during evolution
            self.eta = self.gc.Ag_extract(Ea, self.E, self.gc.xb1,  self.xb)
            
        self.mature=True ### mark that antigen extraction is finished
        
        if self.eta is None:
            self.eta = 0
            
        ## convert to fitness
        if self.gc.no_selection:
            self.w = self.gc.Tcell_selection(self.E-Ea)
        else:
            self.w = self.gc.Tcell_selection(self.eta)  
        return self.eta
    
    def apoptosis(self, threshold=-1):
        '''
        B cell dies if eta = 0 or with a uniform probability
        '''
        if self.eta == 0:
            self.goDie()
            return

        if self.gc.random_death: ### die randomely
            if np.random.uniform()<1-self.pa:
                self.goDie()
        else: ### only low-affinity B cells die
            if self.E < threshold:
                self.goDie()
        return
            
    