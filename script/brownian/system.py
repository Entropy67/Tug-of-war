"""
    12/09/2019, simulate single antigen extraction
"""


################ import module ############
import numpy as np
from . import bonds as bd
from . import utilities as utl
from . import force as force
import json

###### constants ##################
PI = 3.1415926
kT = 300*1.38E-23



#### default parameter #####
import os
path = os.path.dirname(__file__)
with open(path +"/prm.json", "r") as fp:
    default_parameter = json.load(fp)
    
class System:
    """
    A class used to represent the tug-of-war system
    ...
    Attributes
    -------
    x1: float
        APC-Ag bond extension
    
    x2: float
        BCR-Ag bond extension

    Methods
    -------
    
    loadprm(prm):
        load parameters from json file
    
    setup():
        setup the system
        
    init():
        init the bond extensions. Set them to zero
    
    run():
        run simulation
    
    """
    
    def __init__(self, prm=None):
        
        
        self.prm = prm #dict, a dict that contains all parameters
        self.force_gen = force.Force_prm() #Force_prm class, a class that models different force schemes
        
        self.noise = True #boolean, true if noise is added, false if the system is deterministic
        self.output = True #boolean, true if we want output log info for debug
        self.manyRun = False #boolean, true if we want to run simulations many times else false
        
        self.frozenBCR = False 
        # boolean, true if we freeze BCR movement, 
        # this means the boundary along xb = xb^ddagger is reflective and BCR-Ag is unbreakable
    
        self.frozenAPC = False
        # boolean, true if we freeze APC movement, 
        # this means the boundary along xa=xb^ddagger is reflective and APC-Ag is unbreakable
    
        self.fixAPC = False # boolean,  true if we fix APC-Ag bond, this mean xa = constant
        
        self.fixBCR = False # boolean, true if we fix BCR-Ag bond, this means xb = constant
        
        
        self.numRun = 20 # int, number of repeated experiments
        self.numSample = 200 # int, number of runs we take to calculate the extraction chance
        
        self.setup()
        
        pass
    
    def loadprm(self, prm):
        self.prm = prm
        self.setup()
        pass
    
    def setup(self):
        self.potential = self.prm["potential"] ### read potential type, either "cusp" or "linear-cubic"
        self.record_time = self.prm["record_time"] ### step intervals to record the particle position
        
        self.tm = self.prm["tm"] ### max simulation duration
        self.dt = self.prm["dt"] ### simulation time step
        self.time_unit = self.prm["time_unit"] ### time unit
        
        self.force_gen.loadprm(self.prm) #### specify the force
        
        self.bd1, self.bd2 = bd.getBonds(
            self.prm["Eb1"], 
            self.prm["Eb2"], 
            self.prm["xb1"], 
            self.prm["xb2"], 
            self.prm["potential"],
            output=False) ### get bond, bd1 is APC-Ag and bd2 is BCR-Ag
        
        
        self.bd1.setup() ### setup APC-Ag bond
        self.bd2.setup() ### setup BCR-Ag bond
        
        if not self.noise: ### turn off the noise
            self.bd1.noiseOff() 
            self.bd2.noiseOff()
            
        self.sqrtdt = np.sqrt(self.dt) ### square root of dt, to save sim time
        self.gma_eff = self.bd1.gma*self.bd2.gma/(self.bd1.gma+self.bd2.gma) ### effective gamma = gma1 * gma2 / (gma1 + gma2)
        return
    
    
    def _updateForce(self, t):
        """
        update the force at time t
        """
        ## return force in nN
        return self.force_gen.get_f(t)*1.0E-3
    
    
    def init(self):
        """
        initialize the trajectory, set the particle to the origin
        """
        self.x1_traj, self.x2_traj, self.f_traj = [], [], []
        self.t = 0
        self.x1, self.x2 = 0, 0
        pass
    
    def run(self):
        """
        run simulation
        """
        if self.manyRun:
            return self._manyRun(self.output)
        else:
            return self._run(self.output)
        
    def _manyRun(self, output=False):
        """
        run the experiment many times
        In each experiment, we simulate numSample trajectories to calculate eta
        """
        n = self.numRun #### number of repeated experiments
        
        self.pList = [] #### list contains simulated eta values
        self.tendList = [] ### rupture time list
        self.fendList = [] #### rupture force list
        for i in range(n):
            p, stdDumn, tend, stdDumn2, fend = self._run(output=output) ### one experiment to calculate eta
            self.pList.append(p)
            self.tendList.append(tend)
            self.fendList.append(fend)
    
        self.eta = np.mean(self.pList) ### ensemble averaged eta
        self.eta_std = np.std(self.pList) 
        self.tend = np.mean(self.tendList)
        self.tend_std = np.std(self.tendList)
        self.fend = np.mean(self.fendList)
        self.fend_std=np.std(self.fendList)
        if output:
            print("\nn\teta\teta_std\ttend\ttend_std")
            print("{0:d}\t{1:.3f}\t{2:.3f}\t{3:.2f}\t{4:.2f}".format(n, self.eta, self.eta_std, self.tend, self.tend_std))
        
        return self.eta, self.eta_std, self.tend, self.tend_std, self.fend  
    
        
    def _run(self, output=True):
        """
        run one experiment (contains numSample rupture simulations) to calculate eta
        """
        self.tList = [] ### rupture time list
        self.fList = [] ### rupture force list
        
        self.fSuList = [] ### record the rupture force when the antigen is extracted
        self.fFaList = [] ### record the rupture force when the antigen is lost
        
        self.tSuList = [] ### record the rupture time when the antigen is extracted
        self.tFaList = [] ### record the rupture time if the antigen is lost
        self.x1List = [] ### record rupture x1 position
        self.x2List = [] ### record rupture x2 position
        
        count = np.zeros(3, dtype=int) ### count the number of breaking events
        #### count[0] = number of cases where no bond breaks
        #### count[1] = number of APC-Ag breaks
        #### count[2] = number of BCR-Ag breaks
        
        
        for j in range(self.numSample):
            if output:
                utl.printProgress(j, self.numSample)
            flag, p, t, f = self.run1(init=True) ### simulate one single trajectory
            if flag: ### if the complex is broken
                count[p] += 1
                self.tList.append(t)
                self.fList.append(f)
                if p>1:
                    self.tSuList.append(t)
                    self.fSuList.append(f)
                    self.x2List.append(self.x2)
                elif p==1:
                    self.tFaList.append(t)
                    self.x1List.append(self.x1)
                    self.fFaList.append(f)
            else:
                print("Warning: simulation is not finished! please addjust tm")
            
            
        self.eta = count[2]/sum(count)
        self.tend = np.mean(self.tList)
        self.tend_std = np.std(self.tList)
        self.fend = np.mean(self.fList)
        self.fend_std=np.std(self.fList)
        if output:
            print("\ncount\tf\teta\tend")
            print("{0:d}\t{1:.2f}\t{2:.3f}\t{3:.2f}".format(sum(count), self.fend, self.eta, self.tend))
        return self.eta, 0, self.tend, np.std(self.tList), self.fend
        
        
        
    def run1(self, init=True):
        '''
        simulate a single rupture trajectory
        
        @parameter:
                init: true if we want to initilize the particle position else false
        @return: 
                flag: boolean, true if the complex is broken else false
                p: int, 
                    0 : apc-ag-bcr, no bond breaks
                    1 : apc-ag bcr
                    2 : apc ag-bcr
                t: float, rupture time
                f: float, rupture force
        
        '''
        if init:
            self.init() ### set bond extensions to zero
        flag = False ## mark break or not
        step = 0
        while (self.t<self.tm):
            f = self._updateForce(self.t) ### get force, f in nN
            self._step(f)  ### two bond movement
            
            flag, p = self._breakOrNot() ### check if the complex is broken or not
            
            if step % self.record_time == 0:
                #### record the position and force
                self.x1_traj.append(self.x1)
                self.x2_traj.append(self.x2)
                self.f_traj.append(1000*f)
                
            if flag:
                break
            self.t += self.dt
            step += 1
        return flag, p, self.t, 1000*f
    
    
    
    def _step(self, f):
        """
        simulate one step of Brownian motion
        
        @arguments:
            f: float, force, in nN
        
        """
        xi1 = np.random.normal(0, self.bd1.std) ### white noise xi1
        xi2 = np.random.normal(0, self.bd2.std) ### white noise xi2
        
        fx1, fx2 = self._drift_force()  #### get drift force
        
        ######## calculate displacement according to the  Langevin equation #########
        #--------------------------------------------------------------------------------------
        dx1 = self.dt*(fx1-fx2+xi1/self.sqrtdt)/self.bd1.gma
        dx2 = self.dt*(fx2/self.gma_eff-fx1/self.bd1.gma+f/self.bd2.gma+xi2/(self.sqrtdt*self.bd2.gma)-xi1/(self.sqrtdt*self.bd1.gma))
         #--------------------------------------------------------------------------------------

        ##### update the coordinates
        if self.fixAPC and not self.fixBCR:
            self.x2 += self.dt*(fx2+f+xi2/self.sqrtdt)/self.bd2.gma
            self.x1 += 0
        
        elif self.fixBCR and not self.fixAPC:
            self.x1 += self.dt*(fx1+f+xi1/self.sqrtdt)/self.bd1.gma
            self.x2 += 0

        else:
            self.x1 += dx1
            self.x2 += dx2
        return

    
    def _drift_force(self):
        '''
        return the drift force exerting on the two bonds
        @return 
            fx1: potential force generated by APC-Ag bond
            fx2: potential drift force generated by BCR-Ag bond
        '''
        unit = 1.0E18
        fx1, fx2 = 0, 0
        if self.potential == "cusp":
            #### cusp harmonic potential
            fx1 = -self.bd1.k1*self.x1
            fx2 = -self.bd2.k1*self.x2
        elif self.potential =="linear-cubic":
            #### linear cubic potential
            fx1 = (-1.5*self.bd1.Eb/(self.bd1.xb)+1.5*self.bd1.Eb*((2*self.x1-self.bd1.xb)/self.bd1.xb)**2/(self.bd1.xb))*unit
            fx2 = (-1.5*self.bd2.Eb/(self.bd2.xb)+1.5*self.bd2.Eb*((2*self.x2-self.bd2.xb)/self.bd2.xb)**2/(self.bd2.xb))*unit
            
        return fx1, fx2
    
    def _breakOrNot(self):
        """
        check if the complex is broken or not
        @return:
            flag: boolean, true if the complex is broken else false
            p: 
                0: complex is not broken
                1: APC-Ag is broken
                2: BCR-Ag is broken
        """
        if self.bd1.broken(self.x1): #### if APC-Ag is broken
            if not self.frozenAPC:
                return True, 2
            else:
                if self.bd2.broken(self.x2):
                    return True, 1
                else:
                    self.x1 = self.bd1.reflect(self.x1) ## reflect
                    return False, 0
        
        elif self.bd2.broken(self.x2): #### if BCR-Ag is broken
            if not self.frozenBCR:
                return True, 1
            else:
                if self.bd1.broken(self.x1):
                    return True, 2
                else:
                    self.x2 = self.bd2.reflect(self.x2) ## reflect
                    return False, 0
        
        else:
            return False, 0
        
        
    def _getTraj(self):
        """
        get time,  x1, x2 trajectories
        """
        
        length = len(self.x1_traj)
        
        self.t_traj = [self.dt*i*self.record_time*self.time_unit for i in range(length)]
        return self.t_traj, self.x1_traj, self.x2_traj
    
    
    def _pot(self, x1, x2, f=0):
        """
        calculate the potential at (x1, x2) under force f, used to visulize the potential function
        """
        ## f in nN
        unit = 1.0E-18/kT
        
        if self.potential == "cusp":
            if x1<self.bd1.x1 and x2<self.bd2.x1:
                return 0.5*self.bd1.k1*x1**2*unit+0.5*self.bd2.k1*x2**2*unit-f*(x1+x2)*unit
            else:
                if x1>self.bd1.x1 and x2<self.bd2.x1 and self.frozenAPC:
                    return 0.5*self.bd1.k1*x1**2*unit+0.5*self.bd2.k1*x2**2*unit-f*(x1+x2)*unit
                elif x2>self.bd2.x1 and x1<self.bd1.x1 and self.frozenBCR:
                    return 0.5*self.bd1.k1*x1**2*unit+0.5*self.bd2.k1*x2**2*unit-f*(x1+x2)*unit
                else:
                    return None
        
        elif self.potential == "linear-cubic":
            part1 = 0.75*self.bd1.Eb*((2*x1-self.bd1.xb)/self.bd1.xb)/kT-0.25*self.bd1.Eb*((2*x1-self.bd1.xb)/self.bd1.xb)**3/kT
            part2 = 0.75*self.bd2.Eb*((2*x2-self.bd2.xb)/self.bd2.xb)/kT-0.25*self.bd2.Eb*((2*x2-self.bd2.xb)/self.bd2.xb)**3/kT
            return part1+part2-f*(x1+x2)*unit
    

  