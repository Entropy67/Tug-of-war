"""
Implement the GC class
"""

from .Bcell import *
from .plasmaCell import *
from .prm import *
from .utilities import *
from . import singleBond as sb
from . import theory as theory

import json
import numpy as np
import matplotlib.pyplot as plt

from time import gmtime, strftime
from collections import defaultdict
from scipy import integrate
import warnings
import os
from mpmath import mp
import heapq
from termcolor import colored
from copy import deepcopy
from datetime import datetime

INFINITE_PLASMA = True ### infinite plasma cell capacity



class GC:
    def __init__(self, f0=10.0, prm=prm_list):
            ## current number of alived B cells
        self.f = f0    ## pulling force applied by B cells
        self.prm = deepcopy(prm)
        
        self.qty_name = ["e",        ### GC B cell affinity mean
                         "eStd",     ### GC B cell affinity std
                         "eVar",     ### GC B cell affinity variance
                         "eMp",      ### most prob affinity
                         "tau",      ### GC B cell force-free bond lifetime
                         "tauVar",   ### lifetime variance
                         "tauF",     ### GC B cell lifetime under force
                         "tau0",     ### initial lifetime
                         "xb",       ### GC B cell bond length
                         "xbStd",    ### GC B cell bond length std
                         "ep",       ### plasma cell affinity, mean
                         "epStd",    ### plasma cell affinity, std
                         "tau_p",    ### plasma cell bond lifetime
                         "tau_pStd", ### plams cell lifetime std
                         "eab",      ### Ab pool mean affinity
                         "eabStd",   ### Ab pool affinity std
                         "eabVar",   ### Ab pool affinity var
                         "xb_ab",    ### Ab pool mean xb
                         "xb_abStd", ### Ab pool xb std
                         "Ec",       ### thresholding feedback affinity
                         "tau_ab",   ### Ab pool quality
                         "tau_abStd",### Ab pool quality std
                         "n",        ### B cell pop size
                         "np",       ###  plasma cell pop size
                         "nab",      ### antibody number
                         "w",        ### B cell fitness, mean
                         "wStd",     ### B cell fitness, std
                         "eta",      ### B cell extraction chance, mean
                         "etaStd",   ### extraction chance, std
                         "cov",      ### covariance between fitness and E
                         "cov_lmd",  ### cov / mean fit
                         "cov2",     ### covariance between fitness and E^2
                         "de",       ### affinity gap between Ab and B cell]
                         "deltaE",   ### affinity diff
                         "sens",     ### sensitivity, d lmd / d Gb
                         "sens_var", ### sens * var
                         "finished", ### finished or not, single value,
                         "mat_rate", ### maturation rate
                         "f",        ### force
                         "Nb",       ### bottle neck pop size
                        ]
            
        self.dataset = {}
                          
        self.eta_record = {} 

        self.no_selection = False
        
        self.setup()
        pass
    
    def load_prm(self, prm):
        ### load parameters
        self.Td = int(prm["Td"])  ### feedback delay time
        self.Tp = int(prm["Tp"])  ### plasma cell lifetime
        self.xb1 = prm["xb1"]     ### APC-Ag bond length
        self.xb2 = prm["xb2"]     ### BCR-Ag bond length
        self.w0 = prm["w0"]       ### max proliferation rate
        
        self.N0 = prm["N0"]       ### initial B cell pop size
        self.cag = prm["cag"]     ### antigen density
        self.dE = prm["dE"]       ### mutation step size for barrier height
        self.dxb = prm["dxb"]     ### mutation step size for bond length
        
        self.pa = prm["pa"]       ### probability for apoptosis
        self.pd = prm["pd"]       ### probability for death
        self.pm = prm["pm"]       ### probability for mutation
        self.pm_xb = prm["pm_xb"]
        
        self.Ea = prm["Ea"]       ### tether affinity, no feedback only
        self.Eb = prm["Eb"]       ### naive B cell barrier height
        
        self.df = prm["df"]       ### force increament per GC cycle
        
        self.f = prm["f"]         ### force 
        
        self.eta0 =prm["eta0"]    ### half saturation in fitness when kernel=2
        self.Nc = prm["Nc"]       ### GC B cell pop capacity
        self.Npc = prm["Npc"]
        
        self.Nab = prm["Nab"]     ### feedback Ab pool size
        
        self.cutoff = prm["cutoff"]
        
        ## for feedback only
        self.feedback = prm["feedback"]     ### bool, if feedback or not
        
        self.randomAb = prm["randomAb"]     ### randomize tether affinity
        
        self.useSim = prm["useSim"]         ### use numerical simulation to obtain eta
        if self.useSim:
            self.sto = sb.System()
            
        self.update_rule = prm["update_rule"]
        
        self.output = prm["output"]         ### write out infomration and warnnings
        
        self.useWeek = prm["useWeek"]       ### False: use GC cycle as time, True: use weeks
        
        self.debug= prm["debug"]
        
        self.useBell = prm["useBell"]       ### use Bell formula
        
        self.goodDiff = prm["goodDiff"]     ### only chose B cell with fitness above the average for differentiation. 
        
        self.useBinomial = prm["useBinomial"]     ### use n_ag instead of eta to compute fitness
        
        self.no_selection = prm["no_selection"]   ### pure diffusion
        
        self.potential = prm["potential"]     ### type of interacing potential
        
        self.dump_eta = prm["dump_eta"]       ### save eta or not
        
        self.random_death = (prm["death"] == "random")   ### random death          
        
        self.eta_record_file = "eta_record/"+self.potential +"/" + prm["eta_file"]
        
        return
 
    
    
    
        
    def setup(self):
        self.load_prm(self.prm)
        
        self.N = 0        ### num of B cells
        self.NP = 0       ### num of plasma cell
        self.Ec = self.Eb ### threshold affinity to enter the feedback pool
        
        self.t = 0        ### current time
        self.tm = 0       ### total time when finished
        self.finished = False   ## True if time reached tm
        
        self._used_sim_already= False 
        self._warned = False
        
        self.agents = []     ### GC B cell list
        self.plasma = []     ### plasma cell list
        self.Ab_pool = []    ### feedback Ab pool
        self.my_Ab_pool = [] ### current Ab pool
        self.Ab_pool_taul_list = [] ### current Ab pool binding quality
        self.tether_aff = self.Eb ### tether affinity
        
        self.division_record = [] ### record number of divisions for all B cells
        for qty in self.qty_name:
            self.dataset[qty] = [] ### store all quantity of interest

        self.etaList = []    ### current eta list
        
        self.fit_mean, self.fit_max = 0, 0 ### mean fitness, max fitness
        
        for i in range(self.N0):
            ### initialize N0 B cells
            bcell = Bcell(self, self.Eb, self.xb2)
            self.addAliveBcell(bcell)
            
        self.checkN() ### double check if population size is syncrhonized
        
        self.N_ave = self.N  ### pop size averaing dark zone step and light zone step
        self.N_min = self.N  ### minimal pop size during evolution
        self.N_last = self.N ### pop size of previous cycle, to smooth the pop curve when it is close to the capacity
        
        self.lightZoneStep() ### initialize extraction
        self.collectData()   ### collect all data
        self.dataset["tau0"].append(theory.lifetime_linear_cubic(self.Eb, self.xb2))
        
        self.alive = True
        return
    
    
    def run(self, tm=10):
        if self.dump_eta:    
            self.load_eta_record()
        while self.t < tm and self.N>0:
            self.step()
            if self.output:
                sb.printProgress(self.t, tm)
            if self.N > 2E4:
                write("warning: population too large %d"%self.N)
        self.tm = self.t + 1
        
        if self.N == 0:
            self.alive = False
        if self.dump_eta and len(self.eta_record)>0:   
            self.dump_eta_record()
        if self.output:
            if self.N == 0:
                write("Info: GC dead")
            else:
                write("Info: finished!")
        if self.N == 0:
            self.dataset["finished"] = False
            return False
        else:
            self.dataset["finished"] = True
            return True
        
    
    def step(self):
        '''
        a complete GC cycle
        '''
        self.N_ave = self.N/2
        
        ######## light zone step: extract antigen, differentiate, die ######
        self.lightZoneStep()
        
        ######## remove dead B cells ######
        self.clean()
        self.N_ave += self.N/2
        
        ######## collecting data #####
        self.collectData()
        
        ######## dark zone step: divide and mutate ####
        self.darkZoneStep()
        
        ######## remove dead B cells #####
        self.clean()
        
        ####### plasma cell step: die #####
        self.plasmaStep()
        
        
        ####### update force if force is dynamic #####
        self.update_force()
        
        ###### update time #####
        self.t += 1
        return
    
    
    
    def lightZoneStep(self, test=False):
        '''
        light zone step
        1. update the feedback Ab pool at different time
        2. get the current feedback Ab pool: my_Ab_pool
        3. extract antigens
        '''
        if self.feedback and not test:
            ### update the feedback Ab pool
            self.update_Ab_pool()
        
        ## intial eta list
        if self.feedback and self.t>self.Td:
            if self.randomAb: ### random sample Ab as tether
                self.my_Ab_pool = self.Ab_pool[self.t-self.Td] ## feedback delay
                np.random.shuffle(self.my_Ab_pool)             ## random shuffle of the feedback Ab pool
            else: ### use the mean property as tether
                self.my_Ab_pool = [np.mean(self.Ab_pool[self.t-self.Td], axis=0).tolist()]
        elif self.feedback:  ### if no feedback Ab has been generated
            self.my_Ab_pool = [(self.Eb, self.xb2)]
        else: ## no feedback
            self.my_Ab_pool = [(self.Ea, self.xb1)]

        ## B cells extract antigens 
        Ab_num = len(self.my_Ab_pool)
        order = 0
        for bcell in self.agents:
            if bcell.alive:
                eta = bcell.extract(*self.my_Ab_pool[order % Ab_num])
                order += 1
        self.etaList = self.getCurrentEtaList()
        return
    
    def darkZoneStep(self): 
        
        '''
        dark zone step, including
        1. apoptosis: this happened in the light zone. 
                      Here we put it in the dark zone step right after antigen extraction.
                      Thus does not change the results at all. 
                      The only difference is when we collect the data
        2. differentiate: this happened in the light zone. Placing it here does not change the results
        3. proliferate and mutate
        '''
        ### apoptosis
        threshold = -1
        if not self.random_death:
            ### select the threshold affinity for survive according to the fraction of survive
            threshold = self.quick_select(self.getCurrentAffinityList()) 
 
        for bcell in self.agents:
            if bcell.alive:
                bcell.apoptosis(threshold) ### if threshold = -1, die randomly
                
        ### differentiate
        if self.goodDiff:
            mean_eta = np.mean(self.etaList)
        else:
            mean_eta = -1
        for bcell in self.agents:
            if bcell.alive:
                bcell.differentiation(mean_eta)
                
                
        ### proliferation and mutation
        daughter_cell_list = []
        self.division_record = []
        
        for bcell in self.agents:
            if bcell.alive:
                ret = bcell.comp_divide()
                if ret:
                    daughter_cell_list += ret
                self.division_record.append([bcell.eta, len(ret)])
        for bcell in daughter_cell_list:
            self.addAliveBcell(bcell)   
        return
    
    
    def plasmaStep(self):
        #### plasma cell dies after Tp
        for cell in self.plasma:
            if cell.alive:
                if self.t - cell.t0 > cell.T and cell.T>0:
                    cell.goDie()
        newPlasmaList = []
        for cell in self.plasma:
            if cell.alive:
                newPlasmaList.append(cell)
            else:
                self.NP -= 1
        self.plasma = newPlasmaList
        return
    
    
    def update_force(self):
        self.f += self.df
        return
    
    def collectData(self):
        '''
         collect data of interest during GC evolution
        '''
        # check if N = number of alive B cell or not
        self.checkN()
        
        
        ## record number of B cells ## averaging over the entire circle
        self.append("n", (self.N_ave + self.N_last)/2)
        
        self.N_min = min(self.N_min, (self.N_ave + self.N_last)/2)
        self.append("Nb", self.N_min)
        
        self.N_last = self.N_ave
        
        ## record affinity distribution
        elist = self.getCurrentAffinityList()
        eb = self._get_mean(elist)
        eb_var = self._get_var(elist)
        self.append("e", eb)
        self.append("eStd", np.sqrt(eb_var))
        self.append("eVar", eb_var)
        self.append("eMp", self._get_most_prob(elist))
        
        if len(self.dataset["e"]) > self.cutoff:
            self.append("mat_rate", np.mean(diff(self.dataset["e"])[-self.cutoff:]))
        else:
            if len(self.dataset["e"]) >1:
                self.append("mat_rate", np.mean(diff(self.dataset["e"])))
            else:
                self.append("mat_rate", 0)
        ### record bond length distribution
        xblist = self.getCurrentXbList()
        self.append("xb", self._get_mean(xblist))
        self.append("xbStd", self._get_std(xblist))
        
        ## record lifetime distribution
        taulist = self.getCurrentLifetimeList()
        tau_mean = self._get_mean(taulist)
        tau_var = self._get_var(taulist)
        self.append("tau", tau_mean)
        self.append("tauVar", tau_var)
        self.append("tauStd", np.sqrt(tau_var))
        
        ### record lifetime under force
        self.append("tauF",  theory.lifetime_linear_cubic(eb, self._get_mean(xblist), f0=self.f))
            
        ## record number of plasma cells
        self.append("np", self.NP)
        
        ## record plasma cell affinity
        eplist = self.getCurrentPlasmaAffinityList()
        self.append("ep", self._get_mean(eplist))
        self.append("epStd", self._get_std(eplist))
        
        ### record plasma cell lifetime distribution
        tau_p_list = self.getCurrentPlasmaLifetimeList()
        self.append("tau_p", self._get_mean(tau_p_list))
        self.append("tau_pStd", np.sqrt(self._get_var(tau_p_list)))
        
        
        ### record Ab pool affinity
        Ab_pool_Eb = [Ab[0] for Ab in self.my_Ab_pool]
        eab = self._get_mean(Ab_pool_Eb)
        eab_std = self._get_std(Ab_pool_Eb)
        self.append("eab", eab)
        self.append("eabStd", eab_std)
        self.append("eabVar", eab_std**2)
        self.append("nab", len(self.my_Ab_pool))
        self.append("Ec", min(Ab_pool_Eb))
        
        Ab_pool_xb = [Ab[1] for Ab in self.my_Ab_pool]
        self.append("xb_ab", self._get_mean(Ab_pool_xb))
        self.append("xb_abStd", self._get_std(Ab_pool_xb))
        
        Ab_pool_tau = [theory.lifetime_linear_cubic(Ab[0], Ab[1]) for Ab in self.my_Ab_pool]
        self.append("tau_ab", self._get_mean(Ab_pool_tau))
        self.append("tau_abStd", self._get_std(Ab_pool_tau))
        
        ### record extraction list
        etaList = self.getCurrentEtaList()
        eta = np.mean(etaList)
        self.append("eta", eta)
        self.append("etaStd", np.std(etaList))
        ## record extract chance, fitness and affinity
        
        fitList = self.getCurrentFitList()
        fit = self._get_mean(fitList)
        self.fit_mean = fit
        self.fit_max = self._get_max(fitList)
        self.append("w", fit)
        self.append("wStd", self._get_std(fitList))
        
        if self.N>0:
            self.append("de", eb-eab)
        else:
            self.append("de", np.nan)
        
        cov = my_cov(elist, fitList)
        self.append("cov", cov)
        self.append("cov_lmd", cov/fit)
        
        e2list = [ei**2 for ei in elist]
        self.append("cov2", my_cov(e2list, fitList))
        
        if self.tm == 0 or self.N<1:
            deltaE = 0
        else:
            deltaE = eb - self.dataset["e"][self.tm-1]
        self.append("deltaE", deltaE)
        
        self.append("f", self.f)
        return
    
    
    
    ################### fuctions related to B cell selection ####################
    def Ag_extract(self, Ea, Eb, xa, xb):
        '''
        this function is called in each B cell: bcell.extract()
        @parameter: Ea: tether affinity; Eb: BCR affinity; xa: tether bond length; xb: BCR bond length
        @return: extraction probability or nag
        '''
        if self.useBell:
            eta0, flag, barrier_height = theory.extRatio_linear_cubic(xa, xb, Ea, Eb, 0, output=False, Eb_min=self.prm["Eb_min"])
            tau0 = 1/eta0 - 1
            eta = 1.0/(1.0+tau0*np.exp(self.f*(xb-xa)/4.012))
        else:
            if self.potential=="cusp":    
                eta, flag = theory.extRatio(xa, xb, Ea, Eb, self.f, output=False)
            elif self.potential=="linear-cubic":
                eta, flag, barrier_height = theory.extRatio_linear_cubic(xa, xb, Ea, Eb, self.f, output=False, Eb_min=self.prm["Eb_min"])
            if not flag and not self._warned and self.output and False:
                print("*GC warning: low barrier! Ea={0:.1f}, Eb={1:.1f}, dE={2:.2f}, xa={5:.2f}, xb={4:.2f}, f={3:.3f}".format(Ea, Eb, barrier_height, self.f, xb, xa))
                self._warned=True

            if not flag and self.useSim:
                if self.output and not self._used_sim_already:
                    if eta is None and self.output and False:
                        print("*start sim @t={0:.0f}, F={1:.2f}, eta_theory=None".format(self.t, self.f), end=", ")
                        pass
                    elif self.output and False:
                        print("*start sim @t={0:.0f}, F={1:.2f}, eta_theory={2:.3f}".format(self.t, self.f, eta), end=", ")
                eta = self.read_eta_from_record(Ea, Eb, xa, xb, self.f)

        if self.useBinomial:
            #### get antigen number
            return np.random.binomial(self.cag, eta)
        else:
            #### use extraction chance directly
            return eta
        
        

        
        
    def Tcell_selection(self, eta):
        '''
        determine B cell fitness according to the extraction chance
        this function is called by a B cell
        @return fitness
        '''
        if self.no_selection:
            ## weak selection
            if eta<3:
                return self.w0*np.exp(0.1*eta)*(1-self.N/self.Nc)/np.exp(0.3)
            else:
                return self.w0*(1-self.N/self.Nc)
        
        if self.useBinomial:
            return self.w0*eta/(self.cag*self.eta0+eta) * (1-self.N/self.Nc)
        else:
            return self.w0*eta/(self.eta0+eta) * (1-self.N/self.Nc)
    
    def update_Ab_pool(self):
        '''
         update Ab pool
        '''
        if not self.feedback:
            print("warning! updating tether pool even without feedback")
            return
        if self.update_rule == "all":
            plasma_aff = self.getCurrentPlasmaAffinityList(include_xb=True)+[(self.Eb, self.xb2)]
            self.Ab_pool.append(plasma_aff)
        elif self.update_rule == "topK":
            self.Ab_pool.append(self.getCurrentTopPlasmaAffinityList(self.Nab))
        else:
            raise Exception("updating rule not found")
        return
    

    #################### useful functions ####################
    def run_sim(self, Ea, Eb,  dE=1, xa=1.5, xb=2.0, f=None):
        '''
        run Brownian simulations to get eta
        '''
        
        if not self.useSim:
            return None
        self.sto.prm["Eb1"] = Ea
        self.sto.prm["Eb2"] = Eb
        if f is None:
            self.sto.prm["f0"] = self.f
        else:
            self.sto.prm["f0"] = f
        self.sto.prm["xb1"] = xa
        self.sto.prm["xb2"] = xb
        
        self.sto.prm["potential"] = self.potential
        if dE>3:
            self.sto.prm["dt"] = 2
            self.sto.numSample=100
        elif dE<=0.1:
            self.sto.prm["dt"] = 0.2
            self.sto.numSample=400
        else:
            self.sto.prm["dt"] = 1
            self.sto.numSample=200
            
        self.sto.setup()
        print("@", end="")
        return self.sto.run()
    
    def interpolation(key):
        if key in self.eta_record:
            return self.eta_record[key]
        return 0
       
    def load_eta_record(self):
        try:
            with open(self.eta_record_file, "r") as myfile:
                while True:
                    line = myfile.readline()
                    if not line:
                        break
                    s = line.split('\t')
                    if self.feedback:
                        key = (int(s[0]), int(s[1]), int(s[2]), int(s[3]), int(s[4]))
                        val = float(s[5])
                    else:
                        key = (int(s[0]), int(s[1]), int(s[2]))
                        val = float(s[3])
                    self.eta_record[key] = val
        except:
            #self.eta_record = {}
            pass
        return
    
    def dump_eta_record(self):
        with open(self.eta_record_file, "w") as myfile:
            for key, val in self.eta_record.items():
                if self.feedback:
                    myfile.write("{0:d}\t{1:d}\t{2:d}\t{3:.4f}\n".format(key[0], key[1], key[2], key[3], key[4], val))
                else:
                    myfile.write("{0:d}\t{1:d}\t{2:.4f}\n".format(key[0], key[1], key[2], val))
        return
    
    def read_eta_from_record(self, Ea, Eb, xa, xb, f):
        if not self._used_sim_already and self.output and False:
            print("using simulation!") 
        key = ( my_round(Ea), my_round(Eb), my_round(xa), my_round(xb), my_round(f) )
        if key in self.eta_record:
            eta = self.eta_record[key]
        else:
            eta = self.run_sim(Ea, Eb, xa=xa, xb=xb, f=f)
            self.eta_record[key] = eta
            self._used_sim_already=True
        return eta
    
        
    def checkN(self):
        assert len(self.agents)==self.N, \
            "at t={0:d} population wrong: len(agents)={1:d} N={2:d}".format(self.t, len(self.agents), self.N)
        
    def clean(self):
        ## remove dead B cells form germinal center
        newAgents = []
        for bcell in self.agents:
            if bcell.alive:
                newAgents.append(bcell) 
        self.agents = newAgents
        self.N = len(self.agents)
        self.checkN()
        return
    
    def addAliveBcell(self, bcell):
        if bcell.alive:
            self.agents.append(bcell)
            self.N += 1
            
        return
    
    def removeDeadBcell(self, bcell):
        if not bcell.alive:
            self.agents.remove(bcell)
            self.N -= 1
        return
        
    def getCurrentPlasmaAffinityList(self, include_xb=False):
        Elist = []
        for plasma in self.plasma:
            if plasma.alive:
                if include_xb:
                    Elist.append( (plasma.E, plasma.xb) )
                else:
                    Elist.append(plasma.E)
        return Elist
    
    def getCurrentPlasmaLifetimeList(self):
        taulist = []
        for plasma in self.plasma:
            if plasma.alive:
                taulist.append(plasma.tau)
        return taulist
    
    def getCurrentFeedbackAbList(self):
        return [Ab[0] for Ab in self.my_Ab_pool]
    
    def getCurrentTopPlasmaAffinityList(self, k=1):
        tau0 = theory.lifetime_linear_cubic( self.Eb, self.xb2)
        Elist = [(tau0, self.Eb, self.xb2)]
        for plasma in self.plasma:
            if plasma.alive: # and plasma.tau > tau0:
                if len(Elist) < k:
                    heapq.heappush(Elist, (plasma.tau, plasma.E, plasma.xb))
                else:
                    heapq.heappushpop(Elist, (plasma.tau, plasma.E, plasma.xb)) ### first push and then pop. Always maintain the largest K values
        
        self.Ec = Elist[0][1]
        self.Ab_pool_taul_list = [ei[0] for ei in Elist]
        return [(ei[1], ei[2]) for ei in Elist]
    
    
        
    def getCurrentAffinityList(self):
        Elist = []
        for bcell in self.agents:
            if bcell.alive:
                Elist.append(bcell.E)
        return Elist
    
    def getCurrentLifetimeList(self):
        taulist = []
        for bcell in self.agents:
            if bcell.alive:
                taulist.append(bcell.tau)
        return taulist
    
    def getCurrentXbList(self):
        xblist = []
        for bcell in self.agents:
            if bcell.alive:
                xblist.append(bcell.xb)
        return xblist

    def getCurrentEtaList(self):
        etaList = []
        for bcell in self.agents:
            if bcell.alive and bcell.mature:
                ### skip the new born B cell
                etaList.append(bcell.eta)
                #raise Exception("some B cells have zero eta")
        return etaList
    
    def getCurrentFitList(self):
        fitList= []
        for bcell in self.agents:
            if bcell.alive and bcell.mature:
                fitList.append((bcell.w*self.pa*(1-self.pd)))
        return fitList
    
    def _get_sens(self, eta, Eb):
        deta = self.eta0*(1-eta)/(self.eta0+eta)
        Ef = self.f*self.xb2/(2*Eb*4.012)
        dtau = 1#(1-1.5/Eb - Ef**2-(Ef/Eb)/(1-Ef))
        dw =deta*dtau
        return dw
    
    def _get_sens2(self, w0, w, Eb0, Eb):
        return (w-w0)/(w*(Eb-Eb0))
    
    def _get_sens3(self, ea, eb):
        de = 0.01
        eta1, flag1, _ = theory.extRatio_linear_cubic(e10=ea, e20=eb+de, xb1=self.xb1, xb2=self.xb2, f0=self.f)
        eta2, flag2, _= theory.extRatio_linear_cubic(e10=ea, e20=eb-de, xb1=self.xb1, xb2=self.xb2, f0=self.f)
        
        if not (flag1 and flag2):
            return np.nan
        
        eta = (eta1+eta2)/2
        eta0 = self.prm["eta0"]
        deta_dGb = (eta1-eta2)/(2*de)
        return eta0/(eta*(eta0+eta))*deta_dGb
        
    def _get_mean(self, alist):
        if len(alist)>0:
            return np.mean(alist)
        else:
            return np.nan
    def _get_std(self, alist):
        if len(alist)>0:
            return np.std(alist)
        else:
            return np.nan
    def _get_max(self, alist):
        if len(alist)>0:
            return max(alist)
        else:
            return np.nan
        
    def _get_var(self, alist):
        if len(alist)>0:
            return np.var(alist)
        else:
            return np.nan
        
    def _get_most_prob(self, array):
        hist, bins = np.histogram(array, bins=30)
        index = np.argmax(hist)
        return (bins[index] +bins[index+1])/2
    
    def quick_select(self, alist):
        if not alist:
            return -1
        
        k = int(len(alist)*self.pa)
        if k>= len(alist):
            return min(alist)
        return alist[self._quick_select(alist, 0, len(alist)-1, k)]
    
    def _quick_select(self, alist, left, right, k):
        if left >= right:
            return k
        
        p = self._partition(alist, left, right)
        if p==k:
            return p
        elif p>k:
            return self._quick_select(alist, left, p-1, k)
        else:
            return self._quick_select(alist, p+1, right, k)
        
    
    def _partition(self, alist, left, right):
        pivot = alist[right]
        a = left
        for i in range(left, right):
            if alist[i]>=pivot:
                alist[i], alist[a] = alist[a], alist[i]
                a +=1
                
        alist[right], alist[a] = alist[a], alist[right]
        return a
        
    def _print(self, s):
        if type(s) == type(""):
            print(s)
        else:
            for si in s:
                print(si, end=", ")
            print("")
        return
        
    def append(self, qty, value):
        
        if qty in self.qty_name:
            
            self.dataset[qty].append(value)
        else:
            self.qty_name.append(qty)
            self.dataset[qty] = [value]
    
    def addPlasma(self, E, xb):
        if INFINITE_PLASMA or self.NP < self.Npc:
            self.plasma.append(PlasmaCell(self, E, xb, self.t))
            self.NP += 1
        else:
            ### if plasma pool is full, randomly delete one plasma cell and add the new one
            rand_idx = np.random.randint(self.NP)
            self.plasma[rand_idx] = PlasmaCell(self, E, xb, self.t)
        return
    
    def _death_process(self):
        ### apoptosis
        for bcell in self.agents:
            if bcell.alive:
                bcell.apoptosis()
        return
    
    def _birth_process(self):
        ### proliferation and mutation
        daughter_cell_list = []
        
        for bcell in self.agents:
            if bcell.alive:
                ret = bcell.comp_divide()
                if ret:
                    daughter_cell_list += ret
        
        for bcell in daughter_cell_list:
            self.addAliveBcell(bcell)
        return
    
    
    def plotQty(self, qty, cr="r", label="", ax=None, filling=False,  **keyargs):
        if not ax:
            fig, ax = plt.subplots(figsize=(4,3), dpi=100)
            ax.set_xlabel("GC cycle", fontsize=12)
            
            ax.set_ylabel(qty, fontsize=12)
            
        if filling:
            mean = np.asarray(self.dataset[qty])
            std = np.asarray(self.dataset[qty+"Std"])
            for i in range(len(mean)):
                if mean[i] is None:
                    continue
                start_id = i
                break
            if self.useWeek:
                tlist = convert_to_week(range(len(mean)))
            else:
                tlist = range(len(mean))
            ax.plot(tlist[start_id:], mean[start_id:],  color=cr, label=label,  **keyargs)
            ax.fill_between(tlist[start_id:], mean[start_id:]-std[start_id:], mean[start_id:]+std[start_id:], color=cr, alpha=.3)
        else:
            if self.useWeek:
                tlist = convert_to_week(range(len(self.dataset[qty])))
            else:
                tlist = range(self.tm)
            ax.plot(tlist, self.dataset[qty], color=cr, label=label)
        return ax