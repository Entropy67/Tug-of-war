"""
functions to scan the parameter space
"""


#### import modules
import matplotlib.pyplot as plt
import numpy as np
import json
import os.path
from mpmath import mp


class Data:
    
    """
    Data class
    """
    
    def __init__(self):
        self.qty_name = [] #### name of quantities to be recorded
        self.dataset = {} ### store all interested data
        pass
    
    def setup(self):
        """
        initialize the dataset
        """
        for qty in self.qty_name:
            self.dataset[qty] = [] 
    
    def copy(self, other_obj):
        """
        copy other data object
        """
        self.qty_name = other_obj.qty_name.copy()
        for qty in other_obj.qty_name:
            try:
                self.dataset[qty] = other_obj.dataset[qty].copy()
            except:
                self.dataset[qty] = []
        self._unzipDataset()
        return
    
    def append(self, qtys, values):
        if isinstance(qtys, list):
            for qty, v in zip(qtys, values):
                self._append(qty, v)
        else:
            self._append(qtys, values)
        return
    
    def _append(self, qty, value):
        """
        append new (qty, value) to the dataset
        arguments:
            qty: string, name of quantity
            value:, value of quantity
        """
        if qty in self.qty_name:
            self.dataset[qty].append(value)
        else:
            self.qty_name.append(qty)
            self.dataset[qty] = [value]
        return
    
    def appends(self, qty, value):
        """
        append a list of qty and value to the dataset
        """
        #assert(isinstance(qty, list), "input should be a list of qtys")
        for q, v in zip(qty, value):
            self.append(q, v)
        return
    
    def get(self, qty):
        """
        get the value of qty
        """
        if qty in self.qty_name:
            return np.asarray(self.dataset[qty])
        else:
            return []
    
    def _unzipDataset(self):
        ### convert from dataset to other lists
        pass
    
    def _zip2Dataset(self):
        ### convert from lists to dataset
        pass
    
    def dump(self, filename):
        """
        save data to file
        """
        self._zip2Dataset()

        ### save to json
        v = 0
        while os.path.exists(filename+'.json'):
            print("file exists! ")
            filename += gen_pi(v)
            v += 1
        
        with open(filename+'.json', 'w') as fp:
            json.dump(self.dataset, fp)
        return
    
    def load(self, filename):
        """
        load data from filename
        """
        ### load the json
        with open(filename, 'r') as fp:
            self.dataset = json.load(fp)
        
        self._unzipDataset()
        
        return
    
    def plotQty(self, qty, ax=None, **keyargs):
        """
        simple visualization of the data
        """
        if not ax:
            fig, ax = plt.subplots(figsize=(5, 4), dpi=100)
        plt.plot(prm_list, qty, **keyargs)
        return ax

    

class Scan_prm(Data):
    
    """
    derived class to scan parameters
    """
    
    def __init__(self, sto):
        super().__init__()
        self.sto = sto ### system class to simulate tug-of-war process
        self.etalist = [] ## eta list
        self.trlist = [] ## rupture time
        self.frlist = [] ## rupture force per bond
        self.prm_list = []  ### record parameter list
        self.prm_name = ""
        
        self.qty_name = ["prm_name", 
                         "prms", 
                         "eta", 
                         "fr", 
                         "tr", 
                         "etastd", 
                         "frstd", 
                         "trstd"
                        ]
        
        
        self.etastdlist, self.trstdlist, self.frstdlist = [], [], []
        self.labelSize= 14
        pass

    def run(self, prm_name, prm_list):
        """
        scan prm_name in prm_list
        arguments:
            prm_name: string, name of parameter to be scanned
            prm_list: list, a list of values of parameter to be scanned
        """
        print(prm_name, "\teta\ttend\tfend")
        self.prm_name = prm_name

        for prm in prm_list:
            self.sto.prm[prm_name] = prm
            self._run() #### run the simulation
            print("{0:.3f}\t{1:.3f}\t{2:.3f}\t{3:.3f}".format(
                prm, self.sto.eta,
                self.sto.tend, self.sto.fend
            ))
            self.prm_list.append(prm)
        return
            
    def _run(self):
        #self.sto.output=False
        self.sto.setup()
        self.sto.run() ### run the simulation
        
        self.etalist.append(self.sto.eta)
        self.trlist.append(self.sto.tend)
        self.frlist.append(self.sto.fend)
        self.frstdlist.append(self.sto.fend_std)
        self.trstdlist.append(self.sto.tend_std)

        assert len(self.frstdlist)>0, "frstd is not recorded"
        if self.sto.manyRun:
            self.etastdlist.append(self.sto.eta_std)
        return
    
    def clean(self):
        self.etalist = []
        self.frlist = []
        self.trlist = []
        self.prm_list = []
        self.etastdlist = []
        self.frstdlist = []
        self.trstdlist = []
        return 
    
    def copy(self, other_obj):
        self.etalist = other_obj.etalist
        self.frlist = other_obj.frlist
        self.trlist = other_obj.trlist
        self.prm_list = other_obj.prm_list
        self.prm_name = other_obj.prm_name
        
        self.etastdlist=other_obj.etastdlist
        self.trstdlist = other_obj.trstdlist
        self.frstdlist = other_obj.frstdlist
        return
    
    def _zip2dataset(self):
        ## zip to dataset
        self.dataset = {}
        dataset = [self.prm_name, self.prm_list, self.etalist,
                  self.frlist, self.trlist,
                  self.etastdlist, self.frstdlist, self.trstdlist]
        for i, qty in enumerate(self.qty_name):
            self.dataset[qty] = dataset[i]
            
        return
    
    def _unzipDataset(self):
        ### unzip to lists
        self.etalist = self.dataset["eta"]
        self.frlist = self.dataset["fr"]
        self.trlist = self.dataset["tr"]
        self.prm_list = self.dataset["prms"]
        self.prm_name = self.dataset["prm_name"]
        
        self.etastdlist=self.dataset["etastd"]
        self.trstdlist = self.dataset["trstd"]
        self.frstdlist = self.dataset["frstd"]
        return
    
    
    
    
    def plot(self, axes=None):
        if not axes:
            fig, axes = plt.subplots(figsize=(9,4), ncols=2)
        
        ax1= axes[0]
        if self.sto.manyRun:
            ax1.errorbar(self.prm_list, self.etalist, yerr=self.etastdlist, capsize=3)
        else:
            ax1.plot(self.prm_list, self.etalist, '-o')
        ax1.set_xlabel(self.prm_name, fontsize=self.labelSize)
        ax1.set_ylabel("extraction chance", fontsize=self.labelSize)
        
        ax2 = axes[1]
        if self.sto.manyRun:
            ax2.errorbar(self.prm_list, self.frlist, yerr=self.frstdlist, capsize=3)
        else:
            ax2.plot(self.prm_list, self.frlist, '-o')
        ax2.set_xlabel(self.prm_name, fontsize=self.labelSize)
        ax2.set_ylabel("rupture force", fontsize=self.labelSize)
        return axes
    
    def plotFr(self, ax=None, errorbar=True, fmt='-o', **keyargs):
        if not ax:
            fig, ax = plt.subplots(figsize=(4,3.3), dpi=100)
        if self.sto.manyRun and errorbar:
            ax.errorbar(self.frlist, self.etalist, yerr=self.etastdlist,fmt=fmt, capsize=5, alpha=0.5, fillstyle='none', **keyargs)
            ax.plot(self.frlist, self.etalist, fmt, fillstyle='none', **keyargs)
        else:
            ax.plot(self.frlist, self.etalist, fmt, fillstyle='none', **keyargs)
        ax.set_xlabel("rupture force, pN", fontsize=self.labelSize)
        ax.set_ylabel("extraction chance", fontsize=self.labelSize)
        
        return ax

    
    
    
    
class Dis_curve(Scan_prm):
    
    def __init__(self, sto):
        super().__init__(sto)
        
    def run(self, prm_list):
        return super().run("Eb2", prm_list)
        

    def plotDisCurve(self, ax=None, errorbar=True, **keywards):
        if not ax:
            fig, ax = plt.subplots(figsize=(4,3.3), dpi=100)
        
            plt.xlabel(r"BCR-Ag potential height, $\Delta G_2^{\dagger}$", fontsize=self.labelSize)
            plt.ylabel("extracted chance, $\eta$", fontsize=self.labelSize)
        
        if self.sto.manyRun and errorbar:
            ax.errorbar(self.prm_list, self.etalist, yerr=self.etastdlist, capsize=3, **keywards)
        else:
            ax.plot(self.prm_list, self.etalist, '-', **keywards)
        return ax



def findEb(sto, emin, emax, target_eta=0.5, tol=0.002, numSample=1000, output=True):
    error = 1
    emid = 0
    
    info=""
    
    def run(e2):
        sto.prm["Eb2"] = e2
        sto.setup()
        sto.run()
        return sto.eta
    
    sto.numSample = numSample
    
    eta_min = run(emin)
    if eta_min > target_eta:
        info_tmp = "emin too large, emin={0:.3f}, eta={1:.3f}\n".format(emin, eta_min)
        info += info_tmp
        if output:
            print(info_tmp)
        return None, info
    eta_max = run(emax)
    if eta_max < target_eta:
        info_tmp= "emax too small, emax={0:.3f}, eta={1:.3f}\n".format(emax, eta_max)
        info += info_tmp
        if output:
            print(info_tmp)
        return None, info
    info_tmp = "emid\teta\t error\n"
    info += info_tmp
    if output:
        print(info_tmp)
    while abs(error) > tol and emax - emin > tol:
        emid = (target_eta - eta_min)*(emax-emin)/(eta_max-eta_min) + emin
        run(emid)
        error = sto.eta - target_eta
        if error > 0:
            eta_max = sto.eta
            emax = emid
        else:
            eta_min = sto.eta
            emin = emid
        info_tmp = "{0:.3f}\t{1:.3f}\t{2:.3f}\n".format(emid, sto.eta, error)
        info += info_tmp
        if output:
            print(info_tmp)
        
    return emid, info


class FindDisRange:
    
    def __init__(self, sto):
        self.sto = sto
        self.log = ""
        self.tol = 0.005
        self.numSample = 2000
        self.output=True
        self.log_file = "log.txt"
    
    def run(self, e1min, e1max, e2min, e2max, e3min, e3max, output=False):
        self.log += "find Eb for eta=0.1\n"
        e_low, info = findEb(self.sto, e1min, e1max, target_eta=0.1,
                     tol=self.tol, numSample=self.numSample,output=output
                      )
        if self.output:
            print(info)
        
        self.log += info
        self.write_log(self.log)
        
        self.log += "find Eb for eta=0.5\n"
        e_mid, info = findEb(self.sto, e2min, e2max, target_eta=0.5,
                     tol=self.tol,numSample=self.numSample, output=output
                      )
        if self.output:
            print(info)
        self.log += info
        self.write_log(self.log)
        
        self.log += "find Eb for eta=0.9\n"
        e_high, info = findEb(self.sto, e3min, e3max, target_eta=0.9,
                     tol=self.tol,numSample=self.numSample, output=output
                      )
        if self.output:
            print(info)
        self.log += info
        self.write_log(self.log)
        
        return e_low, e_mid, e_high
    
    def write_log(self, info):
        filename= self.log_file
        if self.output:
            print(info)
        with open(filename, "a") as myfile:
            myfile.write(info)
        self.log = ""
        return
    
class Scan_disRange(Data):
    
    def __init__(self, sto):
        super().__init__()
        self.sto = sto
        self.qty_name=[
            "prm_name",
            "prms",
            "elow",
            "emid",
            "ehigh",
            "de"
        ]
        self.prm_list = []
        self.finder = FindDisRange(sto)
        self.finder.output=False
        self.elow_list = []
        self.emid_list = []
        self.ehigh_list = []
        self.de_list = []
        return
        
    def run(self, prm_name, prm_list, logfile="dogR.txt"):
        print(prm_name, "\tEb1\tEmid\tEb2\tdEb")
        
        self.finder.log_file = "log/" + logfile
        self.prm_name= prm_name
        
        for prm in prm_list:
            self.sto.prm[prm_name] = prm
            self.sto.setup()
            self.finder.write_log("run for "+ prm_name + "="+str(prm))
            
            e_low, e_mid, e_high = self.finder.run(0, 50, 15, 60, 10, 50)
            if e_low is None or e_high is None:
                print("elow or e_high is none", end=', ')
                print(self.finder.log)
                return
            self.elow_list.append(e_low)
            self.emid_list.append(e_mid)
            self.ehigh_list.append(e_high)
            self.de_list.append(e_high-e_low)
            self.prm_list.append(prm)
            print("{0:.3f}\t{1:.3f}\t{2:.3f}\t{3:.3f}\t{4:.4f}".format(prm, e_low, e_mid, e_high, e_high-e_low))
            self.finder.write_log("log/"+logfile)
        self._zip2Dataset()
        return

        
    def _zip2Dataset(self):
        self.dataset= {}
        dataset = [self.prm_name, self.prm_list, self.elow_list, self.emid_list, self.ehigh_list, self.de_list]

        for i, qty in enumerate(self.qty_name):
            self.dataset[qty] = dataset[i]
        return
    
    def _unzipDataset(self):
        self.prm_name = self.dataset["prm_name"]
        self.prm_list = self.dataset["prms"]
        self.elow_list = self.dataset["elow"]
        self.emid_list = self.dataset["emid"]
        self.ehigh_list = self.dataset["ehigh"]
        self.de_list = self.dataset["de"]
        return
    
    def plot(self):
        fig, ax = plt.subplot(figsize=(5,4), dpi=100)
        plt.plot(self.prm_list, self.de_list, '-o')
        plt.show()
        return
    
    def plot_disRange(self,ax=None, twinx=False, color='steelblue', filename=""):

        if not ax:
            fig,ax = plt.subplots(figsize=(4, 3.3), dpi=100)
            ax.set_xlabel("Pulling force, f, pN", fontsize=14)
            ax.set_ylabel("BCR affinity range", fontsize=14)
        #plt.plot(dog.prm_list, dog.de_list, '-o')
        ax.plot(self.prm_list, self.elow_list, '-o', color=color, fillstyle='none', ms=4)
        ax.plot(self.prm_list, self.ehigh_list,'-o', color=color, fillstyle='none',ms=4)
        ax.fill_between(self.prm_list, self.elow_list, self.ehigh_list, where=self.ehigh_list>self.elow_list, alpha=0.4,color=color)


        #fx = np.linspace(0, 10, 100)
        #ax.plot(fx, prm["Eb1"]-np.log(0.2*9)-fx*0.5/4.012, '--', color='k', alpha=0.6)
        #ax.plot(fx, prm["Eb1"]-np.log(0.1/9)-fx*0.5/4.012, '--', color='k', alpha=0.6)


        if twinx:
            tx=ax.twinx()
            tx.plot(self.prm_list, self.de_list, '-or', ms=4)
        if filename:
            plt.tight_layout()
            plt.savefig(filename, format='pdf')
        return ax
    
def gen_pi(n):
    """
    generate n-th digit of pi
    """
    if n==0:
        return '_3'
    elif n==1:
        return '.1'
    else:
        return str(mp.pi)[n+1]