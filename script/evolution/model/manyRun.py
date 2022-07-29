
"""
simulate many GC runs to get statistics
"""

from copy import deepcopy
from .GC import *
import matplotlib.pyplot as plt

################################
## get statistical property
class manyRun:
    def __init__(self, prm=prm_list):
        ## control penal
        self.prm = deepcopy(prm)
        
        self.output = True
        self.fix_success = False
        self.save_traj = True
        
        self.traj_qty_name = ["e",
                              "eStd",
                             
                              "tau", 
                              "tauStd", 
                             
                              "xb", 
                              "xbStd", 
                             
                              "eab", 
                              "eabStd",
                             
                              "xb_ab", 
                              "xb_abStd",
                         
                              "tau_ab", 
                              "tau_abStd",
                         
                              "ep", 
                              "epStd", 
                              
                              
                              "Ec",
                         
                              "tau_p", 
                              "tau_pStd",
                         
                              "n", 
                              "np", 
                              #"nab",
                         
                              "w", 
                              "wStd",
                         
                              "eta", 
                              "etaStd"
                             ]
        
        self.saved_qty = ["e", ### affinity at tmax
                          "eStd", ### affinity std
                          "xb", ### bond length at tmax
                          "xbStd", ### bond length std
                          "tau",  ### intrinsic lifetime
                          "tauStd",  ### lifetime std
                          "eab",  ### ab pool mean affinity
                          "eabStd", ### ab pool affinity std
                          "tau_ab", ### ab pool intrinsic lifetime
                          "n",  ### b cell pop size
                          "np", ### plasma cell pop size
                          "w",  ### fitness average
                          "wStd", ### fitness std
                          "eta",  ### average eta
                          "etaStd",
                          "tau0", ### initial lifetime] 
                          "mat_rate", ### maturation rate
                          "surv_prob", ### survival probability
                          "tend", ### ending time 
                          "Ec", 
                          "de",
                          "Nb"
                         ]
        self.trajdata = {}
        
        self.dataset = {}
        self.eta_record = {}
        self.mean_dataset = {}
        self.setup()
        self.sample_rate = 5
        self.num_run = 20
        return
    
    def set_prm(self, qty, value):
        self.prm[qty] = value
        self.setup()
        return 
    
    def setup(self):
        self.tm = self.prm["tm"]
        self.f=self.prm["f"]
        self.eta_record = {}
        self.gc = GC(f0=self.f, prm=self.prm)
        self.gc.output=False
        self.surv_prob = None
        if self.save_traj:
            for qty in self.traj_qty_name:
                self.trajdata[qty] = []
            
        for qty in self.saved_qty:
            self.dataset[qty] = []
        pass
    

    
    
    def dump_all(self, filename=""):
        dump(self.dataset, filename + "_dataset")
        dump(self.trajdata, filename + "_trajdata")
        return
        
    
    def load_eta_record(self):
        self.gc.load_eta_record()
        self.eta_record = self.gc.eta_record
        return
    
    def dump_eta_record(self):
        self.gc.eta_record = self.eta_record
        self.gc.dump_eta_record()
        return
    
    def run(self, output=False):
        
        self.setup()
        
        fix_success=self.fix_success
        success_count = 0
        dead_count = 0
        tot_count = 0

        mat_rate, var_tmp, eta_tmp = 0, 0, 0
        if self.output and output:
            print("manyRun: check parameter: f=%f" % self.prm["f"], ", xa= %f" % self.prm["xb1"])
        while success_count < self.num_run:
            if output and self.output:
                print(">>>> starting GC {0:d}, last run: mr= {1:.3f}, std={2:.3f}, eta={3:.3f}, tm={4:.1f}, success={5:d}".format(tot_count, mat_rate, var_tmp, eta_tmp, self.gc.t, success_count))
            
            self.gc.setup()
            self.gc.eta_record = self.eta_record
            self.gc.output = (output and self.output)
            success = self.gc.run(tm=self.tm)
            self.eta_record = self.gc.eta_record
            tot_count += 1
            
            if self.save_traj:
                self.append_traj()
            self.append_qty()
            
            self.dataset["tend"].append(self.gc.t)
            if not success:
                dead_count += 1
                if tot_count>2000:
                    if success_count == 0:
                        print("WARNING:---- All GC died!")
                        self.surv_prob = 0
                        return False
                    else:
                        break
            else:  
                success_count += 1
                var_tmp = np.mean(self.gc.dataset["eStd"][-20:])
                mat_rate = np.mean(diff(self.gc.dataset["e"])[-20:])
                eta_tmp = np.mean(self.gc.dataset["eta"][-20:])
            
            if (not fix_success) and tot_count>self.num_run:
                break
        self.surv_prob = success_count/(success_count+dead_count)
        self.dataset["surv_prob"] = [self.surv_prob]
        if output and self.output:
            print(colored("finished! success={0:d}, dead={1:d}, percen={2:.3f}\t".format(success_count, dead_count, self.surv_prob), "grey"))
        if success_count == 0:
            print("WARNING: ----- All GC died!")
            return False
        else:
            if self.output:
                print("many run: summary: surv_prob = {2:.4f}\t, Eb = {3:.3f}{5:.3f}\t, tau={4:.3f}".format(success_count, dead_count, self.surv_prob, np.mean(self.dataset["e"]), np.std(self.dataset["e"]), np.log10(np.mean(self.dataset["tau"])/self.gc.dataset["tau"][0])))
        
        #self.get_mean()
        return True
    
    def append_qty(self):
        for qty in self.saved_qty:
            if qty in self.gc.dataset.keys() and self.gc.alive:
                self.dataset[qty].append(self.gc.dataset[qty][-1])
        
    def append_traj(self):
        for qty in self.traj_qty_name:
            traj_tmp = self.gc.dataset[qty]
            self.trajdata[qty].append(traj_tmp[::self.sample_rate])
            ### include the last one
            if len(traj_tmp) % self.sample_rate != 1:
                self.trajdata[qty][-1].append(traj_tmp[-1])
        return
    
    def append(self, qty, value):
        self.dataset[qty].append(value)
        pass
    
    def get_data(self):
        return self.dataset
    
    def get_mean_traj(self, qty):
        return np.mean(self.trajdata[qty], axis=0)
    
    def plotQty(self, qty, cr="r", label="", ax=None, filling=False, multiDim=True, **keyargs):
        if not ax:
            fig, ax = plt.subplots(figsize=(4,3), dpi=100)
            ax.set_xlabel("time", fontsize=12)
            ax.set_ylabel(qty, fontsize=12)
        
        if multiDim:
            mean = np.mean(self.trajdata[qty], axis=0)
            std = np.std(self.trajdata[qty], axis=0)
            for i in range(len(mean)):
                if mean[i] is None:
                    continue
                start_id = i
                break
        else:
            mean = self.trajdata[qty]
            
        if self.prm["useWeek"]:
            tlist = convert_to_week(range(len(mean)))
        else:
            tlist = range(len(mean))
        ax.plot(tlist[start_id:], mean[start_id:],  color=cr, label=label,  **keyargs)
        
        if filling:
            ax.fill_between(tlist[start_id:], mean[start_id:]-std[start_id:], mean[start_id:]+std[start_id:], color=cr, alpha=.3)
        return ax
    
    def plot_combined(self, axes=None, cr='r', show_ab=True, filling=True):
        qty_name = ["B cell population",  "Affinity", "eta", "affinity gap", "variance", "sensitivity"]
        if axes is None:
            fig, axes = plt.subplots(figsize=(6, 8), dpi=100, ncols=2, nrows=3)
            plt.subplots_adjust(wspace=0.3)
            i=0
            for ax in axes:
                for axi in ax:
                    axi.set_ylabel(qty_name[i])
                    i += 1

        ax1= self.plotQty("n", ax=axes[0, 0], cr=cr, filling=filling)
        ax1.set_ylim(0, 1200)
        ax2= self.plotQty("e", ax=axes[0, 1], cr=cr, filling=filling)
        if show_ab:
            self.plotQty("eab", ax=axes[0, 1], filling=filling, cr=cr, linestyle='--')
        ax3=self.plotQty("eta", ax=axes[1, 0], cr=cr, filling=filling)
        ax3.set_ylim(0, 1)
        ax4=self.plotQty("de", ax=axes[1, 1], cr=cr, filling=filling)
        ax4.set_ylim(0, 10)
        
        ax5=self.plotQty("e_var", ax=axes[2, 0], cr=cr, filling=filling)
        ax5.set_ylim(0, 4)
        
        ax5=self.plotQty("sens", ax=axes[2, 1], cr=cr, filling=filling)
        return axes
    
    
        
def gen_pi(n):
    if n==0:
        return '_3'
    elif n==1:
        return '.1'
    else:
        return str(mp.pi)[n+1]
    
def save_prm(prm, filename):
    v = 0
    while os.path.exists(filename):
        print("file exists! ")
        filename += gen_pi(v)
        v += 1
        
    with open(filename+".txt", "w") as myfile:
        for item, amount in prm.items():  # dct.iteritems() in Python 2
            myfile.write("{}: {}".format(item, amount))
    return


def dump(self, dataset, filename, unique=True, mod='w'):
    """
    save data to file
    """

    if unique:
        ### create uniq filename
        v = 0
        while os.path.exists(filename+'.json'):
            print("file exists! ")
            filename += gen_pi(v)
            v += 1

    with open(filename+'.json', mod) as fp:
        json.dump(dataset, fp)
    return filename
    