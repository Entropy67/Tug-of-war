# Tug-of-war

---

Tug-of-war is a python project to simulate the dynamics of antigen extraction process and its influence on GC evolution. The paper is now available on Arxiv:  Hongda Jiang and Shenshen Wang, [**Immune cells use active tugging forces to distinguish affinity and accelerate evolution**](https://arxiv.org/abs/2209.13991).

![new tug-of-war image](/tug-of-war-demo.png)

Basic examples are given in the jupyter notebook.


## Modules

The program contains two main parts: **Brownian motion** simulation and GC **evolution** simulation. 

- Brownian motion

- > We simulate the Langevin equation and calculate the probability to extract an antigen. 

- GC evolution

- >  we perform agent-based simulation to simulate the interactions between GC B cells, antigen presenting cells and plasma cells. Below we introduce main classes. 



### brownian

#### What does the package contain?

This module contains the following files:

- *bonds.py* 

  > implement the **Bond** class. 

- *force.py* 

- > implement the **Force_prm** class. The main function is **get_f(time)**, which returns the force magnitude at any time according to a given force scheme. 

- *system.py* : 

- >  contains the **System** class, which represents the tug-of-war system. It has two main attributes: x<sub>1</sub> and x<sub>2</sub>. The method **run1()** simulates a single antigen extraction process. Method **run()** simulate many realizations of antigen extraction.  

- *scan_parameter.py* 

- > this library makes it easier to scan different parameters. It implememts two main classes:  **Data** and **Scan_prm**. The Data class utlizes a dict to store all the quantities of interest, which can be saved to file or loaded from file conveniencely. The Scan_prm class runs simulations for a list of parameters and store the results of each parameter automatically. 

- *plot.py* 

- > This library helps to visulize the Brownian motion trajectories. 

- *theory.py* 

- > Analytical results.



#### How to use the package?

To simulate a single antigen extraction event, run the following code:

```python
import script.brownian.utilities as utl
import script.brownian.system as system
## prepare the parameter dict
prm = utilities.getDefaultPrm()
## create a System object
sys = system.System(prm)
sys.manyRun = False
### run a single event
flag, p, tend, fend = sys.run1()
print("simulation finished")
print(">> bond breaking? : ", flag)
print(">> which bond was broken? (1: BCR-Ag; 2: APC-Ag): ", p)
print(">> rupture time: ", tend, "us")
print(">> rupture force: ", fend, "pN")
```

To simulate many runs of antigen extraction and get statistics:

```python
sys.output=True
sys.manyRun = False ## change it to True if you want to know uncertainty of eta
sys.numSample = 1000
sys.run() ### will run Brownian simulation 1000 times
print("simulation finished")
print(">> antigen extraction chance: ", sys.eta)
print(">> average complex lifetime: ", sys.tend, "us")

```

To scan parameters:

```python
import script.brownian.scan_parameter as scan_parameter
#### create a Scan_prm object
agent = scan_parameter.Scan_prm(sto = sys)
agent.run("f0", [0, 5, 10, 15, 20, 25, 30]) ### will run simulations for f0=0, 5, 10, ..., 30

#### plot eta vs force:
fig, ax = plt.subplots(figsize=(3, 2), dpi=100)
plt.plot(agent.dataset["prms"], agent.dataset["eta"])
plt.show()
```



### evolution

#### What does this package contain?

- *controller.py*

- > Implement the **Controller** class, which provides basic functions to control a simulator. The main atrribute is **dataset**, which contains all the simulation results. It also makes it easier to save the results.

- *data.py*

- > **MyData** class provides methods to add, change and save the data. 

- *scan.py*

- > Similar to scan_parameter.py. The **Scanner** class inherits properties of **Controller**. It manages GC simulations under different parameter setting. 

- *utilities.py*

- > Utility functions

- *model*

- > GC simulation model

  - *bonds.py*

  - *singleBond.py*

  - *theory.py*

    > The above three packages are the same as what in the brownian package, providing Brownian simulations to get the extraction chance. 

  - *Bcell.py*

  - > Implement the **Bcell** class. A B cell can divide, mutate, differentiate, die and extract antigens. 

  - *plasmaCell.py*

  - > Implement the **PlasmaCell** class. A plasma cell will die after **T<sub>p</sub>** cycles.

  - *GC.py*

  - > Contains the **GC** class. A GC has a collection of B cells and plasma cells. In each GC cycle, it performs **lightZoneStep** and **darkZoneStep**, followed by data collection.

  - *manyRun.py*

  - > Simulate many realizations of GC evolution to get statistical properties. 

#### How to use the package?

To run a single GC simulation

```python
import script.evolution.model.GC as gc
import script.evolution.model.prm as prm

prm0 = prm.prm_list.copy()
## with Ab feedback
prm0["feedback"] = True
## without Ab feedback
prm0["feedback"] = False

gc0 = gc.GC(prm=prm0)
gc0.run(tm=100)
ax = gc0.plotQty("e", filling = True)
```

To run many GC simulations

```python
import evolution.model.manyRun as mr
mr0 = mr.manyRun(prm=prm0)
mr0.num_run = 5 ### number of GCs
mr0.run(output=True)
```

To scan parameters

```python
import evolution.scan_parameter as scan
prm1 = prm.prm_list.copy()

#### prepare the parameter scan info
prm1["prm_name"] = "f" ## the parameter to be scanned
prm1["prm_list_prm"] = [0, 41, 2] ### np.arange(0, 41, 2)
prm1["unit"] = 1
prm1["prm_scale"] = "linear" # "log10" etc.


### scan parameters
scanner = scan.Scaner(agent = mr0)
scanner.load_config(mr0.prm)
scanner.run()


### make plot
#### Make figure
def make_plot(scanner, qty="e", twinx=False):
    fig, ax = plt.subplots(figsize=(3.5,2.8), dpi=150)
    plt.errorbar(scanner.dataset.get("prms"), scanner.dataset.get(qty), yerr=scanner.dataset.get(qty+"_std"), capsize=3, fmt='ok', fillstyle='none', ms=5)
    #plt.ylim(14, 22)
    ax.set(xlabel="Force, F (pN)", ylabel="Output affinity (kT)")
    if twinx:
        tx= ax.twinx()
        tx.plot(scanner.dataset.get("prms"), scanner.dataset.get("surv_prob"), '--sr', fillstyle='none')
        tx.set_ylabel("Fraction of surviving GCs", color='r')
    return ax

ax = make_plot(scanner, qty="e", twinx=True)
ax.set_ylim(12, 24)
ax.patch.set_facecolor('none')
plt.tight_layout()
plt.show()
```
