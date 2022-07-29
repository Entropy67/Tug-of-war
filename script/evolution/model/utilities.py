"""
utitity functions
"""
import numpy as np
from time import gmtime, strftime
import matplotlib.pyplot as plt
from termcolor import colored

week_per_cycle = 0.5/7

def convert_to_week(tm):
    try:
        return tm*week_per_cycle
    except:
        tmp = []
        for ti in tm:
            tmp.append(ti*week_per_cycle)
        return tmp

def convert_to_affinity(ei):
    return np.exp(np.asarray(ei))/np.exp(14)
#################################################################


def write(s, path=""):
    if path == "":
        print(s)
    else:
        t = "\n" + strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ":  " + s
        with open(path, 'a') as log:
            log.write(t)
    return

    
def diff(a, step=1):
    return [(a[i+step]-a[i])/step for i in range(len(a)-step)]

def my_cov(a, b):
    if len(a) != len(b):
        write("Error, cov array length doesn't match")
        raise ValueError()
    if len(a)>0:
        return np.inner(a, b)/len(a)-np.mean(a)*np.mean(b)
    else:
        return None

def get_mean(a):
    b = []
    lengthlist = [len(ai) for ai in a]
    if len(lengthlist)==0:
        return 0
    length = min(lengthlist)
    n = len(a)
    for j in range(length):
        s = 0
        for i in range(n):
            s += a[i][j]
        b.append(s/n)
    return np.asarray(b)

def get_std(a):
    b = []
    lengthlist = [len(ai) for ai in a]
    if len(lengthlist)==0:
        return 0
    length = min(lengthlist)
    n = len(a)
    for j in range(length):
        s = 0
        b.append(np.std([a[k][j] for k in range(n)]))
    return np.asarray(b)

def setAxWidth(ax, frame_width):
    ax.spines['top'].set_linewidth(frame_width)
    ax.spines['right'].set_linewidth(frame_width)
    ax.spines['bottom'].set_linewidth(frame_width)
    ax.spines['left'].set_linewidth(frame_width)
    return
            
     
def getAx(xlabel="time", ylabel="quantity", fontsize=15, figsize=(6, 5)):
    fig, ax = plt.subplots(figsize=figsize, dpi=100)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    setAxWidth(ax, 1.5)
    return ax

def print_dict(dct, ncol=1):
    print("Items held:")
    if ncol==1:
        for item, amount in dct.items():  # dct.iteritems() in Python 2
            print("{}:\t {}".format(item, amount))
    else:
        idx = 0
        end = ["\t\t"]*(ncol-1) + ["\n"]
        max_length = max([len(k) for k in dct.keys()])
        fmt = "%"+str(max_length)+"s"
        for item, val in dct.items():
            print(fmt % item, end=": \t")
            print(colored("{}".format(val), "magenta"), end=end[idx % ncol])
            idx += 1

    return

        
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


def dump(dataset, filename, unique=True, mod='w'):
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


def my_round(x, tol=0.1):
    return int(round(x/tol))
