"""
a class to visualize the trajectory in a 2D potential
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.cm as cm
from . import bonds as bd

PI = 3.1415926
kT = 300*1.38E-23
labelSize = 16
tickSize = 14

class Plot:
    
    def __init__(self, sys):
        
        self.model = sys
        
        self.coords = "x1x2" ## x1x2 or xy
        pass
    
    
    def convertToXY(self):
        t_traj, x1_traj, x2_traj = self.model._getTraj()
        self.x_traj = []
        self.y_traj = []
        for i in range(len(t_traj)):
            self.x_traj.append(x1_traj[i])
            self.y_traj.append(x1_traj[i]+x2_traj[i])
        
        return self.x_traj, self.y_traj
        
        
        
    
    def plotPotTraj(self, cr='r', alpha=0.4):
        tlist, potlist = self.model._getPotTraj()
        ax = self._plotPotTraj(tlist, potlist, cr=cr, alpha=alpha)
        plt.show()
        return
    
    def _plotPotTraj(self, tlist, potlist, ax=None, cr='r', alpha=0.4, lb=""):
        if not ax:
            fig, ax = plt.subplots(figsize=(8,3), dpi=100)
            setAxWidth(ax, 1.5)
            plt.xlabel("time, ms", fontsize=20)
            plt.ylabel("potential, kT", fontsize=20)
            
        if lb:
            ax.plot(tlist, potlist, '-', color=cr, alpha=alpha, label=lb)
        else:
            ax.plot(tlist, potlist, '-', color=cr, alpha=alpha)
        plt.scatter(tlist[0], potlist[0], color=cr, alpha=0.8)
        return ax
        
    def plotTimeDis(self, tList = None, ax=None, cr='b', alpha=0.5):
        if not ax:
            fig, ax = plt.subplots(figsize=(6,5))
            setAxWidth(ax, 1.5)
            plt.xlabel("lifetime, ms", fontsize=labelSize)
            plt.ylabel("counts", fontsize=labelSize)
        if not isinstance(tList, list):
            tList = self.model.tList.copy()
            for i in range(len(tList)):
                tList[i] = tList[i]*self.model.time_unit
        n, bins, patches = ax.hist(tList, 20, facecolor=cr, density=False,ec='white', alpha=alpha)
        
        return ax
    
    
    def plotForceDis(self, fList=None, ax=None, cr='b', alpha=1.0):
        if not ax:
            fig, ax = plt.subplots(figsize=(6,4))
            setAxWidth(ax, 1.5)
            plt.xlabel("rupture force, pN", fontsize=labelSize)
            plt.ylabel("count", fontsize=labelSize)
        if not isinstance(fList, list):
            fList = self.model.fList
        n, bins, patches = ax.hist(fList, 20, facecolor=cr, density=True,ec='white', alpha=alpha)
        
        return ax
    
    
    def plotPotential(self, x_traj=None, y_traj=None, overlay=True, plotTraj=True, title="", vector_field=False, pot_surface=True, colorful=True, traj_color="lime", sample_rate = 10, label_point=True, return_ax=False, figsize=(4,4)):
        
        
        
        ns = 100
        fig, ax = plt.subplots(figsize=figsize, dpi=100)
        
        setAxWidth(ax,  1.3)
        if self.coords == "x1x2":
            #plt.xlabel("$x_1-x_{10}$", fontsize=labelSize)
            #plt.ylabel("$x_2-x_{20}$", fontsize=labelSize)
            plt.xlabel("$x_a$", fontsize=labelSize)
            plt.ylabel("$x_b$", fontsize=labelSize)
        elif self.coords == "xy":
            plt.xlabel("$x-x_0$", fontsize=labelSize)
            plt.ylabel("$y-y_0$", fontsize=labelSize)
        xmin = -1.0*self.model.bd1.x1
        xmax = 1.5*self.model.bd1.x1
        
        if self.coords == "x1x2":
            ymin = -1.0*self.model.bd2.x1
            ymax = 1.5*self.model.bd2.x1
        elif self.coords == "xy":
            ymin = -1.0*(self.model.bd2.x1+self.model.bd1.x1)
            ymax = 1.5*(self.model.bd2.x1+self.model.bd1.x1)
        
        #ymin, ymax = -1, 2
        #xmin, xmax = -1, 2
        
        x1 = np.linspace(xmin, xmax, ns)
        x2 = np.linspace(ymin, ymax, ns)
        
        #plt.xlim(-2, 8)
        #plt.ylim(-4, 8)
        
        
        pot = np.zeros((ns, ns))
        for i in range(ns):
            for j in range(ns):
                curr_t = self.model._updateForce(self.model.t)
                if self.coords == "x1x2":
                    tmp = self.model._pot(x1[i], x2[j], curr_t)
                elif self.coords == "xy":
                    tmp = self.model._pot(x1[i], x2[j]-x1[i], curr_t)
                    
                if not tmp or tmp>21:
                    pot[j, i] = None
                elif tmp<-100:
                    pot[j, i] = None
                else:
                    pot[j, i] = tmp
        #
        if vector_field:
            nf = 20
            dx = 0.01
            x1f = np.linspace(xmin, xmax, nf)
            x2f = np.linspace(ymin, ymax, nf)
            uvec = np.zeros((nf, nf))
            vvec = np.zeros((nf, nf))
            for i in range(nf):
                for j in range(nf):
                    curr_t = self.model._updateForce(self.model.t)
                    if self.coords == "x1x2":
                        tmp1 = self.model._pot(x1f[i], x2f[j], curr_t)
                        tmp2 = self.model._pot(x1f[i]+dx, x2f[j], curr_t)
                        if type(tmp1)!=type(None) and type(tmp2)!=type(None):
                            u = -(tmp2-tmp1)/dx
                            tmp3 = self.model._pot(x1f[i], x2f[j]+dx, curr_t)
                        if type(tmp3)!=type(None) and type(tmp1)!=type(None):
                            v = -(tmp3-tmp1)/dx
                            r = np.sqrt(u**2+v**2)
                        uvec[j, i] = u/r
                        vvec[j, i] = v/r
                    elif self.coords == "xy":
                        tmp1 = self.model._pot(x1f[i], x2f[j]-x1f[i], curr_t)
                        tmp2 = self.model._pot(x1f[i]+dx, x2f[j]-x1f[i]-dx, curr_t)
                        if type(tmp1)!=type(None) and type(tmp2)!=type(None):
                            u = -(tmp2-tmp1)/dx
                            tmp3 = self.model._pot(x1f[i], x2f[j]+dx-x1f[i], curr_t)
                        if type(tmp3)!=type(None) and type(tmp1)!=type(None):
                            v = -(tmp3-tmp1)/dx
                            r = np.sqrt(u**2+v**2)
                        uvec[j, i] = u/r
                        vvec[j, i] = v/r
        #hmap=plt.pcolor(x1,x2,pot,cmap="bwr")
        #hmap=plt.contour(x1, x2, pot, levels=18, cmap="coolwarm")
        if pot_surface:
            hmap=plt.contourf(x1, x2, pot, levels=16, cmap="coolwarm")
            #cbar = plt.colorbar(hmap, ticks=[-])
            #cbar.ax.set_yticklabels(['< -1', '0', '> 1'])
        if self.coords == "x1x2":
            plt.axvline(x=self.model.bd1.x1, linestyle='--', color='k', lw=1.2)
            plt.axhline(y=self.model.bd2.x1, linestyle='--', color='k', lw=1.2)
        
        
            #plt.xticks([-2,-1, 0, 1, 2, 3], fontsize=13)
            #plt.yticks([-2, -1, 0, 1, 2, 3], fontsize=13)
            
        if vector_field:
            q = ax.quiver(x1f, x2f, uvec, vvec)
        #xb = [xi+self.bd2.x1 for xi in x1]
        #plt.plot(x1, xb, '-r', lw=2.0)
        #plt.xlim(xmin, xmax)
        #plt.ylim(ymin, ymax)
        
        
        if overlay:
            colors = cm.rainbow(np.linspace(0, 1, 11))
            #plt.scatter(self.model.x1, self.model.x2, color=colors, alpha=0.8)
            
        if label_point:
            ### label the position of attractor and saddle point
            xa_A, xb_A = self.model.get_attractor()
            xa_Sa, xb_Sa = self.model.get_saddle_a()
            xa_Sb, xb_Sb = self.model.get_saddle_b()
            
            if self.coords=="x1x2":
                plt.scatter([xa_A], [xb_A])
                plt.scatter([xa_Sa, xa_Sb], [xb_Sa, xb_Sb])
            else:
                plt.scatter([xa_A], [xa_A+xb_A])
                plt.scatter([xa_Sa, xa_Sb], [xa_Sa + xb_Sa, xa_Sb + xb_Sb])
            
        crange = 0.5
        
        if plotTraj:
            #ax.plot(self.model.x1_traj, self.model.x2_traj, alpha=0.7, color="green")
            if isinstance(x_traj, type(None)) or isinstance(y_traj, type(None)):
                if self.coords == "x1x2":
                    
                    x_traj = self.model.x1_traj
                    y_traj = self.model.x2_traj
                elif self.coords == "xy":
                    self.convertToXY()
                    x_traj = self.x_traj
                    y_traj = self.y_traj
            N = len(x_traj)//sample_rate
            #plt.scatter(self.model.x1, self.model.x2, color=colors, alpha=0.8)
            x_traj_sample =  [x_traj[0+i*sample_rate] for i in range(N)] + x_traj[-sample_rate:]
            y_traj_sample = [y_traj[0+i*sample_rate] for i in range(N)] + y_traj[-sample_rate:]
            if colorful:
                ## plot colorful traj
                for i in range(N-1):
                    #ax.plot(x_traj[i:i+2], y_traj[i:i+2],alpha=0.5, color='k')#  color=plt.cm.jet(i/N))
                    ## colorful plot
                    ax.plot(x_traj_sample[i:i+2], y_traj_sample[i:i+2],alpha=0.8,  color=plt.cm.jet(0.3+crange*i/N))
                    #ax.plot(x_traj[i:i+2], y_traj[i:i+2],  color='k', alpha=i/N)
            else:
                ax.plot(x_traj_sample, y_traj_sample, alpha=0.7, color=traj_color, linewidth=0.4)

            #### plotting deterministic trajectory
            if self.coords == "x1x2":
                #plt.plot(self.model.x1_traj, self.model.x2_traj,'--', color='k')
                pass
            else:
                y_traj_det = [self.model.x1_traj[i]+self.model.x2_traj[i] for i in range(len(self.model.x1_traj))]
                #plt.plot(self.model.x1_traj, y_traj_det,'--', color='k')
        if title != "":
            plt.savefig(title, fmt='pdf')
            plt.close()
        else:
            if return_ax:
                return ax, hmap
            else:
                plt.show()
        
        return None, None
    
    

        
    def lifetimeDis(self, ax=None, sufa=False,sufa_ratio=False, norm=False, cr='r', lb="", xvariable="time"):
        
        tend = np.mean(self.model.tList)
        tList = []
        
        tau_mean = np.mean(self.model.tList)

        if sufa_ratio:
            tList =[]
            print("sufa ratio mode")
            for i in range(len(self.model.tSuList)):
                tList.append(self.model.tSuList[i]/self.model.tFaList[i])
        else:
            if xvariable=="time":
                for tau in self.model.tList:
                    if norm:
                        tList.append(tau/tau_mean)
                    else:
                        tList.append(tau/100)
            elif xvariable=="force":
                tau_mean = np.mean(self.model.fList)
                print("fmean =", tau_mean)
                for tau in self.model.fList:
                    
                    if norm:
                        tList.append(tau/tau_mean)
                    else:
                        tList.append(tau)
            
            
        if sufa:
            print("tend\ttSu\ttFa")
            tSu = np.mean(self.model.tSuList)
            tFa = np.mean(self.model.tFaList)
            print("{0:.2f}\t{1:.2f}\t{2:.2f}".format(tend, tSu, tFa))
        else:
            print("{0:.2f}".format(tend))
            
            
        if not ax: 
            fig, ax = plt.subplots(dpi=100)
            if norm:
                if xvariable=="time":
                    plt.xlabel(r"$\tau/<\tau>$", fontsize=labelSize)
                elif xvariable=="force":
                    plt.xlabel(r"$F/<F>$", fontsize=labelSize)
                plt.ylabel("distribution", fontsize=labelSize)
            else:
                plt.xlabel("time, ms", fontsize=labelSize)
                #plt.xscale('log')
                plt.ylabel("counts", fontsize=labelSize)
            if sufa_ratio:
                plt.xlabel("tSu/tFa", fontsize=labelSize)
        
        
        
        if not sufa:
            if norm:
                bins = np.linspace(0, 5, 30)
                
                count,bins, = np.histogram(tList, bins=bins, density=True)
                plt.plot(bins[:-1], count,'-o', ms=5.0, fillstyle='none', color=cr, lw=1.5, label=lb)
                n, bins, patches = ax.hist(tList, bins=bins, color=cr,ec='white', alpha=0.4, label=lb, density=True)
            else:
                bins = np.logspace(-3, 3, 70)
                #n, bins, patches = ax.hist(tList, bins=bins, color=cr,ec='white', alpha=0.4, label=lb, density=False)
                n, bins, patches = ax.hist(tList, 40, color=cr,ec='white', alpha=0.4, label=lb, normed=False)
        else:
            n,bins, = np.histogram(self.model.tList, 30)
            plt.plot(bins[:-1], n, '-k', lw=3.0, label="rupture time")
            n, bins, patches = ax.hist(self.model.tSuList, 30, color='navy',ec='white', alpha=0.4, label="sucess time")
            n, bins, patches = ax.hist(self.model.tFaList, 30, color='r',ec='white', alpha=0.4, label="fail time")
        #plt.legend(loc=0, frameon=False, fontsize=labelSize)
        #plt.show()
        return ax
        
        
    
    def getPop(self, manyRun=True, sync=False, max_length=1000):
        apc_ag_bcr = []
        apc_ag = []
        ag_bcr = []
        f_traj = []
        if manyRun:
            ## get sample number
            num_sample = len(self.model.mol_rec_many)
            t_end_list = []
            for reci in self.model.mol_rec_many:
                t_end_list.append(len(reci))
            
            length = max(t_end_list)
            print("num_sample:{0:d}\tmax_length:{1:d}".format(num_sample, length))
            if length>max_length:
                step = length//max_length
            else:
                step = 1
            def get_sequence(k):
                lst = []
                for i in range(0, length, step):
                    tmp_list = []
                    for j in range(num_sample):
                        if i < t_end_list[j]:
                            tmp_list.append(self.model.mol_rec_many[j][i][k])
                        else:
                            if sync:
                                tmp_list.append(self.model.mol_rec_many[j][t_end_list[j]-1][k])
                            else:
                                tmp_list.append(None)
                    lst.append(tmp_list)
                return lst
            
            apc_ag_bcr = get_sequence(0)
            apc_ag = get_sequence(1)
            ag_bcr = get_sequence(2)
            t = [self.model.record_time*i*self.model.dt*step/1000 for i in range(len(apc_ag_bcr))]
            for ti in t:
                f_traj.append(1000*self.model._updateForce(1000*ti))
            
            return t, apc_ag_bcr, apc_ag, ag_bcr, f_traj
        else:
            length = len(self.model.mol_rec)
            if length>max_length:
                step = length//max_length
            else:
                step = 1
            for i in range(0, length, step):
                apc_ag_bcr.append(self.model.mol_rec[i][0])
                apc_ag.append(self.model.mol_rec[i][1])
                ag_bcr.append(self.model.mol_rec[i][2])
                f_traj.append(self.model.f_traj[i])
            t = [self.model.record_time*i*self.model.dt*step/1000 for i in range(len(apc_ag_bcr))]
            return t, apc_ag_bcr, apc_ag, ag_bcr, f_traj
    
    def plotPop(self, ax=None, manyRun=True):
        if ax==None:
            fig, ax = plt.subplots(dpi=100)
            plt.xlabel("time, s", fontsize=14)
            plt.ylabel("num of bonds", fontsize=14)
            
        t, apc_ag_bcr, apc_ag, ag_bcr, f_traj = self.getPop(manyRun)
        if manyRun:
            alpha0 = 0.3
        else:
            alpha0 = 1.0
        plt.plot(t, apc_ag_bcr, color='r', alpha=alpha0)
        plt.plot(t, apc_ag, color='b', alpha=alpha0)
        plt.plot(t, ag_bcr, color='darkgreen', alpha=alpha0)
        #plt.legend(["APC=Ag=BCR", "APC=Ag  BCR", "APC  Ag=BCR"], frameon=False)
        if not manyRun:
            ax2 = ax.twinx()
            ax2.plot(t, f_traj, 'k-', alpha=0.4)
        plt.show()
        return
    
    def plotMean(self, ax=None, fill_between=True, max_length=1000):
        if ax==None:
            fig, ax = plt.subplots(dpi=100)
            plt.xlabel("time, s", fontsize=14)
            plt.ylabel("num of bonds", fontsize=14)
        t, lst1, lst2, lst3, f = self.getPop(sync=True, max_length=max_length)
        length = len(t)
        lst1_mean = np.zeros(length)
        lst2_mean = np.zeros(length)
        lst3_mean = np.zeros(length)
        lst1_std = np.zeros(length)
        lst2_std = np.zeros(length)
        lst3_std = np.zeros(length)
        
        
        for i, item in enumerate(lst1):
            lst1_mean[i] = np.mean(item)
            lst1_std[i] = np.std(item)
        for i, item in enumerate(lst2):
            lst2_mean[i] = np.mean(item)
            lst2_std[i] = np.std(item)
        for i, item in enumerate(lst3):
            lst3_mean[i] = np.mean(item)
            lst3_std[i] = np.std(item)
            
        #plt.yscale('log')
        plt.plot(t, lst1_mean, color='r', lw=2)
        plt.plot(t, lst2_mean, color='b', lw=2)
        plt.plot(t, lst3_mean, color='darkgreen', lw=2)
        
        if fill_between:
            plt.fill_between(t, lst1_mean-lst1_std, lst1_mean+lst1_std, color="red", alpha=0.1)
            plt.fill_between(t, lst2_mean-lst2_std, lst2_mean+lst2_std, color="blue", alpha=0.1)
            plt.fill_between(t, lst3_mean-lst3_std, lst3_mean+lst3_std, color="green", alpha=0.1)
        #plt.legend(["APC=Ag=BCR", "APC=Ag  BCR", "APC  Ag=BCR"], frameon=False)
        
        ax2 = ax.twinx()
        ax2.plot(t, f, '-k', alpha=0.5)
        plt.show()
        return
    
    
    
    
def printProgress(n, N):
    percent = int(100.0*n/N)
    toPrint = "progress: "
    for i in range(percent//5):
        toPrint += '|'
    toPrint += "{:d}%".format(percent)
    print(toPrint, end='\r')
    return


def setAxWidth(ax, frame_width):
    ax.spines['top'].set_linewidth(frame_width)
    ax.spines['right'].set_linewidth(frame_width)
    ax.spines['bottom'].set_linewidth(frame_width)
    ax.spines['left'].set_linewidth(frame_width)
    return

def getAx(xlabel="time", ylabel="quantity", fontsize=15):
    fig, ax = plt.subplots(figsize=(6,5))
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    setAxWidth(ax, 2.0)
    return ax