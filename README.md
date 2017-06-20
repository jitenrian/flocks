# flocks
#Asymmetric exchange of information in flocks
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 02:15:42 2017

@author: jitendrian
"""
import csv
import numpy as np
import matplotlib.pyplot as plt
import os

"*****************************************************************************"
#distance between ath bird and bth bird
def dist_ab(rax,ray,rbx,rby):
        return np.sqrt(np.power(rbx-rax,2)+np.power(rby-ray,2))
    
"*****************************************************************************"    
#projection of distance between ath bird and bth bird on x-axis
def xy_ab(rax,rbx):
    return (rax-rbx)

"*****************************************************************************"    
#flock behaviour
def exchange(J, A, chi, eta, dt, T, N, L, IR):
    G = int(T/dt)
    
    # J,A = exchange interaction
    # chi = spin inertia
    # eta = 
    # dt = grid size of time
    # T = total time
    # N = flock size
    # L = length of a square box
    # IR = interaction radius
    # G =  total grids in T time
    
        
    ## make arrays of size N.
    vxi = np.zeros(N)
    vyi = np.zeros(N)
    vxf = np.zeros(N)
    vyf = np.zeros(N)
    rxi = np.zeros(N)
    ryi = np.zeros(N)
    rxf = np.zeros(N)
    ryf = np.zeros(N)
    thei = np.zeros(N)
    thef = np.zeros(N)
    si = np.zeros(N)
    sf = np.zeros(N)
    # vxi = x-component of initial velocity
    # vyi = y-component of initial velocity
    # vxf = x-component of next velocity
    # vyf = y-component of next velocity
    # rxi = x-component of initial position
    # ryi = y-component of initial position
    # rxf = x-component of next position
    # ryf = y-component of next position
    # thi = initial theta
    # thf = next theta
    # si = initial spin
    # sf = spin in next instant
    #images
    "***********************************************************************"
    #exchange interaction>summation over all b's in neighbour of a (J_ab*(v_aXv_b))
    def interaction(a):
        sJ_ab = 0
        for i in range(N):
            if i!=a:
                d = dist_ab(rxi[a],ryi[a],rxi[i],ryi[i])
                if d < IR:
                    phi = ((xy_ab(rxi[a],rxi[i])*np.cos(thei[a])) - (xy_ab(ryi[a],ryi[i])*np.sin(thei[a])))/d
                    sJ_ab = sJ_ab + (J - phi*A)*np.sin(thei[i]-thei[a])
        return sJ_ab
    
    "*************************************************************************"

    "*************************************************************************"
    "             initiation of velocity, spin and position                   "
    "*************************************************************************"
    for t in range(N):
        thei[t] = np.pi/2 + (np.sin(t*np.pi/10)/20)
        vxi[t] = np.cos(thei[t])
        vyi[t] = np.sin(thei[t])
        si[t] = np.cos(t*np.pi/10)/20
    
    q = 0
    for t in range(L):
        for w in range(L):
            rxi[q], ryi[q] = t, w
            q = q+1
            
    "*************************************************************************"
    "                              saving location                            "
    "*************************************************************************"         
    Z = eta*np.sqrt(J/chi)
    p_path = os.getcwd()
    folder = 'J_' + str(J) + '_A_' + str(A) + '_chi_' + str(chi) + '_eta_' + str(eta) + '_mod_A_' + str(Z)
    os.mkdir(folder)
    n_path = os.path.join(p_path,folder)
    os.chdir(n_path)
    "*************************************************************************"
    "                             time evolution                              "
    "*************************************************************************"

    for g in range(G):
        
        for t in range(N):
            sf[t] = si[t] + dt*(interaction(t) - ((eta/chi)*si[t]))
            thef[t] = thei[t] + (si[t]/chi)*dt
            vxf[t] = np.cos(thef[t])
            vyf[t] = np.sin(thef[t])
            rxf[t] = rxi[t] + vxf[t]*dt
            ryf[t] = ryi[t] + vyf[t]*dt
        "*********************************************************************"
        "                               periodicity                           "
        "*********************************************************************"
        for t in range(N):
            if rxf[t]<0:
                rxf[t] = rxf[t] + L
            if ryf[t]<0:
                ryf[t] = ryf[t] + L
            if rxf[t]>L:
                rxf[t] = rxf[t] - L
            if ryf[t]>L:
                ryf[t] = ryf[t] - L
               
        for t in range(N):
            vxi[t] = vxf[t]
            vyi[t] = vyf[t]
            rxi[t] = rxf[t]
            ryi[t] = ryf[t]
            si[t] = sf[t]
            thei[t] = thef[t]
            
        dr = zip(rxi,ryi,vxi,vyi)
        
        name = 'test'+str(g)+'.csv'

        with open(name, 'wb') as fp:
            a = csv.writer(fp)
            a.writerow(['rxi','ryi','vxi','vyi'])
            for x in dr:
                a.writerow(x)
        fp.close()
        G = '%04d' %g 
        plt.figure()    
        plt.plot(rxf,ryf,'.')
        plt.savefig('image'+str(G)+'.png')
        plt.close()
    
    os.chdir(p_path)

if __name__=="__main__":
    J = 2.49
    chi = 4.01
    eta = 0.81
    dt = 0.05
    T = 10.0
    N = 900
    L = 30
    IR = 2.0
    A_c = eta*np.sqrt(J/chi)
    AC = [A_c*2.5,A_c*3.5,A_c*4.5]
    
    for A in AC:
        exchange(J, A, chi, eta, dt, T, N, L, IR)       
