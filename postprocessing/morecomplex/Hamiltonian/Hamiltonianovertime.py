# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 16:53:37 2015

@author: jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from matplotlib2tikz import save as tikz_save
import os

def interpquarticval(aj,bj,cj,dj,ej,xj,x):
    
    return aj*(x -xj)*(x -xj)*(x -xj)*(x -xj) + bj*(x -xj)*(x -xj)*(x -xj) \
    + cj*(x -xj)*(x -xj) + dj*(x -xj)+ ej
    
def interpquarticgrad(aj,bj,cj,dj,ej,xj,x):
    
    return 4*aj*(x -xj)*(x -xj)*(x -xj) + 3*bj*(x -xj)*(x -xj) \
    + 2*cj*(x -xj) + dj
    
def interpquartcoeff(q,j,dx):
    i24 = 1.0 / 24.0
    i12 = 1.0 / 12.0
    idx = 1.0/dx
    aj = i24*idx*idx*idx*idx*(q[j+2] - 4*q[j+1] + 6*q[j] - 4*q[j-1] + q[j-2])
    bj = i12*idx*idx*idx*(q[j+2] - 2*q[j+1] + 2*q[j-1] - q[j-2])
    cj = i24*idx*idx*(-q[j+2] + 16*q[j+1] - 30*q[j] + 16*q[j-1] - q[j-2])
    dj = i12*idx*(-q[j+2] + 8*q[j+1] - 8*q[j-1] + q[j-2])
    ej = q[j]
    
    return aj,bj,cj,dj,ej

def Hamiltonianacrosscell(x,h,u,epsilon,sigma,j,dx):
    #so we have h,u at midpoints
    #epsilon and sigma are everywhere
    i3 = 1.0/3.0
    ie = 1.0 /epsilon

    #jth cell
    uaj,ubj,ucj,udj,uej = interpquartcoeff(u,j,dx)
    haj,hbj,hcj,hdj,hej = interpquartcoeff(h,j,dx)
    
    #first gauss point
    fgp = 0.5*dx*sqrt(3.0/5.0) + x[j]
    fgph = interpquarticval(haj,hbj,hcj,hdj,hej,x[j],fgp)
    fgpu = interpquarticval(uaj,ubj,ucj,udj,uej,x[j],fgp)
    fgpux = interpquarticgrad(uaj,ubj,ucj,udj,uej,x[j],fgp)
    
    fgpe = epsilon*epsilon*fgph*fgpu*fgpu + i3*(epsilon*epsilon*sigma*sigma)*fgph*fgph*fgph*fgpux*fgpux \
        + (fgph-1)*(fgph-1)
        
    #second gauss point
    sgp = x[j]
    sgph = interpquarticval(haj,hbj,hcj,hdj,hej,x[j],sgp)
    sgpu = interpquarticval(uaj,ubj,ucj,udj,uej,x[j],sgp)
    sgpux = interpquarticgrad(uaj,ubj,ucj,udj,uej,x[j],sgp)
    
    sgpe = epsilon*epsilon*sgph*sgpu*sgpu + i3*(epsilon*epsilon*sigma*sigma)*sgph*sgph*sgph*sgpux*sgpux \
        + (sgph-1)*(sgph-1)

    #third gauss point
    tgp = -0.5*dx*sqrt(3.0/5.0) + x[j]
    tgph = interpquarticval(haj,hbj,hcj,hdj,hej,x[j],tgp)
    tgpu = interpquarticval(uaj,ubj,ucj,udj,uej,x[j],tgp)
    tgpux = interpquarticgrad(uaj,ubj,ucj,udj,uej,x[j],tgp)
    
    tgpe = epsilon*epsilon*tgph*tgpu*tgpu + i3*(epsilon*epsilon*sigma*sigma)*tgph*tgph*tgph*tgpux*tgpux \
        + (tgph-1)*(tgph-1) 
        
    Hamilcell = 0.5*dx*( (5.0/9.0)*fgpe + (8.0/9.0)*sgpe + (5.0/9.0)*tgpe)
    
    return Hamilcell
    
def HankEnergyacrosscell(x,h,u,g,j,dx):
    #so we have h,u at midpoints
    #epsilon and sigma are everywhere
    i3 = 1.0/3.0

    #jth cell
    uaj,ubj,ucj,udj,uej = interpquartcoeff(u,j,dx)
    haj,hbj,hcj,hdj,hej = interpquartcoeff(h,j,dx)
    
    #first gauss point
    fgp = 0.5*dx*sqrt(3.0/5.0) + x[j]
    fgph = interpquarticval(haj,hbj,hcj,hdj,hej,x[j],fgp)
    fgpu = interpquarticval(uaj,ubj,ucj,udj,uej,x[j],fgp)
    fgpux = interpquarticgrad(uaj,ubj,ucj,udj,uej,x[j],fgp)
    
    fgpe = fgph*fgpu*fgpu + g*fgph*fgph + i3*(fgph*fgph*fgph)*fgpux*fgpux
        
    #second gauss point
    sgp = x[j]
    sgph = interpquarticval(haj,hbj,hcj,hdj,hej,x[j],sgp)
    sgpu = interpquarticval(uaj,ubj,ucj,udj,uej,x[j],sgp)
    sgpux = interpquarticgrad(uaj,ubj,ucj,udj,uej,x[j],sgp)
    
    sgpe = sgph*sgpu*sgpu + g*sgph*sgph + i3*(sgph*sgph*sgph)*sgpux*sgpux

    #third gauss point
    tgp = -0.5*dx*sqrt(3.0/5.0) + x[j]
    tgph = interpquarticval(haj,hbj,hcj,hdj,hej,x[j],tgp)
    tgpu = interpquarticval(uaj,ubj,ucj,udj,uej,x[j],tgp)
    tgpux = interpquarticgrad(uaj,ubj,ucj,udj,uej,x[j],tgp)
    
    tgpe = tgph*tgpu*tgpu + g*tgph*tgph + i3*(tgph*tgph*tgph)*tgpux*tgpux
        
    Hamilcell = 0.5*dx*( (5.0/9.0)*fgpe + (8.0/9.0)*sgpe + (5.0/9.0)*tgpe)
    
    return Hamilcell
    
def Hamiltonianall(x,h,u,epsilon,sigma,nBC,dx):

    n = len(x)
    sum1 = 0.0
    Hamilbycell = []
    for i in range(nBC,n - nBC):
       sum1 = sum1 + Hamiltonianacrosscell(x,h,u,epsilon,sigma,i,dx)
       #Hamilbycell.append(Hamiltonianacrosscell(x,h,u,epsilon,sigma,i,dx))
    return 0.5*sum1#, Hamilbycell
    
def HankEnergyall(x,h,u,g,nBC,dx):

    n = len(x)
    sum1 = 0.0
    #Hamilbycell = []
    for i in range(nBC,n - nBC):
       sum1 = sum1 + HankEnergyacrosscell(x,h,u,g,i,dx)
       #Hamilbycell.append(HankEnergyacrosscell(x,h,u,g,i,dx))
    return 0.5*sum1#, Hamilbycell

def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 

def dambreaksmooth(x,x0,base,eta0,diffuse,dx):
    from numpy import tanh
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    
    for i in range(n):
        h[i] = base + 0.5*eta0*(1 + tanh(diffuse*(x0 - abs(x[i]))))

    return h,u  

def dambreaksmoothHamilcomp(x,eta0,dx):
    from numpy import tanh
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    
    for i in range(n):
        h[i] =0.5*eta0*(1 + tanh(250 - abs(x[i])))
    return h,u      


# Smooth Dambreak

startx = 0.0
endx = 1000

g = 9.81
epsilon = 1.0
sigma = 1.0

    
wdir = "../../../../data/raw/Cserre/dbE/o3/hf1p8hl1p0/"
sdir = "../../../../data/postprocessing/dbE/o3/hf1p8hl1p0/"
if not os.path.exists(sdir):
        os.makedirs(sdir)
        
filesi = os.listdir(wdir)

times = []
Hvals = []
Evals = []

for filen in filesi:
    
    s = wdir + filen
    with open(s,'r') as file1:
         readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
         h = []
         u = []
         x = []
         j = -1
         for row in readfile:       
             if (j >= 0):
                dx = float(row[0])
                dt = float(row[1])
                t = float(row[2])
                x.append(float(row[3]))
                h.append(float(row[4]))
                u.append(float(row[6]))
             j = j + 1
         x = array(x)
         h = array(h)     
  
    #plot(x,h)
    #xlabel("$x$ ($m$)")
    #ylabel("$h$ ($m$)")
    
    
    #BC
    nBC = 2
    hbeg = h[0]*ones(nBC)
    hend = h[-1]*ones(nBC)
    ubeg = u[0]*ones(nBC)
    uend = u[-1]*ones(nBC)
    xbeg = arange(startx - nBC*dx,startx,dx)
    xend = arange(endx + dx,endx + (nBC+1)*dx)
    
    xbc = concatenate([xbeg,x,xend])
    hbc = concatenate([hbeg,h,hend])
    ubc = concatenate([ubeg,u,uend])
    
    Hval = Hamiltonianall(xbc,hbc,ubc,epsilon,sigma,nBC,dx)
    Eval = HankEnergyall(xbc,hbc,ubc,g,nBC,dx)
    
    times.append(t)
    Hvals.append(Hval)
    Evals.append(Eval)

mix = zip(times,Evals,Hvals)
mix = sorted(mix)
timess = [p[0] for p in mix]
Evalss = [p[1] for p in mix]
Hvalss = [p[2] for p in mix]


s = sdir + "timeenergy.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     writefile2.writerow(['dx' ,'dt','time','Evals', 'Hvals'])        
                   
     for j in range(n):
         writefile2.writerow([str(dx),str(dt),str(times[i]), str(Evalss[j]), str(Hvals[j])]) 