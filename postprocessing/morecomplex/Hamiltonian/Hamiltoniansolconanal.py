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
    
    fgpe = epsilon*fgph*fgpu*fgpu + i3*(epsilon*sigma*sigma)*fgph*fgph*fgph*fgpux*fgpux \
        + ie*fgph*fgph
        
    #second gauss point
    sgp = x[j]
    sgph = interpquarticval(haj,hbj,hcj,hdj,hej,x[j],sgp)
    sgpu = interpquarticval(uaj,ubj,ucj,udj,uej,x[j],sgp)
    sgpux = interpquarticgrad(uaj,ubj,ucj,udj,uej,x[j],sgp)
    
    sgpe = epsilon*sgph*sgpu*sgpu + i3*(epsilon*sigma*sigma)*sgph*sgph*sgph*sgpux*sgpux \
        + ie*sgph*sgph

    #third gauss point
    tgp = -0.5*dx*sqrt(3.0/5.0) + x[j]
    tgph = interpquarticval(haj,hbj,hcj,hdj,hej,x[j],tgp)
    tgpu = interpquarticval(uaj,ubj,ucj,udj,uej,x[j],tgp)
    tgpux = interpquarticgrad(uaj,ubj,ucj,udj,uej,x[j],tgp)
    
    tgpe = epsilon*tgph*tgpu*tgpu + i3*(epsilon*sigma*sigma)*tgph*tgph*tgph*tgpux*tgpux \
        + ie*tgph*tgph  
        
    Hamilcell = 0.5*dx*( (5.0/9.0)*fgpe + (8.0/9.0)*sgpe + (5.0/9.0)*tgpe)
    
    return Hamilcell
    
def Hamiltonianall(x,h,u,epsilon,sigma,nBC,dx):

    n = len(x)
    sum1 = 0.0
    for i in range(nBC,n - nBC):
       sum1 = sum1 + Hamiltonianacrosscell(x,h,u,epsilon,sigma,i,dx)
    return 0.5*sum1
        
        

dxw = "15"
wdirord = "o3"

stx = -50
etx = 250

    
#wdir = "../../../../data/raw/Joesolconnon/"  +wdirord +"/" + dxw + "/"
#sdir = "../../../../data/postprocessing/Joesolconnon/"  +wdirord +"/" + dxw + "/"
wdir = "../../../../data/raw/hinonling10/"  +wdirord +"/" + dxw + "/"
sdir = "../../../../data/postprocessing/hinonling10/"  +wdirord +"/" + dxw + "/"
if not os.path.exists(sdir):
        os.makedirs(sdir)
     
s = wdir + "saveoutputtslast.txt"
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
            h.append(float(row[3]))
            u.append(float(row[5]))
         j = j + 1
     x = arange(stx,etx+dx,dx)#array(xt)
     h = array(h)     
  
plot(x,h)
xlabel("$x$ ($m$)")
ylabel("$h$ ($m$)")

#BC
nBC = 2
hbeg = h[0]*ones(nBC)
hend = h[-1]*ones(nBC)
ubeg = u[0]*ones(nBC)
uend = u[-1]*ones(nBC)

xbc = arange(stx-nBC*dx,etx+(nBC+1)*dx,dx)
hbc = concatenate([hbeg,h,hend])
ubc = concatenate([ubeg,u,uend])

Hend = Hamiltonianall(xbc,hbc,ubc,1.0,9.81,nBC,dx)

s = wdir + "saveoutputts1.txt"
with open(s,'r') as file1:
     readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     h = []
     u = []
     x = []
     j = 0
     for row in readfile:       
         if (j >= 0):
            dx = float(row[0])
            dt = float(row[1])
            t = float(row[2])
            h.append(float(row[3]))
            u.append(float(row[5]))
         j = j + 1
     x = arange(stx,etx+dx,dx)#array(xt)
     h = array(h)     
  
plot(x,h)
xlabel("$x$ ($m$)")
ylabel("$h$ ($m$)")

#BC
nBC = 2
hbeg = h[0]*ones(nBC)
hend = h[-1]*ones(nBC)
ubeg = u[0]*ones(nBC)
uend = u[-1]*ones(nBC)

xbc = arange(stx-nBC*dx,etx+(nBC+1)*dx,dx)
hbc = concatenate([hbeg,h,hend])
ubc = concatenate([ubeg,u,uend])

Hbeg = Hamiltonianall(xbc,hbc,ubc,1.0,9.81,nBC,dx)

#legend()


    

