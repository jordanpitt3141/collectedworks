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
       Hamilbycell.append(Hamiltonianacrosscell(x,h,u,epsilon,sigma,i,dx))
    return 0.5*sum1, Hamilbycell
    
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
    
    
#Matching Paper, we also see that this gives very close to analytic solution by Wolfram Alpha
"""
dxw = "12"
wdirord = "o3"

ll = int(dxw)
diffuse = 1

#dx = 10.0 / (2**ll)
dx = 0.1
l = 0.1
dt = l*dx
startx = -700
endx = 700.0 + dx
startt = 0.0
endt = 30.0+(dt*0.9)  
        
szoomx = startx
ezoomx = endx
        
g = 9.81
    
gap = int(0.5/dt)
        
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
hf = 1.8
hl = 1.0
base = hl
eta0 = hf - hl
x0 = 500
#hm,um = dambreaksmooth(x,x0,base,eta0,diffuse,dx) 

eta0 = 0.4182
hm,um = dambreaksmoothHamilcomp(x,eta0,dx)  

epsilon = hf/hl
sigma = g 

#BC
nBC = 2
hbeg = hm[0]*ones(nBC)
hend = hm[-1]*ones(nBC)
ubeg = um[0]*ones(nBC)
uend = um[-1]*ones(nBC)

xbeg = arange(startx - nBC*dx,startx,dx)
xend = arange(endx + dx,endx + (nBC+1)*dx)

xbc = concatenate([xbeg,x,xend])
hbc = concatenate([hbeg,hm,hend])
ubc = concatenate([ubeg,um,uend])     


Hbeg = Hamiltonianall(xbc,hbc,ubc,1.0,1.0,nBC,dx)
"""



#Nonlinear Soliton

# We conserve the Hamiltonian for analytic solution up to 10**-11 accuracy which is floating point for this many points
stx = -600
etx = 600

g = 1.0
dxw = "10"
wdirord = "o3"

    
wdir = "../../../../data/raw/Joesolconnon/"  +wdirord +"/" + dxw + "/"
sdir = "../../../../data/postprocessing/Joesolconnon/"  +wdirord +"/" + dxw + "/"
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
            #h.append(float(row[3]))
            #u.append(float(row[5]))
            h.append(float(row[6]))
            u.append(float(row[7]))
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
Eend = HankEnergyall(xbc,hbc,ubc,g,nBC,dx)

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
            #h.append(float(row[3]))
            #u.append(float(row[5]))
            h.append(float(row[6]))
            u.append(float(row[7]))
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
Ebeg = HankEnergyall(xbc,hbc,ubc,g,nBC,dx)

#legend()


"""
# Soliton Interaction

startx = -200
endx = 200

dxw = "15"
wdirord = "o3"

epsilon = 1.0
sigma = 1.0

    
wdir = "../../../../data/raw/Cserre/solitonothers/collDMcopy/"  +wdirord +"/"
#sdir = "../../../../data/postprocessing/Joesolconnon/"  +wdirord +"/" + dxw + "/"
#if not os.path.exists(sdir):
#        os.makedirs(sdir)
     
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
            x.append(float(row[3]))
            h.append(float(row[4]))
            u.append(float(row[6]))
         j = j + 1
     x = array(x)
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
xbeg = arange(startx - nBC*dx,startx,dx)
xend = arange(endx + dx,endx + (nBC+1)*dx)

xbc = concatenate([xbeg,x,xend])
hbc = concatenate([hbeg,h,hend])
ubc = concatenate([ubeg,u,uend])

Hend = Hamiltonianall(xbc,hbc,ubc,epsilon,sigma,nBC,dx)

s = wdir + "saveoutputts1.txt"
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
  
plot(x,h)
xlabel("$x$ ($m$)")
ylabel("$h$ ($m$)")

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

Hbeg = Hamiltonianall(xbc,hbc,ubc,epsilon,sigma,nBC,dx)
"""


# Smooth Dambreak
"""
startx = 0.0
endx = 1000

dxw = "10"
wdirord = "o3"
diffusew = "20"
diffuses = [0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1.0,2.5,5.0,7.5,10.0,25.0,50.0,75.0,100.0,250.0,500.0,750.0,1000.0]

epsilon = 0.8
sigma = 9.81
g = 9.81

    
wdir = "../../../../data/raw/Joebigsmooth/"  +wdirord +"/" + dxw + "/" + diffusew + "/"
#sdir = "../../../../data/postprocessing/Joesolconnon/"  +wdirord +"/" + dxw + "/"
#if not os.path.exists(sdir):
#        os.makedirs(sdir)
    
s = wdir + "outlast.txt"
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
            alpha = float(row[7])
         j = j + 1
     x = array(x)
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
xbeg = arange(startx - nBC*dx,startx,dx)
xend = arange(endx + dx,endx + (nBC+1)*dx)

xbc = concatenate([xbeg,x,xend])
hbc = concatenate([hbeg,h,hend])
ubc = concatenate([ubeg,u,uend])

Hend,Henda = Hamiltonianall(xbc,hbc,ubc,epsilon,sigma,nBC,dx)
Eend = HankEnergyall(xbc,hbc,ubc,g,nBC,dx)

#Have to start it up initial conditions
ll = int(dxw)
dx = 10.0 / (2**ll)
l = 0.01
dt = l*dx
startx = 0.0
endx = 1000.0 + dx
startt = 0.0
endt = 30.0+(dt*0.9)  
        
szoomx = startx
ezoomx = endx
                
g = 9.81
       
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
hf = 1.8
hl = 1.0
base = hl
eta0 = hf - hl
x0 = 500
hi,ui = dambreaksmooth(x,x0,base,eta0,alpha,dx)    
  
plot(x,hi)
xlabel("$x$ ($m$)")
ylabel("$h$ ($m$)")

#BC
nBC = 2
hbeg = hi[0]*ones(nBC)
hend = hi[-1]*ones(nBC)
ubeg = ui[0]*ones(nBC)
uend = ui[-1]*ones(nBC)

xbeg = arange(startx - nBC*dx,startx,dx)
xend = arange(endx + dx,endx + (nBC+1)*dx)

xbc = concatenate([xbeg,x,xend])

hibc = concatenate([hbeg,hi,hend])
uibc = concatenate([ubeg,ui,uend])

Hbeg, Hbega = Hamiltonianall(xbc,hibc,uibc,epsilon,sigma,nBC,dx)
Ebeg = HankEnergyall(xbc,hibc,uibc,g,nBC,dx)
"""