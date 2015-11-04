# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 23:36:31 2015

@author: jordan
"""
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

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
   
dx = 10.0 / (2**12)
l = 0.01
dt = l*dx
theta = 1.2
startx = 0.0
endx = 1000.0 + dx
startt = 0.0
endt = dt  
    
g = 9.81
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
base = 1.0
eta0 = 0.8
x0 = 500

diffuses = [0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1.0,2.5,5.0,7.5,10.0,25.0,50.0,75.0,100.0,250.0,500.0,750.0,1000.0]

for i in range(12,len(diffuses)):
    diffuse = diffuses[i]
    hh,uu = dambreaksmooth(x,x0,base,eta0,diffuse,dx)
    s = "diff = " + str(diffuse)
    plot(x,hh,label=s)   

s = "Initial Conditions Dam Break: dx = " + str(dx)
title(s)
xlabel("x (m)")
ylabel("h (m)")
legend()


"""
sdir = "../../results/show/steve2/initial/"

for ll in range(3,16):
    dx = 10.0 / (2**ll)
    l = 0.01
    dt = l*dx
    theta = 1.2
    startx = 0.0
    endx = 1000.0 + dx
    startt = 0.0
    endt = dt  
        
    g = 9.81
        
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
    base = 1.0
    eta0 = 0.8
    x0 = 500
    
    diffuses = [0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1.0,2.5,5.0,7.5,10.0,25.0,50.0,75.0,100.0,250.0,500.0,750.0,1000.0]
    
    for i in range(len(diffuses)):
        diffuse = diffuses[i]
        hh,uu = dambreaksmooth(x,x0,base,eta0,diffuse,dx)
        s = "diff = " + str(diffuse)
        plot(x,hh,label=s)   
    
    s = "Initial Conditions Dam Break: dx = " + str(dx)
    title(s)
    xlabel("x (m)")
    ylabel("h (m)")
    legend()
    s = sdir +"diffalldx"+ str(dx) +".png"       
    savefig(s, bbox_inches='tight')        
    clf()
"""    