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
from numpy import ones

from numpy import tanh

def Minitial(xbeg,xend):
    A = xbeg + 2*sqrt(2.0/3.0)*tanh(sqrt(3.0/8.0)*xbeg)
    B = xend + 2*sqrt(2.0/3.0)*tanh(sqrt(3.0/8.0)*xend)
    
    return B - A   
    
    
def Pinitial(xbeg,xend):
    A = 7.23326*tanh(0.612372*xbeg)
    B = 7.23326*tanh(0.612372*xend)
    
    return B - A  
    
def uhacrosscell(x,h,u,j,dx): 
    return dx*h[j]*u[j]
    
def Hamiltonianall(x,h,u,dx):
    n = len(x)
    sum1 = 0.0
    for i in range(n):
        sum1 = sum1 + uhacrosscell(x,h,u,i,dx)
    return  sum1
    


def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 


def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a

def soliton (x,t,g,a0,a1):
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)
  
def solitoninit(n,a0,a1,g,x,t0,dx):
    h = zeros(n)
    u = zeros(n)
    c = sqrt(g*(a0 + a1))
    for i in range(n):
        h[i] = soliton(x[i],t0,g,a0,a1)
        u[i] =  c* ((h[i] - a0) / h[i])
    
    return h,u

"""
#Nonlinear Soliton
dxw = "12"
wdirord = "FDcent"

    
wdir = "../../../../data/raw/solconnonsmallg10FDall/"  +wdirord +"/" + dxw + "/"
sdir = "../../../../data/postprocessing/solconsmallg10FD/"  +wdirord +"/" + dxw + "/"
if not os.path.exists(sdir):
        os.makedirs(sdir)
     
s = wdir + "outlast.txt"
with open(s,'r') as file1:
     readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     h = []
     u = []
     he = []
     ue = []
     x = []
     j = -1
     for row in readfile:       
         if (j >= 0):
            dx = float(row[0])
            dt = float(row[1])
            t = float(row[2])
            #h.append(float(row[3]))
            x.append(float(row[6]))
            h.append(float(row[7]))
            u.append(float(row[8]))
            he.append(float(row[9]))
            ue.append(float(row[10]))
         j = j + 1    
  
plot(x,h)
plot(x,he)
n = len(x)
o = ones(n)

hi,ui = solitoninit(n,1,1,9.81,x,0,dx)

Ma = Hamiltonianall(x,he,o,dx)
Pa = Hamiltonianall(x,he,ue,dx)

Mn = Hamiltonianall(x,h,o,dx)
Pn = Hamiltonianall(x,h,u,dx)

Mi = Hamiltonianall(x,hi,o,dx)
Pi = Hamiltonianall(x,hi,ui,dx)
"""

wdirord = "grim"
#wdir = "../../../../data/raw/NEWdata/FDredo/grim/"
#sdir = "../../../../data/postprocessing/scFDallAE/grim/"

sdir = "../../../../data/postprocessing/scFDallAC/"+wdirord+"/"
if not os.path.exists(sdir):
        os.makedirs(sdir)
Mns = []
Pns = []
Mis = []
Pis = []
dxs = []
#range(6,19)
for ki in range(6,21):
    
    #Nonlinear Soliton
    dxw = str(ki)
    
        
    wdir = "../../../../data/raw/NEWdata/FDredo/"  +wdirord +"/" + dxw + "/"
         
    s = wdir + "outlast.txt"
    with open(s,'r') as file1:
         readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
         h = []
         u = []
         he = []
         ue = []
         x = []
         j = -1
         for row in readfile:       
             if (j >= 0):
                dx = float(row[0])
                dt = float(row[1])
                t = float(row[2])
                #h.append(float(row[3]))
                x.append(float(row[6]))
                h.append(float(row[7]))
                u.append(float(row[8]))
                he.append(float(row[9]))
                ue.append(float(row[10]))
             j = j + 1    
      
    #plot(x,h)
    #plot(x,he)
    n = len(x)
    o = ones(n)
    
    #hi,ui = solitoninit(n,1,1,9.81,x,0,dx)
    
    Mn = Hamiltonianall(x,h,o,dx)
    Pn = Hamiltonianall(x,h,u,dx)
    
    xbeg = -50 - 0.5*dx
    xend = 250 + 0.5*dx

    Pi = Pinitial(xbeg,xend)
    Mi = Minitial(xbeg,xend)
    
    #Mi = Hamiltonianall(x,hi,o,dx)
    #Pi = Hamiltonianall(x,hi,ui,dx)
    
    Pns.append(Pn)
    Pis.append(Pi)
    Mns.append(Mn)
    Mis.append(Mi)
    dxs.append(dx)

    
Mns = array(Mns)
Pns = array(Pns)
Mis = array(Mis)
Pis = array(Pis)

relerrP = abs(Pis - Pns)/ abs(Pis)
relerrM = abs(Mis - Mns)/ abs(Mis)

n= len(dxs)
s = sdir + "con.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
     writefile2.writerow(['dx','Pis','Pns' ,'relerrP','Mis','Mns','relerrM'])        
                       
     for j in range(n):
         writefile2.writerow([str(dxs[j]),str(Pis[j]),str(Pns[j]),str(relerrP[j]),str(Mis[j]),str(Mns[j]),str(relerrM[j])])  

s = sdir + "conh.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",relerrM[i])
        file1.write(s)         
    
s = sdir + "conuh.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",relerrP[i])
        file1.write(s)   
