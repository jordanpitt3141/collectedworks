# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 16:53:37 2015

@author: jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
import os
from numpy import ones

from numpy import tanh
from Hamil import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

def copyarraytoC(a):
    n = len(a)
    b = mallocPy(n)
    for i in range(n):
        writetomem(b,i,a[i])
    return b
    
def copyarrayfromC(a,n):
    b = [0]*n
    for i in range(n):
        b[i] = readfrommem(a,i)
        
    return b

def totalmassa(xb,xe,x0,h0,h1,alpha):
    return 0.5*(h1 + h0)*(xe - xb)

def totalmomentum():
    return 0

def Hamilana(xb,xe,x0,h0,h1,alpha,g):
    p1 = (g/8.0)*(xe - xb)*( ((h0+h1)**2) + ((h1 - h0)**2))
    p2 = alpha*(g/4.0)*(h1 - h0)**2*(tanh(0.5*(xb - xe)/alpha))
    
    return p1 + p2
        

def midpointtoca(h,dx):
    n = len(h)
    b = zeros(n)
    c = zeros(n)
    a = zeros(n)
    i24 = 1.0/24

    for i in range(n): 
        a[i-1] = -i24
        b[i] = 26*i24
        c[i] = -i24
        
    #i =0
    i = 0;
    b[i] = 1.0;
    c[i] = 0.0;

    #i=n-1
    i = n-1;
    a[i-1] = 0.0;
    b[i] = 1.0; 
    
    return TDMApy(a,b,c,h)
    
def TDMApy(a,b,c,d):
    n = len(d)
    alpha = []
    beta = []
    x = [0]*n
    
    alpha.append((1.0*c[0])/b[0])
    beta.append((1.0*d[0])/b[0] )  
 
    for i in range(1,n-1):
        m = 1.0 / (b[i] - a[i-1]*alpha[i-1])
        alpha.append(c[i]* m)
        beta.append((d[i] - a[i-1]*beta[i-1]) * m)
        
    m = 1.0 / (b[n-1] - a[n-2]*alpha[n-2])
    beta.append((d[n-1] - a[n-2]*beta[n-2]) * m)  

    x[n-1] = beta[n-1]
    
    for i in range(n-2,-1,-1):
        x[i] = beta[i] - alpha[i]*x[i+1]
 
    return array(x)
    
def uhacrosscell(x,h,u,j,dx): 
    return dx*h[j]*u[j]
    
def hacrosscell(x,h,j,dx): 
    return dx*h[j]
    
def uhall(x,h,u,dx):
    n = len(x)
    sum1 = 0.0
    for i in range(n):
        sum1 = sum1 + uhacrosscell(x,h,u,i,dx)
    return  sum1
    
def hall(x,h,dx):
    n = len(x)
    sum1 = 0.0
    for i in range(n):
        sum1 = sum1 + hacrosscell(x,h,i,dx)
    return  sum1

#have to convert to cell averages    


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

#Dam Break

wdirord = "o3"
#wdir = "../../../../data/raw/NEWdata/FDredo/grim/"
#sdir = "../../../../data/postprocessing/scFDallAE/grim/"

wdirb = "../../../../data/raw/bigsmoothtargetted/"+wdirord+"/"
dxs = []
alphas = []
TEes = []
TMes = []
TPes = []
TEas = []
TMas = []
TPas = []
RelEs = []
RelMs = []
RelPs = []
#range(6,19)
diff = "20"
g = 9.81

sdir = "../../../../data/postprocessing/CONuhHNT/"+ diff + "/"
if not os.path.exists(sdir):
        os.makedirs(sdir)
for ki in range(1,2):
    
    #Nonlinear Soliton
    dxw = str(2**ki)
    
        
    wdir = wdirb + dxw + "/" + diff + "/"
         
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
                x.append(float(row[3]))
                h.append(float(row[4]))
                u.append(float(row[6]))
                ai = float(float(row[7]))
             j = j + 1    
    n = len(x)
    o = ones(n)
    niBC = 4
    startx = x[0]
    endx = x[0]
    
    u0 = u[0]*ones(niBC)
    u1 = u[-1]*ones(niBC)   
    h0 = h[0]*ones(niBC)
    h1 = h[-1]*ones(niBC)
    
    xbeg = arange(startx - niBC*dx,startx,dx)
    xend = arange(endx + dx,endx + (niBC+1)*dx) 
    
    xbc =  concatenate([xbeg,x,xend])
    hbc =  concatenate([h0,h,h1])
    ubc =  concatenate([u0,u,u1])
    
    xbc_c = copyarraytoC(xbc)
    hbc_c = copyarraytoC(hbc)
    ubc_c = copyarraytoC(ubc)
         
    Eval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)
    Mome = uhall(x,h,u,dx) + 0.5*g*30*(1**2 - 1.8**2)
    Mass = hall(x,h,dx) 
    
    TE = Hamilana(x[0] - 0.5*dx,x[-1] + 0.5*dx,500,1.0,1.8,1.0/ai,g)
    
    TM = totalmassa(x[0] - 0.5*dx,x[-1] + 0.5*dx,500,1.0,1.8,1.0/ai)
    TP = totalmomentum() 
    
    RelE = abs(Eval - TE) / abs(TE)
    RelM = abs(Mass - TM) / abs(TM)
    RelP = abs(Mome - TP)
    
    h = array(h)
    u = array(u)
    
    alpha = 1.0/ ai
    dxs.append(dx)
    alphas.append(alpha)
    TEes.append(Eval)
    TMes.append(Mass)
    TPes.append(Mome)
    TEas.append(TE)
    TMas.append(TM)
    TPas.append(TP)
    RelEs.append(RelE)
    RelMs.append(RelM)
    RelPs.append(RelP)
    



n= len(dxs)
s = sdir + "con.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
     writefile2.writerow(['dx','alpha','TEes' ,'TMes','TPes','TEas' ,'TMas','TPas','RelEs' ,'RelMs','RelPs'])        
                       
     for j in range(n):
         writefile2.writerow([str(dxs[j]),str(alphas[j]),str(TEes[j]),str(TMes[j]), \
                              str(TPes[j]),str(TEas[j]),str(TMas[j]),str(TPas[j]),str(RelEs[j]),str(RelMs[j]),str(RelPs[j])])  

s = sdir + "conh.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",RelMs[i])
        file1.write(s)         
    
s = sdir + "conuh.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",RelPs[i])
        file1.write(s)  
        
s = sdir + "conH.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",RelEs[i])
        file1.write(s)  
