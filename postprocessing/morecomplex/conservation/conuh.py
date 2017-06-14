# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 16:53:37 2015

@author: jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import os
from numpy import ones
from Hamil import *

from numpy import tanh

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

#eta = 1
def M0p7initial(xbeg,xend):
    A = xbeg + 0.7*sqrt(68.0/21.0)*tanh(sqrt(21.0/68.0)*xbeg)
    B = xend + 0.7*sqrt(68.0/21.0)*tanh(sqrt(21.0/68.0)*xend)
    
    return B - A   
    
    
def P0p7initial(xbeg,xend):
    A = sqrt(16.677)*0.7*sqrt(68.0/21.0)*tanh(sqrt(21.0/68.0)*xbeg)
    B = sqrt(16.677)*0.7*sqrt(68.0/21.0)*tanh(sqrt(21.0/68.0)*xend)
    
    return B - A  
    
def M1initial(xbeg,xend):
    A = xbeg + 2*sqrt(2.0/3.0)*tanh(sqrt(3.0/8.0)*xbeg)
    B = xend + 2*sqrt(2.0/3.0)*tanh(sqrt(3.0/8.0)*xend)
    
    return B - A   
    
    
def P1initial(xbeg,xend):
    A = 7.23326*tanh(0.612372*xbeg)
    B = 7.23326*tanh(0.612372*xend)
    
    return B - A  
    
    
def SolE(xbeg,xend):
    AB = 21.0068*tanh(0.555719*xbeg) - 19.2569*arctanh(0.641689*tanh(0.555719*xbeg))
    AE = 21.0068*tanh(0.555719*xend) - 19.2569*arctanh(0.641689*tanh(0.555719*xend))
    
    BB = 9.81*(xbeg) + tanh(0.555719*xbeg)*(2.88329*sech(0.555719*xbeg)**2 + 30.4805)
    BE =9.81*(xend) + tanh(0.555719*xend)*(2.88329*sech(0.555719*xend)**2 + 30.4805)
    
    CB = 307.641*(tanh(0.555719*xbeg)*(0.049539 - 0.00937224*sech(0.555719*xbeg)**2) -0.0625954*arctanh(0.641689*(tanh(0.555719*xbeg))))
    CE = 307.641*(tanh(0.555719*xend)*(0.049539 - 0.00937224*sech(0.555719*xend)**2) -0.0625954*arctanh(0.641689*(tanh(0.555719*xend))))

    
    A = AE - AB  
    B = BE - BB    
    C = CE - CB


    #1527.68293
    return 0.5*(A + B + C)

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
    


def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 

def sech(x):
  a = 2./(exp(x) + exp(-x))
  return a

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
wdirord = "FDcent"
#wdir = "../../../../data/raw/NEWdata/FDredo/grim/"
#sdir = "../../../../data/postprocessing/scFDallAE/grim/"

sdir = "../../../../data/postprocessing/FDREREDON/"+wdirord+"/"
if not os.path.exists(sdir):
        os.makedirs(sdir)
Mns = []
Pns = []
Mis = []
Pis = []
dxs = []
Ens = []
Eis = []
#range(6,19)
g = 9.81
for ki in range(6,20):
    
    #Nonlinear Soliton
    dxw = str(ki)
    
        
    wdir = "../../../../data/raw/FDreredo/"  +wdirord +"/" + dxw + "/"
         
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
    
    #hi,ui = solitoninit(n,1,1,9.81,x,0,dx)

    En = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)
    Pn = uhall(xbc_c,hbc_c,ubc_c,n + 2*niBC,niBC,dx)
    Mn = hall(xbc_c,hbc_c,n + 2*niBC,niBC,dx) 
    
    xbeg = -50 - 0.5*dx
    xend = 250 + 0.5*dx

    Pi = P0p7initial(xbeg,xend)
    Mi = M0p7initial(xbeg,xend)
    Ei = SolE(xbeg,xend)
    
    #Mi = Hamiltonianall(x,hi,o,dx)
    #Pi = Hamiltonianall(x,hi,ui,dx)
    
    Pns.append(Pn)
    Pis.append(Pi)
    Mns.append(Mn)
    Mis.append(Mi)
    dxs.append(dx)
    Ens.append(En)
    Eis.append(Ei)

Ens = array(Ens)
Eis = array(Eis)    
Mns = array(Mns)
Pns = array(Pns)
Mis = array(Mis)
Pis = array(Pis)

relerrP = abs(Pis - Pns)/ abs(Pis)
relerrM = abs(Mis - Mns)/ abs(Mis)
relerrE = abs(Eis - Ens)/ abs(Eis)

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
        
s = sdir + "conH.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",relerrE[i])
        file1.write(s) 
"""

wdirord = "o3"
#wdir = "../../../../data/raw/NEWdata/FDredo/grim/"
#sdir = "../../../../data/postprocessing/scFDallAE/grim/"

wdirb = "../../../../data/raw/bigsmoothtargetted/"+wdirord+"/"
Mns = []
Pns = []
Mis = []
Pis = []
dxs = []
Ens = []
Eis = []
alphas = []
#range(6,19)
diffi= "12"
g = 9.81

sdir = "../../../../data/postprocessing/CONuhHNAT/"+ diffi + "/"

if not os.path.exists(sdir):
        os.makedirs(sdir)
for ki in range(1,2):
    
    #Nonlinear Soliton
    dxw = str(2**ki)
    
        
    wdir = wdirb + "/" + dxw + "/" + diffi + "/"
         
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
      
    #plot(x,h)
    #plot(x,he)
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
    
    #hi,ui = solitoninit(n,1,1,9.81,x,0,dx)

    En = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)
    Pn = uhall(xbc_c,hbc_c,ubc_c,n + 2*niBC,niBC,dx) + 0.5*g*30*(h[-1]**2 - h[0]**2)
    Mn = hall(xbc_c,hbc_c,n + 2*niBC,niBC,dx)
    
    xbeg = 0 - 0.5*dx
    xend = 1000 + 0.5*dx

    
    Ei = Hamilana(x[0] - 0.5*dx,x[-1] + 0.5*dx,500,1.0,1.8,1.0/ai,g)
    
    Mi = totalmassa(x[0] - 0.5*dx,x[-1] + 0.5*dx,500,1.0,1.8,1.0/ai)
    Pi = totalmomentum() 
    
    #Mi = Hamiltonianall(x,hi,o,dx)
    #Pi = Hamiltonianall(x,hi,ui,dx)
    
    alpha = 1.0/ ai
    alphas.append(alpha)
    
    Pns.append(Pn)
    Pis.append(Pi)
    Mns.append(Mn)
    Mis.append(Mi)
    dxs.append(dx)
    Ens.append(En)
    Eis.append(Ei)

Ens = array(Ens)
Eis = array(Eis)    
Mns = array(Mns)
Pns = array(Pns)
Mis = array(Mis)
Pis = array(Pis)

relerrP = abs(Pis - Pns)
relerrM = abs(Mis - Mns)/ abs(Mis)
relerrE = abs(Eis - Ens)/ abs(Eis)


n= len(dxs)
s = sdir + "con.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
     writefile2.writerow(['dx','alphas','Pis','Pns' ,'relerrP','Mis','Mns','relerrM'])        
                       
     for j in range(n):
         writefile2.writerow([str(dxs[j]),str(alphas[j]),str(Pis[j]),str(Pns[j]),str(relerrP[j]),str(Mis[j]),str(Mns[j]),str(relerrM[j])])  

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
        
s = sdir + "conH.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",relerrE[i])
        file1.write(s) 
