# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:09:21 2015

@author: jordan
"""
from CoupAvec import *
from scipy import *
import csv
import os
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy.linalg import norm  

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
    
def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 
    
def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a

  return a0 + a1*sech2(k*phi)
    
def exp1Sol(n,a0,a1,b0,b1,x,dx):
    v = zeros(n)
    u = zeros(n)
    for i in range(n):
        u[i] = a0 + a1*sech2(x[i])
        v[i] = b0 + b1*sech2(x[i])

    return u,v 
    
def exp1Sho(n,a0,a1,b0,b1,x,dx):
    v = zeros(n)
    u = zeros(n)
    for i in range(n):
        
        if(x[i] > -50 and x[i]< 50):
            u[i] = a1
            v[i] = b1
        else:
            u[i] = a0
            v[i] = b0

    return u,v 
        

###Chris Theta Dambreak with Hamiltonian
wdir = "../../../data/CoupAdvec/Exp1/"
if not os.path.exists(wdir):
    os.makedirs(wdir)
        
theta = 1.2    
a0 = 1.0
a1 = 0.1
b0 = 1.0
b1 = 0.1

a = 1
b = 10

dx = 0.01
Cr = 0.5
lam = 1.0 /(a*(b0 + b1) + b*(a0 + a1))
dt = Cr*lam*dx
startx = -200
endx = 100.0 + dx
startt = 0.0
endt = 5 + dt 
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)

u,v = exp1Sol(n,a0,a1,b0,b1,x,dx)  
pu,pv = exp1Sol(n,a0,a1,b0,b1,x,dx)  
  
nBC = 3
nBCs = 4
niBC = nBCs
u0 = u[0]*ones(nBCs)
u1 = u[-1]*ones(nBCs)  
v0 = v[0]*ones(nBCs)
v1 = v[-1]*ones(nBCs)

pu0 = pu[0]*ones(nBCs)
pu1 = pu[-1]*ones(nBCs)  
pv0 = pv[0]*ones(nBCs)
pv1 = pv[-1]*ones(nBCs)
    
u_c = copyarraytoC(u)
v_c = copyarraytoC(v)
FDu_c = copyarraytoC(u)
FDv_c = copyarraytoC(v)
u0_c  = copyarraytoC(u0)
u1_c  = copyarraytoC(u1)
v0_c  = copyarraytoC(v0)
v1_c  = copyarraytoC(v1)

pu_c = copyarraytoC(pu)
pv_c = copyarraytoC(pv)
pu0_c  = copyarraytoC(pu0)
pu1_c  = copyarraytoC(pu1)
pv0_c  = copyarraytoC(pv0)
pv1_c  = copyarraytoC(pv1)

xbeg = arange(startx - niBC*dx,startx,dx)
xend = arange(endx + dx,endx + (niBC+1)*dx,dx) 

xbc =  concatenate([xbeg,x,xend])  

ubc =  concatenate([u0,u,u1]) 
vbc =  concatenate([v0,v,v1]) 

pubc =  concatenate([pu0,pu,pu1]) 
pvbc =  concatenate([pv0,pv,pv1])
xbc_c = copyarraytoC(xbc)
pubc_c = copyarraytoC(pubc)
pvbc_c = copyarraytoC(pvbc)
ubc_c = copyarraytoC(ubc)
vbc_c = copyarraytoC(vbc)
FDubc_c = copyarraytoC(ubc)
FDvbc_c = copyarraytoC(vbc)

TUb = HankEnergyall(xbc_c,ubc_c,n,nBC,dx)
TVb = HankEnergyall(xbc_c,vbc_c,n,nBC,dx)
   

for i in range(1,len(t)):  
    evolvewrap(u_c,v_c,u0_c,u1_c,v0_c,v1_c,a,b,dx,dt,nBC,n,nBCs,theta)
    evolvewrapFD(FDu_c, FDv_c, pubc_c, pvbc_c, u0_c,u1_c, v0_c, v1_c, a,b,dx,dt,nBC,n,nBCs);
    print t[i]     

conc(v0_c , v_c,v1_c,niBC,n ,niBC , vbc_c)
conc(u0_c , u_c,u1_c,niBC,n ,niBC , ubc_c) 
conc(v0_c , FDv_c,v1_c,niBC,n ,niBC , FDvbc_c)
conc(u0_c , FDu_c,u1_c,niBC,n ,niBC , FDubc_c) 
TUe = HankEnergyall(xbc_c,ubc_c,n,nBC,dx)
TVe = HankEnergyall(xbc_c,vbc_c,n,nBC,dx)

TFDUe = HankEnergyall(xbc_c,FDubc_c,n,nBC,dx)
TFDVe = HankEnergyall(xbc_c,FDvbc_c,n,nBC,dx)

uend = copyarrayfromC(u_c,n)
vend = copyarrayfromC(v_c,n)      
FDuend = copyarrayfromC(FDu_c,n)
FDvend = copyarrayfromC(FDv_c,n)  

plot(x,u,label="initial u")
plot(x,v,label="initial v")
plot(x,FDuend,label="FD sol for u")
plot(x,FDvend,label="FD sol for v")
plot(x,uend,label="FVM sol for u")
plot(x,vend,label="FVM sol for v")
#xlim([-50,50])
legend()

reluerr = abs(TUe - TUb) / abs(TUb)
relverr = abs(TVe - TVb) / abs(TVb)

relFDuerr = abs(TFDUe - TUb) / abs(TUb)
relFDverr = abs(TFDVe - TVb) / abs(TVb)

print("Conservation for FVM")
print("for u : %e" %(reluerr)) 
print("for v : %e" %(relverr)) 

print("Conservation for FD")
print("for u : %e" %(relFDuerr)) 
print("for v : %e" %(relFDverr)) 


deallocPy(u_c)   
deallocPy(v_c)
deallocPy(v0_c)
deallocPy(v1_c)
deallocPy(u0_c)
deallocPy(u1_c)

