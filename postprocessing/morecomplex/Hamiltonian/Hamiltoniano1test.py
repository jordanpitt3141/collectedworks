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
    
#Tested from CFD forum post
def pentadiagsolve(e,a,d,c,f,B):
    n = len(d)
    X = zeros(n)
    
    for i in range(1,n-1):
        xmult = float(a[i-1]) / d[i-1]
        
        d[i] = d[i] - xmult*c[i-1]
        c[i] = c[i] - xmult*f[i-1]
        B[i] = B[i] - xmult*B[i-1]
        
        xmult = float(e[i-1]) /d[i-1]
        a[i] = a[i] - xmult*c[i-1]
        d[i+1] = d[i+1] - xmult*f[i-1]
        B[i+1] = B[i+1] - xmult*B[i-1]
        
    xmult = float(a[n-2]) / d[n-2]
    d[n-1] = d[n-1] - xmult*c[n-2]
    X[n-1] = (B[n-1] - xmult*B[n-2]) / float(d[n-1])
    X[n-2] = (B[n-2] - c[n-2]*X[n-1]) / float(d[n-2])
    
    for i in range(n-3,-1,-1):
        X[i] = (B[i] - f[i]*X[i+2] - c[i]*X[i+1])/float(d[i])
        
    return X    
    
def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 

def midpointtocellaverages(mq,dx):
    #no BC required, assumes that the averages and midpoints at the boundaries are the same
    idx = 1.0/dx
    i24 = 1.0 / 24.0
    n = len(mq)
    
    a = zeros(n-1)
    b = zeros(n)
    c = zeros(n-1)
    for i in range(1,n-1):
        ai = -i24
        bi = 26*i24
        ci = -i24

        a[i-1] = ai
        b[i] = bi
        c[i] = ci
    
    #i = 0
    i = 0
    ai =0.0 #-i24
    bi =1.0 #26*i24
    ci =0.0 #-i24

    b[i] = bi
    c[i] = ci
    
    #mq[i] = mq[i] - ai*qbeg[0]
    
    #i = 0
    i = n-1
    ai =0.0# -i24
    bi =1.0# 26*i24
    ci =0.0# -i24

    a[i-1] = ai
    b[i] = bi
    
    #mq[i] = mq[i] - ci*qend[0]
    
    q = TDMApy(a,b,c,mq)
    
    return q
    
def cellaveragestomidpoints(q,dx):
    #no BC required, assumes that the averages and midpoints at the boundaries are the same
    i24 = 1.0 / 24.0
    n = len(q)
    mq = zeros(n)
    for i in range(1,n-1):
        #iterate  over the cell midpoints, there are 2 edge values for each (except the first and last cell)
        
        #variables
        #ai = (q[i+1] - 2*q[i] + q[i-1])*0.5*idx*idx
        #bi = (q[i+1] - q[i-1])*0.5*idx
        ci = i24*(-q[i+1] + 26*q[i]  -q[i-1])
        mq[i] = ci
    
    #i = 0
    i = 0
    ci = q[i] #i24*(-q[i+1] + 26*q[i] - qbeg[0])
    mq[i] = ci
    
    #i = n-1
    i = n-1
    ci = q[i]#i24*(-qend[0] + 26*q[i] - q[i-1])
    mq[i] = ci 
    
    return mq

def solveufromGh(G,h,hbeg,hend,ubeg,uend,dx):
    #takes midpoint values of G,h and gives midpoint values of u
    idx = 1.0 / dx
    i12 = 1.0 / 12.0
    i3 = 1.0 / 3.0
    n = len(G)
    
    a = zeros(n-2)
    b = zeros(n-1)
    c = zeros(n)
    d = zeros(n-1)
    e = zeros(n-2)
    
    for i in range(2,n-2):
        th = h[i]
        thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*h[i-1] + h[i-2] )
        
        ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
        bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
        ci = th + (30*i12*idx*idx)*(i3*th*th*th)
        di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
        ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2
        
        
        
        a[i-2] = ai
        b[i-1] =  bi
        c[i] = ci
        d[i] = di
        e[i] = ei
        
    #boundary    
    #i=0
    i=0
    th = h[i]
    thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*hbeg[-1] + hbeg[-2] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2

 
    c[i] = ci
    d[i] = di
    e[i] = ei
    
    G[i] = G[i] - ubeg[-1]*bi - ubeg[-2]*ai
    
    #i=1
    i=1
    th = h[i]
    thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*h[i-1] + hbeg[-1] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2

 
    c[i] = ci
    d[i] = di
    e[i] = ei
    b[i-1] = bi 
    
    G[i] = G[i] - ubeg[-1]*ai
    
    #boundary    
    #i=n-2
    i=n-2
    th = h[i]
    thx = i12*idx*(-hend[0] + 8*h[i+1] - 8*h[i-1] + h[i-2] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2

    a[i-2] = ai
    b[i-1] = bi 
    c[i] = ci
    d[i] = di
    
    G[i] = G[i]- uend[0]*ei
    
    #i=n-1
    i=n-1
    th = h[i]
    thx = i12*idx*(-hend[1] + 8*hend[0] - 8*h[i-1] + h[i-2] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2

    a[i-2] = ai
    b[i-1] = bi 
    c[i] = ci
    
    G[i] = G[i] -uend[0]*di - uend[1]*ei
    
    u = pentadiagsolve(a,b,c,d,e,G)
    return u
    
def solveGfromuh(u,h,hbeg,hend,ubeg,uend,dx):
    #takes midpoint values of u,h and gives midpoint values of G    
    
    idx = 1.0 / dx
    i12 = 1.0 / 12.0
    i3 = 1.0 / 3.0
    n = len(u)
    
    G = zeros(n)
    
    for i in range(2,n-2):
        th = h[i]
        thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*h[i-1] + h[i-2] )
        
        ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
        bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
        ci = th + (30*i12*idx*idx)*(i3*th*th*th)
        di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
        ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2
        
        G[i] = ai*u[i-2] + bi*u[i-1] + ci*u[i] + di*u[i+1] + ei*u[i+2]
        
        
    #boundary    
    #i=0
    i=0
    th = h[i]
    thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*hbeg[-1] + hbeg[-2] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2
    
    G[i] = ai*ubeg[-2] + bi*ubeg[-1] + ci*u[i] + di*u[i+1] + ei*u[i+2]

    
    #i=1
    i=1
    th = h[i]
    thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*h[i-1] + hbeg[-1] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2

    G[i] = ai*ubeg[-1] + bi*u[i-1] + ci*u[i] + di*u[i+1] + ei*u[i+2]
    
    #boundary    
    #i=n-2
    i=n-2
    th = h[i]
    thx = i12*idx*(-hend[0] + 8*h[i+1] - 8*h[i-1] + h[i-2] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2

    G[i] = ai*u[i-2] + bi*u[i-1] + ci*u[i] + di*u[i+1] + ei*uend[0]
    
    #i=n-1
    i=n-1
    th = h[i]
    thx = i12*idx*(-hend[1] + 8*hend[0] - 8*h[i-1] + h[i-2] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2

    G[i] = ai*u[i-2] + bi*u[i-1] + ci*u[i] + di*uend[0] + ei*uend[1]

    return G

def quadraticrecon(q,dx,i):
    i24 = 1.0/24.0
    i12 = 1.0/12.0
    idx = 1.0/dx
    ai = (i24*idx*idx*idx*idx)*(q[i+2] - 4*q[i+1] + 6*q[i] - 4*q[i-1] + q[i-2])
    bi = (i12*idx*idx*idx)*(q[i+2] - 2*q[i+1] + 2*q[i-1] - q[i-2])
    ci = (i24*idx*idx)*(-q[i+2] + 16*q[i+1] - 30*q[i] + 16*q[i-1] - q[i-2])
    di = (i12*idx)*(-q[i+2] + 8*q[i+1] - 8*q[i-1] + q[i-2])
    ei = q[i]   
    
    return ai,bi,ci,di,ei
    
def quadraticcell(q,i,dx,x,ex):
    ai,bi,ci,di,ei = quadraticrecon(q,dx,i)
    return ai*(ex - x[i])**4 + bi*(ex - x[i])**3 + ci*(ex - x[i])**2 + di*(ex - x[i]) + ei
    
def quadraticderivcell(q,i,dx,x,ex):
    ai,bi,ci,di,ei = quadraticrecon(q,dx,i)
    return 4*ai*(ex - x[i])**3 + 3*bi*(ex - x[i])**2 + 2*ci*(ex - x[i]) + di

def hamiltonianincell(u,h,epsilon,sigma,i,x,dx):
    sqrt35 = sqrt(3.0/5.0)
    iepsilon = 1.0 / epsilon
    i9 = 1.0 / 9.0
    
    firstpth = epsilon*quadraticcell(h,i,dx,x,-0.5*dx*sqrt35 + x[i])
    firstptu = epsilon*quadraticcell(u,i,dx,x,-0.5*dx*sqrt35 + x[i])
    firstptdu = epsilon*quadraticderivcell(u,i,dx,x,-0.5*dx*sqrt35 + x[i])
    
    secondpth = epsilon*quadraticcell(h,i,dx,x,x[i])
    secondptu = epsilon*quadraticcell(u,i,dx,x,x[i])
    secondptdu = epsilon*quadraticderivcell(u,i,dx,x,x[i])
    
    thirdpth = epsilon*quadraticcell(h,i,dx,x,0.5*dx*sqrt35 + x[i])
    thirdptu = epsilon*quadraticcell(u,i,dx,x,0.5*dx*sqrt35 + x[i]) 
    thirdptdu = epsilon*quadraticderivcell(u,i,dx,x,0.5*dx*sqrt35 + x[i])
    
    firstptham = epsilon*firstpth*firstptu*firstptu + 0.5*epsilon*sigma*sigma*firstpth*firstpth*firstpth*firstptdu*firstptdu \
    + iepsilon*firstpth*firstpth
    
    secondptham = epsilon*secondpth*secondptu*secondptu + 0.5*epsilon*sigma*sigma*secondpth*secondpth*secondpth*secondptdu*secondptdu \
    + iepsilon*secondpth*secondpth
    
    thirdptham = epsilon*thirdpth*thirdptu*thirdptu + 0.5*epsilon*sigma*sigma*thirdpth*thirdpth*thirdpth*thirdptdu*thirdptdu \
    + iepsilon*thirdpth*thirdpth
    
    return 0.5*dx*( 5*i9*firstptham + 8*i9*secondptham + 5*i9*thirdptham )

def hamiltonianallcells(ubc,hbc,epsilon,sigma,xbc,nBC,dx):
    n = len(xbc)
    for i in range(nBC,n-nBC):
        hamiltonianincell(u,h,epsilon,sigma,i,x,dx)
        
    
    return 0.5*dx*( 5*i9*firstptham + 8*i9*secondptham + 5*i9*thirdptham )    
    
    
def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a

def soliton (x,t,g,a0,a1):
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)
  
def solitoninit(n,a0,a1,g,x,t0,dx):
    #t0 = 0
    h = zeros(n)
    u = zeros(n)
    G = zeros(n)
    c = sqrt(g*(a0 + a1))
    k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
    for i in range(n):
        h[i] = soliton(x[i],t0,g,a0,a1)
        u[i] =  c* ((h[i] - a0) / h[i])
        endG = 4*a1*a1*k*k*tanh(k*x[i])*tanh(k*x[i])*sech2(k*x[i])*sech2(k*x[i]) - 2*h[i]*a1*k*k*sech2(k*x[i])*(1 - 3*tanh(k*x[i])*tanh(k*x[i]))
        G[i] = u[i]*h[i] - (c*a0/3.0)*endG
    
    return h,u,G
 
## Setup   
dx = 0.1
hdx = 0.01
dx2hdx = 0.1
hdx2dx = 10
a0 = 1.0
a1 = 1.0
Cr = 0.5
g = 9.81
l = Cr / (sqrt(g*(a0 + a1)))
dt = l*dx
startx = -100.0
endx = 500.0
startt = 0.0
endt = 100 + dt
    
szoomx = startx
ezoomx = endx
    
#number of boundary conditions (one side)
nfcBC = 4 #for flux calculation
nGsBC = 4 #for solving G from u,h
niBC = nGsBC + nfcBC #total
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
xbc,t = makevar(startx - nGsBC*dx,endx+ nGsBC*dx,dx,startt,endt,dt)
hx,ht = makevar(startx,endx,hdx,startt,endt,dt)
n = len(x)

t0 = 0   
hm,um,Gm = solitoninit(n,a0,a1,g,x,t0,dx)

umbeg = um[0]*ones(niBC)
umend = um[-1]*ones(niBC)
hmbeg = hm[0]*ones(niBC)
hmend = hm[-1]*ones(niBC)    
Gmbeg = Gm[0]*ones(niBC)
Gmend = Gm[-1]*ones(niBC) 
        
#calculate G midpoint
cnBC = niBC - nGsBC
        
umbc = concatenate([umbeg[-cnBC:],um,umend[0:cnBC]]) 
hmbc = concatenate([hmbeg[-cnBC:],hm,hmend[0:cnBC]]) 
Gmbc = concatenate([Gmbeg[-cnBC:],Gm,Gmend[0:cnBC]])  

#using third order FD       
Gmbcn = solveGfromuh(umbc,hmbc,hmbeg[0:-cnBC],hmend[-cnBC:],umbeg[0:-cnBC],umend[-cnBC:],dx)  
Gmn = Gmbc[cnBC:-cnBC]

#calculate averages
Gabc = midpointtocellaverages(Gmbc,dx)
habc = midpointtocellaverages(hmbc,dx)
uabc = midpointtocellaverages(umbc,dx)
        
#so we can just go from here with Ga ang ha?
Gabeg = Gabc[0:cnBC]
Ga = Gabc[cnBC:-cnBC]
Gaend = Gabc[-cnBC:]
habeg = habc[0:cnBC]
ha = habc[cnBC:-cnBC]
haend = habc[-cnBC:]
uabeg = uabc[0:cnBC]
ua = uabc[cnBC:-cnBC]
uaend = uabc[-cnBC:]



### Given cell averages get midpoints and solve for u
hp = cellaveragestomidpoints(ha,dx)
Gp = cellaveragestomidpoints(Ga,dx)
up = solveufromGh(Gp,hp,hmbeg,hmend,umbeg,umend,dx)

upmbc = concatenate([umbeg[-cnBC:],up,umend[0:cnBC]]) 
hpmbc = concatenate([hmbeg[-cnBC:],hp,hmend[0:cnBC]]) 
Gpmbc = concatenate([Gmbeg[-cnBC:],Gp,Gmend[0:cnBC]])  


##### TEST Reconstruction

"""
#test specific 
xbeg = (100.0)/dx
ht1  = quadraticcell(hp,xbeg,dx,x,x[xbeg])
ht2  = quadraticcell(hp,xbeg,dx,x,x[xbeg+1])
ht3  = quadraticcell(hp,xbeg,dx,x,x[xbeg+2])
ht4  = quadraticcell(hp,xbeg,dx,x,x[xbeg-1])
ht5  = quadraticcell(hp,xbeg,dx,x,x[xbeg-2])

print(hp[xbeg],ht1,hp[xbeg]-ht1)
print(hp[xbeg+1],ht2,hp[xbeg+1]-ht2 )
print(hp[xbeg+2],ht3,hp[xbeg+2]-ht3)
print(hp[xbeg-1],ht4,hp[xbeg-1]-ht4)
print(hp[xbeg-2],ht5,hp[xbeg-2]-ht5)
"""


#test all
"""
hn = len(hx)
nbc = len(upmbc)
hpr = zeros(hn)
Gpr = zeros(hn)
upr = zeros(hn)
for i in range(cnBC,nbc - cnBC ):
    k = i - cnBC
    for j in range(hdx2dx*k,hdx2dx*k +hdx2dx):
        ex = hx[j]
        #print(x[k],hx[j])
        print(i,k,j)
        hpr[j] = quadraticcell(hpmbc,i,dx,xbc,ex)
        Gpr[j] = quadraticcell(Gpmbc,i,dx,xbc,ex)
        upr[j] = quadraticcell(upmbc,i,dx,xbc,ex)
        
"""
    

#####get hamiltonian numerically

#ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)