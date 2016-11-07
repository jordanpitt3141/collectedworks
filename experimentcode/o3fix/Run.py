from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv
import os
from Serre3 import *
from numpy.linalg import norm
from time import time

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
    

def interpquarticvalPy(aj,bj,cj,dj,ej,xj,x):
    
    return aj*(x -xj)*(x -xj)*(x -xj)*(x -xj) + bj*(x -xj)*(x -xj)*(x -xj) \
    + cj*(x -xj)*(x -xj) + dj*(x -xj)+ ej
    
def interpquarticgradPy(aj,bj,cj,dj,ej,xj,x):
    
    return 4*aj*(x -xj)*(x -xj)*(x -xj) + 3*bj*(x -xj)*(x -xj) \
    + 2*cj*(x -xj) + dj
    
def interpquartcoeffPy(q,j,dx):
    i24 = 1.0 / 24.0
    i12 = 1.0 / 12.0
    idx = 1.0/dx
    aj = i24*idx*idx*idx*idx*(q[j+2] - 4*q[j+1] + 6*q[j] - 4*q[j-1] + q[j-2])
    bj = i12*idx*idx*idx*(q[j+2] - 2*q[j+1] + 2*q[j-1] - q[j-2])
    cj = i24*idx*idx*(-q[j+2] + 16*q[j+1] - 30*q[j] + 16*q[j-1] - q[j-2])
    dj = i12*idx*(-q[j+2] + 8*q[j+1] - 8*q[j-1] + q[j-2])
    ej = q[j]
    
    return aj,bj,cj,dj,ej
    
def HankEnergyacrosscellPy(x,h,u,g,j,dx):
    #so we have h,u at midpoints
    #epsilon and sigma are everywhere
    i3 = 1.0/3.0

    #jth cell
    uaj,ubj,ucj,udj,uej = interpquartcoeffPy(u,j,dx)
    haj,hbj,hcj,hdj,hej = interpquartcoeffPy(h,j,dx)
    
    #first gauss point
    fgp = 0.5*dx*sqrt(3.0/5.0) + x[j]
    fgph = interpquarticvalPy(haj,hbj,hcj,hdj,hej,x[j],fgp)
    fgpu = interpquarticvalPy(uaj,ubj,ucj,udj,uej,x[j],fgp)
    fgpux = interpquarticgradPy(uaj,ubj,ucj,udj,uej,x[j],fgp)
    
    fgpe = fgph*fgpu*fgpu + g*fgph*fgph + i3*(fgph*fgph*fgph)*fgpux*fgpux
        
    #second gauss point
    sgp = x[j]
    sgph = interpquarticvalPy(haj,hbj,hcj,hdj,hej,x[j],sgp)
    sgpu = interpquarticvalPy(uaj,ubj,ucj,udj,uej,x[j],sgp)
    sgpux = interpquarticgradPy(uaj,ubj,ucj,udj,uej,x[j],sgp)
    
    sgpe = sgph*sgpu*sgpu + g*sgph*sgph + i3*(sgph*sgph*sgph)*sgpux*sgpux

    #third gauss point
    tgp = -0.5*dx*sqrt(3.0/5.0) + x[j]
    tgph = interpquarticvalPy(haj,hbj,hcj,hdj,hej,x[j],tgp)
    tgpu = interpquarticvalPy(uaj,ubj,ucj,udj,uej,x[j],tgp)
    tgpux = interpquarticgradPy(uaj,ubj,ucj,udj,uej,x[j],tgp)
    
    tgpe = tgph*tgpu*tgpu + g*tgph*tgph + i3*(tgph*tgph*tgph)*tgpux*tgpux
        
    Hamilcell = 0.5*dx*( (5.0/9.0)*fgpe + (8.0/9.0)*sgpe + (5.0/9.0)*tgpe)
    
    return Hamilcell
    
def HankEnergyallPy(x,h,u,g,nBC,dx):

    n = len(x)
    sum1 = 0.0
    #Hamilbycell = []
    for i in range(nBC,n - nBC):
       sum1 = sum1 + HankEnergyacrosscellPy(x,h,u,g,i,dx)
       #Hamilbycell.append(HankEnergyacrosscell(x,h,u,g,i,dx))
    return 0.5*sum1#, Hamilbycell

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
    

def dambreak(x,xc,hf,hl):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    
    for i in range(n):
       
        if(x[i] >= xc):
            h[i] = hl
        else:
            h[i]= hf
        
    return u,h
    

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

def soliton2interactinit(n,a0,a11,solbeg1,solend1,direction1,a12,solbeg2,solend2,direction2,g,x,t0,dx):
    h = zeros(n)
    u = zeros(n)
    c1 = sqrt(g*(a0 + a11))
    c2 = sqrt(g*(a0 + a12))
    for i in range(n):
        if (x[i] > solbeg1 and x[i] < solend1):
            h[i] = soliton(abs(x[i] - 0.5*(solbeg1 + solend1)),t0,g,a0,a11)
            u[i] = direction1*c1*( (h[i] - a0) / h[i] )
        elif (x[i] > solbeg2 and x[i] < solend2):
            h[i] = soliton(abs(x[i] - 0.5*(solbeg2 + solend2)),t0,g,a0,a12)
            u[i] =  direction2*c2* ((h[i] - a0) / h[i])
        else:
            h[i] = a0
            u[i] = 0.0
    
    return h,u        

def experiment1(x,b,h0,h1,dx):
    n = len(x)
    u = zeros(n)
    h = ones(n)*h1
    for i in range(n):
        if (x[i] <0 and x[i] > -2*b):
            h[i] = h0

    return h,u
    
def experiment2(x,h0,u0,h1,u1,dx):
    n = len(x)
    u = ones(n)*u0
    h = ones(n)*h0
    for i in range(n):
        if (i > n/2):
            h[i] = h1
            u[i] = u1

    return h,u
    
def dambreaksmooth(x,x0,base,eta0,diffuse,dx):
    from numpy import tanh
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    dh = zeros(n)
    du = zeros(n)
    ddu = zeros(n)
    G = zeros(n)
    
    for i in range(n):
        h[i] = base + 0.5*eta0*(1 + tanh(diffuse*(x0 - x[i])))
        u[i] = 0.0

    return h,u
    
def dambreaksmoothG(x,x0,base,eta0,diffuse,dx):
    from numpy import tanh
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    dh = zeros(n)
    du = zeros(n)
    ddu = zeros(n)
    G = zeros(n)
    
    for i in range(n):
        h[i] = base + 0.5*eta0*(1 + tanh(diffuse*(x0 - x[i])))
        u[i] = 0.0
        dh[i] = -0.5*diffuse*etah*sech2(diffuse*(x0 - x[i]))
        du[i] = 0.0
        ddu[i] = 0.0
        G[i] = u[i]*h[i] - h[i]*h[i]*dh[i]*du[i] - (1.0/3.0)*h[i]*h[i]*h[i]*ddu[i]

    return h,u,G

def DSWsmooth(x,x0,base,eta0,diffuse,dx):
    from numpy import tanh
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    
    for i in range(n):
        h[i] = base + 0.5*eta0*(1 + tanh(diffuse*(x0 - abs(x[i]))))
        u[i] = 2*(sqrt(g*h[i]) -1)

    return h,u    
    
def DBSsmooth(x,x0,baseh,etah,baseu,etau,diffuse,dx):
    from numpy import tanh
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    dh = zeros(n)
    du = zeros(n)
    ddu = zeros(n)
    G = zeros(n)
    
    for i in range(n):
        h[i] = baseh + 0.5*etah*(1 + tanh(diffuse*(x0 - x[i])))
        u[i] = baseu + 0.5*etau*(1 + tanh(diffuse*(x0 - x[i])))
        dh[i] = -0.5*diffuse*etah*sech2(diffuse*(x0 - x[i]))
        du[i] = -0.5*diffuse*etau*sech2(diffuse*(x0 - x[i]))
        ddu[i] = -diffuse*diffuse*etau*tanh(diffuse*(x0 - x[i]))*sech2(diffuse*(x0 - x[i]))
        G[i] = u[i]*h[i] - h[i]*h[i]*dh[i]*du[i] - (1.0/3.0)*h[i]*h[i]*h[i]*ddu[i]
        
        
    print(x[102],ddu[102],diffuse*diffuse*tanh(diffuse*(x0 - x[102]))*sech2(diffuse*(x0 - x[102])))

    return h,u,G 
    
def DBSREVsmooth(x,x0,baseh,etah,baseu,etau,diffuse,dx):
    from numpy import tanh
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    dh = zeros(n)
    du = zeros(n)
    ddu = zeros(n)
    G = zeros(n)
    
    for i in range(n):
        h[i] = baseh + 0.5*etah*(1 + tanh(diffuse*(x0 - x[i])))
        u[i] = baseu + 0.5*etau*(1 + tanh(diffuse*(x0 - x[i])))
        dh[i] = -0.5*diffuse*etah*sech2(diffuse*(x0 - x[i]))
        du[i] = -0.5*diffuse*etau*sech2(diffuse*(x0 - x[i]))
        ddu[i] = -diffuse*diffuse*etau*tanh(diffuse*(x0 - x[i]))*sech2(diffuse*(x0 - x[i]))
        G[i] = u[i]*h[i] - h[i]*h[i]*dh[i]*du[i] - (1.0/3.0)*h[i]*h[i]*h[i]*ddu[i]
        
        
    print(x[102],ddu[102],diffuse*diffuse*tanh(diffuse*(x0 - x[102]))*sech2(diffuse*(x0 - x[102])))

    return h,u,G 
"""
#DB Smooth ASPECT
#diffuses = [0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1.0,2.5,5.0,7.5,10.0,25.0,50.0,75.0,100.0,250.0,500.0,750.0,1000.0]
diffuses = [10]
#diffuses = [0.025]
dxlist = [10]

#diffuses = [2.5]
#dxlist = [5]
wdirb = "../../../data/raw/DBASPECTRAT/o3/"
hf = 8.0
#for ll in range(3,16):
for ll in dxlist:
    for k in range(len(diffuses)):
        wdir = wdirb + str(ll) + "/" + str(diffuses[k]) + "/" +str(hf)+ "/"
        if not os.path.exists(wdir):
            os.makedirs(wdir) 
        dx = 10.0 / (2**ll)
        l = 0.01
        dt = l*dx
        startx = -600.0
        endx = 1400.0 + dx
        startt = 0.0
        endt = 100 + (dt*0.9)  
                
        szoomx = startx
        ezoomx = endx
                
        #number of boundary conditions (one side)
        nfcBC = 4 #for flux calculation
        nGsBC = 2 #for solving G from u,h
        niBC = nGsBC + nfcBC #total
                
        g = 9.81
        
        gap = int(1.0/dt)
                
        x,t = makevar(startx,endx,dx,startt,endt,dt)
        n = len(x)
        hl = 1.0
        baseh = hl
        etah = hf - hl
        x0 = 500
        diffuse = diffuses[k]   
        hm,um,Gm = dambreaksmoothG(x,x0,baseh,etah,diffuse,dx)
        #plot(x,hm)
        #plot(x,um)
        #show()
                        
        umbegi = um[0]*ones(niBC)
        umendi = um[-1]*ones(niBC)
        hmbegi = hm[0]*ones(niBC)
        hmendi = hm[-1]*ones(niBC)   
        Gmbegi = Gm[0]*ones(niBC)
        Gmendi = Gm[-1]*ones(niBC) 
                
        umbeg = umbegi
        umend = umendi
        hmbeg = hmbegi
        hmend = hmendi
        Gmbeg = Gmbegi
        Gmend = Gmendi
                
        #calculate G midpoint
        cnBC = niBC - nGsBC
                
        umbc = concatenate([umbeg[-cnBC:],um,umend[0:cnBC]]) 
        hmbc = concatenate([hmbeg[-cnBC:],hm,hmend[0:cnBC]])       
        Gmbc = concatenate([Gmbeg[-cnBC:],Gm,Gmend[0:cnBC]])  
                
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
            
        Ga_c = copyarraytoC(Ga)
        Gabeg_c = copyarraytoC(Gabeg)
        Gaend_c = copyarraytoC(Gaend)
        ha_c = copyarraytoC(ha)
          
        habeg_c = copyarraytoC(habeg)
        haend_c = copyarraytoC(haend)
          
        uabeg_c = copyarraytoC(uabeg)
        uaend_c = copyarraytoC(uaend)
           
        hmbeg_c = copyarraytoC(hmbeg)
        hmend_c = copyarraytoC(hmend)
            
        umbeg_c = copyarraytoC(umbeg)
        umend_c = copyarraytoC(umend)
           
        u_c= mallocPy(n)
        G_c = mallocPy(n)
        h_c = mallocPy(n)
        
        xbeg = arange(startx - niBC*dx,startx,dx)
        xend = arange(endx + dx,endx + (niBC+1)*dx) 
        
        xbc =  concatenate([xbeg,x,xend])  
        
        xbc_c = copyarraytoC(xbc)
        hbc_c = mallocPy(n + 2*niBC)
        ubc_c = mallocPy(n + 2*niBC)
        Evals = []
            
        for i in range(1,len(t)):  
            
            if(i ==1 or i %gap == 0):
                ca2midpt(ha_c,dx,n,h_c)
                ca2midpt(Ga_c,dx,n,G_c)
                ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
                
                conc(hmbeg_c , h_c,hmend_c,niBC,n ,niBC , hbc_c)
                conc(umbeg_c , u_c,umend_c,niBC,n ,niBC , ubc_c)        
                Eval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)
                
                Evals.append(Eval)
                u = copyarrayfromC(u_c,n)
                G = copyarrayfromC(G_c,n)
                h = copyarrayfromC(h_c,n)
                s = wdir + "out" + str(i) +  ".txt"
                with open(s,'a') as file2:
                     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                    
                     writefile2.writerow(['dx' ,'dt','time','Eval','cell midpoint', 'height(m)', 'G' , 'u(m/s)',"diffuse"])        
                                   
                     for j in range(n):
                         writefile2.writerow([str(dx),str(dt),str(t[i]),str(Eval), str(x[j]), str(h[j]) , str(G[j]) , str(u[j]), str(diffuse)])   
            evolvewrap(Ga_c,ha_c,Gabeg_c,Gaend_c,habeg_c,haend_c,hmbeg_c,hmend_c,uabeg_c,uaend_c,umbeg_c,umend_c,nfcBC,nGsBC,g,dx,dt,n,cnBC,niBC)
            print (t[i])
            
                
        ca2midpt(ha_c,dx,n,h_c)
        ca2midpt(Ga_c,dx,n,G_c)
        ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
        
        conc(hmbeg_c , h_c,hmend_c,niBC,n ,niBC , hbc_c)
        conc(umbeg_c , u_c,umend_c,niBC,n ,niBC , ubc_c)        
        Eval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)
        
        Evals.append(Eval)
        u = copyarrayfromC(u_c,n)
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        s = wdir + "outlast.txt"
        with open(s,'a') as file2:
             writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
             writefile2.writerow(['dx' ,'dt','time','Eval','cell midpoint', 'height(m)', 'G' , 'u(m/s)',"diffuse"])        
                           
             for j in range(n):
                 writefile2.writerow([str(dx),str(dt),str(t[i]),str(Eval), str(x[j]), str(h[j]) , str(G[j]) , str(u[j]), str(diffuse)])    
        deallocPy(u_c)   
        deallocPy(h_c)
        deallocPy(G_c)
        deallocPy(ha_c)
        deallocPy(Ga_c)
        deallocPy(habeg_c)
        deallocPy(Gabeg_c)
        deallocPy(haend_c)
        deallocPy(Gaend_c)
        deallocPy(uabeg_c)
        deallocPy(uaend_c)
        deallocPy(umbeg_c)
        deallocPy(umend_c)
        deallocPy(hmbeg_c)
        deallocPy(hmend_c)

"""


"""
#DSW Reverse
#diffuses = [0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1.0,2.5,5.0,7.5,10.0,25.0,50.0,75.0,100.0,250.0,500.0,750.0,1000.0]
diffuses = [10]
#diffuses = [0.025]
dxlist = [10]

#diffuses = [2.5]
#dxlist = [5]
wdirb = "../../../data/raw/DBSWREVanew100s/o3/"
#for ll in range(3,16):
for ll in dxlist:
    for k in range(len(diffuses)):
        wdir = wdirb + str(ll) + "/" + str(diffuses[k]) + "/"
        if not os.path.exists(wdir):
            os.makedirs(wdir) 
        dx = 10.0 / (2**ll)
        l = 0.01
        dt = l*dx
        startx = 0.0
        endx = 1100.0 + dx
        startt = 0.0
        endt = 100 + (dt*0.9)  
                
        szoomx = startx
        ezoomx = endx
                
        #number of boundary conditions (one side)
        nfcBC = 4 #for flux calculation
        nGsBC = 2 #for solving G from u,h
        niBC = nGsBC + nfcBC #total
                
        g = 9.81
        
        gap = int(1.0/dt)
                
        x,t = makevar(startx,endx,dx,startt,endt,dt)
        n = len(x)
        hf = 1.8#(sqrt(1.8) +1)**2/4.0
        hl = 1.36898
        baseh = hl
        etah = hf - hl
        baseu = 1.07498
        etau = 0.0 - 1.07498
        x0 = 500
        diffuse = diffuses[k]   
        hm,um,Gm = DBSREVsmooth(x,x0,baseh,etah,baseu,etau,diffuse,dx)
        #plot(x,hm)
        #plot(x,um)
        #show()
                        
        umbegi = um[0]*ones(niBC)
        umendi = um[-1]*ones(niBC)
        hmbegi = hm[0]*ones(niBC)
        hmendi = hm[-1]*ones(niBC)   
        Gmbegi = Gm[0]*ones(niBC)
        Gmendi = Gm[-1]*ones(niBC) 
                
        umbeg = umbegi
        umend = umendi
        hmbeg = hmbegi
        hmend = hmendi
        Gmbeg = Gmbegi
        Gmend = Gmendi
                
        #calculate G midpoint
        cnBC = niBC - nGsBC
                
        umbc = concatenate([umbeg[-cnBC:],um,umend[0:cnBC]]) 
        hmbc = concatenate([hmbeg[-cnBC:],hm,hmend[0:cnBC]])       
        Gmbc = concatenate([Gmbeg[-cnBC:],Gm,Gmend[0:cnBC]])  
                
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
            
        Ga_c = copyarraytoC(Ga)
        Gabeg_c = copyarraytoC(Gabeg)
        Gaend_c = copyarraytoC(Gaend)
        ha_c = copyarraytoC(ha)
          
        habeg_c = copyarraytoC(habeg)
        haend_c = copyarraytoC(haend)
          
        uabeg_c = copyarraytoC(uabeg)
        uaend_c = copyarraytoC(uaend)
           
        hmbeg_c = copyarraytoC(hmbeg)
        hmend_c = copyarraytoC(hmend)
            
        umbeg_c = copyarraytoC(umbeg)
        umend_c = copyarraytoC(umend)
           
        u_c= mallocPy(n)
        G_c = mallocPy(n)
        h_c = mallocPy(n)
        
        xbeg = arange(startx - niBC*dx,startx,dx)
        xend = arange(endx + dx,endx + (niBC+1)*dx) 
        
        xbc =  concatenate([xbeg,x,xend])  
        
        xbc_c = copyarraytoC(xbc)
        hbc_c = mallocPy(n + 2*niBC)
        ubc_c = mallocPy(n + 2*niBC)
        Evals = []
            
        for i in range(1,len(t)):  
            
            if(i ==1 or i %gap == 0):
                ca2midpt(ha_c,dx,n,h_c)
                ca2midpt(Ga_c,dx,n,G_c)
                ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
                
                conc(hmbeg_c , h_c,hmend_c,niBC,n ,niBC , hbc_c)
                conc(umbeg_c , u_c,umend_c,niBC,n ,niBC , ubc_c)        
                Eval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)
                
                Evals.append(Eval)
                u = copyarrayfromC(u_c,n)
                G = copyarrayfromC(G_c,n)
                h = copyarrayfromC(h_c,n)
                s = wdir + "out" + str(i) +  ".txt"
                with open(s,'a') as file2:
                     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                    
                     writefile2.writerow(['dx' ,'dt','time','Eval','cell midpoint', 'height(m)', 'G' , 'u(m/s)',"diffuse"])        
                                   
                     for j in range(n):
                         writefile2.writerow([str(dx),str(dt),str(t[i]),str(Eval), str(x[j]), str(h[j]) , str(G[j]) , str(u[j]), str(diffuse)])   
            evolvewrap(Ga_c,ha_c,Gabeg_c,Gaend_c,habeg_c,haend_c,hmbeg_c,hmend_c,uabeg_c,uaend_c,umbeg_c,umend_c,nfcBC,nGsBC,g,dx,dt,n,cnBC,niBC)
            print (t[i])
            
                
        ca2midpt(ha_c,dx,n,h_c)
        ca2midpt(Ga_c,dx,n,G_c)
        ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
        
        conc(hmbeg_c , h_c,hmend_c,niBC,n ,niBC , hbc_c)
        conc(umbeg_c , u_c,umend_c,niBC,n ,niBC , ubc_c)        
        Eval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)
        
        Evals.append(Eval)
        u = copyarrayfromC(u_c,n)
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        s = wdir + "outlast.txt"
        with open(s,'a') as file2:
             writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
             writefile2.writerow(['dx' ,'dt','time','Eval','cell midpoint', 'height(m)', 'G' , 'u(m/s)',"diffuse"])        
                           
             for j in range(n):
                 writefile2.writerow([str(dx),str(dt),str(t[i]),str(Eval), str(x[j]), str(h[j]) , str(G[j]) , str(u[j]), str(diffuse)])    
        deallocPy(u_c)   
        deallocPy(h_c)
        deallocPy(G_c)
        deallocPy(ha_c)
        deallocPy(Ga_c)
        deallocPy(habeg_c)
        deallocPy(Gabeg_c)
        deallocPy(haend_c)
        deallocPy(Gaend_c)
        deallocPy(uabeg_c)
        deallocPy(uaend_c)
        deallocPy(umbeg_c)
        deallocPy(umend_c)
        deallocPy(hmbeg_c)
        deallocPy(hmend_c)
"""

#Dispersive Shockwaves
"""
#diffuses = [0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1.0,2.5,5.0,7.5,10.0,25.0,50.0,75.0,100.0,250.0,500.0,750.0,1000.0]
diffuses = [10]
#diffuses = [0.025]
dxlist = [6]

#diffuses = [2.5]
#dxlist = [5]
wdirb = "../../../data/raw/DBSWanew100s/o3/"
#for ll in range(3,16):
for ll in dxlist:
    for k in range(len(diffuses)):
        wdir = wdirb + str(ll) + "/" + str(diffuses[k]) + "/"
        if not os.path.exists(wdir):
            os.makedirs(wdir) 
        dx = 10.0 / (2**ll)
        l = 0.01
        dt = l*dx
        startx = 0.0
        endx = 800.0 + dx
        startt = 0.0
        endt = 100.0 + (dt*0.9)  
                
        szoomx = startx
        ezoomx = endx
                
        #number of boundary conditions (one side)
        nfcBC = 4 #for flux calculation
        nGsBC = 2 #for solving G from u,h
        niBC = nGsBC + nfcBC #total
                
        g = 9.81
        
        gap = int(1.0/dt)
                
        x,t = makevar(startx,endx,dx,startt,endt,dt)
        n = len(x)
        hf = 1.36898#(sqrt(1.8) +1)**2/4.0
        hl = 1.0
        baseh = hl
        etah = hf - hl
        baseu = 0.0
        etau = 1.07498 - 0.0
        x0 = 300
        diffuse = diffuses[k]   
        hm,um,Gm = DBSsmooth(x,x0,baseh,etah,baseu,etau,diffuse,dx)
        #plot(x,hm)
        #plot(x,um)
        #show()
                        
        umbegi = um[0]*ones(niBC)
        umendi = um[-1]*ones(niBC)
        hmbegi = hm[0]*ones(niBC)
        hmendi = hm[-1]*ones(niBC)   
        Gmbegi = Gm[0]*ones(niBC)
        Gmendi = Gm[-1]*ones(niBC) 
                
        umbeg = umbegi
        umend = umendi
        hmbeg = hmbegi
        hmend = hmendi
        Gmbeg = Gmbegi
        Gmend = Gmendi
                
        #calculate G midpoint
        cnBC = niBC - nGsBC
                
        umbc = concatenate([umbeg[-cnBC:],um,umend[0:cnBC]]) 
        hmbc = concatenate([hmbeg[-cnBC:],hm,hmend[0:cnBC]])       
        Gmbc = concatenate([Gmbeg[-cnBC:],Gm,Gmend[0:cnBC]])  
                
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
            
        Ga_c = copyarraytoC(Ga)
        Gabeg_c = copyarraytoC(Gabeg)
        Gaend_c = copyarraytoC(Gaend)
        ha_c = copyarraytoC(ha)
          
        habeg_c = copyarraytoC(habeg)
        haend_c = copyarraytoC(haend)
          
        uabeg_c = copyarraytoC(uabeg)
        uaend_c = copyarraytoC(uaend)
           
        hmbeg_c = copyarraytoC(hmbeg)
        hmend_c = copyarraytoC(hmend)
            
        umbeg_c = copyarraytoC(umbeg)
        umend_c = copyarraytoC(umend)
           
        u_c= mallocPy(n)
        G_c = mallocPy(n)
        h_c = mallocPy(n)
        
        xbeg = arange(startx - niBC*dx,startx,dx)
        xend = arange(endx + dx,endx + (niBC+1)*dx) 
        
        xbc =  concatenate([xbeg,x,xend])  
        
        xbc_c = copyarraytoC(xbc)
        hbc_c = mallocPy(n + 2*niBC)
        ubc_c = mallocPy(n + 2*niBC)
        Evals = []
            
        for i in range(1,len(t)):  
            
            if(i ==1 or i %gap == 0):
                ca2midpt(ha_c,dx,n,h_c)
                ca2midpt(Ga_c,dx,n,G_c)
                ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
                
                conc(hmbeg_c , h_c,hmend_c,niBC,n ,niBC , hbc_c)
                conc(umbeg_c , u_c,umend_c,niBC,n ,niBC , ubc_c)        
                Eval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)
                
                Evals.append(Eval)
                u = copyarrayfromC(u_c,n)
                G = copyarrayfromC(G_c,n)
                h = copyarrayfromC(h_c,n)
                s = wdir + "out" + str(i) +  ".txt"
                with open(s,'a') as file2:
                     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                    
                     writefile2.writerow(['dx' ,'dt','time','Eval','cell midpoint', 'height(m)', 'G' , 'u(m/s)',"diffuse"])        
                                   
                     for j in range(n):
                         writefile2.writerow([str(dx),str(dt),str(t[i]),str(Eval), str(x[j]), str(h[j]) , str(G[j]) , str(u[j]), str(diffuse)])   
            evolvewrap(Ga_c,ha_c,Gabeg_c,Gaend_c,habeg_c,haend_c,hmbeg_c,hmend_c,uabeg_c,uaend_c,umbeg_c,umend_c,nfcBC,nGsBC,g,dx,dt,n,cnBC,niBC)
            print (t[i])
            
                
        ca2midpt(ha_c,dx,n,h_c)
        ca2midpt(Ga_c,dx,n,G_c)
        ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
        
        conc(hmbeg_c , h_c,hmend_c,niBC,n ,niBC , hbc_c)
        conc(umbeg_c , u_c,umend_c,niBC,n ,niBC , ubc_c)        
        Eval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)
        
        Evals.append(Eval)
        u = copyarrayfromC(u_c,n)
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        s = wdir + "outlast.txt"
        with open(s,'a') as file2:
             writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
             writefile2.writerow(['dx' ,'dt','time','Eval','cell midpoint', 'height(m)', 'G' , 'u(m/s)',"diffuse"])        
                           
             for j in range(n):
                 writefile2.writerow([str(dx),str(dt),str(t[i]),str(Eval), str(x[j]), str(h[j]) , str(G[j]) , str(u[j]), str(diffuse)])    
        deallocPy(u_c)   
        deallocPy(h_c)
        deallocPy(G_c)
        deallocPy(ha_c)
        deallocPy(Ga_c)
        deallocPy(habeg_c)
        deallocPy(Gabeg_c)
        deallocPy(haend_c)
        deallocPy(Gaend_c)
        deallocPy(uabeg_c)
        deallocPy(uaend_c)
        deallocPy(umbeg_c)
        deallocPy(umend_c)
        deallocPy(hmbeg_c)
        deallocPy(hmend_c)
"""



"""
#Dam Break Much Longer
import os
#diffuses = [0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1.0,2.5,5.0,7.5,10.0,25.0,50.0,75.0,100.0,250.0,500.0,750.0,1000.0]
diffuses = [10.0]
dxlist = [8]
wdirb = "../../../data/raw/longcontactdiscdx8diff10fileio/o3/"
#for ll in range(3,16):
for ll in dxlist:
    for k in range(len(diffuses)):
        wdir = wdirb + str(ll) + "/" + str(k) + "/"
        if not os.path.exists(wdir):
            os.makedirs(wdir) 
        dx = 10.0 / (2**ll)
        Cr = 0.5
        g = 9.81
        hf = 1.8
        l = 1.0 / sqrt(g*hf)
        dt = 0.01*dx#Cr*l*dx
        startx = -900
        endx = 1800.0 + dx
        startt = 0.0
        endt = 300.0+(dt*0.9)  
                
        szoomx = startx
        ezoomx = endx
                
        #number of boundary conditions (one side)
        nfcBC = 4 #for flux calculation
        nGsBC = 2 #for solving G from u,h
        niBC = nGsBC + nfcBC #total
                
            
        gap = int(0.1/dt)
                
        x,t = makevar(startx,endx,dx,startt,endt,dt)
        n = len(x)
        hl = 1.0
        base = hl
        eta0 = hf - hl
        x0 = 500
        diffuse = diffuses[k]   
        hm,um = dambreaksmooth(x,x0,base,eta0,diffuse,dx)
        #plot(x,hm)
        #show()

                        
        umbegi = zeros(niBC)
        umendi = zeros(niBC)
        hmbegi = ones(niBC)
        hmendi = ones(niBC)    
                
        for i in range(niBC):
            umbegi[i] = um[0]
            umendi[i] = um[-1]
            hmbegi[i] = hm[0]
            hmendi[i] = hm[-1]
                
        umbeg = umbegi
        umend = umendi
        hmbeg = hmbegi
        hmend = hmendi
                
        #calculate G midpoint
        cnBC = niBC - nGsBC
                
        umbc = concatenate([umbeg[-cnBC:],um,umend[0:cnBC]]) 
        hmbc = concatenate([hmbeg[-cnBC:],hm,hmend[0:cnBC]])       
        Gmbc = solveGfromuh(umbc,hmbc,hmbeg[0:-cnBC],hmend[-cnBC:],umbeg[0:-cnBC],umend[-cnBC:],dx)  
                
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
            
        Ga_c = copyarraytoC(Ga)
        Gabeg_c = copyarraytoC(Gabeg)
        Gaend_c = copyarraytoC(Gaend)
        ha_c = copyarraytoC(ha)
          
        habeg_c = copyarraytoC(habeg)
        haend_c = copyarraytoC(haend)
          
        uabeg_c = copyarraytoC(uabeg)
        uaend_c = copyarraytoC(uaend)
           
        hmbeg_c = copyarraytoC(hmbeg)
        hmend_c = copyarraytoC(hmend)
            
        umbeg_c = copyarraytoC(umbeg)
        umend_c = copyarraytoC(umend)
            
        u_c= mallocPy(n)
        G_c = mallocPy(n)
        h_c = mallocPy(n)
        jtrack = 0    
        for i in range(1,len(t)):  
            if(i == 1 or i%gap ==0):
                ca2midpt(ha_c,dx,n,h_c)
                ca2midpt(Ga_c,dx,n,G_c)
                ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
                u = copyarrayfromC(u_c,n)
                G = copyarrayfromC(G_c,n)
                h = copyarrayfromC(h_c,n)
                s = wdir + "out" +str(jtrack)+".txt"
                with open(s,'a') as file2:
                     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                    
                     writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'G' , 'u(m/s)',"diffuse"])        
                                   
                     for j in range(n):
                         writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(G[j]) , str(u[j]), str(diffuse)])  
                jtrack = jtrack + 1     
                
            evolvewrap(Ga_c,ha_c,Gabeg_c,Gaend_c,habeg_c,haend_c,hmbeg_c,hmend_c,uabeg_c,uaend_c,umbeg_c,umend_c,nfcBC,nGsBC,g,dx,dt,n,cnBC,niBC)
            print (t[i])
            
                
        ca2midpt(ha_c,dx,n,h_c)
        ca2midpt(Ga_c,dx,n,G_c)
        ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
        u = copyarrayfromC(u_c,n)
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        s = wdir + "outlast.txt"
        with open(s,'a') as file2:
             writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
             writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'G' , 'u(m/s)',"diffuse"])        
                           
             for j in range(n):
                 writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(G[j]) , str(u[j]), str(diffuse)])    
        deallocPy(u_c)   
        deallocPy(h_c)
        deallocPy(G_c)
        deallocPy(ha_c)
        deallocPy(Ga_c)
        deallocPy(habeg_c)
        deallocPy(Gabeg_c)
        deallocPy(haend_c)
        deallocPy(Gaend_c)
        deallocPy(uabeg_c)
        deallocPy(uaend_c)
        deallocPy(umbeg_c)
        deallocPy(umend_c)
        deallocPy(hmbeg_c)
        deallocPy(hmend_c)
"""

"""
#Dam Break targetted
difflist = [20]

deltaxa = [1]
dxlist = [deltaxa]

import os
diffuses = [0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1.0,2.5,5.0,7.5,10.0,25.0,50.0,75.0,100.0,250.0,500.0,750.0,1000.0]
wdirb = "../../data/raw/bigsmoothTEST/o3/"
for lk in range(len(difflist)):
    for ll in dxlist[lk]:
        wdir = wdirb + str(ll) + "/" + str(difflist[lk]) + "/"
        dx = ll*(10.0 / (2**10))
        if not os.path.exists(wdir):
            os.makedirs(wdir) 
        l = 0.01
        dt = l*dx
        startx = -100
        endx = 1100.0 + dx
        startt = 0.0
        endt = 100.0+(dt*0.9)  
                
        szoomx = startx
        ezoomx = endx
                
        #number of boundary conditions (one side)
        nfcBC = 4 #for flux calculation
        nGsBC = 2 #for solving G from u,h
        niBC = nGsBC + nfcBC #total
                
        g = 9.81
            
        gap = int(5.0/dt)
                
        x,t = makevar(startx,endx,dx,startt,endt,dt)
        n = len(x)
        hf = 1.8
        hl = 1.0
        base = hl
        eta0 = hf - hl
        x0 = 500
        diffuse = diffuses[difflist[lk]]  
        hm,um = dambreaksmooth(x,x0,base,eta0,diffuse,dx)
                        
        umbegi = zeros(niBC)
        umendi = zeros(niBC)
        hmbegi = ones(niBC)
        hmendi = ones(niBC)    
                
        for i in range(niBC):
            umbegi[i] = um[0]
            umendi[i] = um[-1]
            hmbegi[i] = hm[0]
            hmendi[i] = hm[-1]
                
        umbeg = umbegi
        umend = umendi
        hmbeg = hmbegi
        hmend = hmendi
                
        #calculate G midpoint
        cnBC = niBC - nGsBC
                
        umbc = concatenate([umbeg[-cnBC:],um,umend[0:cnBC]]) 
        hmbc = concatenate([hmbeg[-cnBC:],hm,hmend[0:cnBC]])       
        Gmbc = solveGfromuh(umbc,hmbc,hmbeg[0:-cnBC],hmend[-cnBC:],umbeg[0:-cnBC],umend[-cnBC:],dx)  
                
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
            
        Ga_c = copyarraytoC(Ga)
        Gabeg_c = copyarraytoC(Gabeg)
        Gaend_c = copyarraytoC(Gaend)
        ha_c = copyarraytoC(ha)
          
        habeg_c = copyarraytoC(habeg)
        haend_c = copyarraytoC(haend)
          
        uabeg_c = copyarraytoC(uabeg)
        uaend_c = copyarraytoC(uaend)
           
        hmbeg_c = copyarraytoC(hmbeg)
        hmend_c = copyarraytoC(hmend)
            
        umbeg_c = copyarraytoC(umbeg)
        umend_c = copyarraytoC(umend)
            
        u_c= mallocPy(n)
        G_c = mallocPy(n)
        h_c = mallocPy(n)
            
        for i in range(1,len(t)):  
            if(i % gap == 0 or i ==1):
                
                ca2midpt(ha_c,dx,n,h_c)
                ca2midpt(Ga_c,dx,n,G_c)
                ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
                u = copyarrayfromC(u_c,n)
                G = copyarrayfromC(G_c,n)
                h = copyarrayfromC(h_c,n)
                s = wdir +"saveout" + str(i/gap) +".txt"
                with open(s,'a') as file2:
                     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                    
                     writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'G' , 'u(m/s)',"diffuse"])        
                                   
                     for j in range(n):
                         writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(G[j]) , str(u[j]), str(diffuse)]) 
            evolvewrap(Ga_c,ha_c,Gabeg_c,Gaend_c,habeg_c,haend_c,hmbeg_c,hmend_c,uabeg_c,uaend_c,umbeg_c,umend_c,nfcBC,nGsBC,g,dx,dt,n,cnBC,niBC)
            print (t[i])
            
                
        ca2midpt(ha_c,dx,n,h_c)
        ca2midpt(Ga_c,dx,n,G_c)
        ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
        u = copyarrayfromC(u_c,n)
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        s = wdir + "outlast.txt"
        with open(s,'a') as file2:
             writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
             writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'G' , 'u(m/s)',"diffuse"])        
                           
             for j in range(n):
                 writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(G[j]) , str(u[j]), str(diffuse)])    
        deallocPy(u_c)   
        deallocPy(h_c)
        deallocPy(G_c)
        deallocPy(ha_c)
        deallocPy(Ga_c)
        deallocPy(habeg_c)
        deallocPy(Gabeg_c)
        deallocPy(haend_c)
        deallocPy(Gaend_c)
        deallocPy(uabeg_c)
        deallocPy(uaend_c)
        deallocPy(umbeg_c)
        deallocPy(umend_c)
        deallocPy(hmbeg_c)
        deallocPy(hmend_c)
"""

"""
#Dam Break 
from time import time
wdir = "../../data/trackleadsol/o3/"
hf = 1.8
hl = 1.0
g = 9.81
dx = 10.0 /(2**10)
l = 0.01
dt = l*dx
startx = 0.0
endx = 1000.0 + dx
startt = 0.0
endt = 100 + dt

if not os.path.exists(wdir):
    os.makedirs(wdir) 
        
szoomx = startx
ezoomx = endx
        
#number of boundary conditions (one side)
nfcBC = 4 #for flux calculation
nGsBC = 2 #for solving G from u,h
niBC = nGsBC + nfcBC #total
        
g = 9.81
    
gap = int(0.5/dt)
        
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
      
base = hl
eta0 = hf - hl
x0 = 500
diffuse = 1000
hm,um = dambreaksmooth(x,x0,base,eta0,diffuse,dx)
                
umbegi = zeros(niBC)
umendi = zeros(niBC)
hmbegi = ones(niBC)
hmendi = ones(niBC)    
        
for i in range(niBC):
    umbegi[i] = um[0]
    umendi[i] = um[-1]
    hmbegi[i] = hm[0]
    hmendi[i] = hm[-1]
        
umbeg = umbegi
umend = umendi
hmbeg = hmbegi
hmend = hmendi
        
#calculate G midpoint
cnBC = niBC - nGsBC
        
umbc = concatenate([umbeg[-cnBC:],um,umend[0:cnBC]]) 
hmbc = concatenate([hmbeg[-cnBC:],hm,hmend[0:cnBC]])       
Gmbc = solveGfromuh(umbc,hmbc,hmbeg[0:-cnBC],hmend[-cnBC:],umbeg[0:-cnBC],umend[-cnBC:],dx)  
        
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
    
Ga_c = copyarraytoC(Ga)
Gabeg_c = copyarraytoC(Gabeg)
Gaend_c = copyarraytoC(Gaend)
ha_c = copyarraytoC(ha)
  
habeg_c = copyarraytoC(habeg)
haend_c = copyarraytoC(haend)
  
uabeg_c = copyarraytoC(uabeg)
uaend_c = copyarraytoC(uaend)
   
hmbeg_c = copyarraytoC(hmbeg)
hmend_c = copyarraytoC(hmend)
    
umbeg_c = copyarraytoC(umbeg)
umend_c = copyarraytoC(umend)
    
u_c= mallocPy(n)
G_c = mallocPy(n)
h_c = mallocPy(n)

hm_c = copyarraytoC(hm)
um_c = copyarraytoC(um)

xbeg = arange(startx - niBC*dx,startx,dx)
xend = arange(endx + dx,endx + (niBC+1)*dx) 

xbc =  concatenate([xbeg,x,xend])  

xbc_c = copyarraytoC(xbc)
hbc_c = mallocPy(n + 2*niBC)
ubc_c = mallocPy(n + 2*niBC)
Evals = []
aplus = []
aplusx = []
aplust = []
    
for i in range(1,len(t)):

    if(i ==1 or i %gap == 0):
        ca2midpt(ha_c,dx,n,h_c)
        ca2midpt(Ga_c,dx,n,G_c)
        conc(hmbeg_c , h_c,hmend_c,niBC,n ,niBC , hbc_c)
        conc(umbeg_c , u_c,umend_c,niBC,n ,niBC , ubc_c)        
        Eval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)
        
        Evals.append(Eval)
        ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
        u = copyarrayfromC(u_c,n)
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        
        mi = n - 2
        for mi in range(n-1,-1,-1):
            if(h[mi -1] < h[mi]) and (h[mi] > 1.1 ):
                break
        aplus.append(h[mi])
        aplusx.append(x[mi])
        aplust.append(t[i])
        s = wdir + "out" + str(i) +  ".txt"
        with open(s,'a') as file2:
             writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
             writefile2.writerow(['dx' ,'dt','time','Eval','cell midpoint', 'height(m)', 'G' , 'u(m/s)'])        
                           
             for j in range(n):
                 writefile2.writerow([str(dx),str(dt),str(t[i]),str(Eval), str(x[j]), str(h[j]) , str(G[j]) , str(u[j])])   
                 
    evolvewrap(Ga_c,ha_c,Gabeg_c,Gaend_c,habeg_c,haend_c,hmbeg_c,hmend_c,uabeg_c,uaend_c,umbeg_c,umend_c,nfcBC,nGsBC,g,dx,dt,n,cnBC,niBC)
    print (t[i])
    
ha1 = copyarrayfromC(ha_c,n)
Ga1 = copyarrayfromC(Ga_c,n)        
ca2midpt(ha_c,dx,n,h_c)
ca2midpt(Ga_c,dx,n,G_c)

conc(hmbeg_c , h_c,hmend_c,niBC,n ,niBC , hbc_c)
conc(umbeg_c , u_c,umend_c,niBC,n ,niBC , ubc_c)        
Eval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)

Evals.append(Eval)

ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
u = copyarrayfromC(u_c,n)
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)

mi = n - 2
for mi in range(n-1,-1,-1):
    if(h[mi -1] < h[mi]) and (h[mi] > 1.1 ):
        break
aplus.append(h[mi])
aplusx.append(x[mi])
aplust.append(t[i])
ha2 = copyarrayfromC(ha_c,n)
Ga2 = copyarrayfromC(Ga_c,n) 
s = wdir + "outlast.txt"
with open(s,'a') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
        writefile2.writerow(['dx' ,'dt','time','Eval','cell midpoint', 'height(m)', 'G' , 'u(m/s)'])        
                           
        for j in range(n):
                writefile2.writerow([str(dx),str(dt),str(t[i]),str(Eval), str(x[j]), str(h[j]) , str(G[j]) , str(u[j])])   
                
s = wdir + "aplus.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(['x' ,'t','aplus' ,"Grim"])        
           
    for j in range(len(aplus)):
        writefile2.writerow([str(aplusx[j]),str(aplust[j]),str(aplus[j]),str(0.739976603390100695296990254)]) 
                 

deallocPy(u_c)   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(ha_c)
deallocPy(Ga_c)
deallocPy(habeg_c)
deallocPy(Gabeg_c)
deallocPy(haend_c)
deallocPy(Gaend_c)
deallocPy(uabeg_c)
deallocPy(uaend_c)
deallocPy(umbeg_c)
deallocPy(umend_c)
deallocPy(hmbeg_c)
deallocPy(hmend_c)
"""



################# DAM BREAK Accuracy
"""
wdir = "../../../data/raw/dbh/o3/"
for k in range(20):
    g = 9.81
    hf = 1.8
    hl = 1.0 
    dx = 100.0 / (2**k)
    Cr = 0.5
    l = Cr / sqrt(g*hf)
    dt = l*dx
    dt = l*dx
    startx = 0.0
    endx = 1000.0 + dx
    startt = 0.0
    endt = 30 + dt
            
    szoomx = startx
    ezoomx = endx
            
    #number of boundary conditions (one side)
    nfcBC = 4 #for flux calculation
    nGsBC = 2 #for solving G from u,h
    niBC = nGsBC + nfcBC #total

    wdatadir = wdir+ str(k) + "/"   
    if not os.path.exists(wdatadir):
        os.makedirs(wdatadir)
        
    gap = max(1,int(0.5/dt))
            
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
        
    xc = 500   
    um,hm = dambreak(x,xc,hf,hl)
                    
    umbegi = zeros(niBC)
    umendi = zeros(niBC)
    hmbegi = ones(niBC)
    hmendi = ones(niBC)    
            
    for i in range(niBC):
        umbegi[i] = um[0]
        umendi[i] = um[-1]
        hmbegi[i] = hm[0]
        hmendi[i] = hm[-1]
            
    umbeg = umbegi
    umend = umendi
    hmbeg = hmbegi
    hmend = hmendi
            
    #calculate G midpoint
    cnBC = niBC - nGsBC
            
    umbc = concatenate([umbeg[-cnBC:],um,umend[0:cnBC]]) 
    hmbc = concatenate([hmbeg[-cnBC:],hm,hmend[0:cnBC]])       
    Gmbc = solveGfromuh(umbc,hmbc,hmbeg[0:-cnBC],hmend[-cnBC:],umbeg[0:-cnBC],umend[-cnBC:],dx)  
            
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
        
    Ga_c = copyarraytoC(Ga)
    Gabeg_c = copyarraytoC(Gabeg)
    Gaend_c = copyarraytoC(Gaend)
    ha_c = copyarraytoC(ha)
      
    habeg_c = copyarraytoC(habeg)
    haend_c = copyarraytoC(haend)
      
    uabeg_c = copyarraytoC(uabeg)
    uaend_c = copyarraytoC(uaend)
       
    hmbeg_c = copyarraytoC(hmbeg)
    hmend_c = copyarraytoC(hmend)
        
    umbeg_c = copyarraytoC(umbeg)
    umend_c = copyarraytoC(umend)
        
    u_c= mallocPy(n)
    G_c = mallocPy(n)
    h_c = mallocPy(n)
        
    for i in range(1,len(t)):
        if(i ==1 or i %gap == 0):
            ca2midpt(ha_c,dx,n,h_c)
            ca2midpt(Ga_c,dx,n,G_c)
            ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
            u = copyarrayfromC(u_c,n)
            G = copyarrayfromC(G_c,n)
            h = copyarrayfromC(h_c,n)
            s = wdatadir + "out" + str(i) +  ".txt"
            with open(s,'a') as file2:
                 writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                
                 writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'G' , 'u(m/s)'])        
                               
                 for j in range(n):
                     writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(G[j]) , str(u[j])])    
        evolvewrap(Ga_c,ha_c,Gabeg_c,Gaend_c,habeg_c,haend_c,hmbeg_c,hmend_c,uabeg_c,uaend_c,umbeg_c,umend_c,nfcBC,nGsBC,g,dx,dt,n,cnBC,niBC)
        print (t[i])
        
            
    ca2midpt(ha_c,dx,n,h_c)
    ca2midpt(Ga_c,dx,n,G_c)
    ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
    u = copyarrayfromC(u_c,n)
    G = copyarrayfromC(G_c,n)
    h = copyarrayfromC(h_c,n)
    s = wdatadir + "outlast.txt"
    with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
         writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'G' , 'u(m/s)'])        
                       
         for j in range(n):
             writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(G[j]) , str(u[j])])    
             
    s = wdir + str(k)+ ".txt"
    with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
         writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'G' , 'u(m/s)'])        
                   
         for j in range(n):
             writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(G[j]) , str(u[j])])    
    deallocPy(u_c)   
    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(ha_c)
    deallocPy(Ga_c)
    deallocPy(habeg_c)
    deallocPy(Gabeg_c)
    deallocPy(haend_c)
    deallocPy(Gaend_c)
    deallocPy(uabeg_c)
    deallocPy(uaend_c)
    deallocPy(umbeg_c)
    deallocPy(umend_c)
    deallocPy(hmbeg_c)
    deallocPy(hmend_c)
"""


######### Solitary Soliton
"""
dx = 0.001
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
nGsBC = 2 #for solving G from u,h
niBC = nGsBC + nfcBC #total
    
    
wdir = "../../../data/raw/Cserre/Evalt/o3/dx0p05"

if not os.path.exists(wdir):
    os.makedirs(wdir)
    
gap = int(10.0/dt)
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)

t0 = 0   
hm,um = solitoninit(n,a0,a1,g,x,t0,dx)
            
umbeg = um[0]*ones(niBC)
umend = um[-1]*ones(niBC)
hmbeg = hm[0]*ones(niBC)
hmend = hm[-1]*ones(niBC)    
            
#calculate G midpoint
cnBC = niBC - nGsBC
    
umbc = concatenate([umbeg[-cnBC:],um,umend[0:cnBC]]) 
hmbc = concatenate([hmbeg[-cnBC:],hm,hmend[0:cnBC]])       
Gmbc = solveGfromuh(umbc,hmbc,hmbeg[0:-cnBC],hmend[-cnBC:],umbeg[0:-cnBC],umend[-cnBC:],dx)  
    
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

Ga_c = copyarraytoC(Ga)
Gabeg_c = copyarraytoC(Gabeg)
Gaend_c = copyarraytoC(Gaend)
ha_c = copyarraytoC(ha)

habeg_c = copyarraytoC(habeg)
haend_c = copyarraytoC(haend)

uabeg_c = copyarraytoC(uabeg)
uaend_c = copyarraytoC(uaend)

hmbeg_c = copyarraytoC(hmbeg)
hmend_c = copyarraytoC(hmend)

umbeg_c = copyarraytoC(umbeg)
umend_c = copyarraytoC(umend)

u_c = mallocPy(n)
G_c = mallocPy(n)
h_c = mallocPy(n)

xbeg = arange(startx - niBC*dx,startx,dx)
xend = arange(endx + dx,endx + (niBC+1)*dx) 

xbc =  concatenate([xbeg,x,xend])  

xbc_c = copyarraytoC(xbc)
hbc_c = mallocPy(n + 2*niBC)
ubc_c = mallocPy(n + 2*niBC)
Evals = []

for i in range(1,len(t)):
    
    if(i % gap == 0 or i ==1):
        
        ca2midpt(ha_c,dx,n,h_c)
        ca2midpt(Ga_c,dx,n,G_c)
        
        conc(hmbeg_c , h_c,hmend_c,niBC,n ,niBC , hbc_c)
        conc(umbeg_c , u_c,umend_c,niBC,n ,niBC , ubc_c)        
        Eval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)
        
        Evals.append(Eval)
        
        ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
        u = copyarrayfromC(u_c,n)
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        
        
        c = sqrt(g*(a0 + a1))
        htrue = zeros(n)
        utrue = zeros(n)
        for j in range(n):             
            he = soliton(x[j],t[i],g,a0,a1)
            htrue[j] = he
            utrue[j] = c* ((he - a0) / he) 
            
        s = wdir + "saveoutputts" + str(i) + ".txt"
        print t[i]
        print(h[30],G[30]) 
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

            writefile2.writerow(['dx' ,'dt','time', 'height(m)', 'G' , 'u(m/s)','true height', 'true velocity'  ])        
               
            for j in range(n):
                writefile2.writerow([str(dx),str(dt),str(t[i]), str(h[j]) , str(G[j]) , str(u[j]), str(htrue[j]), str(utrue[j])])  
             
        
    evolvewrap(Ga_c,ha_c,Gabeg_c,Gaend_c,habeg_c,haend_c,hmbeg_c,hmend_c,uabeg_c,uaend_c,umbeg_c,umend_c,nfcBC,nGsBC,g,dx,dt,n,cnBC,niBC)
    print t[i]
    print(h[30],G[30]) 

    
ca2midpt(ha_c,dx,n,h_c)
ca2midpt(Ga_c,dx,n,G_c)
ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
u = copyarrayfromC(u_c,n)
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)

c = sqrt(g*(a0 + a1))
htrue = zeros(n)
utrue = zeros(n)
for j in range(n):             
    he = soliton(x[j],t[i],g,a0,a1)
    htrue[j] = he
    utrue[j] = c* ((he - a0) / he) 
    
s = wdir + "saveoutputtslast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['dx' ,'dt','time', 'height(m)', 'G' , 'u(m/s)','true height', 'true velocity' ])        
               
     for j in range(n):
         writefile2.writerow([str(dx),str(dt),str(t[i]), str(h[j]) , str(G[j]) , str(u[j]), str(htrue[j]), str(utrue[j])])  

for i in range(1,len(Evals)):
    print(Evals[0] - Evals[i])
"""
"""
######### Colliding Soliton Experiments

wdir = "../../../data/raw/Cserre/solitonint/collnonlindx0p1N/o3/"

if not os.path.exists(wdir):
    os.makedirs(wdir)

dx = 0.1

a0 = 1.0
a11 = 1.0
solbeg1 = 100.0
solend1 = 200.0
direction1 = 1.0
a12 = 1.0
solbeg2 = 200.0
solend2 = 300.0
direction2 = -1.0

Cr = 0.5
g = 9.81
#g = 1.0
l = Cr / (sqrt(g*1.5*(a0 + a11 + a12)))
dt = l*dx
startx = -100.0
endx = 500.0 + dx
startt = 0.0
endt = 50 + dt

    
#number of boundary conditions (one side)
nfcBC = 4 #for flux calculation
nGsBC = 2 #for solving G from u,h
niBC = nGsBC + nfcBC #total
    
gap = int(0.5/dt)
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)

t0 = 0   
hm,um = soliton2interactinit(n,a0,a11,solbeg1,solend1,direction1,a12,solbeg2,solend2,direction2,g,x,t0,dx)
          
umbeg = um[0]*ones(niBC)
umend = um[-1]*ones(niBC)
hmbeg = hm[0]*ones(niBC)
hmend = hm[-1]*ones(niBC)    
            
#calculate G midpoint
cnBC = niBC - nGsBC
    
umbc = concatenate([umbeg[-cnBC:],um,umend[0:cnBC]]) 
hmbc = concatenate([hmbeg[-cnBC:],hm,hmend[0:cnBC]])       
Gmbc = solveGfromuh(umbc,hmbc,hmbeg[0:-cnBC],hmend[-cnBC:],umbeg[0:-cnBC],umend[-cnBC:],dx)  
    
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

Ga_c = copyarraytoC(Ga)
Gabeg_c = copyarraytoC(Gabeg)
Gaend_c = copyarraytoC(Gaend)
ha_c = copyarraytoC(ha)

habeg_c = copyarraytoC(habeg)
haend_c = copyarraytoC(haend)

uabeg_c = copyarraytoC(uabeg)
uaend_c = copyarraytoC(uaend)

hmbeg_c = copyarraytoC(hmbeg)
hmend_c = copyarraytoC(hmend)

umbeg_c = copyarraytoC(umbeg)
umend_c = copyarraytoC(umend)

u_c = mallocPy(n)
G_c = mallocPy(n)
h_c = mallocPy(n)

hm_c = copyarraytoC(hm)
um_c = copyarraytoC(um)

xbeg = arange(startx - niBC*dx,startx,dx)
xend = arange(endx + dx,endx + (niBC+1)*dx) 

xbc =  concatenate([xbeg,x,xend])  

xbc_c = copyarraytoC(xbc)
hbc_c = mallocPy(n + 2*niBC)
ubc_c = mallocPy(n + 2*niBC)
Evals = []

for i in range(1,len(t)):
    
    if(i % gap == 0 or i ==1):
        
        ca2midpt(ha_c,dx,n,h_c)
        ca2midpt(Ga_c,dx,n,G_c)
        
        ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
        
        conc(hmbeg_c , h_c,hmend_c,niBC,n ,niBC , hbc_c)
        conc(umbeg_c , u_c,umend_c,niBC,n ,niBC , ubc_c)        
        Eval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)
        
        Evals.append(Eval)
        
        u = copyarrayfromC(u_c,n)
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
            
        s = wdir + "saveoutputts" + str(i) + ".txt"
        print t[i]
        print(h[30],G[30]) 
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

            writefile2.writerow(['dx' ,'dt','time','Eval','x value', 'height(m)', 'G' , 'u(m/s)','true height', 'true velocity' ])        
                   
            for j in range(n):
                writefile2.writerow([str(dx),str(dt),str(t[i]),str(Eval), str(x[j]), str(h[j]) , str(G[j]) , str(u[j])]) 
             
        
    evolvewrap(Ga_c,ha_c,Gabeg_c,Gaend_c,habeg_c,haend_c,hmbeg_c,hmend_c,uabeg_c,uaend_c,umbeg_c,umend_c,nfcBC,nGsBC,g,dx,dt,n,cnBC,niBC)
    print t[i]
    print(h[30],G[30]) 

    
ca2midpt(ha_c,dx,n,h_c)
ca2midpt(Ga_c,dx,n,G_c)

ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)

conc(hmbeg_c , h_c,hmend_c,niBC,n ,niBC , hbc_c)
conc(umbeg_c , u_c,umend_c,niBC,n ,niBC , ubc_c)        
Eval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)

Evals.append(Eval)

u = copyarrayfromC(u_c,n)
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)
    
s = wdir + "saveoutputtslast.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(['dx' ,'dt','time','Eval','x value', 'height(m)', 'G' , 'u(m/s)','true height', 'true velocity' ])        
           
    for j in range(n):
        writefile2.writerow([str(dx),str(dt),str(t[i]),str(Eval), str(x[j]), str(h[j]) , str(G[j]) , str(u[j])]) 
"""

"""
#Soliton Time
from time import time
import os
from scipy.linalg import norm
a0 = 1.0
a1 = 10.0
dx = 0.1#100.0 / (2.0**5)
Cr = 0.2
g = 9.81
l = 1 / (sqrt(g*(a0 + a1)))
dt = Cr*l*dx
startx = -200.0
endx = 800.0
startt = 0.0
endt = 50.0 + dt
    
szoomx = startx
ezoomx = endx
    
#number of boundary conditions (one side)
nfcBC = 4 #for flux calculation
nGsBC = 2 #for solving G from u,h
niBC = nGsBC + nfcBC #total
    
    
wdir = "../../../data/raw/timecomplong/o3/"

if not os.path.exists(wdir):
    os.makedirs(wdir) 

gap = int(0.5/dt)
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)

t0 = 0   
hm,um = solitoninit(n,a0,a1,g,x,t0,dx)
            
umbeg = um[0]*ones(niBC)
umend = um[-1]*ones(niBC)
hmbeg = hm[0]*ones(niBC)
hmend = hm[-1]*ones(niBC)    
            
#calculate G midpoint
cnBC = niBC - nGsBC
    
umbc = concatenate([umbeg[-cnBC:],um,umend[0:cnBC]]) 
hmbc = concatenate([hmbeg[-cnBC:],hm,hmend[0:cnBC]])       
Gmbc = solveGfromuh(umbc,hmbc,hmbeg[0:-cnBC],hmend[-cnBC:],umbeg[0:-cnBC],umend[-cnBC:],dx)  
    
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

Ga_c = copyarraytoC(Ga)
Gabeg_c = copyarraytoC(Gabeg)
Gaend_c = copyarraytoC(Gaend)
ha_c = copyarraytoC(ha)

habeg_c = copyarraytoC(habeg)
haend_c = copyarraytoC(haend)

uabeg_c = copyarraytoC(uabeg)
uaend_c = copyarraytoC(uaend)

hmbeg_c = copyarraytoC(hmbeg)
hmend_c = copyarraytoC(hmend)

umbeg_c = copyarraytoC(umbeg)
umend_c = copyarraytoC(umend)

u_c = mallocPy(n)
G_c = mallocPy(n)
h_c = mallocPy(n)

t0 = time()
for i in range(1,len(t)):    
    evolvewrap(Ga_c,ha_c,Gabeg_c,Gaend_c,habeg_c,haend_c,hmbeg_c,hmend_c,uabeg_c,uaend_c,umbeg_c,umend_c,nfcBC,nGsBC,g,dx,dt,n,cnBC,niBC)
    print t[i] 
t1 = time()
timeelapse = t1 - t0
    
ca2midpt(ha_c,dx,n,h_c)
ca2midpt(Ga_c,dx,n,G_c)
ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
u = copyarrayfromC(u_c,n)
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)

c = sqrt(g*(a0 + a1))
htrue = zeros(n)
utrue = zeros(n)
for j in range(n):             
    he = soliton(x[j],t[i],g,a0,a1)
    htrue[j] = he
    utrue[j] = c* ((he - a0) / he) 
normh = norm(h - htrue,ord=1) / norm(htrue,ord=1)    
normu = norm(u - utrue,ord=1) / norm(utrue,ord=1) 
s = wdir + "saveoutputtslast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['dx' ,'dt','time', 'height(m)', 'G' , 'u(m/s)','true height', 'true velocity',"hnorm","unorm", "time taken", "number of steps", "average time per step" ])        
               
     for j in range(n):
         writefile2.writerow([str(dx),str(dt),str(t[i]), str(h[j]) , str(G[j]) , str(u[j]), str(htrue[j]), str(utrue[j]),str(normh),str(normu), str(timeelapse),str(len(t) - 1) , str((1.0*timeelapse)/ (len(t) - 1) )])
"""

"""
##Accuracy Test
### Soliton Accuracy ################
wdir = "../../../data/raw/Solnon0p7/o3/"

if not os.path.exists(wdir):
    os.makedirs(wdir)

s = wdir + "savenorms.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(['dx','Normalised L1-norm Difference Height', ' Normalised L1-norm Difference Velocity','Hamiltonian Difference'])
    
for k in range(6,21):
    dx = 100.0 / (2**k)
    a0 = 1.0
    a1 = 0.7
    g = 9.81
    #g = 1.0
    Cr = 0.5
    l = 1.0 / (sqrt(g*(a0 + a1)))
    dt = Cr*l*dx
    startx = -250.0
    endx = 250.0 + dx
    startt = 0.0
    endt = 50 + dt
    
    print(dx,dt)
        
    szoomx = startx
    ezoomx = endx
        
    #number of boundary conditions (one side)
    nfcBC = 4 #for flux calculation
    nGsBC = 2 #for solving G from u,h
    niBC = nGsBC + nfcBC #total
        
        
    wdatadir = wdir+ str(k) + "/" 
    if not os.path.exists(wdatadir):
        os.makedirs(wdatadir)
    
    gap =int(10.0/dt)
        
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
    
    t0 = 0.0   
    hm,um = solitoninit(n,a0,a1,g,x,t0,dx)
                
    umbeg = um[0]*ones(niBC)
    umend = um[-1]*ones(niBC)
    hmbeg = hm[0]*ones(niBC)
    hmend = hm[-1]*ones(niBC)    
                
    #calculate G midpoint
    cnBC = niBC - nGsBC
        
    umbc = concatenate([umbeg[-cnBC:],um,umend[0:cnBC]]) 
    hmbc = concatenate([hmbeg[-cnBC:],hm,hmend[0:cnBC]])       
    Gmbc = solveGfromuh(umbc,hmbc,hmbeg[0:-cnBC],hmend[-cnBC:],umbeg[0:-cnBC],umend[-cnBC:],dx)  
        
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
    
    Ga_c = copyarraytoC(Ga)
    Gabeg_c = copyarraytoC(Gabeg)
    Gaend_c = copyarraytoC(Gaend)
    ha_c = copyarraytoC(ha)
    
    habeg_c = copyarraytoC(habeg)
    haend_c = copyarraytoC(haend)
    
    uabeg_c = copyarraytoC(uabeg)
    uaend_c = copyarraytoC(uaend)
    
    hmbeg_c = copyarraytoC(hmbeg)
    hmend_c = copyarraytoC(hmend)
    
    umbeg_c = copyarraytoC(umbeg)
    umend_c = copyarraytoC(umend)
    
    u_c = mallocPy(n)
    G_c = mallocPy(n)
    h_c = mallocPy(n)
    
    xbeg = arange(startx - niBC*dx,startx,dx)
    xend = arange(endx + dx,endx + (niBC+1)*dx) 
    
    xbc =  concatenate([xbeg,x,xend])  
    
    xbc_c = copyarraytoC(xbc)
    hbc_c = mallocPy(n + 2*niBC)
    ubc_c = mallocPy(n + 2*niBC)
    Evals = []
    
    for i in range(1,len(t)):
        
        if(i % gap == 0 or i ==1):
            
            ca2midpt(ha_c,dx,n,h_c)
            ca2midpt(Ga_c,dx,n,G_c)
            ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
            
            conc(hmbeg_c , h_c,hmend_c,niBC,n ,niBC , hbc_c)
            conc(umbeg_c , u_c,umend_c,niBC,n ,niBC , ubc_c)        
            Eval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)
            
            Evals.append(Eval)
            u = copyarrayfromC(u_c,n)
            G = copyarrayfromC(G_c,n)
            h = copyarrayfromC(h_c,n)
            
            
            c = sqrt(g*(a0 + a1))
            htrue = zeros(n)
            utrue = zeros(n)
            for j in range(n):             
                he = soliton(x[j],t[i],g,a0,a1)
                htrue[j] = he
                utrue[j] = c* ((he - a0) / he) 
                
            s = wdatadir + "saveoutputts" + str(i) + ".txt"
            print t[i]
            print(h[3],G[3]) 
            with open(s,'a') as file2:
                writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
                writefile2.writerow(['dx' ,'dt','time','xi','Eval', 'height(m)', 'G' , 'u(m/s)','true height', 'true velocity'  ])        
                   
                for j in range(n):
                    writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]),str(Eval), str(h[j]) , str(G[j]) , str(u[j]), str(htrue[j]), str(utrue[j])])  
                 
            
        evolvewrap(Ga_c,ha_c,Gabeg_c,Gaend_c,habeg_c,haend_c,hmbeg_c,hmend_c,uabeg_c,uaend_c,umbeg_c,umend_c,nfcBC,nGsBC,g,dx,dt,n,cnBC,niBC)
        print t[i]
        print(h[3],G[3]) 
    
        
    ca2midpt(ha_c,dx,n,h_c)
    ca2midpt(Ga_c,dx,n,G_c)
    ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
    
    conc(hmbeg_c , h_c,hmend_c,niBC,n ,niBC , hbc_c)
    conc(umbeg_c , u_c,umend_c,niBC,n ,niBC , ubc_c)        
    Eval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)
    
    Evals.append(Eval)
    u = copyarrayfromC(u_c,n)
    G = copyarrayfromC(G_c,n)
    h = copyarrayfromC(h_c,n)
    
    c = sqrt(g*(a0 + a1))
    htrue = zeros(n)
    utrue = zeros(n)
    for j in range(n):             
        he = soliton(x[j],t[i],g,a0,a1)
        htrue[j] = he
        utrue[j] = c* ((he - a0) / he) 
        
    s = wdatadir + "saveoutputtslast.txt"
    with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
         writefile2.writerow(['dx' ,'dt','time','xi','Eval', 'height(m)', 'G' , 'u(m/s)','true height', 'true velocity'  ])        
                   
         for j in range(n):
             writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]),str(Eval), str(h[j]) , str(G[j]) , str(u[j]), str(htrue[j]), str(utrue[j])])  
     
    normhdiffi = norm(h - htrue,ord=1) / norm(htrue,ord=1)
    normudiffi = norm(u -utrue,ord=1) / norm(utrue,ord=1)
    normHamdiff = (Evals[-1] - Evals[0])/ Evals[0]

    s = wdir + "savenorms.txt"
    with open(s,'a') as file1:
        writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

        writefile.writerow([str(dx),str(normhdiffi), str(normudiffi),str(normHamdiff)])
"""


"""
##### SEGUR AND HAMMACK EXPERIMENT 1
tl = 60.0
b = 0.61
h0 = 0.09
h1 = 0.1
g = 9.81

dx = 0.01
Cr = 0.5
l = Cr / sqrt(g*h1)
dt = l*dx
startx = -tl
endx = tl + dx
startt = 0
endt = 50 + dt  

wdir = "../../../data/raw/segur/o3/"

if not os.path.exists(wdir):
    os.makedirs(wdir)


#number of boundary conditions (one side)
nfcBC = 4 #for flux calculation
nGsBC = 2 #for solving G from u,h
niBC = nGsBC + nfcBC #total
        
    
gap = 1
        
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
    
  
hm,um = experiment1(x,b,h0,h1,dx)
plot(x,hm)
xlim([-2,2])
               
umbegi = zeros(niBC)
umendi = zeros(niBC)
hmbegi = ones(niBC)
hmendi = ones(niBC)    
        
for i in range(niBC):
    umbegi[i] = um[0]
    umendi[i] = um[-1]
    hmbegi[i] = hm[0]
    hmendi[i] = hm[-1]
        
umbeg = umbegi
umend = umendi
hmbeg = hmbegi
hmend = hmendi
        
#calculate G midpoint
cnBC = niBC - nGsBC
        
umbc = concatenate([umbeg[-cnBC:],um,umend[0:cnBC]]) 
hmbc = concatenate([hmbeg[-cnBC:],hm,hmend[0:cnBC]])       
Gmbc = solveGfromuh(umbc,hmbc,hmbeg[0:-cnBC],hmend[-cnBC:],umbeg[0:-cnBC],umend[-cnBC:],dx)  
        
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
    
Ga_c = copyarraytoC(Ga)
Gabeg_c = copyarraytoC(Gabeg)
Gaend_c = copyarraytoC(Gaend)
ha_c = copyarraytoC(ha)
  
habeg_c = copyarraytoC(habeg)
haend_c = copyarraytoC(haend)
  
uabeg_c = copyarraytoC(uabeg)
uaend_c = copyarraytoC(uaend)
   
hmbeg_c = copyarraytoC(hmbeg)
hmend_c = copyarraytoC(hmend)
    
umbeg_c = copyarraytoC(umbeg)
umend_c = copyarraytoC(umend)
    
u_c= mallocPy(n)
G_c = mallocPy(n)
h_c = mallocPy(n)
    
for i in range(1,len(t)):
    if(i ==1 or i %gap == 0):
        ca2midpt(ha_c,dx,n,h_c)
        ca2midpt(Ga_c,dx,n,G_c)
        ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
        u = copyarrayfromC(u_c,n)
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        s = wdir + "out" + str(i) +  ".txt"
        with open(s,'a') as file2:
             writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
             writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'G' , 'u(m/s)'])        
                           
             for j in range(n):
                 writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(G[j]) , str(u[j])])    
    evolvewrap(Ga_c,ha_c,Gabeg_c,Gaend_c,habeg_c,haend_c,hmbeg_c,hmend_c,uabeg_c,uaend_c,umbeg_c,umend_c,nfcBC,nGsBC,g,dx,dt,n,cnBC,niBC)
    print (t[i])
    
        
ca2midpt(ha_c,dx,n,h_c)
ca2midpt(Ga_c,dx,n,G_c)
ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
u = copyarrayfromC(u_c,n)
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)
s = wdir + "outlast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'G' , 'u(m/s)'])        
                   
     for j in range(n):
         writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(G[j]) , str(u[j])])    
deallocPy(u_c)   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(ha_c)
deallocPy(Ga_c)
deallocPy(habeg_c)
deallocPy(Gabeg_c)
deallocPy(haend_c)
deallocPy(Gaend_c)
deallocPy(uabeg_c)
deallocPy(uaend_c)
deallocPy(umbeg_c)
deallocPy(umend_c)
deallocPy(hmbeg_c)
deallocPy(hmend_c)       
"""


##### Chanson
tl = 300
b = 0.61
h0 = 0.09
h1 = 0.1
g = 9.81

h0 = 0.192
u0 = 0.199
h1 = 0.22
u1 = 0.0  

dx = 0.01115
Cr = 0.5
l = Cr / sqrt(g*h0)
dt = l*dx
startx = -tl
endx = tl + dx
startt = 0
endt = 5 + dt  

wdir = "../../../data/raw/Canstest/o3/"

if not os.path.exists(wdir):
    os.makedirs(wdir)


#number of boundary conditions (one side)
nfcBC = 4 #for flux calculation
nGsBC = 2 #for solving G from u,h
niBC = nGsBC + nfcBC #total
        
    
gap = int(0.5/dt)
        
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)

#u1 = -0.199

#etm = 1.0
#um = 0.0
#etp = 0.22/0.192
#up = 
  
hm,um =  experiment2(x,h0,u0,h1,u1,dx)
plot(x,hm)
xlim([-2,2])
               
umbegi = zeros(niBC)
umendi = zeros(niBC)
hmbegi = ones(niBC)
hmendi = ones(niBC)    
        
for i in range(niBC):
    umbegi[i] = um[0]
    umendi[i] = um[-1]
    hmbegi[i] = hm[0]
    hmendi[i] = hm[-1]
        
umbeg = umbegi
umend = umendi
hmbeg = hmbegi
hmend = hmendi
        
#calculate G midpoint
cnBC = niBC - nGsBC
        
umbc = concatenate([umbeg[-cnBC:],um,umend[0:cnBC]]) 
hmbc = concatenate([hmbeg[-cnBC:],hm,hmend[0:cnBC]])       
Gmbc = solveGfromuh(umbc,hmbc,hmbeg[0:-cnBC],hmend[-cnBC:],umbeg[0:-cnBC],umend[-cnBC:],dx)  
        
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
    
Ga_c = copyarraytoC(Ga)
Gabeg_c = copyarraytoC(Gabeg)
Gaend_c = copyarraytoC(Gaend)
ha_c = copyarraytoC(ha)
  
habeg_c = copyarraytoC(habeg)
haend_c = copyarraytoC(haend)
  
uabeg_c = copyarraytoC(uabeg)
uaend_c = copyarraytoC(uaend)
   
hmbeg_c = copyarraytoC(hmbeg)
hmend_c = copyarraytoC(hmend)
    
umbeg_c = copyarraytoC(umbeg)
umend_c = copyarraytoC(umend)
    
u_c= mallocPy(n)
G_c = mallocPy(n)
h_c = mallocPy(n)
    
for i in range(1,len(t)):
    if(i ==1 or i %gap == 0):
        ca2midpt(ha_c,dx,n,h_c)
        ca2midpt(Ga_c,dx,n,G_c)
        ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
        u = copyarrayfromC(u_c,n)
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        s = wdir + "out" + str(i) +  ".txt"
        with open(s,'a') as file2:
             writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
             writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'G' , 'u(m/s)'])        
                           
             for j in range(n):
                 writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(G[j]) , str(u[j])])    
    evolvewrap(Ga_c,ha_c,Gabeg_c,Gaend_c,habeg_c,haend_c,hmbeg_c,hmend_c,uabeg_c,uaend_c,umbeg_c,umend_c,nfcBC,nGsBC,g,dx,dt,n,cnBC,niBC)
    print (t[i])
    
        
ca2midpt(ha_c,dx,n,h_c)
ca2midpt(Ga_c,dx,n,G_c)
ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
u = copyarrayfromC(u_c,n)
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)
s = wdir + "outlast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'G' , 'u(m/s)'])        
                   
     for j in range(n):
         writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(G[j]) , str(u[j])])    
deallocPy(u_c)   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(ha_c)
deallocPy(Ga_c)
deallocPy(habeg_c)
deallocPy(Gabeg_c)
deallocPy(haend_c)
deallocPy(Gaend_c)
deallocPy(uabeg_c)
deallocPy(uaend_c)
deallocPy(umbeg_c)
deallocPy(umend_c)
deallocPy(hmbeg_c)
deallocPy(hmend_c) 