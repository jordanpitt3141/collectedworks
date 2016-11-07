# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:09:21 2015

@author: jordan
"""
from Serre2 import *
from scipy import *
import csv
import os
from numpy.linalg import norm  
from matplotlib.pyplot import plot,ylim
from scipy.special import ellipj,ellipk,ellipe
    
def copyarraytoC(a):
    n = len(a)
    b = mallocPy(n)
    for i in range(n):
        writetomem(b,i,a[i])
    return b
    
def copywritearraytoC(a,b):
    n = len(a)
    for i in range(n):
        writetomem(b,i,a[i])
    
def copyarrayfromC(a,n):
    b = [0]*n
    for i in range(n):
        b[i] = readfrommem(a,i)
        
    return b
    
def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 

        
#gives exact up to linears, so is second order accurate huzzah    
def getGfromupy(h,u,bed,u0,u1,h0,h1,b0,b1,dx):
    idx = 1.0 / dx
    ithree = 1.0 / 3.0
        
    n = len(h)

    G = zeros(n)
        
    for i in range(1,n-1):
        th = h[i]
        thx = 0.5*idx*(h[i+1] - h[i-1])
        tbx = 0.5*idx*(bed[i+1] - bed[i-1])
        tbxx = idx*idx*(bed[i+1] -2*bed[i] + bed[i-1])
        
        D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
        
        ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
        bi = D + 2.0*ithree*idx*idx*th*th*th
        ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
        
        G[i] = ai*u[i-1] + bi*u[i] + ci*u[i+1]
        
    #boundary    
    #i=0
    i=0
    th = h[i]
    thx = 0.5*idx*(h[i+1] - h0)
    tbx = 0.5*idx*(bed[i+1] - b0)
    tbxx = idx*idx*(bed[i+1] -2*bed[i] + b0)
            
    D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
            
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
 
    G[i] = ai*u0 + bi*u[i] + ci*u[i+1]
    
    #i = n-1
    i = n-1

    th = h[i]
    thx = 0.5*idx*(h1 - h[i-1])
    tbx = 0.5*idx*(b1 - bed[i-1])
    tbxx = idx*idx*(b1 -2*bed[i] + bed[i-1])
        
    D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
        
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx

    G[i] = ai*u[i-1] + bi*u[i] + ci*u1
            
    return G  
    

def dambreak(x,hf,hc,hl,bot,dx):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    bed = bot*ones(n)
    
    for i in range(n):
        if (x[i] < hc):
            h[i] = hf
        else:
            h[i] = hl
    
    G = getGfromupy(h,u,bed,0.0,0.0,h[0],h[-1],0.0,0.0,dx)
    return h,G,bed

       

def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a

def soliton (x,t,g,a0,a1):
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)
  
def solitoninit(n,a0,a1,g,x,t0,bot,dx):
    h = zeros(n)
    bx = zeros(n)
    u = zeros(n)
    c = sqrt(g*(a0 + a1))
    for i in range(n):
        bx[i] = bot
        h[i] = soliton(x[i],t0,g,a0,a1)
        u[i] =  c* ((h[i] - a0) / h[i])
         
    G = getGfromupy(h,u,bx,0.0,0.0,a0,a0,0.0,0.0,dx)
    
    return h,u,G,bx 
    
def experiment1(x,b,h0,h1,dx):
    n = len(x)
    u = zeros(n)
    bx = zeros(n)
    h = ones(n)*h1
    for i in range(n):
        if (x[i] <0 and x[i] > -2*b):
            h[i] = h0
    G = getGfromupy(h,u,bx,0.0,0.0,h1,h1,0.0,0.0,dx)

    return h,G,bx 

def dambreaksmooth(x,x0,base,eta0,diffuse,bot,dx):
    from numpy import tanh
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    bed = bot*ones(n)
    
    for i in range(n):
        h[i] = base + 0.5*eta0*(1 + tanh(diffuse*(x0 - abs(x[i]))))
    
    G = getGfromupy(h,u,bed,0.0,0.0,h[0],h[-1],0.0,0.0,dx)
    return h,G,bed

def powerfunction(r,n):
    if ( r >= 0 and r<= 1):
        return (1-r)**n
    else:
        return 0
    
    
    
def flowoverbump(x,stage,center,width,height,vel,l):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    bed = zeros(n)
    
    for i in range(n): 
        r = abs(x[i] - center) / width
        bed[i] = height*(powerfunction(r,l + 2)*((l*l + 4*l + 3)*r*r*(1.0/3) + (l + 2)*r  + 1))
        h[i] = stage - bed[i]
        u[i] = vel
        
    G = getGfromupy(h,u,bed,u[0],u[-1],h[0],h[-1],bed[0],bed[-1],dx)    
    
    
    return h,G,bed
    
def soloverbump(x,a0,a1,solbeg,solend,t0,g,stage,center,width,height,vel,l):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    bed = zeros(n)
    
    c = sqrt(g*(a0 + a1))
    for i in range(n):
        
        if (x[i] > solbeg and x[i] < solend):
            r = abs(x[i] - center) / width
            bed[i] = height*(powerfunction(r,l + 2)*((l*l + 4*l + 3)*r*r*(1.0/3) + (l + 2)*r  + 1))
            h[i] = soliton(x[i],t0,g,a0,a1) - bed[i]
            u[i] =  c* ((h[i] - a0) / h[i])
        else:
            r = abs(x[i] - center) / width
            bed[i] = height*(powerfunction(r,l + 2)*((l*l + 4*l + 3)*r*r*(1.0/3) + (l + 2)*r  + 1))
            h[i] = stage - bed[i]
            u[i] =  vel
            
        
    G = getGfromupy(h,u,bed,u[0],u[-1],h[0],h[-1],bed[0],bed[-1],dx)    
    
    
    return h,G,bed
    
def DBsmoothoverbump(x,hadd,dbstart, diffuse,g,stage,center,width,height,l):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    bed = zeros(n)

    for i in range(n):
        
        if (x[i] < dbstart):
            r = abs(x[i] - center) / width
            bed[i] = height*(powerfunction(r,l + 2)*((l*l + 4*l + 3)*r*r*(1.0/3) + (l + 2)*r  + 1))
            h[i] = stage - bed[i] + 0.5*hadd*(1 + tanh(diffuse*(dbstart - x[i])))
        else:
            r = abs(x[i] - center) / width
            bed[i] = height*(powerfunction(r,l + 2)*((l*l + 4*l + 3)*r*r*(1.0/3) + (l + 2)*r  + 1))
            h[i] = stage - bed[i]
            
        
    G = getGfromupy(h,u,bed,u[0],u[-1],h[0],h[-1],bed[0],bed[-1],dx)    
    
    
    return h,G,bed
    
def flatlake(x,dx,x0,stage):
    
    n = len(x)
    u = zeros(n)
    h = zeros(n)
    bed = zeros(n)
    for i in range(n):
        if(x[i] > x0):
            bed[i] = 5.0
            h[i] = stage - bed[i]
        else:
            h[i] = stage
    
    G = getGfromupy(h,u,bed,u[0],u[-1],h[0],h[-1],bed[0],bed[-1],dx)    
    
    
    return h,G,bed   

def dnsq(eta,m):
    
    sn,cn,dn,ph = ellipj(eta,m)
    
    return dn*dn
    

def cnoidalwaves(x,t,dx,a0,a1,g,k):
    
    n = len(x)
    u = zeros(n)
    h = zeros(n)
    bed = zeros(n)
    
    m = k*k
    
    h0 = a0 + a1*(float(ellipe(m)) / ellipk(m))    
    
    c = sqrt((g*a0*(a0 + a1)*(a0 + (1 - k*k)*a1))) / float(h0)
    
    Kc = sqrt(float(3*a1) / (4*a0*(a0 + a1)*(a0 + (1-k*k)*a1)))

    
    for i in range(n):
        h[i] = a0 + a1*dnsq(Kc*(x[i] - c*t),m)
        u[i] = c *(1 - float(h0)/h[i])
        
        
    h0i = a0 + a1*dnsq(Kc*(x[0] - dx - c*t),m)
    u0i = c *(1 - float(h0)/h0i) 
    
    h1i = a0 + a1*dnsq(Kc*(x[-1] + dx - c*t),m)
    u1i = c *(1 - float(h0)/h1i)
    
    G = getGfromupy(h,u,bed,u0i,u1i,h0i,h1i,bed[0],bed[-1],dx)   
    
    return h,u,G,bed
    
def eigenvectors(x,t,g,k,h0,u0,h1,u1):
    n = len(x)
    w = u0*k + k*sqrt(g*h0)*sqrt(3.0/ (h0*h0*k*k + 3))
    imu = sqrt(-1)
    bed = zeros(n)
    h = zeros(n)
    u = zeros(n)
    for i in range(n):
        h[i] = h0 + h1*exp(imu*(w*t + k*x[i]))
        u[i] = u0 + u1*exp(imu*(w*t + k*x[i]))
    
    cx = x[0] - dx     
    h0i = h0 + h1*exp(imu*(w*t + k*cx))
    u0i = u0 + u1*exp(imu*(w*t + k*cx))
    
    cx = x[-1] + dx 
    h1i = h0 + h1*exp(imu*(w*t + k*cx))
    u1i = h0 + h1*exp(imu*(w*t + k*cx))
    
    G = getGfromupy(h,u,bed,u0i,u1i,h0i,h1i,bed[0],bed[-1],dx) 
    
    return h,u,G,bed
    
def DingFlume(x,dx):
    n = len(x)
    bed = zeros(n)
    h = zeros(n)
    u = zeros(n)
    for i in range(n):

        if(0 <= x[i] <= 6):
            h[i] = 0.4            
        elif(6 < x[i] <= 12):
            bed[i] = 0.05*(x[i] - 6)
            h[i] = 0.4 - bed[i]
        elif(12 < x[i] <= 14):
            bed[i] = 0.3
            h[i] = 0.1
        elif(14 < x[i] <= 17):
            bed[i] = 0.3 - 0.1*(x[i] - 14)
            h[i] = 0.4 - bed[i]
        else:
            bed[i] = 0
            h[i] = 0.4      
    G = getGfromupy(h,u,bed,0,0,0.4,0.4,0,0,dx)
    return h,u,G,bed
    
def DingEdge(x,u,u0,u1,a0,T0,l0,t):
    il0 = 1.0 / l0
    iT0 = 1.0 / T0
    n = len(x)
    v = zeros(n)
    eta = zeros(n)
    omega = k*sqrt(9.81*hb) * sqrt(3.0 / (k*k*hb*hb*hb + 3))
    for i in range(n):
        ft = a0*cos(2*pi*(x[i]*il0 - t*iT0))
        st = pi*a0*a0*il0*cos(4*pi*(x[i]*il0 - t*iT0))
        tt = -pi*pi*0.5*a0*a0*a0*il0*il0*(cos(2*pi*(x[i]*il0 - t*iT0)) - cos(6*pi*(x[i]*il0 - t*iT0)))
        eta[i] = 0.4 + ft + st + tt
        ut = (omega/(k*hb))*(ft + st + tt)
        v[i] = ut
    
    cx = x[0] - dx    
    ft = a0*cos(2*pi*(cx*il0 - t*iT0))
    st = pi*a0*a0*il0*cos(4*pi*(cx*il0 - t*iT0))
    tt = -pi*pi*0.5*a0*a0*a0*il0*il0*(cos(2*pi*(cx*il0 - t*iT0)) - cos(6*pi*(cx*il0 - t*iT0)))
    h0 = 0.4 + ft + st + tt
    u0 = (omega/(k*hb))*(ft + st + tt)
    
    cx = x[-1] + dx    
    ft = a0*cos(2*pi*(cx*il0 - t*iT0))
    st = pi*a0*a0*il0*cos(4*pi*(cx*il0 - t*iT0))
    tt = -pi*pi*0.5*a0*a0*a0*il0*il0*(cos(2*pi*(cx*il0 - t*iT0)) - cos(6*pi*(cx*il0 - t*iT0)))
    h1 = 0.4 + ft + st + tt
    u1 = (omega/(k*hb))*(ft + st + tt)
    G = getGfromupy(eta,u,bed,u0,u1,h0,h1,0,0,dx)    
    return eta,G

def BejiEdge(x,t0):
    n = len(x)
    eta = zeros(n)
    bed = zeros(n)
    v = zeros(n)
    hb = 0.4
    omega = k*sqrt(9.81*hb) * sqrt(3.0 / (k*k*hb*hb*hb + 3))
    for i in range(n):
        ft = a0*sin(k*x[i] - omega*t )
        ut = (omega/(k*hb))*ft
        eta[i] = hb + ft
        v[i] = ut
    
    cx = x[0] - dx    
    ft = a0*sin(k*cx - omega*t )
    ut = (omega/(k*hb))*ft
    h0 = hb + ft
    v0 = ut
    
    cx = x[-1] + dx    
    ft = a0*sin(k*cx - omega*t )
    ut = (omega/(k*hb))*ft
    h1 = hb + ft
    v1 = ut
    G = getGfromupy(eta,v,bed,v0,v1,h0,h1,0,0,dx)  
    return eta,v,G
    
def CosineEdge(x,a0,k,t):
    n = len(x)
    eta = zeros(n)
    bed = zeros(n)
    v = zeros(n)
    hb = 0.4
    omega = k*sqrt(9.81*hb) * sqrt(3.0 / (k*k*hb*hb*hb + 3))
    for i in range(n):
        ft = a0*sin(k*x[i] - omega*t )
        ut = (omega/(k*hb))*ft
        eta[i] = hb + ft
        v[i] = ut
    
    cx = x[0] - dx    
    ft = a0*sin(k*cx - omega*t )
    ut = (omega/(k*hb))*ft
    h0 = hb + ft
    v0 = ut
    
    cx = x[-1] + dx    
    ft = a0*sin(k*cx - omega*t )
    ut = (omega/(k*hb))*ft
    h1 = hb + ft
    v1 = ut
    G = getGfromupy(eta,v,bed,v0,v1,h0,h1,0,0,dx)  
    return eta,v,G
    
def SolitonEdge(x,g,d1,a0,c0,t0):
    n = len(x)
    eta = zeros(n)
    v = zeros(n)
    bed = zeros(n)
    cs = sqrt(g*(0.4 + d1))
    for i in range(n):
        eta[i] = soliton(x[i],t0 - 10,g,0.4,d1)
        v[i] =  cs* ((h[i] - 0.4) / h[i])
    
    cx = x[0] - dx    
    eta0 = soliton(cx,t0 - 10,g,0.4,d1)
    v0 =  cs* ((h[i] - 0.4) / h[i])
    
    cx = x[-1] + dx    
    eta1 = soliton(cx,t0 - 10,g,0.4,d1)
    v1 =  cs* ((h[i] - 0.4) / h[i])
    G = getGfromupy(eta,v,bed,v0,v1,eta0,eta1,0,0,dx)  
    return eta,v,G
    
    
   
"""
### Cnoidal wave with BC
wdatadir = "../../../data/raw/cnoidaltestfixlong/o2/"

if not os.path.exists(wdatadir):
    os.makedirs(wdatadir)
    
s = wdatadir + "savenorms.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(["dx",'theta','l1h', 'l1u'])  
for ij in range(13,14):
    a0 = 1.0
    a1 = 0.1
    k = 0.99
    
    ### WAVE LENGTH
        
    m = k*k
    Kc = sqrt(float(3*a1) / (4*a0*(a0 + a1)*(a0 + (1-m)*a1)))
    
    lamb = 2*ellipk(m) / Kc
    
    g = 9.81
    dx = 100.0 / 2**ij
    Cr = 0.5
    l = Cr / (sqrt(g*(a0 +a1) ))
    dt = l*dx
    theta = 2
    startx = 0.0
    endx = 9*lamb
    startt = 0.0
    endt = 100 + dt  
    
    wdir = wdatadir + str(ij) + "/"
    
    if not os.path.exists(wdir):
        os.makedirs(wdir)
    
    
    nBCn = 3
    nBC = 6
        
    xbc,t = makevar(startx - nBC*dx,endx + nBC*dx,dx,startt,endt,dt)
    
    x = xbc[nBC: -nBC]
    xbeg = xbc[:nBC]
    xend = xbc[-nBC:] 
    
    n = len(x)
    m = len(t)
    
    gap = int(10.0/dt)
    
    t0 = 0.0
        
    
    #initial conditions for time steps
    tij = 0.0
    hBC,uBC,GBC,bedBC = cnoidalwaves(xbc,tij,dx,a0,a1,g,k)
    h = hBC[nBC:-nBC]
    h0 = hBC[:nBC]
    h1 = hBC[-nBC:]
    u = uBC[nBC:-nBC]
    u0 = uBC[:nBC]
    u1 = uBC[-nBC:]
    G = GBC[nBC:-nBC]
    G0 = GBC[:nBC]
    G1 = GBC[-nBC:]
    bed = bedBC[nBC:-nBC]
    b0 = bedBC[:nBC]
    b1 = bedBC[-nBC:]
       
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    bed_c = copyarraytoC(bed)
    x_c = copyarraytoC(x)
    u_c = mallocPy(n)

    
    hBCh,uBCh,GBCh,bedBCh = cnoidalwaves(xbc,tij + dt,dx,a0,a1,g,k)
    h0h = hBCh[:nBC]
    h1h = hBCh[-nBC:]
    u0h = uBCh[:nBC]
    u1h = uBCh[-nBC:]
    G0h = GBCh[:nBC]
    G1h = GBCh[-nBC:]
    b0h = bedBCh[:nBC]
    b1h = bedBCh[-nBC:]
    
    un_c = mallocPy(n+2*nBCn)
    Gn_c = mallocPy(n+2*nBCn)
    hn_c = mallocPy(n+2*nBCn)
    
    h0_c = mallocPy(nBC)
    h1_c = mallocPy(nBC)
    u0_c = mallocPy(nBC)
    u1_c = mallocPy(nBC)
    G0_c = mallocPy(nBC)
    G1_c = mallocPy(nBC)
    b0_c = mallocPy(nBC)
    b1_c = mallocPy(nBC)
    
    h0h_c = mallocPy(nBC)
    h1h_c = mallocPy(nBC)
    u0h_c = mallocPy(nBC)
    u1h_c = mallocPy(nBC)
    G0h_c = mallocPy(nBC)
    G1h_c = mallocPy(nBC)
    b0h_c = mallocPy(nBC)
    b1h_c = mallocPy(nBC)
    
    hi,ui,Gi,bedi = cnoidalwaves(x,tij,dx,a0,a1,g,k)
    
    copywritearraytoC(h0,h0_c)
    copywritearraytoC(h1,h1_c)
    copywritearraytoC(u0,u0_c)
    copywritearraytoC(u1,u1_c)
    copywritearraytoC(G0,G0_c)
    copywritearraytoC(G1,G1_c)
    copywritearraytoC(b0,b0_c)
    copywritearraytoC(b1,b1_c)
    
    copywritearraytoC(h0h,h0h_c)
    copywritearraytoC(h1h,h1h_c)
    copywritearraytoC(u0h,u0h_c)
    copywritearraytoC(u1h,u1h_c)
    copywritearraytoC(G0h,G0h_c)
    copywritearraytoC(G1h,G1h_c)
    copywritearraytoC(b0h,b0h_c)
    copywritearraytoC(b1h,b1h_c) 
    
        
    for i in range(1,len(t)):  
         
        
        evolvewrapBC(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,G0_c,G1_c,h0h_c,h1h_c,u0h_c,u1h_c,G0h_c,G1h_c,b0_c,b1_c,g,dx,dt, n, nBC, nBCn,theta, hn_c,Gn_c,un_c)
        
        #evolvewrapperiodic(G_c,h_c,bed_c,g,dx,dt,n,nBCn,theta,hn_c, Gn_c,un_c);    
        print (t[i])
        
        copywritearraytoC(h0h,h0_c)
        copywritearraytoC(h1h,h1_c)
        copywritearraytoC(u0h,u0_c)
        copywritearraytoC(u1h,u1_c)
        copywritearraytoC(G0h,G0_c)
        copywritearraytoC(G1h,G1_c)
        copywritearraytoC(b0h,b0_c)
        copywritearraytoC(b1h,b1_c)
        
        hBCh,uBCh,GBCh,bedBCh = cnoidalwaves(xbc,t[i] + dt,dx,a0,a1,g,k)
        h0h = hBCh[:nBC]
        h1h = hBCh[-nBC:]
        u0h = uBCh[:nBC]
        u1h = uBCh[-nBC:]
        G0h = GBCh[:nBC]
        G1h = GBCh[-nBC:]
        b0h = bedBCh[:nBC]
        b1h = bedBCh[-nBC:]
        
        copywritearraytoC(h0h,h0h_c)
        copywritearraytoC(h1h,h1h_c)
        copywritearraytoC(u0h,u0h_c)
        copywritearraytoC(u1h,u1h_c)
        copywritearraytoC(G0h,G0h_c)
        copywritearraytoC(G1h,G1h_c)
        copywritearraytoC(b0h,b0h_c)
        copywritearraytoC(b1h,b1h_c) 
            
        tij = t[i]
    
    #getufromGperiodic(h_c,G_c,bed_c, dx ,n,u_c)
    
    #something weird with u at boundaires
    haBC,uaBC,GaBC,bedaBC = cnoidalwaves(xbc,tij,dx,a0,a1,g,k)
    ha = haBC[nBC:-nBC]
    h0 = haBC[:nBC]
    h1 = haBC[-nBC:]
    ua = uaBC[nBC:-nBC]
    u0 = uaBC[:nBC]
    u1 = uaBC[-nBC:]
    Ga = GaBC[nBC:-nBC]
    G0 = GaBC[:nBC]
    G1 = GaBC[-nBC:]
    beda = bedaBC[nBC:-nBC]
    b0 = bedaBC[:nBC]
    b1 = bedaBC[-nBC:]
       
    getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], 0.0, 0.0, dx ,n,u_c)
    u = copyarrayfromC(u_c,n)
    G = copyarrayfromC(G_c,n)
    h = copyarrayfromC(h_c,n)
    
    un = copyarrayfromC(un_c,n+2*nBCn)
    Gn = copyarrayfromC(Gn_c,n+2*nBCn)
    hn = copyarrayfromC(hn_c,n+2*nBCn)
    
    ha,ua,Ga,beda = cnoidalwaves(x,t[-1],dx,a0,a1,g,k)  
    
    s = wdir + "outlast.txt"
    with open(s,'a') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow(['dx' ,'dt','time' ,"cell midpoint", 'height(m)', 'G' , 'u(m/s)' ,'ha', 'u'])        
               
        for j in range(n):
            writefile2.writerow([str(dx),str(dt),str(t[-1]),str(x[j]) ,str(h[j]) , str(G[j]) , str(u[j]), str(ha[j]), str(ua[j])])     
    
    normhdiffi = norm(h - ha,ord=1) / norm(ha,ord=1)
    normudiffi = norm(u -ua,ord=1) / norm(ua,ord=1) 
    
    s = wdatadir + "savenorms.txt"
    with open(s,'a') as file1:
        writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile.writerow([str(dx),str(theta),str(normhdiffi), str(normudiffi)])  
        
    deallocPy(u_c)
    deallocPy(G_c)
    deallocPy(h_c)
    deallocPy(h0_c)
    deallocPy(h1_c)
    deallocPy(u0_c)
    deallocPy(u1_c)
    deallocPy(G0_c)
    deallocPy(G1_c)
    deallocPy(b0_c)
    deallocPy(b1_c)
    deallocPy(h0h_c)
    deallocPy(h1h_c)
    deallocPy(u0h_c)
    deallocPy(u1h_c)
    deallocPy(G0h_c)
    deallocPy(G1h_c)
    deallocPy(b0h_c)
    deallocPy(b1h_c) 
#### Cnoidal Waves
"""

"""
#incoming wave cosine test
wdatadir = "../../../data/raw/CosT/o2/"

if not os.path.exists(wdatadir):
    os.makedirs(wdatadir)
    
#a0 = 0.01
#T0 = 2.02
#l0 = 3.73
    
a0 = 0.01
c0 = 0.1
### WAVE LENGTH

g = 9.81
dx = 0.01
Cr = 0.5
l = Cr / (sqrt(g*(0.5) ))
dt = l*dx
theta = 2
startx = 0
endx = 1000
startt = 0.0
endt = 400 + dt  

wdir = wdatadir

if not os.path.exists(wdir):
    os.makedirs(wdir)


nBCn = 3
nBC = 6
    
xbc,t = makevar(startx - nBC*dx,endx + nBC*dx,dx,startt,endt,dt)

x = xbc[nBC: -nBC]
xbn = xbc[nBC - nBCn:-(nBC - nBCn)]
xbeg = xbc[:nBC]
xend = xbc[-nBC:] 

n = len(x)
m = len(t)

gap = int(10.0/dt)

t0 = 0.0
    

#initial conditions for time steps
tij = 0.0
h,u,G,bed = DingFlume(x,a0,c0,dx)

b0 = zeros(nBC)
b1 = zeros(nBC)
u0 = zeros(nBC)
u1 = zeros(nBC)

G1 = G[-1]*ones(nBC)

#first do just h at edges
h0,u0,G0 =  CosineEdge(xbeg,a0,c0,tij)

#h0,u0,G0 = SolitonEdge(xbeg,g,0.01,a0,c0,tij)
h0h,u0h,G0h =  CosineEdge(xbeg,a0,c0,tij+dt)
#h0h,u0h,G0h  = SolitonEdge(xbeg,g,0.01,a0,c0,tij +dt)

#hc,Gc = CosineEdge(x,u,u[0],u[-1],a0,c0,0)


h1 = 0.4*ones(nBC)

   
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
bed_c = copyarraytoC(bed)
x_c = copyarraytoC(x)
u_c = mallocPy(n)

un_c = mallocPy(n+2*nBCn)
Gn_c = mallocPy(n+2*nBCn)
hn_c = mallocPy(n+2*nBCn)

h1_c = mallocPy(nBC)
u1_c = mallocPy(nBC)
G0_c = mallocPy(nBC)
G1_c = mallocPy(nBC)
b0_c = mallocPy(nBC)
b1_c = mallocPy(nBC)

h0h_c = mallocPy(nBC)
h0_c = mallocPy(nBC)
u0h_c = mallocPy(nBC)
u0_c = mallocPy(nBC)
G0h_c = mallocPy(nBC)
G0_c = mallocPy(nBC)

copywritearraytoC(h0,h0_c)
copywritearraytoC(u0,u0_c)
copywritearraytoC(h1,h1_c)
copywritearraytoC(u1,u1_c)
copywritearraytoC(G1,G1_c)
copywritearraytoC(b0,b0_c)
copywritearraytoC(b1,b1_c)

copywritearraytoC(h0h,h0h_c)
copywritearraytoC(u0h,u0h_c)
copywritearraytoC(G0h,G0h_c)

    
for i in range(1,len(t)):  
    
    #getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], 0.0, 0.0, dx ,n,u_c)
    #u = copyarrayfromC(u_c,n)
    #copywritearraytoC(u[:nBC],u0_c)
    
    #evolvewrapBCwavetank(G_c,h_c,bed_c,h0_c,u0_c,G0_c,h1_c,u1_c,G1_c,h0h_c,u0h_c,G0h_c,h1_c,u1_c,G1_c,b0_c,b1_c,g,dx,dt, n, nBC, nBCn,theta, hn_c,Gn_c,un_c)
    evolvewrapBC(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,G0_c,G1_c,h0h_c,h1_c,u0h_c,u1_c,G0h_c,G1_c,b0_c,b1_c,g,dx,dt, n, nBC, nBCn,theta, hn_c,Gn_c,un_c)    
    tij = tij + dt
    copywritearraytoC(h0h,h0_c)
    copywritearraytoC(u0h,u0_c)
    copywritearraytoC(G0h,G0_c)
    h0h,u0h,G0h =  CosineEdge(xbeg,a0,c0,tij+dt)
    #h0h,u0h,G0h = SolitonEdge(xbeg,g,0.01,a0,c0,tij +dt)
    copywritearraytoC(h0h,h0h_c)
    copywritearraytoC(u0h,u0h_c)
    copywritearraytoC(G0h,G0h_c)
    
 
    #evolvewrapperiodic(G_c,h_c,bed_c,g,dx,dt,n,nBCn,theta,hn_c, Gn_c,un_c);    
    print (t[i])

    
    #print(h0h)
    

getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], 0.0, 0.0, dx ,n,u_c)
u = copyarrayfromC(u_c,n)
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)

un = copyarrayfromC(un_c,n+2*nBCn)
Gn = copyarrayfromC(Gn_c,n+2*nBCn)
hn = copyarrayfromC(hn_c,n+2*nBCn)

deallocPy(u_c)
deallocPy(G_c)
deallocPy(h_c)
deallocPy(h0_c)
deallocPy(h1_c)
deallocPy(u1_c)
deallocPy(G1_c)
deallocPy(b0_c)
deallocPy(b1_c) 
"""




### Dingemans results
wdatadir = "../../../data/raw/DingT/o2/"

if not os.path.exists(wdatadir):
    os.makedirs(wdatadir)
    
#a0 = 0.01
#T0 = 2.02
#l0 = 3.73
    
#a0 = 0.01
#T0 = 2.02
#l0 = 3.73
a0 = 0.01
k = 0.1

### WAVE LENGTH

g = 9.81
dx = 0.1
Cr = 0.5
l = Cr / (sqrt(g*(0.41) ))
dt = 0.039312
dx = dt / l

theta = 2
startx = -50
endx = 300
startt = 0
endt = 1*dt  

wdir = wdatadir

if not os.path.exists(wdir):
    os.makedirs(wdir)


nBCn = 3
nBC = 6
    
xbc,t = makevar(startx - nBC*dx,endx + nBC*dx,dx,startt,endt,dt)

x = xbc[nBC: -nBC]
xbeg = xbc[:nBC]
xend = xbc[-nBC:] 

n = len(x)
m = len(t)

gap = int(10.0/dt)

t0 = 0.0
    

#initial conditions for time steps
tij = 0.0
h,u,G,bed = DingFlume(x,dx)

b0 = zeros(nBC)
b1 = zeros(nBC)
u0 = zeros(nBC)
u1 = zeros(nBC)

G1 = G[-1]*ones(nBC)

#first do just h at edges
h0,u0,G0 =  CosineEdge(xbeg,a0,k,tij)#DingEdge(xbeg,u0,0,0,a0,T0,l0,tij)
h0h,u0h,G0h =  CosineEdge(xbeg,a0,k,tij + dt)#DingEdge(xbeg,u0,0,0,a0,T0,l0,tij+dt)


h1 = 0.4*ones(nBC)

   
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
bed_c = copyarraytoC(bed)
x_c = copyarraytoC(x)
u_c = mallocPy(n)

un_c = mallocPy(n+2*nBCn)
Gn_c = mallocPy(n+2*nBCn)
hn_c = mallocPy(n+2*nBCn)

h0_c = mallocPy(nBC)
h1_c = mallocPy(nBC)
u0_c = mallocPy(nBC)
u1_c = mallocPy(nBC)
G0_c = mallocPy(nBC)
G1_c = mallocPy(nBC)
b0_c = mallocPy(nBC)
b1_c = mallocPy(nBC)

h0h_c = mallocPy(nBC)
u0h_c = mallocPy(nBC)
G0h_c = mallocPy(nBC)

copywritearraytoC(h0,h0_c)
copywritearraytoC(h1,h1_c)
copywritearraytoC(u0,u0_c)
copywritearraytoC(u1,u1_c)
copywritearraytoC(G0,G0_c)
copywritearraytoC(G1,G1_c)
copywritearraytoC(b0,b0_c)
copywritearraytoC(b1,b1_c)

copywritearraytoC(h0h,h0h_c)
copywritearraytoC(u0h,u0h_c)
copywritearraytoC(G0h,G0h_c)

    
for i in range(1,len(t)):  
    
    getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], 0.0, 0.0, dx ,n,u_c)
    u = copyarrayfromC(u_c,n)
    copywritearraytoC(u[:nBC],u0_c)
    
    evolvewrapBC(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,G0_c,G1_c,h0h_c,h1_c,u0_c,u1_c,G0_c,G1_c,b0_c,b1_c,g,dx,dt, n, nBC, nBCn,theta, hn_c,Gn_c,un_c)
    
    #evolvewrapperiodic(G_c,h_c,bed_c,g,dx,dt,n,nBCn,theta,hn_c, Gn_c,un_c);    
    print (t[i])
    
    copywritearraytoC(h0h,h0_c)
    copywritearraytoC(G0h,G0_c)
    copywritearraytoC(u0h,u0_c)
    h0h,u0h,G0h =  CosineEdge(xbeg,a0,k,t[i] + dt)#DingEdge(xbeg,u[:nBC],u[0],u[nBC],a0,T0,l0,t[i]+dt)
    copywritearraytoC(h0h,h0h_c)
    copywritearraytoC(G0h,G0h_c)
    copywritearraytoC(u0h,u0h_c)
    
    #print(h0h)
    

getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], 0.0, 0.0, dx ,n,u_c)
u = copyarrayfromC(u_c,n)
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)

un = copyarrayfromC(un_c,n+2*nBCn)
Gn = copyarrayfromC(Gn_c,n+2*nBCn)
hn = copyarrayfromC(hn_c,n+2*nBCn)

  
deallocPy(u_c)
deallocPy(G_c)
deallocPy(h_c)
deallocPy(h0_c)
deallocPy(h1_c)
deallocPy(u0_c)
deallocPy(u1_c)
deallocPy(G0_c)
deallocPy(G1_c)
deallocPy(b0_c)
deallocPy(b1_c) 


    
####cnoidal waves periodic
"""
wdatadir = "../../../data/raw/cnoidaltest/o2/"

if not os.path.exists(wdatadir):
    os.makedirs(wdatadir)
    
s = wdatadir + "savenorms.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(["dx",'theta','l1h', 'l1u'])  
for ij in range(8,9):
    a0 = 1.0
    a1 = 0.1
    k = 0.99
    
    ### WAVE LENGTH
        
    m = k*k
    Kc = sqrt(float(3*a1) / (4*a0*(a0 + a1)*(a0 + (1-m)*a1)))
    
    lamb = 2*ellipk(m) / Kc
    
    g = 9.81
    dx = (lamb) / (3**ij)
    Cr = 0.5
    l = Cr / (sqrt(g*(a0 +a1) ))
    dt = l*dx
    theta = 2
    startx = dx
    endx = 9*lamb + 0.5*dx
    startt = 0.0
    endt = 1 + dt  
    
    wdir = wdatadir + str(ij) + "/"
    
    if not os.path.exists(wdir):
        os.makedirs(wdir)
    
    
    nBCn = 3
    nBC = 3
        
    xbc,t = makevar(startx - nBCn*dx,endx + nBCn*dx,dx,startt,endt,dt)
    
    x = xbc[nBCn: -nBCn]
    xbeg = x[:nBCn]
    xend = x[-nBCn:] 
    
    n = len(x)
    m = len(t)
    
    gap = int(10.0/dt)
    
    t0 = 0.0
        
    
    #initial conditions for time steps
    tij = 0.0
    hBC,uBC,GBC,bedBC = cnoidalwaves(xbc,tij,dx,a0,a1,g,k)
    h = hBC[nBCn:-nBCn]
    u = uBC[nBCn:-nBCn]
    G = GBC[nBCn:-nBCn]
    bed = bedBC[nBCn:-nBCn]
       
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    bed_c = copyarraytoC(bed)
    x_c = copyarraytoC(x)
    u_c = mallocPy(n)
    
    un_c = mallocPy(n+2*nBCn)
    Gn_c = mallocPy(n+2*nBCn)
    hn_c = mallocPy(n+2*nBCn)
    
    hi,ui,Gi,bedi = cnoidalwaves(x,tij,dx,a0,a1,g,k)
    
        
    for i in range(1,len(t)):  
        evolvewrapperiodic(G_c,h_c,bed_c,g,dx,dt,n,nBCn,theta,hn_c, Gn_c,un_c);    
        print (t[i])
            
        tij = t[i]
    
    getufromGperiodic(h_c,G_c,bed_c, dx ,n,u_c)
    #getufromG(h_c,G_c,bed_c,ubeg[-1],uend[0],hbeg[-1],hend[0], 0.0, 0.0, dx ,n,u_c)
    u = copyarrayfromC(u_c,n)
    G = copyarrayfromC(G_c,n)
    h = copyarrayfromC(h_c,n)
    
    un = copyarrayfromC(un_c,n+2*nBCn)
    Gn = copyarrayfromC(Gn_c,n+2*nBCn)
    hn = copyarrayfromC(hn_c,n+2*nBCn)
    
    ha,ua,Ga,beda = cnoidalwaves(x,t[-1],dx,a0,a1,g,k)  
    
    s = wdir + "outlast.txt"
    with open(s,'a') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow(['dx' ,'dt','time' ,"cell midpoint", 'height(m)', 'G' , 'u(m/s)' ,'ha', 'u'])        
               
        for j in range(n):
            writefile2.writerow([str(dx),str(dt),str(t[-1]),str(x[j]) ,str(h[j]) , str(G[j]) , str(u[j]), str(ha[j]), str(ua[j])])     
    
    normhdiffi = norm(h - ha,ord=1) / norm(ha,ord=1)
    normudiffi = norm(u -ua,ord=1) / norm(ua,ord=1) 
    
    s = wdatadir + "savenorms.txt"
    with open(s,'a') as file1:
        writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile.writerow([str(dx),str(theta),str(normhdiffi), str(normudiffi)])    
        
    deallocPy(u_c)   
    deallocPy(h_c)
    deallocPy(G_c)
"""

"""
#incoming wave cosine test
wdatadir = "../../../data/raw/CosT/o2/"

if not os.path.exists(wdatadir):
    os.makedirs(wdatadir)
    
#a0 = 0.01
#T0 = 2.02
#l0 = 3.73
    
hbase = 0.4
    
a0 = 0.01
c0 = 0.1
an = 4
### WAVE LENGTH

g = 9.81
dx = 0.1
Cr = 0.5
l = Cr / (sqrt(g*(0.5) ))
dt = l*dx
theta = 2
startx = dx
endx = 300
startt = 0.0
endt = 200 + dt  

xsl = endx  - 100
xsu = endx

wdir = wdatadir

if not os.path.exists(wdir):
    os.makedirs(wdir)


nBCn = 3
nBC = 6
    
xbc,t = makevar(startx - nBC*dx,endx + nBC*dx,dx,startt,endt,dt)

x = xbc[nBC: -nBC]
xbn = xbc[nBC - nBCn:-(nBC - nBCn)]
xbeg = xbc[:nBC]
xend = xbc[-nBC:] 

n = len(x)
m = len(t)

gap = int(10.0/dt)

t0 = 0.0
    

#initial conditions for time steps
tij = 0.0
h,u,G,bed = DingFlume(x,a0,c0,dx)

b0 = zeros(nBC)
u0 = zeros(nBC)
u1 = zeros(nBC)

G1 = G[-1]*ones(nBC)

#first do just h at edges
h0,u0,G0 =  CosineEdge(xbeg,a0,c0,tij)

#h0,u0,G0 = SolitonEdge(xbeg,g,0.01,a0,c0,tij)
h0h,u0h,G0h =  CosineEdge(xbeg,a0,c0,tij+dt)
#h0h,u0h,G0h  = SolitonEdge(xbeg,g,0.01,a0,c0,tij +dt)

#hc,Gc = CosineEdge(x,u,u[0],u[-1],a0,c0,0)


h1 = h[-1]*ones(nBC)
b1 = bed[-1]*ones(nBC)


   
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
bed_c = copyarraytoC(bed)
x_c = copyarraytoC(x)
u_c = mallocPy(n)

xbc_c = copyarraytoC(xbc)

un_c = mallocPy(n+2*nBCn)
Gn_c = mallocPy(n+2*nBCn)
hn_c = mallocPy(n+2*nBCn)

h1_c = mallocPy(nBC)
u1_c = mallocPy(nBC)
G0_c = mallocPy(nBC)
G1_c = mallocPy(nBC)
b0_c = mallocPy(nBC)
b1_c = mallocPy(nBC)

h0h_c = mallocPy(nBC)
h0_c = mallocPy(nBC)
u0h_c = mallocPy(nBC)
u0_c = mallocPy(nBC)
G0h_c = mallocPy(nBC)
G0_c = mallocPy(nBC)

copywritearraytoC(h0,h0_c)
copywritearraytoC(u0,u0_c)
copywritearraytoC(h1,h1_c)
copywritearraytoC(u1,u1_c)
copywritearraytoC(G1,G1_c)
copywritearraytoC(b0,b0_c)
copywritearraytoC(b1,b1_c)

copywritearraytoC(h0h,h0h_c)
copywritearraytoC(u0h,u0h_c)
copywritearraytoC(G0h,G0h_c)

    
for i in range(1,len(t)):  
    
    #getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], 0.0, 0.0, dx ,n,u_c)
    #u = copyarrayfromC(u_c,n)
    #copywritearraytoC(u[:nBC],u0_c)
    
    #evolvewrapBCwavetank(G_c,h_c,bed_c,h0_c,u0_c,G0_c,h1_c,u1_c,G1_c,h0h_c,u0h_c,G0h_c,h1_c,u1_c,G1_c,b0_c,b1_c,g,dx,dt, n, nBC, nBCn,theta, hn_c,Gn_c,un_c)
    #evolvewrapBC(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,G0_c,G1_c,h0h_c,h1_c,u0h_c,u1_c,G0h_c,G1_c,b0_c,b1_c,g,dx,dt, n, nBC, nBCn,theta, hn_c,Gn_c,un_c)    
    
    evolvewrapBCSponge(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,G0_c,G1_c,h0h_c,h1_c,u0h_c,u1_c,G0h_c,G1_c,b0_c,b1_c,g,dx,dt,n,nBC,nBCn,theta,hn_c, Gn_c,un_c, xbc_c,xsl,xsu,an,hbase,hbase+a0)

    tij = tij + dt
    copywritearraytoC(h0h,h0_c)
    copywritearraytoC(u0h,u0_c)
    copywritearraytoC(G0h,G0_c)
    h0h,u0h,G0h =  CosineEdge(xbeg,a0,c0,tij+dt)
    #h0h,u0h,G0h = SolitonEdge(xbeg,g,0.01,a0,c0,tij +dt)
    copywritearraytoC(h0h,h0h_c)
    copywritearraytoC(u0h,u0h_c)
    copywritearraytoC(G0h,G0h_c)
    
 
    #evolvewrapperiodic(G_c,h_c,bed_c,g,dx,dt,n,nBCn,theta,hn_c, Gn_c,un_c);    
    print (t[i])

    
    #print(h0h)
    

getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], 0.0, 0.0, dx ,n,u_c)
u = copyarrayfromC(u_c,n)
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)

un = copyarrayfromC(un_c,n+2*nBCn)
Gn = copyarrayfromC(Gn_c,n+2*nBCn)
hn = copyarrayfromC(hn_c,n+2*nBCn)

deallocPy(u_c)
deallocPy(G_c)
deallocPy(h_c)
deallocPy(h0_c)
deallocPy(h1_c)
deallocPy(u1_c)
deallocPy(G1_c)
deallocPy(b0_c)
deallocPy(b1_c) 

"""


#### soliton
"""

wdatadir = "../../../data/raw/soltest/o2/"

if not os.path.exists(wdatadir):
    os.makedirs(wdatadir)
    
s = wdatadir + "savenorms.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(["dx",'theta','l1h', 'l1u'])  

if not os.path.exists(wdatadir):
    os.makedirs(wdatadir)

for ij in range(7,8):
    a0 = 1.0
    a1 = 0.1
    k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
    xsl = 180
    xsu = 200.0
    an = 2
    
    ### WAVE LENGTH
        
    #m = k*k
    #Kc = sqrt(float(3*a1) / (4*a0*(a0 + a1)*(a0 + (1-m)*a1)))
    
    #lamb = 2*ellipk(m) / Kc
    
    g = 9.81
    dx = 0.1
    Cr = 0.5
    l = Cr / (sqrt(g*(a0 +a1) ))
    dt = l*dx
    theta = 2
    startx = -50
    endx = 200
    startt = 0.0
    endt = 200 + dt  
    
    wdir = wdatadir + str(ij) + "/"
    
    if not os.path.exists(wdir):
        os.makedirs(wdir)
    
    
    nBCn = 3
    nBC = 3
        
    xbc,t = makevar(startx - nBCn*dx,endx + nBCn*dx,dx,startt,endt,dt)
    
    x = xbc[nBCn: -nBCn]
    xbeg = x[:nBCn]
    xend = x[-nBCn:] 
    
    n = len(x)
    m = len(t)
    
    gap = int(10.0/dt)
    
    t0 = 0.0
        
    
    #initial conditions for time steps
    tij = 0.0
    hBC,uBC,GBC,bedBC = solitoninit(n + 2*nBCn,a0,a1,g,xbc,tij,0,dx)
    hbeg = hBC[:nBCn]
    hend = hBC[-nBCn:]
    h = hBC[nBCn:-nBCn]
    ubeg = uBC[:nBCn]
    uend = uBC[-nBCn:]
    u = uBC[nBCn:-nBCn]
    Gbeg = GBC[:nBCn]
    Gend = GBC[-nBCn:]
    G = GBC[nBCn:-nBCn]
    bedbeg = bedBC[:nBCn]
    bedend = bedBC[-nBCn:]
    bed = bedBC[nBCn:-nBCn]
       
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    bed_c = copyarraytoC(bed)
    x_c = copyarraytoC(x)
    xbc_c = copyarraytoC(xbc)
    u_c = mallocPy(n)
    
    un_c = mallocPy(n+2*nBCn)
    Gn_c = mallocPy(n+2*nBCn)
    hn_c = mallocPy(n+2*nBCn)
    
    b0_c = copyarraytoC(zeros(nBCn))
    b1_c = copyarraytoC(zeros(nBCn)) 
    #Gbcn = concatenate([Gbeg[:nBCn],G,Gend[-nBCn:]])
    
    hi,ui,Gi,bedi = solitoninit(n,a0,a1,g,x,tij,0,dx)
    
   
    h1 = a0*ones(nBC)
    h0 = a0*ones(nBC)
    u1 = 0.0*ones(nBC)
    u0 = 0.0*ones(nBC)   
    G1 = 0.0*ones(nBC)
    G0 = 0.0*ones(nBC)
    b1 = 0.0*ones(nBC)
    b0 = 0.0*ones(nBC)
    
    h0_c = mallocPy(nBC)
    h1_c = mallocPy(nBC)
    u0_c = mallocPy(nBC)
    u1_c = mallocPy(nBC)
    G0_c = mallocPy(nBC)
    G1_c = mallocPy(nBC)
    b0_c = mallocPy(nBC)
    b1_c = mallocPy(nBC)
    
    copywritearraytoC(h0,h0_c)
    copywritearraytoC(h1,h1_c)
    copywritearraytoC(u0,u0_c)
    copywritearraytoC(u1,u1_c)
    copywritearraytoC(G0,G0_c)
    copywritearraytoC(G1,G1_c)
    copywritearraytoC(b0,b0_c)
    copywritearraytoC(b1,b1_c)
    
    h0h_c = mallocPy(nBC)
    h1h_c = mallocPy(nBC)
    u0h_c = mallocPy(nBC)
    u1h_c = mallocPy(nBC)
    G0h_c = mallocPy(nBC)
    G1h_c = mallocPy(nBC)
    b0h_c = mallocPy(nBC)
    b1h_c = mallocPy(nBC)
    
    copywritearraytoC(h0,h0h_c)
    copywritearraytoC(h1,h1h_c)
    copywritearraytoC(u0,u0h_c)
    copywritearraytoC(u1,u1h_c)
    copywritearraytoC(G0,G0h_c)
    copywritearraytoC(G1,G1h_c)
    copywritearraytoC(b0,b0h_c)
    copywritearraytoC(b1,b1h_c)
        
    for i in range(1,len(t)): 
        
        
        #getufromG(h_c,G_c,bed_c,ubeg[-1],uend[0],hbeg[-1],hend[0], 0.0, 0.0, dx ,n,u_c)
        #u = copyarrayfromC(u_c,n)
        #G = copyarrayfromC(G_c,n)
        #h = copyarrayfromC(h_c,n) 
          
        
        
        
        #evolvewrap(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,G0_c,G1_c,h0h_c,h1h_c,u0h_c,u1h_c,G0h_c,G1h_c,b0_c,b1_c,g,dx,dt,n,nBC,nBCn,theta,hn_c, Gn_c,un_c)
        #evolvewrapperiodic(G_c,h_c,bed_c,g,dx,dt,n,nBCn,theta,hn_c, Gn_c,un_c);    
        evolvewrapBCSponge(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,G0_c,G1_c,h0h_c,h1h_c,u0h_c,u1h_c,G0h_c,G1h_c,b0_c,b1_c,g,dx,dt,n,nBC,nBCn,theta,hn_c, Gn_c,un_c, xbc_c,xsl,xsu,an,a0,a0+a1)
        print (t[i])
            
        tij = t[i]
    
        
        
        
    
    #getufromGperiodic(h_c,G_c,bed_c, dx ,n,u_c)
    getufromG(h_c,G_c,bed_c,ubeg[-1],uend[0],hbeg[-1],hend[0], 0.0, 0.0, dx ,n,u_c)
    u = copyarrayfromC(u_c,n)
    G = copyarrayfromC(G_c,n)
    h = copyarrayfromC(h_c,n)
    
    un = copyarrayfromC(un_c,n+2*nBCn)
    Gn = copyarrayfromC(Gn_c,n+2*nBCn)
    hn = copyarrayfromC(hn_c,n+2*nBCn)
    
    ha,ua,Ga,beda = solitoninit(n,a0,a1,g,x,t[-1],0,dx)
    
    haBC,uaBC,GaBC,bedBC = solitoninit(n + 2*nBCn,a0,a1,g,xbc,0,0,dx)
    
    s = wdir + "outlast.txt"
    with open(s,'a') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow(['dx' ,'dt','time' ,"cell midpoint", 'height(m)', 'G' , 'u(m/s)' ,'ha', 'u'])        
               
        for j in range(n):
            writefile2.writerow([str(dx),str(dt),str(t[-1]),str(x[j]) ,str(h[j]) , str(G[j]) , str(u[j]), str(ha[j]), str(ua[j])])     
    
    normhdiffi = norm(h - ha,ord=1) / norm(ha,ord=1)
    normudiffi = norm(u -ua,ord=1) / norm(ua,ord=1) 
    
    s = wdatadir + "savenorms.txt"
    with open(s,'a') as file1:
        writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile.writerow([str(dx),str(theta),str(normhdiffi), str(normudiffi)]) 
        
    deallocPy(u_c)   
    deallocPy(h_c)
    deallocPy(G_c)
"""




"""
###Soliton over bump scenario
### Energy Test ####################

wdatadir = "../../../data/raw/solbumpEnergucontfullfix0m3/o2/"

if not os.path.exists(wdatadir):
    os.makedirs(wdatadir)


s = wdatadir +  "savenorms.txt"
with open(s,'a') as file3:
    writefile3 = csv.writer(file3, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile3.writerow(['dx' ,'rel Eval'])



for k in range(6,19,3):
    stage = 1
    center = 50.0
    width = 25
    height = 0.5
    el = 4.0
    vel = 0
    
    
    a0 = 1.0
    a1 = 1.0
    
    g = 9.81
    dx = 100.0 / (2**k)
    Cr = 0.5
    l = Cr / (2.5 + sqrt(g*(a0 +a1) ))
    dt = l*dx
    theta = 1.0
    startx = -200
    endx = 250.0 + dx
    startt = 0.0
    endt = 50 + dt  
    
    wdir = wdatadir + str(k) + "/"
    
    if not os.path.exists(wdir):
        os.makedirs(wdir)
        
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
    m = len(t)
    
    gap = int(10.0/dt)
        
    #h,G,bed = flowoverbump(x,stage,center,width,height,vel,el)
    h,G,bed = soloverbump(x,a0,a1,-20,20,0.0,g,stage,center,width,height,vel,el)
        
    nBC = 3
    nBCs = 4
    b0 = bed[0]*ones(nBCs)
    b1 = bed[-1]*ones(nBCs)
    u0 = vel*ones(nBCs)
    u1 = vel*ones(nBCs)   
    h0 = h[0]*ones(nBCs)
    h1 = h[-1]*ones(nBCs)
        
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    bed_c = copyarraytoC(bed)
    h0_c  = copyarraytoC(h0)
    h1_c  = copyarraytoC(h1)
    u0_c  = copyarraytoC(u0)
    u1_c  = copyarraytoC(u1)
    b0_c  = copyarraytoC(b0)
    b1_c  = copyarraytoC(b1)
    u_c = mallocPy(n)
    
    xbeg = arange(startx - nBC*dx,startx,dx)
    xend = arange(endx + dx,endx + (nBC+1)*dx) 
    
    xbc =  concatenate([xbeg,x,xend]) 
    
    xbc_c = copyarraytoC(xbc)
    hbc_c = mallocPy(n + 2*nBC)
    ubc_c = mallocPy(n + 2*nBC)
    bedbc_c = mallocPy(n + 2*nBC)
    conc(b0_c , bed_c,b1_c,nBC,n ,nBC , bedbc_c)
    #HEvals = []
    GNEvals = []
    Etimes = []
    
    tBC = 2
    #initial conditions for time steps
        
    for i in range(1,len(t)): 
    
        if (i == 1 or i%gap == 0):
            ki = i
            getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
            u = copyarrayfromC(u_c,n)
            G = copyarrayfromC(G_c,n)
            h = copyarrayfromC(h_c,n)
            
            conc(h0_c , h_c,h1_c,nBC,n ,nBC , hbc_c)
            conc(u0_c , u_c,u1_c,nBC,n ,nBC , ubc_c)  
            #HEval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*nBC,nBC,dx)
            GNEval = GNall(xbc_c,hbc_c,ubc_c,bedbc_c,g,n + 2*nBC,nBC,dx)
            
            
            
            #HEvals.append(HEval)
            GNEvals.append(GNEval)
            Etimes.append(t[i-1])
            
            s = wdir +  "out" + str(i)+".txt"
            with open(s,'a') as file2:
                writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
                writefile2.writerow(['dx' ,'dt','time','Eval',"cell midpoint" ,'height(m)', 'G' , 'u(m/s)','bed' ])        
                           
                for j in range(n):
                    writefile2.writerow([str(dx),str(dt),str(t[i]), str(GNEval),str(x[j]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j])])
            
        evolvewrap(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,b0_c,b1_c,g,dx,dt,n,nBCs,theta)
        print (t[i])
            
    getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
    u = copyarrayfromC(u_c,n)
    G = copyarrayfromC(G_c,n)
    h = copyarrayfromC(h_c,n)
    
    conc(h0_c , h_c,h1_c,nBC,n ,nBC , hbc_c)
    conc(u0_c , u_c,u1_c,nBC,n ,nBC , ubc_c)        
    #HEval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*nBC,nBC,dx)
    GNEval = GNall(xbc_c,hbc_c,ubc_c,bedbc_c,g,n + 2*nBC,nBC,dx)
    
    
    
    #HEvals.append(HEval)
    GNEvals.append(GNEval)
    Etimes.append(t[-1])
    s = wdir +  "outlast.txt"
    with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
         writefile2.writerow(['dx' ,'dt','time','Eval',"cell midpoint" ,'height(m)', 'G' , 'u(m/s)','bed' ])        
                       
         for j in range(n):
             writefile2.writerow([str(dx),str(dt),str(t[i]),str(GNEval),str(x[j]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j])])
    
    finE = GNEvals[-1]
    tfinE = Etimes[-1]
    begE = GNEvals[0]
    tbegE = Etimes[0]
    relerr = abs(finE - begE) / abs(begE)
    
    print(str(tbegE) + " || " + str(begE))
    print(str(tfinE) + " || " + str(finE))
    s = wdatadir +  "savenorms.txt"
    with open(s,'a') as file3:
        writefile3 = csv.writer(file3, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile3.writerow([str(dx) ,str(relerr)])
    
    deallocPy(u_c)   
    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(h0_c)
    deallocPy(h1_c)
    deallocPy(u0_c)
    deallocPy(u1_c) 
"""