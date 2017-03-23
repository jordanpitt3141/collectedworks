# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 07:07:18 2016

@author: jordan

This code calculates the dispersion error for our FDVM
"""
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

#given a matrix [[a,b],[c,d]] calculates the eigenvalues of it
def eigenvaluecalc2by2(a,b,c,d):
    from scipy import sqrt
    T = a + d
    D = a*d - b*c
    
    l1 = T/2.0 + sqrt((T*T)/4.0 - D)
    l2 = T/2.0 - sqrt((T*T)/4.0 - D)
    
    return l1,l2
 
#the elliptic solver for u o2 (first,second) and o3(third)  
def G2(h0,k,dx):
    from scipy import cos
    i3 = 1.0 / 3.0
    idx = 1.0 / dx
    G1 = h0  - (h0*h0*h0*i3)*(idx*idx)*(2*(cos(k*dx) - 1))
    
    return G1
    
def GFEM2(h0,k,Rp,Rm,dx):
    from scipy import cos
    i3 = 1.0 / 3.0
    im = sqrt(-1)
    idx = 1.0 / dx
    Recon = 1.0 / (e**(-im*k*dx) + 4 +  e**(im*k*dx))
    #Recon = 1.0 / (e**(-im*k*dx)*Rp + 2*Rm + 2*Rp +  e**(im*k*dx)*Rm)
    Hmult = e**(-im*k*dx) + 4 + e**(im*k*dx)
    H3mult = 2*idx*idx*(-e**(-im*k*dx) + 2 - e**(im*k*dx))
    G1 = (h0*Hmult + h0*h0*h0*H3mult)*Recon
    
    return G1
    
def GFEM4(h0,k,Rp,Rm,dx):
    from scipy import cos
    i3 = 1.0 / 3.0
    im = sqrt(-1)
    idx = 1.0 / dx
    Recon = 2.0 / (e**(-im*k*dx)*Rp + 8 + Rm)
    G1 = Recon*( h0*(4 + cos(k*dx/2.0)) + ((40*h0*h0*h0*idx*idx)/3.0)*((1 - cos(k*dx/2.0))) )
    
    return G1
    
def G4(h0,k,dx):
    from scipy import cos
    i3 = 1.0 / 3.0
    i12 = 1.0 / 12.0
    idx = 1.0 / dx
    G1 = h0  - (h0*h0*h0*i3)*(idx*idx)*i12*(32*cos(k*dx) - 2*cos(2*k*dx) - 30)
    
    return G1

#the midpoint to cell average conversion contribution
def M3(k,dx):
    return 24.0/(26 - 2*cos(k*dx))

#the reconstruction contributions for o1 (first), o2(second) and o3(third order)
def Rpo3(k,M,dx):
    from scipy import sqrt,e,sin
    im = sqrt(-1)
    Rp1 = M*e**(im*k*dx)*(5 - e**(im*k*dx) + 2*e**(-im*k*dx)) /6.0
    return Rp1

def Rmo3(k,M,dx):
    im = sqrt(-1)
    Rm1 = M*(5 - e**(-im*k*dx) + 2*e**(im*k*dx)) /6.0
    return Rm1
    
def Rpo2(k,dx):
    from scipy import sqrt,e,sin
    im = sqrt(-1)
    Rp1 = e**(im*k*dx)*(1 - (im*sin(k*dx))/(2.0))
    return Rp1

def Rmo2(k,dx):
    im = sqrt(-1)
    Rm1 = 1 + (im*sin(k*dx))/(2.0)
    return Rm1
    
def Rpo1(k,dx):
    from scipy import sqrt,e,sin
    im = sqrt(-1)
    Rp1 = e**(im*k*dx)
    return Rp1

def Rmo1(k,dx):
    im = sqrt(-1)
    Rm1 = 1
    return Rm1
 
#the smooth reconstruction for u o2 (first,second) and o3(third)   
def Ruo2(k,dx):
    im = sqrt(-1)
    return 0.5*(1 + e**(im*k*dx))
    
def Ruo3(k,dx):
    im = sqrt(-1)
    return (1.0/48.0)*(-3*e**(2*im*k*dx) + 27*e**(im*k*dx) + 27 -3*e**(-im*k*dx) )

#The difference-shift operator (comes out of difference of fluxes)    
def D(k,dx):
    im = sqrt(-1)
    return 1 - e**(-im*k*dx)
 

#calculation of the flux contributions   
def Fun(g,h,Rm,Rp):
    return 0.5*g*h*(Rm+Rp)
    
def Fuu(g,h,G,Rm,Rp):
    return -0.5*sqrt(g*h)*G*(Rp-Rm)
    
def Fnu(h,Ru):
    return h*Ru
    
def Fnn(g,h,Rm,Rp):
    return -0.5*sqrt(g*h)*(Rp-Rm)
    



def o1(k,g,h,dx,dt):
    from scipy import log
    im = sqrt(-1)


    #calculate the elementary contribution factors
    M =   1 
    
    Rp1 = Rpo1(k,dx)
    Rm1 = Rmo1(k,dx)
    Ru = Ruo2(k,dx)
    
    G1 = G2(h,k,dx)  
    
    Fgn1 = Fun(g,h,Rm1,Rp1)
    Fgu1 = Fuu(g,h,G1,Rm1,Rp1)
    Fnn1 = Fnn(g,h,Rm1,Rp1)
    Fnu1 = Fnu(h,Ru)
    
    D1 = D(k,dx)
    
    #get all terms of the matrix F
    Fun1 = D1/(M*dx*G1)*Fgn1
    Fuu1 = D1/(M*dx*G1)*Fgu1
    Fnn1 = D1/(M*dx)*Fnn1
    Fnu1 = D1/(M*dx)*Fnu1
    
    
    #calculate eigenvalues
    l1,l2 = eigenvaluecalc2by2(Fnn1,Fnu1,Fun1,Fuu1)    
    
    #use formula from RK steps to get the numerical dispersion
    wn1 = 1.0/(im*dt)*log(1- dt*l1)
    wn2 = 1.0/(im*dt)*log(1- dt*l2)
    
    return wn1,wn2
    
def o2(k,g,h,dx,dt):
    from scipy import log
    im = sqrt(-1)

    #calculate the elementary contribution factors
    M =   1 
    
    Rp1 = Rpo2(k,dx)
    Rm1 = Rmo2(k,dx)
    Ru = Ruo2(k,dx)
    
    G1 = G2(h,k,dx) 
    
    Fgn1 = Fun(g,h,Rm1,Rp1)
    Fgu1 = Fuu(g,h,G1,Rm1,Rp1)
    Fnn1 = Fnn(g,h,Rm1,Rp1)
    Fnu1 = Fnu(h,Ru)
    
    D1 = D(k,dx)
    
    #get all terms of the matrix F
    Fun1 = D1/(M*dx*G1)*Fgn1
    Fuu1 = D1/(M*dx*G1)*Fgu1
    Fnn1 = D1/(M*dx)*Fnn1
    Fnu1 = D1/(M*dx)*Fnu1
    
    #calculate eigenvalues
    l1,l2 = eigenvaluecalc2by2(Fnn1,Fnu1,Fun1,Fuu1)
    
    #use formula from RK steps to get the numerical dispersion
    wn1 = 1.0/(im*dt)*log(1- dt*l1 + 0.5*dt*dt*l1*l1)
    wn2 = 1.0/(im*dt)*log(1- dt*l2+ 0.5*dt*dt*l2*l2)
    
    return wn1,wn2

def oFEM2(k,g,h,dx,dt):
    from scipy import log
    im = sqrt(-1)

    #calculate the elementary contribution factors
    M =   1 
    
    Rp1 = Rpo2(k,dx)
    Rm1 = Rmo2(k,dx)
    Ru = Ruo2(k,dx)
    
    G1 = GFEM2(h,k,Rp1,Rm1,dx) 
    
    Fgn1 = Fun(g,h,Rm1,Rp1)
    Fgu1 = Fuu(g,h,G1,Rm1,Rp1)
    Fnn1 = Fnn(g,h,Rm1,Rp1)
    Fnu1 = Fnu(h,Ru)
    
    D1 = D(k,dx)
    
    #get all terms of the matrix F
    Fun1 = D1/(M*dx*G1)*Fgn1
    Fuu1 = D1/(M*dx*G1)*Fgu1
    Fnn1 = D1/(M*dx)*Fnn1
    Fnu1 = D1/(M*dx)*Fnu1
    
    #calculate eigenvalues
    l1,l2 = eigenvaluecalc2by2(Fnn1,Fnu1,Fun1,Fuu1)
    
    #use formula from RK steps to get the numerical dispersion
    wn1 = 1.0/(im*dt)*log(1- dt*l1 + 0.5*dt*dt*l1*l1)
    wn2 = 1.0/(im*dt)*log(1- dt*l2+ 0.5*dt*dt*l2*l2)
    
    return wn1,wn2
    
def o3(k,g,h,dx,dt):
    from scipy import log
    im = sqrt(-1)

    #calculate the elementary contribution factors
    M =   M3(k,dx)  
    
    Rp1 = Rpo3(k,M,dx)
    Rm1 = Rmo3(k,M,dx)
    Ru = Ruo3(k,dx)
    
    G1 = G4(h,k,dx)
    
    Fgn1 = Fun(g,h,Rm1,Rp1)
    Fgu1 = Fuu(g,h,G1,Rm1,Rp1)
    Fnn1 = Fnn(g,h,Rm1,Rp1)
    Fnu1 = Fnu(h,Ru)
    
    D1 = D(k,dx)
    
    #get all terms of the matrix F
    Fun1 = D1/(M*dx*G1)*Fgn1
    Fuu1 = D1/(M*dx*G1)*Fgu1
    Fnn1 = D1/(M*dx)*Fnn1
    Fnu1 = D1/(M*dx)*Fnu1
    
    #calculate eigenvalues
    l1,l2 = eigenvaluecalc2by2(Fnn1,Fnu1,Fun1,Fuu1)
    
    #use formula from RK steps to get the numerical dispersion
    wn1 = 1.0/(im*dt)*log(1- dt*l1 + 0.5*dt*dt*l1*l1 - dt*dt*dt*l1*l1*l1/6.0)
    wn2 = 1.0/(im*dt)*log(1- dt*l2+ 0.5*dt*dt*l2*l2- dt*dt*dt*l2*l2*l2/6.0)
    
    return wn1,wn2


def oFEM3(k,g,h,dx,dt):
    from scipy import log
    im = sqrt(-1)

    #calculate the elementary contribution factors
    M =   M3(k,dx)  
    
    Rp1 = Rpo3(k,M,dx)
    Rm1 = Rmo3(k,M,dx)
    Ru = Ruo3(k,dx)
    
    G1 = GFEM4(h,k,Rp1,Rm1,dx)
    
    Fgn1 = Fun(g,h,Rm1,Rp1)
    Fgu1 = Fuu(g,h,G1,Rm1,Rp1)
    Fnn1 = Fnn(g,h,Rm1,Rp1)
    Fnu1 = Fnu(h,Ru)
    
    D1 = D(k,dx)
    
    #get all terms of the matrix F
    Fun1 = D1/(M*dx*G1)*Fgn1
    Fuu1 = D1/(M*dx*G1)*Fgu1
    Fnn1 = D1/(M*dx)*Fnn1
    Fnu1 = D1/(M*dx)*Fnu1
    
    #calculate eigenvalues
    l1,l2 = eigenvaluecalc2by2(Fnn1,Fnu1,Fun1,Fuu1)
    
    #use formula from RK steps to get the numerical dispersion
    wn1 = 1.0/(im*dt)*log(1- dt*l1 + 0.5*dt*dt*l1*l1 - dt*dt*dt*l1*l1*l1/6.0)
    wn2 = 1.0/(im*dt)*log(1- dt*l2+ 0.5*dt*dt*l2*l2- dt*dt*dt*l2*l2*l2/6.0)
    
    return wn1,wn2
    
def wactual(k,g,h0):
    from scipy import sqrt
    w1 = k*sqrt(g*h0)*sqrt(3.0 / (k*k*h0*h0 + 3))
    return w1,-w1

from scipy import pi


h = 1.0
k = 5.0
g = 9.81

w1,w2 = wactual(k,g,h)

#calculate the phase speeds by the relation omega/k
v1 = w1/k
v2 = w2/k

#this then measures the dispersion relation from the highest resolution soliton problem until 0.5
dxs = linspace(100.0/ (2**20),0.5,num=100)   
n = len(dxs)
w1s = zeros(n)
w2s = zeros(n) 
o21s = zeros(n,dtype=complex)
o22s = zeros(n,dtype=complex)
of21s = zeros(n,dtype=complex)
of22s = zeros(n,dtype=complex)
o11s = zeros(n,dtype=complex)
o12s = zeros(n,dtype=complex)
o31s = zeros(n,dtype=complex)
o32s = zeros(n,dtype=complex)
of31s = zeros(n,dtype=complex)
of32s = zeros(n,dtype=complex)
erro21 = zeros(n)
erro22 = zeros(n)
errof21 = zeros(n)
errof22 = zeros(n)
erro11 = zeros(n)
erro12 = zeros(n)
erro31 = zeros(n)
erro32 = zeros(n)
errof31 = zeros(n)
errof32 = zeros(n)

for i in range(n):
    dx = dxs[i]
    Cr = 0.5
    l = Cr/v1
    dt = l*dx
    w1,w2 = wactual(k,g,h)
    
    o11,o12 = o1(k,g,h,dx,dt)
    o21,o22 = o2(k,g,h,dx,dt)
    of21,of22 = oFEM2(k,g,h,dx,dt)
    o31,o32 = o3(k,g,h,dx,dt)
    of31,of32 = oFEM3(k,g,h,dx,dt)


    #Somehow this process does not guarantee that o11 and o12 are either always w1 or w2
    #so we check that o11 is positive or not, and then assign o11s and o12s accordingly so that
    # o11s[i] corresponds to the positive analytic dispersion w1
    if(o11.real > 0):
        o11s[i] = o11
        o12s[i] = o12
    else:
        o12s[i] = o11
        o11s[i] = o12
        
    if(o21.real > 0):
        o21s[i] = o21
        o22s[i] = o22
    else:
        o22s[i] = o21
        o21s[i] = o22
        
    if(of21.real > 0):
        of21s[i] = of21
        of22s[i] = of22
    else:
        of22s[i] = of21
        of21s[i] = of22

    if(o31.real > 0):
        o31s[i] = o31
        o32s[i] = o32
    else:
        o32s[i] = o31
        o31s[i] = o32
        
    if(of31.real > 0):
        of31s[i] = of31
        of32s[i] = of32
    else:
        of32s[i] = of31
        of31s[i] = of32
    
    erro11[i] = abs(w1 - o11s[i]) / abs(w1)    
    erro12[i] = abs(w2 - o12s[i]) / abs(w2)
    
    erro21[i] = abs(w1 - o21s[i]) / abs(w1)    
    erro22[i] = abs(w2 - o22s[i]) / abs(w2)
    
    errof21[i] = abs(w1 - of21s[i]) / abs(w1)    
    errof22[i] = abs(w2 - of22s[i]) / abs(w2)
    
    erro31[i] = abs(w1 - o31s[i]) / abs(w1)    
    erro32[i] = abs(w2 - o32s[i]) / abs(w2)
    
    errof31[i] = abs(w1 - of31s[i]) / abs(w1)    
    errof32[i] = abs(w2 - of32s[i]) / abs(w2)
    
    

    