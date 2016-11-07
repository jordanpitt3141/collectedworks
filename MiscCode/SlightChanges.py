# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 07:07:18 2016

@author: jordan
"""
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

def quadsolver(a,b,c):
    from scipy import sqrt
    v1 = (-b + sqrt(b**2 - 4*a*c))/(2.0*a)
    v2 = (-b - sqrt(b**2 - 4*a*c))/(2.0*a)
    return v1,v2
 
def v(dx,g,h0,Fgu,Fgn,Fnn,Fnu,G,D,M):
    from scipy import sqrt
    im = sqrt(-1)
    idx = 1.0 / dx
    a = -G*M**2
    b=im*idx*D*M*( Fgu + G*Fnn )
    
    c = idx*idx*D*D*(Fnn*Fgu - Fnu*Fgn)
    
    v1,v2 = quadsolver(a,b,c )
    
    return v1, v2
   
def G2(h0,k,dx):
    from scipy import cos
    i3 = 1.0 / 3.0
    idx = 1.0 / dx
    G1 = h0  - (h0*h0*h0*i3)*(idx*idx)*(2*(cos(k*dx) - 1))
    
    return G1
    
def G4(h0,k,dx):
    from scipy import cos
    i3 = 1.0 / 3.0
    i12 = 1.0 / 12.0
    idx = 1.0 / dx
    G1 = h0  - (h0*h0*h0*i3)*(idx*idx)*i12*(32*cos(k*dx) - 2*cos(2*k*dx) - 30)
    
    return G1



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
    
def Ruo2(k,dx):
    im = sqrt(-1)
    return 0.5*(1 + e**(im*k*dx))
    
def Ruo3(k,dx):
    im = sqrt(-1)
    return (1.0/48.0)*(-3*e**(2*im*k*dx) + 27*e**(im*k*dx) + 27 -3*e**(-im*k*dx) )
    
def D(k,dx):
    im = sqrt(-1)
    return 1 - e**(-im*k*dx)
    
def Fgn(g,h,Rm,Rp):
    return 0.5*g*h*(Rm+Rp)
    
def Fgu(g,h,G,Rm,Rp):
    return -0.5*sqrt(g*h)*G*(Rp-Rm)
    
def Fnu(h,Ru):
    return h*Ru
    
def Fnn(g,h,Rm,Rp):
    return -0.5*sqrt(g*h)*(Rp-Rm)
    
def M3(k,dx):
    
    return 24.0/(26 - 2*cos(k*dx))
    #return 1


def o3(dx,k,g,h0):
    
    
    M =   M3(k,dx)  
    
    Rp1 = Rpo3(k,M,dx)
    Rm1 = Rmo3(k,M,dx)
    Ru = Ruo3(k,dx)
    #print("Edges")
    #print(Rp1,Rm1)
    
    G1 = G4(h0,k,dx)
    #print("G")
    #print(G1)    
    
    Fgn1 = Fgn(g,h0,Rm1,Rp1)
    Fgu1 = Fgu(g,h0,G1,Rm1,Rp1)
    Fnn1 = Fnn(g,h0,Rm1,Rp1)
    Fnu1 = Fnu(h0,Ru)
    
    #print("F")
    #print(Fgn1, Fgu1)
    #print(Fnn1, Fnu1) 
    
    D1 = D(k,dx)
    
    #print("D")
    #print(D1)
    
    return v(dx,g,h0,Fgu1,Fgn1,Fnn1,Fnu1,G1,D1,M)

def o2(dx,k,g,h0):
    Rp1 = Rpo2(k,dx)
    Rm1 = Rmo2(k,dx)
    Ru = Ruo2(k,dx)
    #print("Edges")
    #print(Rp1,Rm1)
    
    G1 = G2(h0,k,dx)
    #print("G")
    #print(G1)    
    
    Fgn1 = Fgn(g,h0,Rm1,Rp1)
    Fgu1 = Fgu(g,h0,G1,Rm1,Rp1)
    Fnn1 = Fnn(g,h0,Rm1,Rp1)
    Fnu1 = Fnu(h0,Ru)
    
    #print("F")
    #print(Fgn1, Fgu1)
    #print(Fnn1, Fnu1) 
    
    D1 = D(k,dx)
    
    #print("D")
    #print(D1)
    
    return v(dx,g,h0,Fgu1,Fgn1,Fnn1,Fnu1,G1,D1,1)
    
def o1(dx,k,g,h0):
    Rp1 = Rpo1(k,dx)
    Rm1 = Rmo1(k,dx)
    Ru = Ruo2(k,dx)
    #print("Edges")
    #print(Rp1,Rm1)
    
    G1 = G2(h0,k,dx)
    #print("G")
    #print(G1)    
    
    Fgn1 = Fgn(g,h0,Rm1,Rp1)
    Fgu1 = Fgu(g,h0,G1,Rm1,Rp1)
    Fnn1 = Fnn(g,h0,Rm1,Rp1)
    Fnu1 = Fnu(h0,Ru)
    
    #print("F")
    #print(Fgn1, Fgu1)
    #print(Fnn1, Fnu1) 
    
    D1 = D(k,dx)
    
    #print("D")
    #print(D1)
    
    return v(dx,g,h0,Fgu1,Fgn1,Fnn1,Fnu1,G1,D1,1)
    
def wactual(k,g,h0):
    from scipy import sqrt
    w1 = k*sqrt(g*h0)*sqrt(3.0 / (k*k*h0*h0 + 3))
    return w1,-w1

from scipy import pi

h0 = 1.0
k = 2.5
g = 9.81

Ns = linspace(0.0001,0.2,num=100)  
dxs = linspace(0.0001,0.2,num=100)  
n = len(dxs)
w1s = zeros(n)
w2s = zeros(n) 
o21s = zeros(n,dtype=complex)
o22s = zeros(n,dtype=complex)
o11s = zeros(n,dtype=complex)
o12s = zeros(n,dtype=complex)
o31s = zeros(n,dtype=complex)
o32s = zeros(n,dtype=complex)
erro21 = zeros(n)
erro22 = zeros(n)
erro11 = zeros(n)
erro12 = zeros(n)
erro31 = zeros(n)
erro32 = zeros(n)
for i in range(n):
    dx = ((2*pi)/k)*Ns[i]
    w1,w2 = wactual(k,g,h0)
    w1s[i] = w1
    w2s[i] = w2
    
    o31,o32 = o3(dx,k,g,h0)
    o21,o22 = o2(dx,k,g,h0)
    o11,o12 = o1(dx,k,g,h0)
    o21s[i] = o22
    o22s[i] = o21
    o11s[i] = o12
    o12s[i] = o11
    o31s[i] = o32
    o32s[i] = o31
    
    erro31[i] = abs(w1 - o32) / abs(w1)
    
    erro32[i] = abs(w2 - o31) / abs(w2)
    
    erro21[i] = abs(w1 - o22) / abs(w1)
    
    erro22[i] = abs(w2 - o21) / abs(w2)
    
    erro11[i] = abs(w1 - o12) / abs(w1)
    
    erro12[i] = abs(w2 - o11) / abs(w2)

    