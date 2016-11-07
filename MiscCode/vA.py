# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 07:07:18 2016

@author: jordan
"""
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

def quadsolver(a,b,c):
    v1 = (-b + sqrt(b**2 - 4*a*c))/(2.0*a)
    v2 = (-b - sqrt(b**2 - 4*a*c))/(2.0*a)
    return v1,v2

def v(dx,g,h0,Fgu,Fgn,Fnn,Fnu,G,D):
    from scipy import sqrt
    im = sqrt(-1)
    idx = 1.0 / dx
    a = -G
    
    b=im*idx*D*(Fgu + G*Fnn)
    
    c = idx*idx*D*D*(Fnn*Fgu - Fnu*Fgn)
    
    v1,v2 = quadsolver(a,b,c )
    
    return v1, v2
    
def D(k,dx):
    im = sqrt(-1)
    return 1 - e**(-im*k*dx)
    
def G2(h0,k,dx):
    from scipy import cos
    i3 = 1.0 / 3.0
    idx = 1.0 / dx
    G1 = h0  - (h0*h0*h0*i3)*(idx*idx)*(2*(cos(k*dx) - 1))
    
    return G1

def Fgn(g,h,Rm,Rp):
    return 0.5*g*h*(Rm-Rp)
    
def Fgu(G,Rm,Rp):
    return -G*(Rp-Rm)
    
def Fnu(h,Rm,Rp):
    return 0.5*h*(Rm-Rp)
    
def Fnn(Rm,Rp):
    return -(Rp-Rm)

def Rm(k,dx):
    im = sqrt(-1)
    Rm1 = 1 + (im*sin(k*dx))/(2.0)
    return Rm1
    
def Rp(k,dx):
    im = sqrt(-1)
    Rp1 = e**(im*k*dx)*(1 - (im*sin(k*dx))/(2.0))
    return Rp1
    
def o2(dx,k,h0,g):
    Rp1 = Rp(k,dx)
    Rm1 = Rm(k,dx)
    print("Edges")
    print(Rp1,Rm1)
    
    G1 = G2(h0,k,dx)
    print("G")
    print(G1)    
    
    Fgn1 = Fgn(g,h0,Rm1,Rp1)
    Fgu1 = Fgu(G1,Rm1,Rp1)
    Fnn1 = Fnn(Rm1,Rp1)
    Fnu1 = Fnu(h0,Rm1,Rp1)
    
    print("F")
    print(Fgn1, Fgu1)
    print(Fnn1, Fnu1) 
    
    D1 = D(k,dx)
    
    print("D")
    print(D1)
    
    return v(dx,g,h0,Fgu1,Fgn1,Fnn1,Fnu1,G1,D1)
    
def wactual(k,g,h0):
    from scipy import sqrt
    w1 = k*sqrt(g*h0)*sqrt(3.0 / (k*k*h0*h0 + 3))
    return w1,-w1

"""
h0 = 1
k = 0.5
g = 9.81

dxs = linspace(0,10,num=100)  
n = len(dxs)
w1s = zeros(n)
w2s = zeros(n) 
o21s = zeros(n)
o22s = zeros(n)
err1 = zeros(n)
err2 = zeros(n)
for i in range(n):
    w1,w2 = wactual(k,g,h0)
    w1s[i] = w1
    w2s[i] = w2
    
    o21,o22 = o2(dxs[i],k,h0,g)
    o21s[i] = o22
    o22s[i] = o21
    
    err1[i] = abs(w1.real - o22.real) / abs(w1.real)
    
    err2[i] = abs(w2.real - o21.real) / abs(w2.real)
"""
    