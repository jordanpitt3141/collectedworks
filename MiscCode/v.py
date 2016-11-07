# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 07:07:18 2016

@author: jordan
"""
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

def v(dx,g,h0,F,G):
    from scipy import sqrt
    idx = 1.0 / dx
    
    v1 = (sqrt(-1)*G*F*idx + sqrt(-1)*sqrt(G*G*idx*idx*F*F + 4*G*F*g*h0*h0*idx)) / (2*G)
    v2 = (sqrt(-1)*G*F*idx  - sqrt(-1)*sqrt(G*G*idx*idx*F*F + 4*G*F*g*h0*h0*idx)) / (2*G)
    
    return v1, v2
    
def G2(h0,k,dx):
    from scipy import cos
    i3 = 1.0 / 3.0
    idx = 1.0 / dx
    G1 = h0  - (h0*h0*h0*i3)*(idx*idx)*(2*(cos(k*dx) - 1))
    
    return G1
    
def FG(k,dx,aip,aim,aim1p,aim1m,sp,sm):
    from scipy import exp,sqrt
    FT = (aip*sm - aim*sp) /(aip -aim )
    ST = exp(sqrt(-1)*k*dx)*(aim1p*sm - aim1m*sp) /(aim1p -aim1m )
    return FT + ST


def Spm2(k,dx):
    from scipy import sqrt,sin,exp
    im = sqrt(-1)
    Sp = exp(im*k*dx)*(1 - 0.5*(im*sin(im*k*dx)))
    Sm = (1 + 0.5*(im*sin(im*k*dx)))
    
    return Sp, Sm
    
def Spm1(k,dx):
    from scipy import sqrt,sin,exp
    im = sqrt(-1)
    Sp = exp(im*k*dx)
    Sm = 1
    
    return Sp, Sm
    
def aipm(g,k,h0):
    from scipy import sqrt
    return wactual(k,g,h0)
    
def o2(dx,k,h0,g):
    aip,aim = aipm(g,k,h0)
    aim1p,aim1m = aipm(g,k,h0)
    print(aip,aim)
    sp,sm = Spm2(k,dx)
    print(sp,sm)
    G = G2(h0,k,dx)
    print(G)
    F = FG(k,dx,aip,aim,aim1p,aim1m,sp,sm)
    print(F)
    return v(dx,g,h0,F,G)
    
def wactual(k,g,h0):
    from scipy import sqrt
    w1 = k*sqrt(g*h0)*sqrt(3.0 / (k*k*h0*h0 + 3))
    return w1,-w1

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
    