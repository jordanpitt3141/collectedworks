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
    
    return (26 - 2*cos(k*dx))/24.0
    #return 1
    
def wactual(k,g,h0):
    from scipy import sqrt
    w1 = k*sqrt(g*h0)*sqrt(3.0 / (k*k*h0*h0 + 3))
    return w1,-w1

from scipy import pi

h0 = 1.0
k = 2.5
g = 9.81

dx = 0.1
idx = 1.0 / dx

i12 = 1.0 / 12.0

im = sqrt(-1)

#second derivative check

SD2 = (idx*idx)*(2*(cos(k*dx) - 1))
SD4 = (idx*idx)*i12*(32*cos(k*dx) - 2*cos(2*k*dx) - 30)
SD = -k*k

print(SD, SD2 , SD4)

#reconstruction check
M = M3(k,dx)
Rpo3v = Rpo3(k,M,dx)
Rmo3v = Rmo3(k,M,dx)
Rpo2v = Rpo2(k,dx)
Rmo2v = Rmo2(k,dx)
Rpo1v = Rpo1(k,dx)
Rmo1v = Rmo1(k,dx)

R = e**(0.5*im*k*dx)
print("First Order")
print(Rpo1v,Rmo1v)
print("Second Order")
print(Rpo2v,Rmo2v)
print("third Order")
print(Rpo3v,Rmo3v)
print("Analytic")
print(R)