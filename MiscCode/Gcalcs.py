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
    
def GFEM2A(h0,k,Rp,Rm,dx):
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
    
def GFEM2R(h0,k,Rp,Rm,dx):
    from scipy import cos
    i3 = 1.0 / 3.0
    im = sqrt(-1)
    idx = 1.0 / dx
    Recon = 1.0 / (e**(-im*k*dx)*Rp + 2*Rm + 2*Rp +  e**(im*k*dx)*Rm)
    Hmult = e**(-im*k*dx) + 4 + e**(im*k*dx)
    H3mult = 2*idx*idx*(-e**(-im*k*dx) + 2 - e**(im*k*dx))
    G1 = (h0*Hmult + h0*h0*h0*H3mult)*Recon*e**(0.5*im*k*dx)
    
    print(Recon*H3mult, 1.0/Recon - Hmult)
    
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





def o1(k,g,h,dx):
    from scipy import log
    im = sqrt(-1)


    #calculate the elementary contribution factors
    M =   1 
    
    Rp1 = Rpo1(k,dx)
    Rm1 = Rmo1(k,dx)
    Ru = Ruo2(k,dx)
    
    G1 = G2(h,k,dx)  

    
    return G1
    
def o2(k,g,h,dx):
    from scipy import log
    im = sqrt(-1)

    #calculate the elementary contribution factors
    M =   1 
    
    Rp1 = Rpo2(k,dx)
    Rm1 = Rmo2(k,dx)
    Ru = Ruo2(k,dx)
    
    G1 = G2(h,k,dx) 

    
    return G1

def oFEM2A(k,g,h,dx):
    from scipy import log
    im = sqrt(-1)

    #calculate the elementary contribution factors
    M =   1 
    
    Rp1 = Rpo2(k,dx)
    Rm1 = Rmo2(k,dx)
    Ru = Ruo2(k,dx)
    
    G1 = GFEM2A(h,k,Rp1,Rm1,dx) 

    
    return G1
    
def oFEM2R(k,g,h,dx):
    from scipy import log
    im = sqrt(-1)

    #calculate the elementary contribution factors
    M =   1 
    
    Rp1 = Rpo2(k,dx)
    Rm1 = Rmo2(k,dx)
    Ru = Ruo2(k,dx)
    
    #print(Rp1,Rm1,Ru, e**(im*k*dx*0.5))
    
    G1 = GFEM2R(h,k,Rp1,Rm1,dx) 

    
    return G1
    
def o3(k,g,h,dx):
    from scipy import log
    im = sqrt(-1)

    #calculate the elementary contribution factors
    M =   M3(k,dx)  
    
    Rp1 = Rpo3(k,M,dx)
    Rm1 = Rmo3(k,M,dx)
    Ru = Ruo3(k,dx)
    
    G1 = G4(h,k,dx)

    
    return G1


def oFEM3(k,g,h,dx):
    from scipy import log
    im = sqrt(-1)

    #calculate the elementary contribution factors
    M =   M3(k,dx)  
    
    Rp1 = Rpo3(k,M,dx)
    Rm1 = Rmo3(k,M,dx)
    Ru = Ruo3(k,dx)
    
    G1 = GFEM4(h,k,Rp1,Rm1,dx)

    
    return G1
    

from scipy import pi


h = 1.0
k = 0.5
g = 9.81

i3 = 1.0 / 3.0


#this then measures the dispersion relation from the highest resolution soliton problem until 0.5
dxs = linspace(100.0/ (2**20),0.5,num=1000)  

#dxs = [10**(-10), 10**(-8),10**(-6),10**(-4),10**(-2),1,10] 
n = len(dxs)
w1s = zeros(n)
w2s = zeros(n) 
G1s = zeros(n,dtype=complex)
G2s = zeros(n,dtype=complex)
G3s = zeros(n,dtype=complex)
Gf2as = zeros(n,dtype=complex)
Gf2rs = zeros(n,dtype=complex)
Gf3s = zeros(n,dtype=complex)
erro1 = zeros(n)
erro2 = zeros(n)
erro3 = zeros(n)
errof2a = zeros(n)
errof2r = zeros(n)
errof3 = zeros(n)

for i in range(n):
    dx = dxs[i]
    Ga = h + h**3*k*k*i3
    
    G1v = o1(k,g,h,dx)
    G2v = o2(k,g,h,dx)
    Gf2av = oFEM2A(k,g,h,dx)
    Gf2rv = oFEM2R(k,g,h,dx)
    G3v = o3(k,g,h,dx)
    Gf3v = oFEM3(k,g,h,dx)
    
    G1s[i] = G1v
    G2s[i] = G2v
    G3s[i] = G3v
    Gf2as[i] = Gf2av
    Gf2rs[i] = Gf2rv
    Gf3s[i] = Gf3v
  
    erro1[i] = abs(Ga - G1v) / abs(Ga)    
    erro2[i] = abs(Ga - G2v) / abs(Ga)    
    errof2a[i] = abs(Ga - Gf2av) / abs(Ga) 
    errof2r[i] = abs(Ga - Gf2rv) / abs(Ga)  
    erro3[i] = abs(Ga - G3v) / abs(Ga)        
    errof3[i] = abs(Ga - Gf3v) / abs(Ga)    
    
    

    