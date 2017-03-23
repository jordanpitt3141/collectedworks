# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 07:07:18 2016

@author: jordan
"""
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from scipy.linalg import norm

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
    
    F = matrix([[Fnn1,Fnu1], [Fun1,Fuu1]])
    
    M = eye(2) - dt*F
    
    return M


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
    
    F = matrix([[Fnn1,Fnu1], [Fun1,Fuu1]])
    
    M = 0.5*(2*eye(2) - 2*dt*F + dt*dt*F*F)
    
    return M

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
    
    F = matrix([[Fnn1,Fnu1], [Fun1,Fuu1]])
    
    M = (eye(2) - dt*F + 0.5*dt*dt*F*F + (1.0/6.0)*dt*dt*dt*F*F*F)
    
    return M

def wactual(k,g,h0):
    from scipy import sqrt
    w1 = k*sqrt(g*h0)*sqrt(3.0 / (k*k*h0*h0 + 3))
    return w1,-w1


def advectionFD(k,dt,dx):
    idx = 1.0 / dx
    im = sqrt(-1)
    p1 = abs(1 - dt*idx*im*sin(k*dx))**2
    p2 = abs(cos(k*dx) - dt*idx*im*sin(k*dx))**2
    
    return p1,p2
h = 1.0
k = 1
g = 9.81

w1,w2 = wactual(k,g,h)

#calculate the phase speeds by the relation omega/k
v1 = w1/k
v2 = w2/k


#this then measures the dispersion relation from the highest resolution soliton problem until 0.5
dxs = arange(1,0.01,-0.01) 
n = len(dxs)


for i in range(n):
    dx = dxs[i]
    Cr = 0.5
    l = Cr/v1
    dt = l*dx
    
    #M = o1(k,g,h,dx,dt)
    #N = o2(k,g,h,dx,dt)
    #P = o3(k,g,h,dx,dt)
    print("dx",dx,"M",norm(M,ord=2), "N",norm(N,ord=2) ,"P",norm(P,ord=2) )


    
    

    