# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:09:21 2015

@author: jordan
"""
from Serre2dc import *
from scipy import *
import csv
import os
from numpy.linalg import norm  

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
    
def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 
        
#gives exact up to linears, so is second order accurate huzzah    
def getGfromupy(h,u,u0,u1,h0,h1,dx):
    idx = 1.0 / dx
    ithree = 1.0 / 3.0
        
    n = len(h)

    G = zeros(n)
        
    for i in range(1,n-1):
        th = h[i]
        thx = 0.5*idx*(h[i+1] - h[i-1])
        
        D = th
        
        ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
        bi = D + 2.0*ithree*idx*idx*th*th*th
        ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
        
        G[i] = ai*u[i-1] + bi*u[i] + ci*u[i+1]
        
    #boundary    
    #i=0
    i=0
    th = h[i]
    thx = 0.5*idx*(h[i+1] - h0)
            
    D = th
            
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
 
    G[i] = ai*u0 + bi*u[i] + ci*u[i+1]
    
    #i = n-1
    i = n-1

    th = h[i]
    thx = 0.5*idx*(h1 - h[i-1])
        
    D = th
        
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx

    G[i] = ai*u[i-1] + bi*u[i] + ci*u1
            
    return G  
    

def dambreak(x,hf,hc,hl,dx):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    
    for i in range(n):
        if (x[i] < hc):
            h[i] = hf
        else:
            h[i] = hl
    
    G = getGfromupy(h,u,0.0,0.0,h[0],h[-1],dx)
    return h,G       

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
         
    G = getGfromupy(h,u,0.0,0.0,a0,a0,dx)
    
    return h,G 
    
def experiment1(x,h0,h1,dx):
    n = len(x)
    u = zeros(n)
    h = ones(n)*h1
    for i in range(n):
        if (x[i] <0 and x[i] > -2*b):
            h[i] = h0
    G = getGfromupy(h,u,0.0,0.0,h1,h1,dx)

    return h,G

   
## DAM BREAK time##########################
wdir = "../../data/time/o2dc/"
dx = 1.0
l = 0.01
dt = l*dx
theta = 1.2
startx = 0.0
endx = 1000.0 + dx
startt = 0.0
endt = 30 + dt  
    
g = 9.81
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
    
bot = 0.0
hf = 1.8
hl = 1.0
gap = max(5,int(0.5/dt))
    
h,G = dambreak(x,hf,500,hl,dx)
    
nBC = 3
nBCs = 4
u0 = zeros(nBCs)
u1 = zeros(nBCs)    
h0 = hf*ones(nBCs)
h1 = hl*ones(nBCs)
    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
h0_c  = copyarraytoC(h0)
h1_c  = copyarraytoC(h1)
u0_c  = copyarraytoC(u0)
u1_c  = copyarraytoC(u1)
u_c = mallocPy(n)
    
       
for i in range(1,len(t)):        
    evolvewrap(G_c,h_c,h0_c,h1_c,u0_c,u1_c,g,dx,dt,nBC,n,nBCs,theta)
    print (t[i])

getufromG(h_c,G_c,u0[-1],u1[0],h0[-1],h1[0], dx ,n,u_c)
u = copyarrayfromC(u_c,n)
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)        

deallocPy(u_c)   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(h0_c)
deallocPy(h1_c)
deallocPy(u0_c)
deallocPy(u1_c)
