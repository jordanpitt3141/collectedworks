# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 23:36:31 2015

@author: jordan
"""
from scipy import *
import matplotlib.pyplot as plt
import os
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 

def dambreaksmooth(x,x0,base,eta0,diffuse,dx):
    from numpy import tanh
    h = base + 0.5*eta0*(1 + tanh((x0 - abs(x))/diffuse))
    return h
    
def dambreaksmoothallalphas(X,Y,x0,base,eta0,dx):
    from numpy import tanh
    xs = X.shape
    Z = zeros(xs)
    for i in range(xs[0]):
        for j in range(xs[1]):
            Z[i][j] = dambreaksmooth(X[i][j],x0,base,eta0,Y[i][j],dx)
    return Z

sdir = "../../../data/postprocessing/dbsmoothcontinuum/"
if not os.path.exists(sdir):
        os.makedirs(sdir)   
dx = 0.1
l = 0.01
dt = l*dx
theta = 1.2
startx = 490
endx = 510.0 + dx
startt = 0.0
endt = dt  
    
g = 9.81
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
base = 1.0
eta0 = 0.8
x0 = 500

alphas = arange(0.1,3,0.1)

Xn, Yn = meshgrid(x, alphas)

Z = dambreaksmoothallalphas(Xn,Yn,x0,base,eta0,dx)

"""
fig = plt.figure()
ax = fig.gca(projection='3d')

# Plot the surface.
surf = ax.plot_surface(Xn, Yn, Z, cmap=cm.Spectral,linewidth=0, antialiased=False)
#ax.plot_wireframe(Xn, Yn, Z,rstride=10, cstride=10)
plt.show()
"""

xs = Xn.shape
s = sdir  + "3D.dat"
with open(s,'w') as file1:
    #s ="%s%5s%s%5s%s\n" %("x"," ","y"," ","z")
    #file1.write(s)
    for i in range(xs[0]):
        for j in range(xs[1]):
            s ="%3.8f%5s%1.15f%5s%1.15f\n" %(Xn[i][j]," ",Yn[i][j]," ",Z[i][j])
            file1.write(s)
        file1.write("\n")
    



  