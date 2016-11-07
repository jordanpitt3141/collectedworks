# -*- coding: utf-8 -*-
"""
Created on Thu Sep 11 16:26:29 2014

@author: Jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from Hamil import *
from os import listdir

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
    
def dambreaksmooth(x,x0,base,eta0,diffuse,dx):
    from numpy import tanh
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    
    for i in range(n):
        h[i] = base + 0.5*eta0*(1 + tanh(diffuse*(x0 - abs(x[i]))))
    return h,u 
    
def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 

dxw = "3"
wdirord = "o1"
wdatadir = "../../../../data/raw/bigsmoothtargetted/"  +wdirord +"/"


dxws = listdir(wdatadir)
dxws.sort(key=int)
Evals = []

sdir1 = "../../../../data/postprocessing/dbEnergy/"
if not os.path.exists(sdir1):
    os.makedirs(sdir1)

s = sdir1 + "Evalsinitial.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['dx' ,'beta','Eval'])        

for dxw in dxws:
    wdatadirn = wdatadir + dxw + "/" 
    nums = listdir(wdatadirn)
    nums.sort(key=int)
    for k in nums:
        wdir = wdatadirn  + str(k) + "/" 
        sdir = "../../../../data/postprocessing/dbEnergy/" + wdirord+  "/" + dxw + "/" + str(k) + "/"
        #if not os.path.exists(sdir):
        #    os.makedirs(sdir) 
             
        s = wdir + "outlast.txt"
        with open(s,'r') as file1:
            lines = file1.readlines()
            rl = lines[1].split(",")
            dx = float(rl[0])
            dt = float(rl[1])
            beta = float(rl[7])
      
        l = 0.01
        dt = l*dx
        startx = 0.0
        endx = 1000.0 + dx
        startt = 0.0
        endt = 30.0+(dt*0.9)   
        g = 9.81
                    
        x,t = makevar(startx,endx,dx,startt,endt,dt)
        n = len(x)
            
        bot = 0.0
        hf = 1.8
        hl = 1.0
        base = hl
        eta0 = hf - hl
        x0 = 500
        diffuse = beta
         
        niBC = 3     
        h,u= dambreaksmooth(x,x0,base,eta0,diffuse,dx)  
        xbeg = arange(startx - niBC*dx,startx,dx)
        xend = arange(endx + dx,endx + (niBC+1)*dx) 
        hbeg = h[0]*ones(niBC)
        hend = h[-1]*ones(niBC)
        ubeg = u[0]*ones(niBC)
        uend = u[-1]*ones(niBC)
        xbc =  concatenate([xbeg,array(x),xend])
        hbc =  concatenate([hbeg,array(h),hend])
        ubc =  concatenate([ubeg,array(u),uend])
        
        xbc_c = copyarraytoC(xbc)
        hbc_c = copyarraytoC(hbc)
        ubc_c = copyarraytoC(ubc)
        
        Eval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)
        Evals.append((Eval,wdirord,dxw,beta))
        
        s = sdir1 + "Evalsinitial.txt"
        with open(s,'a') as file2:
             writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
             writefile2.writerow([str(dx) ,str(beta),str(Eval)])   
        
        