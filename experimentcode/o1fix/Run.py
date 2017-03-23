# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:09:21 2015

@author: jordan
"""
from Serre1 import *
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
    

def dambreak(x,hf,hc,hl,bot,dx):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    bed = bot*ones(n)
    
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

def soliton2interactinit(n,a0,a11,solbeg1,solend1,direction1,a12,solbeg2,solend2,direction2,g,x,t0,dx):
    h = zeros(n)
    u = zeros(n)
    c1 = sqrt(g*(a0 + a11))
    c2 = sqrt(g*(a0 + a11))
    for i in range(n):
        if (x[i] > solbeg1 and x[i] < solend1):
            h[i] = soliton(abs(x[i] - 0.5*(solbeg1 + solend1)),t0,g,a0,a11)
            u[i] = direction1*c1*( (h[i] - a0) / h[i] )
        elif (x[i] > solbeg2 and x[i] < solend2):
            h[i] = soliton(abs(x[i] - 0.5*(solbeg2 + solend2)),t0,g,a0,a12)
            u[i] =  direction2*c2* ((h[i] - a0) / h[i])
        else:
            h[i] = a0
            u[i] = 0.0
    G = getGfromupy(h,u,0.0,0.0,a0,a0,dx)
    return h,G  

def experiment1(x,b,h0,h1,dx):
    n = len(x)
    u = zeros(n)
    h = ones(n)*h1
    for i in range(n):
        if (x[i] <0 and x[i] > -2*b):
            h[i] = h0
    G = getGfromupy(h,u,0.0,0.0,h1,h1,dx)

    return h,G  

def dambreaksmooth(x,x0,base,eta0,diffuse,dx):
    from numpy import tanh
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    
    for i in range(n):
        h[i] = base + 0.5*eta0*(1 + tanh(diffuse*(x0 - abs(x[i]))))
    
    G = getGfromupy(h,u,0.0,0.0,h[0],h[-1],dx)
    return h,G       
    
            
"""
## DB
from time import time

wdir = "../../data/trackleadsol/o1/"

if not os.path.exists(wdir):
    os.makedirs(wdir)
    
hf = 1.8
hl = 1.0
g = 9.81

dx = 10.0 /(2**10)
Cr = 0.5
l = 0.01

dt = l*dx

startx = 0.0
endx = 1000.0 + dx
startt = 0.0
endt = 100 + dt  
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
    
bot = 0.0
    
diffuse = 1000
base = hl
eta0 = hf - hl
x0 = 500
h,G = dambreaksmooth(x,x0,base,eta0,diffuse,dx)     
   
nBC = 3
nBCs = 4
u0 = zeros(nBCs)
u1 = zeros(nBCs)    
h0 = hf*ones(nBCs)
h1 = hl*ones(nBCs)

gap = int(0.5/dt)
    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
h0_c  = copyarraytoC(h0)
h1_c  = copyarraytoC(h1)
u0_c  = copyarraytoC(u0)
u1_c  = copyarraytoC(u1)
u_c = mallocPy(n)

aplus = []
aplusx = []
aplust = [] 
    
tim1 =time()    
for i in range(1,len(t)):
    if(i == 1 or i%gap ==0):
        getufromG(h_c,G_c,u0[-1],u1[0],h0[-1],h1[0], dx ,n,u_c)
        u = copyarrayfromC(u_c,n)
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        
        mi = n - 2
        for mi in range(n-1,-1,-1):
            if(h[mi -1] < h[mi]) and (h[mi] > 1.1 ):
                break
        aplus.append(h[mi])
        aplusx.append(x[mi])
        aplust.append(t[i])
        s = wdir + "out" + str(i) + ".txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
            writefile2.writerow(['dx' ,'dt','time' ,"cell midpoint", 'height(m)', 'G' , 'u(m/s)'])        
                   
            for j in range(n):
                writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]) ,str(h[j]) , str(G[j]) , str(u[j])])   
    evolvewrap(G_c,h_c,h0_c,h1_c,u0_c,u1_c,g,dx,dt,nBC,n,nBCs)
    print t[i]
        
getufromG(h_c,G_c,u0[-1],u1[0],h0[-1],h1[0], dx ,n,u_c)
u = copyarrayfromC(u_c,n)
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)

mi = n - 2
for mi in range(n-1,-1,-1):
    if(h[mi -1] < h[mi]) and (h[mi] > 1.1 ):
        break
aplus.append(h[mi])
aplusx.append(x[mi])
aplust.append(t[i])
s = wdir + "outlast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     writefile2.writerow(['dx' ,'dt','time' ,"cell midpoint", 'height(m)', 'G' , 'u(m/s)'])        
                   
     for j in range(n):
         writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]) ,str(h[j]) , str(G[j]) , str(u[j])])   

s = wdir + "aplus.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(['x' ,'t','aplus' ,"Grim"])        
           
    for j in range(len(aplus)):
        writefile2.writerow([str(aplusx[j]),str(aplust[j]),str(aplus[j]),str(0.739976603390100695296990254)])              
deallocPy(u_c)   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(h0_c)
deallocPy(h1_c)
deallocPy(u0_c)
deallocPy(u1_c)
"""


"""
###################### DAM BREAK ACCURACY ################################
wdir = "../../../data/raw/dbh/o1/"
for k in range(20):
    
    hf = 1.8
    hl = 1.0
    g = 9.81
    
    dx = 100.0 / (2**k)
    Cr = 0.5
    l = Cr / sqrt(g*hf)
    dt = l*dx
    startx = 0.0
    endx = 1000.0 + dx
    startt = 0.0
    endt = 30 + dt   
    g = 9.81
    
    wdatadir = wdir+ str(k) + "/"
    if not os.path.exists(wdatadir):
        os.makedirs(wdatadir)
        
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
        
    bot = 0.0
    hf = 1.8
    hl = 1.0
        
    h,G= dambreak(x,hf,500,hl,bot,dx)
       
    nBC = 3
    nBCs = 4
    u0 = zeros(nBCs)
    u1 = zeros(nBCs)    
    h0 = hf*ones(nBCs)
    h1 = hl*ones(nBCs)
    
    gap = max(1,int(0.5/dt))
        
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    h0_c  = copyarraytoC(h0)
    h1_c  = copyarraytoC(h1)
    u0_c  = copyarraytoC(u0)
    u1_c  = copyarraytoC(u1)
    u_c = mallocPy(n)
        
        
        
    for i in range(1,len(t)):
        if(i == 1 or i%gap ==0):
            getufromG(h_c,G_c,u0[-1],u1[0],h0[-1],h1[0], dx ,n,u_c)
            u = copyarrayfromC(u_c,n)
            G = copyarrayfromC(G_c,n)
            h = copyarrayfromC(h_c,n)
            s = wdatadir + "out" + str(i) + ".txt"
            with open(s,'a') as file2:
                writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
                writefile2.writerow(['dx' ,'dt','time' ,"cell midpoint", 'height(m)', 'G' , 'u(m/s)'])        
                       
                for j in range(n):
                    writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]) ,str(h[j]) , str(G[j]) , str(u[j])])   
        evolvewrap(G_c,h_c,h0_c,h1_c,u0_c,u1_c,g,dx,dt,nBC,n,nBCs)
        print t[i]
            
    getufromG(h_c,G_c,u0[-1],u1[0],h0[-1],h1[0], dx ,n,u_c)
    u = copyarrayfromC(u_c,n)
    G = copyarrayfromC(G_c,n)
    h = copyarrayfromC(h_c,n)
    s = wdatadir + "outlast.txt"
    with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
         writefile2.writerow(['dx' ,'dt','time' ,"cell midpoint", 'height(m)', 'G' , 'u(m/s)'])        
                       
         for j in range(n):
             writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]) ,str(h[j]) , str(G[j]) , str(u[j])])   
                 

    s = wdir +  str(k)+ ".txt"
    with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
         writefile2.writerow(['dx' ,'dt','time' ,"cell midpoint", 'height(m)', 'G' , 'u(m/s)'])        
                   
         for j in range(n):
             writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]) ,str(h[j]) , str(G[j]) , str(u[j])])   
             
    deallocPy(u_c)   
    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(h0_c)
    deallocPy(h1_c)
    deallocPy(u0_c)
    deallocPy(u1_c)
"""

#big smooth
"""
wdirb = "../../data/bigsmoothalphainf/o1/"
import os
#diffuses = [0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1.0,2.5,5.0,7.5,10.0,25.0,50.0,75.0,100.0,250.0,500.0,750.0,1000.0]

diffuses = [2.5,10.0,1000.0]
dxs = []
#diffuses = [10.0,1000.0]

#for ll in range(3,16):
for ll in range(14,15):
    for k in range(len(diffuses)):
        wdir = wdirb + str(ll) + "/" + str(k) + "/"
        dx = 10.0 / (2**ll)
        l = 0.01
        dt = l*dx
        startx = 0.0
        endx = 1000.0 + dx
        startt = 0.0
        endt = 30.0+(dt*0.9)   
        g = 9.81
        
        gap = max(1,int(0.02/dt))
                    
        x,t = makevar(startx,endx,dx,startt,endt,dt)
        n = len(x)
            
        bot = 0.0
        hf = 1.8
        hl = 1.0
        base = hl
        eta0 = hf - hl
        x0 = 500
        diffuse = diffuses[k]
        
        if not os.path.exists(wdir):
            os.makedirs(wdir) 
            
        h,G= dambreaksmooth(x,x0,base,eta0,diffuse,dx)
           
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
            evolvewrap(G_c,h_c,h0_c,h1_c,u0_c,u1_c,g,dx,dt,nBC,n,nBCs)
            print t[i]
                
        getufromG(h_c,G_c,u0[-1],u1[0],h0[-1],h1[0], dx ,n,u_c)
        u = copyarrayfromC(u_c,n)
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        s = wdir + "outlast.txt"
        with open(s,'a') as file2:
             writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
             writefile2.writerow(['dx' ,'dt','time' ,"cell midpoint", 'height(m)', 'G' , 'u(m/s)', 'diffuse'])        
                           
             for j in range(n):
                 writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]) ,str(h[j]) , str(G[j]) , str(u[j]), str(diffuse)])   
                 
        deallocPy(u_c)   
        deallocPy(h_c)
        deallocPy(G_c)
        deallocPy(h0_c)
        deallocPy(h1_c)
        deallocPy(u0_c)
        deallocPy(u1_c)
"""


#big smooth targetted
"""
difflist = [1]

deltaxa = [5,6,7,9,10,11,12,13,14,15]
dxlist = [deltaxa]

import os
diffuses = [0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1.0,2.5,5.0,7.5,10.0,25.0,50.0,75.0,100.0,250.0,500.0,750.0,1000.0]
wdirb = "../../data/bigsmoothtargetted/o1/"
for lk in range(len(difflist)):
    for ll in dxlist[lk]:
        wdir = wdirb + str(ll) + "/" + str(difflist[lk]) + "/"
        dx = ll*(10.0 / (2**12))
        l = 0.01
        dt = l*dx
        startx = 0.0
        endx = 1000.0 + dx
        startt = 0.0
        endt = 30.0+(dt*0.9)   
        g = 9.81
        
        gap = max(1,int(0.02/dt))
                    
        x,t = makevar(startx,endx,dx,startt,endt,dt)
        n = len(x)
            
        bot = 0.0
        hf = 1.8
        hl = 1.0
        base = hl
        eta0 = hf - hl
        x0 = 500
        diffuse = diffuses[difflist[lk]]
        
        if not os.path.exists(wdir):
            os.makedirs(wdir) 
            
        h,G= dambreaksmooth(x,x0,base,eta0,diffuse,dx)
           
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
            evolvewrap(G_c,h_c,h0_c,h1_c,u0_c,u1_c,g,dx,dt,nBC,n,nBCs)
            print t[i]
                
        getufromG(h_c,G_c,u0[-1],u1[0],h0[-1],h1[0], dx ,n,u_c)
        u = copyarrayfromC(u_c,n)
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        s = wdir + "outlast.txt"
        with open(s,'a') as file2:
             writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
             writefile2.writerow(['dx' ,'dt','time' ,"cell midpoint", 'height(m)', 'G' , 'u(m/s)', 'diffuse'])        
                           
             for j in range(n):
                 writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]) ,str(h[j]) , str(G[j]) , str(u[j]), str(diffuse)])   
                 
        deallocPy(u_c)   
        deallocPy(h_c)
        deallocPy(G_c)
        deallocPy(h0_c)
        deallocPy(h1_c)
        deallocPy(u0_c)
        deallocPy(u1_c)
"""

"""
################################# SOLITON Accuracy ####################3
#dxs = [100.0,90.0,80.0,70.0,60.0,50.0,40.0,30.0,20.0,10.0,9.0,8.0,7.0,6.0,5.0,4.0,3.0,2.0,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.09,0.08,0.07 \
#,0.06,0.05,0.04,0.03,0.02,0.01]
wdir = "../../../data/raw/Solnon0p7/o1/"

if not os.path.exists(wdir):
    os.makedirs(wdir)

s = wdir + "savenorms.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(['dx','Normalised L1-norm Difference Height', ' Normalised L1-norm Difference Velocity','Hamiltonian Difference'])
    
for k in range(6,21):
    dx = 100.0 / (2**k)
    a0 = 1.0
    a1 = 0.7
    g = 9.81
    #g = 1
    Cr = 0.5
    l = 1.0 / (sqrt(g*(a0 + a1)))
    dt = Cr*l*dx
    startx = -250.0
    endx = 250.0 + dx
    startt = 0
    endt = 50 + dt
    
    wdatadir = wdir+ str(k) + "/"
    if not os.path.exists(wdatadir):
        os.makedirs(wdatadir)
    
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
    
    t0 = 0
    bot = 0
    gap =int(10.0/dt)
    
    h,G= solitoninit(n,a0,a1,g,x,t0,dx)
    
    nBC = 3
    nBCs = 4
    niBC = nBCs
    u0 = zeros(nBCs)
    u1 = zeros(nBCs)    
    h0 = h[0]*ones(nBCs)
    h1 = h[-1]*ones(nBCs)
    
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    h0_c  = copyarraytoC(h0)
    h1_c  = copyarraytoC(h1)
    u0_c  = copyarraytoC(u0)
    u1_c  = copyarraytoC(u1)
    u_c = mallocPy(n)
    
    xbeg = arange(startx - niBC*dx,startx,dx)
    xend = arange(endx + dx,endx + (niBC+1)*dx) 
    
    xbc =  concatenate([xbeg,x,xend])  
    
    xbc_c = copyarraytoC(xbc)
    hbc_c = mallocPy(n + 2*niBC)
    ubc_c = mallocPy(n + 2*niBC)
    Evals = []
    
    
    
    for i in range(1,len(t)):
        
        if(i % gap == 0 or i ==1):
            getufromG(h_c,G_c,u0[-1],u1[0],h0[-1],h1[0], dx ,n,u_c)
            
            conc(h0_c , h_c,h1_c,niBC,n ,niBC , hbc_c)
            conc(u0_c , u_c,u1_c,niBC,n ,niBC , ubc_c)         
            Eval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)
            
            Evals.append(Eval)
            u = copyarrayfromC(u_c,n)
            G = copyarrayfromC(G_c,n)
            h = copyarrayfromC(h_c,n)
            
            c = sqrt(g*(a0 + a1))
            htrue = zeros(n)
            utrue = zeros(n)
            for j in range(n):             
                he = soliton(x[j],t[i],g,a0,a1)
                htrue[j] = he
                utrue[j] = c* ((he - a0) / he) 
                
            s = wdatadir + "saveoutputts" + str(i) + ".txt"

            with open(s,'a') as file2:
                writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
                writefile2.writerow(['dx' ,'dt','time','xi','Eval', 'height(m)', 'G' , 'u(m/s)','true height', 'true velocity'  ])        
                   
                for j in range(n):
                    writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]),str(Eval), str(h[j]) , str(G[j]) , str(u[j]), str(htrue[j]), str(utrue[j])])  
                
        print t[i]
        print(h[1],G[1])     
        evolvewrap(G_c,h_c,h0_c,h1_c,u0_c,u1_c,g,dx,dt,nBC,n,nBCs)
        
    getufromG(h_c,G_c,u0[-1],u1[0],h0[-1],h1[0], dx ,n,u_c)
    
    conc(h0_c , h_c,h1_c,niBC,n ,niBC , hbc_c)
    conc(u0_c , u_c,u1_c,niBC,n ,niBC , ubc_c)        
    Eval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)
    
    Evals.append(Eval)
    u = copyarrayfromC(u_c,n)
    G = copyarrayfromC(G_c,n)
    h = copyarrayfromC(h_c,n)
    
    
    c = sqrt(g*(a0 + a1))
    htrue = zeros(n)
    utrue = zeros(n)
    for j in range(n):             
        he = soliton(x[j],t[-1],g,a0,a1)
        htrue[j] = he
        utrue[j] = c* ((he - a0) / he) 
    
    s = wdatadir + "saveoutputtslast.txt"
    with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
         writefile2.writerow(['dx' ,'dt','time','xi','Eval', 'height(m)', 'G' , 'u(m/s)','true height', 'true velocity'  ])        
                   
         for j in range(n):
             writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]),str(Eval), str(h[j]) , str(G[j]) , str(u[j]), str(htrue[j]), str(utrue[j])])  
                
    normhdiffi = norm(h - htrue,ord=1) / norm(htrue,ord=1)
    normudiffi = norm(u -utrue,ord=1) / norm(utrue,ord=1)  
    normHamdiff = (Evals[-1] - Evals[0])/ Evals[0]

    s = wdir + "savenorms.txt"
    with open(s,'a') as file1:
        writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

        writefile.writerow([str(dx),str(normhdiffi), str(normudiffi),str(normHamdiff)])      
        
    deallocPy(u_c)   
    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(h0_c)
    deallocPy(h1_c)
    deallocPy(u0_c)
    deallocPy(u1_c)
"""

""" 
##### SOLITON
dx = 0.05
a0 = 1.0
a1 = 1.0
g = 9.81
Cr = 0.5
l = Cr / (sqrt(g*(a0 + a1)))
dt = l*dx
startx = -100.0
endx = 500.0
startt = 0.0
endt = 100 + dt

    
wdir = "../../../data/raw/Cserre/solitonothers/highnonlinear/order1/dx0p05"
if not os.path.exists(wdir):
    os.makedirs(wdir)
    
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
    
t0 = 0
gap = int(10.0/dt)
    
h,G= solitoninit(n,a0,a1,g,x,t0,dx)
    
nBC = 3
nBCs = 4
u0 = zeros(nBCs)
u1 = zeros(nBCs)    
h0 = h[0]*ones(nBCs)
h1 = h[-1]*ones(nBCs)
    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
h0_c  = copyarraytoC(h0)
h1_c  = copyarraytoC(h1)
u0_c  = copyarraytoC(u0)
u1_c  = copyarraytoC(u1)
u_c = mallocPy(n)
    
    
    
for i in range(1,len(t)):
       
    if(i % gap == 0 or i ==1):
       getufromG(h_c,G_c,u0[-1],u1[0],h0[-1],h1[0],dx,n,u_c)
       u = copyarrayfromC(u_c,n)
       G = copyarrayfromC(G_c,n)
       h = copyarrayfromC(h_c,n)
 
       c = sqrt(g*(a0 + a1))
       htrue = zeros(n)
       utrue = zeros(n)
       for j in range(n):             
           he = soliton(x[j],t[i],g,a0,a1)
           htrue[j] = he
           utrue[j] = c* ((he - a0) / he) 
                
       s = wdir + "saveoutputts" + str(i) + ".txt"
       print t[i]
       print(h[1],G[1]) 
       with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
            writefile2.writerow(['dx' ,'dt','time', 'height(m)', 'G' , 'u(m/s)','true height', 'true velocity' ])        
                   
            for j in range(n):
                writefile2.writerow([str(dx),str(dt),str(t[i]), str(h[j]) , str(G[j]) , str(u[j]) , str(htrue[j]), str(utrue[j])])  
                 
            
    evolvewrap(G_c,h_c,h0_c,h1_c,u0_c,u1_c,g,dx,dt,nBC,n,nBCs)
        
getufromG(h_c,G_c,u0[-1],u1[0],h0[-1],h1[0], dx ,n,u_c)
u = copyarrayfromC(u_c,n)
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)
    
    
c = sqrt(g*(a0 + a1))
htrue = zeros(n)
utrue = zeros(n)
for j in range(n):             
    he = soliton(x[j],t[-1],g,a0,a1)
    htrue[j] = he
    utrue[j] = c* ((he - a0) / he) 
    
s = wdir + "saveoutputtslast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
     writefile2.writerow(['dx' ,'dt','time', 'height(m)', 'G' , 'u(m/s)','true height', 'true velocity'  ])        
                   
     for j in range(n):
         writefile2.writerow([str(dx),str(dt),str(t[-1]), str(h[j]) , str(G[j]) , str(u[j]), str(htrue[j]), str(utrue[j])])       
    

deallocPy(u_c)   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(h0_c)
deallocPy(h1_c)
deallocPy(u0_c)
deallocPy(u1_c)  
"""


##### SOLITON INTERACTION
wdir = "../../../data/raw/SolColAGNH/o1/"
if not os.path.exists(wdir):
    os.makedirs(wdir)
    
dx = 100.0/ (2**11)

a0 = 1.0
a11 = 0.7
solbeg1 = 100.0
solend1 = 200.0
direction1 = 1.0
a12 = 0.7
solbeg2 = 200.0
solend2 = 300.0
direction2 = -1.0

Cr = 0.5
g = 9.81
#g = 1.0
l = Cr / (sqrt(g*1.5*(a0 + a11 + a12)))
dt = l*dx
startx = -100.0
endx = 500.0 + dx
startt = 0.0
endt = 50 + dt

    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
    
t0 = 0
gap = int(0.5/dt)
    
h,G = soliton2interactinit(n,a0,a11,solbeg1,solend1,direction1,a12,solbeg2,solend2,direction2,g,x,t0,dx)
    
nBC = 3
nBCs = 4
niBC = nBC
u0 = zeros(nBCs)
u1 = zeros(nBCs)    
h0 = h[0]*ones(nBCs)
h1 = h[-1]*ones(nBCs)
    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
h0_c  = copyarraytoC(h0)
h1_c  = copyarraytoC(h1)
u0_c  = copyarraytoC(u0)
u1_c  = copyarraytoC(u1)
u_c = mallocPy(n)

xbeg = arange(startx - niBC*dx,startx,dx)
xend = arange(endx + dx,endx + (niBC+1)*dx) 

xbc =  concatenate([xbeg,x,xend])  

xbc_c = copyarraytoC(xbc)
hbc_c = mallocPy(n + 2*niBC)
ubc_c = mallocPy(n + 2*niBC)
Evals = []
    
    
    
for i in range(1,len(t)):
       
    if(i % gap == 0 or i ==1):
       getufromG(h_c,G_c,u0[-1],u1[0],h0[-1],h1[0],dx,n,u_c)
       u = copyarrayfromC(u_c,n)
       G = copyarrayfromC(G_c,n)
       h = copyarrayfromC(h_c,n)

       conc(h0_c , h_c,h1_c,niBC,n ,niBC , hbc_c)
       conc(u0_c , u_c,u1_c,niBC,n ,niBC , ubc_c)         
       Eval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)
        
       Evals.append(Eval)
                 
       s = wdir + "saveoutputts" + str(i) + ".txt"
       print t[i]
       print(h[1],G[1]) 
       with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
            writefile2.writerow(['dx' ,'dt','time','Eval','x value', 'height(m)', 'G' , 'u(m/s)','true height', 'true velocity' ])        
                   
            for j in range(n):
                writefile2.writerow([str(dx),str(dt),str(t[i]),str(Eval), str(x[j]), str(h[j]) , str(G[j]) , str(u[j])])  
                 
            
    evolvewrap(G_c,h_c,h0_c,h1_c,u0_c,u1_c,g,dx,dt,nBC,n,nBCs)
        
getufromG(h_c,G_c,u0[-1],u1[0],h0[-1],h1[0], dx ,n,u_c)
u = copyarrayfromC(u_c,n)
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)

conc(h0_c , h_c,h1_c,niBC,n ,niBC , hbc_c)
conc(u0_c , u_c,u1_c,niBC,n ,niBC , ubc_c)        
Eval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)

Evals.append(Eval)
    
s = wdir + "saveoutputtslast.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    writefile2.writerow(['dx' ,'dt','time','Eval','x value', 'height(m)', 'G' , 'u(m/s)','true height', 'true velocity' ])        
       
    for j in range(n):
        writefile2.writerow([str(dx),str(dt),str(t[i]),str(Eval), str(x[j]), str(h[j]) , str(G[j]) , str(u[j])])  
    

deallocPy(u_c)   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(h0_c)
deallocPy(h1_c)
deallocPy(u0_c)
deallocPy(u1_c)  


"""
##### SOLITON time
from time import time
import os
from scipy.linalg import norm
dx = 100.0 / (2.0**10)
a0 = 10.0
a1 = 1.0
g = 9.81
Cr = 0.5
l = 1 / (sqrt(g*(a0 + a1)))
dt = Cr*l*dx
startx = -500.0
endx = 1000.0 + dx
startt = 0
endt = 50.0 + dt

theta = 1.2
    
wdir = "../../../data/raw/timecomplong/o1/"
if not os.path.exists(wdir):
    os.makedirs(wdir) 
    
g = 9.81
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
t0 = 0
bot = 0
gap = max(1,int(0.5/dt))
    
h,G= solitoninit(n,a0,a1,g,x,t0,dx)
    
nBC = 3
nBCs = 4
u0 = zeros(nBCs)
u1 = zeros(nBCs)    
h0 = h[0]*ones(nBCs)
h1 = h[-1]*ones(nBCs)
    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
h0_c  = copyarraytoC(h0)
h1_c  = copyarraytoC(h1)
u0_c  = copyarraytoC(u0)
u1_c  = copyarraytoC(u1)
u_c = mallocPy(n)
    
    
t0 = time()      
for i in range(1,len(t)):      
    evolvewrap(G_c,h_c,h0_c,h1_c,u0_c,u1_c,g,dx,dt,nBC,n,nBCs)
    print(t[i])
t1 = time()
timeelapse = t1 - t0
        
getufromG(h_c,G_c,u0[-1],u1[0],h0[-1],h1[0], dx ,n,u_c)
u = copyarrayfromC(u_c,n)
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)
    
    
c = sqrt(g*(a0 + a1))
htrue = zeros(n)
utrue = zeros(n)
for j in range(n):             
    he = soliton(x[j],t[-1],g,a0,a1)
    htrue[j] = he
    utrue[j] = c* ((he - a0) / he) 

normh = norm(h - htrue,ord=1) / norm(htrue,ord=1)    
normu = norm(u - utrue,ord=1) / norm(utrue,ord=1)     
s = wdir + "saveoutputtslast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['dx' ,'dt','time', 'height(m)', 'G' , 'u(m/s)','true height', 'true velocity',"hnorm","unorm", "time taken", "number of steps", "average time per step" ])        
               
     for j in range(n):
         writefile2.writerow([str(dx),str(dt),str(t[i]), str(h[j]) , str(G[j]) , str(u[j]), str(htrue[j]), str(utrue[j]),str(normh),str(normu), str(timeelapse),str(len(t) - 1) , str((1.0*timeelapse)/ (len(t) - 1) )])

deallocPy(u_c)   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(h0_c)
deallocPy(h1_c)
deallocPy(u0_c)
deallocPy(u1_c)  
"""


"""
################ SEGUR AND HAMMACK RECTANGULAR WAVE ###################################
tl = 60.0
b = 0.61
h0 = 0.09
h1 = 0.1
g = 9.81

dx = 0.01
Cr = 0.5
l = Cr / sqrt(g*h1)

dt = l*dx
startx = -tl
endx = tl + dx
startt = 0
endt = 50.0 + dt
    
wdir = "../../../data/raw/segur/o1af/"

if not os.path.exists(wdir):
    os.makedirs(wdir)
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
    
bot = 0
gap = 1
    
h,G = experiment1(x,b,h0,h1,dx)
  
nBC = 3
nBCs = 4
u0 = zeros(nBCs)
u1 = zeros(nBCs)    
h0 = h[0]*ones(nBCs)
h1 = h[-1]*ones(nBCs)
    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
h0_c  = copyarraytoC(h0)
h1_c  = copyarraytoC(h1)
u0_c  = copyarraytoC(u0)
u1_c  = copyarraytoC(u1)
u_c = mallocPy(n)

for i in range(1,len(t)):
    if(i == 1 or i%gap ==0):
        getufromG(h_c,G_c,u0[-1],u1[0],h0[-1],h1[0], dx ,n,u_c)
        u = copyarrayfromC(u_c,n)
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        s = wdir + "out" + str(i) + ".txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
            writefile2.writerow(['dx' ,'dt','time' ,"cell midpoint", 'height(m)', 'G' , 'u(m/s)'])        
                   
            for j in range(n):
                writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]) ,str(h[j]) , str(G[j]) , str(u[j])])   
    evolvewrap(G_c,h_c,h0_c,h1_c,u0_c,u1_c,g,dx,dt,nBC,n,nBCs)
    print t[i]
        
getufromG(h_c,G_c,u0[-1],u1[0],h0[-1],h1[0], dx ,n,u_c)
u = copyarrayfromC(u_c,n)
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)
s = wdir + "outlast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     writefile2.writerow(['dx' ,'dt','time' ,"cell midpoint", 'height(m)', 'G' , 'u(m/s)'])        
                   
     for j in range(n):
         writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]) ,str(h[j]) , str(G[j]) , str(u[j])])   
             
deallocPy(u_c)   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(h0_c)
deallocPy(h1_c)
deallocPy(u0_c)
deallocPy(u1_c)
"""