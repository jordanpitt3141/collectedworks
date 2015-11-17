# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:09:21 2015

@author: jordan
"""
from Serre2 import *
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
def getGfromupy(h,u,bed,u0,u1,h0,h1,b0,b1,dx):
    idx = 1.0 / dx
    ithree = 1.0 / 3.0
        
    n = len(h)

    G = zeros(n)
        
    for i in range(1,n-1):
        th = h[i]
        thx = 0.5*idx*(h[i+1] - h[i-1])
        tbx = 0.5*idx*(bed[i+1] - bed[i-1])
        tbxx = idx*idx*(bed[i+1] -2*bed[i] + bed[i-1])
        
        D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
        
        ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
        bi = D + 2.0*ithree*idx*idx*th*th*th
        ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
        
        G[i] = ai*u[i-1] + bi*u[i] + ci*u[i+1]
        
    #boundary    
    #i=0
    i=0
    th = h[i]
    thx = 0.5*idx*(h[i+1] - h0)
    tbx = 0.5*idx*(bed[i+1] - b0)
    tbxx = idx*idx*(bed[i+1] -2*bed[i] + b0)
            
    D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
            
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
 
    G[i] = ai*u0 + bi*u[i] + ci*u[i+1]
    
    #i = n-1
    i = n-1

    th = h[i]
    thx = 0.5*idx*(h1 - h[i-1])
    tbx = 0.5*idx*(b1 - bed[i-1])
    tbxx = idx*idx*(b1 -2*bed[i] + bed[i-1])
        
    D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
        
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
    
    G = getGfromupy(h,u,bed,0.0,0.0,h[0],h[-1],0.0,0.0,dx)
    return h,G,bed

       

def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a

def soliton (x,t,g,a0,a1):
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)
  
def solitoninit(n,a0,a1,g,x,t0,bot,dx):
    h = zeros(n)
    bx = zeros(n)
    u = zeros(n)
    c = sqrt(g*(a0 + a1))
    for i in range(n):
        bx[i] = bot
        h[i] = soliton(x[i],t0,g,a0,a1)
        u[i] =  c* ((h[i] - a0) / h[i])
         
    G = getGfromupy(h,u,bx,0.0,0.0,a0,a0,0.0,0.0,dx)
    
    return h,G,bx 
    
def experiment1(x,b,h0,h1,dx):
    n = len(x)
    u = zeros(n)
    bx = zeros(n)
    h = ones(n)*h1
    for i in range(n):
        if (x[i] <0 and x[i] > -2*b):
            h[i] = h0
    G = getGfromupy(h,u,bx,0.0,0.0,h1,h1,0.0,0.0,dx)

    return h,G,bx 

def dambreaksmooth(x,x0,base,eta0,diffuse,bot,dx):
    from numpy import tanh
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    bed = bot*ones(n)
    
    for i in range(n):
        h[i] = base + 0.5*eta0*(1 + tanh(diffuse*(x0 - abs(x[i]))))
    
    G = getGfromupy(h,u,bed,0.0,0.0,h[0],h[-1],0.0,0.0,dx)
    return h,G,bed   
"""
dx = 100.0 / (2**13)
l = 0.01
dt = l*dx
theta = 1.2
startx = 0.0
endx = 1000.0 + dx
startt = 0.0
endt = dt  
    
g = 9.81
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
bot = 0.0
base = 1.0
eta0 = 0.8
x0 = 500

for i in range(5):
    diffuse = (1000) / (10.0**i)
    hh,G,bed = dambreaksmooth(x,x0,base,eta0,diffuse,bot,dx)
    s = "diff = " + str(diffuse)
    plot(x,hh,label=s)
    
h,G,bed = dambreak(x,1.8,x0,1.0,0.0,dx)
s = "actual"
plot(x,h,label=s)
legend()
"""
"""
## smooth DAM BREAK time##########################
from time import time
import os
diffuses = [0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1.0,2.5,5.0,7.5,10.0,25.0,50.0,75.0,100.0,250.0,500.0,750.0,1000.0]
wdirb = "../../data/bigsmooth/o2/"
for ll in range(3,16):
    for k in range(len(diffuses)):
        wdir = wdirb + str(ll) + "/" + str(k) + "/"
        dx = 10.0 / (2**ll)
        l = 0.01
        dt = l*dx
        theta = 1.2
        startx = 0.0
        endx = 1000.0 + dx
        startt = 0.0
        endt = 30.0+(dt*0.9) 
        hf = 1.8
        hl = 1.0
        
        bot = 0.0
        base = hl
        eta0 = hf - hl
        x0 = 500
        diffuse = diffuses[k]
        
        if not os.path.exists(wdir):
            os.makedirs(wdir)  
            
        g = 9.81
            
        x,t = makevar(startx,endx,dx,startt,endt,dt)
        n = len(x)
            
        bot = 0.0
        gap = max(5,int(0.5/dt))
            
        h,G,bed = dambreaksmooth(x,x0,base,eta0,diffuse,bot,dx)
            
        nBC = 3
        nBCs = 4
        b0 = bot*ones(nBCs)
        b1 = bot*ones(nBCs)
        u0 = zeros(nBCs)
        u1 = zeros(nBCs)    
        h0 = hf*ones(nBCs)
        h1 = hl*ones(nBCs)
            
        h_c = copyarraytoC(h)
        G_c = copyarraytoC(G)
        bed_c = copyarraytoC(bed)
        h0_c  = copyarraytoC(h0)
        h1_c  = copyarraytoC(h1)
        u0_c  = copyarraytoC(u0)
        u1_c  = copyarraytoC(u1)
        b0_c  = copyarraytoC(b0)
        b1_c  = copyarraytoC(b1)
        u_c = mallocPy(n)
            
            
        for i in range(1,len(t)):            
            evolvewrap(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,b0_c,b1_c,g,dx,dt,nBC,n,nBCs,theta)
            print (t[i])
                
        getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
        u = copyarrayfromC(u_c,n)
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        s = wdir +  "outlast.txt"
        with open(s,'a') as file2:
             writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
             writefile2.writerow(['dx' ,'dt','time',"cell midpoint" ,'height(m)', 'G' , 'u(m/s)','bed',"diffuse" ])        
                           
             for j in range(n):
                 writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j]),str(diffuse)])
             
        deallocPy(u_c)   
        deallocPy(h_c)
        deallocPy(G_c)
        deallocPy(h0_c)
        deallocPy(h1_c)
        deallocPy(u0_c)
        deallocPy(u1_c) 
"""
"""
## smooth DAM BREAK targetted##########################
difflist = [1,6,8,9,12]

deltaxa = [5,6,7,9,10,11,12,13,14,15]
dxlist = [deltaxa,deltaxa,deltaxa,deltaxa,deltaxa]

diffuses = [0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1.0,2.5,5.0,7.5,10.0,25.0,50.0,75.0,100.0,250.0,500.0,750.0,1000.0]
wdirb = "../../data/bigsmoothtargetted/o2/"
for lk in range(len(difflist)):
    for ll in dxlist[lk]:
        wdir = wdirb + str(ll) + "/" + str(difflist[lk]) + "/"
        dx = ll*(10.0 / (2**11))
        l = 0.01
        dt = l*dx
        theta = 1.2
        startx = 0.0
        endx = 1000.0 + dx
        startt = 0.0
        endt = 30.0+(dt*0.9) 
        hf = 1.8
        hl = 1.0
        
        bot = 0.0
        base = hl
        eta0 = hf - hl
        x0 = 500
        diffuse = diffuses[difflist[lk]]
        
        if not os.path.exists(wdir):
            os.makedirs(wdir)  
            
        g = 9.81
            
        x,t = makevar(startx,endx,dx,startt,endt,dt)
        n = len(x)
            
        bot = 0.0
        gap = max(5,int(0.5/dt))
            
        h,G,bed = dambreaksmooth(x,x0,base,eta0,diffuse,bot,dx)
            
        nBC = 3
        nBCs = 4
        b0 = bot*ones(nBCs)
        b1 = bot*ones(nBCs)
        u0 = zeros(nBCs)
        u1 = zeros(nBCs)    
        h0 = hf*ones(nBCs)
        h1 = hl*ones(nBCs)
            
        h_c = copyarraytoC(h)
        G_c = copyarraytoC(G)
        bed_c = copyarraytoC(bed)
        h0_c  = copyarraytoC(h0)
        h1_c  = copyarraytoC(h1)
        u0_c  = copyarraytoC(u0)
        u1_c  = copyarraytoC(u1)
        b0_c  = copyarraytoC(b0)
        b1_c  = copyarraytoC(b1)
        u_c = mallocPy(n)
            
            
        for i in range(1,len(t)):            
            evolvewrap(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,b0_c,b1_c,g,dx,dt,nBC,n,nBCs,theta)
            print (t[i])
                
        getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
        u = copyarrayfromC(u_c,n)
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        s = wdir +  "outlast.txt"
        with open(s,'a') as file2:
             writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
             writefile2.writerow(['dx' ,'dt','time',"cell midpoint" ,'height(m)', 'G' , 'u(m/s)','bed',"diffuse" ])        
                           
             for j in range(n):
                 writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j]),str(diffuse)])
             
        deallocPy(u_c)   
        deallocPy(h_c)
        deallocPy(G_c)
        deallocPy(h0_c)
        deallocPy(h1_c)
        deallocPy(u0_c)
        deallocPy(u1_c) 
"""
""" 
## DAM BREAK time##########################
from time import time
wdir = "../../data/time/o2/"
samp = 100
dx = 0.001
l = 0.01
dt = l*dx
theta = 1.2
startx = 0.0
endx = 1000.0 + dx
startt = 0.0
endt = samp*dt  
    
g = 9.81
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
    
bot = 0.0
hf = 1.8
hl = 1.0
gap = max(5,int(0.5/dt))
    
h,G,bed = dambreak(x,hf,500,hl,bot,dx)
    
nBC = 3
nBCs = 4
b0 = bot*ones(nBCs)
b1 = bot*ones(nBCs)
u0 = zeros(nBCs)
u1 = zeros(nBCs)    
h0 = hf*ones(nBCs)
h1 = hl*ones(nBCs)
    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
bed_c = copyarraytoC(bed)
h0_c  = copyarraytoC(h0)
h1_c  = copyarraytoC(h1)
u0_c  = copyarraytoC(u0)
u1_c  = copyarraytoC(u1)
b0_c  = copyarraytoC(b0)
b1_c  = copyarraytoC(b1)
u_c = mallocPy(n)
    
    
tim1 =time()     
for i in range(1,len(t)):        
    evolvewrap(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,b0_c,b1_c,g,dx,dt,nBC,n,nBCs,theta)
    print (t[i])
tim2 =time()
tt = tim2 - tim1  
        
s = wdir + "outlast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     writefile2.writerow(['samp' ,'total time(s)','average time(s)'])        
     writefile2.writerow([str(samp), str(tt), str(tt/float(samp))]) 
deallocPy(u_c)   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(h0_c)
deallocPy(h1_c)
deallocPy(u0_c)
deallocPy(u1_c)  
"""

## DAM BREAK##########################
from time import time
wdir = "../../../data/dbChris/o2/"
samp = 100
g = 9.81
hf = 10.0
hl = 1.0
dx = 0.1
Cr = 0.2
l = Cr / sqrt(g*hf) 
dt = l*dx
theta = 1.2
startx = 0.0
endx = 1000.0 + dx
startt = 0.0
endt = 30.0 + dt  

if not os.path.exists(wdir):
    os.makedirs(wdir)
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
    
bot = 0.0
gap = max(5,int(0.1/dt))
    
h,G,bed = dambreak(x,hf,500,hl,bot,dx)
    
nBC = 3
nBCs = 4
b0 = bot*ones(nBCs)
b1 = bot*ones(nBCs)
u0 = zeros(nBCs)
u1 = zeros(nBCs)    
h0 = hf*ones(nBCs)
h1 = hl*ones(nBCs)
    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
bed_c = copyarraytoC(bed)
h0_c  = copyarraytoC(h0)
h1_c  = copyarraytoC(h1)
u0_c  = copyarraytoC(u0)
u1_c  = copyarraytoC(u1)
b0_c  = copyarraytoC(b0)
b1_c  = copyarraytoC(b1)
u_c = mallocPy(n)
    
    
    
for i in range(1,len(t)): 
    if (i == 1 or i%gap == 0):
        getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
        u = copyarrayfromC(u_c,n)
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        s = wdir +  "out" + str(i)+".txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
            writefile2.writerow(['dx' ,'dt','time',"cell midpoint" ,'height(m)', 'G' , 'u(m/s)','bed' ])        
                       
            for j in range(n):
                writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j])])
        
    evolvewrap(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,b0_c,b1_c,g,dx,dt,nBC,n,nBCs,theta)
    print (t[i])
        
getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
u = copyarrayfromC(u_c,n)
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)
s = wdir +  "outlast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     writefile2.writerow(['dx' ,'dt','time',"cell midpoint" ,'height(m)', 'G' , 'u(m/s)','bed' ])        
                   
     for j in range(n):
         writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j])])
deallocPy(u_c)   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(h0_c)
deallocPy(h1_c)
deallocPy(u0_c)
deallocPy(u1_c)  


"""
##### DAMBREAK ACCURACY
g = 9.81
Cr = 0.2
hf = 1.8
hl = 1.0
wdir = "../../data/dbh/o2/"
for k in range(16,20):
    dx = 100.0 / (2**k)
    #l = Cr / sqrt(g*hf)
    l = 0.01
    dt = l*dx
    theta = 1.2
    startx = 0.0
    endx = 1000.0 + dx
    startt = 0.0
    endt = 30 + dt
    
    wdatadir = wdir+ str(k) + "/"
    if not os.path.exists(wdatadir):
        os.makedirs(wdatadir)
        
    g = 9.81
        
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
        
    bot = 0.0
    gap = max(1,int(0.5/dt))
        
    h,G,bed = dambreak(x,hf,500,hl,bot,dx)
        
    nBC = 3
    nBCs = 4
    b0 = bot*ones(nBCs)
    b1 = bot*ones(nBCs)
    u0 = zeros(nBCs)
    u1 = zeros(nBCs)    
    h0 = hf*ones(nBCs)
    h1 = hl*ones(nBCs)
        
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    bed_c = copyarraytoC(bed)
    h0_c  = copyarraytoC(h0)
    h1_c  = copyarraytoC(h1)
    u0_c  = copyarraytoC(u0)
    u1_c  = copyarraytoC(u1)
    b0_c  = copyarraytoC(b0)
    b1_c  = copyarraytoC(b1)
    u_c = mallocPy(n)
        
        
        
    for i in range(1,len(t)): 
        if (i == 1 or i%gap == 0):
            getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
            u = copyarrayfromC(u_c,n)
            G = copyarrayfromC(G_c,n)
            h = copyarrayfromC(h_c,n)
            s = wdatadir +  "out" + str(i)+".txt"
            with open(s,'a') as file2:
                writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
                writefile2.writerow(['dx' ,'dt','time',"cell midpoint" ,'height(m)', 'G' , 'u(m/s)','bed' ])        
                           
                for j in range(n):
                    writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j])])
            
        evolvewrap(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,b0_c,b1_c,g,dx,dt,nBC,n,nBCs,theta)
        print (t[i])
            
    getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
    u = copyarrayfromC(u_c,n)
    G = copyarrayfromC(G_c,n)
    h = copyarrayfromC(h_c,n)
    s = wdatadir +  "outlast.txt"
    with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
         writefile2.writerow(['dx' ,'dt','time',"cell midpoint" ,'height(m)', 'G' , 'u(m/s)','bed' ])        
                       
         for j in range(n):
             writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j])])
             
    s = wdir +  str(k)+ ".txt"
    with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
         writefile2.writerow(['dx' ,'dt','time',"cell midpoint" ,'height(m)', 'G' , 'u(m/s)','bed' ])        
                   
         for j in range(n):
             writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j])])
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
wdir = "../../../data/solcononesec/o2/"

if not os.path.exists(wdir):
    os.makedirs(wdir)

s = wdir + "savenorms.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(['dx','Normalised L1-norm Difference Height', ' Normalised L1-norm Difference Velocity'])
    
for k in range(20):
    dx = 100.0 / (2**k)
    l = 0.01
    dt = l*dx
    startx = -500.0
    endx = 1500.0 + dx
    startt = 0
    endt = 1 + dt
    theta = 1.2
    
    wdatadir = wdir+ str(k) + "/"
    if not os.path.exists(wdatadir):
        os.makedirs(wdatadir)
    
    g = 9.81
    
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
    
    a0 = 10.0
    a1 = 1.0
    t0 = 0
    bot = 0
    gap = max(1,int(0.5/dt))
    
    h,G,bed = solitoninit(n,a0,a1,g,x,t0,bot,dx)
    
    nBC = 3
    nBCs = 4
    b0 = bot*ones(nBCs)
    b1 = bot*ones(nBCs)
    u0 = zeros(nBCs)
    u1 = zeros(nBCs)    
    h0 = h[0]*ones(nBCs)
    h1 = h[-1]*ones(nBCs)
    
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    bed_c = copyarraytoC(bed)
    h0_c  = copyarraytoC(h0)
    h1_c  = copyarraytoC(h1)
    u0_c  = copyarraytoC(u0)
    u1_c  = copyarraytoC(u1)
    b0_c  = copyarraytoC(b0)
    b1_c  = copyarraytoC(b1)
    u_c = mallocPy(n)
    
    
    
    for i in range(1,len(t)):
        
        if(i % gap == 0 or i ==1):
            getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
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
            
            print t[i]
            print(h[1],G[1]) 
            with open(s,'a') as file2:
                writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
                writefile2.writerow(['dx' ,'dt','time', 'height(m)', 'G' , 'u(m/s)','bed','true height', 'true velocity' ])        
                   
                for j in range(n):
                    writefile2.writerow([str(dx),str(dt),str(t[i]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j]) , str(htrue[j]), str(utrue[j])])  
                 
            
        evolvewrap(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,b0_c,b1_c,g,dx,dt,nBC,n,nBCs,theta)
        
    getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
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
    
         writefile2.writerow(['dx' ,'dt','time', 'height(m)', 'G' , 'u(m/s)','bed','true height', 'true velocity'  ])        
                   
         for j in range(n):
             writefile2.writerow([str(dx),str(dt),str(t[-1]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j]), str(htrue[j]), str(utrue[j])])       
    
    normhdiffi = norm(h - htrue,ord=1) / norm(htrue,ord=1)
    normudiffi = norm(u -utrue,ord=1) / norm(utrue,ord=1)  

    s = wdir + "savenorms.txt"
    with open(s,'a') as file1:
        writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

        writefile.writerow([str(dx),str(normhdiffi), str(normudiffi)]) 

    deallocPy(u_c)   
    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(bed_c)
    deallocPy(h0_c)
    deallocPy(h1_c)
    deallocPy(u0_c)
    deallocPy(u1_c)
    deallocPy(b0_c)
    deallocPy(b1_c)    
"""
"""
##### SOLITON
dx = 0.1
l = 0.01
dt = l*dx
startx = -500.0
endx = 1500.0 + dx
startt = 0
endt = 100.0 + dt

theta = 1.2
    
wdir = "../../data/Cserre/soliton/order2/t/"
    
g = 9.81
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
    
a0 = 10.0
a1 = 1.0
t0 = 0
bot = 0
gap = max(1,int(0.5/dt))
    
h,G,bed = solitoninit(n,a0,a1,g,x,t0,bot,dx)
    
nBC = 3
nBCs = 4
b0 = bot*ones(nBCs)
b1 = bot*ones(nBCs)
u0 = zeros(nBCs)
u1 = zeros(nBCs)    
h0 = h[0]*ones(nBCs)
h1 = h[-1]*ones(nBCs)
    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
bed_c = copyarraytoC(bed)
h0_c  = copyarraytoC(h0)
h1_c  = copyarraytoC(h1)
u0_c  = copyarraytoC(u0)
u1_c  = copyarraytoC(u1)
b0_c  = copyarraytoC(b0)
b1_c  = copyarraytoC(b1)
u_c = mallocPy(n)
    
    
    
for i in range(1,len(t)):
       
    if(i % gap == 0 or i ==1):
       getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
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
    
            writefile2.writerow(['dx' ,'dt','time', 'height(m)', 'G' , 'u(m/s)','bed','true height', 'true velocity' ])        
                   
            for j in range(n):
                writefile2.writerow([str(dx),str(dt),str(t[i]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j]) , str(htrue[j]), str(utrue[j])])  
                 
            
    evolvewrap(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,b0_c,b1_c,g,dx,dt,nBC,n,nBCs,theta)
        
getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
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
     writefile2.writerow(['dx' ,'dt','time', 'height(m)', 'G' , 'u(m/s)','bed','true height', 'true velocity'  ])        
                   
     for j in range(n):
         writefile2.writerow([str(dx),str(dt),str(t[-1]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j]), str(htrue[j]), str(utrue[j])])       
    

deallocPy(u_c)   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(bed_c)
deallocPy(h0_c)
deallocPy(h1_c)
deallocPy(u0_c)
deallocPy(u1_c)
deallocPy(b0_c)
deallocPy(b1_c)  
"""
"""
from time import time
import os
from scipy.linalg import norm
dx = 100.0 / (2.0**6)
l = 0.01
dt = l*dx
startx = -500.0
endx = 1500.0 + dx
startt = 0
endt = 100.0 + dt

theta = 1.2
    
wdir = "../../data/timecomp2/o2/"
if not os.path.exists(wdir):
    os.makedirs(wdir) 
    
g = 9.81
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
    
a0 = 10.0
a1 = 1.0
t0 = 0
bot = 0
gap = max(1,int(0.5/dt))
    
h,G,bed = solitoninit(n,a0,a1,g,x,t0,bot,dx)
    
nBC = 3
nBCs = 4
b0 = bot*ones(nBCs)
b1 = bot*ones(nBCs)
u0 = zeros(nBCs)
u1 = zeros(nBCs)    
h0 = h[0]*ones(nBCs)
h1 = h[-1]*ones(nBCs)
    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
bed_c = copyarraytoC(bed)
h0_c  = copyarraytoC(h0)
h1_c  = copyarraytoC(h1)
u0_c  = copyarraytoC(u0)
u1_c  = copyarraytoC(u1)
b0_c  = copyarraytoC(b0)
b1_c  = copyarraytoC(b1)
u_c = mallocPy(n)
    
    
t0 = time()    
for i in range(1,len(t)):
    evolvewrap(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,b0_c,b1_c,g,dx,dt,nBC,n,nBCs,theta)
    print t[i]
t1 = time()
timeelapse = t1 - t0
        
getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
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
deallocPy(bed_c)
deallocPy(h0_c)
deallocPy(h1_c)
deallocPy(u0_c)
deallocPy(u1_c)
deallocPy(b0_c)
deallocPy(b1_c)  
"""
"""
################ SEGUR AND HAMMACK RECTANGULAR WAVE ###################################
tl = 60.0
b = 0.61
h0 = 0.09
h1 = 0.1
g = 9.81

dx = 0.01
Cr = 0.2
l = Cr / sqrt(g*h1)
theta = 1.2
dt = l*dx
startx = -tl
endx = tl + dx
startt = 0
endt = 50 + dt  

wdir = "../../data/segur2/o2/"
   
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
    
bot = 0.0

gap = 1
    
h,G,bed = experiment1(x,b,h0,h1,dx)

nBC = 3
nBCs = 4
b0 = bot*ones(nBCs)
b1 = bot*ones(nBCs)
u0 = zeros(nBCs)
u1 = zeros(nBCs)    
h0 = h1*ones(nBCs)
h1 = h1*ones(nBCs)
    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
bed_c = copyarraytoC(bed)
h0_c  = copyarraytoC(h0)
h1_c  = copyarraytoC(h1)
u0_c  = copyarraytoC(u0)
u1_c  = copyarraytoC(u1)
b0_c  = copyarraytoC(b0)
b1_c  = copyarraytoC(b1)
u_c = mallocPy(n)
    
    
    
for i in range(1,len(t)): 
    if (i == 1 or i%gap == 0):
        getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
        u = copyarrayfromC(u_c,n)
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        s = wdir +  "out" + str(i)+".txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
            writefile2.writerow(['dx' ,'dt','time',"cell midpoint" ,'height(m)', 'G' , 'u(m/s)','bed' ])        
                       
            for j in range(n):
                writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j])])
        
    evolvewrap(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,b0_c,b1_c,g,dx,dt,nBC,n,nBCs,theta)
    print (t[i])
        
getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
u = copyarrayfromC(u_c,n)
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)
s = wdir +  "outlast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     writefile2.writerow(['dx' ,'dt','time',"cell midpoint" ,'height(m)', 'G' , 'u(m/s)','bed' ])        
                   
     for j in range(n):
         writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j])])
deallocPy(u_c)   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(h0_c)
deallocPy(h1_c)
deallocPy(u0_c)
deallocPy(u1_c) 
"""