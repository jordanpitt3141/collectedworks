# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:09:21 2015

@author: jordan
"""
from Serre2FDC import *
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
        

def dambreak(x,hf,hc,hl,dx):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    
    for i in range(n):
        if (x[i] < hc):
            h[i] = hf
        else:
            h[i] = hl
    return h,u 

def dambreaksmooth(x,x0,base,eta0,diffuse,dx):
    from numpy import tanh
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    
    for i in range(n):
        h[i] = base + 0.5*eta0*(1 + tanh(diffuse*(x0 - abs(x[i]))))

    return h,u       

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
    
    return h,u
    
def experiment1(x,h0,h1,dx):
    n = len(x)
    u = zeros(n)
    h = ones(n)*h1
    for i in range(n):
        if (x[i] <0 and x[i] > -2*b):
            h[i] = h0

    return h,u
"""
#solitonaccuracy
wdir = "../../../data/solcononesec/FDcent/"

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
    startt = 0.0
    endt = 1 + dt  
        
    g = 9.81
        
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
        
    
    gap = max(5,int(0.5/dt))
    
    a0 = 10.0
    a1 = 1.0
    
    
    h,u = solitoninit(n,a0,a1,g,x,0.0,dx)
    ph,pu = solitoninit(n,a0,a1,g,x,-dt,dx)
        
    nBC = 3
    nBCs = 4
    u0 = zeros(nBCs)
    u1 = zeros(nBCs)    
    h0 = a0*ones(nBCs)
    h1 = a0*ones(nBCs)
        
    h_c = copyarraytoC(h)
    u_c = copyarraytoC(u)
    pubc_c = copyarraytoC(concatenate([u0[-nBC:],pu,u1[:nBC]]))
    phbc_c = copyarraytoC(concatenate([h0[-nBC:],ph,h1[:nBC]]))
    h0_c  = copyarraytoC(h0)
    h1_c  = copyarraytoC(h1)
    u0_c  = copyarraytoC(u0)
    u1_c  = copyarraytoC(u1)
        
          
    for i in range(1,len(t)):            
        evolvewrap(u_c, h_c, pubc_c,phbc_c, h0_c, h1_c,u0_c, u1_c,g,dx,dt,nBC, n,nBCs)    
        print (t[i])
    
    u = copyarrayfromC(u_c,n)
    h = copyarrayfromC(h_c,n)  
    he,ue = solitoninit(n,a0,a1,g,x,t[i],dx)
    
    if not os.path.exists(wdir+ str(k) + "/"):
        os.makedirs(wdir+ str(k) + "/") 
    
    s = wdir+ str(k) + "/"  + "outlast.txt"
    with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
         writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'u(m/s)',"he","ue"])        
                           
         for j in range(n):
             writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(u[j]), str(he[j]),str(ue[j])]) 
             
    normhdiffi = norm(h - he,ord=1) / norm(he,ord=1)
    normudiffi = norm(u -ue,ord=1) / norm(ue,ord=1)  
    
    s = wdir + "savenorms.txt"
    with open(s,'a') as file1:
        writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

        writefile.writerow([str(dx),str(normhdiffi), str(normudiffi)])
    
    deallocPy(u_c)   
    deallocPy(h_c)
    deallocPy(h0_c)
    deallocPy(h1_c)
    deallocPy(u0_c)
    deallocPy(u1_c)
"""

"""
##Soliton
wdir = "../../data/t/"
dx = 0.1
l = 0.01
dt = l*dx
startx = -500.0
endx = 1500.0 + dx
startt = 0.0
endt = 10 + dt  
    
g = 9.81
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
    

gap = 10**100#max(5,int(0.5/dt))

a0 = 10.0
a1 = 1.0


h,u = solitoninit(n,a0,a1,g,x,0.0,dx)
ph,pu = solitoninit(n,a0,a1,g,x,-dt,dx)
    
nBC = 3
nBCs = 4
u0 = zeros(nBCs)
u1 = zeros(nBCs)    
h0 = a0*ones(nBCs)
h1 = a0*ones(nBCs)
    
h_c = copyarraytoC(h)
u_c = copyarraytoC(u)
pubc_c = copyarraytoC(concatenate([u0[-nBC:],pu,u1[:nBC]]))
phbc_c = copyarraytoC(concatenate([h0[-nBC:],ph,h1[:nBC]]))
h0_c  = copyarraytoC(h0)
h1_c  = copyarraytoC(h1)
u0_c  = copyarraytoC(u0)
u1_c  = copyarraytoC(u1)
    
      
for i in range(1,len(t)): 
    if(i % gap == 0 or i ==1):
        u = copyarrayfromC(u_c,n)
        h = copyarrayfromC(h_c,n)  
        he,ue = solitoninit(n,10.0,1.0,g,x,t[i],dx)
        s = wdir + "saveoutputts" + str(i) + ".txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
            writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'u(m/s)',"he","ue"])        
                       
            for j in range(n):
                writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(u[j]), str(he[j]),str(ue[j])])  
        
    evolvewrap(u_c, h_c, pubc_c,phbc_c , h0_c, h1_c,u0_c, u1_c,g,dx,dt,nBC, n,nBCs)    
    print (t[i])

u = copyarrayfromC(u_c,n)
h = copyarrayfromC(h_c,n)  
he,ue = solitoninit(n,a0,a1,g,x,t[i],dx)
s = wdir + "outlast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
     writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'u(m/s)',"he","ue"])        
                       
     for j in range(n):
         writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(u[j]), str(he[j]),str(ue[j])])   
   

deallocPy(u_c)   
deallocPy(h_c)
deallocPy(h0_c)
deallocPy(h1_c)
deallocPy(u0_c)
deallocPy(u1_c)
"""
"""   
## DAM BREAK Smooth ##########################
wdir = "../../data/t/"
dx = 0.1
l = 0.01
dt = l*dx
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

diffuse = 0.5
base = hl
eta0 = hf - hl
x0 = 500
h,u = dambreaksmooth(x,x0,base,eta0,diffuse,dx)   
    
nBC = 3
nBCs = 4
u0 = zeros(nBCs)
u1 = zeros(nBCs)    
h0 = hf*ones(nBCs)
h1 = hl*ones(nBCs)
    
h_c = copyarraytoC(h)
u_c = copyarraytoC(u)
pubc_c = copyarraytoC(concatenate([u0[-nBC:],u,u1[:nBC]]))
phbc_c = copyarraytoC(concatenate([u0[-nBC:],h,u1[:nBC]]))
h0_c  = copyarraytoC(h0)
h1_c  = copyarraytoC(h1)
u0_c  = copyarraytoC(u0)
u1_c  = copyarraytoC(u1)
    
      
for i in range(1,len(t)):
    if(i % gap == 0 or i ==1):
        u = copyarrayfromC(u_c,n)
        h = copyarrayfromC(h_c,n)
        s = wdir + "saveoutputts" + str(i) + ".txt"
        with open(s,'a') as file2:
             writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                
             writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'u(m/s)'])        
                               
             for j in range(n):
                 writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(u[j])])  
    evolvewrap(u_c, h_c, pubc_c,phbc_c, h0_c, h1_c,u0_c, u1_c,g,dx,dt,nBC, n,nBCs)    
    print (t[i])

u = copyarrayfromC(u_c,n)
h = copyarrayfromC(h_c,n)
s = wdir + "outlast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
     writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'u(m/s)'])        
                       
     for j in range(n):
         writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(u[j])])     

deallocPy(u_c)   
deallocPy(h_c)
deallocPy(h0_c)
deallocPy(h1_c)
deallocPy(u0_c)
deallocPy(u1_c)
"""
#big smooth
"""
diffuses = [0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1.0,2.5,5.0,7.5,10.0,25.0,50.0,75.0,100.0,250.0,500.0,750.0,1000.0]
wdirb = "../../data/bigsmooth/grim/"
for ll in range(3,16):
    for k in range(len(diffuses)):
        wdir = wdirb + str(ll) + "/" + str(k) + "/"
        if not os.path.exists(wdir):
            os.makedirs(wdir) 
        dx = 10.0 / (2**ll)
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
        gap = max(1,int(0.02/dt))
        
        diffuse = diffuses[k]
        base = hl
        eta0 = hf - hl
        x0 = 500
        h,u = dambreaksmooth(x,x0,base,eta0,diffuse,dx)   
            
        nBC = 3
        nBCs = 4
        u0 = zeros(nBCs)
        u1 = zeros(nBCs)    
        h0 = hf*ones(nBCs)
        h1 = hl*ones(nBCs)
            
        h_c = copyarraytoC(h)
        u_c = copyarraytoC(u)
        pubc_c = copyarraytoC(concatenate([u0[-nBC:],u,u1[:nBC]]))
        phbc_c = copyarraytoC(concatenate([h0[-nBC:],h,h1[:nBC]]))
        h0_c  = copyarraytoC(h0)
        h1_c  = copyarraytoC(h1)
        u0_c  = copyarraytoC(u0)
        u1_c  = copyarraytoC(u1)
            
              
        for i in range(1,len(t)):
            evolvewrap(u_c, h_c, pubc_c,phbc_c, h0_c, h1_c,u0_c, u1_c,g,dx,dt,nBC, n,nBCs)    
            print (t[i])
        
        u = copyarrayfromC(u_c,n)
        h = copyarrayfromC(h_c,n)
        s = wdir + "outlast.txt"
        with open(s,'a') as file2:
             writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                
             writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'u(m/s)', 'diffuse'])        
                               
             for j in range(n):
                 writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(u[j]), str(diffuse)])     
        
        deallocPy(u_c)   
        deallocPy(h_c)
        deallocPy(h0_c)
        deallocPy(h1_c)
        deallocPy(u0_c)
        deallocPy(u1_c)
"""

#big smooth targeted

difflist = [1,6,8,9,12]

deltaxa = [5,6,7,9,10,11,12,13,14,15]
dxlist = [deltaxa,deltaxa,deltaxa,deltaxa,deltaxa]

diffuses = [0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1.0,2.5,5.0,7.5,10.0,25.0,50.0,75.0,100.0,250.0,500.0,750.0,1000.0]
wdirb = "../../data/bigsmoothtargetted/FDcent/"
for lk in range(len(difflist)):
    for ll in dxlist[lk]:
        wdir = wdirb + str(ll) + "/" + str(difflist[lk]) + "/"
        if not os.path.exists(wdir):
            os.makedirs(wdir) 
        dx = ll*(10.0 / (2**12))
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
        gap = max(1,int(0.02/dt))
        
        diffuse = diffuses[difflist[lk]]
        base = hl
        eta0 = hf - hl
        x0 = 500
        h,u = dambreaksmooth(x,x0,base,eta0,diffuse,dx)   
            
        nBC = 3
        nBCs = 4
        u0 = zeros(nBCs)
        u1 = zeros(nBCs)    
        h0 = hf*ones(nBCs)
        h1 = hl*ones(nBCs)
            
        h_c = copyarraytoC(h)
        u_c = copyarraytoC(u)
        pubc_c = copyarraytoC(concatenate([u0[-nBC:],u,u1[:nBC]]))
        phbc_c = copyarraytoC(concatenate([h0[-nBC:],h,h1[:nBC]]))
        h0_c  = copyarraytoC(h0)
        h1_c  = copyarraytoC(h1)
        u0_c  = copyarraytoC(u0)
        u1_c  = copyarraytoC(u1)
            
              
        for i in range(1,len(t)):
            evolvewrap(u_c, h_c, pubc_c,phbc_c, h0_c, h1_c,u0_c, u1_c,g,dx,dt,nBC, n,nBCs)    
            print (t[i])
        
        u = copyarrayfromC(u_c,n)
        h = copyarrayfromC(h_c,n)
        s = wdir + "outlast.txt"
        with open(s,'a') as file2:
             writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                
             writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'u(m/s)', 'diffuse'])        
                               
             for j in range(n):
                 writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(u[j]), str(diffuse)])     
        
        deallocPy(u_c)   
        deallocPy(h_c)
        deallocPy(h0_c)
        deallocPy(h1_c)
        deallocPy(u0_c)
        deallocPy(u1_c)
