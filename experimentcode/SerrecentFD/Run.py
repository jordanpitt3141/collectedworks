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
        h[i] = base + 0.5*eta0*(1 + tanh(diffuse*(x0 -x[i])))

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
wdir = "../../../data/raw/FDreredo/FDcent/"

if not os.path.exists(wdir):
    os.makedirs(wdir) 
    
s = wdir + "savenorms.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(['dx','Normalised L1-norm Difference Height', ' Normalised L1-norm Difference Velocity', 'Eval Error'])

for k in range(6,21):
    dx = 100.0 / (2**k)
    Cr = 0.5
    
    g = 9.81
    a0 = 1.0
    a1 = 0.7
    
    l = 1.0 / (sqrt(g*(a0 + a1)))
    dt = Cr*l*dx
    startx = -50.0
    endx = 250.0 + dx
    startt = 0
    endt = 50 + dt
        
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
        
    
    gap = max(5,int(0.5/dt))
    
    
    h,u = solitoninit(n,a0,a1,g,x,0.0,dx)
    ph,pu = solitoninit(n,a0,a1,g,x,-dt,dx)
        
    nBC = 3
    niBC = nBC
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
    
    xbeg = arange(startx - niBC*dx,startx,dx)
    xend = arange(endx + dx,endx + (niBC+1)*dx) 
    
    xbc =  concatenate([xbeg,x,xend])  
    
    xbc_c = copyarraytoC(xbc)
    hbc_c = mallocPy(n + 2*niBC)
    ubc_c = mallocPy(n + 2*niBC)
    Evals = []
   
    conc(h0_c , h_c,h1_c,niBC,n ,niBC , hbc_c)
    conc(u0_c , u_c,u1_c,niBC,n ,niBC , ubc_c)         
    Evali = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)     
          
    for i in range(1,len(t)):            
        evolvewrap(u_c, h_c, pubc_c,phbc_c, h0_c, h1_c,u0_c, u1_c,g,dx,dt,nBC, n,nBCs)    
        print (t[i])
    
    conc(h0_c , h_c,h1_c,niBC,n ,niBC , hbc_c)
    conc(u0_c , u_c,u1_c,niBC,n ,niBC , ubc_c)         
    Evalf = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)
    u = copyarrayfromC(u_c,n)
    h = copyarrayfromC(h_c,n)  
    he,ue = solitoninit(n,a0,a1,g,x,t[i],dx)
    
    if not os.path.exists(wdir+ str(k) + "/"):
        os.makedirs(wdir+ str(k) + "/") 
    
    s = wdir+ str(k) + "/"  + "outlast.txt"
    with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
         writefile2.writerow(['dx' ,'dt','time','Evali','Evalf','Eval Error','cell midpoint', 'height(m)', 'u(m/s)',"he","ue"])        
                           
         for j in range(n):
             writefile2.writerow([str(dx),str(dt),str(t[i]), str(Evali), str(Evalf), str(abs(Evali - Evalf)/ abs(Evali)), str(x[j]), str(h[j]) , str(u[j]), str(he[j]),str(ue[j])]) 
             
    normhdiffi = norm(h - he,ord=1) / norm(he,ord=1)
    normudiffi = norm(u -ue,ord=1) / norm(ue,ord=1)  
    
    s = wdir + "savenorms.txt"
    with open(s,'a') as file1:
        writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

        writefile.writerow([str(dx),str(normhdiffi), str(normudiffi), str(abs(Evali - Evalf)/ abs(Evali))])
    
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
##Soliton Collision
wdir = "../../../data/raw/Cserre/solitonothers/collDMcopyhh/FDc/"
if not os.path.exists(wdir):
    os.makedirs(wdir)
    
dx = 0.01

a0 = 1.0
a11 = 0.96
solbeg1 = 100.0
solend1 = 200.0
direction1 = 1.0
a12 = 0.96
solbeg2 = 200.0
solend2 = 300.0
direction2 = -1.0

Cr = 0.5
#g = 9.81
g = 1.0
#l = Cr / (sqrt(g*1.5*(a0 + a11 + a12)))
dt = 0.1*dx
startx = 0.0
endx = 400.0
startt = 0.0
endt = 150 + dt

    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
    
t0 = 0
gap = int(0.5/dt)


h,u = soliton2interactinit(n,a0,a11,solbeg1,solend1,direction1,a12,solbeg2,solend2,direction2,g,x,t0,dx)
ph,pu = soliton2interactinit(n,a0,a11,solbeg1,solend1,direction1,a12,solbeg2,solend2,direction2,g,x,-dt,dx)
    
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
        s = wdir + "saveoutputts" + str(i) + ".txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
            writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'u(m/s)'])        
                       
            for j in range(n):
                writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(u[j])])  
        
    evolvewrap(u_c, h_c, pubc_c,phbc_c , h0_c, h1_c,u0_c, u1_c,g,dx,dt,nBC, n,nBCs)    
    print (t[i])

u = copyarrayfromC(u_c,n)
h = copyarrayfromC(h_c,n)  
s = wdir + "saveoutputtslast.txt"
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
"""
## DAM BREAK Smooth ##########################
wdir = "../../data/raw/longtimedambreakDTfix/FDc/"
if not os.path.exists(wdir):
    os.makedirs(wdir) 
dx = 10.0 /(2**9)
l = 0.01
dt = l*dx*dx
startx = -900
endx = 1800.0 + dx
startt = 0.0
endt = 300 + dt
    
g = 9.81
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
    
bot = 0.0
hf = 1.8
hl = 1.0
gap = int(0.5/dt)

diffuse = 10
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
   
aplus = []
aplusx = []
aplust = []   
      
for i in range(1,len(t)):
    if(i % gap == 0 or i ==1):
        u = copyarrayfromC(u_c,n)
        h = copyarrayfromC(h_c,n)
        
        mi = n - 2
        for mi in range(n-1,-1,-1):
            if(h[mi -1] < h[mi]) and (h[mi] > 1.1 ):
                break
        aplus.append(h[mi])
        aplusx.append(x[mi])
        aplust.append(t[i])
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
        
     writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'u(m/s)'])        
                       
     for j in range(n):
         writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(u[j])])     

s = wdir + "aplus.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(['x' ,'t','aplus' ,"Grim"])        
           
    for j in range(len(aplus)):
        writefile2.writerow([str(aplusx[j]),str(aplust[j]),str(aplus[j]),str(0.739976603390100695296990254)]) 

deallocPy(u_c)   
deallocPy(h_c)
deallocPy(h0_c)
deallocPy(h1_c)
deallocPy(u0_c)
deallocPy(u1_c)
"""

"""
#big smooth NEW


#TEST

#diffuses = [0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1.0,2.5,5.0,7.5,10.0,25.0,50.0,75.0,100.0,250.0,500.0,750.0,1000.0]
diffuses = [2]
wdirb = "../../data/bigsmoothTEST/FDc/"
for ll in range(9,10):
    for k in range(len(diffuses)):
        wdir = wdirb + str(ll) + "/" + str(k) + "/"
        if not os.path.exists(wdir):
            os.makedirs(wdir) 
        dx = 10.0 / (2**ll)
        Cr = 0.5
        g = 9.81
        hf = 1.3
        l = 1.0 / sqrt(g*hf)
        dt = l*dx
        startx = 0.0
        endx = 1200.0 + dx
        startt = 0.0
        endt = 150.0+(dt*0.9) 
        
        
        x,t = makevar(startx,endx,dx,startt,endt,dt)
        n = len(x)
            
        bot = 0.0
        hl = 1.0
        gap = max(1,int(0.02/dt))
        
        diffuse = diffuses[k]
        base = hl
        eta0 = hf - hl
        x0 = 600
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

difflist = [12]

deltaxa = [12,13,14]
dxlist = [deltaxa,deltaxa,deltaxa,deltaxa,deltaxa]

diffuses = [0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1.0,2.5,5.0,7.5,10.0,25.0,50.0,75.0,100.0,250.0,500.0,750.0,1000.0]
wdirb = "../../data/bigsmoothtargettedNEW1/FDcent/"
for lk in range(len(difflist)):
    for ll in dxlist[lk]:
        wdir = wdirb + str(ll) + "/" + str(difflist[lk]) + "/"
        if not os.path.exists(wdir):
            os.makedirs(wdir) 
        dx = (10.0 / (2**ll))
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
