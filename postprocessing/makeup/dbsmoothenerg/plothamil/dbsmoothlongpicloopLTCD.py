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
import os

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
        h[i] = base + 0.5*eta0*(1 + tanh(diffuse*(x0 - x[i])))

    return h,u

wdirord = "o3"    

wdir = "../../../../../data/raw/longcontactdiscdx9diff10fileio/o3/9/0/"
sdir = "../../../../../data/postprocessing/RLHamilInit/"
if not os.path.exists(sdir):
    os.makedirs(sdir) 

for ts in range(0,3005,5):  
    
    sdirf = sdir + str(float(ts/10.0)) + "secs/"
    if not os.path.exists(sdirf):
        os.makedirs(sdirf) 
    
    s = wdir +"out"+ str(int(ts)) + ".txt"
    with open(s,'r') as file1:
         readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
         h = []
         u = []
         x = []
         j = -1
         for row in readfile:       
             if (j >= 0):
                dx = float(row[0])
                dt = float(row[1])
                t = float(row[2])
                x.append(float(row[3]))
                h.append(float(row[4]))
                u.append(float(row[6])) 
                beta = 10.0
                
             j = j + 1
         x = array(x)
         h = array(h)     
      
    startx = x[0]
    endx = x[-1]
    g = 9.81
    n = len(x)
    niBC = 3    
    xbeg = arange(startx - niBC*dx,startx,dx)
    xend = arange(endx + dx,endx + (niBC+1)*dx) 
    hbeg = h[0]*ones(niBC)
    hend = h[-1]*ones(niBC)
    ubeg = u[0]*ones(niBC)
    uend = u[-1]*ones(niBC)
    xbc =  concatenate([xbeg,array(x),xend])
    hbc =  concatenate([hbeg,array(h),hend])
    ubc =  concatenate([ubeg,array(u),uend])
    
    Ham_c = mallocPy(n)
    HamFT_c = mallocPy(n)
    HamST_c = mallocPy(n)
    HamTT_c = mallocPy(n)
    
    
    xbc_c = copyarraytoC(xbc)
    hbc_c = copyarraytoC(hbc)
    ubc_c = copyarraytoC(ubc)
    
    HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx,Ham_c)
    HankEnergyallPT(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx,HamFT_c,HamST_c,HamTT_c)
    Ham = copyarrayfromC(Ham_c,n)
    HamFT = copyarrayfromC(HamFT_c,n)
    HamST = copyarrayfromC(HamST_c,n)
    HamTT = copyarrayfromC(HamTT_c,n)
    
    hi,ui = dambreaksmooth(x,500,1.0,0.8,beta,dx)
    
    hibc =  concatenate([hbeg,array(hi),hend])
    uibc =  concatenate([ubeg,array(ui),uend])
    
    hibc_c = copyarraytoC(hibc)
    uibc_c = copyarraytoC(uibc)
    Hami_c = mallocPy(n)
    HamiFT_c = mallocPy(n)
    HamiST_c = mallocPy(n)
    HamiTT_c = mallocPy(n)
    
    HankEnergyall(xbc_c,hibc_c,uibc_c,g,n + 2*niBC,niBC,dx,Hami_c)
    HankEnergyallPT(xbc_c,hibc_c,uibc_c,g,n + 2*niBC,niBC,dx,HamiFT_c,HamiST_c,HamiTT_c)
    Hami = copyarrayfromC(Hami_c,n)
    HamiFT = copyarrayfromC(HamiFT_c,n)
    HamiST = copyarrayfromC(HamiST_c,n)
    HamiTT = copyarrayfromC(HamiTT_c,n)
    
    #Write all terms to files
    s = sdirf + "allout.txt"
    with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
         writefile2.writerow(['dx' ,'dt','time','beta','cell midpoint', 'FT', 'ST' , 'TT',"FTi","STi","TTi"])        
                       
         for j in range(n):
             writefile2.writerow([str(dx),str(dt),str(t),str(beta), str(x[j]), str(HamFT[j]) , str(HamST[j]) , str(HamTT[j]), str(HamiFT[j]), str(HamiST[j]), str(HamiTT[j])])   
    
    
    #Write the sums and percentages to files
    
    Ham = sum(HamFT) + sum(HamST) + sum(HamTT)
    Hami = sum(HamiFT) + sum(HamiST) + sum(HamiTT)
    
    relFT = (1.0*sum(HamFT)) / Hami
    relST = (1.0*sum(HamST)) / Hami
    relTT = (1.0*sum(HamTT)) / Hami
    
    reliFT = (1.0*sum(HamiFT)) / Hami
    reliST = (1.0*sum(HamiST)) / Hami
    reliTT = (1.0*sum(HamiTT)) / Hami
    
    s = sdirf + "relperc.txt"
    with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
         writefile2.writerow(['dx' ,'dt','time','beta','Ham', 'FT', 'ST' , 'TT','Hami',"FTi","STi","TTi"])        
                       
         writefile2.writerow([str(dx),str(dt),str(t),str(beta), str(Ham), str(relFT) , str(relST) , str(relTT), str(Hami), str(reliFT), str(reliST), str(reliTT)])   
    
    deallocPy(Ham_c)
    deallocPy(HamFT_c)
    deallocPy(HamST_c)
    deallocPy(HamTT_c)
    
    deallocPy(xbc_c)
    deallocPy(hbc_c)
    deallocPy(ubc_c)
    
    deallocPy(hibc_c)
    deallocPy(uibc_c)
    
    deallocPy(Hami_c)
    deallocPy(HamiFT_c)
    deallocPy(HamiST_c)
    deallocPy(HamiTT_c)       
        
plot(x,HamiFT,"--",label="FTi")
plot(x,HamiST,"--",label="STi")
plot(x,HamiTT,"--",label="TTi")

plot(x,HamFT,label="FT")
plot(x,HamST,label="ST")
plot(x,HamTT,label="TT")

legend()

        
"""
s = sdir1 + "Evalsfinal.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow([str(dx) ,str(beta),str(Eval)]) 
"""
        
        