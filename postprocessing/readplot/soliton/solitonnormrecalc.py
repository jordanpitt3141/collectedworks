import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

from Hamil import *
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
    
order = "o3"


wdir = "../../../../data/raw/FDreredo/grim/"
sdir = "../../../../data/postprocessing/L1/FDreredo/grim/"

if not os.path.exists(sdir):
    os.makedirs(sdir) 

gap = 1
        
s = wdir + "savenorms.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    oh = []
    ou = []
    oE = []
    odxs = []
    j = -1
    for row in readfile:       
        if (j >= 0):
            odxs.append(float(row[0]))
            oh.append(float(row[1]))
            ou.append(float(row[2]))
            oE.append(abs(float(row[3])))
                
        j = j + 1
    ou = array(ou)
    oh = array(oh)
    oE = array(oE)
    odxs = array(odxs)

#eta = 1
#Hi = 1527.68293

#eta = 0.7
Hi =  1508.917011


ndxs = []
nEs = []
relE = []
# Go through old files and read them to get energies from outputlast
for i in range(6,20):
    
    
    s = wdir+"/"+str(i)+ "/" + "outlast.txt"
    with open(s,'r') as file1:    
        readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
        h = []
        u = []
        ht = []
        ut = []
        x = []
        j = -1
        for row in readfile:       
            if (j >= 0):
                dx =float(row[0])
                dt =float(row[1])
                t =float(row[2])
                x.append(float(row[6]))
                h.append(float(row[7]))
                u.append(float(row[8]))
                   
            j = j + 1
    startx = -50 + dx
    endx = 250
    a1 = 0.7
    a0 = 1.0
    g = 9.81
    #g= 1.0
    t0 = 0.0
    x = x[1:]
    n = len(x) 
    ht,ut = solitoninit(n,a0,a1,g,x,t,dx)   
    hti,uti = solitoninit(n,a0,a1,g,x,t0,dx) 
    
    niBC = 3    
    xbeg = arange(startx - niBC*dx,startx,dx)
    xend = arange(endx + dx,endx + (niBC+1)*dx,dx) 
    hbeg = h[0]*ones(niBC)
    hend = h[-1]*ones(niBC)
    ubeg = u[0]*ones(niBC)
    uend = u[-1]*ones(niBC)   
    
    htbeg = ht[0]*ones(niBC)
    htend = ht[-1]*ones(niBC)
    utbeg = ut[0]*ones(niBC)
    utend = ut[-1]*ones(niBC)  
    
    htibeg = hti[0]*ones(niBC)
    htiend = hti[-1]*ones(niBC)
    utibeg = uti[0]*ones(niBC)
    utiend = uti[-1]*ones(niBC)   
            
    xbc =  concatenate([xbeg,array(x),xend])
    hbc =  concatenate([hbeg,array(h),hend])
    ubc =  concatenate([ubeg,array(u),uend])
    
    htbc =  concatenate([htbeg,array(ht),htend])
    utbc =  concatenate([utbeg,array(ut),utend])
    
    htibc =  concatenate([htibeg,array(hti),htiend])
    utibc =  concatenate([utibeg,array(uti),utiend])
    
    xbc_c = copyarraytoC(xbc)
    hbc_c = copyarraytoC(hbc)
    ubc_c = copyarraytoC(ubc)
    htbc_c = copyarraytoC(htbc)
    utbc_c = copyarraytoC(utbc)
    htibc_c = copyarraytoC(htibc)
    utibc_c = copyarraytoC(utibc)
    
    Eval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)    

    #Evaltrel.append(abs(Evalt - Eval)/abs(Evalt))
    #Evaltabs.append(abs(Evalt - Eval))
    relE.append(abs(Hi - Eval)/abs(Hi))
    ndxs.append(dx)


    
n = len(relE)

s = sdir + "h.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(odxs[i]," ",oh[i])
        file1.write(s)
s = sdir + "u.dat"
with open(s,'w') as file2:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(odxs[i]," ",ou[i])
        file2.write(s)
s = sdir + "nE.dat"
with open(s,'w') as file2:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(ndxs[i]," ",relE[i])
        file2.write(s)
  

#Energy Only

"""  
s = wdir + "savenorms.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    h = []
    u = []
    E = []
    dxs = []
    j = -1
    for row in readfile:       
        if (j >= 0):
            dxs.append(float(row[0]))
            E.append(abs(float(row[1])))
                
        j = j + 1
    E = array(E)
    dxs = array(dxs)
    
         
#scaling
ldx = dxs
lE = E


    
n = len(lE)
s = sdir + "E.dat"
with open(s,'w') as file2:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(ldx[i]," ",lE[i])
        file2.write(s)
"""    
