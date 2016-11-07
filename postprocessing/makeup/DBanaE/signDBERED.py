
import csv
from numpy.linalg import norm
from scipy import *
import os
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
order = "o3"
wdir = "../../../../data/postprocessing/dbEnergyRED/"+order+"/"
sdir = "../../../../data/postprocessing/DBEC/"+order+"/"


def HamilDB(alpha,dx):
    return 10.3986*(1000 + dx) - 0.7848*((2.0/alpha)*tanh(alpha * (500.0+ 0.5*dx)))

if not os.path.exists(sdir):
    os.makedirs(sdir)

         
s = wdir + "Evals.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    Evalnumis = []
    Evalnumfs = []
    betas = []
    dxs = []
    relEvalnums = []
    j = -1
    for row in readfile:       
        if (j >= 0):
            dxs.append(float(row[0]))
            betas.append(float(row[1]))
            Evalnumfs.append(float(row[2]))
                
        j = j + 1
    Evalnumfs = array(Evalnumfs)
    betas = array(betas)
    dxs = array(dxs)

Evalis = []    
n = len(dxs)
for i in range(n):
    Evalis.append(HamilDB(betas[i],dxs[i]))
    
relE = []    
for i in range(n):
   relE.append(abs(Evalis[i] - Evalnumfs[i])/abs(Evalis[i]))
   
# want to do a couple things, break up into a to inf and x to 0
   
#dx fixed   

ldxs = dxs.tolist()
indlistdx = []
for i in range(3,13):
    num = 10.0/(2**i)
    indlistdx.append(ldxs.index(num))

iln = len(indlistdx) 
for i in range(iln):
    if i == 0:
     sdirn = sdir +"/dxn" + str(i)+ "/" 
     if not os.path.exists(sdirn):
         os.makedirs(sdirn) 
     s = sdirn + "outlast.txt"   
     with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
         writefile2.writerow(['dx' ,'betas','H inital','H num end', 'Relative Error'])        

         for k in range(indlistdx[i],n):
            writefile2.writerow([str(dxs[k]),str(betas[k]),str(Evalis[k]),str(Evalnumfs[k]), str(relE[k])])
     file2.close()
     
     s = sdirn + "relE.dat"     
     with open(s,'a') as file3:
         for k in range(indlistdx[i],n):
             s ="%3.8f%5s%1.15f\n" %(betas[k]," ",relE[k])
             file3.write(s)
        
    else:
     sdirn = sdir +"/dxn" + str(i)+ "/" 
     if not os.path.exists(sdirn):
         os.makedirs(sdirn) 
     s = sdirn + "outlast.txt"   
     with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
         writefile2.writerow(['dx' ,'betas','H inital','H num end', 'Relative Error'])           

         for k in range(indlistdx[i],indlistdx[i-1]):
            writefile2.writerow([str(dxs[k]),str(betas[k]),str(Evalis[k]),str(Evalnumfs[k]), str(relE[k])])
     file2.close()
     
     s = sdirn + "relE.dat"     
     with open(s,'a') as file3:
         for k in range(indlistdx[i],indlistdx[i-1]):
             s ="%3.8f%5s%1.15f\n" %(betas[k]," ",relE[k])
             file3.write(s)

           
# fix beta, first we have to re order everything
           
nbetas = [x for (x,y,z,a,b) in sorted(zip(betas,dxs,Evalis,Evalnumfs,relE), key=lambda pair: pair[0])] 
ndxs = [y for (x,y,z,a,b) in sorted(zip(betas,dxs,Evalis,Evalnumfs,relE), key=lambda pair: pair[0])] 
nEvalis = [z for (x,y,z,a,b) in sorted(zip(betas,dxs,Evalis,Evalnumfs,relE), key=lambda pair: pair[0])] 
nEvalnumfs = [a for (x,y,z,a,b) in sorted(zip(betas,dxs,Evalis,Evalnumfs,relE), key=lambda pair: pair[0])] 
nrelE = [b for (x,y,z,a,b) in sorted(zip(betas,dxs,Evalis,Evalnumfs,relE), key=lambda pair: pair[0])] 

nbetas = array(nbetas)
ndxs = array(ndxs)
nEvalis = array(nEvalis)
nEvalnumfs = array(nEvalnumfs )
nrelE = array(nrelE)


diffs = [0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1.0,2.5,5.0,7.5,10.0,25.0,50.0,75.0,100.0,250.0,500.0,750.0,1000.0]
nd = len(diffs)
lbetas = nbetas.tolist()
indlistbeta = []
for i in range(nd):
    num = diffs[i]
    indlistbeta.append(lbetas.index(num))

ilbn = len(indlistbeta) 
            
for i in range(ilbn):
    if i == ilbn - 1:
     sdirn = sdir +"/betan" + str(diffs[i])+ "/" 
     if not os.path.exists(sdirn):
         os.makedirs(sdirn) 
     s = sdirn + "outlast.txt"   
     with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
         writefile2.writerow(['betas','dx','H inital','H num end', 'Relative Error'])        

         for k in range(indlistbeta[i],n):
            writefile2.writerow([str(nbetas[k]),str(ndxs[k]),str(nEvalis[k]),str(nEvalnumfs[k]), str(nrelE[k])])
     file2.close()
     
     s = sdirn + "relE.dat"     
     with open(s,'a') as file3:
         for k in range(indlistbeta[i],n):
             s ="%3.8f%5s%1.15f\n" %(ndxs[k]," ",nrelE[k])
             file3.write(s)
        
    else:
     sdirn = sdir +"/betan" + str(diffs[i])+ "/" 
     if not os.path.exists(sdirn):
         os.makedirs(sdirn) 
     s = sdirn + "outlast.txt"   
     with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
         writefile2.writerow(['betas','dx','H inital','H num inital','H num end', 'Relative Error'])           

         for k in range(indlistbeta[i],indlistbeta[i+1]):
            writefile2.writerow([str(nbetas[k]),str(ndxs[k]),str(nEvalis[k]),str(nEvalnumfs[k]), str(nrelE[k])])
     file2.close()
     
     s = sdirn + "relE.dat"     
     with open(s,'a') as file3:
         for k in range(indlistbeta[i],indlistbeta[i+1]):
             s ="%3.8f%5s%1.15f\n" %(ndxs[k]," ",nrelE[k])
             file3.write(s)
    