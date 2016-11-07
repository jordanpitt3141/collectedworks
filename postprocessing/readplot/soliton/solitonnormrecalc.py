import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

order = "o3"


wdir = "../../../../data/raw/solconnonsmallg10Hamil/o3/"
sdir = "../../../../data/postprocessing/hinonling10recalc/o3/"

if not os.path.exists(sdir):
    os.makedirs(sdir) 

gap = 1
        
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
            h.append(float(row[1]))
            u.append(float(row[2]))
            E.append(abs(float(row[3])))
                
        j = j + 1
    u = array(u)
    h = array(h)
    E = array(E)
    dxs = array(dxs)
    
         
#scaling
ldx = dxs
lh = h
lu = u
lE = E

#LE are old energy relative errors, nee dot calculate with new initial conditions

Hi = 1527.68293

ndxs = []
nEs = []

# Go through old files and read them to get energies from outputlast
for i in range(6,20):
    
    
    s = wdir+"/"+str(i)+ "/" + "saveoutputtslast.txt"
    with open(s,'r') as file1:
        readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

        j = -1
        for row in readfile:       
            if (j >= 0):
                dx= float(row[0])
                nE = float(row[4])
                    
            j = j + 1

    
    ndxs.append(dx)
    nEs.append(nE)


relE = []
for i in range(len(nEs)):
    relE.append(abs(Hi -nEs[i] ) / abs(Hi))
    

    
n = len(nEs)

s = sdir + "h.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(ldx[i]," ",lh[i])
        file1.write(s)
s = sdir + "u.dat"
with open(s,'w') as file2:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(ldx[i]," ",lu[i])
        file2.write(s)
s = sdir + "nE.dat"
with open(s,'w') as file2:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(ldx[i]," ",relE[i])
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
