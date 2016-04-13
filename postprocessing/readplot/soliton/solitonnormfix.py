import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

order = "3"
wdir = "../../../../data/raw/Chrisnew/solcon/o"+order+"/"
sdir = "../../../../data/postprocessing/Chrisnew/solcon/o"+order+"/"

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
            E.append(float(row[3]))
                
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

    
n = len(lh)

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
s = sdir + "E.dat"
with open(s,'w') as file2:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(ldx[i]," ",lE[i])
        file2.write(s)
        
s = sdir + "ch.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(ldx[i]," ",lh[i]/(ldx[i]**int(order)))
        file1.write(s)
s = sdir + "cu.dat"
with open(s,'w') as file2:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(ldx[i]," ",lu[i]/(ldx[i]**int(order)))
        file2.write(s)
s = sdir + "cE.dat"
with open(s,'w') as file2:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(ldx[i]," ",lE[i]/(ldx[i]**int(order)))
        file2.write(s)
