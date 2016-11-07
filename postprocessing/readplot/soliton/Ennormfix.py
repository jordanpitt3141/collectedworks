import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

order = "o3"
#wdir = "../../../../data/raw/Chrisnew/solcon/o"+order+"/"
#sdir = "../../../../data/postprocessing/Chrisnew/solcon/o"+order+"/"

#wdir = "../../../../data/raw/ChrisDatn/solen/"+order+"/"
#sdir = "../../../../data/postprocessing/ChrisDatn/sendsolen/"+order+"/"

wdir = "../../../../data/raw/solbumpn/"
sdir = "../../../../data/postprocessing/solbumpn/"

if not os.path.exists(sdir):
    os.makedirs(sdir) 

gap = 1
        
s = wdir + "savenorms.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
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
