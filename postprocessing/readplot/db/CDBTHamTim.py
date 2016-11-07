import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy import ones



wdir = "../../../../data/raw/CDBThetaHAMILYTIME/o2/1.1/"
sdir = "../../../../data/postprocessing/CDBC/"


s = wdir + "Hamil.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    t = []
    H = []
    j = -1
    for row in readfile:       
        if (j >= 0):
           
            
            dx =float(row[0])
            theta =float(row[1])
            t.append(float(row[2]))
            H.append(float(row[3]))
            

            
            
                
        j = j + 1

n = len(t)        
t = array(t)
H = array(H)

Hdiff = abs(H[:] - H[0]) / H[0]

   