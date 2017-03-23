import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy import ones


wdir = "../../../../data/raw/bigsmoothtargetted/o3/4/4/"

timeinsecs = 30

gap = 1
g = 9.81
         
#filen = 200*2560
#s = wdir + "saveoutputts" + str(int(filen)) + ".txt"
 
#s = wdir + "saveoutputtslast.txt"
s = wdir +"outlast.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    h = []
    bed = []
    u = []
    he = []
    ue = []
    x = []
    ht = []
    ut = []
    j = -1
    for row in readfile:       
        if (j >= 0):
           
            
            dx =float(row[0])
            dt =float(row[1])
            t =float(row[2])
            x.append(float(row[3]))
            h.append(float(row[4]))
            u.append(float(row[6]))
            #he.append(float(row[8]))
            #ue.append(float(row[9]))
            #bed.append(float(row[7]))
            #diffuse = float(row[7])
            
            """
            dx =float(row[0])
            dt =float(row[1])
            t =float(row[2])
            x.append(float(row[5]))
            h.append(float(row[6]))
            u.append(float(row[8]))
            #diffuse = float(row[7])
            """
            
            
            """
            
            dx =float(row[0])
            dt =float(row[1])
            t =float(row[2])
            x.append(float(row[4]))
            h.append(float(row[5]))
            u.append(float(row[7]))
            """
            

            
            
                
        j = j + 1

    n = len(x)        
    x = array(x)
    u = array(u)
    h = array(h)
    bed = array(bed)
    
s = "hS.dat"
with open(s,'w') as file1:
    for i in range(51200,66560,4):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",h[i])
        file1.write(s)