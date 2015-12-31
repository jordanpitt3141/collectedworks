import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog


#wdir = "../../../../data/raw/Joebigsmooth/o1/13/11/"
#wdir = "../../../../../../PhDold/project/code/data/Joe/alldb/o3/"

#wdir = "../../../../data/dbchris/o3femh110h01testsmallt/"


#wdir = "../../../../data/raw/Joebigsmooth/o1/13/9/"

wdir = "../../../../data/raw/Cserre/solitonothers/collDMcopy/o3/"

gap = 1
g = 9.81
         
#s = wdir + "outlast.txt"
 
s = wdir + "saveoutputtslast.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    h = []
    u = []
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
                
        j = j + 1
    x = array(x)
    u = array(u)
    h = array(h)
   