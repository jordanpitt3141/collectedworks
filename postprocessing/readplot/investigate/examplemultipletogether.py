import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog


#wdir = "../../../../data/raw/Joebigsmooth/o1/13/11/"
#wdir = "../../../../../../PhDold/project/code/data/Joe/alldb/o3/"

#wdir = "../../../../data/dbchris/o3femh110h01testsmallt/"


#wdir = "../../../../data/raw/Joebigsmooth/o1/13/9/"

#wdir = "../../../../data/raw/Cserre/solitonothers/collDMcopy/o3/"
#wdir = "../../../../data/raw/bigsmoothalphainf/o3/10/0/"

#wdir = "../../../../data/raw/DSWalpha/o3/9/0/"
#wdir = "../../../../data/raw/DSWalpha/o3/9/1/"
wdir = "../../../../data/raw/DSWalpha/o3/9/3/"

gap = 1
g = 9.81

dx = 0.0
dt = 0.0
t = 0.0
diffuse = 0.0
         
s = wdir + "outlast.txt"
 
#s = wdir + "saveoutputtslast.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    h = []
    u = []
    x = []    
    hs = []
    us = []
    xs = []
    dxs = []
    dts = []
    ts = []
    diffuses = []
    j = 0
    for row in readfile:       
        if (row[0] == 'dx'):
            print(1,j,row[0])
            x = array(x)
            u = array(u)
            h = array(h)
            hs.append(h)
            us.append(u)
            xs.append(x)
            dxs.append(dx)
            dts.append(dt)
            ts.append(t)
            diffuses.append(diffuse)
            h = []
            u = []
            x = []
            
        else:
            dx =float(row[0])
            dt =float(row[1])
            t =float(row[2])
            x.append(float(row[3]))
            h.append(float(row[4]))
            u.append(float(row[6]))
            diffuse = float(row[7])
        j = j + 1
   