import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy import ones


#wdir = "../../../../data/raw/Joebigsmooth/o1/13/11/"
#wdir = "../../../../../../PhDold/project/code/data/Joe/alldb/o3/"

#wdir = "../../../../data/dbchris/o3femh110h01testsmallt/"


wdir = "../../../../data/raw/Joebigsmooth/o2/10/12/"

#wdir = "../../../../data/raw/Cserre/solitonothers/collDMcopy/o3/"
#wdir = "../../../../data/raw/bigsmoothalphainf/o3/10/0/"

#wdir = "../../../../data/raw/trackleadsola10new/o3/"

#wdir = "../../../../data/raw/DBASPECTRAT/o3/10/10/8.0/"

#wdir = "../../../../data/raw/longcontactdiscdx9diff10fileio/o3/9/0/"
#wdir = "../../../../data/raw/longcontactdiscdx8diff10fileio/o3/8/0/"

#wdir = "../../../../data/raw/LCD/o3/dx=9/diff=1000.0/"


#wdir = "../../../../data/raw/longcontactdisc/o3/dx=9/diff=10.0/"
#wdir = "../../../../data/raw/Joebigsmooth/o3/9/6/"

#wdir = "../../../../gitproj/data/raw/longcontactdiscdx9diff10fileio/o3/6/0/"

#wdir = "../../../../data/raw/DBSTEEPNoLimiter/o2/0/"

#wdir = "../../../../data/raw/longtimedambreakNEW/FDc/"

timeinsecs = 30

gap = 1
g = 9.81
         
filen = 200*2560
#s = wdir + "saveoutputts" + str(int(filen)) + ".txt"
 
#s = wdir + "saveoutputtslast.txt"
s = wdir + "outlast.txt"
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
            u.append(float(row[5]))
            #he.append(float(row[8]))
            #ue.append(float(row[9]))
            #bed.append(float(row[7]))
            diffuse = float(row[6])
            
            
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
   