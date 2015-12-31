# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 16:53:37 2015

@author: jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from matplotlib2tikz import save as tikz_save
import os

dxw = "15"
wdirord = "o1"

    
wdir = "../../../../data/raw/Joesolconnon/"  +wdirord +"/" + dxw + "/"
sdir = "../../../../data/postprocessing/Joesolconnon/"  +wdirord +"/" + dxw + "/"
if not os.path.exists(sdir):
        os.makedirs(sdir)
     
s = wdir + "saveoutputtslast.txt"
with open(s,'r') as file1:
     readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     h = []
     u = []
     x = []
     j = -1
     for row in readfile:       
         if (j >= 0):
            dx = float(row[0])
            dt = float(row[1])
            t = float(row[2])
            h.append(float(row[3]))
            u.append(float(row[5]))
         j = j + 1
     x = arange(-600,600+dx,dx)#array(xt)
     h = array(h)     
  
plot(x,h)
xlabel("$x$ ($m$)")
ylabel("$h$ ($m$)")
#legend()

    

