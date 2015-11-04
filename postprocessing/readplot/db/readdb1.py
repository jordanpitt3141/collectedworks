# -*- coding: utf-8 -*-
"""
Created on Thu Sep 11 16:26:29 2014

@author: Jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

wdir1 = "../../data/dbt/sup0p5/"
sdir = "../db2/"
     
     
s = wdir1 + "saveoutputtslast.txt"
with open(s,'r') as file1:
     readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     h1 = []
     u1 = []
     bed1 = []
     j = -1
     for row in readfile:       
         if (j >= 0):
            dx1 = float(row[0])
            dt1 = float(row[1])
            t1 = float(row[2])
            h1.append(float(row[3]))
            u1.append(float(row[5]))
            bed1.append(float(row[6]))     
            
         j = j + 1
     x1 = arange(0.0, 1000, dx1) 
     h1 = array(h1)
     bed1 = array(bed1)  
   


s = "Hybrid 0.5 (no lim)"
plot(x1,h1 + bed1,'-r', label=s)


s = "Bottom" 
plot(x1,bed1,'-g', label=s)

xlim([0,1000])
ylim([-0.1,12.0])
title("Dam Break")
xlabel("Distance (m)")
ylabel("Height (m)")
legend(loc='best')
s = sdir + "dballh.png"       
savefig(s, bbox_inches='tight')        
clf()


s = "Hybrid 0.5 (no lim)"
plot(x1,h1 + bed1,'-r', label=s)

s = "Bottom" 
plot(x1,bed1,'-g', label=s)

xlim([400,800])
ylim([1.0,12.0])
title("Dam Break")
xlabel("Distance (m)")
ylabel("Height (m)")
legend(loc='best')
s = sdir + "dballhz.png"       
savefig(s, bbox_inches='tight')        
clf()



s = "Hybrid 0.5 (no lim)"
plot(x1,u1,'-r', label=s)

xlim([0,1000])
ylim([-0.1,12.0])
title("Dam Break")
xlabel("Distance (m)")
ylabel("Velocity (m/s)")
legend(loc='best')
s = sdir + "dballv.png"       
savefig(s, bbox_inches='tight')        
clf()



