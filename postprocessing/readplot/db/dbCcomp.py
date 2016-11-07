# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 15:08:15 2015

@author: Jordan
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Sep 11 16:26:29 2014

@author: Jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

wdir1 = "../../data/ct/cdb/"
wdir2 = "../../data/dambreak/centdiff0p5/"
sdir = "../Cdb/"

s = wdir1 + "saveoutputtslast.txt"
with open(s,'r') as file1:
     readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     hv = []
     uv = []
     bedv = []
     j = -1
     for row in readfile:       
         if (j >= 0):
            dxv = float(row[0])
            dtv = float(row[1])
            tv = float(row[2])
            hv.append(float(row[3]))
            uv.append(float(row[5]))
            bedv.append(float(row[6]))     
            
         j = j + 1
     xv = arange(0.0, 1000, dxv) 
     hv = array(hv)
     bedv = array(bedv)
     
s = wdir2 + "saveoutputtslast.txt"
with open(s,'r') as file1:
     readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     hd = []
     ud = []
     bedd = []
     j = -1
     for row in readfile:       
         if (j >= 0):
            dxd = float(row[0])
            dtd = float(row[1])
            td = float(row[2])
            hd.append(float(row[3]))
            ud.append(float(row[5]))
            
         j = j + 1
     xd = arange(0.0, 1000, dxd) 
     hd = array(hd)
     bedd = array(bedd)



s = "PY" 
plot(xd,hd,'-k', label=s)

s = "C"
plot(xv,hv + bedv,'--b', label=s)


xlim([0,1000])
ylim([-0.1,2.0])
title("Dam Break")
xlabel("Distance (m)")
ylabel("Height (m)")
legend(loc='best')
s = sdir + "dballh.png"       
savefig(s, bbox_inches='tight')        
clf()


s = "PY" 
plot(xd,hd,'-k', label=s)

s = "C"
plot(xv,hv + bedv,'--b', label=s)

xlim([400,800])
ylim([1.0,2.0])
title("Dam Break")
xlabel("Distance (m)")
ylabel("Height (m)")
legend(loc='best')
s = sdir + "dballhz.png"       
savefig(s, bbox_inches='tight')        
clf()


s = "PY" 
plot(xd,ud,'-k', label=s)

s = "C"
plot(xv,uv,'--b', label=s)
xlim([0,1000])
ylim([-0.1,2.0])
title("Dam Break")
xlabel("Distance (m)")
ylabel("Velocity (m/s)")
legend(loc='best')
s = sdir + "dballv.png"       
savefig(s, bbox_inches='tight')        
clf()



