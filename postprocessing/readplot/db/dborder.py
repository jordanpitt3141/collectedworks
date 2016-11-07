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

wdir1 = "../../data/dambreak/nolim/"
wdir2 = "../../data/ppm/dba/"
wdir3 = "../../data/limit/tdb/"


sdir = "../order/dba/"

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
     
s = wdir2 + "savelast.txt"
with open(s,'r') as file1:
     readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     hd = []
     ud = []
     bedd = []
     j = -1
     for row in readfile:       
         if (j >= 0):
            dxd = float(row[1])
            td = float(row[0])
            hd.append(float(row[2]))
            ud.append(float(row[3]))
            
         j = j + 1
     xd = arange(0.0, 1000, dxd) 
     hd = array(hd)
     bedd = array(bedd)
     
s = wdir3 + "savelast.txt"
with open(s,'r') as file1:
     readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     hs = []
     us = []
     beds = []
     j = -1
     for row in readfile:       
         if (j >= 0):
            dxs = float(row[1])
            ts = float(row[0])
            hs.append(float(row[2]))
            us.append(float(row[3]))
            
         j = j + 1
     xs = arange(0.0, 1000, dxs) 
     hs = array(hs)
     beds = array(beds)

s = "3rd Order PPM" 
plot(xd,hd,'-k', label=s)

s = "3rd Order Koren" 
plot(xs,hs,'-c', label=s)

s = "2nd Order"
plot(xv,hv + bedv,'-r', label=s)


xlim([0,1000])
ylim([-0.1,2.0])
title("Dam Break")
xlabel("Distance (m)")
ylabel("Height (m)")
legend(loc='best')
s = sdir + "dballh.png"       
savefig(s, bbox_inches='tight')        
clf()

s = "3rd Order PPM" 
plot(xd,hd,'-k', label=s)

s = "3rd Order Koren" 
plot(xs,hs,'-c', label=s)

s = "2nd Order"
plot(xv,hv + bedv,'-r', label=s)

xlim([400,800])
ylim([1.0,2.0])
title("Dam Break")
xlabel("Distance (m)")
ylabel("Height (m)")
legend(loc='best')
s = sdir + "dballhz.png"       
savefig(s, bbox_inches='tight')        
clf()


s = "Difference Between PPM and 2nd"
plot(xv,hd - hv,'-k', label=s)

s = "Difference Between Koren and 2nd"
plot(xv,hs - hv,'-r', label=s)

xlim([400,800])
title("Dam Break")
xlabel("Distance (m)")
ylabel("Height (m)")
legend(loc='best')
s = sdir + "diffhz.png"       
savefig(s, bbox_inches='tight')        
clf()

s = "3rd Order PPM"  
plot(xd,ud,'-k', label=s)

s = "3rd Order Koren"  
plot(xs,us,'-c', label=s)

s = "2nd Order"
plot(xv,uv,'-r', label=s)


xlim([0,1000])
ylim([-0.1,2.0])
title("Dam Break")
xlabel("Distance (m)")
ylabel("Velocity (m/s)")
legend(loc='best')
s = sdir + "dballv.png"       
savefig(s, bbox_inches='tight')        
clf()


s = "3rd Order PPM"  
plot(xd,ud,'-k', label=s)

s = "3rd Order Koren"  
plot(xs,us,'-c', label=s)

s = "2nd Order"
plot(xv,uv,'-r', label=s)


xlim([400,800])
title("Dam Break")
xlabel("Distance (m)")
ylabel("Velocity (m/s)")
legend(loc='best')
s = sdir + "dballvz.png"       
savefig(s, bbox_inches='tight')        
clf()



