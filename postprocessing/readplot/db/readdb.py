# -*- coding: utf-8 -*-
"""
Created on Thu Sep 11 16:26:29 2014

@author: Jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

wdirf = "../../data/fddambreak/dx0p5l0p001/"
wdir1 = "../../data/dbt/dx0p5/"
wdir2 = "../../data/dbt/dx0p05/"
sdir = "../db2/"
     
s = wdirf + "saveoutputtslast.txt"
with open(s,'r') as file1:
     readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     hf = []
     uf = []
     bedf = []
     j = -1
     for row in readfile:       
         if (j >= 0):
            dxf = float(row[0])
            dtf = float(row[1])
            tf = float(row[2])
            hf.append(float(row[3]))
            uf.append(float(row[4]))
            bedf.append(float(row[5]))     
            
         j = j + 1
     xf = arange(0.0, 1000, dxf) 
     hf = array(hf)
     bedf = array(bedf)    
     
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

s = wdir2 + "saveoutputtslast.txt"
with open(s,'r') as file1:
     readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     h2 = []
     u2 = []
     bed2 = []
     j = -1
     for row in readfile:       
         if (j >= 0):
            dx2 = float(row[0])
            dt2 = float(row[1])
            t2 = float(row[2])
            h2.append(float(row[3]))
            u2.append(float(row[5]))
            bed2.append(float(row[6]))     
            
         j = j + 1
     x2 = arange(0.0, 1000, dx2) 
     h2 = array(h2)
     bed2 = array(bed2)      


s = "finite difference 0.5"
plot(xf,hf + bedf,'-k', label=s)

s = "Hybrid 0.05 (no lim)"
plot(x2,h2 + bed2,'-b', label=s)

s = "Hybrid 0.5 (no lim)"
plot(x1,h1 + bed1,'-r', label=s)


s = "Bottom" 
plot(xf,bedf,'-g', label=s)

xlim([0,1000])
ylim([-0.1,2.0])
title("Dam Break")
xlabel("Distance (m)")
ylabel("Height (m)")
legend(loc='best')
s = sdir + "dballh.png"       
savefig(s, bbox_inches='tight')        
clf()


s = "finite difference 0.5"
plot(xf,hf + bedf,'-k', label=s)

s = "Hybrid 0.05 (no lim)"
plot(x2,h2 + bed2,'-b', label=s)

s = "Hybrid 0.5 (no lim)"
plot(x1,h1 + bed1,'-r', label=s)

s = "Bottom" 
plot(xf,bedf,'-g', label=s)

xlim([400,800])
ylim([1.0,2.0])
title("Dam Break")
xlabel("Distance (m)")
ylabel("Height (m)")
legend(loc='best')
s = sdir + "dballhz.png"       
savefig(s, bbox_inches='tight')        
clf()

s = "finite difference 0.5"
plot(xf,hf + bedf,'-k', label=s)
s = "Bottom" 
plot(xf,bedf,'-g', label=s)

xlim([400,800])
ylim([1.0,2.0])
title("Dam Break")
xlabel("Distance (m)")
ylabel("Height (m)")
legend(loc='best')
s = sdir + "dbFDhz.png"       
savefig(s, bbox_inches='tight')        
clf()


s = "Hybrid 0.05 (no lim)"
plot(x2,h2 + bed2,'-b', label=s)


s = "Bottom" 
plot(xf,bedf,'-g', label=s)

xlim([400,800])
ylim([1.0,2.0])
title("Dam Break")
xlabel("Distance (m)")
ylabel("Height (m)")
legend(loc='best')
s = sdir + "db0p05hz.png"       
savefig(s, bbox_inches='tight')        
clf()


s = "finite difference 0.5"
plot(xf,uf,'-k', label=s)

s = "Hybrid 0.5 (no lim)"
plot(x2,u2,'-b', label=s)

s = "Hybrid 0.5 (no lim)"
plot(x1,u1,'-r', label=s)

xlim([0,1000])
ylim([-0.1,2.0])
title("Dam Break")
xlabel("Distance (m)")
ylabel("Velocity (m/s)")
legend(loc='best')
s = sdir + "dballv.png"       
savefig(s, bbox_inches='tight')        
clf()



