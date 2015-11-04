# -*- coding: utf-8 -*-
"""
Created on Thu Sep 11 16:26:29 2014

@author: Jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

betanum = "10"

wdir = "../../../data/dc/db/loop/1/"+betanum+"/"
wdirf = "../../../data/Joe/alldb/o2/"
sdir = "../../results/show/steve1/correct/"
     
s = wdirf + "12.txt"
with open(s,'r') as file1:
     readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     hf = []
     uf = []
     xf = []
     j = -1
     for row in readfile:       
         if (j >= 0):
            dxf = float(row[0])
            dtf = float(row[1])
            tf = float(row[2])
            xf.append(float(row[3]))
            hf.append(float(row[4]))
            uf.append(float(row[6])) 
            #beta = float(row[7])
            
         j = j + 1
     xf = array(xf)
     hf = array(hf)
     
s = wdir + "outlast.txt"
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
            x.append(float(row[3]))
            h.append(float(row[4]))
            u.append(float(row[6])) 
            beta = float(row[7])
            
         j = j + 1
     x = array(x)
     h = array(h)     

s = str(0)
plot(xf,hf ,'-b', label=s)
s = str(beta)
plot(x,h ,'-r', label=s)

ylim([-0.1,2])
xlim([0,1000])
title("Dam Break: dx =" + str(dx))
xlabel("x (m)")
ylabel("h (m)")
legend(loc='best')
s = sdir + betanum+  "dbh.png"       
savefig(s, bbox_inches='tight')        
clf()  


s = str(0)
plot(xf,hf ,'-b', label=s)
s = str(beta)
plot(x,h ,'-r', label=s)
xlim([300,700])
ylim([0.9,1.9])
title("Dam Break: dx =" + str(dx))
xlabel("x (m)")
ylabel("h (m)")
legend(loc='best')
s = sdir + betanum+  "dbhz.png"       
savefig(s, bbox_inches='tight')        
clf()

s = str(0)
plot(xf,hf ,'-b', label=s)
s = str(beta)
plot(x,h ,'-r', label=s)
xlim([500,600])
ylim([1.0,1.6])
title("Dam Break: dx =" + str(dx))
xlabel("x (m)")
ylabel("h (m)")
legend(loc='best')
s = sdir + betanum+   "dbhzz.png"       
savefig(s, bbox_inches='tight')        
clf()

s = str(0)
plot(xf,hf ,'-b', label=s)
s = str(beta)
plot(x,h ,'-r', label=s)
xlim([380,390])
ylim([1.6,1.8])
title("Dam Break: dx =" + str(dx))
xlabel("x (m)")
ylabel("h (m)")
legend(loc='best')
s = sdir + betanum+   "dbhzzslope.png"       
savefig(s, bbox_inches='tight')        
clf()


s = str(0)
plot(xf,uf ,'-b', label=s)
s = str(beta)
plot(x,u ,'-r', label=s)
xlim([0.0,1000])
ylim([-0.1,1.9])
title("Dam Break: dx =" + str(dx))
xlabel("x (m)")
ylabel("u (m/s)")
legend(loc='best')
s = sdir + betanum+   "dbu.png"       
savefig(s, bbox_inches='tight')        
clf()

s = str(0)
plot(xf,uf ,'-b', label=s)
s = str(beta)
plot(x,u ,'-r', label=s)
xlim([300,700])
ylim([0.1,1.8])
title("Dam Break: dx =" + str(dx))
xlabel("x (m)")
ylabel("u (m/s)")
legend(loc='best')
s = sdir + betanum+   "dbuz.png"       
savefig(s, bbox_inches='tight')        
clf()  

s = str(0)
plot(xf,uf ,'-b', label=s)
s = str(beta)
plot(x,u ,'-r', label=s)
xlim([500,600])
ylim([0.8,1.2])
title("Dam Break: dx =" + str(dx))
xlabel("x (m)")
ylabel("u (m/s)")
legend(loc='best')
s = sdir + betanum+   "dbuzz.png"       
savefig(s, bbox_inches='tight')        
clf()            
    
