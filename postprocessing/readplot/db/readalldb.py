# -*- coding: utf-8 -*-
"""
Created on Thu Sep 11 16:26:29 2014

@author: Jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

wdir1 = "../../data/dambreak/centdiffprop0p5/"
wdir3 = "../../data/dambreak/centdiffa0p5/"
wdir2 = "../../data/dambreak/nfinezero/"
wdirf = "../../data/dambreak/fd0p5/"
sdir = "../db2/"

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
            bedd.append(float(row[6]))      
            
         j = j + 1
     xd = arange(0.0, 1000, dxd) 
     hd = array(hd)
     bedd = array(bedd)
     
s = wdir3 + "saveoutputtslast.txt"
with open(s,'r') as file1:
     readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     hs = []
     us = []
     beds = []
     j = -1
     for row in readfile:       
         if (j >= 0):
            dxs = float(row[0])
            dts = float(row[1])
            ts = float(row[2])
            hs.append(float(row[3]))
            us.append(float(row[5]))
            beds.append(float(row[6]))     
            
         j = j + 1
     xs = arange(0.0, 1000, dxs) 
     hs = array(hs)
     beds = array(beds)
     
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
     
s = wdir3 + "saveoutputts1.txt"
with open(s,'r') as file1:
     readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     ihs = []
     ius = []
     ibeds = []
     j = -1
     for row in readfile:       
         if (j >= 0):
            idxs = float(row[0])
            idts = float(row[1])
            its = float(row[2])
            ihs.append(float(row[3]))
            ius.append(float(row[5]))
            ibeds.append(float(row[6]))     
            
         j = j + 1
     ixs = arange(0.0, 1000, dxs) 
     ihs = array(ihs)
     ibeds = array(ibeds)
     



s = "Initial Conditions"
plot(ixs,ihs + ibeds,'-b', label=s)

s = "Bottom" 
plot(ixs,ibeds,'-g', label=s)

xlim([0,1000])
ylim([-0.1,2.0])
title("Dam Break")
xlabel("Distance (m)")
ylabel("Height (m)")
legend(loc='best')
s = sdir + "initdballh.png"       
savefig(s, bbox_inches='tight')        
clf()

s = "Initial Conditions"
plot(ixs,ius,'-b', label=s)

xlim([0,1000])
ylim([-2,2])
title("Dam Break")
xlabel("Distance (m)")
ylabel("Velocity (m/s)")
legend(loc='best')
s = sdir + "initdballv.png"       
savefig(s, bbox_inches='tight')        
clf() 


s = "FD"
plot(xf,hf + bedf,'-k', label=s)

s = "dx = 0.05"
plot(xv,hv + bedv,'-c', label=s)

s = "dx = 0.1"
plot(xs,hs + beds,'-r', label=s)

s = "dx = 0.5"
plot(xd,hd +bedd ,'-b', label=s)


s = "Bottom" 
plot(xv,bedv,'-g', label=s)

xlim([0,1000])
ylim([-0.1,2.0])
title("Dam Break")
xlabel("Distance (m)")
ylabel("Height (m)")
legend(loc='best')
s = sdir + "dballh.png"       
savefig(s, bbox_inches='tight')        
clf()


s = "FD"
plot(xf,hf + bedf,'-k', label=s)

s = "dx=0.05"
plot(xv,hv + bedv,'-c', label=s)

s = "dx=0.1"
plot(xs,hs + beds,'-r', label=s)

s = "dx=0.5"
plot(xd,hd +bedd ,'-b', label=s)


s = "Bottom" 
plot(xv,bedv,'-g', label=s)

xlim([400,800])
ylim([1.0,2.0])
title("Dam Break")
xlabel("Distance (m)")
ylabel("Height (m)")
legend(loc='best')
s = sdir + "dballhz.png"       
savefig(s, bbox_inches='tight')        
clf()


s = "FD"
plot(xf,uf,'-k', label=s)

s = "dx=0.05"
plot(xv,uv,'-c', label=s)

s = "dx=0.1"
plot(xs,us,'-r', label=s)

s = "dx=0.5"
plot(xd,ud,'-b', label=s)

xlim([0,1000])
ylim([-0.1,2.0])
title("Dam Break")
xlabel("Distance (m)")
ylabel("Velocity (m/s)")
legend(loc='best')
s = sdir + "dballv.png"       
savefig(s, bbox_inches='tight')        
clf()



