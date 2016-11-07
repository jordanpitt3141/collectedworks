# -*- coding: utf-8 -*-
"""
Created on Sun Mar 01 17:51:43 2015

@author: Jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from os import listdir

wdir1 = "../../data/dambreak/nfineone/"
wdir2 = "../../data/dambreak/nfinetwo/"
wdir3 = "../../data/dambreak/nfinezero/"
wdirf = "../../data/dambreak/fd0p5/"
sdir = "../db/alltime/"

for i in listdir(wdir1):
    
    s = wdir1 + i
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
     
    s = wdir2 + i
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
     
    s = wdir3 + i
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
     
    s = wdirf + i
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
            
    s = "FD"
    plot(xf,hf + bedf,'-k', label=s)
    
    s = "Zero"
    plot(xs,hs + beds,'-c', label=s)
    
    s = "Two"
    plot(xd,hd +bedd ,'-r', label=s)
    
    s = "One"
    plot(xv,hv + bedv,'-b', label=s)
    
    
    s = "Bottom" 
    plot(xv,bedv,'-g', label=s)
    
    xlim([0,1000])
    ylim([-0.1,2.0])
    title("Dam Break at " + str(tf) + "s" )
    xlabel("Distance (m)")
    ylabel("Height (m)")
    legend(loc='best')
    s = sdir +"height/" + str(tf) + ".png"       
    savefig(s, bbox_inches='tight')        
    clf()
    
    
    s = "FD"
    plot(xf,hf + bedf,'-k', label=s)
    
    s = "Zero"
    plot(xs,hs + beds,'-c', label=s)
    
    s = "Two"
    plot(xd,hd +bedd ,'-r', label=s)
    
    s = "One"
    plot(xv,hv + bedv,'-b', label=s)
    
    
    s = "Bottom" 
    plot(xv,bedv,'-g', label=s)
    
    xlim([400,800])
    ylim([1.0,2.0])
    title("Dam Break at " + str(tf) + "s" )
    xlabel("Distance (m)")
    ylabel("Height (m)")
    legend(loc='best')
    s = sdir +"zoomheight/" + str(tf) + ".png"       
    savefig(s, bbox_inches='tight')        
    clf()
    
    
    s = "FD"
    plot(xf,uf,'-k', label=s)
    
    s = "Zero"
    plot(xs,us,'-c', label=s)
    
    s = "Two"
    plot(xd,ud,'-r', label=s)
    
    s = "One"
    plot(xv,uv,'-b', label=s)
    
    xlim([0,1000])
    ylim([-0.1,2.0])
    title("Dam Break at " + str(tf) + "s" )
    xlabel("Distance (m)")
    ylabel("Velocity (m/s)")
    legend(loc='best')
    s = sdir +"velocity/" + str(tf) + ".png"        
    savefig(s, bbox_inches='tight')        
    clf()
     
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
