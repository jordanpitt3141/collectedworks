# -*- coding: utf-8 -*-
"""
Created on Thu Sep 11 16:26:29 2014

@author: Jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog


wdirfords = ["o1","o2","o3","o3fem"]
wdirfnums = ["11","12","13","14"]

for wdirford in wdirfords:
    for wdirfnum in wdirfnums: 
        wdirf = "../../../data/Joe/alldb/" +wdirford + "/" 
        sdir = "../../results/show/steve1/nocorrect/"
             
        s = wdirf + wdirfnum+ ".txt"
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
        
        plot(xf,hf ,'-b')
        
        #xlim([300,700])
        #ylim([0.9,1.9])
        #xlim([0,1000])
        ylim([-0.1,2])
        xlim([0,1000])
        s = "Dam Break: " + wdirford + " dx = " + str(dxf)
        title(s)
        xlabel("x (m)")
        ylabel("h (m)")
        s = sdir + wdirford +"n" + wdirfnum+  "dbh.png"       
        savefig(s, bbox_inches='tight')        
        clf()  
        
        
        plot(xf,hf ,'-b')
        xlim([300,700])
        ylim([0.9,1.9])
        s = "Dam Break: " + wdirford + " dx = " + str(dxf)
        title(s)
        xlabel("x (m)")
        ylabel("h (m)")
        s = sdir + wdirford +"n" + wdirfnum+  "dbhz.png"       
        savefig(s, bbox_inches='tight')        
        clf()
        
        plot(xf,hf ,'-b')
        xlim([500,600])
        ylim([1.0,1.6])
        s = "Dam Break: " + wdirford + " dx = " + str(dxf)
        title(s)
        xlabel("x (m)")
        ylabel("h (m)")
        s = sdir + wdirford +"n" + wdirfnum+  "dbhzz.png"       
        savefig(s, bbox_inches='tight')        
        clf()
        
        
        plot(xf,uf ,'-b')
        xlim([0.0,1000])
        ylim([-0.1,1.9])
        s = "Dam Break: " + wdirford + " dx = " + str(dxf)
        title(s)
        xlabel("x (m)")
        ylabel("u (m/s)")
        s = sdir + wdirford +"n" + wdirfnum+  "dbu.png"       
        savefig(s, bbox_inches='tight')        
        clf()
        
        plot(xf,uf ,'-b')
        xlim([300,700])
        ylim([0.1,1.8])
        s = "Dam Break: " + wdirford + " dx = " + str(dxf)
        title(s)
        xlabel("x (m)")
        ylabel("u (m/s)")
        s = sdir + wdirford +"n" + wdirfnum+  "dbuz.png"       
        savefig(s, bbox_inches='tight')        
        clf()  
        
        plot(xf,uf ,'-b')
        xlim([500,600])
        ylim([0.8,1.2])
        s = "Dam Break: " + wdirford + " dx = " + str(dxf)
        title(s)
        xlabel("x (m)")
        ylabel("u (m/s)")
        s = sdir + wdirford +"n" + wdirfnum+  "dbuzz.png"       
        savefig(s, bbox_inches='tight')        
        clf()            
        
