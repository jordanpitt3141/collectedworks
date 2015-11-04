# -*- coding: utf-8 -*-
"""
Created on Thu Sep 11 16:26:29 2014

@author: Jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

wdir = "../../data/Cserre/segur1/0p005/o1/"
sdir = "../segur1/0p005/o1/"

inum = 5000
     
s = wdir + "out"+str(inum) + ".txt"
with open(s,'r') as file1:
     readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     h = []
     u = []
     x = []
     bed = []
     j = -1
     for row in readfile:       
         if (j >= 0):
            dx = float(row[0])
            dt = float(row[1])
            t = float(row[2])
            x.append(float(row[3]))
            h.append(float(row[4]))
            u.append(float(row[6]))  
            
         j = j + 1
     x = array(x)
     h = array(h)
     


plot(x,h ,'.')

xlim([-5,5])
s = "Segur at t="+str(t)
title(s)
xlabel("x (m)")
ylabel("Height (m)")
s = sdir +"o1segurtime"+str(t)+ ".png"       
savefig(s, bbox_inches='tight')        
clf()

