# -*- coding: utf-8 -*-
"""
Created on Thu Sep 11 16:26:29 2014

@author: Jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from matplotlib2tikz import save as tikz_save

exp1 = "test"
gap = 50

wdirords = ["o213m","o213o","o213p"]

#nums = [0,3,6,9,12]
#nums = [2,5,8,11,14]
#nums = [1,4,7,10,13]
nums = range(16,17)
for wdirord in wdirords:
    for k in nums:
        wdir = "../../../data/Joe2/smoothbd/" +wdirord+"/"+str(k) + "/" 
             
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
                    try :    
                        beta = float(row[8]) #could be 8 as well for o2
                    except :
                        beta = float(row[7])
                    
                 j = j + 1
             x = array(x)
             h = array(h)     
      
             xt =x[::gap]
             ht =h[::gap]
                
             plot(xt,ht ,'-')

sdir = "../../results/smoothdb/" + exp1 + "/"
        
if not os.path.exists(sdir):
    os.makedirs(sdir) 
xlim([0,1000])
ylim([0.8,2.0])
xlabel("x(m)")
ylabel("h(m)")
s = sdir + "db"+".tikz" 
tikz_save(s);           
clf()
