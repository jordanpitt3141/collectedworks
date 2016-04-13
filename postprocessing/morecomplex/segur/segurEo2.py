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
def minmod(a, b,  c):
    if((a > 0) and (b>0) and  (c>0)):
        return min(a,b,c);
    elif((a < 0) and  (b<0) and  (c<0)):
        return max(a,b,c);
    else:
        return 0.0;


import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from os import listdir
import os

wdir = "../../../../data/raw/segur/o2E/"
sdir = "../../../../data/postprocessing/segur/o2E/"

if not os.path.exists(sdir):
    os.makedirs(sdir)

g = 9.81
h1 = 0.1
h0 = 0.09
theta = 1.2

#xs = [0,5.0,10.0,15.0,20.0]
#poss = [0,5.0,10.0,15.0,20.0]

time = []
Energ = []    

filelist = listdir(wdir)
n = len(filelist)
for i in range(n):
    hs = []
    jp = 0
    s = wdir + filelist[i]
    with open(s,'r') as file1:
        lines = file1.readlines()
        rl = lines[1].split(",")
        dx = float(rl[0])
        dt = float(rl[1])
        t = float(rl[2])
        E = float(rl[3])
         

    ##get cell
    Energ.append(E)
    time.append(t)
         
###### sort them in order ######
mix = zip(time,Energ)
mix = sorted(mix)
time = [p[0] for p in mix]
Energ = [p[1] for p in mix]

s = sdir + "outE.txt"
with open(s,'a') as file2:
        writefile1 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
        writefile1.writerow(['time' ,'Energy'])        
                       
        for i in range(n):
             writefile1.writerow([str(time[i]),str(Energ[i])]) 

time = array(time)
Energ = array(Energ)

s = "E vs t"
plot(time,Energ,'.')
xlim([0.0,50])
ylim([5.8,5.9])
title(s)
xlabel("t")
ylabel("Energ")
s = sdir + "Evt.png"       
savefig(s, bbox_inches='tight')        
clf()
    




