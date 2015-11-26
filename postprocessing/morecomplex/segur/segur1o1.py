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
from os import listdir

wdir = "../../../../data/raw/segur/o1/"
sdir = "../../../../data/postprocessing/segur/o1/"


g = 9.81
dx = 0.1
h1 = 0.1
h0 = 0.09

#xs = [0.0+dx]
poss = [0,5.0,10.0,15.0,20.0]
for pos in poss:
    pos = pos + 0.61
    time = []
    height = []
    
    filelist = listdir(wdir)
    n = len(filelist)
    for i in range(n):
        s = wdir + filelist[i]
        with open(s,'r') as file1:
             readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
             j = -1
             for row in readfile:       
                 if (j >= 0):
                    dx = float(row[0])
                    dt = float(row[1])
                     
                    if(pos >= float(row[3]) - 0.5*dx and pos <= float(row[3]) + 0.5*dx):
                        #print(pos,float(row[3]),float(row[3]) - 0.5*dx,float(row[3]) + 0.5*dx)
                        height.append(float(row[4]))
                        time.append(float(row[2]))    
                    
                 j = j + 1
         
    ###### sort them in order ######
    mix = zip(time,height)
    mix = sorted(mix)
    time = [p[0] for p in mix]
    height = [p[1] for p in mix]
    
    s = sdir + "out" + str(int(pos)) +  ".txt"
    with open(s,'a') as file2:
            writefile1 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
            writefile1.writerow(['dx' ,'dt','time','height'])        
                           
            for i in range(n):
                 writefile1.writerow([str(dx),str(dt),str(time[i]), str(height[i])]) 
    
    time = array(time)
    height = array(height)
    time = sqrt(g/h1)*time - pos/h1 
    height = (height - h1)/h1
    
    s = "Segur x/h1="+str(pos/h1)
    plot(time,height,'.')
    xlim([-10,250])
    ylim([-0.1,0.1])
    title(s)
    xlabel("Scaled t,x")
    ylabel("Scaled Height")
    s = sdir +"o1" + str(int(pos)) + ".png"       
    savefig(s, bbox_inches='tight')        
    clf()
    
    print pos



