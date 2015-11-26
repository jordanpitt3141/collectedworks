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

wdir = "../../../../data/raw/segur/o2/"
sdir = "../../../../data/postprocessing/segur/o2/"


g = 9.81
h1 = 0.1
h0 = 0.09
theta = 1.0

#xs = [0,5.0,10.0,15.0,20.0]
poss = [0,5.0,10.0,15.0,20.0]
for pos in poss:
    pos = pos + 0.61

    time = []
    height = []    
    
    filelist = listdir(wdir)
    n = len(filelist)
    for i in range(n):
        hs = []
        jp = 0
        s = wdir + filelist[i]
        with open(s,'r') as file1:
             readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
             j = -1
             for row in readfile:       
                 if (j >= 0):
                    dx = float(row[0])
                    dt = float(row[1])
                    hs.append(float(row[4]))
                    if(pos >= float(row[3]) - 0.5*dx and pos <= float(row[3]) + 0.5*dx):
                        time.append(float(row[2])) 
                        ix = float(row[3])
                        jp = j                    
                 j = j + 1
        ##get cell
        ih = hs[jp]
        ihm1 = hs[jp-1]
        ihp1 = hs[jp+1]
        dhif = ihp1 - ih
        dhim = 0.5*(ihp1 - ihm1)
        dhib = ih - ihm1
        dhi = minmod(theta*dhif,dhim,theta*dhib)
        nih = ih + (pos - ix)*(dhi/dx)
        height.append(nih)
         
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
    s = sdir + str(int(pos)) + ".png"       
    savefig(s, bbox_inches='tight')        
    clf()
    
    s = sdir + "scaledout" + str(int(pos)) +  ".txt"
    with open(s,'a') as file2:
            writefile1 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
            writefile1.writerow(['dx' ,'dt','time','height'])        
                           
            for i in range(n):
                 writefile1.writerow([str(dx),str(dt),str(time[i]), str(height[i])]) 
    print(pos)




