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
        
def reconstructppm(u,i,dx):
    ithr = 1.0 /3.0
    
    daip1 = 0.5*(u[i+2] - u[i])
    dip1 =((u[i+2] - u[i+1])*(u[i+1] - u[i]) >0)
    dip1 = dip1*min([abs(daip1),2*abs(u[i+2] - u[i+1]), 2*abs(u[i+1] - u[i])])*sign(daip1) 
          
    dai = 0.5*(u[i+1] - u[i-1])
    di =((u[i+1] - u[i])*(u[i] - u[i-1]) >0)
    di = di*min([abs(dai),2*abs(u[i+1] - u[i]), 2*abs(u[i] - u[i-1])])*sign(dai) 
    
    uir = u[i] + 0.5*(u[i+1] - u[i]) + (0.5*ithr)*(di - dip1)
    
    daim1 = 0.5*(u[i] - u[i-2])
    dim1 =((u[i] - u[i-1])*(u[i-1] - u[i-2]) >0)
    dim1 = dim1*min([abs(daim1),2*abs(u[i] - u[i-1]), 2*abs(u[i-1] - u[i-2])])*sign(daim1) 
    uil = u[i-1] + 0.5*(u[i] - u[i-1]) + (0.5*ithr)*(dim1 - di)
    
    #local extrema
    lce = (uir - u[i])*(u[i] - uil)   
    
    uir = u[i]*(lce <= 0) + uir*(lce > 0)
    uil = u[i]*(lce <= 0) + uil*(lce > 0)
        
    #monotonicity
    toclosellhs = (uir- uil)*(u[i] - 0.5*(uil + uir))
    tocloselrhs = (uir- uil)*(uir- uil)/6.0
    
    tocloserlhs = (uir- uil)*(u[i] - 0.5*(uil + uir))
    tocloserrhs = -(uir- uil)*(uir- uil)/6.0
    
    uil = (3*u[i] - 2*uir)*(toclosellhs > tocloselrhs) + uil*(toclosellhs <= tocloselrhs)
    uir = (3*u[i] - 2*uil)*(tocloserlhs < tocloserrhs) + uir*(tocloserlhs >= tocloserrhs)
    
    return uil,uir


import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from os import listdir

wdir = "../../../../data/segur2/o3p/"
sdir = "../../../results/segurf/o3p/"


g = 9.81
h1 = 0.1
h0 = 0.09
theta = 1.0
i3 = 1.0/3.0

#poss= [5.0,10.0,15.0,20.0]
poss =[0.1]
for pos in poss:
    pos = pos
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
        ihl,ihr = reconstructppm(hs,jp,dx)
        ai = (3*ihl + 3*ihr - 6*ih)/(dx*dx)
        bi = (ihr - ihl)/ dx
        ci  = 0.25*(6*ih - ihl - ihr)
        dfx = (pos - ix)
        nih = ai*dfx*dfx + bi*dfx + ci
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
    ylim([-0.3,0.2])
    title(s)
    xlabel("Scaled t,x")
    ylabel("Scaled Height")
    s = sdir +"o1" + str(int(pos)) + ".png"       
    savefig(s, bbox_inches='tight')        
    clf()
    
    s = sdir + "scaledout" + str(int(pos)) +  ".txt"
    with open(s,'a') as file2:
            writefile1 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
            writefile1.writerow(['dx' ,'dt','time','height'])        
                           
            for i in range(n):
                 writefile1.writerow([str(dx),str(dt),str(time[i]), str(height[i])]) 
    print(pos)




