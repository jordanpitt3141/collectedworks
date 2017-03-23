# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 16:53:37 2015

@author: jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
#from matplotlib2tikz import save as tikz_save
import os

dxw = "5"
wdirord = "o1"

#wdirords = ["o3","FDcent","grim","o2","o1"]
#ordsubtup = [[6,7],[5,6],[5,6], [6,8], [6,7]]
#removeint = [[515,545],[515,545],[515,545],[515,545],[515,545]]
#numss = [[1,6,8,9,12],[1,6,8,9,12],[1,6,8,9,12],[1,6,8,9,12],[1]]
#dxwss = [range(9,34),range(9,37),range(9,37),range(9,34),range(9,37)] 

#wdirords = ["o1"]
#ordsubtup = [[6,7]]
#removeint = [[515,545]]

#nums = [1]
#dxws = range(1,513)

#wdirords = ["o3","FDcent","grim","o2"]
#ordsubtup = [[6,7],[5,6],[5,6], [6,8]]
#removeint = [[515,545],[515,545],[515,545],[515,545]]

#nums = [1]
#dxws = range(1,513)

#wdirords = ["o3"]
#ordsubtup = [[6,7],[5,6],[5,6], [6,8]]
#removeint = [[515,545],[515,545],[515,545],[515,545]]
"""
wdirords = ["FDcent","grim"]
ordsubtup = [[5,6],[5,6]]
removeint = [[515,545],[515,545]]
nums = [1,6,8,9,12]
dxws = range(1,513)
"""



#wdirords = ["o3","o2","o1"]
#ordsubtup = [[6,7], [6,8],[6,7]]
wdirords = ["o3"]
ordsubtup = [[6,7]]
#removeint = [[515,545],[515,545],[515,545]]
#removeint = [[515,545]]
removeint = [[0,0]]
nums = [1]
dxws =  range(2,6) + range(6,14)[::2] +   range(16,513)

"""

wdirords = ["FDcent","grim"]
ordsubtup = [[5,6],[5,6]]
removeint = [[515,545],[515,545]]
nums = [1,6,8,9,12]
dxws =  range(1,6) + range(6,14)[::2] +   range(16,513)
"""

for ip in range(len(wdirords)):
    for jp in nums:
        hdx = dxws[0]*(10.0/(2**12))
        hs = []
        us = []
        xs = []
        ts = []
        dxs = []
        dts = []
        diffs = []
        normhs =[]
        normus =[]
        
        diff = str(jp)
        wdirord = wdirords[ip]
                    
        for k in dxws:
            wdir = "../../../../../data/raw/bigsmoothtargetted/"  +wdirord +"/" + str(k)+ "/" + diff + "/"
            sdir = "../../../../../data/postprocessing/smoothdbtarget/WITHMID/1diffmdxcomreal/" + wdirord + "/" +diff+ "/"
            if not os.path.exists(sdir):
                   os.makedirs(sdir)
            

            s = wdir + "outlast.txt"
            if os.path.exists(s):
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
                            u.append(float(row[ordsubtup[ip][0]])) #5 for FDcent, 6 otherwise 
                            beta = float(row[ordsubtup[ip][1]]) #could be 8 as well for o2
                        j = j + 1
                    
                    x = array(x)
                    h = array(h)
                    u = array(u)
                    hs.append(h)
                    us.append(u)
                    xs.append(x)
                    dxs.append(dx)
                    dts.append(dt)
                    ts.append(t)
                    diffs.append(beta)

        s = sdir + "norms.txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                    
            writefile2.writerow(["dx","time","diff" ,"height norm","velocity norm"])           
        # Comparison
        for i in range(len(dxs)):
            dx = dxs[i]
            gap = int(dx/hdx)
            h1 = hs[i].tolist()
            x1 = xs[i].tolist()
            u1 = us[i].tolist()
            h2 = hs[0][::gap].tolist()
            u2 = us[0][::gap].tolist()
            x2 = xs[0][::gap].tolist()
            zbeg = int(removeint[ip][0]/dx)
            zend = int(removeint[ip][1]/dx)
            
            if(len(h1) != len(h2)):
                h1 = h1[:-1]
                u1 = u1[:-1]
                #print(len(h1),len(h2))
                #print("gap")
                #print(dx/hdx,gap)
                

            print(zbeg,zend,dx)
            h1 = h1[0:zbeg] + zeros(zend-zbeg).tolist()  + h1[zend:]
            u1 = u1[0:zbeg] + zeros(zend-zbeg).tolist()  + u1[zend:]
            h2 = h2[0:zbeg] + zeros(zend-zbeg).tolist()  + h2[zend:]
            u2 = u2[0:zbeg] + zeros(zend-zbeg).tolist()  + u2[zend:]
            print(wdirord,len(h1),len(h2))
            h1 = array(h1)
            u1 = array(u1)
            h2 = array(h2)
            u2 = array(u2)
            normh = norm(h1 - h2,ord=1) / norm(h2,ord=1)
            normu = norm(u1 - u2,ord=1) / norm(u2,ord=1)
            
            #plot(xs[i], h1,'--',label=str(dx))
            #plot(xs[i], h2,'-',label=str(dx))
            #legend()
            
            s = sdir + "norms.txt"
            with open(s,'a') as file1:
                writefile1 = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                    
                writefile1.writerow([str(dxs[i]),str(ts[i]),str(diffs[i]) ,str(normh),str(normu)]) 
                
            s = sdir + "h.dat"
            with open(s,'a') as file1:
                 s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",normh)
                 file1.write(s)
                
            s = sdir + "u.dat"
            with open(s,'a') as file1:
                s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",normu)
                file1.write(s) 
            
            normus.append(normu)
            normhs.append(normh)