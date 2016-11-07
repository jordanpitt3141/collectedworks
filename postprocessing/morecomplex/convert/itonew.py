# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 16:53:37 2015

@author: jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from matplotlib2tikz import save as tikz_save
import os
from subprocess import call

wdirords = ["o3","FDcent","grim","o2","o1"]

#wdirords = ["o1"]

wdirbase = "../../../../data/raw/Joebigsmooth/"

wdirlist = os.listdir(wdirbase)

for wdir in wdirlist:
    wdir1 = wdirbase + wdir + "/"
    dxws = os.listdir(wdir1)
    #print(wdir,dxws)
    for dx in dxws:
        idx = int(dx)
        wdir2 = wdirbase + wdir + "/" + dx + "/"
        diffs = os.listdir(wdir2)
        
        for diff in diffs:
            idiff = int(diff)
            #print(wdir,dx,diffs)
            wdir3 = wdirbase + wdir + "/" + dx + "/" + diff + "/"
            #print(wdir3)

            sdir = "../../../../data/raw/bigsmoothtargetted/"  +wdir +"/" + str(2**(12 - idx)) + "/" + diff + "/"
            
            
            if not os.path.exists(sdir):
                    os.makedirs(sdir)
                 
            w = wdir3 + "outlast.txt"
            s = sdir + "outlast.txt"
            call(['cp',w,s]) 
                    
