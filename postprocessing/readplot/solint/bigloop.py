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

dxw = "5"
wdirord = "o1"


wdirords = ["o3","FDc","grim","o2","o1"]
ordsubtup = [[6,7],[5,6],[5,6], [6,8], [6,7]]

#wdirords = ["o3","o2","o1"]
#ordsubtup = [[6,7], [6,8], [6,7]]
for ip in range(len(wdirords)):
    wdirord = wdirords[ip]
    
    ylims = [[0.8,2],[0.8,2]]
    xlims = [[0,400],[0,200]]            
    
    ## for dx 8
    gaps = [4,2]
    zoom = [False,False]
    zoomints = [[500,560],[500,560]]
    zoomgap = [1,1]
    
    
    
    
    for l in range(len(ylims)):
        
        cylim = ylims[l]
        cxlim = xlims[l]
        
        gap = gaps[l]
        wdir = "../../../../data/raw/Cserre/solitonothers/collDMcopy/"  +wdirord +"/"
        sdir = "../../../../data/postprocessing/Cserre/solitonothersfix/collDMcopy/"  +wdirord +"/"
        if not os.path.exists(sdir):
                os.makedirs(sdir)
             
        s = wdir + "saveoutputtslast.txt"
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
                 j = j + 1
             igap = int(gap)
             if zoom[l]:
                 xbeg = int(cxlim[0]/dx)
                 xend = int(cxlim[1]/dx)
                 zbeg = int(zoomints[l][0]/dx)
                 zend = int(zoomints[l][1]/dx)
                 xt = x[xbeg:zbeg:igap] + x[zbeg:zend:zoomgap[l]] + x[zend:xend:igap]
                 ht = h[xbeg:zbeg:igap] + h[zbeg:zend:zoomgap[l]] + h[zend:xend:igap]
             else:
                 xbeg = int(cxlim[0]/dx)
                 xend = int(cxlim[1]/dx)
                 xt = x[xbeg:xend:igap]
                 ht = h[xbeg:xend:igap]
             x = array(xt)
             h = array(ht)      
        plot(x,h)
        ylim(cylim)
        xlim(cxlim)
        xlabel("$x$ ($m$)")
        ylabel("$h$ ($m$)")
        #legend()

            
        """
        plot(x,u ,'-b')
        ylim([-0.1,2])
        xlim([0,1000])
        s = "Dam Break: " + wdirord + " dx = " + str(dx)
        title(s)
        xlabel("x (m)")
        ylabel("u (m/s)")
        """
    
        
        stikz = sdir +str(l)+ ".tikz" 
        tikz_save(stikz);
        s = "Soliton Collision: " + wdirord + " dx = " + str(dx)
        title(s)
        s = sdir +str(l)+".png"       
        savefig(s, bbox_inches='tight')
        legend()
        s = sdir +str(l)+ "leg.png"       
        savefig(s, bbox_inches='tight')        
        clf()
        
        n = len(x)
        s = sdir + "h" +str(l)+ ".dat"
        with open(s,'w') as file1:
            for i in range(n):
                s ="%3.8f%5s%1.15f\n" %(x[i]," ",h[i])
                file1.write(s)
        s = sdir + "u" +str(l)+ ".dat"
        with open(s,'w') as file2:
            for i in range(n):
                s ="%3.8f%5s%1.15f\n" %(x[i]," ",u[i])
                file2.write(s)
    
        """
        \pgfplotsset{compat = newest,
    	every x tick label/.append style={font=\scriptsize, yshift=0.5ex},
    	every y tick label/.append style={font=\scriptsize, xshift=0.5ex}}
        """
    
        #make the tex file to create the document
        s = "\\documentclass[]{article}\n\\usepackage{tikz} \n\\usepackage{amsmath} \n" \
        + "\\usepackage{pgfplots}\n\\usepgfplotslibrary{external}\n\\tikzexternalize\n" \
        + "\\pgfplotsset{compat=newest,every x tick label/.append style={font=\scriptsize}," \
        + "every y tick label/.append style={font=\scriptsize}, every axis plot post/.style={line join=round}}\n"\
        + "\\begin{document}\n   \\input{"+str(l)+ ".tikz" +"}\n\\end{document}"
        
    
        #print(s)
        filen = sdir +str(l)+ ".tex" 
        """
        if not os.path.exists(filen):
            with open(filen,'a') as file1:
                print("created" +filen)
        """
        file1 = open(filen, 'w')
        file1.write(s)
        file1.close() 

    #make makefile
    s = "LAT = pdflatex \nLATFLAGS = -shell-escape\n\n"
    for i in range(len(ylims)):
        newl = str(i) +":\n\t $(LAT) $(LATFLAGS) " +str(i)+".tex\n\n"
        s = s + newl
    newl = "clean:\n\t rm -f  *~ ./*.log ./*.aux ./*.auxlock ./*.dep ./*.dpth ./*.pdf ./*.gz\n\n" 
    s = s + newl
    
    newl = "all: "
    for i in range(len(ylims)):
        newl = newl + str(i) + " "
    s = s + newl
       
    filen =sdir + "Makefile" 
    
    file1 = open(filen, 'w')
    file1.write(s)
    file1.close()
    call(['make','-C',sdir,'all'])           
