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

diffw = "11"
wdirord = "o3"

#wdirords = ["o3","FDcent","grim","o2","o1"]
#ordsubtup = [[6,7],[5,6],[5,6], [6,8], [6,7]]

wdirords = ["o3","o2","o1"]
ordsubtup = [[6,7],[6,8], [6,7]]

#wdirords = ["FDcent","grim"]
#ordsubtup = [[5,6],[5,6]]


diffws = [6,9,10,20]

for ip in range(len(wdirords)):
    for jp in diffws:
        
        diffw = str(jp)
        wdirord = wdirords[ip]

        #nums = [0,3,6,9,12]
        #nums = [2,5,8,11,14]
        #nums = [1,4,7,10,13]
        numbeg = 4
        numend = 11
        nums = range(numbeg,numend)
        
        ylims = [[0.0,2],[1.0,1.8],[1.3,1.45],[1.3,1.45],[1.3,1.45]]
        xlims = [[0,1000],[300,700],[500,560],[520,540],[528,536]]
        gaps = [[2,4,6,8,10,12,14],[1,1,2,6,10,14,14],[1,1,1,1,1,1,1],[1,1,1,1,1,1,1],[1,1,1,1,1,1,1]]
        zoom = [True,True,False,False,False]
        zoomints = [[500,560],[500,560],[500,560],[500,560],[530,540]]
        zoomgap = [2,1,1,1,1]
        
        
        for l in range(len(ylims)):
            
            cylim = ylims[l]
            cxlim = xlims[l]
            
            
            for k in nums:
                gap = gaps[l][k-numbeg]
                
                wdir = "../../../../../../data/Joesmooth/bigsmooth/"  +wdirord +"/" + str(k) + "/" + diffw + "/"
                sdirend = "nb" + str(numbeg) + "ne" + str(numend) + "/"
                sdir = "../../../../../results/smoothdb/FVMpics/1diffmdx/" + wdirord + "/" +str(diffw)+ "/" + sdirend
                if not os.path.exists(sdir):
                        os.makedirs(sdir)
                     
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
                            u.append(float(row[ordsubtup[ip][0]])) #5 for FDcent, 6 otherwise 
                            beta = float(row[ordsubtup[ip][1]]) #could be 8 as well for o2
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
              
                s = str(dx) 
                plot(x,h ,label=s)
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
            s = "Dam Break: " + wdirord + " diff = " + str(beta)
            title(s)
            s = sdir +str(l)+".png"       
            savefig(s, bbox_inches='tight')
            legend()
            s = sdir +str(l)+ "leg.png"       
            savefig(s, bbox_inches='tight')        
            clf()
            
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
