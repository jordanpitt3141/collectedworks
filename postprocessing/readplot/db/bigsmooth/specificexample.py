# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 09:50:51 2015

@author: jordan
"""

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

#wdirords = ["o3","o2","o1"]
#ordsubtup = [[6,7],[6,8], [6,7]]

wdirords = ["FDcent","grim"]
ordsubtup = [[5,6],[5,6]]

#wdirords = ["o1"]
#ordsubtup = [[6,7]]


#diffws = [6,9,10,12]
diffws = [12]

for ip in range(len(wdirords)):
    for jp in diffws:
        
        diffw = str(jp)
        wdirord = wdirords[ip]

        #reverse
        k = 5

        ylims = [[0.0,2],[1.0,1.8],[1.3,1.45],[1.3,1.45],[1.3,1.45]]
        xlims = [[0,1000],[300,700],[500,560],[520,540],[528,536]]        
        gaps = [1, 1, 1, 1 ,1]
        zoom = [True,True,True,False,False]
        zoomints = [[350,650],[350,650],[350,650],[350,650],[350,650]]
        zoomgap = [1, 1, 1,1,1]  
        
        
        for l in range(len(ylims)):
        #for l in range(1):
            
            cylim = ylims[l]
            cxlim = xlims[l]
            
            gap = gaps[l]
            #gap = gaps[l][numbeg - k]
            
            wdir = "../../../../../data/raw/Joebigsmooth/"  +wdirord +"/" + str(k) + "/" + diffw + "/"
            sdirend = "nbs" + str(k) + "/"
            sdir = "../../../../../data/postprocessing/smoothdball/5/Joebigsmoothspec/" + wdirord + "/" +str(diffw)+ "/" + sdirend
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
                     print(l,k,zend,xend,igap)
                     zoomgapi = int(zoomgap[l])
                     xt = x[xbeg:zbeg:igap] + x[zbeg:zend:zoomgapi] + x[zend:xend:igap]
                     ht = h[xbeg:zbeg:igap] + h[zbeg:zend:zoomgapi] + h[zend:xend:igap]
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
