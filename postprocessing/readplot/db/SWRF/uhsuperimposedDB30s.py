# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 16:53:37 2015

@author: jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog, xticks,yticks
import os
from subprocess import call


diffs = [10.0]
sdirbase = "../../../../../data/postprocessing/uhcomp30s/DB/o3/"
cwdirbase = "../../../../../data/raw/longcontactdiscdx9diff10fileio/o3/9/0/"


for diff in diffs:

    #full crop front front back middez middle zz
    ylims = [[0.2,2.2],[1.25,1.5],[0.2,2.2] ,[1.31,1.43]]
    xlims = [[350,650],[425,525],[550,625],[528,536]]
    
    
    gaps = [2,2,2,1]
    
    
    
    
    #read the DBSW file
    for l in range(len(ylims)):
    #for l in range(3):
        
        cylim = ylims[l]
        cxlim = xlims[l]
        
        sdir = sdirbase +str(diff)+ "/"
        if not os.path.exists(sdir):
            os.makedirs(sdir)
        
        sdirf = sdir +str(l) + "/"
        if not os.path.exists(sdirf):
            os.makedirs(sdirf)
        
        
        gap = gaps[l]

        h2 = 1.36898
        u2 = 1.074975
        x2 = 500 + u2*30

        s = cwdirbase + "out300.txt"
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
                 j = j + 1
                 
        igap = int(gap)
        xbeg = int((cxlim[0] + 900)/dx) - 1
        xend = int((cxlim[1] + 900)/dx) + 1
        if(xbeg < 0):
            xbeg = 0
        if(xend > len(x)):
            xend = len(x)
        xt = x[xbeg:xend:igap]
        ht = h[xbeg:xend:igap]
        ut = u[xbeg:xend:igap]
        xDB = array(xt)
        hDB = array(ht)  
        uDB = array(ut) + (h2 - u2)
        

        m = len(x)
        s = str(dx) 
        plot(xDB,hDB,label=s)
        plot(xDB,uDB,label=s)

        
        ylim(cylim)
        xlim(cxlim)
        #eyticks = [h2]
        #yticks(list(yticks()[0]) + eyticks)               
        xlabel("$x$ ($m$)")
        ylabel("$h$ ($m$)")
       
        n = len(xDB)
        s = sdirf + "hDBSW.dat"
        with open(s,'w') as file1:
            for i in range(n):
                s ="%3.8f%5s%1.15f\n" %(xDB[i]," ",hDB[i])
                file1.write(s)
                
        n = len(xDB)
        s = sdirf + "hDB.dat"
        with open(s,'w') as file1:
            for i in range(n):
                s ="%3.8f%5s%1.15f\n" %(xDB[i]," ",uDB[i])
                file1.write(s)
        s = "Dam Break: diff = "
        title(s)
        s = sdir +str(l)+".png"       
        savefig(s, bbox_inches='tight')
        legend()
        s = sdir +str(l)+ "leg.png"       
        savefig(s, bbox_inches='tight')        
        clf()
        
            #legend()
        if(l==0):    
            s = "\\documentclass[]{article} \n\\usepackage{pgfplots} \n\\usepgfplotslibrary{external} \n\\tikzexternalize \n" \
            + "\\usepackage{tikz} \n \\usepackage{amsmath} \n  \\usepackage{pgfplots} \n \\usetikzlibrary{calc} \n" \
            + "\\pgfplotsset{compat = newest,every x tick label/.append style={font=\scriptsize},every y tick label/.append style={font=\\scriptsize,color=white}, every axis plot post/.style={line join=round}} \n" \
            + "\\begin{document} \n \\begin{tikzpicture} \n \\begin{axis}[ \n ylabel near ticks,\n xlabel near ticks, \n" \
            + "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5,}, \n"\
            + 	"extra y tick style={yticklabel style={font=\scriptsize,color=black}}, \n" \
            + "xtick={350,400,450,500,550,600,650}, \n"\
            + "ytick={0.4,0.8,1.2,1.6,2.0}, \n"\
            + "extra y ticks = {0.2,0.6,1.0,1.4,1.8,2.2}, \n"\
            + "scaled y ticks=false, \n clip mode=individual,\n  "\
            + "xmin=" + str(cxlim[0] ) + ", \n xmax=" + str(cxlim[1]) + ", \n ymin =" + str(cylim[0] ) + ", \n  ymax = " + str(cylim[1] ) + ",\n"\
            + "xlabel=$x$ ($m$), \n ylabel=value]"\
            + "\\addplot [blue] table {hDBSW.dat}; \n"\
            + "\\addplot [green!80!black] table {hDB.dat}; \n"\
            + "\\end{axis} \n \\end{tikzpicture} \n \\end{document} \n "
            #print(s)  
        
        elif(l==1): 
            s = "\\documentclass[]{article} \n\\usepackage{pgfplots} \n\\usepgfplotslibrary{external} \n\\tikzexternalize \n" \
            + "\\usepackage{tikz} \n \\usepackage{amsmath} \n  \\usepackage{pgfplots} \n \\usetikzlibrary{calc} \n" \
            + "\\pgfplotsset{compat = newest,every x tick label/.append style={font=\scriptsize},every y tick label/.append style={font=\\scriptsize}, every axis plot post/.style={line join=round}} \n" \
            + "\\begin{document} \n \\begin{tikzpicture} \n \\begin{axis}[ \n ylabel near ticks,\n xlabel near ticks, \n" \
            + "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5,}, \n"\
            + "xtick={425,450,475,500,525}, \n"\
            + "ytick={1.25,1.3,1.35,1.4,1.45,1.5}, \n"\
            + "scaled y ticks=false, \n clip mode=individual,\n  "\
            + "xmin=" + str(cxlim[0] ) + ", \n xmax=" + str(cxlim[1]) + ", \n ymin =" + str(cylim[0] ) + ", \n  ymax = " + str(cylim[1] ) + ",\n"\
            + "xlabel=$x$ ($m$), \n ylabel=value]"\
            + "\\addplot [blue] table {hDBSW.dat}; \n"\
            + "\\addplot [green!80!black] table {hDB.dat}; \n"\
            + "\\end{axis} \n \\end{tikzpicture} \n \\end{document} \n "               
            #print(s)
        elif(l==2):
            s = "\\documentclass[]{article} \n\\usepackage{pgfplots} \n\\usepgfplotslibrary{external} \n\\tikzexternalize \n" \
            + "\\usepackage{tikz} \n \\usepackage{amsmath} \n  \\usepackage{pgfplots} \n \\usetikzlibrary{calc} \n" \
            + "\\pgfplotsset{compat = newest,every x tick label/.append style={font=\scriptsize},every y tick label/.append style={font=\\scriptsize,color=white}, every axis plot post/.style={line join=round}} \n" \
            + "\\begin{document} \n \\begin{tikzpicture} \n \\begin{axis}[ \n ylabel near ticks,\n xlabel near ticks, \n" \
            + "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5,}, \n"\
            + 	"extra y tick style={yticklabel style={font=\scriptsize,color=black}}, \n" \
            + "xtick={550,575,600,625}, \n"\
            + "ytick={0.4,0.8,1.2,1.6,2.0}, \n"\
            + "extra y ticks = {0.2,0.6,1.0,1.4,1.8,2.2}, \n"\
            + "scaled y ticks=false, \n clip mode=individual,\n  "\
            + "xmin=" + str(cxlim[0] ) + ", \n xmax=" + str(cxlim[1]) + ", \n ymin =" + str(cylim[0] ) + ", \n  ymax = " + str(cylim[1] ) + ",\n"\
            + "xlabel=$x$ ($m$), \n ylabel=value]"\
            + "\\addplot [blue] table {hDBSW.dat}; \n"\
            + "\\addplot [green!80!black] table {hDB.dat}; \n"\
            + "\\end{axis} \n \\end{tikzpicture} \n \\end{document} \n "
            #print(s)
        else:     
            s = "\\documentclass[]{article} \n\\usepackage{pgfplots} \n\\usepgfplotslibrary{external} \n\\tikzexternalize \n" \
            + "\\usepackage{tikz} \n \\usepackage{amsmath} \n  \\usepackage{pgfplots} \n \\usetikzlibrary{calc} \n" \
            + "\\pgfplotsset{compat = newest,every x tick label/.append style={font=\scriptsize},every y tick label/.append style={font=\\scriptsize,color=white}, every axis plot post/.style={line join=round}} \n" \
            + "\\begin{document} \n \\begin{tikzpicture} \n \\begin{axis}[ \n ylabel near ticks,\n xlabel near ticks, \n" \
            + "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5,}, \n"\
            + 	"extra y tick style={yticklabel style={font=\scriptsize,color=black}}, \n" \
            + "xtick={528,529,530,531,532,533,534,535,536}, \n"\
            + "ytick={1.32,1.34,1.36,1.38,1.4,1.42}, \n"\
            + "extra y ticks = {1.31,1.33,1.35,1.37,1.39,1.41,1.43}, \n"\
            + "scaled y ticks=false, \n clip mode=individual,\n  "\
            + "xmin=" + str(cxlim[0] ) + ", \n xmax=" + str(cxlim[1]) + ", \n ymin =" + str(cylim[0] ) + ", \n  ymax = " + str(cylim[1] ) + ",\n"\
            + "xlabel=$x$ ($m$), \n ylabel=value]"\
            + "\\addplot [blue] table {hDBSW.dat}; \n"\
            + "\\addplot [green!80!black] table {hDB.dat}; \n" \
            + "\\addplot [black, dashed] coordinates {(500," + str(h2) + " ) (550," + str(h2) + " )};\n"\
            + "\\addplot [black, dotted] coordinates {("+ str(x2) + ",0.0 ) ("+ str(x2) + ",3 )};\n"\
            + "\\end{axis} \n \\end{tikzpicture} \n \\end{document} \n "

            #print(s)
            
        filen = sdirf +str(l)+ ".tex" 
        file1 = open(filen, 'w')
        file1.write(s)
        file1.close() 

        #make makefile
        s = "LAT = pdflatex \nLATFLAGS = -shell-escape\n\n"
        newl = str(l) +":\n\t $(LAT) $(LATFLAGS) " +str(l)+".tex\n\n"
        s = s + newl
        newl = "clean:\n\t rm -f  *~ ./*.log ./*.aux ./*.auxlock ./*.dep ./*.dpth ./*.pdf ./*.gz\n\n" 
        s = s + newl
        
        newl = "all: " +str(l)
        s = s + newl
           
        filen =sdirf + "Makefile" 
        
        file1 = open(filen, 'w')
        file1.write(s)
        file1.close()  
        call(['make','-C',sdirf,'all'])