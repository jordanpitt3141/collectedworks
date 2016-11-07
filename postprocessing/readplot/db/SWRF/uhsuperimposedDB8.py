# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 16:53:37 2015

@author: jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog, xticks,yticks
from matplotlib2tikz import save as tikz_save
import os
from subprocess import call


diffs = [10.0]
wdirbase = "../../../../../data/raw/DBASPECTRAT/o3/10/10/"
sdirbase = "../../../../../data/postprocessing/uhcomp/DBS/o3/10/"

hf = 8.0


for diff in diffs:

    #full crop front front back middez middle zz
    ylims = [[-3,9],[-3,8.5],[3.1,3.8],[-3,8.5],[2,6],[3,3.8]]
    xlims = [[-600-1,1400+1],[500-1,1400+1],[500-1,1100+1],[1150-1,1400+1],[1100-1,1150+1],[1110-1,1120+1]]
    
    
    gaps = [20,15,10,10,2,2]
    
    
    
    
    #read the DBSW file
    #for l in range(len(ylims)):
    for l in range(4,5):
        
        cylim = ylims[l]
        cxlim = xlims[l]
        
        sdir = sdirbase +str(hf)+ "/"
        if not os.path.exists(sdir):
            os.makedirs(sdir)
        
        sdirf = sdir +str(l) + "/"
        if not os.path.exists(sdirf):
            os.makedirs(sdirf)
        
        
        gap = gaps[l]
        
        wdir = wdirbase +str(hf)+ "/"


        s =  wdir + "outlast.txt"
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
                    x.append(float(row[4]))
                    h.append(float(row[5]))
                    u.append(float(row[7]))
                 j = j + 1
                 
        igap = int(gap)
        xbeg = int((cxlim[0]+600)/dx)
        xend = int((cxlim[1]+600)/dx)
        if(xbeg < 0):
            xbeg = 0
        if(xend > len(x)):
            xend = len(x)
        xt = x[xbeg:xend:igap]
        ht = h[xbeg:xend:igap]
        ut = u[xbeg:xend:igap]
        xDB = array(xt)
        hDB = array(ht)  
        uDB = array(ut) -2.686249079425793747541513835072039548079969821372709357108

        
        
        
        
        
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
        s = "Dam Break: diff = " + str(beta)
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
            + "xtick={-600,-200,200,600,1000,1400}, \n"\
            + "ytick={-3,-2,-1,0.0,1,2,3,4,5,6,7,8,9}, \n"\
            + "extra y ticks = {-2.5,-1.5,-0.5,0.5,1.5,2.5,3.43004,3.5,4.5,5.5,6.5,7.5,8.5}, \n"\
            + "scaled y ticks=false, \n clip mode=individual,\n  "\
            + "xmin=-600, \n xmax=1400, \n ymin = -3.0, \n  ymax =9,\n"\
            + "xlabel=$x$ ($m$), \n ylabel=value]\n"\
            + "\\addplot [blue] table {hDBSW.dat}; \n"\
            + "\\addplot [green!80!black] table {hDB.dat}; \n"\
            + "\\end{axis} \n \\end{tikzpicture} \n \\end{document} \n "
            #print(s)  
        
        elif(l==1): 
            s = "\\documentclass[]{article} \n\\usepackage{pgfplots} \n\\usepgfplotslibrary{external} \n\\tikzexternalize \n" \
            + "\\usepackage{tikz} \n \\usepackage{amsmath} \n  \\usepackage{pgfplots} \n \\usetikzlibrary{calc} \n" \
            + "\\pgfplotsset{compat = newest,every x tick label/.append style={font=\scriptsize},every y tick label/.append style={font=\\scriptsize,color=white}, every axis plot post/.style={line join=round}} \n" \
            + "\\begin{document} \n \\begin{tikzpicture} \n \\begin{axis}[ \n ylabel near ticks,\n xlabel near ticks, \n" \
            + "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5,}, \n"\
            + 	"extra y tick style={yticklabel style={font=\scriptsize,color=black}}, \n" \
            + "xtick={400,600,800,1000,1200,1400}, \n"\
            + "ytick= {-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5}, \n"\
            + "extra y ticks ={-3,-2,-1,0.0,1,2,3,3.43004,4,5,6,7,8,9}, \n"\
            + "scaled y ticks=false, \n clip mode=individual,\n  "\
            + "xmin=500, \n xmax=1400, \n ymin = -3.0, \n  ymax =9,\n"\
            + "xlabel=$x$ ($m$), \n ylabel=value]\n"\
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
            + "xtick={250,300,350,400,450,500,550,600}, \n"\
            + "ytick={3.15,3.25,3.35,3.45,3.55,3.65,3.75,3.85}, \n"\
            + "extra y ticks = {3.1,3.2,3.3,3.4,3.43004,3.5,3.6,3.7,3.8}, \n" \
            + "scaled y ticks=false, \n clip mode=individual,\n  "\
            + "xmin=500, \n xmax=1100, \n ymin = 3.1, \n  ymax = 3.8,\n"\
            + "xlabel=$x$ ($m$), \n ylabel=value]\n"\
            + "\\addplot [blue] table {hDBSW.dat}; \n"\
            + "\\addplot [green!80!black] table {hDB.dat}; \n"\
            + "\\end{axis} \n \\end{tikzpicture} \n \\end{document} \n "
            #print(s)
        elif(l==3):     
            s = "\\documentclass[]{article} \n\\usepackage{pgfplots} \n\\usepgfplotslibrary{external} \n\\tikzexternalize \n" \
            + "\\usepackage{tikz} \n \\usepackage{amsmath} \n  \\usepackage{pgfplots} \n \\usetikzlibrary{calc} \n" \
            + "\\pgfplotsset{compat = newest,every x tick label/.append style={font=\scriptsize},every y tick label/.append style={font=\\scriptsize,color=white}, every axis plot post/.style={line join=round}} \n" \
            + "\\begin{document} \n \\begin{tikzpicture} \n \\begin{axis}[ \n ylabel near ticks,\n xlabel near ticks, \n" \
            + "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5,}, \n"\
            + 	"extra y tick style={yticklabel style={font=\scriptsize,color=black}}, \n" \
            + "xtick={1150,1200,1250,1300,1350,1400}, \n"\
            + "ytick= {-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5}, \n"\
            + "extra y ticks ={-3,-2,-1,0.0,1,2,3,3.43004,4,5,6,7,8,9}, \n"\
            + "scaled y ticks=false, \n clip mode=individual,\n  "\
            + "xmin=1150, \n xmax=1400, \n ymin = -3.0, \n  ymax =9,\n"\
            + "xlabel=$x$ ($m$), \n ylabel=value]\n"\
            + "\\addplot [blue] table {hDBSW.dat}; \n"\
            + "\\addplot [green!80!black] table {hDB.dat}; \n"\
            + "\\end{axis} \n \\end{tikzpicture} \n \\end{document} \n "
        elif(l==4):                
            s = "\\documentclass[]{article} \n\\usepackage{pgfplots} \n\\usepgfplotslibrary{external} \n\\tikzexternalize \n" \
            + "\\usepackage{tikz} \n \\usepackage{amsmath} \n  \\usepackage{pgfplots} \n \\usetikzlibrary{calc} \n" \
            + "\\pgfplotsset{compat = newest,every x tick label/.append style={font=\scriptsize},every y tick label/.append style={font=\\scriptsize,color=white}, every axis plot post/.style={line join=round}} \n" \
            + "\\begin{document} \n \\begin{tikzpicture} \n \\begin{axis}[ \n ylabel near ticks,\n xlabel near ticks, \n" \
            + "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5,}, \n"\
            + 	"extra y tick style={yticklabel style={font=\scriptsize,color=black}}, \n" \
            + "xtick={1100,1110,1120,1130,1140,1150}, \n"\
            + "ytick={1.5,2,2.5,3,3.5,4,4.5,5,5.5,6}, \n"\
            + "extra y ticks = {1.75,2.25,2.75,3.25,3.43004,3.75,4.25,4.75,5.25,5.75}, \n" \
            + "scaled y ticks=false, \n clip mode=individual,\n  "\
            + "xmin=1100, \n xmax=1150, \n ymin = 1.75, \n  ymax = 5.75,\n"\
            + "xlabel=$x$ ($m$), \n ylabel=value]\n"\
            + "\\addplot [blue] table {hDBSW.dat}; \n"\
            + "\\addplot [green!80!black] table {hDB.dat}; \n"\
            + "\\addplot [black, dashed] coordinates {(1112,1) (1112,6)} \n;" \
            + "\\end{axis} \n \\end{tikzpicture} \n \\end{document} \n "
        else:                
            s = "\\documentclass[]{article} \n\\usepackage{pgfplots} \n\\usepgfplotslibrary{external} \n\\tikzexternalize \n" \
            + "\\usepackage{tikz} \n \\usepackage{amsmath} \n  \\usepackage{pgfplots} \n \\usetikzlibrary{calc} \n" \
            + "\\pgfplotsset{compat = newest,every x tick label/.append style={font=\scriptsize,color=white},every y tick label/.append style={font=\\scriptsize,color=white}, every axis plot post/.style={line join=round}} \n" \
            + "\\begin{document} \n \\begin{tikzpicture} \n \\begin{axis}[ \n ylabel near ticks,\n xlabel near ticks, \n" \
            + "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5,}, \n"\
            + "xticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5,}, \n"\
            + 	"extra y tick style={yticklabel style={font=\scriptsize,color=black}}, \n" \
            + 	"extra x tick style={xticklabel style={font=\scriptsize,color=black}}, \n" \
            + "xtick={1110.5,1111.5,1112.5,1113.5,1114.5,1115.5,1116.5,1117.5,1118.5,1119.5}, \n"\
            + "extra x ticks = {1110,1111,1112,1113,1114,1115,1116,1117,1118,1119,1120}, \n" \
            + "ytick={3.05,3.15,3.25,3.35,3.45,3.55,3.65,3.75}, \n"\
            + "extra y ticks = {3,3.1,3.2,3.3,3.4,3.43004,3.5,3.6,3.7,3.8}, \n" \
            + "scaled y ticks=false, \n clip mode=individual,\n  "\
            + "xmin=1110, \n xmax=1120, \n ymin = 3.0, \n  ymax = 3.8,\n"\
            + "xlabel=$x$ ($m$), \n ylabel=value]\n"\
            + "\\addplot [blue] table {hDBSW.dat}; \n"\
            + "\\addplot [green!80!black] table {hDB.dat}; \n"\
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