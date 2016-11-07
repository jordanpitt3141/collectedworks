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


diffs = [0.5,2.5,4.5,8.5,10.0]
wdirbase = "../../../../../data/raw/DBSWanew100s/o3/10/"
sdirbase = "../../../../../data/postprocessing/DBSWfrontsingular/o3/10/"


for diff in diffs:

    #full crop front front back middez middle zz
    ylims = [[1.0,1.8],[1.0,1.8],[1,1.8],[1.0,1.8],[1.32,1.42],[1.32,1.42]]
    xlims = [[0-1,1000+1],[200-1,1000+1],[770-1,920+1],[600-1,950+1],[550-1,650+1],[590-1,620+1]]
    
    
    gaps = [10,10,10,10,5,5]
    
    
    
    
    #read the DBSW file
    for l in range(len(ylims)):
    #for l in range(1):
        
        cylim = ylims[l]
        cxlim = xlims[l]
        
        sdir = sdirbase +str(diff)+ "/"
        if not os.path.exists(sdir):
            os.makedirs(sdir)
        
        sdirf = sdir +str(l) + "/"
        if not os.path.exists(sdirf):
            os.makedirs(sdirf)
        
        
        gap = gaps[l]
        
        wdir = wdirbase +str(diff)+ "/"
             
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
                    x.append(float(row[4])+200)
                    h.append(float(row[5]))
                    u.append(float(row[6])) 
                    beta = float(row[7]) 
                 j = j + 1
                 
        igap = int(gap)
        xbeg = int((cxlim[0]-200)/dx)
        xend = int((cxlim[1]-200)/dx)
        if(xbeg < 0):
            xbeg = 0
        if(xend > len(x)):
            xend = len(x)
        xt = x[xbeg:xend:igap]
        ht = h[xbeg:xend:igap]
        xDBSW = array(xt)
        hDBSW = array(ht)  
        m = len(x)
        ap = 1.736397786   
        #h2 = 1.36898
        const = ap*ones(m)
        s = str(dx) 
        plot(xDBSW,hDBSW ,label=s)
        plot(x,const,"--k" ,label="a ref")
        
        ylim(cylim)
        xlim(cxlim)
        #eyticks = [h2]
        #yticks(list(yticks()[0]) + eyticks)               
        xlabel("$x$ ($m$)")
        ylabel("$h$ ($m$)")
        
        n = len(xDBSW)
        s = sdirf + "hDBSW.dat"
        with open(s,'w') as file1:
            for i in range(n):
                s ="%3.8f%5s%1.15f\n" %(xDBSW[i]," ",hDBSW[i])
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
            + "\\pgfplotsset{compat = newest,every x tick label/.append style={font=\scriptsize},every y tick label/.append style={font=\\scriptsize}, every axis plot post/.style={line join=round}} \n" \
            + "\\begin{document} \n \\begin{tikzpicture} \n \\begin{axis}[ \n ylabel near ticks,\n xlabel near ticks, \n" \
            + "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5,}, \n"\
            + "xtick={0,200,400,600,800,1000}, \n"\
            + "ytick={1.0,1.1,1.2,1.3,1.36898,1.4,1.5,1.6,1.7,1.8}, \n"\
            + "scaled y ticks=false, \n clip mode=individual,\n  "\
            + "xmin=0, \n xmax=1000, \n ymin = 1.0, \n  ymax = 1.8,\n"\
            + "xlabel=$x$ ($m$), \n ylabel=$h$ ($m$)]"\
            + "\\addplot [blue] table {hDBSW.dat}; \n"\
            + "\\addplot [black, dashed] coordinates {(-100.0,1.736397786) (1100.0,1.736397786)};\n"\
            + "\\end{axis} \n \\end{tikzpicture} \n \\end{document} \n "
            #print(s)  
        
        elif(l==1): 
            s = "\\documentclass[]{article} \n\\usepackage{pgfplots} \n\\usepgfplotslibrary{external} \n\\tikzexternalize \n" \
            + "\\usepackage{tikz} \n \\usepackage{amsmath} \n  \\usepackage{pgfplots} \n \\usetikzlibrary{calc} \n" \
            + "\\pgfplotsset{compat = newest,every x tick label/.append style={font=\scriptsize},every y tick label/.append style={font=\\scriptsize}, every axis plot post/.style={line join=round}} \n" \
            + "\\begin{document} \n \\begin{tikzpicture} \n \\begin{axis}[ \n ylabel near ticks,\n xlabel near ticks, \n" \
            + "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5,}, \n"\
            + "xtick={0,200,400,600,800,1000}, \n"\
            + "ytick={1.0,1.1,1.2,1.3,1.36898,1.4,1.5,1.6,1.7,1.8}, \n"\
            + "scaled y ticks=false, \n clip mode=individual,\n  "\
            + "xmin=200, \n xmax=1000, \n ymin = 1.0, \n  ymax = 1.8,\n"\
            + "xlabel=$x$ ($m$), \n ylabel=$h$ ($m$)]"\
            + "\\addplot [blue] table {hDBSW.dat}; \n"\
            + "\\addplot [black, dashed] coordinates {(-100.0,1.736397786) (1100.0,1.736397786)};\n"\
            + "\\end{axis} \n \\end{tikzpicture} \n \\end{document} \n "               
            #print(s)
        elif(l==2):
            s = "\\documentclass[]{article} \n\\usepackage{pgfplots} \n\\usepgfplotslibrary{external} \n\\tikzexternalize \n" \
            + "\\usepackage{tikz} \n \\usepackage{amsmath} \n  \\usepackage{pgfplots} \n \\usetikzlibrary{calc} \n" \
            + "\\pgfplotsset{compat = newest,every x tick label/.append style={font=\scriptsize},every y tick label/.append style={font=\\scriptsize}, every axis plot post/.style={line join=round}} \n" \
            + "\\begin{document} \n \\begin{tikzpicture} \n \\begin{axis}[ \n ylabel near ticks,\n xlabel near ticks, \n" \
            + "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5,}, \n"\
            + "xtick={770,800,830,860,890,920}, \n"\
            + "ytick={1.0,1.1,1.2,1.3,1.36898,1.4,1.5,1.6,1.7,1.8}, \n"\
            + "scaled y ticks=false, \n clip mode=individual,\n  "\
            + "xmin=770, \n xmax=920, \n ymin = 1.0, \n  ymax = 1.8,\n"\
            + "xlabel=$x$ ($m$), \n ylabel=$h$ ($m$)]"\
            + "\\addplot [blue] table {hDBSW.dat}; \n"\
            + "\\addplot [black, dashed] coordinates {(-100.0,1.736397786) (1100.0,1.736397786)};\n"\
            + "\\end{axis} \n \\end{tikzpicture} \n \\end{document} \n "
            #print(s)
        elif(l==3):                
            s = "\\documentclass[]{article} \n\\usepackage{pgfplots} \n\\usepgfplotslibrary{external} \n\\tikzexternalize \n" \
            + "\\usepackage{tikz} \n \\usepackage{amsmath} \n  \\usepackage{pgfplots} \n \\usetikzlibrary{calc} \n" \
            + "\\pgfplotsset{compat = newest,every x tick label/.append style={font=\scriptsize},every y tick label/.append style={font=\\scriptsize}, every axis plot post/.style={line join=round}} \n" \
            + "\\begin{document} \n \\begin{tikzpicture} \n \\begin{axis}[ \n ylabel near ticks,\n xlabel near ticks, \n" \
            + "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5,}, \n"\
            + "xtick={600,650,700,750,800,850,900,950}, \n"\
            + "ytick={1.0,1.1,1.2,1.3,1.36898,1.4,1.5,1.6,1.7,1.8}, \n"\
            + "scaled y ticks=false, \n clip mode=individual,\n  "\
            + "xmin=600, \n xmax=950, \n ymin = 1.0, \n  ymax = 1.8,\n"\
            + "xlabel=$x$ ($m$), \n ylabel=$h$ ($m$)]"\
            + "\\addplot [blue] table {hDBSW.dat}; \n"\
            + "\\addplot [black, dashed] coordinates {(-100.0,1.736397786) (1100.0,1.736397786)};\n"\
            + "\\end{axis} \n \\end{tikzpicture} \n \\end{document} \n "
        elif(l==4):                
            s = "\\documentclass[]{article} \n\\usepackage{pgfplots} \n\\usepgfplotslibrary{external} \n\\tikzexternalize \n" \
            + "\\usepackage{tikz} \n \\usepackage{amsmath} \n  \\usepackage{pgfplots} \n \\usetikzlibrary{calc} \n" \
            + "\\pgfplotsset{compat = newest,every x tick label/.append style={font=\scriptsize},every y tick label/.append style={font=\\scriptsize,color=white}, every axis plot post/.style={line join=round}} \n" \
            + "\\begin{document} \n \\begin{tikzpicture} \n \\begin{axis}[ \n ylabel near ticks,\n xlabel near ticks, \n" \
            + "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5,}, \n"\
            + 	"extra y tick style={yticklabel style={font=\scriptsize,color=black}}, \n" \
            + "xtick={550,560,570,580,590,600,610,620,630,640,650}, \n"\
            + "ytick={1.33,1.35,1.37,1.39,1.41,1.43}, \n"\
            + "extra y ticks = {1.32,1.34,1.36,1.36898,1.38,1.40,1.42}, \n" \
            + "scaled y ticks=false, \n clip mode=individual,\n  "\
            + "xmin=550, \n xmax=650, \n ymin = 1.32, \n  ymax = 1.42,\n"\
            + "xlabel=$x$ ($m$), \n ylabel=$h$ ($m$)]"\
            + "\\addplot [blue] table {hDBSW.dat}; \n"\
            + "\\addplot [black, dashed] coordinates {(-100.0,1.736397786) (1100.0,1.736397786)};\n"\
            + "\\end{axis} \n \\end{tikzpicture} \n \\end{document} \n "
        else:                
            s = "\\documentclass[]{article} \n\\usepackage{pgfplots} \n\\usepgfplotslibrary{external} \n\\tikzexternalize \n" \
            + "\\usepackage{tikz} \n \\usepackage{amsmath} \n  \\usepackage{pgfplots} \n \\usetikzlibrary{calc} \n" \
            + "\\pgfplotsset{compat = newest,every x tick label/.append style={font=\scriptsize},every y tick label/.append style={font=\\scriptsize,color=white}, every axis plot post/.style={line join=round}} \n" \
            + "\\begin{document} \n \\begin{tikzpicture} \n \\begin{axis}[ \n ylabel near ticks,\n xlabel near ticks, \n" \
            + "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5,}, \n"\
            + 	"extra y tick style={yticklabel style={font=\scriptsize,color=black}}, \n" \
            + "xtick={590,595,600,605,610,615,620}, \n"\
            + "ytick={1.33,1.35,1.37,1.39,1.41,1.43}, \n"\
            + "extra y ticks = {1.32,1.34,1.36,1.36898,1.38,1.40,1.42}, \n" \
            + "scaled y ticks=false, \n clip mode=individual,\n  "\
            + "xmin=590, \n xmax=620, \n ymin = 1.32, \n  ymax = 1.42,\n"\
            + "xlabel=$x$ ($m$), \n ylabel=$h$ ($m$)]"\
            + "\\addplot [blue] table {hDBSW.dat}; \n"\
            + "\\addplot [black, dashed] coordinates {(-100.0,1.736397786) (1100.0,1.736397786)};\n"\
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