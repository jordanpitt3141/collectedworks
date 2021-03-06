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

wdirbase = "../../../../../data/raw/longcontactdiscdx9diff10fileio/o3/9/0/"
sdirbase = "../../../../../data/postprocessing/RealLongTimeSEP/o3/9/all/"

ylims = [[1.0,1.8],[1.3,1.45]]
xlims = [[300,700],[450,550]]

gaps = [8,1]

#timefixes = [300,200,100,30]
timefixes = [30,100,200,300]

for l in range(len(ylims)):    
    
    for timefix in timefixes:
        x2 = 1.074975295 * timefix + 500
        cylim = ylims[l]
        cxlim = [x2 - 60,x2 + 60]
                
        sdirf = sdirbase +str(l) + "/"
        sdir = sdirf
        if not os.path.exists(sdirf):
            os.makedirs(sdirf)
        
        
        gap = gaps[l]
        
        wdir = wdirbase
        
        filen = 10*timefix
        
             
        s = wdir + "out" + str(filen)  + ".txt"
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
                    beta = float(row[7]) 
                 j = j + 1
                 
        igap = int(gap)
        
        xbeg = int((cxlim[0] + 900)/dx) -1
        xend = int((cxlim[1] + 900)/dx) + 1
        
        if(xbeg < 0):
            xbeg = 0
        if(xend > len(x)):
            xend = len(x)
            
        xt = array(x[xbeg:xend:igap])
        ht = h[xbeg:xend:igap]
        xDBSW = array(xt)
        hDBSW = array(ht)  
        m = len(x)
        ap = 1.7645579  
        h2 = 1.36898
        u2 = 1.07498
        s = str(dx) 
        plot(xDBSW,hDBSW ,label=s)
        
        ylim(cylim)
        xlim(cxlim)
        #eyticks = [h2]
        #yticks(list(yticks()[0]) + eyticks)               
        xlabel("$x$ ($m$)")
        ylabel("$h$ ($m$)")
        
        n = len(xDBSW)
        s = sdirf + "ht="+str(timefix)+"s.dat"
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

    if(l==0):                
        s = "\\documentclass[]{article} \n\\usepackage{pgfplots} \n\\usepgfplotslibrary{external} \n\\tikzexternalize \n" \
        + "\\usepackage{tikz} \n \\usepackage{amsmath} \n  \\usepackage{pgfplots} \n \\usetikzlibrary{calc} \n" \
        + "\\pgfplotsset{compat = newest,every x tick label/.append style={font=\scriptsize},every y tick label/.append style={font=\\scriptsize}, every axis plot post/.style={line join=round}} \n" \
        + "\\begin{document} \n \\begin{tikzpicture} \n \\begin{axis}[ \n ylabel near ticks,\n xlabel near ticks, \n" \
        + "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5,}, \n"\
        + "xtick={-500,0,500,1000,1500}, \n"\
        + "ytick={1.0,1.1,1.2,1.3,1.36898,1.4,1.5,1.6,1.7,1.8}, \n"\
        + "scaled y ticks=false, \n clip mode=individual,\n  "\
        + "xmin=" + str(cxlim[0] ) + ", \n xmax=" + str(cxlim[1]) + ", \n ymin =" + str(cylim[0] ) + ", \n  ymax = " + str(cylim[1] ) + ",\n"\
        + "xlabel=$x$ ($m$), \n ylabel=$h$ ($m$)] \n"\
        + "\\addplot [blue] table {ht="+str(timefixes[0])+"s.dat}; \n"\
        + "\\addplot [green!80!black] table {ht="+str(timefixes[1])+"s.dat}; \n"\
        + "\\addplot [red] table {ht="+str(timefixes[2])+"s.dat}; \n"\
        + "\\addplot [cyan!70!white] table {ht="+str(timefixes[3])+"s.dat}; \n"\
        + "\\addplot [black, dashed] coordinates {(-1000.0,1.7645579 ) (1900.0,1.7645579 )};\n"\
        + "\\end{axis} \n \\end{tikzpicture} \n \\end{document} \n "
        #print(s)  

    else:                  
        s = "\\documentclass[]{article} \n\\usepackage{pgfplots} \n\\usepgfplotslibrary{external} \n\\tikzexternalize \n" \
        + "\\usepackage{tikz} \n \\usepackage{amsmath} \n  \\usepackage{pgfplots} \n \\usetikzlibrary{calc} \n" \
        + "\\pgfplotsset{compat = newest,every x tick label/.append style={font=\scriptsize},every y tick label/.append style={font=\\scriptsize,color=white}, every axis plot post/.style={line join=round}} \n" \
        + "\\begin{document} \n \\begin{tikzpicture} \n \\begin{axis}[ \n ylabel near ticks,\n xlabel near ticks, \n" \
        + "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5,}, \n"\
        + 	"extra y tick style={yticklabel style={font=\scriptsize,color=black}}, \n" \
        + "xtick={450,460,470,480,490,500,510,520,530,540,550}, \n"\
        + "extra y ticks = {1.3,1.325,1.35,1.375,1.4,1.425,1.45}, \n" \
        + "scaled y ticks=false, \n clip mode=individual, \n  "\
        + "xmin=" + str(cxlim[0] ) + ", \n xmax=" + str(cxlim[1]) + ", \n ymin = " + str(cylim[0] ) + ", \n  ymax = " + str(cylim[1] ) + ",\n"\
        + "xlabel=$x$ ($m$), \n ylabel=$h$ ($m$)] \n"\
        + "\\addplot [blue] table {ht="+str(timefixes[0])+"s.dat}; \n"\
        + "\\addplot [green!80!black] table {ht="+str(timefixes[1])+"s.dat}; \n"\
        + "\\addplot [red] table {ht="+str(timefixes[2])+"s.dat}; \n"\
        + "\\addplot [cyan!70!white] table {ht="+str(timefixes[3])+"s.dat}; \n"\
        + "\\addplot [black, dashed] coordinates {(-1000," + str(h2) + " ) (1900.0," + str(h2) + " )};\n"\
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
