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
sdirbase = "../../../../../data/postprocessing/RLTPRes/o3/9/300s/"

ylims = [[1.0,1.8], [1,1.8],[1.35,1.4],[1.355,1.385],[1,1.8]]
xlims = [[-900,1800], [900,1800],[-400,700],[760,880],[1650,1800]]

gaps = [10,5,5,1,2]
timefix = 300

for l in range(len(ylims)):   
#for l in range(2,3):  
    
    cylim = ylims[l]
    cxlim = xlims[l]
            
    sdirf = sdirbase +str(l) + "/"
    sdir = sdirf
    if not os.path.exists(sdirf):
        os.makedirs(sdirf)
    
    
    gap = gaps[l]
    
    wdir = wdirbase
    
    filen = 10*timefix
    
    x2 = 500 + 1.074975295 * timefix
         
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
    xend = int((cxlim[1] + 900 )/dx) + 1
    
    if(xbeg < 0):
        xbeg = 0
    if(xend > len(x)):
        xend = len(x)
        
    xt = array(x[xbeg:xend:igap])
    ht = h[xbeg:xend:igap]
    xDBSW = array(xt)
    hDBSW = array(ht)  
    m = len(x)
    ap = 1.7399758 
    h2 = 1.36898
    u2 = 1.07498
    xf = 500 + 4.131484* timefix
    s = str(dx) 
    plot(xDBSW,hDBSW ,label=s)
    
    ylim(cylim)
    xlim(cxlim)
    #eyticks = [h2]
    #yticks(list(yticks()[0]) + eyticks)               
    xlabel("$x$ ($m$)")
    ylabel("$h$ ($m$)")
    
    n = len(xDBSW)
    s = sdirf + "h.dat"
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
        + "xtick={-900,-450,0,450,900,1350,1800}, \n"\
        + "ytick={1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8}, \n"\
        + "scaled y ticks=false, \n clip mode=individual,\n  "\
        + "xmin=" + str(cxlim[0] ) + ", \n xmax=" + str(cxlim[1]) + ", \n ymin =" + str(cylim[0] ) + ", \n  ymax = " + str(cylim[1] ) + ",\n"\
        + "xlabel=$x$ ($m$), \n ylabel=$h$ ($m$)] \n"\
        + "\\addplot [blue] table {h.dat}; \n"\
        + "\\addplot [black, dashed] coordinates {(-1000.0,1.7399758 ) (1900.0,1.7399758 )};\n"\
        + "\\end{axis} \n \\end{tikzpicture} \n \\end{document} \n "
        #print(s)  
    elif(l==1):                
        s = "\\documentclass[]{article} \n\\usepackage{pgfplots} \n\\usepgfplotslibrary{external} \n\\tikzexternalize \n" \
        + "\\usepackage{tikz} \n \\usepackage{amsmath} \n  \\usepackage{pgfplots} \n \\usetikzlibrary{calc} \n" \
        + "\\pgfplotsset{compat = newest,every x tick label/.append style={font=\scriptsize},every y tick label/.append style={font=\\scriptsize}, every axis plot post/.style={line join=round}} \n" \
        + "\\begin{document} \n \\begin{tikzpicture} \n \\begin{axis}[ \n ylabel near ticks,\n xlabel near ticks, \n" \
        + "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5,}, \n"\
        + "xtick={900,1200,1500,1800}, \n"\
        + "ytick={1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8}, \n"\
        + "scaled y ticks=false, \n clip mode=individual,\n  "\
        + "xmin=" + str(cxlim[0] ) + ", \n xmax=" + str(cxlim[1]) + ", \n ymin =" + str(cylim[0] ) + ", \n  ymax = " + str(cylim[1] ) + ",\n"\
        + "xlabel=$x$ ($m$), \n ylabel=$h$ ($m$)] \n"\
        + "\\addplot [blue] table {h.dat}; \n"\
        + "\\addplot [black, dashed] coordinates {(-1000.0,1.7399758 ) (1900.0,1.7399758)};\n"\
        + "\\end{axis} \n \\end{tikzpicture} \n \\end{document} \n "
    elif(l==2):                
        s = "\\documentclass[]{article} \n\\usepackage{pgfplots} \n\\usepgfplotslibrary{external} \n\\tikzexternalize \n" \
        + "\\usepackage{tikz} \n \\usepackage{amsmath} \n  \\usepackage{pgfplots} \n \\usetikzlibrary{calc} \n" \
        + "\\pgfplotsset{compat = newest,every x tick label/.append style={font=\scriptsize},every y tick label/.append style={font=\\scriptsize}, every axis plot post/.style={line join=round}} \n" \
        + "\\begin{document} \n \\begin{tikzpicture} \n \\begin{axis}[ \n ylabel near ticks,\n xlabel near ticks, \n" \
        + "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5,}, \n"\
        + "xtick={-400,-300,-200,-100,0,100,200,300,400,500,600,700}, \n"\
        + "extra y ticks = {1.35,1.36,1.37,1.38,1.39,1.4}, \n"\
        + "scaled y ticks=false, \n clip mode=individual,\n  "\
        + "xmin=" + str(cxlim[0] ) + ", \n xmax=" + str(cxlim[1]) + ", \n ymin =" + str(cylim[0] ) + ", \n  ymax = " + str(cylim[1] ) + ",\n"\
        + "xlabel=$x$ ($m$), \n ylabel=$h$ ($m$)] \n"\
        + "\\addplot [blue] table {h.dat}; \n"\
        + "\\addplot [black, dashed] coordinates {(-1000.0,1.7399758 ) (1900.0,1.7399758 )};\n"\
        + "\\end{axis} \n \\end{tikzpicture} \n \\end{document} \n "
        
    elif(l==3): 
        s = "\\documentclass[]{article} \n\\usepackage{pgfplots} \n\\usepgfplotslibrary{external} \n\\tikzexternalize \n" \
        + "\\usepackage{tikz} \n \\usepackage{amsmath} \n  \\usepackage{pgfplots} \n \\usetikzlibrary{calc} \n" \
        + "\\pgfplotsset{compat = newest,every x tick label/.append style={font=\scriptsize},every y tick label/.append style={font=\\scriptsize,color=white}, every axis plot post/.style={line join=round}} \n" \
        + "\\begin{document} \n \\begin{tikzpicture} \n \\begin{axis}[ \n ylabel near ticks,\n xlabel near ticks, \n" \
        + "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5,}, \n"\
        + 	"extra y tick style={yticklabel style={font=\scriptsize,color=black}}, \n" \
        + "xtick={760,780,800,820,840,860,880}, \n"\
        + "extra y ticks = {1.355, 1.36, 1.365, 1.37, 1.375,1.38,1.385}, \n" \
        + "scaled y ticks=false, \n clip mode=individual, \n  "\
        + "xmin=" + str(cxlim[0] ) + ", \n xmax=" + str(cxlim[1]) + ", \n ymin = " + str(cylim[0] ) + ", \n  ymax = " + str(cylim[1] ) + ",\n"\
        + "xlabel=$x$ ($m$), \n ylabel=$h$ ($m$)] \n"\
        + "\\addplot [blue] table {h.dat}; \n"\
        + "\\addplot [black, dashed] coordinates {(-1000," + str(h2) + " ) (1900.0," + str(h2) + " )};\n"\
        + "\\addplot [black, dotted] coordinates {("+ str(x2) + ",0.0 ) ("+ str(x2) + ",3 )};\n"\
        + "\\end{axis} \n \\end{tikzpicture} \n \\end{document} \n "
    else:                  
        s = "\\documentclass[]{article} \n\\usepackage{pgfplots} \n\\usepgfplotslibrary{external} \n\\tikzexternalize \n" \
        + "\\usepackage{tikz} \n \\usepackage{amsmath} \n  \\usepackage{pgfplots} \n \\usetikzlibrary{calc} \n" \
        + "\\pgfplotsset{compat = newest,every x tick label/.append style={font=\scriptsize},every y tick label/.append style={font=\\scriptsize,color=white}, every axis plot post/.style={line join=round}} \n" \
        + "\\begin{document} \n \\begin{tikzpicture} \n \\begin{axis}[ \n ylabel near ticks,\n xlabel near ticks, \n" \
        + "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5,}, \n"\
        + 	"extra y tick style={yticklabel style={font=\scriptsize,color=black}}, \n" \
        + "xtick={1650,1670,1690,1710,1730,1750}, \n"\
        + "extra y ticks = {1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8}, \n" \
        + "scaled y ticks=false, \n clip mode=individual, \n  "\
        + "xmin=" + str(cxlim[0] ) + ", \n xmax=" + str(cxlim[1]) + ", \n ymin = " + str(cylim[0] ) + ", \n  ymax = " + str(cylim[1] ) + ",\n"\
        + "xlabel=$x$ ($m$), \n ylabel=$h$ ($m$)] \n"\
        + "\\addplot [blue] table {h.dat}; \n"\
        + "\\addplot [black, dashed] coordinates {(-1000.0,1.7399758 ) (1900.0,1.7399758 )};\n"\
        + "\\addplot [black, dotted] coordinates {("+ str(xf) + ",0.0 ) ("+ str(xf) + ",3 )};\n"\
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
    """
