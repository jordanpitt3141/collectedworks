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

hf = 2.0


for diff in diffs:

    #full crop front front back middez middle zz
    ylims = [[0.0,2.4],[0.0,2.4],[1.38,1.52],[0.0,2.4],[1.41,1.5],[1.41,1.5]]
    xlims = [[-600-1,1400+1],[0-1,1000+1],[250-1,600+1],[650-1,950+1],[600-1,650+1],[620-1,640+1]]
    
    
    gaps = [20,15,10,10,2,2]
    
    
    
    
    #read the DBSW file
    for l in range(len(ylims)):
    #for l in range(1):
        
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
        uDB = array(ut) + 0.14800392876550815771841437274865638707859757004709966416

        
        
        
        
        
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
            + "ytick={0.2,0.6,1.0,1.4,1.8,2.2}, \n"\
            + "extra y ticks = {0.0,0.4,0.8,1.2,1.45384,1.6,2.0,2.4}, \n"\
            + "scaled y ticks=false, \n clip mode=individual,\n  "\
            + "xmin=-600, \n xmax=1400, \n ymin = 0.0, \n  ymax = 2.4,\n"\
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
            + "xtick={0,200,400,600,800,1000}, \n"\
            + "ytick={1.1,1.3,1.5,1.7,1.9,2.1,2.3}, \n"\
            + "extra y ticks = {1.0,1.2,1.4,1.45384,1.6,1.8,2.0,2.2,2.4}, \n"\
            + "scaled y ticks=false, \n clip mode=individual,\n  "\
            + "xmin=0.0, \n xmax=1000, \n ymin = 1.0, \n  ymax = 2.4,\n"\
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
            + "ytick={1.39,1.41,1.43,1.45,1.47,1.49,1.51}, \n"\
            + "extra y ticks = {1.38,1.4,1.42,1.44,1.45384,1.46,1.48,1.5,1.52}, \n" \
            + "scaled y ticks=false, \n clip mode=individual,\n  "\
            + "xmin=250, \n xmax=600, \n ymin = 1.38, \n  ymax = 1.52,\n"\
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
            + "xtick={650,700,750,800,850,900,950}, \n"\
            + "ytick={0.2,0.6,1.0,1.4,1.8,2.2}, \n"\
            + "extra y ticks = {0.0,0.4,0.8,1.2,1.45384,1.6,2.0,2.4}, \n"\
            + "scaled y ticks=false, \n clip mode=individual,\n  "\
            + "xmin=650, \n xmax=950, \n ymin = 0.0, \n  ymax = 2.4,\n"\
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
            + "xtick={600,610,620,630,640,650}, \n"\
            + "ytick={1.41,1.43,1.45,1.47,1.49}, \n"\
            + "extra y ticks = {1.4,1.42,1.44,1.45384,1.46,1.48,1.5}, \n" \
            + "scaled y ticks=false, \n clip mode=individual,\n  "\
            + "xmin=600, \n xmax=650, \n ymin = 1.40, \n  ymax = 1.5,\n"\
            + "xlabel=$x$ ($m$), \n ylabel=value]\n"\
            + "\\addplot [blue] table {hDBSW.dat}; \n"\
            + "\\addplot [green!80!black] table {hDB.dat}; \n"\
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
            + "xtick={621,623,625,627,629,631,633,635,637,639}, \n"\
            + "extra x ticks = {620,622,624,626,628,630,632,634,636,638,640}, \n" \
            + "ytick={1.41,1.43,1.45,1.47,1.49}, \n"\
            + "extra y ticks = {1.4,1.42,1.44,1.45384,1.46,1.48,1.5}, \n" \
            + "scaled y ticks=false, \n clip mode=individual,\n  "\
            + "xmin=620, \n xmax=640, \n ymin = 1.40, \n  ymax = 1.5,\n"\
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