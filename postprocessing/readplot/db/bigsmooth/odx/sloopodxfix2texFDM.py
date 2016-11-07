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

diffw = "11"
wdirord = "o3"

#wdirords = ["o3","FDcent","grim","o2","o1"]
#ordsubtup = [[6,7],[5,6],[5,6], [6,8], [6,7]]

#wdirords = ["o3","o2","o1"]
#ordsubtup = [[6,7],[6,8], [6,7]]

wdirords = ["FDcent"]
ordsubtup = [[5,6]]

#wdirords = ["o3"]
#ordsubtup = [[6,7]]

#wdirords = ["o3","FDcent","grim","o2","o1"]
#ordsubtup = [[6,7],[5,6],[5,6], [6,8], [6,7]]

diffs = [0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1.0,2.5,5.0,7.5,10.0,25.0,50.0,75.0,100.0,250.0,500.0,750.0,1000.0]


diffws = [12]
#diffws = [6]

for ip in range(len(wdirords)):
    for jp in diffws:
        
        diffw = str(jp)
        wdirord = wdirords[ip]

        #reverse
        numbeg = 10
        numend = 9
        incr = -1
        nums = range(numbeg,numend,incr)

        ylims = [[0.0,2],[1.0,1.8],[1.3,1.45],[1.3,1.45],[1.3,1.45]]
        xlims = [[0,1000],[299,701],[499,561],[519,541],[527,537]]
        #gaps = [[2,4,6,8,10,12,14],[1,1,2,6,10,14,14],[1,1,1,1,1,1,1],[1,1,1,1,1,1,1],[1,1,1,1,1,1,1]]
        #zoom = [True,True,False,False,False]
        #zoomints = [[500,560],[500,560],[500,560],[500,560],[530,540]]
        #zoomgap = [2,1,1,1,1] 
        
        gaps = [[8*2**6,8*2**5,8*2**4,8*2**3,8*2**2,8*2,8], \
                [2*2**5,2*2**4,2*2**3,2*2**2,2*2,2,2], \
                [4,4,2,2,1,1,1], \
                [2,2,1,1,1,1,1], \
                [1,1,1,1,1,1,1]]
        zoom = [True,True,False,False,False]
        zoomints = [[350,650],[350,650],[350,650],[350,650],[350,650]]
        zoomgap = [[2**6,2**5,2**4,2**3,4,2,1], \
                    [2**5,2**4,2**3,2**2,2,1,1], \
                    [1,1,1,1,1,1,1],[1,1,1,1,1,1,1],[1,1,1,1,1,1,1]]  
        
        
        for l in range(len(ylims)):
        #for l in range(1):
            
            cylim = ylims[l]
            cxlim = xlims[l]
            
            sdirend = "nb" + str(numbeg) + "ne" + str(numend) + "/"
            sdir = "../../../../../../data/postprocessing/smoothdball/paperpicsN/" + wdirord + "/" +str(diffw)+ "/" + sdirend
            if not os.path.exists(sdir):
                os.makedirs(sdir)
            
            sdirf = sdir +str(l) + "/"
            if not os.path.exists(sdirf):
                os.makedirs(sdirf)
            
            
            for k in nums:
                #gap = gaps[l][k-numbeg]
                gap = gaps[l][numbeg - k]
                
                wdir = "../../../../../../data/raw/Joebigsmooth/"  +wdirord +"/" + str(k) + "/" + diffw + "/"
                     
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
                         #zoomgapi = int(zoomgap[l][k-numbeg])
                         zoomgapi = int(zoomgap[l][numbeg - k])
                         xt = x[xbeg:zbeg:igap] + x[zbeg:zend:zoomgapi] + x[zend:xend:igap]
                         ht = h[xbeg:zbeg:igap] + h[zbeg:zend:zoomgapi] + h[zend:xend:igap]
                     else:
                         xbeg = int(cxlim[0]/dx)
                         xend = int(cxlim[1]/dx)
                         xt = x[xbeg:xend:igap]
                         ht = h[xbeg:xend:igap]
                     x = array(xt)
                     h = array(ht)     
                m = len(x)
                ap = 1.736397786
                #h2 = 1.36898
                const = ap*ones(m)
                s = str(dx) 
                plot(x,h ,label=s)
                plot(x,const,"--k" ,label="a ref")
                
                ylim(cylim)
                xlim(cxlim)
                #eyticks = [h2]
                #yticks(list(yticks()[0]) + eyticks)               
                xlabel("$x$ ($m$)")
                ylabel("$h$ ($m$)")
                
                n = len(xt)
                s = sdirf  + "dx" + str(k) + "diff"+ str(beta) + "h.dat"
                with open(s,'w') as file1:
                    for i in range(n):
                        s ="%3.8f%5s%1.15f\n" %(x[i]," ",h[i])
                        file1.write(s)
            s = "Dam Break: " + wdirord + " diff = " + str(beta)
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
                + "\\pgfplotsset{compat = newest,every x tick label/.append style={font=\scriptsize},every y tick label/.append style={font=\\scriptsize}, every axis plot post/.style={line join=round}}" \
                + "\\begin{document} \n \\begin{tikzpicture} \n \\begin{axis}[ \n ylabel near ticks,\n xlabel near ticks, \n" \
                + "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5,}, \n"\
                + "xtick={300,350,400,450,500,550,600,650,700}, \n"\
                + "ytick={1.0,1.1,1.2,1.3,1.36898,1.4,1.5,1.6,1.7,1.8}, \n"\
                + "scaled y ticks=false, \n clip mode=individual,\n  "\
                + "xmin=300, \n xmax=700, \n ymin = 1.0, \n  ymax = 1.8,\n"\
                + "xlabel=$x$ ($m$), \n ylabel=$h$ ($m$)]"\
                + "\\addplot [blue] table {dx4diff"+ str(beta) +"h.dat}; \n"\
                + "\\addplot [green!80!black] table {dx5diff"+ str(beta) +"h.dat}; \n"\
                + "\\addplot [red] table {dx6diff"+ str(beta) +"h.dat}; \n"\
                + "\\addplot [cyan!70!white] table {dx7diff"+ str(beta) +"h.dat}; \n"\
                + "\\addplot [violet!70!white] table {dx8diff"+ str(beta) +"h.dat}; \n"\
                + "\\addplot [yellow!70!black] table {dx9diff"+ str(beta) +"h.dat};\n"\
                + "\\addplot [black] table {dx10diff"+ str(beta) +"h.dat};\n"\
                + "\\addplot [black, dashed] coordinates {(200.0,1.736397786) (800.0,1.736397786)};\n"\
                + "\\end{axis} \n \\end{tikzpicture} \n \\end{document} \n "
                #print(s)    
            elif(l==1):                
                s = "\\documentclass[]{article} \n\\usepackage{pgfplots} \n\\usepgfplotslibrary{external} \n\\tikzexternalize \n" \
                + "\\usepackage{tikz} \n \\usepackage{amsmath} \n  \\usepackage{pgfplots} \n \\usetikzlibrary{calc} \n" \
                + "\\pgfplotsset{compat = newest,every x tick label/.append style={font=\scriptsize},every y tick label/.append style={font=\\scriptsize}, every axis plot post/.style={line join=round}}" \
                + "\\begin{document} \n \\begin{tikzpicture} \n \\begin{axis}[ \n ylabel near ticks,\n xlabel near ticks, \n" \
                + "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5,}, \n"\
                + "xtick={300,350,400,450,500,550,600,650,700}, \n"\
                + "ytick={1.0,1.1,1.2,1.3,1.36898,1.4,1.5,1.6,1.7,1.8}, \n"\
                + "scaled y ticks=false, \n clip mode=individual,\n  "\
                + "xmin=300, \n xmax=700, \n ymin = 1.0, \n  ymax = 1.8,\n"\
                + "xlabel=$x$ ($m$), \n ylabel=$h$ ($m$)]"\
                + "\\addplot [blue] table {dx4diff"+ str(beta) +"h.dat}; \n"\
                + "\\addplot [green!80!black] table {dx5diff"+ str(beta) +"h.dat}; \n"\
                + "\\addplot [red] table {dx6diff"+ str(beta) +"h.dat}; \n"\
                + "\\addplot [cyan!70!white] table {dx7diff"+ str(beta) +"h.dat}; \n"\
                + "\\addplot [violet!70!white] table {dx8diff"+ str(beta) +"h.dat}; \n"\
                + "\\addplot [yellow!70!black] table {dx9diff"+ str(beta) +"h.dat};\n"\
                + "\\addplot [black] table {dx10diff"+ str(beta) +"h.dat};\n"\
                + "\\addplot [black, dashed] coordinates {(200.0,1.736397786) (800.0,1.736397786)};\n"\
                + "\\end{axis} \n \\end{tikzpicture} \n \\end{document} \n "
                #print(s)
            elif(l==2):
                s = "\\documentclass[]{article} \n\\usepackage{pgfplots} \n\\usepgfplotslibrary{external} \n\\tikzexternalize \n" \
                + "\\usepackage{tikz} \n \\usepackage{amsmath} \n  \\usepackage{pgfplots} \n \\usetikzlibrary{calc} \n" \
                + "\\pgfplotsset{compat = newest,every x tick label/.append style={font=\scriptsize},every y tick label/.append style={font=\\scriptsize}, every axis plot post/.style={line join=round}}" \
                + "\\begin{document} \n \\begin{tikzpicture} \n \\begin{axis}[ \n ylabel near ticks,\n xlabel near ticks, \n" \
                + "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5,}, \n"\
                + "xtick={500,510,520,530,540,550,560}, \n"\
                + "ytick={1.3,1.32,1.34,1.36,1.36898,1.38,1.4,1.42,1.44}, \n"\
                + "scaled y ticks=false, \n clip mode=individual,\n  "\
                + "xmin=500, \n xmax=560, \n ymin = 1.3, \n  ymax = 1.44,\n"\
                + "xlabel=$x$ ($m$), \n ylabel=$h$ ($m$)]"\
                + "\\addplot [blue] table {dx4diff"+ str(beta) +"h.dat}; \n"\
                + "\\addplot [green!80!black] table {dx5diff"+ str(beta) +"h.dat}; \n"\
                + "\\addplot [red] table {dx6diff"+ str(beta) +"h.dat}; \n"\
                + "\\addplot [cyan!70!white] table {dx7diff"+ str(beta) +"h.dat}; \n"\
                + "\\addplot [violet!70!white] table {dx8diff"+ str(beta) +"h.dat}; \n"\
                + "\\addplot [yellow!70!black] table {dx9diff"+ str(beta) +"h.dat};\n"\
                + "\\addplot [black] table {dx10diff"+ str(beta) +"h.dat};\n"\
                + "\\addplot [black, dashed] coordinates {(200.0,1.736397786) (800.0,1.736397786)};\n"\
                + "\\end{axis} \n \\end{tikzpicture} \n \\end{document} \n "
                #print(s)
            elif(l==3):                
                s = "\\documentclass[]{article} \n\\usepackage{pgfplots} \n\\usepgfplotslibrary{external} \n\\tikzexternalize \n" \
                + "\\usepackage{tikz} \n \\usepackage{amsmath} \n  \\usepackage{pgfplots} \n \\usetikzlibrary{calc} \n" \
                + "\\pgfplotsset{compat = newest,every x tick label/.append style={font=\scriptsize},every y tick label/.append style={font=\\scriptsize}, every axis plot post/.style={line join=round}}" \
                + "\\begin{document} \n \\begin{tikzpicture} \n \\begin{axis}[ \n ylabel near ticks,\n xlabel near ticks, \n" \
                + "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5,}, \n"\
                + "xtick={520,525,530,535,540}, \n"\
                + "ytick={1.3,1.32,1.34,1.36,1.36898,1.38,1.4,1.42,1.44}, \n"\
                + "scaled y ticks=false, \n clip mode=individual,\n  "\
                + "xmin=520, \n xmax=540, \n ymin = 1.3, \n  ymax = 1.44,\n"\
                + "xlabel=$x$ ($m$), \n ylabel=$h$ ($m$)]"\
                + "\\addplot [blue] table {dx4diff"+ str(beta) +"h.dat}; \n"\
                + "\\addplot [green!80!black] table {dx5diff"+ str(beta) +"h.dat}; \n"\
                + "\\addplot [red] table {dx6diff"+ str(beta) +"h.dat}; \n"\
                + "\\addplot [cyan!70!white] table {dx7diff"+ str(beta) +"h.dat}; \n"\
                + "\\addplot [violet!70!white] table {dx8diff"+ str(beta) +"h.dat}; \n"\
                + "\\addplot [yellow!70!black] table {dx9diff"+ str(beta) +"h.dat};\n"\
                + "\\addplot [black] table {dx10diff"+ str(beta) +"h.dat};\n"\
                + "\\addplot [black, dashed] coordinates {(200.0,1.736397786) (800.0,1.736397786)};\n"\
                + "\\end{axis} \n \\end{tikzpicture} \n \\end{document} \n "
            else:                
                s = "\\documentclass[]{article} \n\\usepackage{pgfplots} \n\\usepgfplotslibrary{external} \n\\tikzexternalize \n" \
                + "\\usepackage{tikz} \n \\usepackage{amsmath} \n  \\usepackage{pgfplots} \n \\usetikzlibrary{calc} \n" \
                + "\\pgfplotsset{compat = newest,every x tick label/.append style={font=\scriptsize},every y tick label/.append style={font=\\scriptsize}, every axis plot post/.style={line join=round}}" \
                + "\\begin{document} \n \\begin{tikzpicture} \n \\begin{axis}[ \n ylabel near ticks,\n xlabel near ticks, \n" \
                + "yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=5,}, \n"\
                + "xtick={528,529,530,531,532,533,534,535,536}, \n"\
                + "ytick={1.3,1.32,1.34,1.36,1.36898,1.38,1.4,1.42,1.44}, \n"\
                + "scaled y ticks=false, \n clip mode=individual,\n  "\
                + "xmin=528, \n xmax=536, \n ymin = 1.3, \n  ymax = 1.44,\n"\
                + "xlabel=$x$ ($m$), \n ylabel=$h$ ($m$)]"\
                + "\\addplot [blue] table {dx4diff"+ str(beta) +"h.dat}; \n"\
                + "\\addplot [green!80!black] table {dx5diff"+ str(beta) +"h.dat}; \n"\
                + "\\addplot [red] table {dx6diff"+ str(beta) +"h.dat}; \n"\
                + "\\addplot [cyan!70!white] table {dx7diff"+ str(beta) +"h.dat}; \n"\
                + "\\addplot [violet!70!white] table {dx8diff"+ str(beta) +"h.dat}; \n"\
                + "\\addplot [yellow!70!black] table {dx9diff"+ str(beta) +"h.dat};\n"\
                + "\\addplot [black] table {dx10diff"+ str(beta) +"h.dat};\n"\
                + "\\addplot [black, dashed] coordinates {(200.0,1.736397786) (800.0,1.736397786)};\n"\
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
