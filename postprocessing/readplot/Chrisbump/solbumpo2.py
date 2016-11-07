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

diffw = "11"
wdirord = "o2"

#dx = 0.1
#filename = "out1.txt"
#filename = "out181651.txt"
#filename = "out363302.txt"
#filename = "out544953.txt"
#filename = "out726604.txt"
#filename = "out908255.txt"
filename = "outlast.txt"

    
cylim = [0.0,2.0]
cxlim = [-150,250]
gap = 1

wdir = "../../../../data/raw/ChrisDatn/solbumpEnerg/17/"
sdir = "../../../../data/postprocessing/ChrisDatn/solbumpEnerg/17/50/"
if not os.path.exists(sdir):
        os.makedirs(sdir)
     
s = wdir + filename
with open(s,'r') as file1:
     readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     h = []
     u = []
     x = []
     bed = []
     j = -1
     for row in readfile:       
         if (j >= 0):
            dx = float(row[0])
            dt = float(row[1])
            t = float(row[2])
            E = float(row[3])
            x.append(float(row[4]))
            h.append(float(row[5]))
            u.append(float(row[7]))
            bed.append(float(row[8]))
         j = j + 1
     igap = int(gap)
     xbeg = int((cxlim[0] + 150)/dx)
     xend = int((cxlim[1] + 150)/dx)
     xt = x[xbeg:xend:igap]
     ht = h[xbeg:xend:igap]
     bedt = bed[xbeg:xend:igap]

     x = array(xt)
     h = array(ht)
     bed = array(bedt)     


s = str(dx) 
plot(x,h + bed,'b' ,label=s)
plot(x, bed,'r' ,label=s)
ylim(cylim)
xlim(cxlim)
xlabel("$x$ ($m$)")
ylabel("$h$ ($m$)")
#legend()

    
stikz = sdir +"bump.tikz" 
tikz_save(stikz);
s = "Bump " 
title(s)
s = sdir +"b.png"       
savefig(s, bbox_inches='tight')
legend()
s = sdir +"bleg.png"       
savefig(s, bbox_inches='tight')        
clf()

n = len(x)
s = sdir + "b.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",bed[i])
        file1.write(s)
s = sdir + "w.dat"
with open(s,'w') as file2:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",h[i] + bed[i])
        file2.write(s)


#make the tex file to create the document
s = "\\documentclass[]{article}\n\\usepackage{tikz} \n\\usepackage{amsmath} \n" \
+ "\\usepackage{pgfplots}\n\\usepgfplotslibrary{external}\n\\tikzexternalize\n" \
+ "\\pgfplotsset{compat=newest,every x tick label/.append style={font=\scriptsize}," \
+ "every y tick label/.append style={font=\scriptsize}, every axis plot post/.style={line join=round}}\n"\
+ "\\begin{document}\n   \\input{bump.tikz" +"}\n\\end{document}"


filen = sdir + "bump.tex" 

file1 = open(filen, 'w')
file1.write(s)
file1.close() 

#make makefile
s = "LAT = pdflatex \nLATFLAGS = -shell-escape\n\n"
newl = "bump:\n\t $(LAT) $(LATFLAGS) bump.tex\n\n"
s = s + newl
newl = "clean:\n\t rm -f  *~ ./*.log ./*.aux ./*.auxlock ./*.dep ./*.dpth ./*.pdf ./*.gz\n\n" 
s = s + newl

newl = "all: "

   
filen =sdir + "Makefile" 

file1 = open(filen, 'w')
file1.write(s)
file1.close()           
