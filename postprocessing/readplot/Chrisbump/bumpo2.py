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
#filename = "out10251.txt"
#filename = "out20553.txt"
#filename = "out30804.txt"
#filename = "out41055.txt"
#filename = "out51306.txt"
#filename = "out61608.txt"
filename = "outlast.txt"

    
cylim = [0.0,2.0]
cxlim = [500,1500]
gap = 1

wdir = "../../../../data/raw/bumpChrislonger/dx0p1/o2/"
sdir = "../../../../data/postprocessing/bumpChrislongerfix/dx0p1/o2/t700/"
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
            x.append(float(row[3]))
            h.append(float(row[4]))
            u.append(float(row[6]))
            bed.append(float(row[7]))
         j = j + 1
     igap = int(gap)
     xbeg = int((cxlim[0] + 1000)/dx)
     xend = int((cxlim[1] + 1000)/dx)
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
