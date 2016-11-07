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

def experiment1(x,b,h0,h1,dx):
    n = len(x)
    u = zeros(n)
    h = ones(n)*h1
    for i in range(n):
        if (x[i] <0 and x[i] > -2*b):
            h[i] = h0

    return h,u
    
def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 



##### SEGUR AND HAMMACK EXPERIMENT 1
tl = 60.0
b = 0.61
h0 = 0.09
h1 = 0.1
g = 9.81

dx = 0.01
Cr = 0.5
l = Cr / sqrt(g*h1)
dt = l*dx
startx = -tl
endx = tl + dx
startt = 0
endt = 50 + dt  

l = "hs"

sdir = "../../../data/postprocessing/HS/"

if not os.path.exists(sdir):
    os.makedirs(sdir)


#number of boundary conditions (one side)
nfcBC = 4 #for flux calculation
nGsBC = 2 #for solving G from u,h
niBC = nGsBC + nfcBC #total
        
    
gap = 1
        
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
    
  
h,u = experiment1(x,b,h0,h1,dx)

cxlim = [-2.0,2.0]
cylim = [0.09,0.11]

ixbeg = int((cxlim[0]+ 60)/dx) 
ixend = int((cxlim[1]+ 60)/dx) + 1
x = x[ixbeg:ixend]
h = h[ixbeg:ixend]
u = u[ixbeg:ixend]
        
plot(x,h)
ylim([0.09,0.11])
xlim([-2.0,2.0])
xlabel("$x$ ($m$)")
ylabel("$h$ ($m$)")

m = len(x)
s = sdir + "h.dat"
with open(s,'a') as file2:
    for i in range(m):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",h[i])
        file2.write(s)


#stikz = sdir +str(l)+ ".tikz" 
#tikz_save(stikz);

"""
#make the tex file to create the document
s = "\\documentclass[]{article}\n\\usepackage{tikz} \n\\usepackage{amsmath} \n" \
+ "\\usepackage{pgfplots}\n\\usepgfplotslibrary{external}\n\\tikzexternalize\n" \
+ "\\pgfplotsset{compat=newest,every x tick label/.append style={font=\scriptsize}," \
+ "every y tick label/.append style={font=\scriptsize}, every axis plot post/.style={line join=round}}\n"\
+ "\\begin{document}\n   \\input{"+str(l)+ ".tikz" +"}\n\\end{document}"


#print(s)
filen = sdir +str(l)+ ".tex" 

file1 = open(filen, 'w')
file1.write(s)
file1.close() 

#make makefile
s = "LAT = pdflatex \nLATFLAGS = -shell-escape\n\n"
newl = str(l) +":\n\t $(LAT) $(LATFLAGS) " +str(l)+".tex\n\n"
s = s + newl
newl = "clean:\n\t rm -f  *~ ./*.log ./*.aux ./*.auxlock ./*.dep ./*.dpth ./*.pdf ./*.gz\n\n" 
s = s + newl

newl = "all: "
newl = newl + str(l) + " "
s = s + newl
   
filen =sdir + "Makefile" 

file1 = open(filen, 'w')
file1.write(s)
file1.close()
call(['make','-C',sdir,'all'])
"""       
