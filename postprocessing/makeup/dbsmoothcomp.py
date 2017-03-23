# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 23:36:31 2015

@author: jordan
"""
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import os

def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 

def dambreaksmooth(x,x0,base,eta0,diffuse,dx):
    from numpy import tanh
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    
    for i in range(n):
        h[i] = base + 0.5*eta0*(1 + tanh(diffuse*(x0 - abs(x[i]))))

    return h,u
   
dx = 10.0 / (2**8)
l = 0.01
dt = l*dx
theta = 1.2
startx = 400
endx = 600.0 + dx
startt = 0.0
endt = dt  
    
g = 9.81
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
base = 1.0
eta0 = 0.8
x0 = 500

sdir = "../../../data/postprocessing/dbsmoothFinAGN/"

if not os.path.exists(sdir):
        os.makedirs(sdir)

diffuses = [0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1.0,2.5,5.0,7.5,10.0,25.0,50.0,75.0,100.0,250.0,500.0,750.0,1000.0]
listalp = [1,6,9,12]
for i in listalp:
    diffuse = diffuses[i]
    hh,uu = dambreaksmooth(x,x0,base,eta0,diffuse,dx)
    s = str(diffuse)
    plot(x,hh,label=s)  
    
    s = sdir + str(diffuse) + "h.dat"
    with open(s,'w') as file1:
        for i in range(n):
            s ="%3.8f%5s%1.15f\n" %(x[i]," ",hh[i])
            file1.write(s)

#s = "Initial Conditions Dam Break: dx = " + str(dx)
#title(s)
ylim([1.0,1.8])
xlim([300,700])
xlabel("$x$ ($m$)")
ylabel("$h$ ($m$)")
#legend()


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



"""
sdir = "../../results/show/steve2/initial/"

for ll in range(3,16):
    dx = 10.0 / (2**ll)
    l = 0.01
    dt = l*dx
    theta = 1.2
    startx = 0.0
    endx = 1000.0 + dx
    startt = 0.0
    endt = dt  
        
    g = 9.81
        
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
    base = 1.0
    eta0 = 0.8
    x0 = 500
    
    diffuses = [0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1.0,2.5,5.0,7.5,10.0,25.0,50.0,75.0,100.0,250.0,500.0,750.0,1000.0]
    
    for i in range(len(diffuses)):
        diffuse = diffuses[i]
        hh,uu = dambreaksmooth(x,x0,base,eta0,diffuse,dx)
        s = "diff = " + str(diffuse)
        plot(x,hh,label=s)   
    
    s = "Initial Conditions Dam Break: dx = " + str(dx)
    title(s)
    xlabel("x (m)")
    ylabel("h (m)")
    legend()
    s = sdir +"diffalldx"+ str(dx) +".png"       
    savefig(s, bbox_inches='tight')        
    clf()
"""    