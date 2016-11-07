import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from matplotlib2tikz import save as tikz_save

wdirb = "../../../../../data/postprocessing/HamilPlotLONG/o3/" 
sdirb = "../../../../../data/postprocessing/HPL2TexLONGzoom/o3/" 

gap = 1

sdir = sdirb
wdir = wdirb

if not os.path.exists(sdir):
    os.makedirs(sdir)   



x = []
HFT = []
HST = []
HTT = []
HiFT = []
HiST = []
HiTT = [] 

         
s = wdir + "allout.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    j = -1
    for row in readfile:
        if (j >= 0):
            dx = float(row[0])
            dt =float(row[1])
            t =float(row[2])
            beta = float(row[3])
            x.append(float(row[4]))
            HFT.append(0.5*float(row[5]))
            HST.append(0.5*float(row[6]))
            HTT.append(0.5*float(row[7]))
            HiFT.append(0.5*float(row[8]))
            HiST.append(0.5*float(row[9]))
            HiTT.append(0.5*float(row[10]))
        j = j + 1

xbeg = int((550)/dx)
xend = int((650)/dx)
xt =x[xbeg:xend:gap]
HFTt = HFT[xbeg:xend:gap]
HSTt = HST[xbeg:xend:gap]
HTTt = HTT[xbeg:xend:gap]
HiFTt = HiFT[xbeg:xend:gap]
HiSTt = HiST[xbeg:xend:gap]
HiTTt = HiTT[xbeg:xend:gap]



n = len(xt)
s = sdir + "HFT.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(xt[i]," ",HFTt[i])
        file1.write(s)
        
s = sdir + "HST.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(xt[i]," ",HSTt[i])
        file1.write(s)
        
s = sdir + "HTT.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(xt[i]," ",HTTt[i])
        file1.write(s)
        
s = sdir + "HiFT.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(xt[i]," ",HiFTt[i])
        file1.write(s)
        
s = sdir + "HiST.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(xt[i]," ",HiSTt[i])
        file1.write(s)
        
s = sdir + "HiTT.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(xt[i]," ",HiTTt[i])
        file1.write(s)

  