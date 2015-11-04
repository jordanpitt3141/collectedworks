import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from matplotlib2tikz import save as tikz_save





wdir = "../../../data/Joe2/smoothbd/" +wdirord+"/"+str(k) + "/"
sdir = "../../results/Tex/smoothdb/"

#gap = 50
gap = 50

h0 = 0.1
g = 9.81

         
s = wdir + num+".txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    h = []
    u = []
    x = []
    j = -1
    for row in readfile:       
        if (j >= 0):
            dx =float(row[0])
            dt =float(row[1])
            t =float(row[2])
            x.append(float(row[3]))
            h.append(float(row[4]))
            u.append(float(row[6]))
                
        j = j + 1
    x = array(x)
    u = array(u)
    h = array(h)

xt =x[::gap]
ht =h[::gap]

plot(xt,ht ,'-b')
xlim([0,1000])
ylim([0.8,2.0])
xlabel("x(m)")
ylabel("h(m)")
s = sdir + "o" + order +".tikz" 
tikz_save(s);      
clf();
    