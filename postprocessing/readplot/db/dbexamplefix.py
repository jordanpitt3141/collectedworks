import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from matplotlib2tikz import save as tikz_save

order = "2af"
num = "17"
wdir = "../../../../data/raw/ndbh/o"+order+"/"
sdir = "../../../../data/postprocessing/ndbhex/o"+order+"/"+num+"/"

if not os.path.exists(sdir):
    os.makedirs(sdir)

#wdir = "../../../data/Joe/dbh/o"+order+"/"
#sdir = "../../../../written/exportpic/dambreak/ex/"

#gap = 50
gap = 30

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

xbeg = int((350)/dx)
xend = int((650)/dx)
xt =x[xbeg:xend:gap]
ht =h[xbeg:xend:gap]

"""
plot(xt,ht ,'-b')
xlim([350,650])
ylim([1.0,1.8])
xlabel("$x$ ($m$)")
ylabel("$h$ ($m$)")
"""
#s = sdir + "o" + order +"n" + num +".tikz" 
#tikz_save(s);      
#clf();

n = len(xt)
s = sdir + "h.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(xt[i]," ",ht[i])
        file1.write(s)
  