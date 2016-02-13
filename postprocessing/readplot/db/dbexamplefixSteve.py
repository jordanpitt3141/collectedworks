import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from matplotlib2tikz import save as tikz_save

order = "2af"
num = "17"
wdir = "../../../../data/raw/SteveDBEhsmaller/"
sdir = "../../../../data/postprocessing/SteveDBEhsmaller/"

if not os.path.exists(sdir):
    os.makedirs(sdir)

#wdir = "../../../data/Joe/dbh/o"+order+"/"
#sdir = "../../../../written/exportpic/dambreak/ex/"

#gap = 50
gap = 1

h0 = 0.1
g = 9.81

#plot results
"""         
s = wdir + "out1.txt"
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
            Eval = float(row[3])
            x.append(float(row[4]))
            h.append(float(row[5]))
            u.append(float(row[7]))
                
        j = j + 1
    x = array(x)
    u = array(u)
    h = array(h)

xbeg = int((520 - 300)/dx)
xend = int((545- 300)/dx)
xt =x[xbeg:xend:gap]
ht =h[xbeg:xend:gap]
"""
filesi = os.listdir(wdir)

times = []
Evals = []

for filen in filesi:
    s = wdir + filen
    with open(s,'r') as file1:
         readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
         for i,row in enumerate(readfile):
             if i == 1:
                t =float(row[2])
                Eval = float(row[3])
                break
    
    times.append(t)
    Evals.append(Eval)
E0 = 4159.44 
mix = zip(times,Evals)
mix = sorted(mix)
timess = [p[0] for p in mix]
Evalss = [p[1] for p in mix]

Evalabsdiff = []
for i in range(len(Evalss)):
    Evalabsdiff.append(4159.44 - Evalss[i])
    
Evalreldiff = []
for i in range(len(Evalss)):
    Evalreldiff.append((4159.44 - Evalss[i])/Evalss[0])


n = len(Evalss)
s = sdir + "Evals.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(timess[i]," ",Evalss[i])
        file1.write(s)
        
n = len(Evalabsdiff)
s = sdir + "Evalabsdiff.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(timess[i]," ",Evalabsdiff[i])
        file1.write(s)

n = len(Evalreldiff)
s = sdir + "Evalreldiff.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(timess[i]," ",Evalreldiff[i])
        file1.write(s)

"""
plot(xt,ht ,'-b')
xlim([350,650])
ylim([1.0,1.8])
xlabel("$x$ ($m$)")
ylabel("$h$ ($m$)")

#s = sdir + "o" + order +"n" + num +".tikz" 
#tikz_save(s);      
#clf();

"""
"""
n = len(xt)
s = sdir + "hz.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(xt[i]," ",ht[i])
        file1.write(s)
""" 