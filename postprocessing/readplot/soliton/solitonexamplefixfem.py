import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

order = "3"
num = "6"
wdir = "../../../../data/femelemdusolcon/10.txt"
sdir = "../../../results/Chrisfem1/sol/"

if not os.path.exists(sdir):
    os.makedirs(sdir)

gapt = 10
gape = 20

h0 = 0.1
g = 9.81

def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a

def soliton (x,t,g,a0,a1):
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)
  
def solitoninit(n,a0,a1,g,x,t0,dx):
    h = zeros(n)
    u = zeros(n)
    c = sqrt(g*(a0 + a1))
    for i in range(n):
        h[i] = soliton(x[i],t0,g,a0,a1)
        u[i] =  c* ((h[i] - a0) / h[i])
    
    return h,u
         
s = wdir
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    h = []
    x = []
    u = []
    ht = []
    ut = []
    j = -1
    for row in readfile:       
        if (j >= 0):
            dx =float(row[0])
            dt =float(row[1])
            t =float(row[2])
            x.append(float(row[3]))
            h.append(float(row[4]))
            u.append(float(row[6]))
            ht.append(float(row[7]))
            ut.append(float(row[8]))
                
        j = j + 1
    u = array(u)
    h = array(h)
    ut = array(ut)
    ht = array(ht)

a1 = 1.0
a0 = 10.0
g = 9.81
t0 = 0.0
n = len(x)  
hti,uti = solitoninit(n,a0,a1,g,x,t0,dx)  
    
    
    
n = len(h)

s = sdir + "h.dat"
with open(s,'a') as file1:
    for i in range(0,n,gape):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",h[i])
        file1.write(s)
s = sdir + "u.dat"
with open(s,'a') as file2:
    for i in range(0,n,gape):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",u[i])
        file2.write(s)
s = sdir + "ht.dat"
with open(s,'a') as file1:
    for i in range(0,n,gapt):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",ht[i])
        file1.write(s)
s = sdir + "ut.dat"
with open(s,'a') as file2:
    for i in range(0,n,gapt):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",ut[i])
        file2.write(s)
s = sdir + "hti.dat"
with open(s,'a') as file1:
    for i in range(0,n,gapt):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",hti[i])
        file1.write(s)
s = sdir + "uti.dat"
with open(s,'a') as file2:
    for i in range(0,n,gapt):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",uti[i])
        file2.write(s)       