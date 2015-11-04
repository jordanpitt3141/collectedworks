import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

order = "3"
num = "6"
wdir = "../../../data/Joe/soliton/o"+order+"/"+num+"/"
sdir = "../../../../written/test/o"+order+"/"

gap = 1

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
         
s = wdir + "saveoutputtslast.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    h = []
    u = []
    ht = []
    ut = []
    j = -1
    for row in readfile:       
        if (j >= 0):
            dx =float(row[0])
            dt =float(row[1])
            t =float(row[2])
            h.append(float(row[3]))
            u.append(float(row[5]))
                
        j = j + 1
    x = arange(-500,1500+dx,dx)
    u = array(u)
    h = array(h)

a1 = 1.0
a0 = 10.0
g = 9.81
t0 = 0.0
n = len(x) 
ht,ut = solitoninit(n,a0,a1,g,x,t,dx)   
hti,uti = solitoninit(n,a0,a1,g,x,t0,dx)  
    
    
    
n = len(h)

s = sdir + "h.dat"
with open(s,'a') as file1:
    for i in range(0,n,gap):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",h[i])
        file1.write(s)
s = sdir + "u.dat"
with open(s,'a') as file2:
    for i in range(0,n,gap):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",u[i])
        file2.write(s)
s = sdir + "ht.dat"
with open(s,'a') as file1:
    for i in range(0,n,gap):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",ht[i])
        file1.write(s)
s = sdir + "ut.dat"
with open(s,'a') as file2:
    for i in range(0,n,gap):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",ut[i])
        file2.write(s)
s = sdir + "hti.dat"
with open(s,'a') as file1:
    for i in range(0,n,gap):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",hti[i])
        file1.write(s)
s = sdir + "uti.dat"
with open(s,'a') as file2:
    for i in range(0,n,gap):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",uti[i])
        file2.write(s)       