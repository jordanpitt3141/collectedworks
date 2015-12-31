import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from matplotlib2tikz import save as tikz_save

order = "3"
dxn = "11"
#wdir = "../../../../data/raw/solconlong/o"+order+"/" +dxn+ "/"
#sdir = "../../../../data/postprocessing/solconlong/o"+order+"/" +dxn+ "/"

#wdir = "../../../../data/raw/solconlong/o"+order+"/"
#sdir = "../../../../data/postprocessing/solconlongex/o"+order+"/"

#wdir = "../../../../data/raw/solconnonsmallg1n/o"+order+"/" +dxn+ "/"
#sdir = "../../../../data/postprocessing/solconnonsmallg1n/o"+order+"/" +dxn+ "/"

wdir = "../../../../data/raw/hinonling1/o"+order+"/" +dxn+ "/"
sdir = "../../../../data/postprocessing/hinonling1ex/o"+order+"/" +dxn+ "/"


if not os.path.exists(sdir):
    os.makedirs(sdir)

gap = 2
gaps = 1

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
    x = arange(-100,100+dx,dx)
#hig10 -250 to 250
#lowg10 -200 to 700 ?
#hig1 -100 to 100 ?

    
a1 = 1.0
a0 = 1.0
#g = 9.81
g= 1.0
t0 = 0.0
n = len(x) 
ht,ut = solitoninit(n,a0,a1,g,x,t,dx)   
hti,uti = solitoninit(n,a0,a1,g,x,t0,dx) 

xbeg = int((-50 - x[0])/dx)
xend = int((250 - x[0])/dx)
xe = x[xbeg:xend:gap]
ue = u[xbeg:xend:gap]
he = h[xbeg:xend:gap]
ute = ut[xbeg:xend:gaps]
hte = ht[xbeg:xend:gaps]
utie = uti[xbeg:xend:gaps]
htie = hti[xbeg:xend:gaps]
xte = x[xbeg:xend:gaps]

ue = array(ue)
he = array(he)


n = len(xe)
s = sdir + "h.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(xe[i]," ",he[i])
        file1.write(s)
s = sdir + "u.dat"
with open(s,'w') as file2:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(xe[i]," ",ue[i])
        file2.write(s)
        
n = len(xte)       
s = sdir + "ht.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(xte[i]," ",hte[i])
        file1.write(s)
s = sdir + "ut.dat"
with open(s,'w') as file2:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(xte[i]," ",ute[i])
        file2.write(s)
        
s = sdir + "hti.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(xte[i]," ",htie[i])
        file1.write(s)
s = sdir + "uti.dat"
with open(s,'w') as file2:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(xte[i]," ",utie[i])
        file2.write(s)


#Didnt like that
"""
plot(xte,htie,'-k') 
plot(xte,hte,'-b') 
plot(xe,he,'.r')
ylim([0.8,2.0])
xlim([-50.0,250.0])
xlabel("$x$ ($m$)")
ylabel("$h$ ($m$)") 
stikz = sdir + order + ".tikz" 
tikz_save(stikz);
clf()

plot(xte,htie,'-k') 
plot(xte,hte,'-b') 
plot(xe,he,'.r')
ylim([0.8,2.0])
xlim([200.0,250.0])
xlabel("$x$ ($m$)")
ylabel("$h$ ($m$)") 
stikz = sdir + order + "z.tikz" 
tikz_save(stikz);
clf()
"""

    
    
    
   