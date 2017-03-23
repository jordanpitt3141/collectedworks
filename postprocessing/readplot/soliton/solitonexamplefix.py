import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from matplotlib2tikz import save as tikz_save

from Hamil import *
import os

def copyarraytoC(a):
    n = len(a)
    b = mallocPy(n)
    for i in range(n):
        writetomem(b,i,a[i])
    return b
    
def copyarrayfromC(a,n):
    b = [0]*n
    for i in range(n):
        b[i] = readfrommem(a,i)
        
    return b

order = "grim"
dxn = "12"
#wdir = "../../../../data/raw/solconlong/o"+order+"/" +dxn+ "/"
#sdir = "../../../../data/postprocessing/solconlong/o"+order+"/" +dxn+ "/"

#wdir = "../../../../data/raw/solconlong/o"+order+"/"
#sdir = "../../../../data/postprocessing/solconlongex/o"+order+"/"

#wdir = "../../../../data/raw/solconnonsmallg1n/o"+order+"/" +dxn+ "/"
#sdir = "../../../../data/postprocessing/solconnonsmallg1n/o"+order+"/" +dxn+ "/"

#wdir = "../../../../data/raw/solconnonsmallg10FDall/"+order+"/" +dxn+ "/"
#sdir = "../../../../data/postprocessing/solconnonsmallg10FDall/ex/"+order+"/" +dxn+ "/"

wdir = "../../../../data/raw/FDreredo/"+order+"/" +dxn+ "/"
sdir = "../../../../data/postprocessing/FDreredo/exsol/"+order+"/" +dxn+ "/"

if not os.path.exists(sdir):
    os.makedirs(sdir)

gap = 4
gaps = 2

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
      
startx = -50
endx = 250      
#s = wdir + "saveoutputtslast.txt"
s = wdir + "outlast.txt"
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
            h.append(float(row[7]))
            u.append(float(row[8]))
               
        j = j + 1
    x = arange(startx,endx+dx,dx)

#hig10 -250 to 250
#lowg10 -200 to 700 ?
#hig1 -100 to 100 ?

a1 = 0.7
a0 = 1.0
g = 9.81
#g= 1.0
t0 = 0.0
n = len(x) 
ht,ut = solitoninit(n,a0,a1,g,x,t,dx)   
hti,uti = solitoninit(n,a0,a1,g,x,t0,dx) 

niBC = 3    
xbeg = arange(startx - niBC*dx,startx,dx)
xend = arange(endx + dx,endx + (niBC+1)*dx) 
hbeg = h[0]*ones(niBC)
hend = h[-1]*ones(niBC)
ubeg = u[0]*ones(niBC)
uend = u[-1]*ones(niBC)   

htbeg = ht[0]*ones(niBC)
htend = ht[-1]*ones(niBC)
utbeg = ut[0]*ones(niBC)
utend = ut[-1]*ones(niBC)  

htibeg = hti[0]*ones(niBC)
htiend = hti[-1]*ones(niBC)
utibeg = uti[0]*ones(niBC)
utiend = uti[-1]*ones(niBC)         

xbc =  concatenate([xbeg,array(x),xend])
hbc =  concatenate([hbeg,array(h),hend])
ubc =  concatenate([ubeg,array(u),uend])

htbc =  concatenate([htbeg,array(ht),htend])
utbc =  concatenate([utbeg,array(ut),utend])

htibc =  concatenate([htibeg,array(hti),htiend])
utibc =  concatenate([utibeg,array(uti),utiend])

xbc_c = copyarraytoC(xbc)
hbc_c = copyarraytoC(hbc)
ubc_c = copyarraytoC(ubc)
htbc_c = copyarraytoC(htbc)
utbc_c = copyarraytoC(utbc)
htibc_c = copyarraytoC(htibc)
utibc_c = copyarraytoC(utibc)

Eval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)    
Evalt = HankEnergyall(xbc_c,htbc_c,utbc_c,g,n + 2*niBC,niBC,dx)
Evalti = HankEnergyall(xbc_c,htibc_c,utibc_c,g,n + 2*niBC,niBC,dx)    

    

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
        
s = sdir + "E.dat"
with open(s,'w') as file1:
    s ="%20s%5s%1.15f\n" %("Eval: "," ",Eval)
    file1.write(s)
    s ="%20s%5s%1.15f\n" %("Evalt: ", " ",Evalt)
    file1.write(s)
    s ="%20s%5s%1.15f\n" %("Evalti: "," ",Evalti)
    file1.write(s)
    s ="%20s%5s%1.15f\n" %("E - Eti: "," ",Eval - Evalti)
    file1.write(s)
    s ="%20s%5s%1.15f\n" %("E - Eti/Eti: "," ",(Eval - Evalti)/Evalti)
    file1.write(s)
    s ="%20s%5s%1.15f\n" %("E - Et: "," ",Eval - Evalt)
    file1.write(s)
    s ="%20s%5s%1.15f\n" %("E - Et/Et: "," ",(Eval - Evalt)/Evalt)
    file1.write(s)

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

    
    
    
   