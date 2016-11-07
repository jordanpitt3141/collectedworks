import csv
from numpy.linalg import norm
from scipy import *
from Hamil import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

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



wdir = "../../../../data/raw/DSWalpha/o3/9/4/"
#wdir = "../../../../data/raw/DSWalpha/o3/9/1/"
#wdir = "../../../../data/raw/DSWalpha/o3/9/3/"


#wdir = "../../../../data/raw/SteveDBEhsmaller/"

gap = 1
g = 9.81

dx = 0.0
dt = 0.0
t = 0.0
diffuse = 0.0

startx = 200
endx = 1000 
g = 9.81
niBC = 3
         
s = wdir + "outlast.txt"
 
#s = wdir + "saveoutputtslast.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    h = []
    u = []
    x = []    
    hs = []
    us = []
    xs = []
    dxs = []
    dts = []
    ts = []
    diffuses = []
    Evals = []
    j = -1
    for row in readfile:
        if (j == -1):
            print('initial')
        elif (row[0] == 'dx' and j >= 0):
            print(1,j,row[0])
            x = array(x)
            u = array(u)
            h = array(h)
            
            n = len(x)
            xbeg = arange(startx - niBC*dx,startx,dx)
            xend = arange(endx + dx,endx + (niBC+1)*dx) 
            hbeg = h[0]*ones(niBC)
            hend = h[-1]*ones(niBC)
            ubeg = u[0]*ones(niBC)
            uend = u[-1]*ones(niBC)            
            
            xbc =  concatenate([xbeg,x,xend])
            hbc =  concatenate([hbeg,h,hend])
            ubc =  concatenate([ubeg,u,uend])
            
            xbc_c = copyarraytoC(xbc)
            hbc_c = copyarraytoC(hbc)
            ubc_c = copyarraytoC(ubc)
            
            Eval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)
        
            Evals.append(Eval)
            
            hs.append(h)
            us.append(u)
            xs.append(x)
            dxs.append(dx)
            dts.append(dt)
            ts.append(t)
            diffuses.append(diffuse)
            h = []
            u = []
            x = []
            
        else:
            dx =float(row[0])
            dt =float(row[1])
            t =float(row[2])
            x.append(float(row[3]))
            h.append(float(row[4]))
            u.append(float(row[6]))
            diffuse = float(row[7])
        j = j + 1
   