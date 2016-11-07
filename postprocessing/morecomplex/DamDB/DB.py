import csv
from numpy.linalg import norm
from scipy import *
from Hamil import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

def HamilDB(alpha):
    return 10398.6 - 0.7848*((2.0/alpha)*tanh(500.0*alpha))

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



wdir = "../../../../data/postprocessing/DBexact/"

if not os.path.exists(wdir):
    os.makedirs(wdir)

Numevals= []
Anaevals =[]


diffuses = [0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1.0,2.5,5.0,7.5,10.0,25.0,50.0,75.0,100.0,250.0,500.0,750.0,1000.0]
    
for diffuse in diffuses:
    dx = (10.0/(2**10))
    h0 = 1.0
    h1 = 1.8
    x0 = 500
    g = 9.81
    Cr = 0.5
    l = Cr / (sqrt(g*(1)))
    dt = l*dx
    startx = 0.0
    endx = 1000.0
    startt = 0.0
    endt = 2*dt
    
    niBC = 3
             
    s = wdir + "outlast.txt"
     
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
        
    eta0 = (h1-h0)
    h,u = dambreaksmooth(x,x0,h0,eta0,diffuse,dx)
    
    
    niBC = 4
    u0 = u[0]*ones(niBC)
    u1 = u[-1]*ones(niBC)   
    h0 = h[0]*ones(niBC)
    h1 = h[-1]*ones(niBC)
    
    
    xbeg = arange(startx - niBC*dx,startx,dx)
    xend = arange(endx + dx,endx + (niBC+1)*dx) 
    
    xbc =  concatenate([xbeg,x,xend])
    hbc =  concatenate([h0,h,h1])
    ubc =  concatenate([u0,u,u1])
    
    xbc_c = copyarraytoC(xbc)
    hbc_c = copyarraytoC(hbc)
    ubc_c = copyarraytoC(ubc)
   
           
    Eval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)
    Anaeval = HamilDB(diffuse)
    
    Numevals.append(Eval)
    Anaevals.append(Anaeval)

   