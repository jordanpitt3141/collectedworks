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
    
def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 
    
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



wdir = "../../../../data/postprocessing/Eanal/"

if not os.path.exists(wdir):
    os.makedirs(wdir)

Evals= []
Evalfs =[]

relerrs = []
relerrfs = []
dxs = []
    
for k in range(5,20):
    dx = 100.0 / (2**k)
    a0 = 1.0
    a1 = 1.0
    g = 9.81
    Cr = 0.5
    l = Cr / (sqrt(g*(a0 + a1)))
    dt = l*dx
    startx = -100.0
    endx = 500.0
    startt = 0.0
    endt = 100 + dt
    
    niBC = 3
             
    s = wdir + "outlast.txt"
     
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
        
    t0 = t[0]    
    h,u = solitoninit(n,a0,a1,g,x,t0,dx)
    
    tf = t[-1]   
    hf,uf = solitoninit(n,a0,a1,g,x,tf,dx)
    
    niBC = 4
    u0 = u[0]*ones(niBC)
    u1 = u[-1]*ones(niBC)   
    h0 = h[0]*ones(niBC)
    h1 = h[-1]*ones(niBC)
    
    
    
    uf0 = uf[0]*ones(niBC)
    uf1 = uf[-1]*ones(niBC)   
    hf0 = hf[0]*ones(niBC)
    hf1 = hf[-1]*ones(niBC)
    
    xbeg = arange(startx - niBC*dx,startx,dx)
    xend = arange(endx + dx,endx + (niBC+1)*dx) 
    
    xbc =  concatenate([xbeg,x,xend])
    hbc =  concatenate([h0,h,h1])
    ubc =  concatenate([u0,u,u1])
    
    xbc_c = copyarraytoC(xbc)
    hbc_c = copyarraytoC(hbc)
    ubc_c = copyarraytoC(ubc)
    
    hfbc =  concatenate([hf0,hf,hf1])
    ufbc =  concatenate([uf0,uf,uf1])
    
    hfbc_c = copyarraytoC(hfbc)
    ufbc_c = copyarraytoC(ufbc)
    
    
           
    Eval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)
    
    Evalf = HankEnergyall(xbc_c,hfbc_c,ufbc_c,g,n + 2*niBC,niBC,dx)
    Evals.append(Eval)
    Evalfs.append(Evalf)
    
    #rel error with Eval as exact
    relerr = abs(Eval - Evalf)/ abs(Eval)
    relerrs.append(relerr)
    
    #rel error with Evalf as exact
    relerrf = abs(Evalf - Eval)/ abs(Evalf)
    relerrfs.append(relerrf)
    
    dxs.append(dx)

   