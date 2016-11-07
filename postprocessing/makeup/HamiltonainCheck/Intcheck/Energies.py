
import csv
from numpy.linalg import norm
from scipy import *
import os
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from Int import *
def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 
    
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

def HamilDB(alpha):
    return 10398.6 - 0.7848*((2.0/alpha)*tanh(500.0*alpha))
    
def Soliton():
    return 1527.68293

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
    
def dambreaksmooth(x,x0,base,eta0,diffuse,dx):
    from numpy import tanh
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    
    for i in range(n):
        h[i] = base + 0.5*eta0*(1 + tanh(diffuse*(x0 - abs(x[i]))))

    return h,u
    
def constants(x,h0,u0,dx):
    n = len(x)
    h = h0*ones(n)
    u = u0*ones(n)
    return h,u


def interpquarticval(aj,bj,cj,dj,ej,xj,x):
    
    return aj*(x -xj)*(x -xj)*(x -xj)*(x -xj) + bj*(x -xj)*(x -xj)*(x -xj) \
    + cj*(x -xj)*(x -xj) + dj*(x -xj)+ ej
    
def interpquarticgrad(aj,bj,cj,dj,ej,xj,x):
    
    return 4*aj*(x -xj)*(x -xj)*(x -xj) + 3*bj*(x -xj)*(x -xj) \
    + 2*cj*(x -xj) + dj
    
def interpquartcoeff(q,j,dx):
    i24 = 1.0 / 24.0
    i12 = 1.0 / 12.0
    idx = 1.0/dx
    aj = i24*idx*idx*idx*idx*(q[j+2] - 4*q[j+1] + 6*q[j] - 4*q[j-1] + q[j-2])
    bj = i12*idx*idx*idx*(q[j+2] - 2*q[j+1] + 2*q[j-1] - q[j-2])
    cj = i24*idx*idx*(-q[j+2] + 16*q[j+1] - 30*q[j] + 16*q[j-1] - q[j-2])
    dj = i12*idx*(-q[j+2] + 8*q[j+1] - 8*q[j-1] + q[j-2])
    ej = q[j]
    
    return aj,bj,cj,dj,ej
    
def interpquartic(u,h,x,xh,nBC):
    n = len(x)
    nu = []
    nh = []
    nx = []
    
    for i in range(nBC,n-nBC):
        aj,bj,cj,dj,ej = interpquartcoeff(h,i,dx)
        nh.append(interpquarticval(aj,bj,cj,dj,ej,x[i],xh[i]))
        aj,bj,cj,dj,ej = interpquartcoeff(u,i,dx)
        nu.append(interpquarticval(aj,bj,cj,dj,ej,x[i],xh[i]))
        nx.append(xh[i])
    return nh, nu,nx

#### Dam  Break
sdir = "../../../../../data/postprocessing/IntCheck//"
if not os.path.exists(sdir):
    os.makedirs(sdir)
    

deltaxa = range(1,22)
#diffuses = [0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1.0,2.5,5.0,7.5,10.0,25.0,50.0,75.0,100.0,250.0,500.0,750.0,1000.0]
Evals = []


s = sdir + "Energy.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     writefile2.writerow(['dx' ,'alpha','Numerical','analytical', "Rel Err"])
 
for ll in deltaxa:
        
    dx = (100.0 / (2**ll))
    l = 0.01
    dt = l*dx
    startx = 0.0
    endx = 1000.0 + dx
    startt = 0.0
    endt = 30.0+(dt*0.9)  
            
    szoomx = startx
    ezoomx = endx
            
    #number of boundary conditions (one side)
    niBC = 2 #total
            
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
    hm,um = constants(x,2,2,dx)
                    
    umbegi = zeros(niBC)
    umendi = zeros(niBC)
    hmbegi = ones(niBC)
    hmendi = ones(niBC)    
            
    for i in range(niBC):
        umbegi[i] = um[0]
        umendi[i] = um[-1]
        hmbegi[i] = hm[0]
        hmendi[i] = hm[-1]
            
    umbeg = umbegi
    umend = umendi
    hmbeg = hmbegi
    hmend = hmendi
            
    xbeg = arange(startx - niBC*dx,startx,dx)
    xend = arange(endx + dx,endx + (niBC+1)*dx,dx) 

    xbc =  concatenate([xbeg,x,xend])          
    ubc = concatenate([umbeg,um,umend]) 
    hbc = concatenate([hmbeg,hm,hmend])

    xbc_c = copyarraytoC(xbc)
    ubc_c = copyarraytoC(ubc)
    hbc_c = copyarraytoC(hbc) 
    
    
    Evali = Intall(xbc_c,hbc_c,n + 2*niBC,niBC,dx) 
    
    H0 = 2*(1000 + dx)
    
    relErr = abs(H0 - Evali) / abs(H0)
    
    s = sdir + "Energy.txt"
    with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)      
                       
         writefile2.writerow([str(dx),str(Evali),str(H0), str(relErr)])  
         
    s = sdir + "relE.dat"     
    with open(s,'a') as file3:
            s ="%3.8f%5s%1.50f\n" %(dx," ",relErr)
            file3.write(s)
