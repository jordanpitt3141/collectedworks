from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv
from Serre3fem import *
from numpy.linalg import norm

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
    
def TDMApy(a,b,c,d):
    n = len(d)
    alpha = []
    beta = []
    x = [0]*n
    
    alpha.append((1.0*c[0])/b[0])
    beta.append((1.0*d[0])/b[0] )  
 
    for i in range(1,n-1):
        m = 1.0 / (b[i] - a[i-1]*alpha[i-1])
        alpha.append(c[i]* m)
        beta.append((d[i] - a[i-1]*beta[i-1]) * m)
        
    m = 1.0 / (b[n-1] - a[n-2]*alpha[n-2])
    beta.append((d[n-1] - a[n-2]*beta[n-2]) * m)  

    x[n-1] = beta[n-1]
    
    for i in range(n-2,-1,-1):
        x[i] = beta[i] - alpha[i]*x[i+1]
 
    return array(x)
    
#Tested from CFD forum post
def pentadiagsolve(e,a,d,c,f,B):
    n = len(d)
    X = zeros(n)
    
    for i in range(1,n-1):
        xmult = float(a[i-1]) / d[i-1]
        
        d[i] = d[i] - xmult*c[i-1]
        c[i] = c[i] - xmult*f[i-1]
        B[i] = B[i] - xmult*B[i-1]
        
        xmult = float(e[i-1]) /d[i-1]
        a[i] = a[i] - xmult*c[i-1]
        d[i+1] = d[i+1] - xmult*f[i-1]
        B[i+1] = B[i+1] - xmult*B[i-1]
        
    xmult = float(a[n-2]) / d[n-2]
    d[n-1] = d[n-1] - xmult*c[n-2]
    X[n-1] = (B[n-1] - xmult*B[n-2]) / float(d[n-1])
    X[n-2] = (B[n-2] - c[n-2]*X[n-1]) / float(d[n-2])
    
    for i in range(n-3,-1,-1):
        X[i] = (B[i] - f[i]*X[i+2] - c[i]*X[i+1])/float(d[i])
        
    return X    
    
def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 

def midpointtocellaverages(mq,dx):
    #no BC required, assumes that the averages and midpoints at the boundaries are the same
    idx = 1.0/dx
    i24 = 1.0 / 24.0
    n = len(mq)
    
    a = zeros(n-1)
    b = zeros(n)
    c = zeros(n-1)
    for i in range(1,n-1):
        ai = -i24
        bi = 26*i24
        ci = -i24

        a[i-1] = ai
        b[i] = bi
        c[i] = ci
    
    #i = 0
    i = 0
    ai =0.0 #-i24
    bi =1.0 #26*i24
    ci =0.0 #-i24

    b[i] = bi
    c[i] = ci
    
    #mq[i] = mq[i] - ai*qbeg[0]
    
    #i = 0
    i = n-1
    ai =0.0# -i24
    bi =1.0# 26*i24
    ci =0.0# -i24

    a[i-1] = ai
    b[i] = bi
    
    #mq[i] = mq[i] - ci*qend[0]
    
    q = TDMApy(a,b,c,mq)
    
    return q
    
def cellaveragestomidpoints(q,dx):
    #no BC required, assumes that the averages and midpoints at the boundaries are the same
    i24 = 1.0 / 24.0
    n = len(q)
    mq = zeros(n)
    for i in range(1,n-1):
        #iterate  over the cell midpoints, there are 2 edge values for each (except the first and last cell)
        
        #variables
        #ai = (q[i+1] - 2*q[i] + q[i-1])*0.5*idx*idx
        #bi = (q[i+1] - q[i-1])*0.5*idx
        ci = i24*(-q[i+1] + 26*q[i]  -q[i-1])
        mq[i] = ci
    
    #i = 0
    i = 0
    ci = q[i] #i24*(-q[i+1] + 26*q[i] - qbeg[0])
    mq[i] = ci
    
    #i = n-1
    i = n-1
    ci = q[i]#i24*(-qend[0] + 26*q[i] - q[i-1])
    mq[i] = ci 
    
    return mq

def solveufromGh(G,h,hbeg,hend,ubeg,uend,dx):
    #takes midpoint values of G,h and gives midpoint values of u
    idx = 1.0 / dx
    i12 = 1.0 / 12.0
    i3 = 1.0 / 3.0
    n = len(G)
    
    a = zeros(n-2)
    b = zeros(n-1)
    c = zeros(n)
    d = zeros(n-1)
    e = zeros(n-2)
    
    for i in range(2,n-2):
        th = h[i]
        thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*h[i-1] + h[i-2] )
        
        ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
        bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
        ci = th + (30*i12*idx*idx)*(i3*th*th*th)
        di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
        ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2
        
        
        
        a[i-2] = ai
        b[i-1] =  bi
        c[i] = ci
        d[i] = di
        e[i] = ei
        
    #boundary    
    #i=0
    i=0
    th = h[i]
    thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*hbeg[-1] + hbeg[-2] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2

 
    c[i] = ci
    d[i] = di
    e[i] = ei
    
    G[i] = G[i] - ubeg[-1]*bi - ubeg[-2]*ai
    
    #i=1
    i=1
    th = h[i]
    thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*h[i-1] + hbeg[-1] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2

 
    c[i] = ci
    d[i] = di
    e[i] = ei
    b[i-1] = bi 
    
    G[i] = G[i] - ubeg[-1]*ai
    
    #boundary    
    #i=n-2
    i=n-2
    th = h[i]
    thx = i12*idx*(-hend[0] + 8*h[i+1] - 8*h[i-1] + h[i-2] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2

    a[i-2] = ai
    b[i-1] = bi 
    c[i] = ci
    d[i] = di
    
    G[i] = G[i]- uend[0]*ei
    
    #i=n-1
    i=n-1
    th = h[i]
    thx = i12*idx*(-hend[1] + 8*hend[0] - 8*h[i-1] + h[i-2] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2

    a[i-2] = ai
    b[i-1] = bi 
    c[i] = ci
    
    G[i] = G[i] -uend[0]*di - uend[1]*ei
    
    u = pentadiagsolve(a,b,c,d,e,G)
    return u
    
def solveGfromuh(u,h,hbeg,hend,ubeg,uend,dx):
    #takes midpoint values of u,h and gives midpoint values of G    
    
    idx = 1.0 / dx
    i12 = 1.0 / 12.0
    i3 = 1.0 / 3.0
    n = len(u)
    
    G = zeros(n)
    
    for i in range(2,n-2):
        th = h[i]
        thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*h[i-1] + h[i-2] )
        
        ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
        bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
        ci = th + (30*i12*idx*idx)*(i3*th*th*th)
        di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
        ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2
        
        G[i] = ai*u[i-2] + bi*u[i-1] + ci*u[i] + di*u[i+1] + ei*u[i+2]
        
        
    #boundary    
    #i=0
    i=0
    th = h[i]
    thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*hbeg[-1] + hbeg[-2] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2
    
    G[i] = ai*ubeg[-2] + bi*ubeg[-1] + ci*u[i] + di*u[i+1] + ei*u[i+2]

    
    #i=1
    i=1
    th = h[i]
    thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*h[i-1] + hbeg[-1] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2

    G[i] = ai*ubeg[-1] + bi*u[i-1] + ci*u[i] + di*u[i+1] + ei*u[i+2]
    
    #boundary    
    #i=n-2
    i=n-2
    th = h[i]
    thx = i12*idx*(-hend[0] + 8*h[i+1] - 8*h[i-1] + h[i-2] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2

    G[i] = ai*u[i-2] + bi*u[i-1] + ci*u[i] + di*u[i+1] + ei*uend[0]
    
    #i=n-1
    i=n-1
    th = h[i]
    thx = i12*idx*(-hend[1] + 8*hend[0] - 8*h[i-1] + h[i-2] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2

    G[i] = ai*u[i-2] + bi*u[i-1] + ci*u[i] + di*uend[0] + ei*uend[1]

    return G
    

def dambreak(x,xc,hf,hl):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    
    for i in range(n):
       
        if(x[i] >= xc):
            h[i] = hl
        else:
            h[i]= hf
        
    return u,h
    
def dambreaksmooth(x,x0,base,eta0,diffuse,dx):
    from numpy import tanh
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    
    for i in range(n):
        h[i] = base + 0.5*eta0*(1 + tanh(diffuse*(x0 - abs(x[i]))))

    return u,h 
    

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

def experiment1(x,b,h0,h1,dx):
    n = len(x)
    u = zeros(n)
    h = ones(n)*h1
    for i in range(n):
        if (x[i] <0 and x[i] > -2*b):
            h[i] = h0

    return h,u
 
""" 
#Dam Break 
import os

wdir = "../../data/dbh/o3fem/"

if not os.path.exists(wdir):
    os.makedirs(wdir)
    
dx = 1.0
l = 0.01
dt = l*dx
startx = 0.0
endx = 1000.0 + dx
startt = 0.0
endt = 30 + dt
        
szoomx = startx
ezoomx = endx
        
#number of boundary conditions (one side)
nfcBC = 4 #for flux calculation
nGsBC = 2 #for solving G from u,h
niBC = nGsBC + nfcBC #total
        
g = 9.81
    
gap = int(0.5/dt)
        
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
    
xc = 500
hf = 1.8
hl = 1.0    
um,hm = dambreak(x,xc,hf,hl)
                
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
        
#calculate G midpoint
cnBC = niBC - nGsBC
        
umbc = concatenate([umbeg[-cnBC:],um,umend[0:cnBC]]) 
hmbc = concatenate([hmbeg[-cnBC:],hm,hmend[0:cnBC]])       
Gmbc = solveGfromuh(umbc,hmbc,hmbeg[0:-cnBC],hmend[-cnBC:],umbeg[0:-cnBC],umend[-cnBC:],dx)  
        
#calculate averages
Gabc = midpointtocellaverages(Gmbc,dx)
habc = midpointtocellaverages(hmbc,dx)
uabc = midpointtocellaverages(umbc,dx)
        
#so we can just go from here with Ga ang ha?
Gabeg = Gabc[0:cnBC]
Ga = Gabc[cnBC:-cnBC]
Gaend = Gabc[-cnBC:]
habeg = habc[0:cnBC]
ha = habc[cnBC:-cnBC]
haend = habc[-cnBC:]
uabeg = uabc[0:cnBC]
ua = uabc[cnBC:-cnBC]
uaend = uabc[-cnBC:]
    
Ga_c = copyarraytoC(Ga)
Gabeg_c = copyarraytoC(Gabeg)
Gaend_c = copyarraytoC(Gaend)
ha_c = copyarraytoC(ha)
  
habeg_c = copyarraytoC(habeg)
haend_c = copyarraytoC(haend)
  
uabeg_c = copyarraytoC(uabeg)
uaend_c = copyarraytoC(uaend)
   
hmbeg_c = copyarraytoC(hmbeg)
hmend_c = copyarraytoC(hmend)
    
umbeg_c = copyarraytoC(umbeg)
umend_c = copyarraytoC(umend)
    
u_c= mallocPy(n)
G_c = mallocPy(n)
h_c = mallocPy(n)
    
for i in range(1,len(t)):
    if(i ==1 or i %gap == 0):
        ca2midpt(ha_c,dx,n,h_c)
        ca2midpt(Ga_c,dx,n,G_c)
            
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        
        Ga = copyarrayfromC(Ga_c,n)
        ha = copyarrayfromC(ha_c,n)
            
        cnBC = nfcBC
        Gabc = concatenate([Gabeg[-cnBC:],Ga,Gaend[0:cnBC]])  
        habc = concatenate([habeg[-cnBC:],ha,haend[0:cnBC]]) 
        Gabc_c = copyarraytoC(Gabc)
        habc_c = copyarraytoC(habc)
        ufromGh(Gabc_c,habc_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,cnBC, u_c)
        u = copyarrayfromC(u_c,n)
        s = wdir + "out" + str(i) +  ".txt"
        with open(s,'a') as file2:
             writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
             writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'G' , 'u(m/s)'])        
                           
             for j in range(n):
                 writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(G[j]) , str(u[j])])    
    evolvewrap(Ga_c,ha_c,Gabeg_c,Gaend_c,habeg_c,haend_c,hmbeg_c,hmend_c,uabeg_c,uaend_c,umbeg_c,umend_c,nfcBC,nGsBC,g,dx,dt,n,cnBC,niBC)
    print (t[i])
    
        
ca2midpt(ha_c,dx,n,h_c)
ca2midpt(Ga_c,dx,n,G_c)
            
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)
        
Ga = copyarrayfromC(Ga_c,n)
ha = copyarrayfromC(ha_c,n)
            
cnBC = nfcBC
Gabc = concatenate([Gabeg[-cnBC:],Ga,Gaend[0:cnBC]])  
habc = concatenate([habeg[-cnBC:],ha,haend[0:cnBC]]) 
Gabc_c = copyarraytoC(Gabc)
habc_c = copyarraytoC(habc)
ufromGh(Gabc_c,habc_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,cnBC, u_c)
u = copyarrayfromC(u_c,n)
s = wdir + "outlast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'G' , 'u(m/s)'])        
                   
     for j in range(n):
         writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(G[j]) , str(u[j])])    
deallocPy(u_c)   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(ha_c)
deallocPy(Ga_c)
deallocPy(habeg_c)
deallocPy(Gabeg_c)
deallocPy(haend_c)
deallocPy(Gaend_c)
deallocPy(uabeg_c)
deallocPy(uaend_c)
deallocPy(umbeg_c)
deallocPy(umend_c)
deallocPy(hmbeg_c)
deallocPy(hmend_c)
"""

"""
#Dam Break Chris Smooth
import os

wdir = "../../../data/dbchris/o3femh110h01testsmooth/"

if not os.path.exists(wdir):
    os.makedirs(wdir)

g = 9.81
hf = 10.0
hl = 1.0
    
dx = 0.1
Cr = 0.2
l = Cr / sqrt(g*hf)
dt = l*dx
startx = 0.0
endx = 1000.0 + dx
startt = 0.0
endt = 2*dt
        
szoomx = startx
ezoomx = endx
        
#number of boundary conditions (one side)
nfcBC = 4 #for flux calculation
nGsBC = 2 #for solving G from u,h
niBC = nGsBC + nfcBC #total
        
g = 9.81
    
gap = int(0.01/dt)
        
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
    
xc = 500   
#um,hm = dambreaksmooth(x,xc,hl,hf-hl,1.0,dx)
um,hm = dambreak(x,xc,hf,hl)
                
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
        
#calculate G midpoint
cnBC = niBC - nGsBC
        
umbc = concatenate([umbeg[-cnBC:],um,umend[0:cnBC]]) 
hmbc = concatenate([hmbeg[-cnBC:],hm,hmend[0:cnBC]])       
Gmbc = solveGfromuh(umbc,hmbc,hmbeg[0:-cnBC],hmend[-cnBC:],umbeg[0:-cnBC],umend[-cnBC:],dx)  
        
#calculate averages
Gabc = midpointtocellaverages(Gmbc,dx)
habc = midpointtocellaverages(hmbc,dx)
uabc = midpointtocellaverages(umbc,dx)
        
#so we can just go from here with Ga ang ha?
Gabeg = Gabc[0:cnBC]
Ga = Gabc[cnBC:-cnBC]
Gaend = Gabc[-cnBC:]
habeg = habc[0:cnBC]
ha = habc[cnBC:-cnBC]
haend = habc[-cnBC:]
uabeg = uabc[0:cnBC]
ua = uabc[cnBC:-cnBC]
uaend = uabc[-cnBC:]
    
Ga_c = copyarraytoC(Ga)
Gabeg_c = copyarraytoC(Gabeg)
Gaend_c = copyarraytoC(Gaend)
ha_c = copyarraytoC(ha)
  
habeg_c = copyarraytoC(habeg)
haend_c = copyarraytoC(haend)
  
uabeg_c = copyarraytoC(uabeg)
uaend_c = copyarraytoC(uaend)
   
hmbeg_c = copyarraytoC(hmbeg)
hmend_c = copyarraytoC(hmend)
    
umbeg_c = copyarraytoC(umbeg)
umend_c = copyarraytoC(umend)
    
u_c= mallocPy(2*n+1)
G_c = mallocPy(n)
h_c = mallocPy(n)
    
for i in range(1,len(t)):
    if(i ==1 or i %gap == 0):
        ca2midpt(ha_c,dx,n,h_c)
        ca2midpt(Ga_c,dx,n,G_c)
            
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        
        Ga = copyarrayfromC(Ga_c,n)
        ha = copyarrayfromC(ha_c,n)
            
        cnBC = nfcBC
        Gabc = concatenate([Gabeg[-cnBC:],Ga,Gaend[0:cnBC]])  
        habc = concatenate([habeg[-cnBC:],ha,haend[0:cnBC]]) 
        Gabc_c = copyarraytoC(Gabc)
        habc_c = copyarraytoC(habc)
        ufromGh(Gabc_c,habc_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,cnBC, u_c)
        u = copyarrayfromC(u_c,2*n+1)
        u = u[1::2]
        s = wdir + "out" + str(i) +  ".txt"
        with open(s,'a') as file2:
             writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
             writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'G' , 'u(m/s)'])        
                           
             for j in range(n):
                 writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(G[j]) , str(u[j])])    
    evolvewrap(Ga_c,ha_c,Gabeg_c,Gaend_c,habeg_c,haend_c,hmbeg_c,hmend_c,uabeg_c,uaend_c,umbeg_c,umend_c,nfcBC,nGsBC,g,dx,dt,n,cnBC,niBC)
    print (t[i])
    
        
ca2midpt(ha_c,dx,n,h_c)
ca2midpt(Ga_c,dx,n,G_c)
            
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)
        
Ga = copyarrayfromC(Ga_c,n)
ha = copyarrayfromC(ha_c,n)
            
cnBC = nfcBC
Gabc = concatenate([Gabeg[-cnBC:],Ga,Gaend[0:cnBC]])  
habc = concatenate([habeg[-cnBC:],ha,haend[0:cnBC]]) 
Gabc_c = copyarraytoC(Gabc)
habc_c = copyarraytoC(habc)
ufromGh(Gabc_c,habc_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,cnBC, u_c)
u = copyarrayfromC(u_c,2*n+1)
u = u[1::2]
s = wdir + "outlast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'G' , 'u(m/s)'])        
                   
     for j in range(n):
         writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(G[j]) , str(u[j])])  

deallocPy(u_c)   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(ha_c)
deallocPy(Ga_c)
deallocPy(habeg_c)
deallocPy(Gabeg_c)
deallocPy(haend_c)
deallocPy(Gaend_c)
deallocPy(uabeg_c)
deallocPy(uaend_c)
deallocPy(umbeg_c)
deallocPy(umend_c)
deallocPy(hmbeg_c)
deallocPy(hmend_c)
"""
"""
#Dam Break Chris
import os

wdir = "../../../data/raw/Cserre/db/o3femhf10hl0p9/"

if not os.path.exists(wdir):
    os.makedirs(wdir)

g = 9.81
hf = 10.0
hl = 1.0
    
dx = 0.1
Cr = 0.2
l = Cr / sqrt(g*hf)
dt = l*dx
startx = 0.0
endx = 1000.0 + dx
startt = 0.0
endt = 30.0 + dt
        
szoomx = startx
ezoomx = endx
        
#number of boundary conditions (one side)
nfcBC = 4 #for flux calculation
nGsBC = 2 #for solving G from u,h
niBC = nGsBC + nfcBC #total
        
g = 9.81
    
gap = int(10.0/dt)
        
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
    
xc = 500   
um,hm = dambreak(x,xc,hf,hl)
                
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
        
#calculate G midpoint
cnBC = niBC - nGsBC
        
umbc = concatenate([umbeg[-cnBC:],um,umend[0:cnBC]]) 
hmbc = concatenate([hmbeg[-cnBC:],hm,hmend[0:cnBC]])       
Gmbc = solveGfromuh(umbc,hmbc,hmbeg[0:-cnBC],hmend[-cnBC:],umbeg[0:-cnBC],umend[-cnBC:],dx)  
        
#calculate averages
Gabc = midpointtocellaverages(Gmbc,dx)
habc = midpointtocellaverages(hmbc,dx)
uabc = midpointtocellaverages(umbc,dx)
        
#so we can just go from here with Ga ang ha?
Gabeg = Gabc[0:cnBC]
Ga = Gabc[cnBC:-cnBC]
Gaend = Gabc[-cnBC:]
habeg = habc[0:cnBC]
ha = habc[cnBC:-cnBC]
haend = habc[-cnBC:]
uabeg = uabc[0:cnBC]
ua = uabc[cnBC:-cnBC]
uaend = uabc[-cnBC:]
    
Ga_c = copyarraytoC(Ga)
Gabeg_c = copyarraytoC(Gabeg)
Gaend_c = copyarraytoC(Gaend)
ha_c = copyarraytoC(ha)
  
habeg_c = copyarraytoC(habeg)
haend_c = copyarraytoC(haend)
  
uabeg_c = copyarraytoC(uabeg)
uaend_c = copyarraytoC(uaend)
   
hmbeg_c = copyarraytoC(hmbeg)
hmend_c = copyarraytoC(hmend)
    
umbeg_c = copyarraytoC(umbeg)
umend_c = copyarraytoC(umend)
    
u_c= mallocPy(2*n+1)
G_c = mallocPy(n)
h_c = mallocPy(n)
    
for i in range(1,len(t)):
    if(i ==1 or i %gap == 0):
        ca2midpt(ha_c,dx,n,h_c)
        ca2midpt(Ga_c,dx,n,G_c)
            
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        
        Ga = copyarrayfromC(Ga_c,n)
        ha = copyarrayfromC(ha_c,n)
            
        cnBC = nfcBC
        Gabc = concatenate([Gabeg[-cnBC:],Ga,Gaend[0:cnBC]])  
        habc = concatenate([habeg[-cnBC:],ha,haend[0:cnBC]]) 
        Gabc_c = copyarraytoC(Gabc)
        habc_c = copyarraytoC(habc)
        ufromGh(Gabc_c,habc_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,cnBC, u_c)
        u = copyarrayfromC(u_c,2*n+1)
        u = u[1::2]
        s = wdir + "out" + str(i) +  ".txt"
        with open(s,'a') as file2:
             writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
             writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'G' , 'u(m/s)'])        
                           
             for j in range(n):
                 writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(G[j]) , str(u[j])])    
    evolvewrap(Ga_c,ha_c,Gabeg_c,Gaend_c,habeg_c,haend_c,hmbeg_c,hmend_c,uabeg_c,uaend_c,umbeg_c,umend_c,nfcBC,nGsBC,g,dx,dt,n,cnBC,niBC)
    print (t[i])
    
        
ca2midpt(ha_c,dx,n,h_c)
ca2midpt(Ga_c,dx,n,G_c)
            
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)
        
Ga = copyarrayfromC(Ga_c,n)
ha = copyarrayfromC(ha_c,n)
            
cnBC = nfcBC
Gabc = concatenate([Gabeg[-cnBC:],Ga,Gaend[0:cnBC]])  
habc = concatenate([habeg[-cnBC:],ha,haend[0:cnBC]]) 
Gabc_c = copyarraytoC(Gabc)
habc_c = copyarraytoC(habc)
ufromGh(Gabc_c,habc_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,cnBC, u_c)
u = copyarrayfromC(u_c,2*n+1)
u = u[1::2]
s = wdir + "outlast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'G' , 'u(m/s)'])        
                   
     for j in range(n):
         writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(G[j]) , str(u[j])])  

deallocPy(u_c)   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(ha_c)
deallocPy(Ga_c)
deallocPy(habeg_c)
deallocPy(Gabeg_c)
deallocPy(haend_c)
deallocPy(Gaend_c)
deallocPy(uabeg_c)
deallocPy(uaend_c)
deallocPy(umbeg_c)
deallocPy(umend_c)
deallocPy(hmbeg_c)
deallocPy(hmend_c)
"""

"""
### DAMBREAK ACCURACY
g = 9.81
Cr = 0.2
hf = 1.8
hl = 1.0
wdir = "../../data/o3fem/db/"
for k in range(16):
    dx = 100.0/(2**k)
    l = 0.01
    dt = l*dx
    startx = 0.0
    endx = 1000.0 + dx
    startt = 0.0
    endt = 30 + dt
            
    szoomx = startx
    ezoomx = endx
            
    #number of boundary conditions (one side)
    nfcBC = 4 #for flux calculation
    nGsBC = 2 #for solving G from u,h
    niBC = nGsBC + nfcBC #total
            
    g = 9.81
        
    gap = int(0.5/dt)
            
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
        
    xc = 500
    hf = 1.8
    hl = 1.0    
    um,hm = dambreak(x,xc,hf,hl)
                    
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
            
    #calculate G midpoint
    cnBC = niBC - nGsBC
            
    umbc = concatenate([umbeg[-cnBC:],um,umend[0:cnBC]]) 
    hmbc = concatenate([hmbeg[-cnBC:],hm,hmend[0:cnBC]])       
    Gmbc = solveGfromuh(umbc,hmbc,hmbeg[0:-cnBC],hmend[-cnBC:],umbeg[0:-cnBC],umend[-cnBC:],dx)  
            
    #calculate averages
    Gabc = midpointtocellaverages(Gmbc,dx)
    habc = midpointtocellaverages(hmbc,dx)
    uabc = midpointtocellaverages(umbc,dx)
            
    #so we can just go from here with Ga ang ha?
    Gabeg = Gabc[0:cnBC]
    Ga = Gabc[cnBC:-cnBC]
    Gaend = Gabc[-cnBC:]
    habeg = habc[0:cnBC]
    ha = habc[cnBC:-cnBC]
    haend = habc[-cnBC:]
    uabeg = uabc[0:cnBC]
    ua = uabc[cnBC:-cnBC]
    uaend = uabc[-cnBC:]
        
    Ga_c = copyarraytoC(Ga)
    Gabeg_c = copyarraytoC(Gabeg)
    Gaend_c = copyarraytoC(Gaend)
    ha_c = copyarraytoC(ha)
      
    habeg_c = copyarraytoC(habeg)
    haend_c = copyarraytoC(haend)
      
    uabeg_c = copyarraytoC(uabeg)
    uaend_c = copyarraytoC(uaend)
       
    hmbeg_c = copyarraytoC(hmbeg)
    hmend_c = copyarraytoC(hmend)
        
    umbeg_c = copyarraytoC(umbeg)
    umend_c = copyarraytoC(umend)
        
    u_c= mallocPy(n)
    G_c = mallocPy(n)
    h_c = mallocPy(n)
        
    for i in range(1,len(t)):
        evolvewrap(Ga_c,ha_c,Gabeg_c,Gaend_c,habeg_c,haend_c,hmbeg_c,hmend_c,uabeg_c,uaend_c,umbeg_c,umend_c,nfcBC,nGsBC,g,dx,dt,n,cnBC,niBC)
        print (t[i])
        
            
    ca2midpt(ha_c,dx,n,h_c)
    ca2midpt(Ga_c,dx,n,G_c)
                
    G = copyarrayfromC(G_c,n)
    h = copyarrayfromC(h_c,n)
            
    Ga = copyarrayfromC(Ga_c,n)
    ha = copyarrayfromC(ha_c,n)
                
    cnBC = nfcBC
    Gabc = concatenate([Gabeg[-cnBC:],Ga,Gaend[0:cnBC]])  
    habc = concatenate([habeg[-cnBC:],ha,haend[0:cnBC]]) 
    Gabc_c = copyarraytoC(Gabc)
    habc_c = copyarraytoC(habc)
    ufromGh(Gabc_c,habc_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,cnBC, u_c)
    u = copyarrayfromC(u_c,n)
    s = wdir + str(k) + ".txt"
    with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
         writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'G' , 'u(m/s)'])        
                       
         for j in range(n):
             writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(G[j]) , str(u[j])])    
    deallocPy(u_c)   
    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(ha_c)
    deallocPy(Ga_c)
    deallocPy(habeg_c)
    deallocPy(Gabeg_c)
    deallocPy(haend_c)
    deallocPy(Gaend_c)
    deallocPy(uabeg_c)
    deallocPy(uaend_c)
    deallocPy(umbeg_c)
    deallocPy(umend_c)
    deallocPy(hmbeg_c)
    deallocPy(hmend_c)
"""    

"""
#soliton 
wdir = "../../data/test11/o3fem/sol/"
if not os.path.exists(wdir):
    os.makedirs(wdir)
dx = 0.5
a0 = 10.0
a1 = 1.0
g = 9.81
l = 0.5/sqrt(g*(a0 + a1))
dt = l*dx
startx = -500.0
endx = 1500.0 + dx
startt = 0.0#19.8065574333 - dt
endt = 20 #822*dt
        
szoomx = startx
ezoomx = endx
        
#number of boundary conditions (one side)
nfcBC = 4 #for flux calculation
nGsBC = 2 #for solving G from u,h
niBC = nGsBC + nfcBC #total
        
g = 9.81
    
gap = int(0.5/dt)
        
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
      
hm,um = solitoninit(n,a0,a1,g,x,0.0,dx)
                
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
        
#calculate G midpoint
cnBC = niBC - nGsBC
        
umbc = concatenate([umbeg[-cnBC:],um,umend[0:cnBC]]) 
hmbc = concatenate([hmbeg[-cnBC:],hm,hmend[0:cnBC]])       
Gmbc = solveGfromuh(umbc,hmbc,hmbeg[0:-cnBC],hmend[-cnBC:],umbeg[0:-cnBC],umend[-cnBC:],dx)  
        
#calculate averages
Gabc = midpointtocellaverages(Gmbc,dx)
habc = midpointtocellaverages(hmbc,dx)
uabc = midpointtocellaverages(umbc,dx)
        
#so we can just go from here with Ga ang ha?
Gabeg = Gabc[0:cnBC]
Ga = Gabc[cnBC:-cnBC]
Gaend = Gabc[-cnBC:]
habeg = habc[0:cnBC]
ha = habc[cnBC:-cnBC]
haend = habc[-cnBC:]
uabeg = uabc[0:cnBC]
ua = uabc[cnBC:-cnBC]
uaend = uabc[-cnBC:]
    
Ga_c = copyarraytoC(Ga)
Gabeg_c = copyarraytoC(Gabeg)
Gaend_c = copyarraytoC(Gaend)
ha_c = copyarraytoC(ha)
  
habeg_c = copyarraytoC(habeg)
haend_c = copyarraytoC(haend)
  
uabeg_c = copyarraytoC(uabeg)
uaend_c = copyarraytoC(uaend)
   
hmbeg_c = copyarraytoC(hmbeg)
hmend_c = copyarraytoC(hmend)
    
umbeg_c = copyarraytoC(umbeg)
umend_c = copyarraytoC(umend)
    
u_c= mallocPy(2*n + 1)
G_c = mallocPy(n)
h_c = mallocPy(n)
    
for i in range(1,len(t)):
    if(i ==1 or i %gap == 0):
        ca2midpt(ha_c,dx,n,h_c)
        ca2midpt(Ga_c,dx,n,G_c)
            
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        
        Ga = copyarrayfromC(Ga_c,n)
        ha = copyarrayfromC(ha_c,n)
            
        cnBC = nfcBC
        Gabc = concatenate([Gabeg[-cnBC:],Ga,Gaend[0:cnBC]])  
        habc = concatenate([habeg[-cnBC:],ha,haend[0:cnBC]]) 
        Gabc_c = copyarraytoC(Gabc)
        habc_c = copyarraytoC(habc)
        ufromGh(Gabc_c,habc_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,cnBC, u_c)
        u = copyarrayfromC(u_c,2*n + 1)
        u = u[1::2]
        
        c = sqrt(g*(a0 + a1))
        htrue = zeros(n)
        utrue = zeros(n)
        for j in range(n):             
            he = soliton(x[j],t[-1],g,a0,a1)
            htrue[j] = he
            utrue[j] = c* ((he - a0) / he)
        s = wdir + "out" + str(i) +  ".txt"
        with open(s,'a') as file2:
             writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
             writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'G' , 'u(m/s)','h(m) exact', 'u exact'])        
                           
             for j in range(n):
                 writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(G[j]) , str(u[j]), str(htrue[j]), str(utrue[j])])    
    evolvewrap(Ga_c,ha_c,Gabeg_c,Gaend_c,habeg_c,haend_c,hmbeg_c,hmend_c,uabeg_c,uaend_c,umbeg_c,umend_c,nfcBC,nGsBC,g,dx,dt,n,cnBC,niBC)
    print (t[i])
    
        
ca2midpt(ha_c,dx,n,h_c)
ca2midpt(Ga_c,dx,n,G_c)
            
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)
        
Ga = copyarrayfromC(Ga_c,n)
ha = copyarrayfromC(ha_c,n)
            
cnBC = nfcBC
Gabc = concatenate([Gabeg[-cnBC:],Ga,Gaend[0:cnBC]])  
habc = concatenate([habeg[-cnBC:],ha,haend[0:cnBC]]) 
Gabc_c = copyarraytoC(Gabc)
habc_c = copyarraytoC(habc)
ufromGh(Gabc_c,habc_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,cnBC, u_c)
u = copyarrayfromC(u_c,2*n + 1)
u = u[1::2]

c = sqrt(g*(a0 + a1))
htrue = zeros(n)
utrue = zeros(n)
for j in range(n):             
    he = soliton(x[j],t[-1],g,a0,a1)
    htrue[j] = he
    utrue[j] = c* ((he - a0) / he)
s = wdir + "outlast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
     writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'G' , 'u(m/s)','h(m) exact', 'u exact'])        
                           
     for j in range(n):
         writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(G[j]) , str(u[j]), str(htrue[j]), str(utrue[j])])   
deallocPy(u_c)   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(ha_c)
deallocPy(Ga_c)
deallocPy(habeg_c)
deallocPy(Gabeg_c)
deallocPy(haend_c)
deallocPy(Gaend_c)
deallocPy(uabeg_c)
deallocPy(uaend_c)
deallocPy(umbeg_c)
deallocPy(umend_c)
deallocPy(hmbeg_c)
deallocPy(hmend_c)
"""

"""
### Soliton accuracy
import os

wdir = "../../../data/femelemdubreak/"

if not os.path.exists(wdir):
    os.makedirs(wdir)

s = wdir + "savenorms.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(['dx','Normalised L1-norm Difference Height', ' Normalised L1-norm Difference Velocity'])
for k in range(20):
    g = 9.81
    dx = 100.0 / (2**k)
    a0 = 10.0
    a1 = 1.0
    Cr = 0.2
    l = 0.01#Cr / sqrt(g*(a0 + a1))
    dt = l*dx
    startx = -500.0
    endx = 1500.0
    startt = 0.0
    endt = 2*dt
            
    szoomx = startx
    ezoomx = endx
            
    #number of boundary conditions (one side)
    nfcBC = 4 #for flux calculation
    nGsBC = 2 #for solving G from u,h
    niBC = nGsBC + nfcBC #total
            
        
    gap = int(0.5/dt)
            
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
          
    hm,um = solitoninit(n,a0,a1,g,x,0.0,dx)
                    
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
            
    #calculate G midpoint
    cnBC = niBC - nGsBC
            
    umbc = concatenate([umbeg[-cnBC:],um,umend[0:cnBC]]) 
    hmbc = concatenate([hmbeg[-cnBC:],hm,hmend[0:cnBC]])       
    Gmbc = solveGfromuh(umbc,hmbc,hmbeg[0:-cnBC],hmend[-cnBC:],umbeg[0:-cnBC],umend[-cnBC:],dx)  
            
    #calculate averages
    Gabc = midpointtocellaverages(Gmbc,dx)
    habc = midpointtocellaverages(hmbc,dx)
            
    #so we can just go from here with Ga ang ha?
    Gabeg = Gabc[0:cnBC]
    Ga = Gabc[cnBC:-cnBC]
    Gaend = Gabc[-cnBC:]
    habeg = habc[0:cnBC]
    ha = habc[cnBC:-cnBC]
    haend = habc[-cnBC:]
    
    
    uabeg = um[0]*ones(2*cnBC)
    uaend = um[-1]*ones(2*cnBC)
        
    Ga_c = copyarraytoC(Ga)
    Gabeg_c = copyarraytoC(Gabeg)
    Gaend_c = copyarraytoC(Gaend)
    ha_c = copyarraytoC(ha)
      
    habeg_c = copyarraytoC(habeg)
    haend_c = copyarraytoC(haend)
      
    uabeg_c = copyarraytoC(uabeg)
    uaend_c = copyarraytoC(uaend)
       
    hmbeg_c = copyarraytoC(hmbeg)
    hmend_c = copyarraytoC(hmend)
        
    umbeg_c = copyarraytoC(umbeg)
    umend_c = copyarraytoC(umend)
        
    ue_c= mallocPy(2*n + 1)
    G_c = mallocPy(n)
    h_c = mallocPy(n)   
    for i in range(1,len(t)):
        evolvewrap(Ga_c,ha_c,Gabeg_c,Gaend_c,habeg_c,haend_c,hmbeg_c,hmend_c,uabeg_c,uaend_c,umbeg_c,umend_c,nfcBC,nGsBC,g,dx,dt,n,cnBC,niBC)
        print (t[i])
        
    #print("end of loop")        
    ca2midpt(ha_c,dx,n,h_c)
    ca2midpt(Ga_c,dx,n,G_c)
                
    G = copyarrayfromC(G_c,n)
    h = copyarrayfromC(h_c,n)
            
    Ga = copyarrayfromC(Ga_c,n)
    ha = copyarrayfromC(ha_c,n)
                
    cnBC = nfcBC
    Gabc = concatenate([Gabeg[-cnBC:],Ga,Gaend[0:cnBC]])  
    habc = concatenate([habeg[-cnBC:],ha,haend[0:cnBC]]) 
    Gabc_c = copyarraytoC(Gabc)
    habc_c = copyarraytoC(habc)
    ufromGh(Gabc_c,habc_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,cnBC, ue_c)
    ue = copyarrayfromC(ue_c,2*n + 1)
    u = ue[1::2]
    
    c = sqrt(g*(a0 + a1))
    htrue = zeros(n)
    utrue = zeros(n)
    for j in range(n):             
        he = soliton(x[j],t[-1],g,a0,a1)
        htrue[j] = he
        utrue[j] = c* ((he - a0) / he)
    s = wdir + str(k) +".txt"
    with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                
         writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'G' , 'u(m/s)','h(m) exact', 'u exact'])        
                               
         for j in range(n):
             writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(G[j]) , str(u[j]), str(htrue[j]), str(utrue[j])])
    normhdiffi = norm(h - htrue,ord=1) / norm(htrue,ord=1)
    normudiffi = norm(u -utrue,ord=1) / norm(utrue,ord=1)  

    s = wdir + "savenorms.txt"
    with open(s,'a') as file1:
        writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

        writefile.writerow([str(dx),str(normhdiffi), str(normudiffi)]) 
    deallocPy(ue_c)   
    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(ha_c)
    deallocPy(Ga_c)
    deallocPy(habeg_c)
    deallocPy(Gabeg_c)
    deallocPy(haend_c)
    deallocPy(Gaend_c)
    deallocPy(uabeg_c)
    deallocPy(uaend_c)
    deallocPy(umbeg_c)
    deallocPy(umend_c)
    deallocPy(hmbeg_c)
    deallocPy(hmend_c)

"""
"""
#Segur and Hammack
import os
wdir = "../../../data/seguro3fem/"

if not os.path.exists(wdir):
    os.makedirs(wdir)

tl = 60.0
b = 0.61
h0 = 0.09
h1 = 0.1
g = 9.81

dx = 0.01
Cr = 0.2
l = Cr / sqrt(g*h1)
dt = l*dx
startx = -tl
endx = tl + dx
startt = 0.0
endt = 50.0 + dt  
        
szoomx = startx
ezoomx = endx
        
#number of boundary conditions (one side)
nfcBC = 4 #for flux calculation
nGsBC = 2 #for solving G from u,h
niBC = nGsBC + nfcBC #total
        
    
gap = 1
        
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)

   
hm,um = experiment1(x,b,h0,h1,dx)
                
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
        
#calculate G midpoint
cnBC = niBC - nGsBC
        
umbc = concatenate([umbeg[-cnBC:],um,umend[0:cnBC]]) 
hmbc = concatenate([hmbeg[-cnBC:],hm,hmend[0:cnBC]])       
Gmbc = solveGfromuh(umbc,hmbc,hmbeg[0:-cnBC],hmend[-cnBC:],umbeg[0:-cnBC],umend[-cnBC:],dx)  
        
#calculate averages
Gabc = midpointtocellaverages(Gmbc,dx)
habc = midpointtocellaverages(hmbc,dx)
uabc = midpointtocellaverages(umbc,dx)
        
#so we can just go from here with Ga ang ha?
Gabeg = Gabc[0:cnBC]
Ga = Gabc[cnBC:-cnBC]
Gaend = Gabc[-cnBC:]
habeg = habc[0:cnBC]
ha = habc[cnBC:-cnBC]
haend = habc[-cnBC:]
uabeg = uabc[0:cnBC]
ua = uabc[cnBC:-cnBC]
uaend = uabc[-cnBC:]
    
Ga_c = copyarraytoC(Ga)
Gabeg_c = copyarraytoC(Gabeg)
Gaend_c = copyarraytoC(Gaend)
ha_c = copyarraytoC(ha)
  
habeg_c = copyarraytoC(habeg)
haend_c = copyarraytoC(haend)
  
uabeg_c = copyarraytoC(uabeg)
uaend_c = copyarraytoC(uaend)
   
hmbeg_c = copyarraytoC(hmbeg)
hmend_c = copyarraytoC(hmend)
    
umbeg_c = copyarraytoC(umbeg)
umend_c = copyarraytoC(umend)
    
u_c= mallocPy(n)
G_c = mallocPy(n)
h_c = mallocPy(n)
    
for i in range(1,len(t)):
    if(i ==1 or i %gap == 0):
        ca2midpt(ha_c,dx,n,h_c)
        ca2midpt(Ga_c,dx,n,G_c)
            
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        
        Ga = copyarrayfromC(Ga_c,n)
        ha = copyarrayfromC(ha_c,n)
            
        cnBC = nfcBC
        Gabc = concatenate([Gabeg[-cnBC:],Ga,Gaend[0:cnBC]])  
        habc = concatenate([habeg[-cnBC:],ha,haend[0:cnBC]]) 
        Gabc_c = copyarraytoC(Gabc)
        habc_c = copyarraytoC(habc)
        ufromGh(Gabc_c,habc_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,cnBC, u_c)
        u = copyarrayfromC(u_c,n)
        s = wdir + "out" + str(i) +  ".txt"
        with open(s,'a') as file2:
             writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
             writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'G' , 'u(m/s)'])        
                           
             for j in range(n):
                 writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(G[j]) , str(u[j])])    
    evolvewrap(Ga_c,ha_c,Gabeg_c,Gaend_c,habeg_c,haend_c,hmbeg_c,hmend_c,uabeg_c,uaend_c,umbeg_c,umend_c,nfcBC,nGsBC,g,dx,dt,n,cnBC,niBC)
    print (t[i])
    
        
ca2midpt(ha_c,dx,n,h_c)
ca2midpt(Ga_c,dx,n,G_c)
            
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)
        
Ga = copyarrayfromC(Ga_c,n)
ha = copyarrayfromC(ha_c,n)
            
cnBC = nfcBC
Gabc = concatenate([Gabeg[-cnBC:],Ga,Gaend[0:cnBC]])  
habc = concatenate([habeg[-cnBC:],ha,haend[0:cnBC]]) 
Gabc_c = copyarraytoC(Gabc)
habc_c = copyarraytoC(habc)
ufromGh(Gabc_c,habc_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,cnBC, u_c)
u = copyarrayfromC(u_c,n)
s = wdir + "outlast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     writefile2.writerow(['dx' ,'dt','time','cell midpoint', 'height(m)', 'G' , 'u(m/s)'])        
                   
     for j in range(n):
         writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(G[j]) , str(u[j])])    
deallocPy(u_c)   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(ha_c)
deallocPy(Ga_c)
deallocPy(habeg_c)
deallocPy(Gabeg_c)
deallocPy(haend_c)
deallocPy(Gaend_c)
deallocPy(uabeg_c)
deallocPy(uaend_c)
deallocPy(umbeg_c)
deallocPy(umend_c)
deallocPy(hmbeg_c)
deallocPy(hmend_c)
""" 