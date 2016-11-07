# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:09:21 2015

@author: jordan
"""
from Serre2 import *
from scipy import *
import csv
import os
from numpy.linalg import norm  

from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv

def minmodpy(a,b,c):
    if (a>0 and b>0 and c>0):
        return min(a,b,c)
    elif(a<0 and b<0 and c<0):
        return max(a,b,c)
    else:
        return 0
          

#algorithm for solving tridiagonal systems from wiki page    
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

#gives exact up to linears, so is second order accurate huzzah
def getufromGpy(con,bed,u0,u1,h0,h1,b0,b1,dx):
    idx = 1.0 / dx
    ithree = 1.0 / 3.0
    n = len(con)
    
    a = zeros(n-1)
    b = zeros(n)
    c = zeros(n-1)
    
    G = zeros(n)
    
    for i in range(n):
        G[i] = con[i][1]
    
    
    for i in range(1,n-1):
        th = con[i][0]
        thx = 0.5*idx*(con[i+1][0] - con[i-1][0])
        tbx = 0.5*idx*(bed[i+1] - bed[i-1])
        tbxx = idx*idx*(bed[i+1] -2*bed[i] + bed[i-1])
        
        D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
        
        ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
        bi = D + 2.0*ithree*idx*idx*th*th*th
        ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
        
        
        a[i-1] = ai
        b[i] =  bi
        c[i] = ci
        
    #boundary    
    #i=0
    i=0
    th = con[i][0]
    thx = 0.5*idx*(con[i+1][0] - h0)
    tbx = 0.5*idx*(bed[i+1] - b0)
    tbxx = idx*idx*(bed[i+1] -2*bed[i] + b0)
            
    D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
            
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
 
    c[i] = ci
    b[i] = bi
    
    G[i] = G[i] - u0*ai
    
    #i = n-1
    i = n-1
    
    
    th = con[i][0]
    thx = 0.5*idx*(h1 - con[i-1][0])
    tbx = 0.5*idx*(b1 - bed[i-1])
    tbxx = idx*idx*(b1 -2*bed[i] + bed[i-1])
        
    D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
        
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx

    a[i-1] = ai
    b[i] = bi
    G[i] = G[i] - u1*ci
    
    u = TDMApy(a,b,c,G)
        
    return u 

#gives exact up to linears, so is second order accurate huzzah    
def getGfromupy(con,u,bed,u0,u1,h0,h1,b0,b1,dx):
    idx = 1.0 / dx
    ithree = 1.0 / 3.0
        
    n = len(con)

    G = zeros(n)
        
    for i in range(1,n-1):
        th = con[i][0]
        thx = 0.5*idx*(con[i+1][0] - con[i-1][0])
        tbx = 0.5*idx*(bed[i+1] - bed[i-1])
        tbxx = idx*idx*(bed[i+1] -2*bed[i] + bed[i-1])
        
        D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
        
        ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
        bi = D + 2.0*ithree*idx*idx*th*th*th
        ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
        
        G[i] = ai*u[i-1] + bi*u[i] + ci*u[i+1]
        
    #boundary    
    #i=0
    i=0
    th = con[i][0]
    thx = 0.5*idx*(con[i+1][0] - h0)
    tbx = 0.5*idx*(bed[i+1] - b0)
    tbxx = idx*idx*(bed[i+1] -2*bed[i] + b0)
            
    D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
            
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
 
    G[i] = ai*u0 + bi*u[i] + ci*u[i+1]
    
    #i = n-1
    i = n-1

    th = con[i][0]
    thx = 0.5*idx*(h1 - con[i-1][0])
    tbx = 0.5*idx*(b1 - bed[i-1])
    tbxx = idx*idx*(b1 -2*bed[i] + bed[i-1])
        
    D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
        
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx

    G[i] = ai*u[i-1] + bi*u[i] + ci*u1
            
    return G 
    


def evolvetwopy(con,bed,g,u0,u1,h0,h1,b0,b1,beta,dx,dt): 
    n = len(con)
    ncon = zeros((n,2))
    #update h' and G'    
    conp = evolvepy(con,bed,g,beta,u0,u1,h0,h1,b0,b1,dx,dt)
        
    #update h'' and G''
    conpp = evolvepy(conp,bed,g,beta,u0,u1,h0,h1,b0,b1,dx,dt)
    
    ncon = 0.5*(conpp + con)
    #ncon = conp
    return ncon
    
    
def evolvepy(con,bed,g,beta,u0,u1,h0,h1,b0,b1,dx,dt):
    #get averages
    idx = 1.0 / dx  
    ithree = 1.0 / 3.0
    n = len(con)
    
    ncon = zeros((n,2))
    
    #calculate u at time step
    u = getufromGpy(con,bed,u0[-1],u1[0],h0[-1],h1[0],b0[-1],b1[0],dx)
    
    
    #boundaries
    #beginning
    #i=-1
    i = -1
    h = h0[-1]
    hx = 0.5*idx*(con[i+1][0] - h0[-2])
    bx = 0.5*idx*(bed[i+1] - b0[-2])
    bxx = idx*idx*(bed[i+1] - 2*b0[-1] + b0[-2])
    
    D = h + h*hx*bx + 0.5*h*h*bxx + h*bx*bx
    
    ai = h*h*hx*(0.5*idx) - ithree*h*h*h*(idx*idx)
    bi = D + 2.0*ithree*h*h*h*(idx*idx)
    ci = -1.0*h*h*hx*(0.5*idx) - ithree*h*h*h*(idx*idx)
    
    gb3 = ai*u0[-2] + bi*u0[-1] +ci*u[i+1]
    
    #i=-2
    h = h0[-2]
    hx = 0.5*idx*(h0[-1] - h0[-3])
    bx = 0.5*idx*(b0[-1] - b0[-3])
    bxx = idx*idx*(b0[-1] - 2*b0[-2] + b0[-3])
    
    D = h + h*hx*bx + 0.5*h*h*bxx + h*bx*bx
    
    ai = h*h*hx*(0.5*idx) - ithree*h*h*h*(idx*idx)
    bi = D + 2.0*ithree*h*h*h*(idx*idx)
    ci = -1.0*h*h*hx*(0.5*idx) - ithree*h*h*h*(idx*idx)
    
    gb2 = ai*u0[-3] + bi*u0[-2] +ci*u0[-1]
    
    #i=-3
    h = h0[-3]
    hx = 0.5*idx*(h0[-2] - h0[-4])
    bx = 0.5*idx*(b0[-2] - b0[-4])
    bxx = idx*idx*(b0[-2] - 2*b0[-3] + b0[-4])
    
    D = h + h*hx*bx + 0.5*h*h*bxx + h*bx*bx
    
    ai = h*h*hx*(0.5*idx) - ithree*h*h*h*(idx*idx)
    bi = D + 2.0*ithree*h*h*h*(idx*idx)
    ci = -1.0*h*h*hx*(0.5*idx) - ithree*h*h*h*(idx*idx)
    
    gb1 = ai*u0[-4] + bi*u0[-3] +ci*u0[-2]
    
    #i = n
    i = n
    h = h1[0]
    hx = 0.5*idx*(h1[1] - con[i-1][0])
    bx = 0.5*idx*(b1[1] - bed[i-1])
    bxx = idx*idx*(b1[1] - 2*b1[0] + bed[i-1])
    
    D = h + h*hx*bx + 0.5*h*h*bxx + h*bx*bx
    
    ai = h*h*hx*(0.5*idx) - ithree*h*h*h*(idx*idx)
    bi = D + 2.0*ithree*h*h*h*(idx*idx)
    ci = -1.0*h*h*hx*(0.5*idx) - ithree*h*h*h*(idx*idx)
    
    ge1 = ai*u[i-1] + bi*u1[0] + ci*u1[1]
    
    #i = n+1
    h = h1[1]
    hx = 0.5*idx*(h1[2] - h1[0])
    bx = 0.5*idx*(b1[2] - b1[0])
    bxx = idx*idx*(b1[2] - 2*b1[1] + b1[0])
    
    D = h + h*hx*bx + 0.5*h*h*bxx + h*bx*bx
    
    ai = h*h*hx*(0.5*idx) - ithree*h*h*h*(idx*idx)
    bi = D + 2.0*ithree*h*h*h*(idx*idx)
    ci = -1.0*h*h*hx*(0.5*idx) - ithree*h*h*h*(idx*idx)
    
    ge2 = ai*u1[0] + bi*u1[1] +ci*u1[2]
    
    #i = n+2    
    h = h1[2]
    hx = 0.5*idx*(h1[3] - h1[1])
    bx = 0.5*idx*(b1[3] - b1[1])
    bxx = idx*idx*(b1[3] - 2*b1[2] + b1[1])
    
    D = h + h*hx*bx + 0.5*h*h*bxx + h*bx*bx
    
    ai = h*h*hx*(0.5*idx) - ithree*h*h*h*(idx*idx)
    bi = D + 2.0*ithree*h*h*h*(idx*idx)
    ci = -1.0*h*h*hx*(0.5*idx) - ithree*h*h*h*(idx*idx)
    
    ge3 = ai*u1[1] + bi*u1[2] +ci*u1[3]
    
    
    conlb = zeros((3,2))
    conrb = zeros((3,2))
    ubeg = u0[1:4]
    uend = u1[0:3]
    bbeg = b0[1:4]
    bend = b1[0:3]
    
    
    conlb[0][0] = h0[-3]
    conlb[1][0] = h0[-2]
    conlb[2][0] = h0[-1]
    conlb[0][1] = gb1
    conlb[1][1] = gb2
    conlb[2][1] = gb3
    
    conrb[0][0] = h1[0]
    conrb[1][0] = h1[1]
    conrb[2][0] = h1[2]
    conrb[0][1] = ge1
    conrb[1][1] = ge2
    conrb[2][1] = ge3
    
    con = concatenate([conlb, con, conrb])
    bed = concatenate([bbeg,bed,bend])
    u = concatenate([ubeg,u,uend])
        
    #do normal stuff 
        
    #i = 2
    i = 2
    #define the stage
    wi = con[i][0] + bed[i]
    wip1 = con[i+1][0] + bed[i+1]
    wip2 = con[i+2][0] + bed[i+2]
    wip3 = con[i+3][0] + bed[i+3]
    wim1 = con[i-1][0] + bed[i-1]
    wim2 = con[i-2][0] + bed[i-2] 
        
    #reconstruct common values first
        
    #i left and right
        
    #gradients
    dwib = (wi - wim1)
    dwif = (wip1 - wi)
    dwim = 0.5*(wip1 - wim1)
    dhib = (con[i][0] - con[i-1][0])
    dhif = (con[i+1][0] - con[i][0])
    dhim = 0.5*(con[i+1][0] - con[i-1][0])
    dGib = (con[i][1] - con[i-1][1])
    dGif = (con[i+1][1] - con[i][1])
    dGim = 0.5*(con[i+1][1] - con[i-1][1])
    duib = (u[i] - u[i-1])
    duif = (u[i+1] - u[i])
    duim = 0.5*(u[i+1] - u[i-1])
        
    #limiting
    dwi = minmodpy(beta*dwib,beta*dwif,dwim)
    dhi = minmodpy(beta*dhib,beta*dhif,dhim)
    dGi = minmodpy(beta*dGib,beta*dGif,dGim)
    dui = minmodpy(beta*duib,beta*duif,duim)
        
    #reconstruct right
    hir = con[i][0] + 0.5*dhi
    wir = wi + 0.5*dwi
    Gir = con[i][1] + 0.5*dGi
    uir = u[i] + 0.5*dui
    bir = wir - hir
        
    #reconstruct left
    hil = con[i][0] - 0.5*dhi
    wil = wi - 0.5*dwi
    Gil = con[i][1] - 0.5*dGi
    uil = u[i] - 0.5*dui
    bil = wil - hil
        
    #only left of i+1 common but do both
        
    #gradients
    dwip1b = (wip1 - wi)
    dwip1f = (wip2 - wip1)
    dwip1m = 0.5*(wip2 - wi)
    dhip1b = (con[i+1][0] - con[i][0])
    dhip1f = (con[i+2][0] - con[i+1][0])
    dhip1m = 0.5*(con[i+2][0] - con[i][0])
    dGip1b = (con[i+1][1] - con[i][1])
    dGip1f = (con[i+2][1] - con[i+1][1])
    dGip1m = 0.5*(con[i+2][1] - con[i][1])
    duip1b = (u[i+1] - u[i])
    duip1f = (u[i+2] - u[i+1])
    duip1m = 0.5*(u[i+2] - u[i])
        
    #limiting
    dwip1 = minmodpy(beta*dwip1b,beta*dwip1f,dwip1m)
    dhip1 = minmodpy(beta*dhip1b,beta*dhip1f,dhip1m)
    dGip1 = minmodpy(beta*dGip1b,beta*dGip1f,dGip1m)
    duip1 = minmodpy(beta*duip1b,beta*duip1f,duip1m)
        
    #reconstruct right
    hip1r = con[i+1][0] + 0.5*dhip1
    wip1r = wip1 + 0.5*dwip1
    Gip1r = con[i+1][1] + 0.5*dGip1
    uip1r = u[i+1] + 0.5*duip1
    bip1r = wip1r - hip1r
        
    #reconstruct left
    hip1l = con[i+1][0] - 0.5*dhip1
    wip1l = wip1 - 0.5*dwip1
    Gip1l = con[i+1][1] - 0.5*dGip1
    uip1l = u[i+1] - 0.5*duip1
    bip1l = wip1l - hip1l
        
        
    #only right of i-1
    #i-1  right
        
    #gradients
    dwim1b = (wim1 - wim2)
    dwim1f = (wi - wim1)
    dwim1m = 0.5*(wi - wim2)
    dhim1b = (con[i-1][0] - con[i-2][0])
    dhim1f = (con[i][0] - con[i-1][0])
    dhim1m = 0.5*(con[i][0] - con[i-2][0])
    dGim1b = (con[i-1][1] - con[i-2][1])
    dGim1f = (con[i][1] - con[i-1][1])
    dGim1m = 0.5*(con[i][1] - con[i-2][1])
    duim1b = (u[i-1] - u[i-2])
    duim1f = (u[i] - u[i-1])
    duim1m = 0.5*(u[i] - u[i-2])
        
    #limiting
    dwim1 = minmodpy(beta*dwim1b,beta*dwim1f,dwim1m)
    dhim1 = minmodpy(beta*dhim1b,beta*dhim1f,dhim1m)
    dGim1 = minmodpy(beta*dGim1b,beta*dGim1f,dGim1m)
    duim1 = minmodpy(beta*duim1b,beta*duim1f,duim1m)
        
    #reconstruct right
    him1r = con[i-1][0] + 0.5*dhim1
    wim1r = wim1 + 0.5*dwim1
    Gim1r = con[i-1][1] + 0.5*dGim1
    uim1r = u[i-1] + 0.5*duim1
    bim1r = wim1r - him1r
        
    #reconstruct i+2 left
    
    #gradients
    dwip2b = (wip2 - wip1)
    dwip2f = (wip3 - wip2)
    dwip2m = 0.5*(wip3 - wip1)
    dhip2b = (con[i+2][0] - con[i+1][0])
    dhip2f = (con[i+3][0] - con[i+2][0])
    dhip2m = 0.5*(con[i+3][0] - con[i+1][0])
    dGip2b = (con[i+2][1] - con[i+1][1])
    dGip2f = (con[i+3][1] - con[i+2][1])
    dGip2m = 0.5*(con[i+3][1] - con[i+1][1])
    duip2b = (u[i+2] - u[i+1])
    duip2f = (u[i+3] - u[i+2])
    duip2m = 0.5*(u[i+3] - u[i+1])
        
    #limiting
    dwip2 = minmodpy(beta*dwip2b,beta*dwip2f,dwip2m)
    dhip2 = minmodpy(beta*dhip2b,beta*dhip2f,dhip2m)
    dGip2 = minmodpy(beta*dGip2b,beta*dGip2f,dGip2m)
    duip2 = minmodpy(beta*duip2b,beta*duip2f,duip2m)
                
    #reconstruct left
    hip2l = con[i+2][0] - 0.5*dhip2
    wip2l = wip2 - 0.5*dwip2
    Gip2l = con[i+2][1] - 0.5*dGip2
    uip2l = u[i+2] - 0.5*duip2
    bip2l = wip2l - hip2l
    
    #calculate forces              
        
    #right force i
    nbi = max(bip1l,bir)
    hihm = max(0,wir-nbi)
    hihp = max(0,wip1l-nbi)

    her = hihp
    Ger = Gip1l
    uer = uip1l
        
    hel = hihm
    Gel = Gir
    uel = uir
        
    duer = idx*(uip2l - uip1l)
    dber = idx*(bip2l - bip1l)
            
    duel = idx*(uir - uim1r)
    dbel = idx*(bir - bim1r)
        
    sqrtghel = sqrt(g*hel)
    sqrtgher = sqrt(g*her)
    sl = min(0,uel - sqrtghel, uer - sqrtgher)
    sr = max(0,uel + sqrtghel, uer + sqrtgher)
        
    felh = uel*hel
    felG = Gel*uel + 0.5*g*hel*hel - 2*ithree*hel*hel*hel*duel*duel + hel*hel*uel*duel*dbel
    ferh = uer*her
    ferG = Ger*uer + 0.5*g*her*her - 2*ithree*her*her*her*duer*duer + her*her*uer*duer*dber
  
          
    if(sr == 0 and sl ==0):
        foh = 0.0
        foG = 0.0
    else:
        isrmsl = 1.0 / (sr - sl)
        foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(her - hel ))
        foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Ger - Gel ))
    
    fih = foh
    fiG = foG
    himhp = hihp
        
    him1r = hir
    wim1r = wir
    bim1r = bir
    Gim1r = Gir
    uim1r = uir
        
    hil = hip1l
    wil = wip1l
    bil = bip1l
    Gil = Gip1l
    uil = uip1l
        
    hir = hip1r
    wir = wip1r
    bir = bip1r
    Gir = Gip1r
    uir = uip1r
    
    hip1l = hip2l
    wip1l = wip2l
    bip1l = bip2l
    Gip1l = Gip2l
    uip1l = uip2l       
        
    
    dhip1 = dhip2
    dwip1 = dwip2
    duip1 = duip2
    dGip1 = dGip2
    
    for i in range(3,len(con)-3):
        #update both forces at same time

        #define the stage
        wi = con[i][0] + bed[i]
        wip1 = con[i+1][0] + bed[i+1]
        wip2 = con[i+2][0] + bed[i+2]
        wip3 = con[i+3][0] + bed[i+3]
        wim1 = con[i-1][0] + bed[i-1]
        wim2 = con[i-2][0] + bed[i-2] 
        
        #reconstruct common values first
                
        #only left of i+1 common but do both        
        
        #reconstruct right
        hip1r = con[i+1][0] + 0.5*dhip1
        wip1r = wip1 + 0.5*dwip1
        Gip1r = con[i+1][1] + 0.5*dGip1
        uip1r = u[i+1] + 0.5*duip1
        bip1r = wip1r - hip1r
        
        
        #reconstruct i+2 left
        
        #gradients
        dwip2b = (wip2 - wip1)
        dwip2f = (wip3 - wip2)
        dwip2m = 0.5*(wip3 - wip1)
        dhip2b = (con[i+2][0] - con[i+1][0])
        dhip2f = (con[i+3][0] - con[i+2][0])
        dhip2m = 0.5*(con[i+3][0] - con[i+1][0])
        dGip2b = (con[i+2][1] - con[i+1][1])
        dGip2f = (con[i+3][1] - con[i+2][1])
        dGip2m = 0.5*(con[i+3][1] - con[i+1][1])
        duip2b = (u[i+2] - u[i+1])
        duip2f = (u[i+3] - u[i+2])
        duip2m = 0.5*(u[i+3] - u[i+1])
        
        #limiting
        dwip2 = minmod(beta*dwip2b,beta*dwip2f,dwip2m)
        dhip2 = minmod(beta*dhip2b,beta*dhip2f,dhip2m)
        dGip2 = minmod(beta*dGip2b,beta*dGip2f,dGip2m)
        duip2 = minmod(beta*duip2b,beta*duip2f,duip2m)
                
        #reconstruct left
        hip2l = con[i+2][0] - 0.5*dhip2
        wip2l = wip2 - 0.5*dwip2
        Gip2l = con[i+2][1] - 0.5*dGip2
        uip2l = u[i+2] - 0.5*duip2
        bip2l = wip2l - hip2l


                
        #calculate forces              
        
        #right force i
        nbi = max(bip1l,bir)
        hihm = max(0,wir-nbi)
        hihp = max(0,wip1l-nbi)

        her = hihp
        Ger = Gip1l
        uer = uip1l
        
        hel = hihm
        Gel = Gir
        uel = uir
        
        duer = idx*(uip2l - uip1l)
        dber = idx*(bip2l - bip1l)
            
        duel = idx*(uir - uim1r)
        dbel = idx*(bir - bim1r)
        
        sqrtghel = sqrt(g*hel)
        sqrtgher = sqrt(g*her)
        sl = min(0,uel - sqrtghel, uer - sqrtgher)
        sr = max(0,uel + sqrtghel, uer + sqrtgher)
        
        felh = uel*hel
        felG = Gel*uel + 0.5*g*hel*hel - 2*ithree*hel*hel*hel*duel*duel + hel*hel*uel*duel*dbel
        ferh = uer*her
        ferG = Ger*uer + 0.5*g*her*her - 2*ithree*her*her*her*duer*duer + her*her*uer*duer*dber
  
        if(sr == 0 and sl ==0):
            foh = 0.0
            foG = 0.0
        else:
            isrmsl = 1.0 / (sr - sl)
            foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(her - hel ))
            foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Ger - Gel ))
        
        
        #calculate the source term
        th = con[i][0]
        tu = u[i]
        tux = (uil - uir)
        tbx = (bil - bir)
        tbxx = idx*idx*(bed[i+1] - 2*bed[i] + bed[i-1])
        
        sourcer = g*0.5*(hihm*hihm - hir*hir)
        sourcec = g*th*tbx +  0.5*th*th*tu*tux*tbxx - th*tu*tu*tbx*tbxx       
        sourcel = g*0.5*(hil*hil - himhp*himhp)
        
        
        ncon[i-3][0] = con[i][0] - dt*idx*(foh - fih)
        ncon[i-3][1] = con[i][1] - dt*idx*(foG - fiG) + dt*idx*(sourcer+sourcel + sourcec)
        
        fih = foh
        fiG = foG
        himhp = hihp
        
        him1r = hir
        wim1r = wir
        bim1r = bir
        Gim1r = Gir
        uim1r = uir
        
        hil = hip1l
        wil = wip1l
        bil = bip1l
        Gil = Gip1l
        uil = uip1l
        
        hir = hip1r
        wir = wip1r
        bir = bip1r
        Gir = Gip1r
        uir = uip1r
    
        hip1l = hip2l
        wip1l = wip2l
        bip1l = bip2l
        Gip1l = Gip2l
        uip1l = uip2l       
        
    
        dhip1 = dhip2
        dwip1 = dwip2
        duip1 = duip2
        dGip1 = dGip2
        

 
    return ncon   

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
        
#gives exact up to linears, so is second order accurate huzzah    
def getGfromupyc(h,u,bed,u0,u1,h0,h1,b0,b1,dx):
    idx = 1.0 / dx
    ithree = 1.0 / 3.0
        
    n = len(h)

    G = zeros(n)
        
    for i in range(1,n-1):
        th = h[i]
        thx = 0.5*idx*(h[i+1] - h[i-1])
        tbx = 0.5*idx*(bed[i+1] - bed[i-1])
        tbxx = idx*idx*(bed[i+1] -2*bed[i] + bed[i-1])
        
        D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
        
        ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
        bi = D + 2.0*ithree*idx*idx*th*th*th
        ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
        
        G[i] = ai*u[i-1] + bi*u[i] + ci*u[i+1]
        
    #boundary    
    #i=0
    i=0
    th = h[i]
    thx = 0.5*idx*(h[i+1] - h0)
    tbx = 0.5*idx*(bed[i+1] - b0)
    tbxx = idx*idx*(bed[i+1] -2*bed[i] + b0)
            
    D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
            
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
 
    G[i] = ai*u0 + bi*u[i] + ci*u[i+1]
    
    #i = n-1
    i = n-1

    th = h[i]
    thx = 0.5*idx*(h1 - h[i-1])
    tbx = 0.5*idx*(b1 - bed[i-1])
    tbxx = idx*idx*(b1 -2*bed[i] + bed[i-1])
        
    D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
        
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx

    G[i] = ai*u[i-1] + bi*u[i] + ci*u1
            
    return G  
    

def dambreak(x,hf,hc,hl,bot,dx):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    bed = bot*ones(n)
    
    for i in range(n):
        if (x[i] < hc):
            h[i] = hf
        else:
            h[i] = hl
    
    G = getGfromupy(h,u,bed,0.0,0.0,h[0],h[-1],0.0,0.0,dx)
    return h,G,bed

       

def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a

def soliton (x,t,g,a0,a1):
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)
  
def solitoninit(n,a0,a1,g,x,t0,bot,dx):
    h = zeros(n)
    bx = zeros(n)
    u = zeros(n)
    c = sqrt(g*(a0 + a1))
    for i in range(n):
        bx[i] = bot
        h[i] = soliton(x[i],t0,g,a0,a1)
        u[i] =  c* ((h[i] - a0) / h[i])
         
    G = getGfromupy(h,u,bx,0.0,0.0,a0,a0,0.0,0.0,dx)
    
    return h,G,bx 

def powerfunction(r,n):
    if ( r >= 0 and r<= 1):
        return (1-r)**n
    else:
        return 0
    
    
    
def flowoverbump(x,stage,center,width,height,vel,l):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    bed = zeros(n)
    
    for i in range(n): 
        r = abs(x[i] - center) / width
        bed[i] = height*(powerfunction(r,l + 2)*((l*l + 4*l + 3)*r*r*(1.0/3) + (l + 2)*r  + 1))
        h[i] = stage - bed[i]
        u[i] = vel
    
    G = getGfromupyc(h,u,bed,u[0],u[-1],h[0],h[-1],bed[0],bed[-1],dx)    
    
    
    return h,G,bed
    
def flowoverbumpChrispy(x,stage,center,width,height,vel,l):
    n = len(x)
    con = zeros((n,2))
    u = zeros(n)
    bed = zeros(n)
    
    for i in range(n): 
        r = abs(x[i] - center) / width
        bed[i] = height*(powerfunction(r,l + 2)*((l*l + 4*l + 3)*r*r*(1.0/3) + (l + 2)*r  + 1))
        con[i][0] = stage - bed[i]
        u[i] = vel   
    G = getGfromupy(con,u,bed,vel,vel,con[0][0],con[-1][0],bed[0],bed[-1],dx)
    for i in range(n):
        con[i][1] = G[i]   
    
    
    return con,bed
    
### FLOW OVER BUMP ####################
wdir = "../../../data/bumpChris/o2/"

stage = 1.0
center = 1000.0
width = 150
height = 0.5
el = 2.0
vel = 2

g = 9.81
dx = 0.1
Cr = 0.5
l = Cr / (2 + sqrt(g*stage) )
dt = l*dx
theta = 1.0
startx = 000.0
endx = 2000.0 + dx
startt = 0.0
endt = 10*dt  

if not os.path.exists(wdir):
    os.makedirs(wdir)
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)

gap = int(0.1/dt)
    
h,G,bed = flowoverbump(x,stage,center,width,height,vel,el)

hCi = h
GCi = G
    
nBC = 3
nBCs = 4
b0 = bed[0]*ones(nBCs)
b1 = bed[-1]*ones(nBCs)
u0 = vel*ones(nBCs)
u1 = vel*ones(nBCs)   
h0 = h[0]*ones(nBCs)
h1 = h[-1]*ones(nBCs)
    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
bed_c = copyarraytoC(bed)
h0_c  = copyarraytoC(h0)
h1_c  = copyarraytoC(h1)
u0_c  = copyarraytoC(u0)
u1_c  = copyarraytoC(u1)
b0_c  = copyarraytoC(b0)
b1_c  = copyarraytoC(b1)
u_c = mallocPy(n)
    
for i in range(1,len(t)): 
    if (i == 1 or i%gap == 0):
        getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
        u = copyarrayfromC(u_c,n)
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        s = wdir +  "out" + str(i)+".txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
            writefile2.writerow(['dx' ,'dt','time',"cell midpoint" ,'height(m)', 'G' , 'u(m/s)','bed' ])        
                       
            for j in range(n):
                writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j])])
        
    evolvewrap(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,b0_c,b1_c,g,dx,dt,nBC,n,nBCs,theta)
    print (t[i])
        
getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
uC = copyarrayfromC(u_c,n)
GC = copyarrayfromC(G_c,n)
hC = copyarrayfromC(h_c,n)
s = wdir +  "outlast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     writefile2.writerow(['dx' ,'dt','time',"cell midpoint" ,'height(m)', 'G' , 'u(m/s)','bed' ])        
                   
     for j in range(n):
         writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j])])
deallocPy(u_c)   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(h0_c)
deallocPy(h1_c)
deallocPy(u0_c)
deallocPy(u1_c) 


############FLOW OVER BUMP
#set it up so its exact floating point

wdir = "../../../data/bumpChris/o2/"

g = 9.81
dx = 0.1
Cr = 0.5
l = Cr / (2 + sqrt(g*stage) )
dt = l*dx
startx = 0.0
endx = 2000.0 + dx
startt = 0.0
endt = 10*dt  

if not os.path.exists(wdir):
    os.makedirs(wdir)

x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)

beta = 1.0

gapbig = int(0.5/dt)

con,bed = flowoverbumpChrispy(x,stage,center,width,height,vel,el)

hpyi = zeros(n)
Gpyi = zeros(n)
for i in range(n):
    hpyi[i] = con[i][0]
    Gpyi[i] = con[i][1]

x0 = array([x[0] - 4*dx,x[0] - 3*dx,x[0] - 2*dx,x[0] - dx])
x1 = array([x[-1] + dx,x[-1] + 2*dx,x[-1] + 3*dx,x[-1] + 4*dx])

b0 = array([bed[0],bed[0],bed[0],bed[0]])
b1 = array([bed[-1],bed[-1],bed[-1],bed[-1]])
ui = vel
u0 = array([ui,ui,ui,ui])
u1 = array([ui,ui,ui,ui])
    
h0 = array([con[0][0],con[0][0],con[0][0],con[0][0]])
h1 = array([con[-1][0],con[-1][0],con[-1][0],con[-1][0]])

for i in range(1,len(t)):
      
    if(i % gapbig == 0 or i ==1):
        u = getufromGpy(con,bed,u0[-1],u1[0],h0[-1],h1[0],b0[-1],b1[0],dx)
        s = wdir + "saveoutputts" + str(i) + ".txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

            writefile2.writerow(['dx' ,'dt','time', 'height(m)', 'G' , 'u(m/s)','bed' ])        
               
            for j in range(n):
                writefile2.writerow([str(dx),str(dt),str(t[i]), str(con[j][0]) , str(con[j][1]) , str(u[j]),str(bed[j])])  
        file2.close()
             
        
    con = evolvetwopy(con,bed,g,u0,u1,h0,h1,b0,b1,beta,dx,dt)
    print t[i]
    print con[200]

upy = getufromGpy(con,bed,u0[-1],u1[0],h0[-1],h1[0],b0[-1],b1[0],dx)    
s = wdir + "saveoutputtslast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['dx','dt','time', 'height(m)', 'G' , 'u(m/s)','bed' ])        
               
     for j in range(n):
        writefile2.writerow([str(dx),str(dt),str(t[i]), str(con[j][0]) , str(con[j][1]) , str(u[j]),str(bed[j])])  
file2.close()

hpy = zeros(n)
Gpy = zeros(n)
for i in range(n):
    hpy[i] = con[i][0]
    Gpy[i] = con[i][1]


"""
dx = 100.0 / (2**13)
l = 0.01
dt = l*dx
theta = 1.2
startx = 0.0
endx = 1000.0 + dx
startt = 0.0
endt = dt  
    
g = 9.81
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
bot = 0.0
base = 1.0
eta0 = 0.8
x0 = 500

for i in range(5):
    diffuse = (1000) / (10.0**i)
    hh,G,bed = dambreaksmooth(x,x0,base,eta0,diffuse,bot,dx)
    s = "diff = " + str(diffuse)
    plot(x,hh,label=s)
    
h,G,bed = dambreak(x,1.8,x0,1.0,0.0,dx)
s = "actual"
plot(x,h,label=s)
legend()
"""
"""
## smooth DAM BREAK time##########################
from time import time
import os
diffuses = [0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1.0,2.5,5.0,7.5,10.0,25.0,50.0,75.0,100.0,250.0,500.0,750.0,1000.0]
wdirb = "../../data/bigsmooth/o2/"
for ll in range(3,16):
    for k in range(len(diffuses)):
        wdir = wdirb + str(ll) + "/" + str(k) + "/"
        dx = 10.0 / (2**ll)
        l = 0.01
        dt = l*dx
        theta = 1.2
        startx = 0.0
        endx = 1000.0 + dx
        startt = 0.0
        endt = 30.0+(dt*0.9) 
        hf = 1.8
        hl = 1.0
        
        bot = 0.0
        base = hl
        eta0 = hf - hl
        x0 = 500
        diffuse = diffuses[k]
        
        if not os.path.exists(wdir):
            os.makedirs(wdir)  
            
        g = 9.81
            
        x,t = makevar(startx,endx,dx,startt,endt,dt)
        n = len(x)
            
        bot = 0.0
        gap = max(5,int(0.5/dt))
            
        h,G,bed = dambreaksmooth(x,x0,base,eta0,diffuse,bot,dx)
            
        nBC = 3
        nBCs = 4
        b0 = bot*ones(nBCs)
        b1 = bot*ones(nBCs)
        u0 = zeros(nBCs)
        u1 = zeros(nBCs)    
        h0 = hf*ones(nBCs)
        h1 = hl*ones(nBCs)
            
        h_c = copyarraytoC(h)
        G_c = copyarraytoC(G)
        bed_c = copyarraytoC(bed)
        h0_c  = copyarraytoC(h0)
        h1_c  = copyarraytoC(h1)
        u0_c  = copyarraytoC(u0)
        u1_c  = copyarraytoC(u1)
        b0_c  = copyarraytoC(b0)
        b1_c  = copyarraytoC(b1)
        u_c = mallocPy(n)
            
            
        for i in range(1,len(t)):            
            evolvewrap(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,b0_c,b1_c,g,dx,dt,nBC,n,nBCs,theta)
            print (t[i])
                
        getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
        u = copyarrayfromC(u_c,n)
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        s = wdir +  "outlast.txt"
        with open(s,'a') as file2:
             writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
             writefile2.writerow(['dx' ,'dt','time',"cell midpoint" ,'height(m)', 'G' , 'u(m/s)','bed',"diffuse" ])        
                           
             for j in range(n):
                 writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j]),str(diffuse)])
             
        deallocPy(u_c)   
        deallocPy(h_c)
        deallocPy(G_c)
        deallocPy(h0_c)
        deallocPy(h1_c)
        deallocPy(u0_c)
        deallocPy(u1_c) 
"""
"""
## smooth DAM BREAK targetted##########################
difflist = [1,6,8,9,12]

deltaxa = [5,6,7,9,10,11,12,13,14,15]
dxlist = [deltaxa,deltaxa,deltaxa,deltaxa,deltaxa]

diffuses = [0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1.0,2.5,5.0,7.5,10.0,25.0,50.0,75.0,100.0,250.0,500.0,750.0,1000.0]
wdirb = "../../data/bigsmoothtargetted/o2/"
for lk in range(len(difflist)):
    for ll in dxlist[lk]:
        wdir = wdirb + str(ll) + "/" + str(difflist[lk]) + "/"
        dx = ll*(10.0 / (2**11))
        l = 0.01
        dt = l*dx
        theta = 1.2
        startx = 0.0
        endx = 1000.0 + dx
        startt = 0.0
        endt = 30.0+(dt*0.9) 
        hf = 1.8
        hl = 1.0
        
        bot = 0.0
        base = hl
        eta0 = hf - hl
        x0 = 500
        diffuse = diffuses[difflist[lk]]
        
        if not os.path.exists(wdir):
            os.makedirs(wdir)  
            
        g = 9.81
            
        x,t = makevar(startx,endx,dx,startt,endt,dt)
        n = len(x)
            
        bot = 0.0
        gap = max(5,int(0.5/dt))
            
        h,G,bed = dambreaksmooth(x,x0,base,eta0,diffuse,bot,dx)
            
        nBC = 3
        nBCs = 4
        b0 = bot*ones(nBCs)
        b1 = bot*ones(nBCs)
        u0 = zeros(nBCs)
        u1 = zeros(nBCs)    
        h0 = hf*ones(nBCs)
        h1 = hl*ones(nBCs)
            
        h_c = copyarraytoC(h)
        G_c = copyarraytoC(G)
        bed_c = copyarraytoC(bed)
        h0_c  = copyarraytoC(h0)
        h1_c  = copyarraytoC(h1)
        u0_c  = copyarraytoC(u0)
        u1_c  = copyarraytoC(u1)
        b0_c  = copyarraytoC(b0)
        b1_c  = copyarraytoC(b1)
        u_c = mallocPy(n)
            
            
        for i in range(1,len(t)):            
            evolvewrap(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,b0_c,b1_c,g,dx,dt,nBC,n,nBCs,theta)
            print (t[i])
                
        getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
        u = copyarrayfromC(u_c,n)
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        s = wdir +  "outlast.txt"
        with open(s,'a') as file2:
             writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
             writefile2.writerow(['dx' ,'dt','time',"cell midpoint" ,'height(m)', 'G' , 'u(m/s)','bed',"diffuse" ])        
                           
             for j in range(n):
                 writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j]),str(diffuse)])
             
        deallocPy(u_c)   
        deallocPy(h_c)
        deallocPy(G_c)
        deallocPy(h0_c)
        deallocPy(h1_c)
        deallocPy(u0_c)
        deallocPy(u1_c) 
"""
""" 
## DAM BREAK time##########################
from time import time
wdir = "../../data/time/o2/"
samp = 100
dx = 0.001
l = 0.01
dt = l*dx
theta = 1.2
startx = 0.0
endx = 1000.0 + dx
startt = 0.0
endt = samp*dt  
    
g = 9.81
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
    
bot = 0.0
hf = 1.8
hl = 1.0
gap = max(5,int(0.5/dt))
    
h,G,bed = dambreak(x,hf,500,hl,bot,dx)
    
nBC = 3
nBCs = 4
b0 = bot*ones(nBCs)
b1 = bot*ones(nBCs)
u0 = zeros(nBCs)
u1 = zeros(nBCs)    
h0 = hf*ones(nBCs)
h1 = hl*ones(nBCs)
    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
bed_c = copyarraytoC(bed)
h0_c  = copyarraytoC(h0)
h1_c  = copyarraytoC(h1)
u0_c  = copyarraytoC(u0)
u1_c  = copyarraytoC(u1)
b0_c  = copyarraytoC(b0)
b1_c  = copyarraytoC(b1)
u_c = mallocPy(n)
    
    
tim1 =time()     
for i in range(1,len(t)):        
    evolvewrap(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,b0_c,b1_c,g,dx,dt,nBC,n,nBCs,theta)
    print (t[i])
tim2 =time()
tt = tim2 - tim1  
        
s = wdir + "outlast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     writefile2.writerow(['samp' ,'total time(s)','average time(s)'])        
     writefile2.writerow([str(samp), str(tt), str(tt/float(samp))]) 
deallocPy(u_c)   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(h0_c)
deallocPy(h1_c)
deallocPy(u0_c)
deallocPy(u1_c)  
"""
"""
## DAM BREAK##########################
from time import time
wdir = "../../../data/dbChris/o2/"
samp = 100
g = 9.81
hf = 10.0
hl = 1.0
dx = 0.1
Cr = 0.2
l = Cr / sqrt(g*hf) 
dt = l*dx
theta = 1.2
startx = 0.0
endx = 1000.0 + dx
startt = 0.0
endt = 30.0 + dt  

if not os.path.exists(wdir):
    os.makedirs(wdir)
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
    
bot = 0.0
gap = max(5,int(0.1/dt))
    
h,G,bed = dambreak(x,hf,500,hl,bot,dx)
    
nBC = 3
nBCs = 4
b0 = bot*ones(nBCs)
b1 = bot*ones(nBCs)
u0 = zeros(nBCs)
u1 = zeros(nBCs)    
h0 = hf*ones(nBCs)
h1 = hl*ones(nBCs)
    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
bed_c = copyarraytoC(bed)
h0_c  = copyarraytoC(h0)
h1_c  = copyarraytoC(h1)
u0_c  = copyarraytoC(u0)
u1_c  = copyarraytoC(u1)
b0_c  = copyarraytoC(b0)
b1_c  = copyarraytoC(b1)
u_c = mallocPy(n)
    
    
    
for i in range(1,len(t)): 
    if (i == 1 or i%gap == 0):
        getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
        u = copyarrayfromC(u_c,n)
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        s = wdir +  "out" + str(i)+".txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
            writefile2.writerow(['dx' ,'dt','time',"cell midpoint" ,'height(m)', 'G' , 'u(m/s)','bed' ])        
                       
            for j in range(n):
                writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j])])
        
    evolvewrap(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,b0_c,b1_c,g,dx,dt,nBC,n,nBCs,theta)
    print (t[i])
        
getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
u = copyarrayfromC(u_c,n)
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)
s = wdir +  "outlast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     writefile2.writerow(['dx' ,'dt','time',"cell midpoint" ,'height(m)', 'G' , 'u(m/s)','bed' ])        
                   
     for j in range(n):
         writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j])])
deallocPy(u_c)   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(h0_c)
deallocPy(h1_c)
deallocPy(u0_c)
deallocPy(u1_c) 
""" 


"""
##### DAMBREAK ACCURACY
g = 9.81
Cr = 0.2
hf = 1.8
hl = 1.0
wdir = "../../data/dbh/o2/"
for k in range(16,20):
    dx = 100.0 / (2**k)
    #l = Cr / sqrt(g*hf)
    l = 0.01
    dt = l*dx
    theta = 1.2
    startx = 0.0
    endx = 1000.0 + dx
    startt = 0.0
    endt = 30 + dt
    
    wdatadir = wdir+ str(k) + "/"
    if not os.path.exists(wdatadir):
        os.makedirs(wdatadir)
        
    g = 9.81
        
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
        
    bot = 0.0
    gap = max(1,int(0.5/dt))
        
    h,G,bed = dambreak(x,hf,500,hl,bot,dx)
        
    nBC = 3
    nBCs = 4
    b0 = bot*ones(nBCs)
    b1 = bot*ones(nBCs)
    u0 = zeros(nBCs)
    u1 = zeros(nBCs)    
    h0 = hf*ones(nBCs)
    h1 = hl*ones(nBCs)
        
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    bed_c = copyarraytoC(bed)
    h0_c  = copyarraytoC(h0)
    h1_c  = copyarraytoC(h1)
    u0_c  = copyarraytoC(u0)
    u1_c  = copyarraytoC(u1)
    b0_c  = copyarraytoC(b0)
    b1_c  = copyarraytoC(b1)
    u_c = mallocPy(n)
        
        
        
    for i in range(1,len(t)): 
        if (i == 1 or i%gap == 0):
            getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
            u = copyarrayfromC(u_c,n)
            G = copyarrayfromC(G_c,n)
            h = copyarrayfromC(h_c,n)
            s = wdatadir +  "out" + str(i)+".txt"
            with open(s,'a') as file2:
                writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
                writefile2.writerow(['dx' ,'dt','time',"cell midpoint" ,'height(m)', 'G' , 'u(m/s)','bed' ])        
                           
                for j in range(n):
                    writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j])])
            
        evolvewrap(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,b0_c,b1_c,g,dx,dt,nBC,n,nBCs,theta)
        print (t[i])
            
    getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
    u = copyarrayfromC(u_c,n)
    G = copyarrayfromC(G_c,n)
    h = copyarrayfromC(h_c,n)
    s = wdatadir +  "outlast.txt"
    with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
         writefile2.writerow(['dx' ,'dt','time',"cell midpoint" ,'height(m)', 'G' , 'u(m/s)','bed' ])        
                       
         for j in range(n):
             writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j])])
             
    s = wdir +  str(k)+ ".txt"
    with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
         writefile2.writerow(['dx' ,'dt','time',"cell midpoint" ,'height(m)', 'G' , 'u(m/s)','bed' ])        
                   
         for j in range(n):
             writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j])])
    deallocPy(u_c)   
    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(h0_c)
    deallocPy(h1_c)
    deallocPy(u0_c)
    deallocPy(u1_c)
"""

"""
################################# SOLITON Accuracy ####################3
wdir = "../../../data/solcononesec/o2/"

if not os.path.exists(wdir):
    os.makedirs(wdir)

s = wdir + "savenorms.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(['dx','Normalised L1-norm Difference Height', ' Normalised L1-norm Difference Velocity'])
    
for k in range(20):
    dx = 100.0 / (2**k)
    l = 0.01
    dt = l*dx
    startx = -500.0
    endx = 1500.0 + dx
    startt = 0
    endt = 1 + dt
    theta = 1.2
    
    wdatadir = wdir+ str(k) + "/"
    if not os.path.exists(wdatadir):
        os.makedirs(wdatadir)
    
    g = 9.81
    
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
    
    a0 = 10.0
    a1 = 1.0
    t0 = 0
    bot = 0
    gap = max(1,int(0.5/dt))
    
    h,G,bed = solitoninit(n,a0,a1,g,x,t0,bot,dx)
    
    nBC = 3
    nBCs = 4
    b0 = bot*ones(nBCs)
    b1 = bot*ones(nBCs)
    u0 = zeros(nBCs)
    u1 = zeros(nBCs)    
    h0 = h[0]*ones(nBCs)
    h1 = h[-1]*ones(nBCs)
    
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    bed_c = copyarraytoC(bed)
    h0_c  = copyarraytoC(h0)
    h1_c  = copyarraytoC(h1)
    u0_c  = copyarraytoC(u0)
    u1_c  = copyarraytoC(u1)
    b0_c  = copyarraytoC(b0)
    b1_c  = copyarraytoC(b1)
    u_c = mallocPy(n)
    
    
    
    for i in range(1,len(t)):
        
        if(i % gap == 0 or i ==1):
            getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
            u = copyarrayfromC(u_c,n)
            G = copyarrayfromC(G_c,n)
            h = copyarrayfromC(h_c,n)
            
            
            
            
            c = sqrt(g*(a0 + a1))
            htrue = zeros(n)
            utrue = zeros(n)
            for j in range(n):             
                he = soliton(x[j],t[i],g,a0,a1)
                htrue[j] = he
                utrue[j] = c* ((he - a0) / he) 
                
            s = wdatadir + "saveoutputts" + str(i) + ".txt"
            
            print t[i]
            print(h[1],G[1]) 
            with open(s,'a') as file2:
                writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
                writefile2.writerow(['dx' ,'dt','time', 'height(m)', 'G' , 'u(m/s)','bed','true height', 'true velocity' ])        
                   
                for j in range(n):
                    writefile2.writerow([str(dx),str(dt),str(t[i]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j]) , str(htrue[j]), str(utrue[j])])  
                 
            
        evolvewrap(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,b0_c,b1_c,g,dx,dt,nBC,n,nBCs,theta)
        
    getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
    u = copyarrayfromC(u_c,n)
    G = copyarrayfromC(G_c,n)
    h = copyarrayfromC(h_c,n)
    
    
    c = sqrt(g*(a0 + a1))
    htrue = zeros(n)
    utrue = zeros(n)
    for j in range(n):             
        he = soliton(x[j],t[-1],g,a0,a1)
        htrue[j] = he
        utrue[j] = c* ((he - a0) / he) 
    
    s = wdatadir + "saveoutputtslast.txt"
    with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
         writefile2.writerow(['dx' ,'dt','time', 'height(m)', 'G' , 'u(m/s)','bed','true height', 'true velocity'  ])        
                   
         for j in range(n):
             writefile2.writerow([str(dx),str(dt),str(t[-1]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j]), str(htrue[j]), str(utrue[j])])       
    
    normhdiffi = norm(h - htrue,ord=1) / norm(htrue,ord=1)
    normudiffi = norm(u -utrue,ord=1) / norm(utrue,ord=1)  

    s = wdir + "savenorms.txt"
    with open(s,'a') as file1:
        writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

        writefile.writerow([str(dx),str(normhdiffi), str(normudiffi)]) 

    deallocPy(u_c)   
    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(bed_c)
    deallocPy(h0_c)
    deallocPy(h1_c)
    deallocPy(u0_c)
    deallocPy(u1_c)
    deallocPy(b0_c)
    deallocPy(b1_c)    
"""
"""
##### SOLITON
dx = 0.1
l = 0.01
dt = l*dx
startx = -500.0
endx = 1500.0 + dx
startt = 0
endt = 100.0 + dt

theta = 1.2
    
wdir = "../../data/Cserre/soliton/order2/t/"
    
g = 9.81
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
    
a0 = 10.0
a1 = 1.0
t0 = 0
bot = 0
gap = max(1,int(0.5/dt))
    
h,G,bed = solitoninit(n,a0,a1,g,x,t0,bot,dx)
    
nBC = 3
nBCs = 4
b0 = bot*ones(nBCs)
b1 = bot*ones(nBCs)
u0 = zeros(nBCs)
u1 = zeros(nBCs)    
h0 = h[0]*ones(nBCs)
h1 = h[-1]*ones(nBCs)
    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
bed_c = copyarraytoC(bed)
h0_c  = copyarraytoC(h0)
h1_c  = copyarraytoC(h1)
u0_c  = copyarraytoC(u0)
u1_c  = copyarraytoC(u1)
b0_c  = copyarraytoC(b0)
b1_c  = copyarraytoC(b1)
u_c = mallocPy(n)
    
    
    
for i in range(1,len(t)):
       
    if(i % gap == 0 or i ==1):
       getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
       u = copyarrayfromC(u_c,n)
       G = copyarrayfromC(G_c,n)
       h = copyarrayfromC(h_c,n)
 
       c = sqrt(g*(a0 + a1))
       htrue = zeros(n)
       utrue = zeros(n)
       for j in range(n):             
           he = soliton(x[j],t[i],g,a0,a1)
           htrue[j] = he
           utrue[j] = c* ((he - a0) / he) 
                
       s = wdir + "saveoutputts" + str(i) + ".txt"
       print t[i]
       print(h[1],G[1]) 
       with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
            writefile2.writerow(['dx' ,'dt','time', 'height(m)', 'G' , 'u(m/s)','bed','true height', 'true velocity' ])        
                   
            for j in range(n):
                writefile2.writerow([str(dx),str(dt),str(t[i]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j]) , str(htrue[j]), str(utrue[j])])  
                 
            
    evolvewrap(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,b0_c,b1_c,g,dx,dt,nBC,n,nBCs,theta)
        
getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
u = copyarrayfromC(u_c,n)
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)
    
    
c = sqrt(g*(a0 + a1))
htrue = zeros(n)
utrue = zeros(n)
for j in range(n):             
    he = soliton(x[j],t[-1],g,a0,a1)
    htrue[j] = he
    utrue[j] = c* ((he - a0) / he) 
    
s = wdir + "saveoutputtslast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
     writefile2.writerow(['dx' ,'dt','time', 'height(m)', 'G' , 'u(m/s)','bed','true height', 'true velocity'  ])        
                   
     for j in range(n):
         writefile2.writerow([str(dx),str(dt),str(t[-1]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j]), str(htrue[j]), str(utrue[j])])       
    

deallocPy(u_c)   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(bed_c)
deallocPy(h0_c)
deallocPy(h1_c)
deallocPy(u0_c)
deallocPy(u1_c)
deallocPy(b0_c)
deallocPy(b1_c)  
"""
"""
from time import time
import os
from scipy.linalg import norm
dx = 100.0 / (2.0**6)
l = 0.01
dt = l*dx
startx = -500.0
endx = 1500.0 + dx
startt = 0
endt = 100.0 + dt

theta = 1.2
    
wdir = "../../data/timecomp2/o2/"
if not os.path.exists(wdir):
    os.makedirs(wdir) 
    
g = 9.81
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
    
a0 = 10.0
a1 = 1.0
t0 = 0
bot = 0
gap = max(1,int(0.5/dt))
    
h,G,bed = solitoninit(n,a0,a1,g,x,t0,bot,dx)
    
nBC = 3
nBCs = 4
b0 = bot*ones(nBCs)
b1 = bot*ones(nBCs)
u0 = zeros(nBCs)
u1 = zeros(nBCs)    
h0 = h[0]*ones(nBCs)
h1 = h[-1]*ones(nBCs)
    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
bed_c = copyarraytoC(bed)
h0_c  = copyarraytoC(h0)
h1_c  = copyarraytoC(h1)
u0_c  = copyarraytoC(u0)
u1_c  = copyarraytoC(u1)
b0_c  = copyarraytoC(b0)
b1_c  = copyarraytoC(b1)
u_c = mallocPy(n)
    
    
t0 = time()    
for i in range(1,len(t)):
    evolvewrap(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,b0_c,b1_c,g,dx,dt,nBC,n,nBCs,theta)
    print t[i]
t1 = time()
timeelapse = t1 - t0
        
getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
u = copyarrayfromC(u_c,n)
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)
    
    
c = sqrt(g*(a0 + a1))
htrue = zeros(n)
utrue = zeros(n)
for j in range(n):             
    he = soliton(x[j],t[-1],g,a0,a1)
    htrue[j] = he
    utrue[j] = c* ((he - a0) / he) 

normh = norm(h - htrue,ord=1) / norm(htrue,ord=1)    
normu = norm(u - utrue,ord=1) / norm(utrue,ord=1)     
s = wdir + "saveoutputtslast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['dx' ,'dt','time', 'height(m)', 'G' , 'u(m/s)','true height', 'true velocity',"hnorm","unorm", "time taken", "number of steps", "average time per step" ])        
               
     for j in range(n):
         writefile2.writerow([str(dx),str(dt),str(t[i]), str(h[j]) , str(G[j]) , str(u[j]), str(htrue[j]), str(utrue[j]),str(normh),str(normu), str(timeelapse),str(len(t) - 1) , str((1.0*timeelapse)/ (len(t) - 1) )])


deallocPy(u_c)   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(bed_c)
deallocPy(h0_c)
deallocPy(h1_c)
deallocPy(u0_c)
deallocPy(u1_c)
deallocPy(b0_c)
deallocPy(b1_c)  
"""
"""
################ SEGUR AND HAMMACK RECTANGULAR WAVE ###################################
tl = 60.0
b = 0.61
h0 = 0.09
h1 = 0.1
g = 9.81

dx = 0.01
Cr = 0.2
l = Cr / sqrt(g*h1)
theta = 1.2
dt = l*dx
startx = -tl
endx = tl + dx
startt = 0
endt = 50 + dt  

wdir = "../../data/segur2/o2/"
   
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
    
bot = 0.0

gap = 1
    
h,G,bed = experiment1(x,b,h0,h1,dx)

nBC = 3
nBCs = 4
b0 = bot*ones(nBCs)
b1 = bot*ones(nBCs)
u0 = zeros(nBCs)
u1 = zeros(nBCs)    
h0 = h1*ones(nBCs)
h1 = h1*ones(nBCs)
    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
bed_c = copyarraytoC(bed)
h0_c  = copyarraytoC(h0)
h1_c  = copyarraytoC(h1)
u0_c  = copyarraytoC(u0)
u1_c  = copyarraytoC(u1)
b0_c  = copyarraytoC(b0)
b1_c  = copyarraytoC(b1)
u_c = mallocPy(n)
    
    
    
for i in range(1,len(t)): 
    if (i == 1 or i%gap == 0):
        getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
        u = copyarrayfromC(u_c,n)
        G = copyarrayfromC(G_c,n)
        h = copyarrayfromC(h_c,n)
        s = wdir +  "out" + str(i)+".txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
            writefile2.writerow(['dx' ,'dt','time',"cell midpoint" ,'height(m)', 'G' , 'u(m/s)','bed' ])        
                       
            for j in range(n):
                writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j])])
        
    evolvewrap(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,b0_c,b1_c,g,dx,dt,nBC,n,nBCs,theta)
    print (t[i])
        
getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
u = copyarrayfromC(u_c,n)
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)
s = wdir +  "outlast.txt"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     writefile2.writerow(['dx' ,'dt','time',"cell midpoint" ,'height(m)', 'G' , 'u(m/s)','bed' ])        
                   
     for j in range(n):
         writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j])])
deallocPy(u_c)   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(h0_c)
deallocPy(h1_c)
deallocPy(u0_c)
deallocPy(u1_c) 
"""