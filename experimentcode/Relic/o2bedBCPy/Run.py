# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:09:21 2015

@author: jordan
"""
#from Serre2 import *
from scipy import *
import csv
import os
from numpy.linalg import norm  
from matplotlib.pyplot import plot
from scipy.special import ellipj,ellipk,ellipe

from numpy import absolute

i3 = 1.0 / 3.0

    
def copyarraytoC(a):
    n = len(a)
    b = mallocPy(n)
    for i in range(n):
        writetomem(b,i,a[i])
    return b
    
def copywritearraytoC(a,b):
    n = len(a)
    for i in range(n):
        writetomem(b,i,a[i])
    
def copyarrayfromC(a,n):
    b = [0]*n
    for i in range(n):
        b[i] = readfrommem(a,i)
        
    return b
    
def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 

#works for complex
def minmodu(x):
    sx =  10**34
    for xi in x:
        if (abs(xi) < abs(sx)):
            sx = xi
    return sx
    
        
    
def minmod(a,b,c):
    #complex arguments? (only looks at real part for > , <)
    #seperates it if in different half planes
    if(sign(a) == sign(b) and sign(b) == sign(c)):
        return minmodu([a,b,c])
    else:
        return 0
        
def TDMAPY(a,b,c,d):
    n = len(b)
    alpha = zeros(n,dtype=complex)
    beta = zeros(n,dtype=complex)
    x = zeros(n,dtype=complex)

    alpha[0] = c[0] / b[0]
    beta[0] = d[0] / b[0]
    
    for i in range(1,n-1):
        m = 1.0 / (b[i] - a[i-1]*alpha[i-1])
        alpha[i] = c[i]*m
        beta[i] = (d[i] - a[i-1]*beta[i-1])*m

    m = 1.0 / (b[n-1] - a[n-2]*alpha[n-2])
    beta[n-1] = (d[n-1] - a[n-2]*beta[n-2])*m

    x[n-1] = beta[n-1];

    for i in range(n-2,-1,-1):
        x[i] = beta[i] - alpha[i]*x[i+1]
    return x
        
def getufromGperiodicPY(h,G,bed,dx):
    idx = 1.0 / dx;
    n = len(h)
    a = zeros(n-1,dtype=complex)
    b = zeros(n,dtype=complex)
    c = zeros(n-1,dtype=complex)
    ap = zeros(n-1,dtype=complex)
    bp = zeros(n,dtype=complex)
    cp = zeros(n-1,dtype=complex)
    u1 = zeros(n,dtype=complex)
    v1 = zeros(n,dtype=complex)
    y1 = zeros(n,dtype=complex)
    q1 = zeros(n,dtype=complex)
    Gp = zeros(n,dtype=complex)
    ublank = zeros(n,dtype=complex)
    for i in range(1,n-1):
        th = h[i]
        thx = 0.5*idx*(h[i+1] - h[i-1])
        tbx = 0.5*idx*(bed[i+1] - bed[i-1])
        tbxx = idx*idx*(bed[i+1] -2*bed[i]+ bed[i-1])

        D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
        
        ai = -i3*idx*idx*th*th*th + 0.5*idx*th*th*thx
        bi = D + 2.0*i3*idx*idx*th*th*th
        ci = -i3*idx*idx*th*th*th - 0.5*idx*th*th*thx

        a[i-1] = ai
        b[i] = bi
        c[i] = ci


    #Boundary 
    i = 0
    j = n-1
    th = h[i]
    thx = 0.5*idx*(h[i+1] - h[j])
    tbx = 0.5*idx*(bed[i+1] - bed[j])
    tbxx = idx*idx*(bed[i+1] -2*bed[i]+ bed[j])

    D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
	
    ai = -i3*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*i3*idx*idx*th*th*th
    ci = -i3*idx*idx*th*th*th - 0.5*idx*th*th*thx

    b[i] = 2*bi
    c[i] = ci
    b1 = bi
    a1 = ai


    i = n-1
    j = 0
    th = h[i]
    thx = 0.5*idx*(h[j] - h[i-1])
    tbx = 0.5*idx*(bed[j] - bed[i-1])
    tbxx = idx*idx*(bed[j] -2*bed[i]+ bed[i-1])

    D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
	
    ai = -i3*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*i3*idx*idx*th*th*th
    ci = -i3*idx*idx*th*th*th - 0.5*idx*th*th*thx

    a[i-1] = ai
    b[i] = bi + ci*(a1/b1)
    cn = ci

    for i in range(n):
        u1[i] = 0
        v1[i] = 0
        bp[i] = b[i]
        Gp[i] = G[i]
        
        if(i<n-1):
            ap[i] = a[i]
            cp[i] = c[i]

    u1[0] = - b1
    u1[n-1] =  cn
    v1[0] = 1.0
    v1[n-1] = - a1/b1


    y1 = TDMAPY(a,b,c,Gp)

    q1 = TDMAPY(ap,bp,cp,u1)

    dotprodvy = 0
    dotprodvq = 0
    for i in range(n):
        dotprodvy = dotprodvy + v1[i]*y1[i]
        dotprodvq = dotprodvq + v1[i]*q1[i]
        
    for i in range(n):
        ublank[i] = y1[i] - (dotprodvy / (1 + dotprodvq))*q1[i]

    return ublank
        

        
#gives exact up to linears, so is second order accurate huzzah    
def getGfromupy(h,u,bed,u0,u1,h0,h1,b0,b1,dx):
    idx = 1.0 / dx
    ithree = 1.0 / 3.0
        
    n = len(h)
    G = zeros(n,dtype=complex)
        
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

def periodicBCPY(h,G,bed,dx,nBC):
    u = getufromGperiodicPY(h,G,bed,dx)
    
    ubeg = zeros(nBC,dtype=complex)
    hbeg = zeros(nBC,dtype=complex)
    gbeg = zeros(nBC,dtype=complex)
    bbeg = zeros(nBC,dtype=complex)
    uend = zeros(nBC,dtype=complex)
    hend = zeros(nBC,dtype=complex)
    gend = zeros(nBC,dtype=complex)
    bend = zeros(nBC,dtype=complex)

    for i in range(nBC):
        ubeg[i] = u[-nBC + i]
        hbeg[i] = h[-nBC + i]
        gbeg[i] = G[-nBC + i]
        bbeg[i] = bed[-nBC + i]
        
        uend[i] = u[i]
        hend[i] = h[i]
        gend[i] = G[i]
        bend[i] = bed[i]
        
    Gbc = concatenate([gbeg,G,gend])
    hbc = concatenate([hbeg,h,hend])
    bbc = concatenate([bbeg,bed,bend])
    ubc = concatenate([ubeg,u,uend])
    
    return hbc,ubc,Gbc,bbc
     
def evolve(hbc,Gbc,bedbc,ubc,g,theta,dx,dt,nBC,n):
    idx = 1.0 / dx
    
    hblank = zeros(n,dtype=complex)
    Gblank = zeros(n,dtype=complex)
    
    i = nBC - 1
    
    #hil and hir
    dhif = idx*(hbc[i+1] - hbc[i])
    dhib = idx*(hbc[i] - hbc[i-1])
    dhim = 0.5*idx*(hbc[i+1] - hbc[i-1])
    dhilim = minmod(theta*dhif,dhim,theta*dhib)
    hil = hbc[i] - 0.5*dx*dhilim
    hir = hbc[i] + 0.5*dx*dhilim

    #Gil and Gir
    dGif = idx*(Gbc[i+1] - Gbc[i])
    dGib = idx*(Gbc[i] - Gbc[i-1])
    dGim = 0.5*idx*(Gbc[i+1] - Gbc[i-1])
    dGilim = minmod(theta*dGif,dGim,theta*dGib)
    Gil = Gbc[i] - 0.5*dx*dGilim
    Gir = Gbc[i] + 0.5*dx*dGilim

    #wil and wir
    wi = hbc[i] + bedbc[i]
    wip1 = hbc[i+1] + bedbc[i+1]
    wim1 = hbc[i-1] + bedbc[i-1]
    dwif = idx*(wip1 - wi)
    dwib = idx*(wi - wim1)
    dwim = 0.5*idx*(wip1 - wim1)
    dwilim = minmod(theta*dwif,dwim,theta*dwib)
    wil = wi - 0.5*dx*dwilim
    wir = wi + 0.5*dx*dwilim

    #bil and bir
    bil = wil - hil
    bir = wir - hir


    #hip1l and hip1r
    dhip1f = idx*(hbc[i+2] - hbc[i+1])
    dhip1b = idx*(hbc[i+1] - hbc[i])
    dhip1m = 0.5*idx*(hbc[i+2] - hbc[i])
    dhip1lim = minmod(theta*dhip1f,dhip1m,theta*dhip1b)
    hip1l = hbc[i+1] - 0.5*dx*dhip1lim
    hip1r = hbc[i+1] + 0.5*dx*dhip1lim
    
    
    #Gip1l and Gip1r
    dGip1f = idx*(Gbc[i+2] - Gbc[i+1])
    dGip1b = idx*(Gbc[i+1] - Gbc[i])
    dGip1m = 0.5*idx*(Gbc[i+2] - Gbc[i])
    dGip1lim = minmod(theta*dGip1f,dGip1m,theta*dGip1b)
    Gip1l = Gbc[i+1] - 0.5*dx*dGip1lim
    Gip1r = Gbc[i+1] + 0.5*dx*dGip1lim

    #wip1l and wip1
    wip1 = hbc[i+1] + bedbc[i+1]
    wip2 = hbc[i+2] + bedbc[i+2]
    wi = hbc[i] + bedbc[i]
    dwip1f = idx*(wip2 - wip1)
    dwip1b = idx*(wip1 - wi)
    dwip1m = 0.5*idx*(wip2 - wi) 
    dwip1lim = minmod(theta*dwip1f,dwip1m,theta*dwip1b)
    wip1l = wip1 - 0.5*dx*dwip1lim
    wip1r = wip1 + 0.5*dx*dwip1lim

    #bip1l and bip1r
    bip1l = wip1l - hip1l
    bip1r = wip1r - hip1r


    #hip2l and hip2r
    dhip2f = idx*(hbc[i+3] - hbc[i+2])
    dhip2b = idx*(hbc[i+2] - hbc[i+1])
    dhip2m = 0.5*idx*(hbc[i+3] - hbc[i+1])
    dhip2lim = minmod(theta*dhip2f,dhip2m,theta*dhip2b)
    hip2l = hbc[i+2] - 0.5*dx*dhip2lim
    hip2r = hbc[i+2] + 0.5*dx*dhip2lim

    #Gip2l and Gip2r
    dGip2f = idx*(Gbc[i+3] - Gbc[i+2])
    dGip2b = idx*(Gbc[i+2] - Gbc[i+1])
    dGip2m = 0.5*idx*(Gbc[i+3] - Gbc[i+1])
    dGip2lim = minmod(theta*dGip2f,dGip2m,theta*dGip2b)
    Gip2l = Gbc[i+2] - 0.5*dx*dGip2lim
    Gip2r = Gbc[i+2] + 0.5*dx*dGip2lim

    #wip2l and wip2r
    wip2 = hbc[i+2] + bedbc[i+2]
    wip3 = hbc[i+3] + bedbc[i+3]
    wip1 = hbc[i+1] + bedbc[i+1]
    dwip2f = idx*(wip3 - wip2)
    dwip2b = idx*(wip2 - wip1)
    dwip2m = 0.5*idx*(wip3 - wip1)
    dwip2lim = minmod(theta*dwip2f,dwip2m,theta*dwip2b)
    wip2l = wip2 - 0.5*dx*dwip2lim
    wip2r = wip2 + 0.5*dx*dwip2lim

    #bip2l and bip2r
    bip2l = wip2l - hip2l
    bip2r = wip2r - hip2r

    #him1l and him1r
    dhim1f = idx*(hbc[i] - hbc[i-1])
    dhim1b = idx*(hbc[i-1] - hbc[i-2])
    dhim1m = 0.5*idx*(hbc[i] - hbc[i-2])
    dhim1lim = minmod(theta*dhim1f,dhim1m,theta*dhim1b)
    him1l = hbc[i-1] - 0.5*dx*dhim1lim
    him1r = hbc[i-1] + 0.5*dx*dhim1lim

    #Gim1l and Gim1r
    dGim1f = idx*(Gbc[i] - Gbc[i-1])
    dGim1b = idx*(Gbc[i-1] - Gbc[i-2])
    dGim1m = 0.5*idx*(Gbc[i] - Gbc[i-2])
    dGim1lim = minmod(theta*dGim1f,dGim1m,theta*dGim1b)
    Gim1l = Gbc[i-1] - 0.5*dx*dGim1lim
    Gim1r = Gbc[i-1] + 0.5*dx*dGim1lim

    #wim1l and wim1r
    wim1 = hbc[i-1] + bedbc[i-1]
    wi = hbc[i] + bedbc[i]
    wim2 = hbc[i-2] + bedbc[i-2]
    dwim1f = idx*(wi - wim1)
    dwim1b = idx*(wim1 - wim2)
    dwim1m = 0.5*idx*(wi - wim2)
    dwim1lim = minmod(theta*dwim1f,dwim1m,theta*dwim1b)
    wim1l = wim1 - 0.5*dx*dwim1lim
    wim1r = wim1 + 0.5*dx*dwim1lim

    #bim1l and bim1r
    bim1l = wim1l - him1l
    bim1r = wim1r - him1r
    
    nbi = fmax(bip1l, bir)
    hihm = fmax(0, wir - nbi)
    hihp = fmax(0, wip1l - nbi)

    her = hihp
    Ger = Gip1l
    ber = bip1l
    dber = idx*(bip2l - bip1l)
    uer  = 0.5*(ubc[i+1] + ubc[i])
    duer = idx*(ubc[i+1] - ubc[i])

    hel = hihm
    Gel = Gir
    bel = bir
    dbel = idx*(bir - bim1r)
    uel  = 0.5*(ubc[i+1] + ubc[i])
    duel = idx*(ubc[i+1] - ubc[i])

    fhel = uel*hel
    fher = uer*her
	
    fGel = Gel*uel + 0.5*g*hel*hel - 2*i3*hel*hel*hel*duel*duel + hel*hel*uel*duel*dbel
    fGer = Ger*uer + 0.5*g*her*her - 2*i3*her*her*her*duer*duer + her*her*uer*duer*dber

    sqrtghel = sqrt(g* hel)
    sqrtgher = sqrt(g* her)

    sl = fmin(0,fmin(uel - sqrtghel, uer - sqrtgher))
    sr = fmax(0,fmax(uel + sqrtghel, uer + sqrtgher))

    isrmsl = 0.0

    if(sr != sl):
        isrmsl = 1.0 / (sr - sl);

    foh =isrmsl*( sr*fhel - sl*fher + sl*sr*(her - hel))
    foG =isrmsl*( sr*fGel - sl*fGer + sl*sr*(Ger - Gel))

    fih = foh
    fiG = foG
    himhp = hihp
    
    for i in range(nBC , n + nBC):
        
        #hil and hir
        dhif = idx*(hbc[i+1] - hbc[i])
        dhib = idx*(hbc[i] - hbc[i-1])
        dhim = 0.5*idx*(hbc[i+1] - hbc[i-1])
        dhilim = minmod(theta*dhif,dhim,theta*dhib)
        hil = hbc[i] - 0.5*dx*dhilim
        hir = hbc[i] + 0.5*dx*dhilim

        #Gil and Gir
        dGif = idx*(Gbc[i+1] - Gbc[i])
        dGib = idx*(Gbc[i] - Gbc[i-1])
        dGim = 0.5*idx*(Gbc[i+1] - Gbc[i-1])
        dGilim = minmod(theta*dGif,dGim,theta*dGib)
        Gil = Gbc[i] - 0.5*dx*dGilim
        Gir = Gbc[i] + 0.5*dx*dGilim

        #wil and wir
        wi = hbc[i] + bedbc[i]
        wip1 = hbc[i+1] + bedbc[i+1]
        wim1 = hbc[i-1] + bedbc[i-1]
        dwif = idx*(wip1 - wi)
        dwib = idx*(wi - wim1)
        dwim = 0.5*idx*(wip1 - wim1)
        dwilim = minmod(theta*dwif,dwim,theta*dwib)
        wil = wi - 0.5*dx*dwilim
        wir = wi + 0.5*dx*dwilim

        #bil and bir
        bil = wil - hil
        bir = wir - hir


        #hip1l and hip1r
        dhip1f = idx*(hbc[i+2] - hbc[i+1])
        dhip1b = idx*(hbc[i+1] - hbc[i])
        dhip1m = 0.5*idx*(hbc[i+2] - hbc[i])
        dhip1lim = minmod(theta*dhip1f,dhip1m,theta*dhip1b)
        hip1l = hbc[i+1] - 0.5*dx*dhip1lim
        hip1r = hbc[i+1] + 0.5*dx*dhip1lim
        
        
        #Gip1l and Gip1r
        dGip1f = idx*(Gbc[i+2] - Gbc[i+1])
        dGip1b = idx*(Gbc[i+1] - Gbc[i])
        dGip1m = 0.5*idx*(Gbc[i+2] - Gbc[i])
        dGip1lim = minmod(theta*dGip1f,dGip1m,theta*dGip1b)
        Gip1l = Gbc[i+1] - 0.5*dx*dGip1lim
        Gip1r = Gbc[i+1] + 0.5*dx*dGip1lim

        #wip1l and wip1
        wip1 = hbc[i+1] + bedbc[i+1]
        wip2 = hbc[i+2] + bedbc[i+2]
        wi = hbc[i] + bedbc[i]
        dwip1f = idx*(wip2 - wip1)
        dwip1b = idx*(wip1 - wi)
        dwip1m = 0.5*idx*(wip2 - wi) 
        dwip1lim = minmod(theta*dwip1f,dwip1m,theta*dwip1b)
        wip1l = wip1 - 0.5*dx*dwip1lim
        wip1r = wip1 + 0.5*dx*dwip1lim

        #bip1l and bip1r
        bip1l = wip1l - hip1l
        bip1r = wip1r - hip1r


        #hip2l and hip2r
        dhip2f = idx*(hbc[i+3] - hbc[i+2])
        dhip2b = idx*(hbc[i+2] - hbc[i+1])
        dhip2m = 0.5*idx*(hbc[i+3] - hbc[i+1])
        dhip2lim = minmod(theta*dhip2f,dhip2m,theta*dhip2b)
        hip2l = hbc[i+2] - 0.5*dx*dhip2lim
        hip2r = hbc[i+2] + 0.5*dx*dhip2lim

        #Gip2l and Gip2r
        dGip2f = idx*(Gbc[i+3] - Gbc[i+2])
        dGip2b = idx*(Gbc[i+2] - Gbc[i+1])
        dGip2m = 0.5*idx*(Gbc[i+3] - Gbc[i+1])
        dGip2lim = minmod(theta*dGip2f,dGip2m,theta*dGip2b)
        Gip2l = Gbc[i+2] - 0.5*dx*dGip2lim
        Gip2r = Gbc[i+2] + 0.5*dx*dGip2lim

        #wip2l and wip2r
        wip2 = hbc[i+2] + bedbc[i+2]
        wip3 = hbc[i+3] + bedbc[i+3]
        wip1 = hbc[i+1] + bedbc[i+1]
        dwip2f = idx*(wip3 - wip2)
        dwip2b = idx*(wip2 - wip1)
        dwip2m = 0.5*idx*(wip3 - wip1)
        dwip2lim = minmod(theta*dwip2f,dwip2m,theta*dwip2b)
        wip2l = wip2 - 0.5*dx*dwip2lim
        wip2r = wip2 + 0.5*dx*dwip2lim

        #bip2l and bip2r
        bip2l = wip2l - hip2l
        bip2r = wip2r - hip2r

        #him1l and him1r
        dhim1f = idx*(hbc[i] - hbc[i-1])
        dhim1b = idx*(hbc[i-1] - hbc[i-2])
        dhim1m = 0.5*idx*(hbc[i] - hbc[i-2])
        dhim1lim = minmod(theta*dhim1f,dhim1m,theta*dhim1b)
        him1l = hbc[i-1] - 0.5*dx*dhim1lim
        him1r = hbc[i-1] + 0.5*dx*dhim1lim

        #Gim1l and Gim1r
        dGim1f = idx*(Gbc[i] - Gbc[i-1])
        dGim1b = idx*(Gbc[i-1] - Gbc[i-2])
        dGim1m = 0.5*idx*(Gbc[i] - Gbc[i-2])
        dGim1lim = minmod(theta*dGim1f,dGim1m,theta*dGim1b)
        Gim1l = Gbc[i-1] - 0.5*dx*dGim1lim
        Gim1r = Gbc[i-1] + 0.5*dx*dGim1lim

        #wim1l and wim1r
        wim1 = hbc[i-1] + bedbc[i-1]
        wi = hbc[i] + bedbc[i]
        wim2 = hbc[i-2] + bedbc[i-2]
        dwim1f = idx*(wi - wim1)
        dwim1b = idx*(wim1 - wim2)
        dwim1m = 0.5*idx*(wi - wim2)
        dwim1lim = minmod(theta*dwim1f,dwim1m,theta*dwim1b)
        wim1l = wim1 - 0.5*dx*dwim1lim
        wim1r = wim1 + 0.5*dx*dwim1lim

        #bim1l and bim1r
        bim1l = wim1l - him1l
        bim1r = wim1r - him1r
        
        nbi = fmax(bip1l, bir)
        hihm = fmax(0, wir - nbi)
        hihp = fmax(0, wip1l - nbi)

        her = hihp
        Ger = Gip1l
        ber = bip1l
        dber = idx*(bip2l - bip1l)
        uer  = 0.5*(ubc[i+1] + ubc[i])
        duer = idx*(ubc[i+1] - ubc[i])

        hel = hihm
        Gel = Gir
        bel = bir
        dbel = idx*(bir - bim1r)
        uel  = 0.5*(ubc[i+1] + ubc[i])
        duel = idx*(ubc[i+1] - ubc[i])

        fhel = uel*hel
        fher = uer*her
	
        fGel = Gel*uel + 0.5*g*hel*hel - 2*i3*hel*hel*hel*duel*duel + hel*hel*uel*duel*dbel
        fGer = Ger*uer + 0.5*g*her*her - 2*i3*her*her*her*duer*duer + her*her*uer*duer*dber

        sqrtghel = sqrt(g* hel)
        sqrtgher = sqrt(g* her)

        sl = fmin(0,fmin(uel - sqrtghel, uer - sqrtgher))
        sr = fmax(0,fmax(uel + sqrtghel, uer + sqrtgher))

        isrmsl = 0.0

        if(sr != sl):
            isrmsl = 1.0 / (sr - sl);

        foh =isrmsl*( sr*fhel - sl*fher + sl*sr*(her - hel))
        foG =isrmsl*( sr*fGel - sl*fGer + sl*sr*(Ger - Gel))

        th = hbc[i]
        tu = ubc[i]
        tux = -0.5*(ubc[i+1] - ubc[i-1])
        tbx = (bil - bir)
        tbxx = idx*idx*(bedbc[i+1] - 2*bedbc[i] + bedbc[i-1])
		
        sourcer = g*0.5*(hihm*hihm - hir*hir)
        sourcec = g*th*tbx +  0.5*th*th*tu*tux*tbxx - th*tu*tu*tbx*tbxx 
        sourcel = g*0.5*(hil*hil - himhp*himhp)

        hblank[i - nBCn] = hbc[i] - dt*idx*(foh - fih)
        Gblank[i - nBCn] = Gbc[i] - dt*idx*(foG -fiG) + dt*idx*(sourcer+sourcel + sourcec)

        fih = foh
        fiG = foG
        himhp = hihp
        
    return hblank,Gblank
        
def evolvewrap(G,h,bed,g,dx,dt,theta):
    n = len(h)
    nBC = 3
    hbc,ubc,Gbc,bbc = periodicBCPY(h,G,bed,dx,nBC)
    hp,Gp =  evolve(hbc,Gbc,bbc,ubc,g,theta,dx,dt,nBC,n)
    hpbc,upbc,Gpbc,bpbc = periodicBCPY(hp,Gp,bed,dx,nBC)
    hpp,Gpp =  evolve(hpbc,Gpbc,bpbc,upbc,g,theta,dx,dt,nBC,n)
    
    return 0.5*(h + hpp) , 0.5*(G + Gpp)
    
    
    
         

def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a

def soliton (x,t,g,a0,a1):
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)
  
def solitoninit(n,a0,a1,g,x,t0,bot,dx):
    h = zeros(n,dtype=complex)
    bx = zeros(n,dtype=complex)
    u = zeros(n,dtype=complex)
    c = sqrt(g*(a0 + a1))
    for i in range(n):
        bx[i] = bot
        h[i] = soliton(x[i],t0,g,a0,a1)
        u[i] =  c* ((h[i] - a0) / h[i])
         
    G = getGfromupy(h,u,bx,0.0,0.0,a0,a0,0.0,0.0,dx)
    
    return h,u,G,bx 

def dnsq(eta,m):
    
    sn,cn,dn,ph = ellipj(eta,m)
    
    return dn*dn
    

def cnoidalwaves(x,t,dx,a0,a1,g,k):
    
    n = len(x)
    u = zeros(n,dtype=complex)
    h = zeros(n,dtype=complex)
    bed = zeros(n,dtype=complex)
    
    m = k*k
    
    h0 = a0 + a1*(float(ellipe(m)) / ellipk(m))    
    
    c = sqrt((g*a0*(a0 + a1)*(a0 + (1 - k*k)*a1))) / float(h0)
    
    Kc = sqrt(float(3*a1) / (4*a0*(a0 + a1)*(a0 + (1-k*k)*a1)))

    
    for i in range(n):
        h[i] = a0 + a1*dnsq(Kc*(x[i] - c*t),m)
        u[i] = c *(1 - float(h0)/h[i])
        
    b0i = 0
    b1i = 0
    
    cx = x[0] - dx
    
    h0i = a0 + a1*dnsq(Kc*(cx - c*t),m)
    u0i = c *(1 - float(h0)/h0i)
    
    cx = x[-1] + dx
    
    h1i = a0 + a1*dnsq(Kc*(cx - c*t),m)
    u1i = c *(1 - float(h0)/h1i)
    
    G = getGfromupy(h,u,bed,u0i,u1i,h0i,h1i,b0i,b1i,dx)   
    
    return h,u,G,bed
    
def eigenvectors(x,t,g,k,h0,u0,h1,u1):
    n = len(x)
    w = u0*k + k*sqrt(g*h0)*sqrt(3.0/ (h0*h0*k*k + 3))
    imu = sqrt(-1)
    bed = zeros(n,dtype=complex)
    h = zeros(n,dtype=complex)
    u = zeros(n,dtype=complex)
    for i in range(n):
        h[i] = h0 + h1*exp(imu*(w*t + k*x[i]))
        u[i] = u0 + u1*exp(imu*(w*t + k*x[i]))
    
    cx = x[0] - dx     
    h0i = h0 + h1*exp(imu*(w*t + k*cx))
    u0i = u0 + u1*exp(imu*(w*t + k*cx))
    
    cx = x[-1] + dx 
    h1i = h0 + h1*exp(imu*(w*t + k*cx))
    u1i = h0 + h1*exp(imu*(w*t + k*cx))
    
    G = getGfromupy(h,u,bed,u0i,u1i,h0i,h1i,bed[0],bed[-1],dx) 
    
    return h,u,G,bed
    
 
####eigenvectors

wdatadir = "../../../data/raw/eigenvectors/o2/"

if not os.path.exists(wdatadir):
    os.makedirs(wdatadir)

h0v = 1.0
u0v = 1.0
h1v = 0.001
u1v = 0.001
k = 0.1

lamb = (2*pi)/k

g = 9.81
dx = (lamb) / (3**5)
Cr = 0.5
l = Cr / (sqrt(g*(h0v + h1v) ))
dt = l*dx
theta = 2
startx = dx
endx = 3*lamb + 0.5*dx
startt = 0.0
endt = 2*dt  

wdir = wdatadir +  "/"

if not os.path.exists(wdir):
    os.makedirs(wdir)


nBCn = 3
nBC = 3
    
xbc,t = makevar(startx - nBCn*dx,endx + nBCn*dx,dx,startt,endt,dt)

x = xbc[nBCn: -nBCn]
xbeg = x[:nBCn]
xend = x[-nBCn:] 

n = len(x)
m = len(t)

gap = int(10.0/dt)

t0 = 0.0
    

#initial conditions for time steps
tij = 0.0
hBC,uBC,GBC,bedBC = eigenvectors(xbc,0,g,k,h0v,u0v,h1v,u1v)
h = hBC[nBCn:-nBCn]
u = uBC[nBCn:-nBCn]
G = GBC[nBCn:-nBCn]
bed = bedBC[nBCn:-nBCn]



#ublank  = getufromGperiodicPY(h,G,bed,dx)

#hbc,ubc,Gbc,bbc = periodicBCPY(h,G,bed,dx,nBC)

for i in range(1,len(t)): 
    h, G = evolvewrap(G,h,bed,g,dx,dt,theta)
    print(t[i])

u = getufromGperiodicPY(h,G,bed,dx)

ha,ua,Ga,beda = eigenvectors(x,t[-1],g,k,h0v,u0v,h1v,u1v)



"""
## Test accuracy
wdatadir = "../../../data/raw/cnoidaltestPY/o2/"
if not os.path.exists(wdatadir):
    os.makedirs(wdatadir)

s = wdatadir + "savenorms.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

for ij in range(7):

    if not os.path.exists(wdatadir):
        os.makedirs(wdatadir)
    a0 = 1.0
    a1 = 0.1
    k = 0.99
    
    ### WAVE LENGTH
        
    m = k*k
    Kc = sqrt(float(3*a1) / (4*a0*(a0 + a1)*(a0 + (1-m)*a1)))
    
    lamb = 2*ellipk(m) / Kc
    
    g = 9.81
    dx = (lamb) / (3**ij)
    Cr = 0.5
    l = Cr / (sqrt(g*(a0 +a1) ))
    dt = l*dx
    theta = 1.5
    startx = dx
    endx = 3*lamb + 0.5*dx
    startt = 0.0
    endt = 1 + dt  
    
    wdir = wdatadir + str(ij) +  "/"
    
    if not os.path.exists(wdir):
        os.makedirs(wdir)
    
    
    nBCn = 3
    nBC = 3
        
    xbc,t = makevar(startx - nBCn*dx,endx + nBCn*dx,dx,startt,endt,dt)
    
    x = xbc[nBCn: -nBCn]
    xbeg = x[:nBCn]
    xend = x[-nBCn:] 
    
    n = len(x)
    m = len(t)
    
    gap = int(10.0/dt)
    
    t0 = 0.0
        
    
    #initial conditions for time steps
    tij = 0.0
    hBC,uBC,GBC,bedBC = cnoidalwaves(xbc,tij,dx,a0,a1,g,k)
    h = hBC[nBCn:-nBCn]
    u = uBC[nBCn:-nBCn]
    G = GBC[nBCn:-nBCn]
    bed = bedBC[nBCn:-nBCn]
    
    
    
    #ublank  = getufromGperiodicPY(h,G,bed,dx)
    
    #hbc,ubc,Gbc,bbc = periodicBCPY(h,G,bed,dx,nBC)
    
    for i in range(1,len(t)): 
        h, G = evolvewrap(G,h,bed,g,dx,dt,theta)
        print(t[i])
    
    u = getufromGperiodicPY(h,G,bed,dx)
    
    ha,ua,Ga,beda = cnoidalwaves(x,t[-1],dx,a0,a1,g,k)
    
    normhdiffi = norm(h - ha,ord=1) / norm(ha,ord=1)
    normudiffi = norm(u -ua,ord=1) / norm(ua,ord=1) 
    
    s = wdatadir + "savenorms.txt"
    with open(s,'a') as file1:
        writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile.writerow([str(dx),str(theta),str(normhdiffi), str(normudiffi)])   
"""

"""
## Single
wdatadir = "../../../data/raw/cnoidaltest/o2/"

if not os.path.exists(wdatadir):
    os.makedirs(wdatadir)
a0 = 1.0
a1 = 0.1
k = 0.99

### WAVE LENGTH
    
m = k*k
Kc = sqrt(float(3*a1) / (4*a0*(a0 + a1)*(a0 + (1-m)*a1)))

lamb = 2*ellipk(m) / Kc

g = 9.81
dx = (lamb) / (3**5)
Cr = 0.5
l = Cr / (sqrt(g*(a0 +a1) ))
dt = l*dx
theta = 2
startx = dx
endx = 3*lamb + 0.5*dx
startt = 0.0
endt = 1 + dt  

wdir = wdatadir +  "/"

if not os.path.exists(wdir):
    os.makedirs(wdir)


nBCn = 3
nBC = 3
    
xbc,t = makevar(startx - nBCn*dx,endx + nBCn*dx,dx,startt,endt,dt)

x = xbc[nBCn: -nBCn]
xbeg = x[:nBCn]
xend = x[-nBCn:] 

n = len(x)
m = len(t)

gap = int(10.0/dt)

t0 = 0.0
    

#initial conditions for time steps
tij = 0.0
hBC,uBC,GBC,bedBC = cnoidalwaves(xbc,tij,dx,a0,a1,g,k)
h = hBC[nBCn:-nBCn]
u = uBC[nBCn:-nBCn]
G = GBC[nBCn:-nBCn]
bed = bedBC[nBCn:-nBCn]



#ublank  = getufromGperiodicPY(h,G,bed,dx)

#hbc,ubc,Gbc,bbc = periodicBCPY(h,G,bed,dx,nBC)

for i in range(1,len(t)): 
    h, G = evolvewrap(G,h,bed,g,dx,dt,theta)
    print(t[i])

u = getufromGperiodicPY(h,G,bed,dx)

ha,ua,Ga,beda = cnoidalwaves(x,t[-1],dx,a0,a1,g,k)
"""
   

    
####cnoidal waves periodic
"""
wdatadir = "../../../data/raw/cnoidaltest/o2/"

if not os.path.exists(wdatadir):
    os.makedirs(wdatadir)
    
s = wdatadir + "savenorms.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(["dx",'theta','l1h', 'l1u'])  
for ij in range(1,9):
    a0 = 1.0
    a1 = 0.1
    k = 0.99
    
    ### WAVE LENGTH
        
    m = k*k
    Kc = sqrt(float(3*a1) / (4*a0*(a0 + a1)*(a0 + (1-m)*a1)))
    
    lamb = 2*ellipk(m) / Kc
    
    g = 9.81
    dx = (lamb) / (3**ij)
    Cr = 0.5
    l = Cr / (sqrt(g*(a0 +a1) ))
    dt = l*dx
    theta = 2
    startx = dx
    endx = 9*lamb + 0.5*dx
    startt = 0.0
    endt = 1 + dt  
    
    wdir = wdatadir + str(ij) + "/"
    
    if not os.path.exists(wdir):
        os.makedirs(wdir)
    
    
    nBCn = 3
    nBC = 3
        
    xbc,t = makevar(startx - nBCn*dx,endx + nBCn*dx,dx,startt,endt,dt)
    
    x = xbc[nBCn: -nBCn]
    xbeg = x[:nBCn]
    xend = x[-nBCn:] 
    
    n = len(x)
    m = len(t)
    
    gap = int(10.0/dt)
    
    t0 = 0.0
        
    
    #initial conditions for time steps
    tij = 0.0
    hBC,uBC,GBC,bedBC = cnoidalwaves(xbc,tij,dx,a0,a1,g,k)
    h = hBC[nBCn:-nBCn]
    u = uBC[nBCn:-nBCn]
    G = GBC[nBCn:-nBCn]
    bed = bedBC[nBCn:-nBCn]
       
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    bed_c = copyarraytoC(bed)
    x_c = copyarraytoC(x)
    u_c = mallocPy(n)
    
    un_c = mallocPy(n+2*nBCn)
    Gn_c = mallocPy(n+2*nBCn)
    hn_c = mallocPy(n+2*nBCn)
    
    hi,ui,Gi,bedi = cnoidalwaves(x,tij,dx,a0,a1,g,k)
    
        
    for i in range(1,len(t)):  
        evolvewrapperiodic(G_c,h_c,bed_c,g,dx,dt,n,nBCn,theta,hn_c, Gn_c,un_c);    
        print (t[i])
            
        tij = t[i]
    
    getufromGperiodic(h_c,G_c,bed_c, dx ,n,u_c)
    #getufromG(h_c,G_c,bed_c,ubeg[-1],uend[0],hbeg[-1],hend[0], 0.0, 0.0, dx ,n,u_c)
    u = copyarrayfromC(u_c,n)
    G = copyarrayfromC(G_c,n)
    h = copyarrayfromC(h_c,n)
    
    un = copyarrayfromC(un_c,n+2*nBCn)
    Gn = copyarrayfromC(Gn_c,n+2*nBCn)
    hn = copyarrayfromC(hn_c,n+2*nBCn)
    
    ha,ua,Ga,beda = cnoidalwaves(x,t[-1],dx,a0,a1,g,k)  
    
    s = wdir + "outlast.txt"
    with open(s,'a') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow(['dx' ,'dt','time' ,"cell midpoint", 'height(m)', 'G' , 'u(m/s)' ,'ha', 'u'])        
               
        for j in range(n):
            writefile2.writerow([str(dx),str(dt),str(t[-1]),str(x[j]) ,str(h[j]) , str(G[j]) , str(u[j]), str(ha[j]), str(ua[j])])     
    
    normhdiffi = norm(h - ha,ord=1) / norm(ha,ord=1)
    normudiffi = norm(u -ua,ord=1) / norm(ua,ord=1) 
    
    s = wdatadir + "savenorms.txt"
    with open(s,'a') as file1:
        writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile.writerow([str(dx),str(theta),str(normhdiffi), str(normudiffi)])    
        
    deallocPy(u_c)   
    deallocPy(h_c)
    deallocPy(G_c)
"""    



#### soliton

"""
wdatadir = "../../../data/raw/soltest/o2/"

if not os.path.exists(wdatadir):
    os.makedirs(wdatadir)
    
s = wdatadir + "savenorms.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(["dx",'theta','l1h', 'l1u'])  

if not os.path.exists(wdatadir):
    os.makedirs(wdatadir)

for ij in range(1,9):
    a0 = 1.0
    a1 = 1.0
    k = 0.9999999
    
    ### WAVE LENGTH
        
    m = k*k
    Kc = sqrt(float(3*a1) / (4*a0*(a0 + a1)*(a0 + (1-m)*a1)))
    
    lamb = 2*ellipk(m) / Kc
    
    g = 9.81
    dx = (lamb) / (3**ij)
    Cr = 0.5
    l = Cr / (sqrt(g*(a0 +a1) ))
    dt = l*dx
    theta = 2
    startx = -20
    endx = 20
    startt = 0.0
    endt = 1 + dt  
    
    wdir = wdatadir + str(ij) + "/"
    
    if not os.path.exists(wdir):
        os.makedirs(wdir)
    
    
    nBCn = 3
    nBC = 3
        
    xbc,t = makevar(startx - nBCn*dx,endx + nBCn*dx,dx,startt,endt,dt)
    
    x = xbc[nBCn: -nBCn]
    xbeg = x[:nBCn]
    xend = x[-nBCn:] 
    
    n = len(x)
    m = len(t)
    
    gap = int(10.0/dt)
    
    t0 = 0.0
        
    
    #initial conditions for time steps
    tij = 0.0
    hBC,uBC,GBC,bedBC = solitoninit(n + 2*nBCn,a0,a1,g,xbc,tij,0,dx)
    hbeg = hBC[:nBCn]
    hend = hBC[-nBCn:]
    h = hBC[nBCn:-nBCn]
    ubeg = uBC[:nBCn]
    uend = uBC[-nBCn:]
    u = uBC[nBCn:-nBCn]
    Gbeg = GBC[:nBCn]
    Gend = GBC[-nBCn:]
    G = GBC[nBCn:-nBCn]
    bedbeg = bedBC[:nBCn]
    bedend = bedBC[-nBCn:]
    bed = bedBC[nBCn:-nBCn]
       
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    bed_c = copyarraytoC(bed)
    x_c = copyarraytoC(x)
    u_c = mallocPy(n)
    
    un_c = mallocPy(n+2*nBCn)
    Gn_c = mallocPy(n+2*nBCn)
    hn_c = mallocPy(n+2*nBCn)
    
    b0_c = copyarraytoC(zeros(nBCn))
    b1_c = copyarraytoC(zeros(nBCn)) 
    #Gbcn = concatenate([Gbeg[:nBCn],G,Gend[-nBCn:]])
    
    hi,ui,Gi,bedi = solitoninit(n,a0,a1,g,x,tij,0,dx)
        
    for i in range(1,len(t)): 
        
        
        #getufromG(h_c,G_c,bed_c,ubeg[-1],uend[0],hbeg[-1],hend[0], 0.0, 0.0, dx ,n,u_c)
        #u = copyarrayfromC(u_c,n)
        #G = copyarrayfromC(G_c,n)
        #h = copyarrayfromC(h_c,n) 
          
        
        
        
        #evolvewrap(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,G0_c,G1_c,h0h_c,h1h_c,u0h_c,u1h_c,G0h_c,G1h_c,b0_c,b1_c,g,dx,dt,n,nBC,nBCn,theta,hn_c, Gn_c,un_c)
        evolvewrapperiodic(G_c,h_c,bed_c,g,dx,dt,n,nBCn,theta,hn_c, Gn_c,un_c);    
        print (t[i])
            
        tij = t[i]
    
        
        
        
    
    getufromGperiodic(h_c,G_c,bed_c, dx ,n,u_c)
    #getufromG(h_c,G_c,bed_c,ubeg[-1],uend[0],hbeg[-1],hend[0], 0.0, 0.0, dx ,n,u_c)
    u = copyarrayfromC(u_c,n)
    G = copyarrayfromC(G_c,n)
    h = copyarrayfromC(h_c,n)
    
    un = copyarrayfromC(un_c,n+2*nBCn)
    Gn = copyarrayfromC(Gn_c,n+2*nBCn)
    hn = copyarrayfromC(hn_c,n+2*nBCn)
    
    ha,ua,Ga,beda = solitoninit(n,a0,a1,g,x,t[-1],0,dx)
    
    haBC,uaBC,GaBC,bedBC = solitoninit(n + 2*nBCn,a0,a1,g,xbc,0,0,dx)
    
    s = wdir + "outlast.txt"
    with open(s,'a') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow(['dx' ,'dt','time' ,"cell midpoint", 'height(m)', 'G' , 'u(m/s)' ,'ha', 'u'])        
               
        for j in range(n):
            writefile2.writerow([str(dx),str(dt),str(t[-1]),str(x[j]) ,str(h[j]) , str(G[j]) , str(u[j]), str(ha[j]), str(ua[j])])     
    
    normhdiffi = norm(h - ha,ord=1) / norm(ha,ord=1)
    normudiffi = norm(u -ua,ord=1) / norm(ua,ord=1) 
    
    s = wdatadir + "savenorms.txt"
    with open(s,'a') as file1:
        writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile.writerow([str(dx),str(theta),str(normhdiffi), str(normudiffi)]) 
        
    deallocPy(u_c)   
    deallocPy(h_c)
    deallocPy(G_c)

"""



"""
###Soliton over bump scenario
### Energy Test ####################

wdatadir = "../../../data/raw/solbumpEnergucontfullfix0m3/o2/"

if not os.path.exists(wdatadir):
    os.makedirs(wdatadir)


s = wdatadir +  "savenorms.txt"
with open(s,'a') as file3:
    writefile3 = csv.writer(file3, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile3.writerow(['dx' ,'rel Eval'])



for k in range(6,19,3):
    stage = 1
    center = 50.0
    width = 25
    height = 0.5
    el = 4.0
    vel = 0
    
    
    a0 = 1.0
    a1 = 1.0
    
    g = 9.81
    dx = 100.0 / (2**k)
    Cr = 0.5
    l = Cr / (2.5 + sqrt(g*(a0 +a1) ))
    dt = l*dx
    theta = 1.0
    startx = -200
    endx = 250.0 + dx
    startt = 0.0
    endt = 50 + dt  
    
    wdir = wdatadir + str(k) + "/"
    
    if not os.path.exists(wdir):
        os.makedirs(wdir)
        
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
    m = len(t)
    
    gap = int(10.0/dt)
        
    #h,G,bed = flowoverbump(x,stage,center,width,height,vel,el)
    h,G,bed = soloverbump(x,a0,a1,-20,20,0.0,g,stage,center,width,height,vel,el)
        
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
    
    xbeg = arange(startx - nBC*dx,startx,dx)
    xend = arange(endx + dx,endx + (nBC+1)*dx) 
    
    xbc =  concatenate([xbeg,x,xend]) 
    
    xbc_c = copyarraytoC(xbc)
    hbc_c = mallocPy(n + 2*nBC)
    ubc_c = mallocPy(n + 2*nBC)
    bedbc_c = mallocPy(n + 2*nBC)
    conc(b0_c , bed_c,b1_c,nBC,n ,nBC , bedbc_c)
    #HEvals = []
    GNEvals = []
    Etimes = []
    
    tBC = 2
    #initial conditions for time steps
        
    for i in range(1,len(t)): 
    
        if (i == 1 or i%gap == 0):
            ki = i
            getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
            u = copyarrayfromC(u_c,n)
            G = copyarrayfromC(G_c,n)
            h = copyarrayfromC(h_c,n)
            
            conc(h0_c , h_c,h1_c,nBC,n ,nBC , hbc_c)
            conc(u0_c , u_c,u1_c,nBC,n ,nBC , ubc_c)  
            #HEval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*nBC,nBC,dx)
            GNEval = GNall(xbc_c,hbc_c,ubc_c,bedbc_c,g,n + 2*nBC,nBC,dx)
            
            
            
            #HEvals.append(HEval)
            GNEvals.append(GNEval)
            Etimes.append(t[i-1])
            
            s = wdir +  "out" + str(i)+".txt"
            with open(s,'a') as file2:
                writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
                writefile2.writerow(['dx' ,'dt','time','Eval',"cell midpoint" ,'height(m)', 'G' , 'u(m/s)','bed' ])        
                           
                for j in range(n):
                    writefile2.writerow([str(dx),str(dt),str(t[i]), str(GNEval),str(x[j]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j])])
            
        evolvewrap(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,b0_c,b1_c,g,dx,dt,n,nBCs,theta)
        print (t[i])
            
    getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], b0[-1], b1[0], dx ,n,u_c)
    u = copyarrayfromC(u_c,n)
    G = copyarrayfromC(G_c,n)
    h = copyarrayfromC(h_c,n)
    
    conc(h0_c , h_c,h1_c,nBC,n ,nBC , hbc_c)
    conc(u0_c , u_c,u1_c,nBC,n ,nBC , ubc_c)        
    #HEval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*nBC,nBC,dx)
    GNEval = GNall(xbc_c,hbc_c,ubc_c,bedbc_c,g,n + 2*nBC,nBC,dx)
    
    
    
    #HEvals.append(HEval)
    GNEvals.append(GNEval)
    Etimes.append(t[-1])
    s = wdir +  "outlast.txt"
    with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
         writefile2.writerow(['dx' ,'dt','time','Eval',"cell midpoint" ,'height(m)', 'G' , 'u(m/s)','bed' ])        
                       
         for j in range(n):
             writefile2.writerow([str(dx),str(dt),str(t[i]),str(GNEval),str(x[j]), str(h[j]) , str(G[j]) , str(u[j]),str(bed[j])])
    
    finE = GNEvals[-1]
    tfinE = Etimes[-1]
    begE = GNEvals[0]
    tbegE = Etimes[0]
    relerr = abs(finE - begE) / abs(begE)
    
    print(str(tbegE) + " || " + str(begE))
    print(str(tfinE) + " || " + str(finE))
    s = wdatadir +  "savenorms.txt"
    with open(s,'a') as file3:
        writefile3 = csv.writer(file3, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile3.writerow([str(dx) ,str(relerr)])
    
    deallocPy(u_c)   
    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(h0_c)
    deallocPy(h1_c)
    deallocPy(u0_c)
    deallocPy(u1_c) 
"""