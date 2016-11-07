
import csv
from numpy.linalg import norm
from scipy import *
import os
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
order = "grim"
wdir = "../../../../data/raw/NEWdata/FDredo/"+order+"/"
sdir = "../../../../data/postprocessing/scFDallA/"+order+"/"

from numpy import tanh, arctanh   
def Soliton(d):
    AB = -22.6552*arctanh(0.707107*tanh(0.612372*(-50 - d))) + 32.0393*tanh(0.612372*(-50 - d))
    AE = -22.6552*arctanh(0.707107*tanh(0.612372*(250 + d))) + 32.0393*tanh(0.612372*(250 + d))
    
    BB = 9.81*(-50 - d) + (42.7191 + 5.33989*sech2(0.612372*(-50 - d)))*tanh(0.612372*(-50 - d))
    BE = 9.81*(250 + d) + (42.7191 + 5.33989*sech2(0.612372*(250 + d)))*tanh(0.612372*(250 + d))
    
    CB = -22.6552*arctanh(0.707107*tanh(0.612372*(-50 - d))) + (21.3595 - 5.33988*sech2(0.612372*(-50 - d)))*tanh(0.612372*(-50 - d))
    CE = -22.6552*arctanh(0.707107*tanh(0.612372*(250 + d))) + (21.3595 - 5.33988*sech2(0.612372*(250 + d)))*tanh(0.612372*(250 + d))
    

    
    A = AE - AB  
    B = BE - BB    
    C = CE - CB


    #1527.68293
    return 0.5*(A + B + C)
    
def SolitonFD():
    return 1527.68293

if not os.path.exists(sdir):
    os.makedirs(sdir)
def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a


for k in range(6,21):
    
    wdirn = wdir + "/" + str(k) +  "/"     
    s = wdirn + "outlast.txt"
    with open(s) as fp:
        for i, line in enumerate(fp):
            if i == 2:
                nllist = line.split(",")
                dx = float(nllist[0])
                Ef = float(nllist[4])
            elif i > 2:
                break
    Ei = Soliton(0.5*dx)
    relErr = abs(Ei - Ef) / abs(Ei)
    
    s = sdir + "relE.dat"     
    with open(s,'a') as file3:
         s ="%3.8f%5s%1.15f\n" %(dx," ",relErr)
         file3.write(s)