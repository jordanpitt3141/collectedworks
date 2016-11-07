import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog



def TV(a):
    return sum(abs(a[1:] - a[:-1]))

nf = 18    
orders = ["1af","2af","3"]     
#orders = ["1"]     
for order in orders:
    wdir = "../../../../data/raw/ndbh/o"+order+"/"
    sdir = "../../../../data/postprocessing/ndbh/o"+order+"/"
    if not os.path.exists(sdir):
        os.makedirs(sdir)
    L1h = []
    L1u = []
    dxs = []
    for i in range(nf):         
        s = wdir + str(i) + ".txt"
        with open(s,'r') as file1:
            readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            xm = []    
            hm = []
            um = []
            j = -1
            for row in readfile:       
                if (j >= 0):
                    dx =float(row[0])
                    dt = float(row[1])
                    time = float(row[2])
                    xm.append(float(row[3]))
                    hm.append(float(row[4]))
                    um.append(float(row[6]))
                        
                j = j + 1
            um = array(um)
            hm = array(hm)
            xm = array(xm)
            
            dxs.append(dx)
            L1h.append(TV(hm))
            L1u.append(TV(um))
    loglog(dxs,L1h,'o',label=order)
    s = sdir + "h.dat"
    n = len(L1h)
    with open(s,'a') as file1:
        for i in range(n):
            s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",L1h[i])
            file1.write(s)
    s = sdir + "u.dat"
    with open(s,'a') as file2:
        for i in range(n):
            s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",L1u[i])
            file2.write(s)
title("Log-Log Plot of TV for Dambreak")
xlabel("dx")    
ylabel("TV") 
legend()