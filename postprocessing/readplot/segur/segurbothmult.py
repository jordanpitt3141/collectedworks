import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog


orders = ["2","3","3p"]
poss = [0,5,10,15,20]

for pos in poss:
    smts = []
    smhs =[]
    for order in orders:
        wdir = "../../results/segur3/o"+order+"/"
        sdir = "../../results/segur3/o"+order+"/"
        h0 = 0.1
        g = 9.81
             
        s = wdir + "out"+str(pos) + ".txt"
        with open(s,'r') as file1:
             readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
             mh = []
             mt = []
             j = -1
             for row in readfile:       
                 if (j >= 0):
                    print(pos,j,row[0])
                    mdx = float(row[0])
                    mdt = float(row[1])
                    mt.append(float(row[2]))
                    mh.append(float(row[3]))
                    
                 j = j + 1
             mt = array(mt)
             mh = array(mh)
        
        s = wdir + "Seguroutxe"+str(int(pos))+".csv"
        with open(s,'r') as file1:
             readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
             eh = []
             et = []
             j = -1
             for row in readfile:       
                 if (j >= 0):
                    et.append(float(row[0]))
                    eh.append(float(row[1]))
                    
                 j = j + 1
             et = array(et)
             eh = array(eh)
             
        #scaling
        smts.append(mt*sqrt(g/h0) - (pos/h0))
        smhs.append(1.5*(mh - h0) / (h0))
    for i in range(len(orders)):
        plot(smts[i],smhs[i] ,'-',label=orders[i])
    plot(et,eh,'-k',label="data")         
    xlim([-10,250])
    s = "Segur at x/h1 =" + str(pos/h0)
    title(s)
    xlabel("Scaled time [(t*sqrt(g/h1)) - x/h1]")
    ylabel("Scaled height [(h- h1)/h1]")
    legend()
    s = sdir +"3hx"+ str(int(pos)) +".png"       
    savefig(s, bbox_inches='tight')        
    clf()
"""
    #plot(et,eh)         
    plot(smt,smh ,'.r')
    xlim([-10,250])
    s = order + "order Segur at x/h1 =" + str(pos/h0)
    title(s)
    xlabel("Scaled time [(t*sqrt(g/h1)) - x/h1]")
    ylabel("Scaled height [(h- h1)/h1]")
    s = sdir +"nmx" + str(int(pos)) + ".png"       
    savefig(s, bbox_inches='tight')        
    clf()
    
    #plot(et,eh)         
    plot(smt,smh ,'.r')
    xlim([-10,250])
    s = order + "order Segur at x/h1 =" + str(pos/h0)
    title(s)
    xlabel("time")
    ylabel("height")
    s = sdir +"nmux"+ str(int(pos)) +".png"       
    savefig(s, bbox_inches='tight')        
    clf()
    
    plot(et,eh)         
    #plot(smt,smh ,'.r')
    xlim([-10,250])
    s = order + "order Segur at x/h1 =" + str(pos/h0)
    title(s)
    xlabel("Scaled time [(t*sqrt(g/h1)) - x/h1]")
    ylabel("Scaled height [(h- h1)/h1]")
    s = sdir +"exx"+ str(int(pos)) +".png"       
    savefig(s, bbox_inches='tight')        
    clf()
"""