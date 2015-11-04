import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

wdir = "../../results/segur2/o2/"
sdir = wdir

poss = [0,5,10,15,20]
for pos in poss:
    h0 = 0.1
    g = 9.81
         
    s = sdir + "out"+str(pos) + ".txt"
    with open(s,'r') as file1:
         readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
         mh = []
         mt = []
         j = -1
         for row in readfile:       
             if (j >= 0):
                print(pos,j)
                mdx = float(row[0])
                mdt = float(row[1])
                mt.append(float(row[2]))
                mh.append(float(row[3]))
                
             j = j + 1
         mt = array(mt)
         mh = array(mh)
    
    s = sdir + "Seguroutxe"+str(int(pos))+".csv"
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
    smt = mt*sqrt(g/h0) - (pos/h0)
    smh = (mh - h0) / (h0)
    
    plot(et,eh)         
    plot(smt,smh ,'.r')
    xlim([-10,250])
    s = "Segur at x/h1 =" + str(pos/h0)
    title(s)
    xlabel("Scaled time [(t*sqrt(g/h1)) - x/h1]")
    ylabel("Scaled height [(h- h1)/h1]")
    s = sdir +"bothx"+ str(int(pos)) +".png"       
    savefig(s, bbox_inches='tight')        
    clf()
    
    #plot(et,eh)         
    plot(smt,smh ,'.r')
    xlim([-10,250])
    s = "Segur at x/h0 =" + str(pos/h0)
    title(s)
    xlabel("Scaled time [(t*sqrt(g/h1)) - x/h1]")
    ylabel("Scaled height [(h- h1)/h1]")
    s = sdir +"nmx" + str(int(pos)) + ".png"       
    savefig(s, bbox_inches='tight')        
    clf()
    
    #plot(et,eh)         
    plot(smt,smh ,'.r')
    xlim([-10,250])
    s = "Segur at x/h0 =" + str(pos/h0)
    title(s)
    xlabel("time")
    ylabel("height")
    s = sdir +"nmux"+ str(int(pos)) +".png"       
    savefig(s, bbox_inches='tight')        
    clf()
    
    plot(et,eh)         
    #plot(smt,smh ,'.r')
    xlim([-10,250])
    s = "Segur at x/h0 =" + str(pos/h0)
    title(s)
    xlabel("Scaled time [(t*sqrt(g/h1)) - x/h1]")
    ylabel("Scaled height [(h- h1)/h1]")
    s = sdir +"exx"+ str(int(pos)) +".png"       
    savefig(s, bbox_inches='tight')        
    clf()