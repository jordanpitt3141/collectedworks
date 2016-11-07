import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

order = "3p"
wdir = "../../results/segurf/o"+order+"/"
expdir = "../../results/segurexp/"
sdir = "../../../exportpic/segurf/o"+order+"/"

#poss = [0,5]
#poss = [0]
poss = [0,5,10,15,20]
for pos in poss:
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
    
    s = expdir + "Seguroutxe"+str(int(pos))+".csv"
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
    smh = 1.5*(mh - h0) / (h0)
    
    plot(smt,smh ,'.r')
    plot(et,eh,'-b')  
    ylim([-0.13,0.06])       
    xlim([-10,250])
    #s = order + "order Segur at x/h1 =" + str(pos/h0)
    #title(s)
    #xlabel("Scaled time [(t*sqrt(g/h1)) - x/h1]")
    #ylabel("Scaled height [(h- h1)/h1]")
    s = sdir +"bothx"+ str(int(pos)) +".png"       
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