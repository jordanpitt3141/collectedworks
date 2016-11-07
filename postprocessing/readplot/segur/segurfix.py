import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

order = "1af0f"
wdir = "../../../../data/postprocessing/segur/o"+order+"/"
expdir = "../../../../data/postprocessing/segurexp/"
sdir = "../../../../data/postprocessing/segurfull/o"+order+"/"

if not os.path.exists(sdir):
    os.makedirs(sdir) 

gap = 1
poss = [0]
#poss = [0,5,10,15,20]
#poss = [5,10,15,20]
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
    
    n = len(smh)
    m = len(eh)
    
    #write to a file
    s = sdir + "msegurp"+ str(int(pos)) +  ".dat"
    with open(s,'a') as file2:
        for i in range(0,n,gap):
            if(smt[i] <= 250 and smt[i] > -50):
                s ="%3.8f%5s%1.15f\n" %(smt[i]," ",smh[i])
                file2.write(s)
    s = sdir + "esegurp"+ str(int(pos)) +  ".dat"
    with open(s,'a') as file2:
        for i in range(m):
            s ="%3.8f%5s%1.15f\n" %(et[i]," ",eh[i])
            file2.write(s)
