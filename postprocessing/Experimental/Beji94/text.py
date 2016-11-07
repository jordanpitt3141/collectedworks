from struct import unpack
from matplotlib.pyplot import plot,ylim
import csv


wdir = "../../../../data/Experimental/Data 1994 Paper/DATA 1994 TEXT/"
sdir = "../../../../data/Experimental/Data 1994 Paper/CSV/"
ts = []
wg1s = []
wg2s = []
wg3s = []
wg4s = []
wg5s = []
wg6s = []
wg7s = []
wg8s= []

exp = "sh"
 
s = wdir + exp+".dat"
with open(s, mode='r') as file1: # b is important -> binary
    j = -2
    for line in file1:
        if(j < -1):
            ls = line.split('  ')
            lsrs = [s for s in ls if s != '']
            lsrscf = map(str, lsrs)
            #print(lsrscf)
        if(j >0):
            ls = line.split('  ')
            lsrs = [s for s in ls if s != '']
            lsrscf = map(float, lsrs)
            #print(lsrscf)
            t,wg1,wg2,wg3,wg4,wg5,wg6,wg7 = lsrscf
            ts.append(t)
            wg1s.append(wg1)
            wg2s.append(wg2)
            wg3s.append(wg3)
            wg4s.append(wg4)
            wg5s.append(wg5)
            wg6s.append(wg6)
            wg7s.append(wg7)
            #print(lsrscf)
        j= j+ 1
    n = len(ts)
    s = sdir + exp+".csv"
    with open(s,'a') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow(['Time(s)' ,'wg1(cm)','wg2(cm)','wg3(cm)', 'wg4(cm)', 'wg5(cm)' , 'wg6(cm)' ,'wg7(cm)'])        
               
        for j in range(n):
            writefile2.writerow([str(ts[j]),str(wg1s[j]),str(wg2s[j]),str(wg3s[j]) ,str(wg4s[j]) , str(wg5s[j]) , str(wg6s[j]), str(wg7s[j])])     
    

   