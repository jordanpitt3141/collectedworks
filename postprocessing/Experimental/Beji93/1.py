from struct import unpack



wdir = "../../../../data/Experimental/Data 1993 Paper/DATA 1993 BINARY/"

 
s = wdir + "JLN"
with open(s, mode='rb') as file1: # b is important -> binary
    fileContent = file1.read()
    h = unpack('@ffffffff',fileContent[30:])



   