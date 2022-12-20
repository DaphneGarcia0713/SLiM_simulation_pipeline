import sys


name = sys.argv[1]       #this is array of command line
filehandle = open("vcf/"+name+".vcf",'r')     #r means read only. filehandle is file object

f = open('RealEffectValueFiles/'+name+".real", 'a',)                      #change this file name each time
for line in filehandle:
    if line[0] == '#':
        continue              #this means "skip"
    data = line.strip().split()     
    INFO = data[7].strip().split()
    y = INFO[0].split(';')
    s=(y[1].lstrip('S='))
    f.write(s+ "\n")
    
    #this script takes in the vcf from command line, and puts the real trait values formatted into file f