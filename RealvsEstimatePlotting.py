import sys

# This file inputs real files and estimated files in terminal.
# gets the first column of real values (only one column), and 
# stores it in real. Gets second line of estimated values (bc 
# first line is p value), stores it in est. Opens new file 
# and for i in est,writes real and est in 2 columns. This gets
# the right real value for its estimate bc it gets the 
# estimates's number of mutation?? double check this??


# get real values

name = sys.argv[1]

print("Usually takes less than a min")
fh = open("RealEffectValueFiles/"+name+".real",'r')   #I think the r means read only   #argv1 should be real 

real = []

for line in fh:

    real.append(float(line.strip()))

fh.close()


# get est values

fh = open("EstimatedEffectValueFiles/"+name+".estim",'r')


est = {}

for line in fh:

    data = line.strip().split()

    est[int(data[2])-1] = float(data[1])

fh.close()



for ind in est:

    with open("RealvsEstimatePlotCoordinates/"+name+".coord", 'a') as f:          ##CHANGE THIS
        print(est[ind], real[ind], file=f)
        
print("done :)")        
        
#how to make this better: change ind to i, get rid of spaces