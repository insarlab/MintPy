import numpy as np

data = np.loadtxt('GSRM_plate_outlines.gmt',dtype=str)
data = np.flip(data,1)
# Locate the starting position of each plate
bnds_index, = np.where(data[:,1] == '>')
n = len(bnds_index)

# Separate the boundaries of each plate and write it in a file
for i in range(n):
    vi = bnds_index[i]
    j1 = vi+1
    if i == n-1:
        np.savetxt(data[vi][0], data[j1:-1], fmt='%10s %9s',header=data[vi][0],comments='')
    else:
        j2 = bnds_index[i+1]
        np.savetxt(data[vi][0], data[j1:j2], fmt='%10s %9s',header=data[vi][0],comments='')
