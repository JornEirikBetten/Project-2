import numpy as np

def eigval(val, N): 
	return 1./N**2*(2-2*np.cos(val*np.pi/7.))

eigvals = []

for i in range(1,7,1): 
	temp = float(i)
	eigvals.append(eigval(temp, 6))

print(eigvals) 
