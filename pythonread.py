import pandas as pd
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import numpy as np

def eigvec(i, N): 
	eigvecN = []
	for j in range(0, N+2, 1) : 
		val = float(j)
		eigvecN.append(np.sin(val*i*np.pi/(float(N)+1)))
	eigvecN = np.array(eigvecN)
	return eigvecN
true_eigvecN101 = eigvec(1.0, 10)
true_eigvecN102 = eigvec(2.0, 10)
true_eigvecN103 = eigvec(3.0, 10)

true_eigvecN1001 = eigvec(1.0, 100)
true_eigvecN1002 = eigvec(2.0, 100)
true_eigvecN1003 = eigvec(3.0, 100)

sN10 = 0
sN100 = 0
for i in range(len(true_eigvecN101)): 
	sN10 += true_eigvecN101[i]*true_eigvecN101[i]



for i in range(len(true_eigvecN1001)): 
	sN100 += true_eigvecN1001[i]*true_eigvecN1001[i]


data = pd.read_csv("eigensolve_N_10.csv")

x = data["x"]
eigvec1 = data["0.0008101"]*np.sqrt(sN10)
eigvec2 = data["0.003175"]*np.sqrt(sN10)
eigvec3 = data["0.006903"]*np.sqrt(sN10)


plt.title("Plot of first three eigenvectors, with N=10", fontsize=20) 
plt.plot(x, eigvec1, label=r'Eigenvector of $\lambda_1=0.0008$')
plt.plot(x, eigvec2, label=r'Eigenvector of $\lambda_2=0.0032$') 
plt.plot(x, eigvec3, label=r'Eigenvector of $\lambda_3=0.0069$') 
plt.plot(x, true_eigvecN101, label=r'Exact eigenvector of $\lambda_1$')
plt.plot(x, true_eigvecN102, label=r'Exact eigenvector of $\lambda_2$') 
plt.plot(x, true_eigvecN103, label=r'Exact eigenvector of $\lambda_3$') 
plt.xlabel(r'$\hat{x}_i$', fontsize=18)
plt.ylabel(r'$v_i$', fontsize=18)
plt.legend(loc=4)
plt.savefig("EigenvectorsN10.pdf")
plt.show()

data1 = pd.read_csv("eigensolve_N_100.csv")

x = data1["x"]
eigvec1 = data1["9.677e-08"]*np.sqrt(sN100)
eigvec2 = data1["3.869e-07"]*np.sqrt(sN100)
eigvec3 = data1["8.702e-07"]*np.sqrt(sN100)


plt.title("Plot of first three eigenvectors, with N=100", fontsize=20) 
plt.plot(x, eigvec1, label=r'Eigenvector with $\lambda_1=9.6775e-08$')
plt.plot(x, eigvec2, label=r'Eigenvector with $\lambda_2=3.8689e-07$') 
plt.plot(x, eigvec3, label=r'Eigenvector with $\lambda_3=8.7015e-07$') 
plt.plot(x, true_eigvecN1001, label=r'Exact eigenvector of $\lambda_1$')
plt.plot(x, true_eigvecN1002, label=r'Exact eigenvector of $\lambda_2$') 
plt.plot(x, true_eigvecN1003, label=r'Exact eigenvector of $\lambda_3$') 
plt.xlabel(r'$\hat{x}_i$', fontsize=18)
plt.ylabel(r'$v_i$', fontsize=18)
plt.legend(loc=4)
plt.savefig("EigenvectorsN100.pdf")
plt.show()
