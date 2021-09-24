import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import numpy as np


N = np.array([2, 4, 8, 16, 32, 64, 128, 256])
tridiag_iter = np.array([1, 6, 77, 319, 1233, 4480, 14858, 40868])
dense_iter = np.array([1, 15, 90, 394, 1686, 7034, 29031, 118626])

log_2_N = np.log2(N)
log_2_tridiag_iter = np.log2(tridiag_iter)
log_2_dense_iter = np.log2(dense_iter)

plt.plot(log_2_N, log_2_tridiag_iter, label='Tridiag')
plt.plot(log_2_N, log_2_dense_iter, label='Dense')
plt.xlabel(r'$Log_2(N)$')
plt.ylabel(r'$Log_2(Iterations)$') 
plt.title(r'The number of transformations versus matrix size')
plt.legend()
plt.savefig("IterationsvsN.pdf")
plt.show()
