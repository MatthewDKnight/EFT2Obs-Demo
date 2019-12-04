N = 10000
A = 10.26
B = 26.3

def N_SM(x):
	return N/(1+A*x+B*x**2)

def N_INT(x):
	return (A*x*N)/(1+A*x+B*x**2)

def N_BSM(x):
	return (B*x**2*N)/(1+A*x+B*x**2)

def f_SM(x):
	return 1/np.sqrt(N_SM(x))

def f_INT(x):
	return 1/np.sqrt(N_INT(x))

def f_BSM(x):
	return 1/np.sqrt(N_BSM(x))

def f_A(x):
	return np.sqrt(f_SM(x)**2 + f_INT(x)**2)

def f_B(x):
	return np.sqrt(f_SM(x)**2 + f_BSM(x)**2)

import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0,1,10000)
y1 = f_A(x)
y2 = f_B(x)

plt.plot(x,y1,label="INT")
plt.plot(x,y2,label="BSM")
plt.legend()
plt.ylim(0,1)

sum = y1+y2
min = sum[0]
min_loc = 0
for i in range(len(x)):
        if sum[i] < min:
                min = sum[i]
                min_loc = x[i]
print(min_loc, min)

plt.plot(x,y1+y2)
plt.plot([min_loc,min_loc],[0,1])
plt.ylim(0,1)
plt.show()
