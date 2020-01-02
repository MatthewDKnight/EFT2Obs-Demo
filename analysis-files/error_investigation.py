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
	return 0.01

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

plt.plot(x,y1,label="f_A")
plt.plot(x,y2,label="f_B")

sum = y1+y2
min = sum[0]
min_loc = 0
for i in range(len(x)):
        if sum[i] < min:
                min = sum[i]
                min_loc = x[i]
print(min_loc, min)

plt.plot(x,y1+y2, label = "f_A+f_B")
plt.plot([min_loc,min_loc],[0,1], label = "cWW=%f"%min_loc)

plt.ylim(0,0.2)
plt.xlabel("cWW")
plt.ylabel("Fractional error")
plt.legend()
plt.savefig("error_investigation_0.png")
plt.show()

def plotN(option):
	global N
	for N in [10,100,1000,10000]:
		x = np.linspace(0,1,100)
		y1 = f_A(x)
		y2 = f_B(x)
			
		if option==1:
			plt.plot(x,y1+y2, label="N=%d"%N)
		elif option==2:
			plt.plot(x,y1, label="N=%d"%N)
		elif option==3:
			plt.plot(x,y2, label="N=%d"%N)

	plt.plot([min_loc, min_loc],[0,2], label = "cWW=%f"%min_loc)
		
	plt.ylim(0,2)
	plt.xlabel("cWW")
	plt.ylabel("Fractional error")
	plt.legend()
	plt.savefig("error_investigation_%d.png"%option)
	plt.show()

plotN(1)
plotN(2)
plotN(3)
