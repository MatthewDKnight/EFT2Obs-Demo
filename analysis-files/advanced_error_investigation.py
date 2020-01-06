N = 1000000
#A = [55.351, 0, 0, 0,0]
#B = [[772.3750,0,173.7206,192.3112,0], [0,773.4795,0,0,220.4026], [0,0,12089.2548,741.6070,0], [0,0,0,1984.2937,0],[0,0,0,0,2045.2175]]

#ordering; cG, tcG, c2G, c3G, tc3G

def N_SM(x, A, B):
	return N/(1+abs(np.dot(A,x))+abs(np.dot(x, np.dot(B,x))))

def N_INT(x, A, B, param_index):
	i=param_index
	if A[i]!=0:
		return abs(A[i]*x[i]*N)/(1+abs(np.dot(A,x))+abs(np.dot(x, np.dot(B,x))))
	else:
		return None

def N_BSM(x, A, B, param_index1, param_index2):
	i=param_index1
	j=param_index2
	if B[i][j]!=0:
		return abs(B[i][j]*x[i]*x[j]*N)/(1+abs(np.dot(A,x))+abs(np.dot(x, np.dot(B,x))))
	else:
		return None

def f_SM(x, A, B):
	return 1/np.sqrt(N_SM(x, A, B))

def f_INT(x, A, B, param_index):
	if N_INT(x, A, B, param_index)!=None:
		return 1/np.sqrt(N_INT(x, A, B, param_index))
	else:
		return None

def f_BSM(x, A, B, param_index1, param_index2):
	if N_BSM(x, A, B, param_index1, param_index2)!=None:
		return 1/np.sqrt(N_BSM(x, A, B, param_index1, param_index2))
	else:
		return None

def f_A(x, A, B, param_index):
	if f_INT(x, A, B, param_index)!=None:
		return np.sqrt(f_SM(x, A, B)**2 + f_INT(x, A, B, param_index)**2)
	else:
		return 0

def f_B(x, A, B, param_index1, param_index2):
	if f_BSM(x, A, B, param_index1, param_index2)!=None:
		return np.sqrt(f_SM(x, A, B)**2 + f_BSM(x, A, B, param_index1, param_index2)**2)
	else:
		return 0

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import sys
sys.path.append('../../eftfitter-1')
from read2 import initialiseScalingFunctions


#x = np.linspace(0,1,10000)
#y1 = f_A(x)
#y2 = f_B(x)

#plt.plot(x,y1,label="INT")
#plt.plot(x,y2,label="BSM")
#plt.legend()
#plt.ylim(0,1)

def f_sum(x, *args):
	A = args[0][0]
        B = args[0][1]
	f_sum=0
	for i in range(len(A)):
		f_sum+=f_A(x, A, B, i)
		for j in range(len(A)):
			f_sum+=f_B(x, A, B, i,j)
	return f_sum
				

#sum = y1+y2
#f_min = f_sum[0]
#min_loc = 0
#for i in range(len(x)):
#        if f_sum[i] < f_min:
#                f_min = f_sum[i]
#                min_loc = x[i]

def minimize_f_sum(process_name, POIs):
	functions, name_ordering, scaling_list = initialiseScalingFunctions(POIs=POIs)
	if process_name in functions.keys():
		A=functions[process_name][0]
        	B=functions[process_name][1]
	else:
		print ("Process name not in list of known processes")
	for a in range(len(A)):
		if A[a]==0:
			A[a]=0.1
		for b in range(len(A)):
			if B[a][b]==0:
				B[a][b]=0.1 
	estimates=[0 for i in range(len(scaling_list))]
	f_min=9999
	for i in range(-2,3):
		guess=10**i
		init=[guess for j in range(len(scaling_list))]
		min_loc = scipy.optimize.minimize(f_sum, init, args=[A, B], tol=1e-8)
		if f_sum(min_loc.x, [A, B])<f_min:
			estimates=min_loc.x
			f_min=f_sum(estimates, [A, B])
			
	#f_min=f_sum(min_loc.x, [A, B])
	scalings=[10**i for i in scaling_list]
	name_ordering=[i.split('_')[0] for i in name_ordering]
	for i in range(len(name_ordering)):
		print(name_ordering[i])
		print(f_INT(min_loc.x, A, B, i))
		for j in range(len(name_ordering)):
			print(name_ordering[i])
			print(name_ordering[j])
			print(f_BSM(min_loc.x, A, B, i, j))
	print(A)
	print(B)	 
	return min_loc.x, f_min, name_ordering

estimates= minimize_f_sum('hww', ['cWWMinuscB', 'cHW', 'tcHW'])
print (estimates)

file1 = open("decay_default_estimates.txt","w")
L = [str(estimates[2][i]) + ": " + str(estimates[0][i]) + "\n"  for i in range(len(estimates[0]))]
print(L)
file1.writelines(L) 
file1.close() 

#plt.plot(x,f_sum)
#plt.plot([min_loc,min_loc],[0,1])
#plt.ylim(0,1)
#for i in range(5):
#                print (f_A(min_loc.x, i))
#                for j in range(5):
#                        print(f_B(min_loc.x,i,j))
#plt.show()
