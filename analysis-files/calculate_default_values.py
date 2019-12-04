from scipy.optimize import minimize


def calculate_default(A, B):
	return minimize(obj_func, x0 = 0.1, args=(A, B))
	


def obj_func(x, *args):
	A, B = args
	return (A*x-1)**2+(B*x**2-1)**2+(B*x**2-A*x)**2

A=10.26
B=26.3

print (calculate_default(A, B))
