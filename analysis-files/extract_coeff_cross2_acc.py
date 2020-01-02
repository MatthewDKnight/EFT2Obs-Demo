import yoda
import matplotlib.pyplot as plt
import numpy as np
import json
from parameter_conversion import params as param_convert #dictionary to switch naming conventions

def getBins(hist):
	bin_values = []
       
	for bin in hist:
		bin_values.append(bin.sumW)
	
	overflow = hist.overflow.sumW
	bin_values.append(overflow)
        bin_values.append(hist.sumW())
	return bin_values
	

def plot(hist):
	bin_values = getBins(hist)
	bin_values = np.array(bin_values)
	norm = bin_values[0]
	bin_values = (bin_values/norm) * 27.45
	
	widths = [] 
	for i in range(len(bin_edges)-1):
		widths.append(bin_edges[i+1]-bin_edges[i])

	left = bin_edges[:len(bin_edges)-1] #remove last bin edge
	height = bin_values[:len(bin_values)-1] #remove >200 bin value

	print(left)
	print(height)
	print(widths)

	plt.bar(left, height, widths, align="edge")
	plt.show()

def extractCoeff(hists_no_acc, hists_acc):
	bin_values_no_acc = [getBins(hist) for hist in hists_no_acc] #values of bins for each reweighting
	bin_values_no_acc = np.array(bin_values_no_acc)
	bin_values_no_acc = bin_values_no_acc/bin_values_no_acc[0] #divide by SM value -> yields
	#bin_values_no_acc = bin_values_no_acc - 1 #now left with INT and BSM terms

	bin_values_acc = [getBins(hist) for hist in hists_acc] #values of bins for each reweighting
	bin_values_acc = np.array(bin_values_acc)
	bin_values_acc = bin_values_acc/bin_values_acc[0] #divide by SM value -> yields
	#bin_values_acc = bin_values_acc - 1 #now left with INT and BSM terms


	#bin_values = [getBins(hist) for hist in hists] #values of bins for each reweighting
	#bin_values = np.array(bin_values)
	#bin_values = bin_values/bin_values[0] #divide by SM value -> yields
	#bin_values = bin_values - 1 #now left with INT and BSM terms

	no_params = len(pars)
	no_bins = len(bin_values_acc[0])
	
	#print(no_params)
	#print(no_bins)
	#print(bin_values_acc)
	#print(bin_values_no_acc)
        acceptance_values=bin_values_acc/bin_values_no_acc
	acceptance_values=acceptance_values/acceptance_values[0]
	acceptance_values=acceptance_values-1
	#print(acceptance_values)

	

	#create a to carry all coefficients of for every bin
	coeffs = [ [[0.0,0.0] for i in range(no_params)] for j in range(no_bins)]
	coeffs = np.array(coeffs)
	
	for i in range(1, no_params+1): #for each parameter
		print("Solving for %s" % pars[i-1]['name'])
		step = pars[i-1]['step']
		for j in range(no_bins):
			#values of INT+BSM terms for bin in two steps of param
			mu1 = acceptance_values[i*2-1][j]
			mu2 = acceptance_values[i*2][j]
			
			A = (4*mu1-mu2)/step
			B = (mu2-2*mu1)/(2*(step/2)**2)
			print(A, B)

			coeffs[j][i-1][0]=A
			coeffs[j][i-1][1]=B

	#extract cross terms
	cross_terms = np.zeros((no_bins, no_params, no_params))

	counter = 1 + no_params*2
	print(counter)
	print(no_bins)
	print(len(acceptance_values))
	for i in range(no_params):
		for j in range(no_params):
			#print(i, j)
			#print(counter)
			if i!=j and j>i:
				step1 = pars[i]['step']
				step2 = pars[j]['step']
				#print(no_bins)
				#print(counter)
				for k in range(no_bins):
					mu = acceptance_values[counter-1][k]
				
					#remove linear and quadratic (not cross) terms
					mu += -coeffs[k][i][0]*step1 -coeffs[k][i][1]*step1**2 -coeffs[k][j][0]*step2 - coeffs[k][j][1]*step2**2
					cross_term = mu / (step1*step2)
 					print(cross_term)	
					cross_terms[k][i][j] = cross_term
				counter += 1

	return coeffs, cross_terms
		

"""
def extractCoeff(hists_no_acc, hists_acc):
	bin_values_no_acc = [getBins(hist) for hist in hists_no_acc] #values of bins for each reweighting
	bin_values_no_acc = np.array(bin_values_no_acc)
	bin_values_no_acc = bin_values_no_acc/bin_values_no_acc[0] #divide by SM value -> yields
	#bin_values_no_acc = bin_values_no_acc - 1 #now left with INT and BSM terms

	bin_values_acc = [getBins(hist) for hist in hists_acc] #values of bins for each reweighting
	bin_values_acc = np.array(bin_values_acc)
	bin_values_acc = bin_values_acc/bin_values_acc[0] #divide by SM value -> yields
	#bin_values_acc = bin_values_acc - 1 #now left with INT and BSM terms


	no_params = len(pars)
	no_bins = len(bin_values_no_acc[0])
	

	print(bin_values_acc)
	print(bin_values_no_acc)
        acceptance_values=bin_values_acc/bin_values_no_acc
	acceptance_values=acceptance_values/acceptance_values[0]
	acceptance_values=acceptance_values-1
	print(acceptance_values)
	#print(no_params)
	#print(no_bins)

	#create a to carry all coefficients of for every bin
	coeffs = [ [[0.0,0.0] for i in range(no_params)] for j in range(no_bins)]
	coeffs = np.array(coeffs)
	
	for i in range(1, no_params+1): #for each parameter
		print("Solving for %s" % pars[i-1]['name'])
		step = pars[i-1]['step']
		for j in range(no_bins):
			#values of INT+BSM terms for bin in two steps of param
			mu1 = acceptance_values[i*2-1][j]
			mu2 = acceptance_values[i*2][j]
			
			A = (4*mu1-mu2)/step
			B = (mu2-2*mu1)/(2*(step/2)**2)
			print(A, B)

			coeffs[j][i-1][0]=A
			coeffs[j][i-1][1]=B

	return coeffs
"""


def writeTextFile(coeffs, cross_terms, filename="equations.txt"):
	with open(filename, "w") as file:
		for j in range(len(coeffs)): #for each bin
			bin = coeffs[j]	     #the A and B coeffs for every parameter for this bin
			string = "%s:1"%bin_names[j]

			for i in range(len(bin)): #for each parameter
				param_name = param_convert[pars[i]['name']] #convert naming convention
				#extract A and B coefficients for this parameter and bin
				A = bin[i][0]
				B = bin[i][1]

				#make exception for cG and cA parameters
				if param_name in ["cG","tcG", "cA"]:
					A = A / (4*np.pi)**2
					B = B / (4*np.pi)**4
				if abs(A) > 10e-5:
					if A>0:
						string += " + %f * %s"%(A, param_name)
					elif A<0:
						string += " %f * %s"%(A, param_name)
				if abs(B) > 10e-5:
					if B>0:
						string += " + %f * %s * %s"%(B, param_name, param_name)
					elif B<0:
						string += " %f * %s * %s"%(B, param_name, param_name)
			
				#lets add the cross terms
				for k in range(len(bin)): #for every parameter
					cross_term = cross_terms[j][i][k]
					if abs(cross_term) > 10e-5:
						param_name2 = param_convert[pars[k]['name']]
						if param_name in ["cG", "tcG", "cA"]:
							cross_term = cross_term / (4*np.pi)**2
						if param_name2 in ["cG", "tcG", "cA"]:
							cross_term = cross_term / (4*np.pi)**2
						if cross_term > 0:
							string += " + %f * %s * %s"%(cross_term, param_name, param_name2) 
						elif cross_term < 0:	
							string += " %f * %s * %s"%(cross_term, param_name, param_name2)	
			string += "\n"		
			file.write(string)

def getHistNumber(hist):
        split_str = hist.name.split("w")[1]
        split_str = split_str[:len(split_str)-1]
        return int(split_str)

def orderHists(hists):
        ordered = False
        while ordered==False:
                ordered=True
                for i in range(len(hists)-1):
                        if getHistNumber(hists[i]) > getHistNumber(hists[i+1]):
                                #swap histograms
                                ordered = False
                                hists[i], hists[i+1] = hists[i+1], hists[i]
        return hists

aos = yoda.read("../acceptance_equations_3.yoda", asdict = False)
H_PT_hists = [h for h in aos if h.path.startswith("/SimpleHiggs/H_PT")]
H_PT_hists = H_PT_hists[1:] #get rid of first histogram

acc_H_PT_hists = [h for h in aos if h.path.startswith("/SimpleHiggs/acc_H_PT")]
acc_H_PT_hists = acc_H_PT_hists[1:] #get rid of first histogram

N_jets_hists = [h for h in aos if h.path.startswith("/SimpleHiggs/N_jets")]
N_jets_hists = N_jets_hists[1:] #get rid of first histogram

acc_N_jets_hists = [h for h in aos if h.path.startswith("/SimpleHiggs/acc_N_jets")]
acc_N_jets_hists = acc_N_jets_hists[1:] #get rid of first histogram


#have to order histograms with rw0 at start
H_PT_hists = orderHists(H_PT_hists)
acc_H_PT_hists = orderHists(acc_H_PT_hists)
N_jets_hists = orderHists(N_jets_hists)
acc_N_jets_hists = orderHists(acc_N_jets_hists)


#Bins for pT are 0,20,45,80,120,200
bin_edges = [0,20,45,80,120,200]
bin_names = ["GG2H_PTH_0_20", "GG2H_PTH_20_45", "GG2H_PTH_45_80", "GG2H_PTH_80_120", "GG2H_PTH_120_200", "GG2H_PTH_GT200"] 
bin_names.append("Total")

with open("../config.json") as jsonfile:
	pars = json.load(jsonfile)

coeffs, cross_terms = extractCoeff(H_PT_hists, acc_H_PT_hists)
writeTextFile(coeffs, cross_terms, "../acceptance_equations_3.txt")

#coeffs, cross_terms = extractCoeff(acc_H_PT_hists)
#writeTextFile(coeffs, cross_terms, "equations_acc.txt")


#plot(hists[0])
