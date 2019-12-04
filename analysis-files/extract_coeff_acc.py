import yoda
import matplotlib.pyplot as plt
import numpy as np
import json
from parameter_conversion import params as param_convert #dictionary to switch naming conventions



#def getBins(hist):
#	bin_values = []
#       
	#merge bins into 0,20,45,80,120,200,400 
#	hist.mergeBins(0,3)
#	hist.mergeBins(1,5)
#	hist.mergeBins(2,8)
#	hist.mergeBins(3,10)
#	hist.mergeBins(4,19)
#	hist.mergeBins(5,len(hist)-1)
#	
#	for bin in hist:
#		bin_values.append(bin.sumW)
#	
#	overflow = hist.overflow.sumW
#	bin_values[-1] += overflow
#
#	return bin_values


def getBins(hist):
	bin_values = []
       
	for bin in hist:
		bin_values.append(bin.sumW)
	
	overflow = hist.overflow.sumW
	bin_values.append(overflow)
        #bin_values.append(hist.sumW())
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


	no_params = len(pars)
	no_bins = len(bin_values_no_acc[0])
	

	#print(bin_values_acc)
	#print(bin_values_no_acc)
        acceptance_values=bin_values_acc/bin_values_no_acc
	acceptance_values=acceptance_values/acceptance_values[0]
	acceptance_values=acceptance_values-1
	#print(acceptance_values)
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
		
def writeTextFile(coeffs):
	with open("equations.txt", "w") as file:
		print(len(coeffs))
		for j in range(len(coeffs)):
			bin = coeffs[j]
			string = "%s:1"%bin_names[j]
			for i in range(len(bin)):
				param_name = param_convert[pars[i]['name']]
				A = bin[i][0]
				B = bin[i][1]
				#make exception for cG and cA parameters
				if param_name in ["cG","tcG", "cA"]:
					A = A / (4*np.pi)**2
					B = B / (4*np.pi)**4
				if A>0:
					string += " + %f * %s"%(A, param_name)
				else:
					string += " %f * %s"%(A, param_name)
				if B>0:
					string += " + %f * %s * %s"%(B, param_name, param_name)
				else:
					string += " %f * %s * %s"%(B, param_name, param_name)
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
		

aos = yoda.read("SimpleHiggs.yoda", asdict = False)
hists_no_acc = [h for h in aos if h.path.startswith("/SimpleHiggs/H_PT")]
hists_no_acc = hists_no_acc[1:] #get rid of first histogram
hists_acc = [h for h in aos if h.path.startswith("/SimpleHiggs/acc_H_PT")]
hists_acc = hists_acc[1:] #get rid of first histogram


#have to order histograms with rw0 at start
hists_no_acc = orderHists(hists_no_acc)
hists_acc = orderHists(hists_acc)

#Bins for pT are 0,20,45,80,120,200
bin_edges = [0,20,45,80,120,200]
bin_names = ["GG2H_PTH_0_20", "GG2H_PTH_20_45", "GG2H_PTH_45_80", "GG2H_PTH_80_120", "GG2H_PTH_120_200", "GG2H_PTH_GT200"] 

with open("config.json") as jsonfile:
	pars = json.load(jsonfile)

coeffs = extractCoeff(hists_no_acc, hists_acc)
writeTextFile(coeffs)

#plot(hists[0])

###################################################################################################
"""

xs_list_no_acc=[]
xs_list_acc=[]
for i in aos:
	#if 'H_PT[rw' in i.name:
	if i.path.startswith("/SimpleHiggs/N_jets"):
		bin_list_no_acc=[]
		num_bins=len(i)
		for j in range(len(i)):
			#print(i[j])
			bin_list_no_acc.append(i[j].sumW)
		xs_list_no_acc.append(bin_list_no_acc)
	if i.path.startswith("/SimpleHiggs/acc_N_jets"):
		bin_list_acc=[]
		for j in range(len(i)):
			#print(i[j])
			bin_list_acc.append(i[j].sumW)
		xs_list_acc.append(bin_list_acc)

		
xs_list_no_acc=xs_list_no_acc[1:4]
xs_list_acc=xs_list_acc[1:4]

cWW_values=[0.0, 0.5, 1.0]

x = np.linspace(0,1,10000)


for k in range(num_bins):
	cWW_xs_values_no_acc=np.array([i[k] for i in xs_list_no_acc])
	cWW_xs_values_acc=np.array([i[k] for i in xs_list_acc])
	#cHW_xs_values=[xs_list[0], xs_list[3], xs_list[4]]
	acceptance=cWW_xs_values_acc/cWW_xs_values_no_acc
	#print(cWW_xs_values_acc)
	acceptance_SM=acceptance[0]
	normed_acceptance=acceptance/acceptance_SM

	p = np.polyfit(cWW_values, normed_acceptance, 2)
	print(p)
	#plt.plot(cWW_values, normed_acceptance, 'bo')
	#plt.plot(x, np.polyval(p,x))
	#plt.show()
"""
