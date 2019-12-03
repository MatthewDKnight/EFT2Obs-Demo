import yoda
import matplotlib.pyplot as plt
import numpy as np
import json
from parameter_conversion import params as param_convert #dictionary to switch naming conventions

def getBins(hist):
	bin_values = []
       
	#merge bins into 0,20,45,80,120,200,400 
	hist.mergeBins(0,3)
	hist.mergeBins(1,5)
	hist.mergeBins(2,8)
	hist.mergeBins(3,10)
	hist.mergeBins(4,19)
	hist.mergeBins(5,len(hist)-1)
	
	for bin in hist:
		bin_values.append(bin.sumW)
	
	overflow = hist.overflow.sumW
	bin_values[-1] += overflow

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

def extractCoeff(hists):
	bin_values = [getBins(hist) for hist in hists] #values of bins for each reweighting
	bin_values = np.array(bin_values)
	bin_values = bin_values/bin_values[0] #divide by SM value -> yields
	bin_values = bin_values - 1 #now left with INT and BSM terms

	no_params = len(pars)
	no_bins = len(bin_values[0])
	
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
			mu1 = bin_values[i*2-1][j]
			mu2 = bin_values[i*2][j]
			
			A = (4*mu1-mu2)/step
			B = (mu2-2*mu1)/(2*(step/2)**2)
			print(A, B)

			coeffs[j][i-1][0]=A
			coeffs[j][i-1][1]=B

	return coeffs
		
def writeTextFile(coeffs):
	with open("equations.txt", "w") as file:
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
		

aos = yoda.read("Rivet.yoda", asdict = False)
hists = [h for h in aos if h.path.startswith("/HiggsTemplateCrossSectionsStage1/pT_Higgs")]
hists = hists[1:] #get rid of first histogram

#have to order histograms with rw0 at start
hists = orderHists(hists)

#Bins for pT are 0,20,45,80,120,200
bin_edges = [0,20,45,80,120,200]
bin_names = ["GG2H_PTH_0_20", "GG2H_PTH_20_45", "GG2H_PTH_45_80", "GG2H_PTH_80_120", "GG2H_PTH_120_200", "GG2H_PTH_GT200"] 

with open("config.json") as jsonfile:
	pars = json.load(jsonfile)

coeffs = extractCoeff(hists)
writeTextFile(coeffs)

#plot(hists[0])
