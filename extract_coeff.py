import yoda
import matplotlib.pyplot as plt
import numpy as np
import json

def getBins(hist):
	bin_centres = []
	bin_values = []
	widths = []
	errors = []	
	
	bin_value = 0
	sum_error2 = 0
	for small_bin in hist:
		small_edges = small_bin.xEdges
		if small_edges[1] in bin_edges:
			bin_value = bin_value + small_bin.area
			sum_error2 = sum_error2 + small_bin.areaErr**2
			bin_values.append(bin_value)
			errors.append(np.sqrt(sum_error2))
			bin_value = 0
			sum_error2 = 0
			
			outer_edge_index = bin_edges.index(small_edges[1])
			edges = [bin_edges[outer_edge_index-1], bin_edges[outer_edge_index]]
			bin_centre = (float(edges[0])+float(edges[1]))/2
			bin_centres.append(bin_centre)

			width = edges[1] - edges[0]
			widths.append(width)
		else:
			bin_value = bin_value + small_bin.area
			sum_error2 = sum_error2 + small_bin.areaErr**2

	return bin_centres, bin_values, widths, errors

def plot(hist):
	bin_centres, bin_values, widths, errors = getBins(hist)
	bin_values = np.array(bin_values)
	norm = bin_values[0]
	bin_values = (bin_values/norm) * 27.45
	
	errors = np.array(errors)
	errors = (errors/norm) * 27.45

	print(bin_values)
	print(errors)
	plt.bar(bin_centres, bin_values, widths, align="center")
	plt.errorbar(bin_centres, bin_values, yerr=errors)
	plt.show()

def extractCoeff(hists):
	bin_values = [getBins(hist)[1] for hist in hists] #values of bins for each reweighting
	bin_values = np.array(bin_values)
	bin_values = bin_values/bin_values[0] #divide by SM value -> yields
	bin_values = bin_values - 1 #now left with INT and BSM terms

	no_params = len(pars)
	no_bins = len(bin_values[0])
	
	print(no_params)
	print(no_bins)

	#create a to carry all coefficients of for every bin
	coeffs = [ [[0.0,0.0] for i in range(no_params)] for j in range(no_bins)]
	coeffs = np.array(coeffs)
	
	for i in range(1, no_params+1): #for each parameter
		print("Solving for %s" % pars[i-1]['name'])
		step = pars[i-1]['step']/2
		for j in range(no_bins):
			#values of INT+BSM terms for bin in two steps of param
			mu1 = bin_values[i*2-1][j]
			mu2 = bin_values[i*2][j]
			
			A = (4*mu1-mu2)/(2*step)
			B = (mu2-2*mu1)/(2*step**2)
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
				param = pars[i]['name']
				A = bin[i][0]
				B = bin[i][1]
				#make exception for cG and cA parameters
				if param in ["cG","tcG", "cA"]:
					A = A / (4*np.pi)**2
					B = B / (4*np.pi)**4
				if A>0:
					string += " + %f * %s"%(A, param)
				else:
					string += " %f * %s"%(A, param)
				if B>0:
					string += " + %f * %s * %s"%(B, param, param)
				else:
					string += " %f * %s * %s"%(B, param, param)
			string += "\n"		
			file.write(string)

aos = yoda.read("Rivet.yoda", asdict = False)
hists = [h for h in aos if h.path.startswith("/HiggsTemplateCrossSectionsStage1/pT_Higgs")]
hists = hists[1:] #get rid of first histogram

#Bins for pT are 0,20,45,80,120,200
bin_edges = [0,20,45,80,120,200]
bin_names = ["GG2H_PTH_0_20", "GG2H_PTH_20_45", "GG2H_PTH_45_80", "GG2H_PTH_80_120", "GG2H_PTH_120_200"] 

with open("config.json") as jsonfile:
	pars = json.load(jsonfile)

coeffs = extractCoeff(hists)
writeTextFile(coeffs)
