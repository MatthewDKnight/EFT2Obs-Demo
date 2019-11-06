print("Importing modules")
import yoda
import matplotlib.pyplot as plt
import numpy as np

print("Reading file")
aos = yoda.read("Rivet.yoda", asdict = False)

print("Rest")
hists = [h for h in aos if h.path.startswith("/HiggsTemplateCrossSectionsStage1/pT_Higgs")]

j_hists = [h for h in aos if h.path.startswith("/HiggsTemplateCrossSectionsStage1/Njets")]

#Bins for pT are 0,20,45,80,120,200

bin_edges = [0,20,45,80,120,200]

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
	
def plotAll():
	for hist in hists:
		bin_centres = []
		bin_values = []
		for bin in hist:
			edges = bin.xEdges
			bin_centres.append((edges[0]+edges[1])/2)
			bin_values.append(bin.area)

		plt.bar(bin_centres, bin_values,width=5,align='center')
		plt.show()

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

plot(hists[1])
