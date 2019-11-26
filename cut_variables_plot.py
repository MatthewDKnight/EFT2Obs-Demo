import yoda
import matplotlib.pyplot as plt
import numpy as np

aos = yoda.read("SimpleHiggs.yoda", asdict=False)

colours = {0:"red", 2:"blue", 4:"blue", 5:"blue"}

def plot(hist, norm = True):
	bin_values = []
	widths = []
	edges = []
	for bin in hist:
		bin_values.append(bin.height)
		widths.append(bin.xWidth)
		edges.append(bin.xMin)	

	if norm:
		bin_values = np.array(bin_values)
		bin_values = bin_values / np.sum(bin_values)

	number = findReweightNumber(hist)
	
	plt.bar(edges, bin_values, widths, align="edge",alpha=0.5, color = colours[number], label = number)

def plotDiff(hist1, hist2):
	bin_values1 = []
	bin_values2 = []
        centres = []
        for bin in hist1:
        	bin_values1.append(bin.height)
        	centres.append((bin.xMin+bin.xMax)/2)	
        
	for bin in hist2:
		bin_values2.append(bin.height)
                                                                                                   
        bin_values1 = np.array(bin_values1)
        bin_values1 = bin_values1 / np.sum(bin_values1)

	bin_values2 = np.array(bin_values2)
	bin_values2 = bin_values2 / np.sum(bin_values2)

        plt.scatter(centres, bin_values1-bin_values2)

def findUniqueNames(hists):
	names = []
	for hist in hists:
		name = hist.name.split("[")[0]
		if name not in names:
			names.append(name)
	return names

def findReweightNumber(hist):
	rw_number = hist.name.split("[")[1]
	number = rw_number[2]
        number = int(number)
	return number

exclude = ["H_PT", "N_jets", "lepton_no", "acceptance"]

hists = []

for hist in aos:
	if hist.path.startswith("/SimpleHiggs"):
                wanted = True
		for name in exclude:
			if name in hist.name:
				wanted = False
		if wanted:
			hists.append(hist)

names = findUniqueNames(hists)
hists = hists[len(names):] #remove weight merging histograms

def main():
        wanted_reweights = [0,4]
	for name in names:
		hists_plotted = []
		for hist in hists:
			if (name in hist.name) and (findReweightNumber(hist) in wanted_reweights):
				plot(hist)
				hists_plotted.append(hist)
		plt.title(name)
		plt.legend()
		#plt.savefig("yodaPlots/%s.png"%name)
		plt.show()
		plt.clf()

		plotDiff(hists_plotted[0], hists_plotted[1])
		plt.show()

