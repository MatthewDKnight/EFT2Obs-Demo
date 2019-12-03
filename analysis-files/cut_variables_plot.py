import yoda
import matplotlib.pyplot as plt
import numpy as np
import matplotlib

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)

aos = yoda.read("SimpleHiggs_old.yoda", asdict=False)

colours = {0:"red", 2:"blue", 4:"green", 5:"blue"}
legend = {0:"SM", 1:"cWW", 2:"cWW", 3:"cHW", 4:"cHW"}

cut_points = {"leading_lep_PT":25,"trailing_lep_PT":13,
	      "leading_lep_eta":2.5,"trailing_lep_eta":2.5,
              "dilep_mass":12,"dilep_PT":30,"trailing_lep_mT":30,
	      "H_mT":60, "delta_phi":0, "cos_dphi":0}

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
	
	plt.bar(edges, bin_values, widths, align="edge",alpha=0.5, color = colours[number], label = legend[number])
	return np.max(bin_values)

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

	diff = -bin_values1+bin_values2
	
        plt.scatter(centres, diff)

	x = np.linspace(centres[0],centres[-1],2)
	y = [0,0]
	plt.plot(x,y,color="red")

	min = abs(np.min(diff))
	max = abs(np.max(diff))
	max = 1.2*np.max([min,max])
	plt.ylim(-max, max)
	
	return max

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
	for i in [2,4]:
		wanted_reweights = [0,i]
		for name in names:
			plt.figure(figsize=(10,7))
			hists_plotted = []
			for hist in hists:
				if (name in hist.name) and (findReweightNumber(hist) in wanted_reweights):
					max = plot(hist)
					hists_plotted.append(hist)
			
			x = [cut_points[name], cut_points[name]]
			y = [0, max]
			plt.plot(x,y,"--",color="black", linewidth=3)

			plt.xlabel(name)
			plt.ylabel("Probability")
			plt.legend()
			plt.savefig("cutVariablesPlots/%s%s.png"%(name,i))
			#plt.show()
			#plt.clf()
			
			plt.figure(figsize=(10,3))
			max = plotDiff(hists_plotted[0], hists_plotted[1])
			
			y = [-max, max]
			plt.plot(x,y,"--",color="black",linewidth=3)

			plt.xlabel(name)
			plt.ylabel("P(%s) - P(SM)"%legend[findReweightNumber(hists_plotted[1])])
			plt.savefig("cutVariablesPlots/%s%s_diff.png"%(name,i))
			#plt.show()

