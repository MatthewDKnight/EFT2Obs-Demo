import yoda
import matplotlib.pyplot as plt
import numpy as np

aos = yoda.read("SimpleHiggs.yoda", asdict=False)

def plot(hist):
	bin_values = []
	widths = []
	edges = []
	for bin in hist:
		bin_values.append(bin.height)
		widths.append(bin.xWidth)
		edges.append(bin.xMin)	
	
	plt.bar(edges, bin_values, widths, align="edge")
	plt.title(hist.name)
        plt.savefig("yodaPlots/%s.png"%hist.name)
	#plt.show()
        plt.clf()

for hist in aos:
	if hist.path.startswith("/SimpleHiggs"):
		plot(hist)
	if hist.path.startswith("/RAW/SimpleHiggs"):
		plot(hist)


