import matplotlib.pyplot as plt
import yoda
import numpy as np
import sys
import scipy.optimize as spo

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

def quad(x, a, b):
	return 1 + a*x + b*x**2

aos=yoda.read(sys.argv[1], asdict=False)

#grabs all H_PT histograms (all the reweightings)
hists = []
for hist in aos:
	if hist.path.startswith("/SimpleHiggs/H_PT"):
		hists.append(hist)	
		
hists = hists[1:] #remove first histogram
hists = orderHists(hists) #orders them in terms of reweighting number

#get cross sections from each histogram and calculate the scaling, mu
xs = []
for hist in hists:
	xs.append(hist.sumW())
xs = np.array(xs)
mu = xs/xs[0]

#default param values to sample
param_vals = [0.1*i for i in range(11)]
param_vals = np.array(param_vals)

#find error in each cross section from total number of events and assume poisson
no_events = hists[0].numEntries()

error = 1/np.sqrt(no_events) * mu

#plot cross section points + a fitted quadratic
plt.errorbar(param_vals, mu, error, fmt = 'o')
plt.xlabel('cWW')
plt.ylabel('Scaling, mu(cWW)')

popt, pcov = spo.curve_fit(quad, param_vals, mu)
x = np.linspace(0,1,10000)
y = quad(x, popt[0], popt[1])

e = quad(param_vals, popt[0], popt[1])
chi2 = np.sum( (mu-e)**2/error**2 )
print chi2

plt.plot(x, y)
plt.savefig("quad_validate.png")
plt.show()
