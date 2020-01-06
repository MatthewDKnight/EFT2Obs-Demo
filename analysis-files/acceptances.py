import matplotlib.pyplot as plt
import yoda
import numpy as np
import sys

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

aos=yoda.read(sys.argv[1], asdict=False)

no_acc_hists = []
acc_hists = []
for i in aos:
	if i.path.startswith("/SimpleHiggs/H_PT"):
		no_acc_hists.append(i)
	elif i.path.startswith("/SimpleHiggs/acc_H_PT"):
		acc_hists.append(i)
no_acc_hists = no_acc_hists[1:]
acc_hists = acc_hists[1:]
orderHists(no_acc_hists)
orderHists(acc_hists)	

pt_bin_edges = [0,20,45,80,120,200]
pt_bin_names = ["GG2H_PTH_0_20", "GG2H_PTH_20_45", "GG2H_PTH_45_80", "GG2H_PTH_80_120", "GG2H_PTH_120_200", "GG2H_PTH_GT200"]

njet_bin_edges = [0,1,2,3,4]
njet_bin_names = ["0j", "1j", "2j", "3j", ">3j"]

bounds = {"cww":[-0.015,0.015], "chw":[-0.012,0.016], "cphl":[-0.05,0.05],"Cross":[-1.0,1.0]}

"""
Plot the total acceptance as a function of the parameters to test for the form of the function
"""
SM_acc = acc_hists[0].sumW()/no_acc_hists[0].sumW()
acceptance_arrays = {}
params = ["cww","chw","cphl","Cross"]
for i in range(len(params)):
	param = params[i]
	
	no_acc_hists_temp = no_acc_hists[(1+(9*i)):(9*(i+1)+1)]
	acc_hists_temp = acc_hists[(1+(9*i)):(9*(i+1)+1)]

	no_acc_xs = []
	acc_xs = []
	for j in range(9):
		no_acc_xs.append(no_acc_hists_temp[j].sumW())
		acc_xs.append(acc_hists_temp[j].sumW())
	no_acc_xs = np.array(no_acc_xs)
	acc_xs = np.array(acc_xs)
	
	acceptance = acc_xs/no_acc_xs
	normed_acceptance = acceptance / SM_acc

	acceptance_arrays[param] = normed_acceptance

	bound = bounds[param]
	interval = bound[1] - bound[0]
	x = [(bound[0] + (interval/8)*i) for i in range(9)]

	plt.scatter(x, normed_acceptance)
	plt.show()

y1 = acceptance_arrays["cww"] + acceptance_arrays["chw"] + acceptance_arrays["cphl"] - 2
y2 = acceptance_arrays["Cross"]
plt.plot(x, y1, label="Total")
plt.plot(x, y2, label="Cross")
plt.legend()
plt.show()

"""
Determine the maximum acceptance for each bin.
"""

pt_bin_edges = [0,20,45,80,120,200]
pt_bin_names = ["GG2H_PTH_0_20", "GG2H_PTH_20_45", "GG2H_PTH_45_80", "GG2H_PTH_80_120", "GG2H_PTH_120_200", "GG2H_PTH_GT200"]
njet_bin_edges = [0,1,2,3,4]
njet_bin_names = ["0j", "1j", "2j", "3j", ">3j"]

bin_names = pt_bin_names + njet_bin_names

aos=yoda.read(sys.argv[1], asdict=False)

H_PT_hists = [h for h in aos if h.path.startswith("/SimpleHiggs/H_PT")]
H_PT_hists = H_PT_hists[1:] #get rid of first histogram
orderHists(H_PT_hists)

acc_H_PT_hists = [h for h in aos if h.path.startswith("/SimpleHiggs/acc_H_PT")]
acc_H_PT_hists = acc_H_PT_hists[1:] #get rid of first histogram
orderHists(acc_H_PT_hists)

N_jets_hists = [h for h in aos if h.path.startswith("/SimpleHiggs/N_jets")]
N_jets_hists = N_jets_hists[1:] #get rid of first histogram
orderHists(N_jets_hists)

acc_N_jets_hists = [h for h in aos if h.path.startswith("/SimpleHiggs/acc_N_jets")]
acc_N_jets_hists = acc_N_jets_hists[1:] #get rid of first histogram
orderHists(acc_N_jets_hists)

params = ["cww","chw","cphl"]
max_acceptances_holder = []
for k in range(len(params)):
        param = params[k]

        H_PT_hists_temp = H_PT_hists[(1+(9*k)):(9*(k+1)+1)]
        acc_H_PT_hists_temp = acc_H_PT_hists[(1+(9*k)):(9*(k+1)+1)]
	N_jets_hists_temp = N_jets_hists[(1+(9*k)):(9*(k+1)+1)]
        acc_N_jets_hists_temp = acc_N_jets_hists[(1+(9*k)):(9*(k+1)+1)]

	max_acceptances = []
	no_pt_bins = len(H_PT_hists_temp[0])
	no_njet_bins = len(N_jets_hists_temp[0])
	no_reweightings = 9

	#find maximum acceptance for each bin
	for i in range(no_pt_bins):
		SM_acceptance = acc_H_PT_hists[0][i].sumW/H_PT_hists[0][i].sumW
		acceptances = []
		for j in range(no_reweightings):
			acceptances.append(acc_H_PT_hists_temp[j][i].sumW/H_PT_hists_temp[j][i].sumW)
		max_acceptances.append(max(acceptances)/SM_acceptance)
	SM_acceptance = acc_H_PT_hists[0].overflow.sumW/H_PT_hists[0].overflow.sumW
        acceptances = []
        for j in range(no_reweightings):
                acceptances.append(acc_H_PT_hists_temp[j].overflow.sumW/H_PT_hists_temp[j].overflow.sumW)
        max_acceptances.append(max(acceptances)/SM_acceptance)

	for i in range(no_njet_bins):
		SM_acceptance = acc_N_jets_hists[0][i].sumW/N_jets_hists[0][i].sumW
		acceptances = []
		for j in range(no_reweightings):
			acceptances.append(acc_N_jets_hists_temp[j][i].sumW/N_jets_hists_temp[j][i].sumW)
		max_acceptances.append(max(acceptances)/SM_acceptance)
	SM_acceptance = acc_N_jets_hists[0].overflow.sumW/N_jets_hists[0].overflow.sumW
        acceptances = []
        for j in range(no_reweightings):
               acceptances.append(acc_N_jets_hists_temp[j].overflow.sumW/N_jets_hists_temp[j].overflow.sumW)
        max_acceptances.append(max(acceptances)/SM_acceptance)
	
	print("------------------------------------")
	print("Parameter: %s"%param)
	for i in range(len(bin_names)):
		print(bin_names[i], max_acceptances[i])

	max_acceptances_holder.append(np.array(max_acceptances))

#sum acceptance effect over all parameters
max_acceptances = max_acceptances_holder[0] + max_acceptances_holder[1] + max_acceptances_holder[2] - 2
print("---------------------------------")
print("Total")
for i in range(len(bin_names)):
	print(bin_names[i], max_acceptances[i])

with open("max_acceptances.txt", "w") as f:
	for i in range(len(bin_names)):
		string = bin_names[i] + ": %f"%max_acceptances[i] + "\n"
		f.write(string)
	
