import matplotlib.pyplot as plt
import yoda
import numpy as np

aos=yoda.read("SimpleHiggs.yoda", asdict=False)

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

plt.xlabel('cWW')
plt.ylabel('Acceptance')

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

