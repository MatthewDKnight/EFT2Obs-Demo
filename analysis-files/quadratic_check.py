import matplotlib.pyplot as plt
import yoda
import numpy as np

aos=yoda.read("SimpleHiggs.yoda", asdict=False)

xs_list=[]
for i in aos:
	#if 'H_PT[rw' in i.name:
	if i.path.startswith("/SimpleHiggs/H_PT"):
		print i.name
		xs_list.append(i.sumW())
		
xs_list=xs_list[1:]


#cWW_values=[10**(i) for i in range(-4, 1)]
#cHW_values=[0., 0.5, 1.]
cWW_xs_values=xs_list
#cHW_xs_values=[xs_list[0], xs_list[3], xs_list[4]]
cWW_values = [0,0.0001,0.001,0.01,0.1,0.2,0.4,0.6,0.8,1]
cWW_xs_values = np.array(cWW_xs_values)

#cWW_xs_values[-1] += 10
#cWW_xs_values[-2] += -5
#cWW_xs_values[-3] += -2

error = 1e-2 * cWW_xs_values
print error

plt.errorbar(cWW_values, cWW_xs_values, error, fmt = 'o')
plt.xlabel('cWW')
plt.ylabel('Cross-section')
#plt.savefig('cWW_cross-section_plot.png')
#plt.show()
#plt.clf()

#plt.plot(cHW_values, cHW_xs_values)
#plt.xlabel('cHW')
#plt.ylabel('Cross-section')
#plt.savefig('cHW_cross-section_plot.png')

x = np.linspace(0,1,10000)

p = np.polyfit(cWW_values[-3:], cWW_xs_values[-3:], 2)
plt.plot(x, np.polyval(p,x))
plt.show()
