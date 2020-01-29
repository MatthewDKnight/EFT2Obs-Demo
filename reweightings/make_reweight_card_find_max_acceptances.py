"""
Generates the reweight card needed for find_max_acceptance.py file. 
hard-coded to work only for cHW, cWW and cpHL. If you want it to work for more 
parameters then the bounds for these parameters must be added.
"""

import json
import sys

bounds = {"cww":[-15e-2, 15e-2],
	  "chw":[-12e-2, 16e-2],
	  "cphl":[-5e-2, 5e-2],
	  "cb":[-15e-2, 15e-2],
	  "ca":[-10e-4, 10e-4],
	  "chb":[-4.5e-2, 7.5e-2]}	  

def PrintBlock(pars, vals, label):
    res = []
    res.append('launch --rwgt_name=%s' % label)
    for par, val in zip(pars, vals):
        res.append('set %s %i %g' % (par['block'], par['index'], val))
    return res

with open(sys.argv[1]) as jsonfile:
    pars = json.load(jsonfile)

print pars

output = []

initvals = [0.] * len(pars)

current_i = 0
output.extend(PrintBlock(pars, initvals, 'rw%s' % current_i))
current_i += 1

#turn on each parameter separately, sampling across its bounds
for i in range(len(pars)):
	vals = list(initvals)
	param_name = pars[i]['name'].strip(" ")
	bound = bounds[param_name]
	step = (bound[1]-bound[0])/8
	for j in range(0,9):
		vals[i] = bound[0] + j*step
		output.extend(PrintBlock(pars, vals, 'rw%s' % current_i))
		current_i += 1

#turn on parameters all at once
for j in range(0,9):
	vals = list(initvals)
	for i in range(len(pars)):
		param_name = pars[i]['name'].strip(" ")
        	bound = bounds[param_name]
	        step = (bound[1]-bound[0])/8
		vals[i] = bound[0] + j*step
	output.extend(PrintBlock(pars, vals, 'rw%s' % current_i))
        current_i += 1

print '\n'.join(output)

if len(sys.argv) > 2:
    with open(sys.argv[2], 'w') as outfile:
            outfile.write('\n'.join(output))
