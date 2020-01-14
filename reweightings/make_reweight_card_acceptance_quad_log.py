"""
Generates the reweight card needed for acceptance_quad_log.py which
checks whether the acceptance follows a quadratic form.
For each parameter c_j, there will be 11 reweightings
with c_j = 0,0.1,0.2,...,1, and all other parameters set to 0.
"""

import json
import sys

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

#for i in xrange(len(pars)):
#    vals = list(initvals)
#    for j in range(5):
#       vals[i] = pars[i]['step']*10**j
#       output.extend(PrintBlock(pars, vals, 'rw%s' % current_i))
#       current_i += 1

for i in range(len(pars)):
	vals = list(initvals)
	for j in range(1,11):
		vals[i] = j*0.1
		output.extend(PrintBlock(pars, vals, 'rw%s' % current_i))
		current_i += 1

print '\n'.join(output)

if len(sys.argv) > 2:
    with open(sys.argv[2], 'w') as outfile:
            outfile.write('\n'.join(output))
