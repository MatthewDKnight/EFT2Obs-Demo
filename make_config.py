import argparse
import os
import sys
import json
import parameter_conversion

#get arguments
parser = argparse.ArgumentParser()
parser.add_argument('--process', '-p', default='ggF')
parser.add_argument('--output', '-o', default='config.json')
parser.add_argument('--block', default='newcoup')
parser.add_argument('--default-value', type=float, default=None)
parser.add_argument('--limit-pars', type=str, default=None)
parser.add_argument('--default-value-inactive', type=float, default=None)
parser.add_argument('--default-value-file', type=str, default=None)
args = parser.parse_args()

#check to see if default values are given
if args.default_value is None and args.default_value_file is None:
	raise Exception("Must provide a default value for all POIs or a file containing the default values")

#get param_card
sys.path.append(os.path.join(os.environ['MG_DIR'], args.process, 'bin', 'internal'))
import check_param_card as param_card_mod
param_card = param_card_mod.ParamCard('%s/%s/Cards/param_card.dat' % (os.environ['MG_DIR'], args.process))
# then to modify an entry of such object you can do
ids = [X[0] for X in param_card[args.block].keys()]
print ids

#find ids of POIs
limit_pars = list(ids)
if args.limit_pars is not None:
    limit_pars = [int(X) for X in args.limit_pars.split(',')]

#read in default values if file is supplied
filename = args.default_value_file
default_values = {}
if filename != None:
	with open(filename, "r") as f:
		for row in f:
			param, default_value = row.split(":")
			default_values[param] = float(default_value)

#set default values of all parameters (POIs and not POIs)
for i in ids:
    par = param_card[args.block].param_dict[(i,)]
    if args.default_value is not None and i in limit_pars:
        par.value = args.default_value
    if args.default_value_file is not None and i in limit_pars:
	param_name = par.comment.strip(" ")
	par.value = default_values[parameter_conversion.params[param_name]]
    if i not in limit_pars and args.default_value_inactive is not None:
        par.value = args.default_value_inactive

param_card.write('param_card.dat')

output = []

for i in limit_pars:
    if i not in limit_pars:
        continue
    par = param_card[args.block].param_dict[(i,)]
    output.append({
        "name": par.comment,
        "block": args.block,
        "index": i,
        "step": par.value
        })

with open(args.output, 'w') as outfile:
        outfile.write(json.dumps(output, sort_keys=True, indent=4))
# param_card[block].param_dict[lhaid].value = value
