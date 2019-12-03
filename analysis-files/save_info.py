#import check_param_card as param_card_mod

CARDS_DIR="cards/ggF/"

import os
import sys

sys.path.append(os.path.join(os.environ['MG_DIR'], "ggF", 'bin', 'internal'))
import check_param_card as param_card_mod

param_card = param_card_mod.ParamCard(CARDS_DIR+"param_card.dat")
#par = param_card[args.block].param_dict[(i,)]

info_text = ""

with open(CARDS_DIR+"run_card.dat", "r") as f:
        for row in f:
                if "nevents" in row:
                        number = row.split("=")[0]
                        number = int(number)
                        break
        info_text += "Number of events: %d\n"%number


info_text += "\nReweighting\n"
info_text += "---------------------------\n"

with open(CARDS_DIR+"reweight_card.dat", "r") as f:
	rw = 0
	for row in f:
		if "launch" in row:
			info_text += "rw%d \n"%rw
			rw +=1
		elif "set" in row:
			numbers = row.split("newcoup")[1]
			numbers = numbers.split(" ")
			
			index = numbers[1]
			value = numbers[2].split("\n")[0]
			
			index = int(index)
			value = float(value)
			
			par = param_card["newcoup"].param_dict[(index,)]
			info_text += "%s = %f \n"%(par.comment, value)			
		

info_text += "\nProcesses\n"
info_text += "---------------------------\n"

with open(CARDS_DIR+"proc_card.dat", "r") as f:
	for row in f:
		if row[0]=="#":
			continue
		elif "generate" in row:
			info_text += row.split("generate")[1]
		elif "add process" in row:
			info_text += row.split("add process")[1]

with open("%s.txt"%sys.argv[1], "w") as f:
	f.write(info_text)

print(info_text)
