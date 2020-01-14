import sys

def splitRow(string):
	"""
	Will split a string into terms that start with '+' or '-'.
	"""
	pieces = []
	i = 0
	for j in range(len(string)):
		if string[j] == "+" or string[j] == "-":
			pieces.append(string[i:j])
			i = j
	pieces.append(string[i:].strip("\n"))
	return pieces


def decomposeTerm(term):
	"""
	Takes in a term like like '-0.682 * cT ' or '+ 10.8796 * cWW * cWW '
	and seperates it into the coefficient and the corresponding EFT
	parameter(s).
	"""
	term = term.replace("*", "") #remove *
	print(term)
	#find where the first parameter by finding first lowercase letter in string
	for i in range(len(term)):
		if ord(term[i])>97 and ord(term[i])<122: #use ascii
			first_half = term[:i] #contains coefficient
			second_half = term[i:] #contains term(s)
			break
	#for first_half, remove spaces and then float to give coefficient
	first_half = first_half.replace(" ", "")
	coefficient = float(first_half) 
	
	#make list of EFT parameters in term
	params = second_half.split(" ", 1)
	#print(params)
	if len(params) == 2:
		if params[1] == "": #if only one term
			del params[1]
		else:
			params[1] = params[1].replace(" ", "")

	#add exception for cG and cA
	#exceptions = ["cG", "cA"]
	#for param in params:
	#	if param in exceptions:
	#		coefficient = coefficient * (4*np.pi)**2

	return coefficient, params


prod_bin_names=[]
decay_bin_names=[]
prod_equations=[]
decay_equations = []
with open(sys.argv[1]) as prod:
	prod_content = prod.readlines()
	
for i in prod_content:
	prod_bin_names.append(i.split(':')[0])
	prod_equation = i.split(':')[1]
	prod_equations.append(prod_equation)

with open(sys.argv[2]) as decay:
	decay_content = decay.readlines()

for i in decay_content:
	decay_bin_names.append(i.split(':')[0])
	decay_equation = i.split(':')[1]
	decay_equations.append(decay_equation)

#print(splitRow('1-0.682 * cT + 10.86 * cWW * cWW')[1:])
prod_split_equations=[]
for i in range(len(prod_equations)):
	split_equation=[]
	for j in splitRow(prod_equations[i])[1:]:
		split_equation.append(decomposeTerm(j))
	prod_split_equations.append(split_equation)
#print(prod_split_equations)
		
decay_split_equations=[]
for i in range(len(decay_equations)):
	split_equation=[]
	for j in splitRow(decay_equations[i])[1:]:
		split_equation.append(decomposeTerm(j))
	decay_split_equations.append(split_equation)
#print(decay_split_equations)

params_list=[]
combined_split_equations=[]
for k in range(len(prod_split_equations)):
	params = []
	prod_equation=prod_split_equations[k]
	decay_equation=decay_split_equations[k]
	combined_equation=prod_equation+decay_equation
	prod_A_terms=[]
	decay_A_terms=[]
	for a in prod_equation:
		if len(a[1])==1:
			prod_A_terms.append(a)
	for b in decay_equation:
		if len(b[1])==1:
			decay_A_terms.append(b)
	for c in prod_A_terms:
		for d in decay_A_terms:
			new_cross_term=(c[0]*d[0], [c[1][0], d[1][0]])
			combined_equation.append(new_cross_term)
	combined_split_equations.append(combined_equation)

print(combined_split_equations)
				
	
with open(sys.argv[3], "w") as file:
		for i in range(len(combined_split_equations)):
			split_equation = combined_split_equations[i]

			string = "%s:1"%(prod_bin_names[i])

			for j in range(len(split_equation)):
				term = split_equation[j]

				coeff = term[0]
				params = term[1]
			
				if abs(coeff) > 10e-5:
					if len(params)==1:	
						if coeff>0:
							string += " + %f * %s"%(coeff, params[0])
						elif coeff<0:
							string += " %f * %s"%(coeff, params[0])
					if len(params)==2:
						if coeff>0:
							string += " + %f * %s * %s"%(coeff, params[0], params[1])
						elif coeff<0:
							string += " %f * %s * %s"%(coeff, params[0], params[1])	
			string += "\n"		
			file.write(string)


