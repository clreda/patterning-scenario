# -*- coding: utf-8 -*-
 
import sys

from globals import *
from utils import *
from z3 import *
 
###################
## TOOLS         ##
###################

def print_error_message(error, msg=""):
	print("ERROR: " + error)
	if (msg):
		print("MSG: " + msg)

## s Z3 Solver
## bit-vector variables in vect_bv have same name length and same size
def enumerate_solutions(s, vect_bv, multi_binary_dict, nb=0, print_model=False):
	print("_____________________________________")
	print("\nMSG: Enumerate solutions from solver!")
	ss = str(s.check())
	i = 0
	while (ss == "sat" and ((nb > 0 and nb > i) or nb==0)):
		print("\n" + ss + " Solution #" + str(i+1))
		if (print_model):
			print("\n-------- MODEL --------")
			print(s.model())
			print("-----------------------\n")
		print(" "*(len(str(vect_bv[0]))+3) + concat(multi_binary_dict["genes"], " | "))
		none_ones, model_ones = [], []
		for x in vect_bv:
			if (s.model()[x]!=None): 
				print(str(x) + ":  " + concat(list(get_binary_dec(s.model()[x], vect_bv[0].size())), "   | "))
				model_ones += [x]
			else:
				print(str(x) + ":  Don't care!!")
				none_ones += [x]
		if (model_ones):
			s.add(map_or([s.model()[x] != x for x in model_ones]))
		if (none_ones):
			print("\nIs None = " + concat(none_ones, ", "))
			s.add(map_and([None != x for x in none_ones]))
		ss = str(s.check())
		i += 1
	if (nb > 0 and nb == i):
		print("\nMSG: Max. number of solutions reached.")
	else:
		print("\nMSG: No " + ifthenelse(i>0, "more ", "") + "solutions.")
	print("_____________________________________")
 
#######################
## MODELS.py test    ##
#######################
 
def test_models(model="model.net", observations="observations.spec"):
	model, observations = [path_to_models + tests_model + e for e in [model, observations]]
	from models import get_model, get_observations, get_from_model, read_files
	print("TEST START ---------")
	max_lines = 5
	if (False):
		print("\n** File:")
		with open(model, "r") as f:
			i = 1
			for line in f:
				if (i<max_lines):
					print("LINE #" + str(i) + ": " + line.split("\n")[0])
				i += 1
		print("** End of file")
	with open(model, "r") as f:
		tmp = get_model(f)
		if (False):
			print("\n** Result:")
			for x in tmp:
				print("-- " + x[0] + " --")
				if (len(x[1])>1):
					try:
						x[1].keys()
						for y in x[1].keys():
							print("- " + y + " " + str(x[1][y]))
					except:
						for y in x[1]:
							print(y)
				else:
					print(x)
			print("** End of result")
	with open(observations, "r") as f:
		genes = map(lambda x:x[0], get_from_model(tmp, "genes"))
		nsteps = get_from_model(tmp, "directives")["nsteps"]
		tmp = get_observations(f, genes, nsteps, verbose=False)
		if (False):
			print("\n** Result:")
			for x in tmp:
				if (len(x) > 1):
					for y in x:
						print(y)
				else:
					print(x)
			print("** End of result")
	tmp = read_files(model, observations, return_phenotypes=False, verbose=True)
	for x in tmp:
		print("")
		print(x)
		print("")
	print("\nTEST END -----------")

###################################
## FORMULA_TO_TRUTHTABLE.py test ##
###################################

def test_formula_to_tt(model="model.net", observations="observations.spec"):
	model, observations = [path_to_models + tests_model + e for e in [model, observations]]
	from formula_to_truthtable import enumerate_vectors, formula_to_tt
	from models import get_from_model, read_files
	print("TEST START ---------")
	print("\n** TESTING BOOLEAN VECTOR ENUMERATION **")
	boolean = True
	imax = 10
	v = 2
	for i in range(1, imax+1):
		args = ([range(v)]*i)
		res = enumerate_vectors(args)
		## Number of vectors in F^n where |F| = m is m^n
		## Length of each vector is m
		boolean = boolean and len(res[0]) == i and len(res) == v**i
		strf = [reduce(lambda x,y: x+y, map(str, outloop)) for outloop in res]
		## Generate unique vectors
		boolean = boolean and len(set(strf)) == v**i
	print(str(boolean) + " == True (boolean vector enumeration)")
	print("\n** TESTING CONVERSION FROM FORMULAS TO TRUTH TABLE **")
	print("\n--- Dummy example")
	verbose = True
	genes = [["A", 1], ["B", 2], ["C", 1]]
	grns = [{"A": {"function_1": {"output": 1, "formula": "( A=1 and B=2 )"}},
		"B": {"function_1": {"output": 1, "formula": "( A=1 )"}, 
			"function_2": {"output": 2, "formula": "( A=1 and C=1 )"}},
		"C": {"function_1": {"output": 1, "formula": "( A=1 )"}}}, 
		{"A": {"function_1": {"output": 1, "formula": "( A=1 and B=2 and C=1 )"}},
		"B": {"function_1": {"output": 1, "formula": "( C=1 )"}, 
			"function_2": {"output": 2, "formula": "( B=1 and C=1 )"}},
		"C": {"function_1": {"output": 1, "formula": "( A=1 )"}}}]
	GRNs = [[i, grns[i]] for i in range(len(grns))] 
	grns = formula_to_tt(GRNs, genes, verbose=verbose)
	print("\n--- Toy model")
	params = read_files(model, observations)
	## params = [genes, patterns, directives, constants, Observations, Fixpoint, GRNs]
	GRNs = params[-1]
	GRNs.setdefault(None, None)
	Observations = params[-3]
	genes = params[0]
	grns_list = [GRNs[x] for x in [oo["GRN"] for oo in Observations]]
	grns_list = [[i, grns_list[i]] for i in range(len(Observations))]
	tt_object = formula_to_tt(grns_list, genes, verbose=verbose)
	print("\nTEST END -----------")

##########################
## PATTERNS             ##
##########################

def aux_diffusion(rate, max_value, source, fields, f):
	import matplotlib.pyplot as plt
	from patterns import default_diffusion, dist
	nfields = len(fields)
	flatten = lambda x : [yy for y in x for yy in y]
	x = range(nfields+1)
	yy = [f({"source": fields[source][1], "field": fields[nf][1], "rate": rate, "max_value": max_value})
		for nf in range(nfields)]
	print(yy)
	y = flatten([[yyy, yyy] for yyy in yy])
	x = flatten([[x[j], x[j+1]] for j in range(nfields)])
	plt.axis([min(x), max(x), 0, max(y)+1])
	plt.plot(x, y, "b")
	plt.title("Diffusion function for rate=" + str(rate) + " max_value=" + str(max_value) 
		+ " from source " + str(source))
	plt.xlabel('Position (' 
			+ concat(["field #" + str(n+1) for n in range(nfields)], ", ") + ')')
	plt.ylabel('Concentration level thresholds')
	## Draw the separation between fields
	for nf in range(nfields+2):
		plt.axvline(x[nf], color="r", linestyle=":")
	## Draw the different concentration level values
	for vg in range(len(y)):
		plt.axhline(vg, color="g", linestyle=":")
	plt.show()
	plt.clf()

def test_diffusion():
	from patterns import default_diffusion, version1_diffusion
	print("TEST START ---------")
	nfields = 3
	## 2D
	fields = [[x, [float(x), min(0, float(x-1))]] for x in range(1, nfields+1)]
	## 1D
	fields = [[x, [float(x)]] for x in range(1, nfields+1)]
	source = 0
	max_value = 3
	f = version1_diffusion #default_diffusion 
	print("\nDiffusion of a product of max concentration level = " + str(max_value)
		+ " across " + str(nfields) + " fields from source " + str(source+1))
	print("\n--- Rate = 0 (no diffusion)")
	rate = 0 
	aux_diffusion(rate, max_value, source, fields, f)
	print("\n--- Rate = 0.5 (~50% diffusion)")
	rate = 0.5
	aux_diffusion(rate, max_value, source, fields, f)
	print("\n--- Rate = 0.75 (~75% diffusion)")
	rate = 0.75
	aux_diffusion(rate, max_value, source, fields, f)
	print("\n--- Rate = 1 (~100% diffusion)")
	rate = 1
	aux_diffusion(rate, max_value, source, fields, f)
	print("\nTEST END   ---------")

def test_diffusion2():
	## For Figure in report
	from patterns import default_diffusion, version1_diffusion
	print("TEST START ---------")
	nfields = 3
	fields = [[x, [float(x)]] for x in range(1, nfields+1)]
	source = 2
	max_value_M2, rate_M2 = 1, 1
	max_value_M1, rate_M1 = 3, 0.5
	f = version1_diffusion #default_diffusion 
	import matplotlib.pyplot as plt
	from patterns import default_diffusion, dist
	nfields = len(fields)
	flatten = lambda x : [yy for y in x for yy in y]
	x = range(nfields+1)
	yy2 = [f({"source": fields[source][1], "field": fields[nf][1], "rate": rate_M2, "max_value": max_value_M2})
		for nf in range(nfields)]
	yy1 = [f({"source": fields[source][1], "field": fields[nf][1], "rate": rate_M1, "max_value": max_value_M1})
		for nf in range(nfields)]
	y1 = flatten([[yyy, yyy] for yyy in yy1])
	y2 = flatten([[yyy, yyy] for yyy in yy2])
	x = flatten([[x[j], x[j+1]] for j in range(nfields)])
	plt.axis([min(x), max(x), 0, max(max(y1), max(y2))+1])
	plt.plot(x, y1, "b")
	plt.plot(x, y2, "r")
	plt.xlabel('Position (' 
			+ concat(["field #" + str(n+1) for n in range(nfields)], ", ") + ')')
	plt.ylabel('Concentration level thresholds')
	## Draw the separation between fields
	for nf in range(nfields+2):
		plt.axvline(x[nf], color="r", linestyle=":")
	## Draw the different concentration level values
	for vg in range(len(y1)):
		plt.axhline(vg, color="g", linestyle=":")
	for vg in range(len(y2)):
		plt.axhline(vg, color="g", linestyle=":")
	plt.fill_between(x, y2, facecolor='red', alpha=0.5)
	plt.fill_between(x, y1, facecolor='blue', alpha=0.5)
	plt.show()
	print("TEST END   ---------")

def test_precompute_levels():
	from patterns import precompute_levels
	print("TEST START ---------")
	genes = [["A", 3], ["B", 4], ["C", 1], ["D", 2]]
	print("Multi-level gene variables: ")
	print(concat([k[0] + " (" + str(k[1]) + " levels)" for k in genes], ", "))
	multi_binary_dict = multi_to_binary(genes)
	nfields = 6
	fields = [[x, [float(x)]] for x in range(1, nfields+1)]
	rate = 0.75
	idx_variable = 2 #gene A, third concentration level
	source = 0 #first field
	[ls, ids] = precompute_levels(source, rate, fields, idx_variable, multi_binary_dict)
	print("\nDiffusion of gene " + str(genes[idx_variable][0]) + " from Field " + str(source)
		+ " at rate=" + str(rate))
	for x in ls:
		print("In Field " + str(x[0]) + " " 
			+ str(multi_binary_dict["genes"][x[2]]) + " = " + str(x[1]))
	print("\nTEST END   ---------")

def test_compute_diffusion():
	from patterns import compute_diffusion, precompute_levels
	from build_conditions import condition_multi_level
	print("TEST START ---------")
	genes = [["A", 3], ["B", 4], ["C", 1], ["D", 2]]
	print("Multi-level gene variables: ")
	print(concat([k[0] + " (" + str(k[1]) + " levels)" for k in genes], ", "))
	multi_binary_dict = multi_to_binary(genes)
	m = len(multi_binary_dict["genes"])
	nfields = 4
	fields = [[x, [float(x)]] for x in range(1, nfields+1)]
	diffusing_genes = [
		## gene A, second concentration level diffusing from field #1
		{"rate": 0.75, "var": 1, "source": [x[0] for x in fields].index(1)},
		## gene B, third concentration level diffusing from field #4
		{"rate": 0.75, "var": 5, "source": [x[0] for x in fields].index(nfields)},
	]
	idx_lists = [precompute_levels(x["source"], x["rate"], 
			fields, x["var"], multi_binary_dict) for x in diffusing_genes]
	for i in range(len(diffusing_genes)):
		x = diffusing_genes[i]
		print("\nDiffusion of gene " + str(genes[multi_binary_dict["map_to_vgenes"][x["var"]][0]][0])
			+ " from Field " + str(x["source"])
			+ " at rate=" + str(x["rate"]))
		for y in idx_lists[i][0]:
			print("In Field " + str(y[0]) + " " 
				+ str(multi_binary_dict["genes"][y[2]]) + " = " + str(y[1]))
	print("")
	## Diffusion vectors 
	vect_bv = [[[BitVec('D_st[' + str(s) + '][' + str(i) + '][' + str(nf) + ']', m) for nf in range(nfields)] 
		for i in range(len(diffusing_genes))] for s in range(nfields)]
	## Aggregated diffusion vectors for each field
	vect_ls = [BitVec('D_st[' + str(nf) + ']' + ' '*6, m) for nf in range(nfields)]
	## Initial states, for each field
	Q_prev = [BitVec('QP_' + str(nf) + ' '*9, m) for nf in range(nfields)]
	## Final states, for each field
	Q_upd = [BitVec('QU_' + str(nf) + ' '*9, m) for nf in range(nfields)]
	s = Solver()
	####################################################################
	### Q_upd[nf] == Q_prev[nf] | D_step[nf] for all field nf
	### -----------------------------------------------
	### D_step[nf] is the vector associated to diffusion of products
	### in neighbouring fields of nf in state Q_prev
	### D_step[nf] = aggregation of all diffusion vectors in field nf
	### D_step[nf] == |_i,j D_step[source_i][diffusing_gene_j][nf] 
	### -----------------------------------------------
	### D_step[source_i][diffusing_gene_j][nf] is the vector of 
	### levels of diffusing_gene_j from field source_i in field nf
	####################################################################
	## Initialize states
	convert_to_binary = lambda lst : [y for x in [[1]*lst[i]+[0]*(genes[i][1]-lst[i]) for i in range(len(lst))] for y in x]
	generate_s = False
	if (generate_s):
		from random import randint
		states = [convert_to_binary([randint(0, x[1]) for x in genes]) for nf in range(nfields)]
	else:
		## 1 solution with no diffusion (Q_prev == Q_upd)
		states = [[0, 0, 1, 2]]*nfields
		states = [[0, 0, 1, 1]]*nfields
		states = [[0, 0, 0, 1]]*nfields
		states = [[1, 0, 0, 1]]*nfields
		## 1 solution with only var1 diffusion
		x1 = 1
		states = [[2, 0, 0, 0]]*x1 + [[0]*len(genes)]*(nfields-x1)
		x11 = 1
		states = [[1, 0, 0, 0]]*x11 + [[0]*len(genes)]*(nfields-x11)
		## 1 solution with only var2 diffusion
		x2 = 1
		states = [[0]*len(genes)]*(nfields-x2) + [[0, 3, 0, 0]]*x2
		## 1 solution with both diffusions
		x3 = 1
		states = [[2, 3, 0, 0]]*x3 + [[0]*len(genes)]*(nfields-2*x3) + [[2, 3, 0, 0]]*x3
		states = map(convert_to_binary, states)
	print(ifthenelse(generate_s, "Generated i", "I") + "nitial states:")
	print(" "*(len(str(Q_prev[nf]))+3) + concat(multi_binary_dict["genes"], " | "))
	for nf in range(nfields):
		print(str(Q_prev[nf]) + ":  " + concat(states[nf], "   | "))
		s.add(cond_equal_value(Q_prev[nf], states[nf]))
	print("")
	## Condition on multi-level variables
	s = condition_multi_level(s, multi_binary_dict, [[Q_prev, Q_upd]])
	## Simple
	if (False):
		## Let's suppose all diffusing genes described in diffusing_genes
		## diffuse at this step of the computation
		for i in range(len(diffusing_genes)):
			## Activate diffusion of the input diffusing genes
			## from the associated sources
			s.add(compute_diffusion(vect_bv[diffusing_genes[i]["source"]][i], 
				nfields, Q_prev, idx_lists[i], 
				verbose=True, debug=True))
	if (False):
		## Let's suppose all diffusing genes described in diffusing_genes
		## diffuse at this step of the computation
		for k in range(len(diffusing_genes)):
			for j in range(nfields):
				if (j == diffusing_genes[k]["source"]):
					print("D_st[" + str(j) + "][" + str(k) + "][nf] == diffusion of gene " + str(k)
						+ " from source " + str(j))
					s.add(compute_diffusion(vect_bv[diffusing_genes[k]["source"]][k], 
						nfields, Q_prev, idx_lists[k],	verbose=True, debug=False))
				else:
					print("D_st[" + str(j) + "][" + str(k) + "][nf] == 0")
					s.add(map_and([vect_bv[j][k][nf] == buildZERO(m) for nf in range(nfields)]))
	if (False):
		## Let's suppose all diffusing genes described in diffusing_genes
		## diffuse at this step of the computation
		for k in range(len(diffusing_genes)):
			for j in range(nfields):
				if (j == diffusing_genes[k]["source"]):
					print("D_st[" + str(j) + "][" + str(k) + "][nf] == diffusion of gene " + str(k)
						+ " from source " + str(j))
					s.add(compute_diffusion(vect_bv[diffusing_genes[k]["source"]][k], 
						nfields, Q_prev, idx_lists[k],	verbose=True, debug=False))
				else:
					print("D_st[" + str(j) + "][" + str(k) + "][nf] == 0")
					s.add(map_and([vect_bv[j][k][nf] == buildZERO(m) for nf in range(nfields)]))
		or_bit = lambda ls : reduce(lambda a,b : a|b, ls)
		f = lambda nf : [vect_bv[j][k][nf] for j in range(nfields) for k in range(len(diffusing_genes))]
		s.add(map_and([vect_ls[nf] == or_bit(f(nf)) for nf in range(nfields)]))
		## Apply this to the next states
		s.add(map_and([Q_upd[nf] == Q_prev[nf]|vect_ls[nf] for nf in range(nfields)]))
	if (True):
		## Apply the diffusion only when the associated gene products appear in previous states
		is_nactive = lambda k, bv : cond_value(Q_prev[diffusing_genes[k]["source"]], diffusing_genes[k]["var"], bv)
		is_diff = lambda j, k : Implies(is_nactive(k, true), compute_diffusion(vect_bv[diffusing_genes[k]["source"]][k], nfields, Q_prev, idx_lists[k]))
		is_not_diff = lambda j, k : Implies(is_nactive(k, false), map_and([vect_bv[j][k][nf] == buildZERO(m) for nf in range(nfields)]))
		diff = lambda j, k : And(is_diff(j, k), is_not_diff(j, k))
		print("Conditions:")
		for k in range(len(diffusing_genes)):
			j = diffusing_genes[k]["source"]
			print("[ Q_prev[nf][" + str(k) + "] == 1 => D_st[" + str(j) 
				+ "][" + str(k) + "][nf] == diffusion of gene " + str(k)
				+ " from source " + str(j) + "]")
			print("AND [ Q_prev[nf][" + str(k) + "] == 0 => D_st[" + str(j) 
				+ "][" + str(k) + "][nf] == 0 for all nf < " + str(nfields) + "]")
			cond = diff(j, k)
			#print("\n-------\n" + str(cond) + "\n*******\n")
			s.add(cond)
		for j in range(nfields):
			for k in range(len(diffusing_genes)):
				if (not(j == diffusing_genes[k]["source"])):
					print("D_st[" + str(j) + "][" + str(k) + "][nf] == 0 for all nf < " + str(nfields))
					cond = map_and([vect_bv[j][k][nf] == buildZERO(m) for nf in range(nfields)])
					#print(str(cond) + "\n**************************************\n")
					s.add(cond)
		or_bit = lambda ls : reduce(lambda a,b : a|b, ls)
		f = lambda nf : [vect_bv[j][k][nf] for j in range(nfields) for k in range(len(diffusing_genes))]
		print("D_st[0] == " + concat(f(0), "|"))
		s.add(map_and([vect_ls[nf] == or_bit(f(nf)) for nf in range(nfields)]))
		## Apply this to the next states
		s.add(map_and([Q_upd[nf] == Q_prev[nf]|vect_ls[nf] for nf in range(nfields)]))
	var = [z for x in vect_bv for y in x for z in y] + vect_ls	
	enumerate_solutions(s, var + Q_prev + Q_upd, multi_binary_dict, nb=3)
	print("\nTEST END   ---------")

##########################
## BUILD CONDITIONS     ##
##########################

def test_cond_observation():
	from build_conditions import condition_observation
	print("TEST START   ---------")
	genes = [["A", 3], ["B", 4], ["C", 1], ["D", 2]]
	print("Multi-level gene variables: ")
	print(concat([k[0] + " (" + str(k[1]) + " levels)" for k in genes], ", "))
	multi_binary_dict = multi_to_binary(genes)
	m = len(multi_binary_dict["genes"])
	nfields = 3
	k = 1
	s = Solver()
	exp_names = ["Trajectory"]
	Observations = [
		{
			"name": exp_names[0], "step": 0, "field": 0, "GRN": "InitialGRN", 
			"phenotype": [['A', '1'], ['B', '0'], ['C', '0']]
		},
		{
			"name": exp_names[0], "step": 0, "field": 1, "GRN": "InitialGRN", 
			"phenotype": [['A', '1'], ['B', '4'], ['C', '0']]
		},
		{
			"name": exp_names[0], "step": 0, "field": 2, "GRN": "InitialGRN", 
			"phenotype": [['A', '2'], ['B', '2'], ['C', '0']]
		}
	]
	states = [[[BitVec('q_f=' + str(nf) + "-t=" + str(t) + "-e=" + str(i), m) 
		for nf in range(nfields)] for t in range(k+1)] for i in range(len(exp_names))]
	s = condition_observation(s, Observations, exp_names, multi_binary_dict, 
		states, verbose=True, debug=True)
	x, y, z = states[0][0][0], states[0][0][1], states[0][0][2]
	s1, s2 = Solver(), Solver()
	print("\nDEBUG == TRUE")
	s1.add(And(And(And(And(And(And(Extract(1, 1, x) == 0,Extract(2, 2, x) == 0), Extract(0, 0, x) == 1), And(And(And(Extract(3, 3, x) == 0, Extract(4, 4, x) == 0), Extract(5, 5, x) == 0), Extract(6, 6, x) == 0)), Extract(7, 7, x) == 0), And(And(And(And(Extract(1, 1, y) == 0, Extract(2, 2, y) == 0), Extract(0, 0, y) == 1), And(And(And(Extract(3, 3, y) == 1, Extract(4, 4, y) == 1), Extract(5, 5, y) == 1), Extract(6, 6, y) == 1)), Extract(7, 7, y) == 0)), And(And(And(And(Extract(2, 2, z) == 0, Extract(0, 0, z) == 1), Extract(1, 1, z) == 1), And(And(And(Extract(5, 5, z) == 0, Extract(6, 6, z) == 0), Extract(3, 3, z) == 1), Extract(4, 4, z) == 1)), Extract(7, 7, z) == 0)))
	enumerate_solutions(s1, [x, y, z], multi_binary_dict)
	s2.add(And(And(And(And(And(And(And(And(And(And(Extract(1,1,x) ==0,Extract(2,2,x) ==0), Extract(0,0,x) == 1), And(And(And(Extract(3,3,x) ==0,Extract(4,4,x) ==0),Extract(5,5,x) ==0), Extract(6,6,x) == 0)), Extract(7, 7, x) == 0), And(And(Extract(1,1,y) == 0, Extract(2,2,y) == 0), Extract(0, 0, y) == 1)), And(And(And(Extract(3,3,y) == 1, Extract(4,4,y) == 1), Extract(5, 5, y) == 1), Extract(6, 6, y) == 1)), Extract(7, 7, y) == 0), And(And(Extract(2, 2, z) == 0, Extract(0, 0, z) == 1), Extract(1, 1, z) == 1)), And(And(And(Extract(5, 5, z) == 0, Extract(6, 6, z) == 0), Extract(3, 3, z) == 1), Extract(4, 4, z) == 1)), Extract(7, 7, z) == 0))
	print("\nDEBUG == FALSE")
	enumerate_solutions(s2, [x, y, z], multi_binary_dict)
	print("\nTEST END   ---------")

def test_cond_multi_level():
	from build_conditions import condition_multi_level
	from formula_to_truthtable import multi_to_boolean, dict_to_tt
	from utils import convert_grn
	print("TEST START   ---------")
	genes = [["A", 1], ["B", 2]]
	print("Multi-level gene variables: ")
	print(concat([k[0] + " (" + str(k[1]) + " levels)" for k in genes], ", "))
	multi_binary_dict = multi_to_binary(genes)
	m = len(multi_binary_dict["genes"])
	nfields = 1
	exp_names = ["Trajectory"]
	t = 0
	i = 0
	Q_updated = [BitVec('q_f=' + str(nf) + "-t=" + str(t) + "-e=" + str(i), m) for nf in range(nfields)]
	Q_next = [BitVec('q_f=' + str(nf) + "-t=" + str(t+1) + "-e=" + str(i), m) for nf in range(nfields)]
	states = [[Q_updated, Q_next] for e in range(len(exp_names))]
	s = Solver()
	## Add multi-level conditions
	s = condition_multi_level(s, multi_binary_dict, states)
	enumerate_solutions(s, Q_updated+Q_next, multi_binary_dict)
	nbsol = reduce(lambda x,y: x*y,[reduce(lambda x,y: x*y, [x[1]+1 for x in genes])]*len(states[0]))
	print("\n#Solutions = " + str(nbsol))
	print("\nTEST END   ---------")

def test_cond_transition():
	from build_conditions import condition_transition, condition_observation, condition_multi_level
	from build_conditions import condition_sync_transition, condition_async_transition
	from formula_to_truthtable import multi_to_boolean, dict_to_tt
	from utils import convert_grn
	print("TEST START   ---------")
	genes = [["A", 3], ["B", 4], ["C", 1], ["D", 2]]
	print("Multi-level gene variables: ")
	print(concat([k[0] + " (" + str(k[1]) + " levels)" for k in genes], ", "))
	multi_binary_dict = multi_to_binary(genes)
	m = len(multi_binary_dict["genes"])
	nfields = 1
	exp_names = ["Trajectory"]
	Observations = [
		{
			"name": exp_names[0], "step": 0, "field": 0, "GRN": "InitialGRN", 
			"phenotype": [['A', '3'], ['B', '2'], ['C', '1'], ['D', '1']]
		}
	]
	## "Monotonous" functions
	GRNs = {"InitialGRN": 
		{
			"A": {
				"function_1": {"output": 1, "formula": "( A=1 )"},
				"function_2": {"output": 2, "formula": "( A=1 and B=2 )"},
				"function_3": {"output": 3, "formula": "( A=1 and B=2 and C=1 )"},
			},
			"B": {
				"function_1": {"output": 1, "formula": "( A=1 )"}, 
				"function_2": {"output": 2, "formula": "( A=1 and C=1 )"},
				"function_3": {"output": 3, "formula": "( A=2 and C=1 )"},
				"function_4": {"output": 4, "formula": "( A=3 and C=1 )"},
			},
			"C": {"function_1": {"output": 1, "formula": "( A=1 )"}}, 
			"D": {
				"function_1": {"output": 1, "formula": "( A=1 )"},
				"function_2": {"output": 2, "formula": "( A=1 and B=1 )"},
			}
		}
	}
	t = 0
	i = 0
	Q_updated = [BitVec('q_f=' + str(nf) + "-t=" + str(t) + "-e=" + str(i), m) for nf in range(nfields)]
	Q_next = [BitVec('q_f=' + str(nf) + "-t=" + str(t+1) + "-e=" + str(i), m) for nf in range(nfields)]
	states = [[Q_updated, Q_next] for e in range(len(exp_names))]
	update_sync = False #True
	[vectors, GRNs, npossibilities] = dict_to_tt(GRNs, multi_binary_dict)
	bvectors = multi_to_boolean(vectors, m, multi_binary_dict)
	## Same constant GRN for each field
	grn_bv = convert_grn(GRNs, m, multi_binary_dict, npossibilities)*nfields
	do_test = map(bool, [1]*4)
	if (do_test[0]):
		s = Solver()
		## Add multi-level conditions
		s = condition_multi_level(s, multi_binary_dict, states)
		## Add initial state
		s = condition_observation(s, Observations, exp_names, multi_binary_dict, states)
		s = condition_transition(s, nfields, Q_next, Q_updated, update_sync, m, grn_bv, bvectors, verbose=False, debug=True)
		enumerate_solutions(s, Q_next+Q_updated, multi_binary_dict, nb=0)
		if (update_sync):
			print("1 expected solution.")
		else:
			nb_vectors = sum([x[1] for x in multi_binary_dict["vgenes"]])
			print("<= " + str(nb_vectors) 
				+ " expected solutions = #bvectors with one and only 1 distinct coordinate and "
				+ "respecting the multi-level constraint")
	if (do_test[1]):
		s1 = Solver()
		i = map(concat, vectors).index(concat([x[1] for x in Observations[0]["phenotype"]]))
		update = ifthenelse(update_sync, condition_sync_transition, condition_async_transition)
		s1 = condition_multi_level(s1, multi_binary_dict, states)
		s1 = condition_observation(s1, Observations, exp_names, multi_binary_dict, states)
		nf = 0
		hyp = cond_equal_value(Q_updated[nf], bvectors[i])
		ccl = update(Q_updated[nf], Q_next[nf], m, grn_bv[nf], i, verbose=True)
		s1.add(Implies(hyp, ccl))
		enumerate_solutions(s1, Q_next+Q_updated, multi_binary_dict, nb=0)
		if (update_sync):
			print("1 expected solution = " + str(states[0][t+1][nf]) + " == " 
				+ concat([simplify(extract1BV(grn_bv[nf][j], i)) for j in range(m)])
				+ ", " + str(states[0][t][nf]) + " == " + concat(bvectors[i]))
		else:
			nb_vectors = sum([x[1] for x in multi_binary_dict["vgenes"]])
			print("<= " + str(nb_vectors) 
				+ " expected solutions = #bvectors with one and only 1 distinct coordinate and "
				+ str(states[0][t+1][nf])
				+ ", " + str(states[0][t][nf]) + " == " + concat(bvectors[i]))
	if (do_test[2]):
		s2 = Solver()
		i = map(concat, vectors).index(concat([x[1] for x in Observations[0]["phenotype"]]))
		update = ifthenelse(update_sync, condition_sync_transition, condition_async_transition)
		s2 = condition_multi_level(s2, multi_binary_dict, states)
		s2 = condition_observation(s2, Observations, exp_names, multi_binary_dict, states)
		nf = 0
		hyp1 = cond_equal_value(Q_updated[nf], bvectors[i])
		ccl1 = update(Q_updated[nf], Q_next[nf], m, grn_bv[nf], i)
		i1 = (i+1)%len(bvectors)
		hyp2 = cond_equal_value(Q_updated[nf], bvectors[i1])
		ccl2 = update(Q_updated[nf], Q_next[nf], m, grn_bv[nf], i1)
		s2.add(And(Implies(hyp1, ccl1), Implies(hyp2, ccl2)))
		enumerate_solutions(s2, Q_next+Q_updated, multi_binary_dict, nb=0)
		if (update_sync):
			print("1 expected solution = " + str(states[0][t+1][nf]) + " == " 
				+ concat([simplify(extract1BV(grn_bv[nf][j], i)) for j in range(m)])
				+ ", " + str(states[0][t][nf]) + " == " + concat(bvectors[i]))
		else:
			nb_vectors = sum([x[1] for x in multi_binary_dict["vgenes"]])
			print("<= " + str(nb_vectors) 
				+ " expected solutions = #bvectors with one and only 1 distinct coordinate and "
				+ str(states[0][t+1][nf])
				+ ", " + str(states[0][t][nf]) + " == " + concat(bvectors[i]))
	if (do_test[3]):
		s3 = Solver()
		i = map(concat, vectors).index(concat([x[1] for x in Observations[0]["phenotype"]]))
		update = ifthenelse(update_sync, condition_sync_transition, condition_async_transition)
		s3 = condition_multi_level(s3, multi_binary_dict, states)
		s3 = condition_observation(s3, Observations, exp_names, multi_binary_dict, states)
		nf = 0
		hyp1 = cond_equal_value(Q_updated[nf], bvectors[i])
		ccl1 = update(Q_updated[nf], Q_next[nf], m, grn_bv[nf], i)
		i1 = (i+1)%len(bvectors)
		hyp2 = cond_equal_value(Q_updated[nf], bvectors[i1])
		ccl2 = update(Q_updated[nf], Q_next[nf], m, grn_bv[nf], i1)
		i2 = (i+2)%len(bvectors)
		hyp3 = cond_equal_value(Q_updated[nf], bvectors[i2])
		ccl3 = update(Q_updated[nf], Q_next[nf], m, grn_bv[nf], i2)
		s3.add(And(Implies(hyp1, ccl1), And(Implies(hyp2, ccl2), Implies(hyp3, ccl3))))
		enumerate_solutions(s3, Q_next+Q_updated, multi_binary_dict, nb=0)
		if (update_sync):
			print("1 expected solution = " + str(states[0][t+1][nf]) + " == " 
				+ concat([simplify(extract1BV(grn_bv[nf][j], i)) for j in range(m)])
				+ ", " + str(states[0][t][nf]) + " == " + concat(bvectors[i]))
		else:
			nb_vectors = sum([x[1] for x in multi_binary_dict["vgenes"]])
			print("<= " + str(nb_vectors) 
				+ " expected solutions = #bvectors with one and only 1 distinct coordinate and "
				+ str(states[0][t+1][nf])
				+ ", " + str(states[0][t][nf]) + " == " + concat(bvectors[i]))
	print("\nTEST END   ---------")

def test_cond_fixpoint():
	from build_conditions import condition_fixpoint
	print("TEST START   ---------")
	genes = [["A", 3], ["B", 4], ["C", 1], ["D", 2]]
	print("Multi-level gene variables: ")
	print(concat([k[0] + " (" + str(k[1]) + " levels)" for k in genes], ", "))
	multi_binary_dict = multi_to_binary(genes)
	m = len(multi_binary_dict["genes"])
	nfields = 2
	exp_names = ["Trajectory"]
	Observations = [
		{
			"name": exp_names[0], "step": 0, "field": 0, "GRN": "InitialGRN", 
			"phenotype": [['A', '3'], ['B', '2'], ['C', '1'], ['D', '1']]
		},
		{
			"name": exp_names[0], "step": 0, "field": 1, "GRN": "InitialGRN", 
			"phenotype": [['A', '2'], ['B', '1'], ['C', '1'], ['D', '1']]
		},
		{
			"name": exp_names[0], "step": 1, "field": 0, "GRN": "InitialGRN", 
			"phenotype": [['A', '0'], ['B', '1'], ['C', '1'], ['D', '1']]
		},
		{
			"name": exp_names[0], "step": 1, "field": 1, "GRN": "InitialGRN", 
			"phenotype": [['A', '0'], ['B', '0'], ['C', '0'], ['D', '0']]
		}
	]
	## experiment name : step : field {'condition', 'GRN'}
	Fixpoint = {'Trajectory': {'1': 
				{'0': {'condition': True},
				'1': {'condition': True}}
				}
			}
	t = 0
	i = 0
	k = 4
	states = [[[BitVec('q_f=' + str(nf) + "-t=" + str(t) + "-e=" + str(i), m) for nf in range(nfields)] 
			for t in range(k+1)]
		for e in range(len(exp_names))]
	print("\n-- With constant GRNs:")
	s = Solver()
	s = condition_fixpoint(s, Observations, exp_names, Fixpoint, states, k, m, grns_bv=None, 
			verbose=True, debug=True)
	print("\n-- With variable GRNs:")
	npossibilities = reduce(lambda x,y: x*y, [x[1]+1 for x in genes])
	grn_name = lambda e,t,nf,j : 'GRN^(' + str(e) + ',t=' + str(t) + ')_(' + str(nf) + ',' + str(j) + ')'
	grns_bv = [[[[
			BitVec(grn_name(e, t, nf, j), npossibilities) 
		for j in range(m)] 
		for nf in range(nfields) ] 
		for t in range(k) ]
		for e in range(len(exp_names)) ]
	## experiment name : step : field {'condition', 'GRN'}
	Fixpoint = {'Trajectory': {'1': 
				{'0': {'condition': True, 'GRN': True},
				'1': {'condition': True}}
				}
			}
	s = condition_fixpoint(s, Observations, exp_names, Fixpoint, states, k, m, grns_bv=grns_bv, 
			verbose=True, debug=True)
	print("\nTEST END   ---------")

def test_cond_utils():
	from build_conditions import reduce_bv, uone_op, morphogen_pattern_ids
	print("\nTEST START ---------")
	print("\n\nAUX reduce_bv")
	tests = ["11101", "1111000010010101010", "110100110010101"]
	for qval in tests:
		s = Solver()
		q = BitVec("q", len(qval))
		for i in range(len(qval)):
			s.add(cond_value(q, i, ifthenelse(int(qval[i]), true, false)))
		[s, bq] = reduce_bv(s, q, "bq")
		if (str(s.check()) == "sat"):
			lst = map(lambda x: get_binary_dec(s.model()[x], len(qval)), [q, bq])
			print("sat: " + "q = " + lst[0] + "; bq = " + lst[1])
	print("\n\nAUX uone_op")
	q = BitVec('q', 5)
	s = Solver()
	coefs_idx = range(3)
	print("Only one among the coefficient indices: " + concat(coefs_idx, ", "))
	coefs = [extract1BV(q, c) for c in coefs_idx]
	bv = uone_op(coefs)
	s.add(bv == true)
	ss = str(s.check())
	while (ss == "sat"):
		print("--")
		print("q = " + get_binary_dec(s.model()[q], 5))
		s.add(q != s.model()[q])
		ss = str(s.check())
	print("\n\nAUX morphogen_pattern_ids")
	patterns = [{'source': 1, 'rate': 1.0, 'morphogen': 'A'}, 
			{'source': 1, 'rate': 0.5, 'morphogen': 'A'},
			{'source': 4, 'rate': 1.0, 'morphogen': 'B'},
			{'source': 1, 'rate': 0.5, 'morphogen': 'B'},
			{'source': 1, 'rate': 1.0, 'morphogen': 'D'},
			{'source': 4, 'rate': 0.5, 'morphogen': 'D'}]
	npatt = len(patterns)
	tests = [["morphogen"], ["morphogen", "source"]]
	expected_results = map(str, [[[0, 1], [2, 3], [4, 5]], [[0, 1]]])
	for i in range(len(tests)):
		print("---")
		patterns_idx = morphogen_pattern_ids(patterns, npatt, eq_feature_list=tests[i])
		for pp in patterns:
			print(pp)
		print(str(patterns_idx) + " == " + expected_results[i])
	print("\nTEST END   ---------")

def test_cond_pattern_selection():
	from build_conditions import condition_pattern_selection
	print("\nTEST START ---------")
	genes = [["A", 3], ["B", 4], ["C", 1], ["D", 2]]
	print("Multi-level gene variables: ")
	print(concat([k[0] + " (" + str(k[1]) + " levels)" for k in genes], ", "))
	multi_binary_dict = multi_to_binary(genes)
	patterns = [{'source': 1, 'rate': 1.0, 'morphogen': 'A'}, 
			{'source': 1, 'rate': 0.5, 'morphogen': 'A'},
			{'source': 4, 'rate': 1.0, 'morphogen': 'B'},
			{'source': 1, 'rate': 0.5, 'morphogen': 'B'},
			{'source': 1, 'rate': 1.0, 'morphogen': 'D'},
			{'source': 4, 'rate': 0.5, 'morphogen': 'D'}]
	npatt = len(patterns)
	nsteps = 8
	patterning = [BitVec("patt_" + str(i), npatt) for i in range(nsteps)]
	## Values
	test_var = {"max_nb_patterns_per_level": 2, 
		"max_nb_patterns": 3, 
		"min_nb_patterns": 3, 
		"max_nb_pattern_times": 1}
	## Conditions
	s = Solver()
	s = condition_pattern_selection(s, patterning, npatt, nsteps, patterns, multi_binary_dict, test_var=test_var, verbose=False, debug=False)
	patterns_dict = {"genes": [concat([y[0]+":"+str(x[y]) for y in x.keys()], ",") for x in patterns]}
	print("\n-- Values:")
	for x in test_var.items():
		print(x[0] + " == " + str(x[1]))
	enumerate_solutions(s, patterning, patterns_dict, nb=2, print_model=False)
	print("\nTEST END   ---------")

def test_diff_model():
	from build_conditions import difference_model
	print("\nTEST START ---------")
	size = 3
	a = BitVec('x', size)
	bv_vars = {"a": [a]}
	varname = "a"
	s = Solver()
	s.add(UGT(a, 0))
	print("Condition: x != 0")
	ss = str(s.check())
	i = 0
	while(ss == "sat"):
		print("\n" + ss + ": MODEL = " + str(s.model()))
		s = difference_model(s, bv_vars, s.model(), varname, verbose=True, debug=True)
		i += 1
		ss = str(s.check())
	print("\n#expected solutions = " + str(2**size-1) + " == " + str(i))
	print("\nTEST END   ---------")

def test_cond_patterning():
	from build_conditions import condition_patterning, condition_multi_level
	from patterns import precompute_levels
	print("\nTEST START ---------")
	genes = [["A", 3], ["B", 4], ["C", 1], ["D", 2]]
	print("Multi-level gene variables: ")
	print(concat([k[0] + " (" + str(k[1]) + " levels)" for k in genes], ", "))
	multi_binary_dict = multi_to_binary(genes)
	m = len(multi_binary_dict["genes"])
	nfields = 3
	fields = [[x, [float(x)]] for x in range(1, nfields+1)]
	nsteps = 1
	patterns = [{'source': 1, 'rate': 1.0, 'morphogen': 'A'}, 
			{'source': 1, 'rate': 0.5, 'morphogen': 'A'},
			{'source': nfields, 'rate': 1.0, 'morphogen': 'B'},
			{'source': 1, 'rate': 0.5, 'morphogen': 'B'},
			{'source': 1, 'rate': 1.0, 'morphogen': 'D'},
			{'source': nfields, 'rate': 0.5, 'morphogen': 'D'}]
	npatt = len(patterns)
	idx_lists = [precompute_levels(patterns[i]["source"]-1, patterns[i]["rate"], fields, 
		multi_binary_dict["mapping"][patterns[i]["morphogen"]], multi_binary_dict) 
			for i in range(npatt)]
	print("\n--- Computation of idx_lists")
	for i in range(len(idx_lists)):
		print("\n-- Pattern #" + str(i+1))
		print("Morphogen " + patterns[i]["morphogen"]
			+ " -- Binary variables " 
			+ concat([multi_binary_dict["genes"][x] for x in idx_lists[i][1]], ", "))
		for nf in range(nfields):
			print("Field #" + str(idx_lists[i][0][nf][0]+1)
				+ " - Level " + str(idx_lists[i][0][nf][1])
				+ " - " + multi_binary_dict["genes"][idx_lists[i][0][nf][2]] + " = "
				+ str(int(idx_lists[i][0][nf][1] > 0)))
	vect_bv = [[BitVec('x_' + str(nf) + "t=" + str(t), m) 
			for nf in range(nfields)] 
		for t in range(nsteps)]
	vect_ls = [[[BitVec('y_' + str(nf) + "t=" + str(t) + "p=" + str(p), m) 
				for nf in range(nfields)] 
			for t in range(nsteps)] 
		for p in range(npatt)]
	Q_prev = [BitVec('QP_field=' + str(nf), m) for nf in range(nfields)]
	Q_upd = [BitVec('QU_field=' + str(nf), m) for nf in range(nfields)]
	patterning = [BitVec('Step^' + str(i), npatt) for i in range(nsteps+1)]
	s = Solver()
	s = condition_patterning(s, 0, nsteps, nfields, Q_prev, Q_upd, idx_lists, m, npatt, patterning, vect_ls, vect_bv, multi_binary_dict, verbose=True, debug=False)
	## Select only one pattern by step
	s.add(map_and([count_1s(x) == bv(0, npatt) for x in patterning]))
	## Multi-Level condition
	s = condition_multi_level(s, multi_binary_dict, [[Q_prev]])
	s = condition_multi_level(s, multi_binary_dict, [[Q_upd]])
	if (True):
		## Select pattern #p at each step
		p = npatt-1
		print("\n-- Select pattern #" + str(p+1) + ":")
		print(patterns[p])
		for x in patterning:
			s.add(x == idxList2BV([p], npatt))
	if (True):
		## No gene products at first
		print("\n-- Initial step: " + concat(Q_prev, ", ") + " == 0")
		for x in Q_prev:
			s.add(x == buildZERO(m))
	enumerate_solutions(s, Q_prev+Q_upd, multi_binary_dict, nb=5, print_model=False)
	print("\nTEST END   ---------")

def test_cond_diffusion():
	from build_conditions import condition_diffusion, condition_multi_level
	from patterns import precompute_diffusion
	print("\nTEST START ---------")
	genes = [["A", 3], ["B", 4], ["C", 1], ["D", 2]]
	print("Multi-level gene variables: ")
	print(concat([k[0] + " (" + str(k[1]) + " levels)" for k in genes], ", "))
	multi_binary_dict = multi_to_binary(genes)
	m = len(multi_binary_dict["genes"])
	nfields = 3
	fields = [[x, [float(x)]] for x in range(1, nfields+1)]
	nsteps = 1
	diffusion_rate = 0.5
	idx_lists = [[precompute_diffusion(fields, multi_binary_dict, g, nf, diffusion_rate) 
		for nf in range(nfields)] for g in range(m)]
	print("\n--- Computation of idx_lists")
	for g in range(m):
		print("\n-- Diffusion of " 
			+ multi_binary_dict["vgenes"][multi_binary_dict["map_to_vgenes"][g][0]][0]
			+ " (" + multi_binary_dict["genes"][g] + ")"
			+ "\n-- Binary variables " 
			+ concat([multi_binary_dict["genes"][x] for x in idx_lists[g][0][1]], ", ")
			+ " with diffusion rate=" + str(diffusion_rate))
		for nf in range(nfields):
			print("--- Source field #" + str(nf+1) + ":")
			for f in range(nfields):
				print("Field #" + str(idx_lists[g][nf][0][f][0]+1)
					+ " - Level " + str(idx_lists[g][nf][0][f][1])
					+ " - " 
					+ multi_binary_dict["genes"][idx_lists[g][nf][0][f][2]] 
					+ " = "
					+ str(int(idx_lists[g][nf][0][f][1] > 0)))
	vect_bv = [[BitVec('x_' + str(nf) + "t=" + str(t), m) 
			for nf in range(nfields)] 
		for t in range(nsteps)]
	vect_ls = [[[[BitVec('y_' + str(nf) + "t=" + str(t) + "g=" + str(g) + "s=" + str(f), m) 
				for nf in range(nfields)] 
			for t in range(nsteps)] 
		for g in range(m)] for f in range(nfields)]
	Q_diff = [BitVec('QP_field=' + str(nf), m) for nf in range(nfields)]
	Q_nxt = [BitVec('QU_field=' + str(nf), m) for nf in range(nfields)]
	s = Solver()
	for step in range(nsteps):
		s = condition_diffusion(s, step, nfields, Q_diff, Q_nxt, m, vect_ls,
			vect_bv, multi_binary_dict, idx_lists, verbose=True, debug=False)
	## Multi-Level condition
	s = condition_multi_level(s, multi_binary_dict, [[Q_diff]])
	s = condition_multi_level(s, multi_binary_dict, [[Q_nxt]])
	if (True):
		## Set initially diffusing products
		lst = [[0, 1, 2], [8, 9], [3, 4, 5, 6]]
		for ii in range(len(lst)):
			if (len(lst[ii])):
				print("\n-- Initial step: " + str(Q_diff[ii]) + ": only " 
				+ concat([multi_binary_dict["genes"][x] for x in lst[ii]], ", ") + " active")
				s.add(Q_diff[ii] == idxList2BV(lst[ii], m))
			else:
				print("\n-- Initial step: " + str(Q_diff[ii]) + " == 0")
				s.add(Q_diff[ii] == buildZERO(m))
	print("\nDiffusion rate = " + str(diffusion_rate))
	enumerate_solutions(s, Q_diff+Q_nxt, multi_binary_dict, nb=0, print_model=False)
	print("\nTEST END   ---------")

#############################
## PATTERN_SOLVER.py test  ##
#############################
 
def test_pattern_solver(model="model.net", observations="observations.spec"):
	model, observations = [path_to_models + tests_model + e for e in [model, observations]]
	from pattern_solver import pattern_solver
	from models import read_files
	print("TEST START ---------")
	params = read_files(model, observations)
	try:
		res = pattern_solver(params, verbose=True, debug=False)
	except ValueError:
		print("ERROR: Something occurred during execution")
	print("\nTEST END -----------")

##########################
## UTILS                ##
##########################

def test_one_selection():
	from build_conditions import one_selection
	from utils import map_and, get_binary_dec, count_1s
	print("TEST START ---------")
	n = 8
	l = 2
	hn = int(n/l)
	patterns_list = [range(hn)]*l
	for i in range(1, l):
		patterns_list[i] = [x+i*hn for x in patterns_list[i]]
	x = BitVec('x', n)
	s = Solver()
	cond_ls = [one_selection(x, np_list) for np_list in patterns_list]
	s.add(map_and(cond_ls))
	scheck = str(s.check())
	print("---")
	print(scheck)
	nsol = 0
	while (scheck == "sat"):
		M = s.model()
		print(M)
		## Enumeration of all possibilities
		s.add(x != M[x])
		print(get_binary_dec(M[x], n))
		for k in range(l):
			i, j = min(patterns_list[k]), max(patterns_list[k])
			vector = Extract(j, i, M[x])
			print("List #" + str(k+1) + ": " 
				+ str(simplify(count_1s(vector))) + " <= 1")
		scheck = str(s.check())
		print("---")
		print(scheck)
		nsol += 1
	## Number of solutions: per list (l lists), hn possible values + 1 value 
	## if no element from the list is selected
	print("Solution number: " + str(nsol) + " == " + str((hn+1)**l))
	print("\nTEST END -----------")

## for all index i, ls[i] == 0, 1 (integer list of size |q|)
def test_cond_value():
	from random import randint
	print("TEST START ---------")
	l = 10
	q = BitVec('q', l)
	ls = [randint(0, 1) for i in range(l)]
	print("Input = " + str(ls))
	vector = idxList2BV(filter_nones([ifthenelse(int(ls[i]), i) for i in range(len(ls))]), q.size())
	ls_vect = concat(map(str, ls))
	print("Vector associated with input: " + ls_vect)
	bv_vect = get_binary_dec(vector, q.size())
	print("Bit-Vector associated with input: " + bv_vect)
	print("Are the two vectors equal? " + str(ls_vect == bv_vect))
	cond_ls = [cond_value(q, i, ifthenelse(int(ls[i]), true, false)) for i in range(len(ls))]
	print("Constraints: " + str(cond_ls))
	print("Condition: " + str(map_and(cond_ls)))
	print("\nTEST END ---------")

def rev(bvv):
	return(concat(list(reversed(list(bvv)))))

def aux_binary_dec(nb, size, bdec):
	print("** test for n=" + str(nb))
	v = BitVecVal(nb, size)
	bd = rev(bdec)
	calc = sum([int(bd[i])*2**i for i in range(len(bd))])
	print(get_binary_dec(v, size) + " == \n" 
			+ rev(bdec) + "0"*(max(0, size-len(bdec))) + "_2 ")
	print(" == " + get_binary_dec(v) + "_10 == " + str(nb) + " -- " + str(calc) + " : " 
			+ str(nb == int(get_binary_dec(v))) + " == True")	

def test_utils_bv():
	print("------ START TEST")
	print(">>> Test get_binary_vec:")
	print("Order of print: 2^0 | 2^1 | ... | 2^i | 2^(i+1) | ...")
	aux_binary_dec(4329327034368, 48, "111111000000000000000000000000000000000000")
	aux_binary_dec(1152921503533105152, 60, "111111111111111111111111111111000000000000000000000000000000")
	print(get_binary_dec(BitVecVal(3,3), 3) + " == 110")
	print(get_binary_dec(BitVecVal(3,3), 5) + " == 11000")
	print(str(get_binary_dec(BitVecVal(3,3))) + " == 3")
	print(">>> Test bv:")
	print(get_binary_dec(bv(0, 3), 3) + " == 100")
	print(get_binary_dec(bv(1, 3), 3) + " == 010")
	print(get_binary_dec(bv(2, 3), 3) + " == 001")
	print(get_binary_dec(bv(3, 3), 3) + " == 000")
	print(">>> Test buildZERO:")
	print(get_binary_dec(buildZERO(3), 3) + " == 000")
	print(get_binary_dec(buildZERO(6), 6) + " == 000000")
	print(">>> Test buildONE:")
	print(get_binary_dec(buildONE(3), 3) + " == 111")
	print(get_binary_dec(buildONE(6), 6) + " == 111111")
	print(">>> Test trim1BV:")
	print(">> First try:")
	print(get_binary_dec(BitVecVal(15, 5), 5) + " == 11110")
	print(get_binary_dec(trim1BV(BitVecVal(15, 5), 0), 4) + " == 1110")
	print(get_binary_dec(trim1BV(BitVecVal(15, 5), 1), 4) + " == 1110")
	print(get_binary_dec(trim1BV(BitVecVal(15, 5), 2), 4) + " == 1110")
	print(get_binary_dec(trim1BV(BitVecVal(15, 5), 3), 4) + " == 1110")
	print(get_binary_dec(trim1BV(BitVecVal(15, 5), 4), 4) + " == 1111")
	print(">> Second try:")
	print(get_binary_dec(BitVecVal(12, 5), 5) + " == 01100")
	print(get_binary_dec(trim1BV(BitVecVal(12, 5), 0), 4) + " == 0110")
	print(get_binary_dec(trim1BV(BitVecVal(12, 5), 1), 4) + " == 0110")
	print(get_binary_dec(trim1BV(BitVecVal(12, 5), 2), 4) + " == 0010")
	print(get_binary_dec(trim1BV(BitVecVal(12, 5), 3), 4) + " == 0010")
	print(get_binary_dec(trim1BV(BitVecVal(12, 5), 4), 4) + " == 0011")
	print(">>> Test extract1BV:")
	print(str(simplify(extract1BV(BitVecVal(3, 3), 0))) + " == 1")
	print(str(simplify(extract1BV(BitVecVal(3, 3), 1))) + " == 1")
	print(str(simplify(extract1BV(BitVecVal(3, 3), 2))) + " == 0")
	print(">>> Test idxList2BV:")
	print(get_binary_dec(idxList2BV([2, 3], 4), 4) + " == 0011")
	print(get_binary_dec(idxList2BV([0], 4), 4) + " == 1000")
	print(get_binary_dec(idxList2BV([], 4), 4) + " == 0000")
	print(">>> Test cond_value:")
	print(str(simplify(cond_value(BitVecVal(3, 3), 0, true))) + " == True")
	print(str(simplify(cond_value(BitVecVal(3, 3), 0, false))) + " == False")
	print(">>> Test equal_except1BV:")
	print(get_binary_dec(BitVecVal(3, 3), 3))
	print(get_binary_dec(BitVecVal(1, 3), 3))
	print("Coordinate 0: ")
	print(str(simplify(equal_except1BV(BitVecVal(3, 3), BitVecVal(1, 3), 0))) + " == False")
	print("Coordinate 1: ")
	print(str(simplify(equal_except1BV(BitVecVal(3, 3), BitVecVal(1, 3), 1))) + " == True")
	print("Coordinate 2: ")
	print(str(simplify(equal_except1BV(BitVecVal(3, 3), BitVecVal(1, 3), 2))) + " == False")
	print(">>> Test cond_value_vec:")
	print("(1): " + get_binary_dec(BitVecVal(3, 3), 3))
	print("(2): " + get_binary_dec(BitVecVal(1, 3), 3))
	print("Coordinate 0: ")
	print(str(simplify(cond_value_vec(BitVecVal(3, 3), BitVecVal(1, 3), 0))) + " == True")
	print("Coordinate 1: ")
	print(str(simplify(cond_value_vec(BitVecVal(3, 3), BitVecVal(1, 3), 1))) + " == False")
	print("Coordinate 2: ")
	print(str(simplify(cond_value_vec(BitVecVal(3, 3), BitVecVal(1, 3), 2))) + " == True")
	print("Coordinate 0 (1) and 2 (2): ")
	print(str(simplify(cond_value_vec(BitVecVal(3, 3), BitVecVal(1, 3), 0, c1=2))) + " == False")
	print("Coordinate 1 (1) and 0 (2): ")
	print(str(simplify(cond_value_vec(BitVecVal(3, 3), BitVecVal(1, 3), 1, c1=0))) + " == True")
	print(">>> Test cond_equal_value:")
	values = [1, 1, 0]
	print("Values: " + str(values))
	print("One bit-vector variable x > 0 of size 3...")
	x = BitVec("x", 3)
	s = Solver()
	s.add(x > buildZERO(3))
	print("Conditions x[0] == Values[0], x[1] == Values[1], x[2] == Values[2]")
	s.add(cond_equal_value(x, values))
	## To see the conditions in native Z3 format       ##
	#print(s.sexpr())
	print(str(s.check()) + " == sat")
	M = s.model()
	print("x = " + get_binary_dec(M[x], 3) + " == " + concat(map(str, values)))
	print(">>> Test count_1s:")
	print(str(simplify(count_1s(BitVecVal(4329327034368, 48)))) + " == " 
		+ str(len(filter(lambda x:bool(int(x)), list("111111000000000000000000000000000000000000")))))
	print(str(simplify(count_1s(BitVecVal(1152921503533105152, 60)))) + " == " 
		+ str(len(filter(lambda x:bool(int(x)), list("111111111111111111111111111111000000000000000000000000000000")))))
	print(">>> Test get_element_between:")
	print(get_element_between("bl#ablal*ba", "#", "*") + " == ablal")
	print(">>> Test get_element_after:")
	print(get_element_after("bl#ablal*ba", "#") + " == ablal*ba")
	print(">>> Test grep:")
	print(str(grep("bl#ablal*ba", "#")) + " == True")
	print(str(grep("bl#ablal*ba", "c")) + " == False")
	print(">>> Test count:")
	print(str(count("bl#ablal*ba", "c")) + " == 0")
	print(str(count("bl#ablal*ba", "b")) + " == 3")
	print(str(count("bl#ablal*ba", "a")) + " == 3")
	print(str(count("bl#ablal*ba", "*")) + " == 1")
	print("------ END TEST")

def test_utils_multi_level():
	print("------ START TEST")
	genes = [["A", 3], ["B", 4], ["C", 1], ["D", 2]]
	print("Multi-level gene variables: ")
	print(concat([k[0] + " (" + str(k[1]) + " levels)" for k in genes], ", "))
	multi_binary_dict = multi_to_binary(genes)
	print("\n-- REAL GENE VARIABLES")
	print(concat([k[0] + " (" + str(k[1]) + " levels)" for k in multi_binary_dict["vgenes"]], ", "))
	print("\n-- BINARY VARIABLES")
	print(multi_binary_dict["genes"])
	print(str(len(multi_binary_dict["genes"])) + " == " + str(sum([k[1] for k in genes])))
	print("\n-- MAPPING BETWEEN REAL GENE AND BINARY VARIABLES")
	for k in range(len(multi_binary_dict["genes"])):
		x = multi_binary_dict["map_to_vgenes"][k]
		print(multi_binary_dict["genes"][k] + " b.v. of " + multi_binary_dict["vgenes"][x[0]][0]
			+ " - Level #" + str(x[1]) + "/" + str(multi_binary_dict["vgenes"][x[0]][1]))
	print("\n-- MAP_LIST")
	print("Integer identifier for each gene: " + str(multi_binary_dict["mapping"].items()))
	for k in range(len(multi_binary_dict["genes"])):
		id_gene, value_gene = multi_binary_dict["map_to_vgenes"][k]
		print(multi_binary_dict["vgenes"][id_gene][0] + "-" + str(value_gene) + " == " + multi_binary_dict["genes"][k])
		print(k in multi_binary_dict["map_list"][multi_binary_dict["mapping"][multi_binary_dict["vgenes"][id_gene][0]]])
	print("\n------------------------------")
	print("-- Conversion from multi-level vectors to Boolean ones (P. Van Ham's mapping)")
	vectors = [[0, 3, 1, 2], [0, 1, 1, 1], [2, 4, 0, 1], [3, 4, 1, 2]]
	m = len(multi_binary_dict["genes"])
	bvs = multi_to_boolean(vectors, m, multi_binary_dict)
	pgenes = concat([x[0] for x in genes], " | ")
	pbvs = concat(multi_binary_dict["genes"], " | ")
	msg = "Multi-Level = "
	for i in range(len(vectors)):
		print(" "*len(msg) + "-"*len(pgenes))
		print(" "*len(msg) + pgenes)
		print(msg + concat(vectors[i], " | "))
		print(" "*len(msg) + "-"*len(pbvs))
		print(" "*len(msg) + pbvs)
		print("Boolean =     " + concat(bvs[i], "   | "))
	print("------ END TEST")

def test_utils_convert_grn(model="model.net", observations="observations.spec"):
	from formula_to_truthtable import dict_to_tt
	print("------ START TEST")
	genes = [["A", 3], ["B", 4], ["C", 1], ["D", 2]]
	print("Multi-level gene variables: ")
	print(concat([k[0] + " (" + str(k[1]) + " levels)" for k in genes], ", "))
	multi_binary_dict = multi_to_binary(genes)
	m = len(multi_binary_dict["genes"])
	nfields = 1
	k = 1
	exp_names = ["Trajectory"]
	nexp = len(exp_names)
	Observations = [
		{
			"name": exp_names[0], "step": 0, "field": 0, "GRN": "InitialGRN", 
			"phenotype": [['A', '3'], ['B', '2'], ['C', '1'], ['D', '1']]
		}
	]
	## "Monotonous" functions
	GRNs = {"InitialGRN": 
		{
			"A": {
				"function_1": {"output": 1, "formula": "( A=1 )"},
				"function_2": {"output": 2, "formula": "( A=1 and B=2 )"},
				"function_3": {"output": 3, "formula": "( A=1 and B=2 and C=1 )"},
			},
			"B": {
				"function_1": {"output": 1, "formula": "( A=1 )"}, 
				"function_2": {"output": 2, "formula": "( A=1 and C=1 )"},
				"function_3": {"output": 3, "formula": "( A=2 and C=1 )"},
				"function_4": {"output": 4, "formula": "( A=3 and C=1 )"},
			},
			"C": {"function_1": {"output": 1, "formula": "( A=1 )"}}, 
			"D": {
				"function_1": {"output": 1, "formula": "( A=1 )"},
				"function_2": {"output": 2, "formula": "( A=1 and B=1 )"},
			}
		}
	}
	npossibilities1 = reduce(lambda x,y: x*y, [x[1]+1 for x in genes])
	s = Solver()
	print("\n-- Convert GRNs into Truth Tables")
	[vectors, GRNs, npossibilities] = dict_to_tt(GRNs, multi_binary_dict)
	print("Checking correct number of multi-level vectors = " 
		+ str(npossibilities) + " == " + str(npossibilities1) 
		+ ": " + str(npossibilities == npossibilities1))
	map_to_vgenes = multi_binary_dict["map_to_vgenes"]
	print("\n-- AUX: colGRN")
	idx_grn = 0 #first GRN
	j = 0 #gene A first level
	idx = map_to_vgenes[j][0]
	grn = col_grn(GRNs, map_to_vgenes, idx_grn, j)
	x = concat(map(lambda x: int(float(x)), list(grn)))
	y = concat(map(lambda x: int(float(x)), list(GRNs["InitialGRN"][0:, idx])))
	print("Variable " + multi_binary_dict["genes"][j] + " from gene " + multi_binary_dict["vgenes"][idx][0])
	print("Extracted GRN vector:\n" + x)
	print("GRN vector:\n" + y)
	print("Equality: " + str(x == y) + " == True")
	j = 8 #gene D first level
	idx = map_to_vgenes[j][0]
	grn = col_grn(GRNs, map_to_vgenes, idx_grn, j)
	x = concat(map(lambda x: int(float(x)), list(grn)))
	y = concat(map(lambda x: int(float(x)), list(GRNs["InitialGRN"][0:, idx])))
	print("Variable " + multi_binary_dict["genes"][j] + " from gene " + multi_binary_dict["vgenes"][idx][0])
	print("Extracted GRN vector:\n" + x)
	print("GRN vector:\n" + y)
	print("Equality: " + str(x == y) + " == True")
	print("\n-- AUX to_binary")
	print("Variable " + multi_binary_dict["genes"][j] + " from gene " + multi_binary_dict["vgenes"][idx][0])
	print("Multi-level GRN vector for gene " + multi_binary_dict["vgenes"][idx][0] + ":\n" 
		+ x)
	idx_list = [to_binary(map_to_vgenes, grn[i], j, i) for i in range(len(grn))]
	bgrn = [ifthenelse(i in idx_list, 1, 0) for i in range(npossibilities)]
	print("Boolean GRN vector for variable " + multi_binary_dict["genes"][j] + ":\n" + concat(bgrn))
	print("\n-- Convert GRNs (into variable Bit-Vectors)")
	[s, grns_bv] = convert_grn_bv(s, Observations, GRNs, m, k, nexp, nfields, multi_binary_dict, 
		exp_names, npossibilities, verbose=True, debug=True)
	flatten = lambda ls : [e for x in ls for y in x for z in y for e in z]
	nexp, t, field = [0]*3
	print("\nlength: " + str(all(map(lambda x:x.size()==npossibilities, flatten(grns_bv)))) + " == True")
	print("\n-- Convert GRNs (into constant Bit-Vectors)")
	grns_bv = convert_grn(GRNs, m, multi_binary_dict, npossibilities, verbose=True, debug=True)
	print("\n-- Bit-Vectors")
	for i in range(len(grns_bv)):
		print("\n---- GRN #" + str(i+1))
		for j in range(m):
			print("GRN vector associated with binary variable #" + str(j+1) + ": " 
				+ multi_binary_dict["genes"][j])
			print(get_binary_dec(grns_bv[i][j], npossibilities))
	print("\nlength: " + str(all(map(lambda x:x.size()==npossibilities, 
		[y for x in grns_bv for y in x]))) + " == True")
	print("------ END TEST")

def test_binary_to_multi():
	print("------ START TEST")
	from numpy import transpose, zeros
	from random import randint
	genes = [["A", 3], ["B", 4], ["C", 1], ["D", 2]]
	possible_values = [[0, g[1]] for g in genes]
	print("Multi-level gene variables: ")
	print(concat([k[0] + " (" + str(k[1]) + " levels)" for k in genes], ", "))
	multi_binary_dict = multi_to_binary(genes)
	nfields = 1
	k = 5
	m = len(multi_binary_dict["genes"])
	nlin = (k+1)*nfields
	trajectory = zeros((nlin, m))
	## Populate the matrix
	lines = [[randint(possible_values[i][0], possible_values[i][1]) for i in range(len(genes))] for j in range(nlin)]
	print("\n -- Generated multi-level columns")
	print(" "*len("Step # : ") + concat([g[0]+" "*max(0, 2-len(g[0])) for g in genes], "| "))
	for j in range(nlin):
		print("Step #" + str(j) + ": " + concat(lines[j], " | "))
	boolean_lines = multi_to_boolean(lines, m, multi_binary_dict)
	print("\n-- Conversion:")
	for i in range(nlin):
		print(concat(lines[i]) + " => " + concat(boolean_lines[i]))
	for j in range(nlin):
		trajectory[j, 0:] = boolean_lines[j]
	print("\n--- Binary trajectory matrix")
	print(concat([g+" "*max(0, 1-len(g)) for g in multi_binary_dict["genes"]], " | "))
	for i in range(nlin):
		print(concat(map(int, trajectory[i, ]), "   | "))
	trajectory_multi = binary_to_multi(trajectory, k, multi_binary_dict)
	print("\n--- Multi-level trajectory matrix")
	print(concat([g[0]+" "*max(0, 2-len(g[0])) for g in genes], " | "))
	for i in range(nlin):
		print(concat(map(int, trajectory_multi[i, ]), "  | "))
	print("------ END TEST")
 
##########################
## VISUALIZATION        ##
##########################

def test_to_igraph():
	from visualize import model_to_igraph
	import numpy as np
	print("TEST START ---------")
	patterns = [{'source': 1, 'rate': 1.0, 'morphogen': 'Bcd'}
		,{'source': 1, 'rate': 0.5, 'morphogen': 'Hb-mat'}
		,{'source': 4, 'rate': 1.0, 'morphogen': 'Cad'}
		,{'source': 1, 'rate': 0.5, 'morphogen': 'Bcd'}
		,{'source': 1, 'rate': 1.0, 'morphogen': 'Hb-mat'}
		,{'source': 4, 'rate': 0.5, 'morphogen': 'Cad'}
		,{'source': 4, 'rate': 1.0, 'morphogen': 'Bcd'}
		,{'source': 4, 'rate': 0.5, 'morphogen': 'Hb-mat'}
		,{'source': 1, 'rate': 1.0, 'morphogen': 'Cad'}
		,{'source': 4, 'rate': 0.5, 'morphogen': 'Bcd'}
		,{'source': 4, 'rate': 1.0, 'morphogen': 'Hb-mat'}
		,{'source': 1, 'rate': 0.5, 'morphogen': 'Cad'}]
	modelname = "gap-gene-sanchez"
	i = 1
	filename = "test_files/result_" + modelname + "_" + str(i+1) + "_pattern_matrix.csv"
	pattern_matrix = np.loadtxt(open(filename, "rb"), delimiter=",", skiprows=1)
	print("\n* Pattern matrix:")
	print("Patterns v")
	if (len(np.shape(pattern_matrix)) > 1):
		print("t= " + concat(map(str, range(np.shape(pattern_matrix)[1])), " | "))
	else:
		print("t = 0")
	print(pattern_matrix)
	model_to_igraph(pattern_matrix, i, modelname, patterns, verbose=True)
	print("\nConversion to igraph done!")
	print("\nTEST END ---------")

def test_to_plot():
	from visualize import model_to_plot
	import numpy as np
	print("TEST START ---------")
	multi_binary_dict = {
		'map_list': [[0, 1, 2], [3], [4, 5, 6], [7, 8, 9], [10, 11], [12, 13], [14]], 
		'vgenes': [['Bcd', 3], ['Gt', 1], ['Hb-zyg', 3], ['Hb-mat', 3], ['Kr', 2], ['Cad', 2], ['Kni', 1]], 
		'genes': ['Bcd-1', 'Bcd-2', 'Bcd-3', 'Gt-1', 'Hb-zyg-1', 'Hb-zyg-2', 
			'Hb-zyg-3', 'Hb-mat-1', 'Hb-mat-2', 'Hb-mat-3', 'Kr-1', 'Kr-2', 'Cad-1', 'Cad-2', 'Kni-1'], 
		'map_to_vgenes': [[0, 1], [0, 2], [0, 3], [1, 1], [2, 1], [2, 2], [2, 3], [3, 1], 
			[3, 2], [3, 3], [4, 1], [4, 2], [5, 1], [5, 2], [6, 1]], 
		'mapping': {'Hb-zyg': 2, 'Gt': 1, 'Kni': 6, 'Bcd': 0, 'Kr': 4, 'Hb-mat': 3, 'Cad': 5}
	}
	k = 8
	nexp = 1
	nfields = 4
	modelname = concat("gap gene sanchez".split(" "), "-")
	i = 1
	filename = "test_files/result_" + modelname + "_" + str(i+1) + "_trajectories.csv"
	trajectory = np.loadtxt(open(filename, "rb"), delimiter=",", skiprows=1)
	nrow, ncol = np.shape(trajectory)
	print("#fields=" + str(nfields) + " #exp=" + str(nexp) + " #steps=" + str(k) + " #bvar=" + str(ncol))
	print("#rows=" + str(nrow) + ", #fields x (#steps+1) = " + str(nfields*(k+1)) 
		+ ": True == " + str(nrow == nfields*(k+1)))
	if (False):
		for idx_i in range(1, nrow+1):
			print(trajectory[idx_i-1, 0:])
			if (idx_i%nfields==0):
				print("")
	model_to_plot(trajectory, k, nexp, i, modelname, multi_binary_dict, verbose=True, debug=True)
	print("\nConversion to plot(s) done!")
	print("\nTEST END ---------")

##########################
## CALL                 ##
##########################
 
x = filter_nones([ifthenelse(x[0:5] == "test_", x) for x in globals().keys()])
tests = [globals()[y] for y in x]
names = [n[5:] for n in x]
i = 0
lenargvC = len(sys.argv) == 2

print("\nTESTING COMMANDS --------------")
print("Available tests for functions: " + str(x))
print("Names: " + str(names))
print("#provided arguments: " + ifthenelse(lenargvC, str(2), "Wrong number"))
print("-------------------------------\n")

## Run all tests in a row
def all_tests():
	print("\n-- START OF ALL TESTS")
	symb = "*"
	for f in x:
		ws = len(str(f))+11
		deco = "\n" + symb*ws + "\n"
		print(deco + symb + " TITLE: " + f + " " + symb + deco)
		eval(f)()
	print("\n-- END OF ALL TESTS")

done = False
while (i < len(names) and lenargvC):
	if (names[i] == sys.argv[1]):
		tests[i]()
		done = True
		break
	i += 1
if (lenargvC and sys.argv[1] == "--all"):
	all_tests()
	done = True
if (lenargvC and sys.argv[1] == "--h" and not done):
	done = True
if (i == len(names) and lenargvC and not done):
	error = "If you wanted to run a test on a file, maybe you mistyped the filename."
	msg = "You may want to use the following command: \'python test.py " + str(names) + "\'."
	print_error_message(error, msg)
