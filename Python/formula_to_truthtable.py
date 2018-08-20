# -*- coding: utf-8 -*-

import numpy as np
import itertools as it

from utils import *
from globals import *

#############################
## TOOLS                   ##
#############################

#' Enumerate all vectors in args[0] x args[1] x ... x args[|args|]
#'
#' @param args list of value lists
#' @return res list of vectors (list of elements in each list of args)
def enumerate_vectors(args):
	return([list(x) for x in list(it.product(*args))])
        
#############################
## FORMULA -> TRUTH TABLE  ##
#############################

## GRN is formatted as follows (see models.py): 
## {"gene1": {"function_1": {"formula": [infix expression as character string with parentheses and white spaces],
##                           "output": [integer level of expression when associated formula is evaluated to True]}
##		}, ... }
## e.g. {'A': {'function_1': {'formula': '( ( B=1 or C=0 ) and D=0 )', 'output': 1}}}

#' Convert formulas (extracted from the observations file) to truth tables
#'
#' @param grn_list list of pairs (integer identifier, GRN (dictionary of formulas))
#' @param genes list of pairs (gene name, #levels)
#' @param verbose logical for printing messages
#' @return res res = [enumeration of possible vector values, 
#' 			gene ordering, 
#'			dictionary of (GRN identifier, truth table)]
def formula_to_tt(grn_list, genes, verbose=False):
	grns = dict()
	## Same keys for every GRN
	## (it is needed to fix the numbering)
	keys = [g[0] for g in genes]
	m = len(keys)
	## Get the maximum concentration level for each gene
	args = [range(g[1]+1) for g in genes]
	## Get the number of possible values for each gene
	nvalues = sum([len(a) for a in args])
	if (verbose):
		print("Genes = " + str(keys))
		print("Values = " + str(args))
		print("#values = " + str(nvalues))
	## Numbering of (multi-level, integer) vectors
	vectors = enumerate_vectors(args)
	if (verbose):
		print("Vectors:")
		print([reduce(lambda x,y:x+y, [str(j) for j in v]) for v in vectors])
	## For every GRN
	for [idf, grn] in grn_list:
		if (verbose):
			print("\n-- GRN #" + str(idf+1))
		if (not grn):
			print("GRN is None!")
			grns.setdefault(idf, None)
			continue
		## Get all formulas
		items = [grn[k] for k in keys]
		## Initialization of truth table
		## of size #possible vectors x #genes 
		tt = np.zeros((len(vectors), m))
		for idx_output in range(len(keys)):
			## Get list of formulas for this GRN
			functions = [item[1] for item in items[idx_output].items()]
			for function in functions:
				## Get the formula (character string)
				f = function["formula"].split(" ")
				## Replace equality litteral by their Boolean value
				## (according to integer vector #v_idx)
				## e.g. if equality litteral = "Gene=i" and that
				## vectors[v_idx][Gene] = i, 
				## then replace equality litteral by "True"
				## otherwise if vectors[v_idx][Gene] != i
				## then replace equality litteral by "False"	
				for v_idx in range(len(vectors)):
					v = vectors[v_idx]
					v = [keys[j]+"="+str(v[j]) for j in range(len(v))]
					fv = [ifthenelse(grep(f[i], "="), str((f[i] in v)), f[i])
					 for i in range(len(f))]
					## After replacement by Boolean values,
					## evaluate the formula to the corresponding Boolean value
					res = eval(concat(fv, sep=" "))
					## If the formula is satisfied with respect to vectors[v_idx]
					## replace the zero value in the truth table by the output
					## concentration level (by the maximum concentration level)
					if (res):
						## Maximum level of gene expression ##
						tt[v_idx, idx_output] = max(function["output"], tt[v_idx, idx_output])
		if (verbose):
			print("** Formulas")
			for gf in grn.items():
				print(gf)
			print("** Truth Table")
			print(concat([g[0]+" "*(2-len(g[0])) for g in genes] + ["Vector"], " | "))
			vs = [reduce(lambda x,y:x+y, [str(j) for j in v]) for v in vectors]
			tt_lines = np.concatenate((tt, np.reshape(vs, (len(vs), 1))), 1)
			for i in range(np.shape(tt_lines)[0]):
				print(concat(map(lambda x: int(float(x)), tt_lines[i, :-1]) 
					+ [tt_lines[i, -1]], "  | "))
		grns.setdefault(idf, tt)
	return([vectors, keys, grns])

#' Enumerate all Boolean vectors corresponding to all possible values
#' of binary variables and turn the dictionaries of GRNs into Truth Tables
#' 
#' @param GRNs dictionary of GRNs (as described in the observations file)
#' @param multi_binary_dict dictionary which contains the full mapping from binary
#'      variables to real genes and the way back
#' @param verbose logical for printing messages
#' @param debug logical for printing (hopefully) useful insights to the code
#' @return grn_vector build (variable) Z3 Bit-Vectors associated with GRN of every gene
def dict_to_tt(GRNs, multi_binary_dict, verbose=False, debug=False):
	import numpy as np
	if (any([len(x[1]) < len(multi_binary_dict["vgenes"]) for x in GRNs.items()])):
		lst = [len(x[1]) < len(multi_binary_dict["vgenes"]) for x in GRNs.items()]
		grn_name = GRNs.keys()[lst.index(True)]
		length = len(GRNs[grn_name])
		print("ERROR: Not all GRFs provided in GRN \'" + grn_name 
			+ "\': #GRFs = " + str(length) + " < " + str(len(multi_binary_dict["vgenes"])))
		raise ValueError
	tt_obj = formula_to_tt(GRNs.items(), multi_binary_dict.get("vgenes"))
	## [x[0] for x in multi_binary_dict['vgenes']] = tt_obj[1]
	vectors = tt_obj[0]
	vs = [reduce(lambda x,y:x+y, [str(j) for j in v]) for v in vectors]
	npossibilities = reduce(lambda x,y : x*y, [(gi[1]+1) for gi in multi_binary_dict.get("vgenes")])
	GRNs_tt = dict()
	for x in tt_obj[2].items():
		tt = np.concatenate((x[1], np.reshape(vs, (len(vs), 1))), 1)
		GRNs_tt.setdefault(x[0], tt)
		if (debug):
			print("\n--- Dimensions of the truth table")
			print("#rows = " 
				+ str(tt.shape[0]) + " == " + str(npossibilities))
			print("#cols = " 
				+ str(tt.shape[1]-1) + " == " 
				+ str(len(multi_binary_dict.get("vgenes"))))			
	if (debug and verbose):
		for j in range(len(GRNs_tt)):
			[idf, tt] = GRNs_tt.items()[j]
			print("\nGRN \'" + str(idf) + "\'")
			for x in GRNs[idf].items():
				print(x)
			print(concat([gi[0]+" "*max(0, 4-len(gi[0])) for gi in multi_binary_dict.get("vgenes")]))
			for i in range(tt.shape[0]):
				print(concat(map(lambda x: int(float(x)), tt[i, :-1]) + [tt[i, -1]], " | "))
	return([vectors, GRNs_tt, npossibilities])
