# -*- coding: utf-8 -*-

from z3 import *

from utils import *
from globals import aggregation_function, application_function, diffusion_function

#########################
# AGGREGATION FUNCTIONS #
#########################

## Dictionary of possible aggregation functions
aggregation = {
	## Sum patterning vectors: independant, uniform morphogen diffusion
	"logical_sum": lambda params : reduce(lambda a, b: a|b, params["arg"]),
	"arithmetic_sum": None
}

#' Aggregate (multiple) patterning vectors at a given time step
#'
#' @param params parameters for aggregation (time step, state, bit-vector, ...)
#' @param by which aggregation function should be used? (modifiable in globals.py)
#' @return res Z3 bit-vector that is the aggregation of patterning vectors
def aggregate_vectors(params, nf, by=aggregation_function):
	params.update({"nf": nf})
	if (params["type"] == "diffusion"):
		params.update({"arg": [params["vect_ls"][so][g][params["step"]][params["nf"]] 
			for so in range(params["nfields"]) for g in range(params["m"])]})
	elif (params["type"] == "patterning"):
		params.update({"arg": [params["vect_ls"][p][params["step"]][params["nf"]] 
			for p in range(params["npatt"])]})
	else:
		raise ValueError
	return(aggregation[by](params))

#########################
# APPLICATION FUNCTIONS #
#########################

## Dictionary of possible functions applying patterning effects to state variable
application = {
	## Sum environmental changes (morphogen diffusion) and current phenotype
	"logical_sum": lambda params : params["Q_previous"][params["nf"]] | params["arg"],
	"arithmetic_sum": None
}

#' Apply patterning result to system state at a given time step
#'
#' @param params parameters for aggregation (time step, state, bit-vector, ...)
#' @param by which application function should be used? (modifiable in globals.py)
#' @return res Z3 bit-vector that encodes the application of patterning changes to input system state
def apply_vectors(params, nf, by=application_function):
	params.update({"nf": nf})
	if (params["type"] in ["diffusion", "patterning"]):
		params.update({"arg": params["vect_ls"][params["step"]][params["nf"]]})
	else:
		raise ValueError
	return(application[by](params))

#########################
# DIFFUSION FUNCTIONS   #
#########################

## Distance between fields
dist = lambda source, field, power : sum([abs(source[i]-field[i])**power for i in range(len(field))])

def default_diffusion(params, power=1):
	d = dist(params["source"], params["field"], power)
	if (d == 0):
		coef = 1
	else:
		coef = params["rate"]/d
	return(int(params["max_value"]*coef))

def version1_diffusion(params, power=1):
	return(int(params["max_value"]*(params["rate"]**(dist(params["source"], params["field"], power)))))

## Dictionary of possible diffusion functions applied to state variable
diffusion = {
	"default": default_diffusion,
	"version1": version1_diffusion
}

#' (Discrete) diffusion of morphogen (as an nonincreasing staircase function)
#'
#' @param max_value maximum concentration level for the considered morphogen
#' @param source source field position of the morphogen diffusion source field
#' @param rate 0 <= rate <= 1 percentage of concentration loss -"default"- (or diffusion rate -"version1") from one field to another
#' @param field position of the considered field for which the morphogen concentration level should be computed
#' @return value concentration level (integer) of morphogen in field #nf
def apply_diffusion(params, field, by=diffusion_function):
	params.update({"field": field})
	return(int(diffusion[by](params)))

#########################
#########################

#' Precomputing gene levels in each field for a given source field, diffusion rate and axis
#' 
#' @param source integer index for source field in list "fields"
#' @param rate diffusion rate
#' @param fields lists of pairs field ID (integer), field position (list of float numbers)
#' @param idx_variable integer identifier for binary variable associated with highest diffusing conc. level
#' @param multi_binary_dict dictionary which contains the full mapping from binary
#'      variables to real genes and the way back
#' @param verbose logical for printing messages
#' @param debug logical for printing (hopefully) useful insights to the code
#' @return idx_lists lists of 
#'             ([field identifier, updated expression value, corresponding binary variable], 
#'               list of same-gene binary variable indices) for each gene
def precompute_levels(source, rate, fields, idx_variable, multi_binary_dict, verbose=False, debug=False):
	nfields = len(fields)
	idx_var, max_value = multi_binary_dict.get("map_to_vgenes")[idx_variable]
	ids = multi_binary_dict.get("map_list")[idx_var]
	genes = multi_binary_dict.get("genes")
	## Compute concentration levels of considered morphogen
	## in every field using this patterning function
	params = {"max_value": max_value, "rate": rate, "source": fields[source][1]}
	idxs = [apply_diffusion(params, field[1]) for field in fields]
	## List of [field no., level of gene product in the corresponding field, index of binary variable]
	## + ids: list of indices of binary variables associated with considered morphogen
	idx_list = [[[nf, idxs[nf], ids[max(0, idxs[nf]-1)]] for nf in range(nfields)], ids]
	if (verbose):
		print("-- " + multi_binary_dict.get("mapping").keys()[idx_var] + " on axis:" + str(axis+1) 
			+ " rate:" + str(rate) + "/fld from fld of index " + str(source))
	if (debug):
		print("- Max concentration lvl = " + str(max_value))
		for nf in range(nfields):
			print("Level in field #" + str(nf+1) + ": " 
				+ genes[idx_list[0][nf][-1]] + " = " 
				+ str(idx_list[0][nf][1]) + "/" 
				+ str(max_value))
	return(idx_list)

#########################
#########################

#' Precomputing gene product levels in each field for every associated binary variable
#' 
#' @param fields lists of pairs field ID (integer), field position (list of float numbers)
#' @param multi_binary_dict dictionary which contains the full mapping from binary
#'      variables to real genes and the way back
#' @param g binary variable integer identifier
#' @param source diffusion source field integer index
#' @param diffusion_rate value for diffusion rate
#' @param verbose logical for printing messages
#' @param debug logical for printing (hopefully) useful insights to the code
#' @return level_idx_lists lists of 
#'             ([gene identifier, updated expression value, corresponding binary variable], 
#'               list of same-gene binary variable indices) for each gene
def precompute_diffusion(fields, multi_binary_dict, g, source, diffusion_rate, verbose=False, debug=False):
	nfields = len(fields)
	diff_idx_list = precompute_levels(source, diffusion_rate, fields, g, 
		multi_binary_dict, verbose=verbose, debug=debug)
	return(diff_idx_list)

#' Precomputing morphogen levels in each field for every patterning function
#' 
#' @param fields lists of pairs field ID (integer), field position (list of float numbers)
#' @param multi_binary_dict dictionary which contains the full mapping from binary
#'      variables to real genes and the way back
#' @param pattern dictionary that describes the patterning function: features (keys) are
#'        "axis", "rate" (step function), "morphogen" (name), "source" (field ID)
#' @param verbose logical for printing messages
#' @param debug logical for printing (hopefully) useful insights to the code
#' @return level_idx_lists lists of 
#'             ([gene identifier, updated expression value, corresponding binary variable], 
#'               list of same-gene binary variable indices) for each gene
def precompute_pattern(fields, multi_binary_dict, pattern, verbose=False, debug=False):
	nfields = len(fields)
	source, rate = [x[0] for x in fields].index(int(pattern["source"])),  float(pattern["rate"])
	morphogen = pattern["morphogen"]
	idx_variable = multi_binary_dict.get("mapping")[morphogen]
	## Index of binary variable associated with highest conc. level
	g = multi_binary_dict.get("map_list")[idx_variable][-1]
	level_idx_list = precompute_levels(source, rate, fields, g, 
		multi_binary_dict, verbose=verbose, debug=debug)
	return(level_idx_list)

#########################
#########################

#' Computing diffusion (from patterning or regulatory action) vector
#' 
#' @param vect_ls list of Z3 bit-vectors of patterning/diffusion vectors (associated with 
#'            each patterning function/each binary variable, of size m = #binary variables)
#'            at each time step, in each field (independent from experiments)
#' @param nfields number of fields
#' @param Q_prev list of previous system states (bit-vectors of size m) of length nfields
#' @param idx_lists lists of 
#'             ([gene identifier, updated expression value, corresponding binary variable], 
#'               list of same-gene binary variable indices) for each gene
#' @param verbose logical for printing messages
#' @param debug logical for printing (hopefully) useful insights to the code
#' @return cond Z3 logical condition that encodes the computation of every patterning/gene product diffusion vector in vect_ls
def compute_diffusion(vect_bv, nfields, Q_prev, idx_list, verbose=False, debug=False):
	## What is the element concentration level in each field?
	if (debug):
		print("\n~~~ Binary variable indices: " + str(idx_list[1]))
		for nf in range(nfields):
			print("- FLD#" + str(nf+1) 
				+ " VAR#" + str(idx_list[0][nf][2]) 
				+ "=" + str(idx_list[0][nf][1]))
	cond_ls = [
			## vect_(field:nf)[idx] == 1 if (computed level > 0) else 0
			## where idx is the index of binary variable associated with element
			cond_value(vect_bv[nf],
				## In the list of computed levels, select the element
				## associated with field nf and get the identifier of binary variable
				## to update 
				idx_list[0][nf][2],
				## If the value to update is > 0 then return true bit-vector
				## else false bit-vector
				ifthenelse(idx_list[0][nf][1]>0, true, false)
			) for nf in range(nfields)
		]
	## Set other concentration levels (i.e. != idx) to 0 
	## (if there are not below the updated level)
	## otherwise to 1
	cond_ls += filter_nones([
		## If the index is not equal to the binary variable index to update
		ifthenelse(not(idx_list[0][nf][2]==idx),
		## For a given binary variable index idx
		## implement condition vect_(field:nf)[idx] == 0 
		## (if the associated binary variable is not active) 
		## or vect_(field:nf)[idx] == 1 
		## (if it is active: i.e. when the previous constraint was about 
		## a higher conc. level being active)
			cond_value(vect_bv[nf], 
				idx, 
				## If idx_list[0][nf][2] == higher level of variable
				## associated with index idx
				## and if this level was set to "active", then
				## binary variable associated with idx is set to 1
				## else it is set to 0
				ifthenelse(idx in idx_list[1] and idx_list[0][nf][2]>idx
					 and idx_list[0][nf][1]>0, 
					true, 
					false)
			)
		)
		## for every field, for every binary variable identifier
		## associated with the diffusing "gene"
		for nf in range(nfields) for idx in range(vect_bv[nf].size())
	])
	if (debug):
		print("\n--- Computation of conditions on patterning functions")
		for e in idx_list[0]:
			print("Var #" + str(e[2]) + "=" + str(e[1]) + " in field #" + str(e[0]))
		for i in range(len(cond_ls)):
			print(cond_ls[i])
	## Aggregate those conditions with a huge "AND" operator
	cond = map_and(cond_ls)
	return(cond)
