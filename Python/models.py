# -*- coding: utf-8 -*-

from globals import *
from copy import deepcopy
from utils import *

##########################
## Read full text file  ##
##########################

#________________#
#   Tools        #
#________________#

#' Get fields in text file
#'
#' @param line character string
#' @param dic Python dictionary to store the result
#' @return x result (with conversion to float or int if possible)
def get_field(line, dic, convert="float"):
	x = line.split(" ")
	try:
		res = ifthenelse(convert=="float", float, int)(x[1])
	except:
		res = x[1]
	dic[x[0]] = res
	return(dic)

#' Get element from model object (list of items of a dictionary)
#'
#' @param tmp list of pairs (key, item)
#' @param name key of element to get
#' @return item item associated with key in tmp
get_from_model = lambda tmp, name : filter(lambda x : x[0]==name, tmp)[0][1]

#' Oneliner (for use in comprehension lists) to update dictionary
#' 
#' @param res dictionary
#' @param new_res new element
#' @param cmd key associated with new element
#' @return res updated dictionary
def modify_result(res, new_res, cmd):
	res[cmd] = new_res
	return(res)

#' Convert list of features describing a patterning function into
#' a dictionary, and add it to the pattern list
#'
#' @param patterns list of dictionaries that describe patterns
#' @param lst list of features associated with a patterning function
#' @param sep delimiter for each string element in lst between feature name and value
#' @param tofloat list of character strings (feature names) which value should be converted to float
#' @return patterns updated list of patterns
def todict(patterns, lst, sep=":", tofloat=["rate"]):
	## dictionary for the considered patterning function
	result = dict()
	## For each pair of (feature name, feature value)
	for x in [x.split(sep) for x in lst]:
		try:
			x[1] = ifthenelse(x[0] in tofloat, float, int)(x[1])
		except:
			## If the feature value can't be converted to another type
			## than string, add it to the dictionary as it is
			None
		result[x[0]] = x[1]
	return(patterns + [result])

#' Update dictionary describing model
#'
#' @param dictionary current dictionary describing model
#' @param params either to create a new element in the dictionary
#'               or to update an element
#' @return dictionary updated dictionary
def update_dictionary(dictionary, params):
	if (len(params) == 2):
		di = dict()
		## Creation of dictionary {typeE : element}
		di.setdefault(params[1], params[0])
		dictionary.update(di)
	else:
		if (dictionary.get(params[0])):
			di = update_dictionary(dictionary.get(params[0]), params[1:])
			dictionary.setdefault(params[0], di)
		else:
			dictionary.setdefault(params[0], update_dictionary(dict(), params[1:]))
	return(dictionary)

#' Dictionary traversal
#'
#' @param di dictionary
#' @param params list of keys
#' @return res element (or dictionary) returned after going through all keys
#'          in params
def call_dictionary(di, params):
	if (len(params)):
		di = filter(lambda x:x[0]==params[0], di.items())[0][1]
		return(call_dictionary(di, params[1:]))
	else:
		return(di)

#________________#
#   Model        #
#________________#

#' Parse model from text file
#'
#' @param f Python file object associated with text file
#' @param verbose logical for printing messages
#' @return res list that contains the description of the abstract model
#' [directives, patterns, genes, constants]
def get_model(f, verbose=False):
	res = { "patterns": [] }
	directives = dict()
	constants = dict()
	execute = lambda command, line : {
		## Keyword in the model file                  ## Keyword in the resulting dictionary
		"directive": lambda line : modify_result(res, 
			get_field(line, directives, convert="int"), "directives"),
		"fields": lambda line : modify_result(res, 
			[[int(get_element_before(x, "[")), map(float, get_element_between(x, "[", "]").split(","))] 
					for x in line.split(", ")], "fields"),
		"constant": lambda line : modify_result(res, 
			get_field(line, constants, convert="float"), "constants"),
		"genes": lambda line : modify_result(res,
			[[w[0], int(w[1])] for w in [[y.split("]")[0] for y in x.split("[")] 
				for x in line.split(";")[0].split("], ")]],
			"genes"),
		"pattern": lambda line : modify_result(res,
					todict(res["patterns"], line.split("\t"))
				, "patterns")
	}[command](line)
	## Parsing
	for line in f:
		command = line.split(" ")[0].split("\t")[0]
		line = line.split("\n")[0] + "\n"
		line = line.split(command + {"pattern": "\t"}.get(command, " "))[1:]
		line = concat(line[:-1] + [line[-1].split(";\n")[0]])
		res = execute(command, line)
	return(res.items())

#________________#
# Observations   #
#________________#

#' Parse conditions
#' 
#' @param x string associated with observations file
#' @return res list [Fixpoint, Phenotypes, GRNs, Summary]
def read_conditions(x, genes, verbose=False):
	Fixpoint = dict()
	Summary = dict()
	## Build summary            
	## - Summary is a dictionary of {[name]: {[step]: {[field]: "condition", "GRN"}}}
	summaries = filter(lambda x : grep(x, "|="), x)
	for y in summaries:
		name = get_element_between(y, "#", "[")
		step = get_element_between(y, "[", "]")
		field = get_element_between(y, "{", "}")
		element = get_element_after(y, "$")
		typeE = {True: "phenotype"}.get(grep(y, "Phenotypes{"), "GRN")
		params = [name, step, field, element, typeE]
		Summary = update_dictionary(Summary, params)
	## Build fixpoint conditions
	## - Fixpoint is a dictionary of {[name]: {[step]: {[field]: "condition", "GRN"}}}
	fixpoints = filter(lambda x : grep(x, "fixpoint("), x)
	for y in fixpoints:
		name = get_element_between(y, "#", "[")
		step = get_element_between(y, "[", "]")
		field = get_element_between(y, "{", "}")
		typeE = {True: "condition"}.get(grep(y, "Phenotypes{"), "GRN")
		params = [name, step, field, True, typeE]
		Fixpoint = update_dictionary(Fixpoint, params)
	x = filter(lambda x : not(grep(x, "#")), x)
	Phenotypes = dict()
	GRNs = dict()
	## Build phenotype & GRN dictionaries
	## - Phenotypes is a dictionary with keys: phenotype names, values: list of (node name, value)
	## in the associated condition
	m = len(x)
	i = 0
	while (i < m):
		## - GRNs is a dictionary with keys: phenotype names, values: GRFs in the associated condition
		if (grep(x[i], "GRNs") and grep(x[i], "$")):
			grn = dict()
			current_gene, current_formula, current_output = [None]*3
			name_GRN = get_element_between(x[i], "$", ":=")
			bracket_count = 1
			while (bracket_count > 0):
				i += 1
				bracket_count += count(x[i], "{")-count(x[i], "}")
				try:
					header = get_element_between(x[i], "$", ":=")
					if (header in genes):
						current_gene = header
					if (header == "formula"):
						i += 1
						current_formula = x[i]
					if (header == "output"):
						current_output = int(get_element_after(x[i], ":="))
					if (current_gene and current_formula and current_output):
						di = dict()
						di.setdefault("function_" + str(current_output), 
							{"formula": current_formula, 
							"output": current_output}
						)
						if (grn.get(current_gene)):
							di.update(grn.get(current_gene))
						grn[current_gene] = di
						current_formula, current_output = [None]*2
				except:
					continue
			GRNs.setdefault(name_GRN, grn)
		if (grep(x[i], "Phenotypes") and grep(x[i], "$")):
			phenotype_name = get_element_between(x[i], "$", ":=")
			j = i
			while (j < m and not (grep(x[j], "}"))):
				j += 1
			y = reduce(lambda a,b: a+" "+b, x[i:(j+1)]).split("{")[1].split("}")[0].split(" and ")
			y = [[concat(yij.split(" ")) for yij in yi.split("=")] for yi in y]
			Phenotypes.setdefault(phenotype_name, y)
			i = j
		i += 1
	if (verbose):
		print(all([len(y[1].items()) == len(genes) for y in GRNs.items()]))
	return([Fixpoint, Phenotypes, GRNs, Summary])

#' Parse observations file
#' 
#' @param f Python file object associated with observations file
#' @param genes set of genes/nodes
#' @param nsteps maximum number of steps
#' @param return_phenotypes logical for returning dictionary of cell phenotypes
#' @param verbose logical for printing messages
#' @return res list that describes the experiments
def get_observations(f, genes, nsteps, return_phenotypes=False, verbose=False):
	x = filter(lambda y:y and not(y[:2]=="//"), 
		[concat(y.split("\n")[0].split("\t")).split(";")[0] for y in f.readlines()])
	[Fixpoint, Phenotypes, GRNs, Summary] = read_conditions(x, genes)
	if (verbose):
		lst = [y + ": " for y in ["FIXPOINT", "PHENOTYPES", "GRNS", "SUMMARY"]]
		for i in range(len(lst)):
			print("\n" + lst[i] + str([Fixpoint, Phenotypes, GRNs, Summary][i]))
		print("")
	## One observation: {name, time step, field, corresponding GRN set, [phenotypes]}
	Observations = []
	for trajectory in call_dictionary(Summary, []).keys():
		for step in call_dictionary(Summary, [trajectory]).keys():
			for field in range(len(call_dictionary(Summary, [trajectory, step]).keys())):
				di = call_dictionary(Summary, [trajectory, step]).get(str(field+1))
				Observations.append({"name": trajectory, 
					"step": int(step), 
					"field": int(field), 
					## name of GRN instead of GRN itself: GRNs.get(di.get("GRN"))
					"GRN": di.get("GRN"),
					"phenotype": Phenotypes.get(di.get("phenotype"))
				})
	if (return_phenotypes):
		return([Observations, Fixpoint, GRNs, Phenotypes])
	return([Observations, Fixpoint, GRNs])

#________________#
#   Both files   #
#________________#

#' Parse both model and observations files (located in the same folder)
#'
#' @param model model file name
#' @param observations experimental file name
#' @param return_phenotypes logical for returning dictionary of cell phenotypes
#' @param verbose logical for printing messages
#' @return res an instance of the GRN inference problem
def read_files(model="model.net", observations="observations.spec", return_phenotypes=False, verbose=False):
	model, observations = [path_to_models + e for e in [model, observations]]
	## Extract model        ##
	with open(model, "r") as f:
		res = get_model(f, verbose=verbose)
	## Extract observations ##
	with open(observations, "r") as f:
		genes = [x[0] for x in get_from_model(res, "genes")]
		nsteps = get_from_model(res, "directives")["nsteps"]
		res_obs = get_observations(f, genes, nsteps, return_phenotypes=return_phenotypes, verbose=verbose)
	lst = ["genes", "patterns", "directives", "constants", "fields"]
	return([get_from_model(res, x) for x in lst] + res_obs)

