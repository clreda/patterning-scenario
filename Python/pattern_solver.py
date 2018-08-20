# -*- coding: utf-8 -*-

from time import time 
from z3 import *
import numpy as np

from formula_to_truthtable import dict_to_tt
from utils import *
from build_conditions import *
from patterns import precompute_pattern, precompute_diffusion
from globals import change_grn

##########################
## SOLUTION PROCESSING  ##
##########################

#' Generate a model from a solution returned by the solver
#' 
#' @param s Z3 solver object
#' @param nexp number of trajectories
#' @param m number of binary variables
#' @param npatt number of possible patterning functions
#' @param nfields number of fields
#' @param k maximum length of trajectories
#' @param nsteps maximum number of selected patterns
#' @param bv_vars dictionary of bit-vectors describing the (potential) solutions
#' @param level_idx_lists lists of 
#'             ([gene identifier, updated expression value, corresponding binary variable], 
#'               list of same-gene binary variable indices) for each gene
#' @param multi_binary_dict dictionary which contains the full mapping from binary
#'      variables to real genes and the way back
#' @param bvectors enumeration of Boolean vectors for the binary variables in the model
#' @param res_list list of solution models so far
#' @param pattern_matrix0 if not None, allows predictions using the patterning matrix 
#'        (of selected patterning functions) provided
#' @param patterns dictionary of patterning functions
#' @param uniqueness character string for solution uniqueness
#' @param verbose logical for printing messages
#' @param debug logical for printing (hopefully) useful insights to the code
#' @return res [updated solution list, updated solver with uniqueness condition, 
#'             logical indicating if the solver should stop looking for conditions]
def get_solutions(s, nexp, m, npatt, nfields, k, nsteps, bv_vars, level_idx_lists, multi_binary_dict, bvectors, patterns, uniqueness, res_list=[], pattern_matrix0=None, verbose=False, debug=False):
	no_more_model = False
	M = None
	print("\n-------------------------\n-- The instance is " + str(s.check()) + "istiable")
	try:
		M = s.model()
	except:
		if (res_list):
			print("\nNo other model found.")
		else:
			print("\nNo model found.")
		no_more_model = True
		return([res_list, s, no_more_model])
	if (M):
		print("\nMODEL #" + str(len(res_list)+1) + " :")
		if (pattern_matrix0 == None):
			## COMPUTATION OF THE PATTERN MATRIX (patterning function selection)
			## pattern_matrix[p, t] = is patterning function #p applied
			## at time step t? == 1 iff. bit-vector Pattern^(t=t,u)_p > 0
			## iff. patterning[t]_p == true
			pattern_matrix = np.zeros((npatt, nsteps))
			for t in range(nsteps):
				vect = get_binary_dec(M[bv_vars.get("patterning")[t]], npatt)
				for p in range(npatt):
					pattern_matrix[p, t] = int(vect[p])
		else:
			pattern_matrix = pattern_matrix0
		if (debug):
			print("\n* Pattern matrix:")
			print("t= " + concat([str(i)+" "*(3-len(str(i))+1) for i in range(nsteps)]))
			print(pattern_matrix)
		if (verbose):
			print("\n** Summary")
			for t in range(nsteps):
				selected = filter(lambda p: pattern_matrix[p, t]==1, range(npatt))
				## Print selected patterning functions
				if (not selected):
					print("t=" + str(t) + ": No pattern function applied")
				else:
					get_vgene = lambda i : multi_binary_dict.get("vgenes")[
						multi_binary_dict.get("map_to_vgenes")[
							level_idx_lists[i][1][0]
						][0]
					][0]
					morphogens = map(get_vgene, selected)
					print("t=" + str(t) + ": morphogen" + ifthenelse(len(selected)>1, "s", "") 
						+ " " + concat(morphogens, ", "))
					## Prints diffusion for each morphogen
					print(concat(["Field #" + str(i) for i in range(nfields)], " | "))
					for i in range(len(selected)):
						print(morphogens[i] + ": " 
							+ concat([str(x[1]) 
						for x in level_idx_lists[selected[i]][0]], " | "))
		## COMPUTATION OF THE TRAJECTORY MATRIX
		## (i.e. evolution of system state through time and position)
		## trajectories[exp, t, n] = system state (of size m = #genes) at step #t, 
		## in field #n, in experiment #exp
		trajectories = np.zeros((nexp, k+1, nfields, m))
		for e in range(nexp):
			for t in range(k+1):
				for n in range(nfields):
					vect = get_binary_dec(M[bv_vars.get("states")[e][t][n]], m)
					trajectories[e, t, n, 0:] = np.reshape(map(int, vect), (1, m))
		if (verbose):
			print("\n* Trajectories")
			for exp in range(len(bv_vars.get('states'))):
				print("\n** Trajectory #" + str(exp+1))
				for t in range(k+1):
					print("\n-- Time step t=" + str(t))
					print(concat(multi_binary_dict.get("genes"), " | "))
					for n in range(nfields):
						print(get_binary_dec(M[bv_vars.get('states')[exp][t][n]], m))
		if (debug):
			print("\n* State transitions")
			for exp in range(len(bv_vars.get('states'))):
				print("\n** Trajectory #" + str(exp+1))
				for n in range(nfields):
					print("\n* Field #" + str(n+1))
					print("-- Time:" + " "*len(str(k)) + concat(multi_binary_dict.get("genes"), " | "))
					for t in range(k+1):
						pr = lambda name : [get_binary_dec(M[bv_vars.get(name)[exp][t][n]], m)[i] + " "*(len(multi_binary_dict.get("genes")[i])-1) for i in range(m)]
						print("-- t = " + str(t) + " "*(len(str(k))-len(str(t))) + ": " + concat(pr("states"), " | "))
						if (t < k):
							if (len(bv_vars.get('patterning')) > t):
								patt_select = get_binary_dec(M[bv_vars.get('patterning')[t]], npatt)
								print("Patterning selection: [" + patt_select + "]")
								lst = filter_nones([ifthenelse(int(patt_select[i]), 
									patterns[i]['morphogen']) for i in range(npatt)])
								if (len(lst)>0):
									print("Morphogen(s): " + concat(lst, ","))
								else:
									print("Morphogen(s): none")
							print("=> (Patt)" + " "*len(str(k)) + concat(pr("updated_states"), " | "))
							print("=> (GRNs)" + " "*len(str(k)) + concat(pr("diffusion_states"), " | "))
							print("=> (Diff) -> ")
		## GET THE RESULTS OF GRNs
		## (if change_grn == True)
		if ("grns_bv" in bv_vars):
			grns = np.zeros((nexp, k, nfields, len(bvectors), len(multi_binary_dict["genes"])))
			all_grns = [ffkq for q in bv_vars.get("grns_bv") for kq in q for fkq in kq for ffkq in fkq]
			#print([M[x] for x in all_grns])
			for x in bv_vars.get("grns_bv"):
				for t in range(k):
					for n in range(nfields):
						for j in range(m):
							vect = get_binary_dec(M.evaluate(x[t][n][j], model_completion=False), len(bvectors))
							grns[e, t, n, 0:, j] = np.reshape(map(int, vect), (1, len(bvectors)))
			if (verbose):
				print("\n* GRNs")
				for exp in range(nexp):
					print("\n** GRN in exp=" + str(exp+1))
					for t in range(k):
						print("\n-- Time step t=" + str(t))
						for n in range(nfields):
							print("\n-- Field #" + str(n+1))
							for j in range(m):
								print("\n" + concat(multi_binary_dict.get("genes"), " | "))
								for b in range(len(bvectors)):
									print(concat(map(str, bvectors[b]), " | ") + ":" 
										+ get_binary_dec(
										M[bv_vars.get('grns_bv')[exp][t][n][j]], 
										len(bvectors))[b]  + " for gene " 
										+ multi_binary_dict["genes"][j])		
		else:
			grns = None
		res_list.append([pattern_matrix, trajectories, grns])
	if ("patterning" in uniqueness and pattern_matrix0==None):
		## Uniqueness of models: selection of patterning         ##
		## Add condition patterning != SOLUTION[patterning]      ##
		print("\n-- Adding uniqueness condition for patterning variables")
		s = difference_model(s, bv_vars, M, "patterning", verbose=verbose, debug=debug)
	if ("grns" in uniqueness and "grns_bv" in bv_vars.keys()):
		## Uniqueness of models: selection of GRNs               ##
		## Add condition {GRNs} != {SOLUTION[GRNs]}              ##
		print("\n-- Adding uniqueness condition for GRN variables")
		all_grns = [ffkq for q in bv_vars.get("grns_bv") for kq in q for fkq in kq for ffkq in fkq]
		bv_vars.update({"all_grns": all_grns})
		s = difference_model(s, bv_vars, M, "all_grns", verbose=verbose, debug=debug)
	if ("states" in uniqueness):
		## Uniqueness of models: state trajectories              ##
		## Add condition {states} != {SOLUTION[states]}          ##
		print("\n-- Adding uniqueness condition for state variables")
		all_states = [fkq for q in bv_vars.get("states") for kq in q for fkq in kq]
		bv_vars.update({"all_states": all_states})
		s = difference_model(s, bv_vars, M, "all_states", verbose=verbose, debug=debug)
	return([res_list, s, no_more_model])

##########################
## PATTERNING INFERENCE ##
##########################

#' Solve an instance of the patterning inference problem
#' 
#' Given a set of putative patterning factors, the number of steps,
#' some info about experimental data, find the order, the morphogen(s),
#' the diffusion type and the source of the patterning agent. 
#' 
#' @param params output of read_files function in models.py
#'               that describes the model and the observations
#' @param pattern_matrix0 if not None, allows predictions using the patterning matrix 
#'        	(of selected patterning functions) provided
#' @param solmax if not None, the solver will enumerate up to solmax solutions
#' @return res_list list of models + list of possible patterns + mapping between binary and multi-level variables
def pattern_solver(params, pattern_matrix0=None, solmax=None, verbose=False, debug=False):
	[genes, patterns, directives, constants, fields, Observations, Fixpoint, GRNs] = params
	multi_binary_dict = multi_to_binary(genes)
	nfields = len(fields)
	m = len(multi_binary_dict.get("genes"))
	k = directives.get("nsteps")
	uniqueness = directives.get("uniqueness")
	diffusion_rate = constants.get("diffusion-rate")
	## Pattern functions should be selected before this step
	patterning_step = directives.get("patterning_step")
	npatt = len(patterns)
	exp_names = list(set([oo["name"] for oo in Observations]))
	nexp = len(exp_names)
	ngrns = len(GRNs.items())
	solmax = ifthenelse(not solmax, directives.get("limit"), solmax)
	if (not(pattern_matrix0==None)):
		## Uniqueness in pattern selection is not required
		## Need uniqueness at least for state trajectories
		uniqueness = filter(lambda x:not(x == "patterning"), list(set(uniqueness+['states'])))
	if (verbose):
		print("\nPARAMETER VALUES:")
		lstt = [multi_binary_dict.get("genes"), ["id " + str(f[0]) + " = " + concat(f[1], ",") for f in fields], nfields, m, 
			k, npatt, solmax, 
			nexp, ngrns, patterning_step,
			aggregation_function, application_function, uniqueness,
			max_nb_patterns_per_level, max_nb_patterns, min_nb_patterns, max_nb_pattern_times, change_grn,
			diffusion_rate]
		lst = ["genes", "fields", "#fields", "#genes", 
			"#steps", "#patterns", "solmax", 
			"#experiments", "#grns", "patterning end step",
			"aggregation", "application", "uniqueness", 
			"max. #patt/step", "max. #patt", "min. #patt", "max. #times/patt", "allowing GRN changes?",
			"diffusion rate"]
		for i in range(len(lst)):
			print(lst[i] + ": " + str(lstt[i]))
	#____________________________________________________#
	#  Initialization of constants and variables         #
	#____________________________________________________#
	s = Solver()
	## Patterning selection matrix  		    ##
	## Same for every experiment                        ##
	## patterning_step should be <= nsteps              ##
	nsteps = ifthenelse(patterning_step==0, k, min(patterning_step, k))
	if (pattern_matrix0 == None):
		patterning = [BitVec('Step^' + str(i), npatt) for i in range(nsteps)]
	else:
		_, nsteps = np.shape(pattern_matrix0)
		pattern_selection = [filter(lambda j: pattern_matrix0[j, i]==1, range(npatt)) for i in range(nsteps)]
		patterning = [idxList2BV(pattern_selection[i], npatt) for i in range(nsteps)]
	## Phenotype for each field at each step	    ##
	state_name = lambda i, e, j : 'State^' + str(i) + '_(' + str(e) + ', ' + str(j) + ')'
	states = [[[BitVec(state_name(i, e, j), m) for j in range(nfields)] for i in range(k+1)] for e in range(nexp)]
	updated_states = []
	diffusion_states = []
	## Pattern vectors for each possible patterning function
	pattern_name = lambda i, j, p : 'Pattern^(t=' + str(i) + ',' + str(j) + ')_' + str(p)
	ei_ls = [[[BitVec(pattern_name(i, j, p), m) for j in range(nfields)] for i in range(nsteps)] for p in range(npatt)]
	patterning_name = lambda i, j : 'Pattern^(t=' + str(i) + ',' + str(j) + ')'
	e_ls = [[BitVec(patterning_name(i, j), m) for j in range(nfields)] for i in range(nsteps)]
	## Diffusion vectors from each possible gene product from each possible source field at each step
	diffusion_name = lambda e, i, j, g, so : 'Diffusion^(e=' + str(e) + ',t=' + str(i) + ',nf=' + str(j) + ')_g=' + str(g) + '_s=' + str(so)
	di_ls = [[[[[BitVec(diffusion_name(e, i, j, g, so), m) for j in range(nfields)] for i in range(k)] for g in range(m)] for so in range(nfields)] for e in range(nexp)]
	## Aggregation of diffusion vectors for each gene product in each field
	diffusing_name = lambda e, i, j : 'Diffusion^(e=' + str(e) + ',t=' + str(i) + ',' + str(j) + ')'
	d_ls = [[[BitVec(diffusing_name(e, i, j), m) for j in range(nfields)] for i in range(k)] for e in range(nexp)]
	## GRN for each field at each step		    ##
	## Every provided GRN must contain functions for all##
	## genes                                            ##
	## Pair: integer identifier (observation) x GRN     ##
	## First convert GRNs to truth tables               ##
	[vectors, GRNs, npossibilities] = dict_to_tt(GRNs, multi_binary_dict, verbose=verbose, debug=debug)
	## Convert enumerated multi-level vectors to        ##
	## Boolean vectors                                  ##
	bvectors = multi_to_boolean(vectors, m, multi_binary_dict, verbose=verbose, debug=debug)
	## STARTING POINT                                   ##
	runtime = time()
	#____________________________________________________#
	#  Conditions on multi-level variables               #
	#____________________________________________________#
	s = condition_multi_level(s, multi_binary_dict, states, verbose=verbose, debug=debug)
	#____________________________________________________#
	#  Conditions on pattern selection                   #
	#____________________________________________________#
	s = condition_pattern_selection(s, patterning, npatt, nsteps, patterns, multi_binary_dict, verbose=verbose, debug=debug)
	#____________________________________________________#
	#  Conditions on GRNs (from observations)            #
	#____________________________________________________#
	if (change_grn):
		## Then some of the bit-vectors are unknown variables
		## grns_bv is a list of length #experiments of lists of length #steps
		## of lists of length #fields of lists of bit-vectors (of length #boolean vectors) of length #genes
		## s is updated with the logical constraints from the observations
		## on GRNs
		s, grns_bv = convert_grn_bv(s, Observations, GRNs, m, k, nexp, nfields, multi_binary_dict, exp_names, npossibilities, 
			verbose=verbose, debug=debug)
	else:
		## Convert truth tables (i.e. matrices) to constant bit-vectors
		## grns_bv is a list of length #GRNs of lists of length #fields of lists 
		## of bit-vectors of length #genes
		grns_bv = convert_grn(GRNs, m, multi_binary_dict, npossibilities, verbose=verbose, debug=debug)
	#____________________________________________________#
	#  Conditions on observations (from observations)    #
	#____________________________________________________#
	s = condition_observation(s, Observations, exp_names, multi_binary_dict, states, verbose=verbose, debug=debug)
	## Fixpoints                                        ##
	if (Fixpoint):
		s = condition_fixpoint(s, Observations, exp_names, Fixpoint, states, k, m,
			grns_bv=ifthenelse(change_grn, grns_bv, None), verbose=verbose, debug=debug)
		fixpoint_steps = [int(min(Fixpoint[key].keys())) for key in Fixpoint.keys()]
	else:
		fixpoint_steps = []
	#____________________________________________________#
	#  Conditions on system state transitions            #
	#____________________________________________________#
	apply_functions = directives.get("apply")
	## Precompute patterns                              ##
	if (verbose):
		print("\n--- Precomputation of patterning functions")
	level_idx_lists = [precompute_pattern(fields, multi_binary_dict, patterns[p], verbose=debug, debug=debug)
			 for p in range(npatt)]
	if (verbose):
		print("\n--- Precomputation of diffusion functions")
	diff_idx_lists = [[precompute_diffusion(fields, multi_binary_dict, g, so, diffusion_rate, verbose=debug, debug=debug) 
		for so in range(nfields) ] for g in range(m)]
	for e in range(nexp):
		updated_states_exp = []
		diffusion_states_exp = []
		for step in range(k):
			## grn_bv is a list of length nfields = #fields 
			## of lists of bit-vectors (of length m = #genes)
			if (not change_grn):
				obs = filter(lambda x:x["name"]==exp_names[e] and not(x["GRN"]==None), Observations)
				## One per field because change_grn == False
				## meaning that the GRN AND INITIAL CONDITIONS *HAVE TO* be specified for each field in the observations file
				grns_names = [[x["field"], x["GRN"]] for x in obs]
				if (len(grns_names) > nfields):
					print("ERROR: Too many values for GRNs!")
					raise ValueError
				tmp = sorted([[x[0], grns_bv[GRNs.keys().index(x[1])]] for x in grns_names], key=lambda x:x[0])
				grn_bv = [x[1] for x in tmp]
			else:
				grn_bv = grns_bv[e][step]
			## Transition (Q_previous -> Q_updated (patterning) -> Q_diffusion -> Q_next) for each step, field
			Q_next = states[e][step+1]
			Q_previous = states[e][step]
			updated_name = lambda nf : "State-updated^" + str(step) + "_(" + str(e) + ',' + str(nf) + ')'
			Q_updated = [BitVec(updated_name(nf), m) for nf in range(nfields)]
			updated_states_exp += [Q_updated]
			diffusion_name = lambda nf : "State-diffusion^" + str(step) + "_(" + str(e) + ',' + str(nf) + ')'
			Q_diffusion = [BitVec(diffusion_name(nf), m) for nf in range(nfields)]
			diffusion_states_exp += [Q_diffusion]
			##---------- PATTERN FUNCTION   -----------##
			s = condition_patterning(s, step, nsteps, nfields, Q_previous, Q_updated, level_idx_lists, m, 
				npatt, patterning, ei_ls, e_ls, multi_binary_dict, verbose=verbose, debug=debug)
			##----------- APPLY FUNCTIONS   -----------##
			#TODO other apply functions
			## --- STATE TRANSITION AT EACH STEP    ---##
			if ("update" in apply_functions):
				if (verbose):
					print("\n* Condition Q(t=" + str(step) 
						+ ") -> Q(t=" + str(step+1) + ") in " + exp_names[e])
				## If fixpoint => synchronous transition condition
				update_sync = directives.get("updates") == "sync" or step in fixpoint_steps
				s = condition_transition(s, nfields, Q_diffusion, Q_updated, update_sync, m, 
						grn_bv, bvectors, verbose=verbose, debug=debug)
				## At least "update" in apply_functions so 
				## no other default condition for Q_updated 
				## (e.g. Q_updated == Q_previous)
			##--------- DIFFUSION FUNCTIONS   ---------##
			## (cell signalling)                       ##
			s = condition_diffusion(s, step, nfields, Q_diffusion, Q_next, m, di_ls[e], 
				d_ls[e], multi_binary_dict, diff_idx_lists, verbose=verbose, debug=debug)
		updated_states += [updated_states_exp]
		diffusion_states += [diffusion_states_exp]
	#____________________________________________________#
	#  Solution processing                               #
	#____________________________________________________#
	bv_vars = {"patterning": patterning, "states": states, "updated_states": updated_states, "diffusion_states": diffusion_states}
	if (change_grn):
		bv_vars.update({"grns_bv": grns_bv})
	solution_no, res_list, no_more_model = 0, [], False
	while (solution_no < solmax and not no_more_model):
		[res_list, s, no_more_model] = get_solutions(s, nexp, m, npatt, nfields, k, nsteps, bv_vars,
			level_idx_lists, multi_binary_dict, bvectors, patterns, uniqueness, res_list=res_list, 
			pattern_matrix0=pattern_matrix0, verbose=False, debug=debug)
		solution_no += 1
	if (solution_no == solmax and verbose):
		print("\n--- Maximum number of solutions (" + str(solmax) + ") reached.")
	if (no_more_model and verbose):
		print("\n--- " + str(solution_no-1) + " solution" + ifthenelse(solution_no-1<2, "", "s") + ".")
	## ENDING POINT                                     ##
	print("\nRUNTIME = " + str(round(time()-runtime, 2)) + " sec.\n")
	return([res_list, patterns, multi_binary_dict, bvectors])
