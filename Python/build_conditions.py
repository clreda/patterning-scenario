# -*- coding: utf-8 -*-

from z3 import *

from utils import *
from patterns import compute_diffusion, aggregate_vectors, apply_vectors

## All of these functions return an updated Z3 solver object

#____________________________________________#
#  CONDITIONS ON multi-level variables       #
#____________________________________________#
#' These are the logical constraints related to multi-level 
#' variables

#' @param s Z3 solver object
#' @param multi_binary_dict dictionary which contains the full mapping from binary
#'      variables to real genes and the way back
#' @param states list of Z3 bit-vectors associated with each state at a given time step, in a given field, ...
#' @param verbose logical for printing messages
#' @param debug logical for printing (hopefully) useful insights to the code
#' @return s Z3 solver object updated with conditions on binary variables associated with
#'  multi-level genes (avoids the selection of inconsistent states such as e.g. A-1=0 and A-4=1)
def condition_multi_level(s, multi_binary_dict, states, verbose=False, debug=False):
	if (verbose):
		print("\n-- Conditions on values of same-gene binary variables")
	## states: bit-vectors of size m = #binary variables
	all_states = [states_e_t_f for states_e in states for states_e_t in states_e for states_e_t_f in states_e_t]
	## get indices of binary variables associated with each "real" gene
	## in increasing order of corresponding concentration level 
	## (only on variables having strictly more than one concentration level)
	all_ids = filter(lambda x:len(x)>1, multi_binary_dict["map_list"])
	coef = lambda state, j, q : extract1BV(state, all_ids[j][q])
	## Condition on (consecutive+transitivity) concentration levels of a same gene
	## value(smaller conc. level) <= value(greater conc. level)
	cond_ml = lambda x, y : UGE(x, y)
	if (not debug):
		cond_ls = [map_and([cond_ml(coef(state, j, q), coef(state, j, q+1)) 
				for state in all_states 
				for q in range(len(all_ids[j])-1)]) 
			for j in range(len(all_ids))]
	else:
		## Same but unfold code
		cond_ls = []
		for j in range(len(all_ids)):
			print("--- gene " + multi_binary_dict["vgenes"][j][0])
			for state in all_states:
				print("------ state " + str(state))
				for q in range(len(all_ids[j])-1):
					print("Condition " + multi_binary_dict["vgenes"][j][0] + "-" + str(q+1) 
						+ " <= " + multi_binary_dict["vgenes"][j][0] + "-" + str(q+1+1))
					cond_ls += [cond_ml(coef(state, j, q), coef(state, j, q+1))]
	cond = map_and(cond_ls)
	s.add(cond)
	return(s)

#___________________________________#
#  CONDITIONS ON observations       #
#___________________________________#
#' These are the logical constraints related to observations 
#' described in the related files (about phenotypes and GRNs)

#' @param s Z3 solver object
#' @param Observations dictionary of observations as described in models.py: name, step, field, phenotype, GRN name
#' @param exp_names list of character strings (observation names, as described in observations file)
#' @param multi_binary_dict dictionary which contains the full mapping from binary
#'      variables to real genes and the way back
#' @param states list of Z3 bit-vectors associated with each state at a given time step, in a given field, ...
#' @param verbose logical for printing messages
#' @param debug logical for printing (hopefully) useful insights to the code
#' @return s Z3 solver object updated with conditions in observations (fixed values for the documented observations)
def condition_observation(s, Observations, exp_names, multi_binary_dict, states, verbose=False, debug=False):
	mapping = multi_binary_dict.get("mapping")
	map_list = multi_binary_dict.get("map_list")
	if (verbose):
		print("\n--- Condition on observations")
	## Build conditions on observations         ##
	## Get state variable associated with a given observation obs
	get_state = lambda obs : states[exp_names.index(obs["name"])][obs["step"]][obs["field"]]
	get_idx_list = lambda phenotype, active : ifthenelse(active, 
						map_list[mapping[phenotype[0]]][:int(phenotype[1])],
						map_list[mapping[phenotype[0]]][int(phenotype[1]):])
	if (not debug):
		cond_ls = [
			map_and([
				## In Observation:exp_names.index(obs["name"])
				## for observation G is above threshold expression level #i-1 
				## at time step:obs["step] and in field:obs["field"]
				map_and([
						## State^(obs["step"])_(obs["field"])[idx] == 0
						## for any idx associated with binary variables for this gene G
						## which correspond to expression levels strictly higher than i
						cond_value(get_state(obs), idx, false) 
						for idx in get_idx_list(phenotype, False)
					]
					+
					[
						## State^(obs["step"])_(obs["field"])[idx] == 1
						## for any idx associated with binary variables for this gene G
						## which correspond to expression levels lower than i
						cond_value(get_state(obs), idx, true) 
					for idx in get_idx_list(phenotype, True)
					]) 
				for phenotype in obs["phenotype"]]) 
			for obs in Observations]
	else:
		## Same code, but unfold
		cond_ls = []
		for obs in Observations:
			field, step, phenotype_obs = obs["field"], obs["step"], obs["phenotype"]
			print("\n-------- OBSERVATION \'" + obs["name"] + "\', step #" 
				+ str(step) + ", field #" + str(field))
			print("Phenotype for field #" + str(field) + ", step #" 
				+ str(step) + ": " + str(phenotype_obs))
			Q = get_state(obs)
			for phenotype in phenotype_obs:
				i = mapping[phenotype[0]]
				i_levels = map_list[i]
				value = int(phenotype[1])
				if (verbose):
					print("Index associated with gene " + phenotype[0] 
						+ " in the list of genes: " + str(i))
					print("Indices in binary variable list associated with " 
						+ phenotype[0] + ": " + str(i_levels))
					print("In observation: " + phenotype[0] + " = " + str(value))
				## Binary variables associated with strictly higher concentration levels are set to 0
				cond_ls_higher = [cond_value(Q, idx, false) for idx in i_levels[value:]]
				if (verbose):
					print("Conditions of inactivation of > concentration levels (if existent):")
					print(cond_ls_higher)
				## Binary variables associated with lower and considered values of GE are set to 1
				cond_ls_lower = [cond_value(Q, idx, true) for idx in i_levels[:value]]
				if (verbose):
					print("Conditions of activation of <= concentration levels (if existent):")
					print(cond_ls_lower)
				cond_ls += [map_and(cond_ls_higher+cond_ls_lower)]
	## Aggregate all these conditions with a huge "And" operator
	cond = map_and(cond_ls)
	if (debug):
		print("Conditions added to solver: " + str(cond))
	s.add(cond)
	return(s)

#_____________________________________#
#  CONDITIONS ON state transitions    #
# (only update function)              #
#_____________________________________#
#' These are the logical constraints related to state updates
#' ruled by GRNs

#' Build transition condition for the next state
#' 
#' @param Q_next future state
#' @param g integer identifier associated with gene
#' @param grn_bv list of Z3 Bit-Vectors associated with current step and field
#' @param i index of the Boolean vector (of binary variable values) equal
#'          to the previous state q_updated
#' @return cond transition condition for a given gene g 
#'               "q_next[g] == GRN(q_updated)[g]" in {0, 1}
def aux_condition_transition(Q_next, g, grn_bv, i):
	return(cond_value_vec(Q_next, grn_bv[g], g, i))

#' Build synchronous transition condition
#'
#' @param Q_updated previous state
#' @param Q_next future state
#' @param m number of binary variables
#' @param grn_bv list of Z3 Bit-Vectors associated with current step and field
#' @param i index of the Boolean vector equal to Q_updated
#' @param verbose logical for printing messages
#' @param debug logical for printing (hopefully) useful insights to the code
#' @return cond condition 
#'             "for all gene g, q_next[g] == GRN(q_updated)[g]" in {0, 1} 
def condition_sync_transition(Q_updated, Q_next, m, grn_bv, i, verbose=False, debug=False):
	cond_ls = [aux_condition_transition(Q_next, g, grn_bv, i) 
			for g in range(m)]
	if (verbose):
		for g in range(m):
			print(cond_ls[g])
			print("GRN(q) = " + str(get_binary_dec(grn_bv[g], grn_bv[g].size())[i]))
	return(map_and(cond_ls))

#' Build asynchronous transition condition
#'
#' @param Q_updated previous state
#' @param Q_next future state
#' @param m number of binary variables
#' @param grn_bv list of Z3 Bit-Vectors associated with current step and field
#' @param i index of the Boolean vector equal to Q_updated
#' @param verbose logical for printing messages
#' @param debug logical for printing (hopefully) useful insights to the code
#' @return cond condition 
#'              "it exists a gene g, [ (q_next[g] == GRN(q_updated)[g]) 
#'              AND (for all gene g'!=g, q_next[g'] == q_updated[g']) ]" 
def condition_async_transition(Q_updated, Q_next, m, grn_bv, i, verbose=False, debug=False):
	## <=> Or on coordinates that can be modified (one and only one is modified)
	cond_ls = [And(equal_except1BV(Q_updated, Q_next, g), aux_condition_transition(Q_next, g, grn_bv, i)) 
			for g in range(m)]
	if (verbose):
		for g in range(m):
			print("\n--- Binary variable #" + str(g) + " is modified:")
			print(cond_ls[g])
			print("GRN(q)[" + str(g) + "] = " + str(get_binary_dec(grn_bv[g], grn_bv[g].size())[i]))
	return(map_or(cond_ls))

#' Build the full transition condition for all genes,
#' for all fields, at a given step of the experiment
#' 
#' @param s Z3 solver object
#' @param nfields number of fields
#' @param Q_next future system state
#' @param Q_updated previous system state
#' @param update_sync Is the update synchronous?
#' @param m number of binary variables
#' @param grn_bv list of Z3 Bit-Vectors associated with current step and field
#' @param bvectors list of Boolean vectors associated with all possible values of the binary variables
#' @param verbose logical for printing messages
#' @param debug logical for printing (hopefully) useful insights to the code
#' @return s updated solver with conditions Q_updated -> Q_next
def condition_transition(s, nfields, Q_next, Q_updated, update_sync, m, grn_bv, bvectors, verbose=False, debug=False):
	update = ifthenelse(update_sync, condition_sync_transition, condition_async_transition)
	if (not debug):
		cond_ls = [
			## <=> Or on the set of Boolean vectors and at least one is chosen
			map_and([
				## If for all variable g, Q_updated[g] == (Boolean Vector #i)[g]
				Implies(cond_equal_value(Q_updated[nf], bvectors[i]), 
				## Then for all variable g, Q_next[g] == GRN_g(Q_updated)[g]
				update(Q_updated[nf], Q_next[nf], m, grn_bv[nf], i))
				for i in range(len(bvectors))
			]) 
			for nf in range(nfields)]
	from copy import deepcopy
	if (debug):
		print("Condition on " + ifthenelse(update_sync, "sync.", "async.") 
			+ " system state transition added to solver:")
		## Same unfold code
		cond_ls = []
		for nf in range(nfields):
			print("- AND_field=" + str(nf))
			cd_ls = []
			for i in range(len(bvectors)):
				print("-- OR_bvector=" + str(i))
				hyp = cond_equal_value(Q_updated[nf], bvectors[i])
				ccl = update(Q_updated[nf], Q_next[nf], m, grn_bv[nf], i, verbose=verbose)
				## GRN vectors associated with current field (for all genes)
				grnq = [get_binary_dec(grn_bv[nf][j], grn_bv[nf][j].size()) for j in range(m)]
				if (update_sync):
					print(" "*4 + "Q_prev == " + concat(bvectors[i]) 
						+ " => Q_next == " + concat([grnq[j][i] for j in range(m)]))
				else:
					for j in range(m):
						bvb, res_bvb = deepcopy(bvectors[i]), deepcopy(bvectors[i])
						print("-"*4 + "OR_updated_gene=" + str(j) 
							+ " (" + str(res_bvb[j]) + " -> " + str(grnq[j][i]) + ")")
						res_bvb[j] = grnq[j][i]
						print(" "*8 + "Q_prev == " + concat(bvb) + " => Q_next == " + concat(res_bvb))
				cd_ls += [Implies(hyp, ccl)]
			cond_ls += [map_and(cd_ls)]
	## Apply this for every field
	cond = map_and(cond_ls)
	s.add(cond)
	return(s)

#_____________________________________#
#  CONDITIONS ON fixpoints            #
#_____________________________________#
#' These are the logical constraints related to fixpoints
#' in observations (related to phenotypes and GRNs)

#' @param s Z3 solver object
#' @param Observations dictionary of observations (as described in models.py)
#' @param exp_names character string list of observation names
#' @param Fixpoint dictionary that summarizes the fixpoint conditions (as described in models.py)
#' @param states list of Z3 bit-vectors associated with system states at each time step, in each field, in each trajectory
#' @param k maximum length of trajectories
#' @param m number of binary variables
#' @param grns_bv if change_grn=True then a list of length #nexp #steps #fields #genes of variable bit-vectors
#' @param verbose logical for printing messages
#' @param debug logical for printing (hopefully) useful insights to the code
#' @return s updated solver with fixpoint conditions
def condition_fixpoint(s, Observations, exp_names, Fixpoint, states, k, m, grns_bv=None, verbose=False, debug=False):
	if (verbose):
		print("\n--- Conditions for fixpoints in observations")
	get_state = lambda key, stime, field : states[exp_names.index(key)][stime][int(field)-1]
	if (grns_bv):
		get_grn = lambda key, stime, field : grns_bv[exp_names.index(key)][stime][int(field)-1]
	else:
		## To comply with Python typing
		get_grn = lambda a, b, c : [None]*m
	get_start_step_fixpoint = lambda key : int(min(Fixpoint[key].keys()))
	get_fields = lambda key : Fixpoint[key][str(get_start_step_fixpoint(key))]
	get_fields_list = lambda key : get_fields(key).keys()
	get_fixpoint = lambda key, field, by : get_fields(key)[field].get(by)
	cond_ls = []
	for key in Fixpoint.keys():
		if (debug):
			print("-- Experiment = " + key)
		## From which step?
		step = get_start_step_fixpoint(key)
		fields = get_fields_list(key)
		for field in fields:
			if (debug):
				print("- Field = " + str(field) + ", Fixpoint start step #" + str(step))
			## Add constraint Q[i] == Q[i+1] for k > i >= time step for fixpoint
			## where Q is a relevant set of states indexed by time steps
			if (get_fixpoint(key, field, "condition")):
				for stime in range(step, k):
					if (debug):
						print("(Phenotype) Step #" + str(stime) + " -> #" + str(stime+1))
					Q_prev = get_state(key, stime, field)
					Q_ne = get_state(key, stime+1, field)
					## State transition condition is ensured by main loop in pattern_solver.py
					cond_ls += [Q_prev == Q_ne]
			## Constraint on GRNs
			if (get_fixpoint(key, field, "GRN")):
				for stime in range(step, k-1):
					if (debug):
						print("(GRN) Step #" + str(stime) + " -> #" + str(stime+1))
					GRN_prev = get_grn(key, stime, field)
					GRN_ne = get_grn(key, stime+1, field)
					for j in range(m):
						cond_ls += [GRN_prev[j] == GRN_ne[j]]
	## cond_ls can be empty if the fix point constraint is on the last time step
	if (cond_ls):
		## Aggregate those conditions with a huge "AND" operator
		cond = map_and(cond_ls)
		if (debug):
			print("\nCondition:\n" + str(cond))
		s.add(cond)
	return(s)

#_____________________________________#
#  CONDITIONS ON patterns             #
#_____________________________________#
#' These are the logical constraints related to patterning
#' and pattern selection (at each step or globally)

#' Replace consecutive nonnegative bits
#' by a single nonnegative bit
#'
#' @param s Z3 solver object
#' @param q Z3 Bit-Vector
#' @param name character string identifier for
#'             resulting Bit-Vector
#' @return res pair (updated solver, bq)
#'         bq is q Z3 Bit-Vector of the same size
#'         where repetitions of 1's were removed
#' number of 1 in bq = number of consecutive 1's in q
def reduce_bv(s, q, name):
	n = q.size()
	bq = BitVec(name, n)
	## First value is kept: bq[0] = q[0]
	cond_ls = [If(cond_value(q, 0, true), 
			cond_value(bq, 0, true), 
			cond_value(bq, 0, false))]
	## For 1 <= i <= n
	## q[i-1] == q[i] == 1 => bq[i] = 0 
	## q[i-1] == q[i] == 0 => bq[i] = 0
	## q[i-1] = 0 and q[i] = 1 => bq[i] = 1
	## q[i-1] = 1 and q[i] = 0 => bq[i] = 0
	## <=> If(q[i-1] == 0 and q[i] == 1, 
	##          bq[i] == 1, bq[i] == 0)
	cond_ls += [If(And(cond_value(q, i-1, false), cond_value(q, i, true)),
		cond_value(bq, i, true),
		cond_value(bq, i, false)
		) for i in range(1, n)]
	cond = map_and(cond_ls)
	s.add(cond)
	return([s, bq])

#' Get indices of similar patterns
#'
#' @param patterns list of dictionaries of features for each patterning function
#' @param npatt number of patterning functions (= length of list patterns)
#' @param eq_feature_list list of feature names which values should be equal for similar patterns
#' @return patterns_list list of lists of indices of similar patterns
def morphogen_pattern_ids(patterns, npatt, eq_feature_list=["morphogen"]): #["morphogen", "source", "rate"]
	done = [False]*npatt
	patterns_list = []
	equality = lambda x, y: map_bool_and([x.get(feature)==y.get(feature) for feature in eq_feature_list])
	while (not all(done)):
		i = done.index(False)
		lst = map(lambda x: ifthenelse(not(x==i),
			equality(patterns[x], patterns[i]),
			True
		), range(npatt))
		done = map(lambda i: lst[i] or done[i], range(npatt))
		lst_ids = filter(lambda i:lst[i], range(npatt))
		if (len(lst_ids) > 1):
			patterns_list += [lst_ids]
	return(patterns_list)

#' UNIQUE-ONE boolean operator to a set of coefficients of a Boolean vector
#' (associativity, commutativity)
#' 
#' @param coefs list of 1-sized bit vectors
#' @return res the equivalent 1-sized bit vector of setting only one
#'    of the input coefficients to 1 (all other coefs are set to 0)
def uone_op(coefs):
	n = len(coefs)
	if (n==0):
		return(true)
	if (n==1):
		return(coefs[0])
	if (n==2):
		## XOR
		return(~coefs[0] & coefs[1] | coefs[0] & ~coefs[1])
	neg_coefs = [map(lambda x:~x, filter_nones([ifthenelse(not(i == j), coefs[j]) for j in range(n)])) for i in range(n)]
	bv_ls = [reduce(lambda x,y: x&y, [coefs[i]] + neg_coefs[i]) for i in range(n)]
	return(reduce(lambda x,y : x | y, bv_ls))

#' UNIQUE-ONE boolean operator to a set of coefficients of a Boolean vector
#' (associativity, commutativity OK) OR all coefficients are zero 
#'
#' @param bv Z3 Bit-Vector
#' @param ids a list of indices in bv such as 
#'          UNIQUE-ONE_{i in ids} bv[i] OR AND_{i in ids} ~bv[i] should be satisfied
#' @return cond condition "UNIQUE-ONE_{i in ids} bv[i] OR AND_{i in ids} ~bv[i]"
def one_selection(bv, ids):
	coefs = [extract1BV(bv, i) for i in ids]
	all_zero = lambda ls : reduce(lambda x,y: x & y, map(lambda x:~x, ls))
	## Either only one single patterning function is selected
	## among those indices, either none is selected
	cond = map_or([x==true for x in [uone_op(coefs), all_zero(coefs)]])
	return(cond)

#' Condition for patterning function selection
#' (relevant parameter values can be modified in globals.py)
#'
#' @param s Z3 solver object
#' @param patterning list of pattern selection vectors (Bit-Vectors of size npatt) ordered by time step
#' @param npatt number of patterning functions
#' @param nsteps maximum number of selected patterns
#' @param patterns list of dictionaries of features for each patterning function
#' @param multi_binary_dict dictionary which contains the full mapping from binary
#'      variables to real genes and the way back
#' @param verbose logical for printing messages
#' @param debug logical for printing (hopefully) useful insights to the code
#' @return s updated solver with patterning conditions
def condition_pattern_selection(s, patterning, npatt, nsteps, patterns, multi_binary_dict, test_var=None, verbose=False, debug=False):
	if (verbose):
		print("\n--- Condition on patterning function selection")
	if (test_var):
		keys = test_var.keys()
		for key in keys:
			globals()[key] = test_var[key]
	cond_ls = []
	## Allow at most one patterning function of a given type 
	## (e.g. same morphogen and source)
	## to be selected at one time step
	patterns_list = morphogen_pattern_ids(patterns, npatt)
	cond_ls += [one_selection(patterning[t], np_list) for t in range(nsteps) for np_list in patterns_list]
	## Allow patterning function selection 
	## only before a given time step
	#if (patterning_step):
	#	cond_ls += [patterning[t]==buildZERO(npatt) for t in range(patterning_step, k)]
	## Maximum number of selected patterning 
	## functions per level (= time step)
	## i.e. limit the number of 1's in each time step
	## max_nb_patterns_per_level <= npatt
	if (max_nb_patterns_per_level and max_nb_patterns_per_level <= npatt):
		bv_n = BitVecVal(max_nb_patterns_per_level, npatt)
		cond_ls += [ULE(count_1s(patterning[i]), bv_n) for i in range(nsteps)]
	## Maximum (resp. minimum) number of selected patterning 
	## functions among those provided
	## i.e. limit the number of 1's in the vector {sum(patterning[t])}_t
	## max_nb_patterns <= npatt and min_nb_patterns <= npatt
	if ((max_nb_patterns or min_nb_patterns) and max_nb_patterns <= npatt and min_nb_patterns <= npatt):
		sum_lines = reduce(lambda a,b : a|b, patterning)
		if (max_nb_patterns):
			bv_n_max = BitVecVal(max_nb_patterns, npatt)
			cond_ls.append(ULE(count_1s(sum_lines), bv_n_max))
		if (min_nb_patterns):
			bv_n_min = BitVecVal(min_nb_patterns, npatt)
			cond_ls.append(UGE(count_1s(sum_lines), bv_n_min))
	## Maximum of times any patterning function 
	## is selected (in a row)
	## that is, if a patterning function is
	## selected in several consecutive time steps, 
	## it only counts as one single time
	## i.e. limit the number of consecutive series of 1 
	## in "lines" of patterning
	## max_nb_pattern_times <= nsteps
	if (max_nb_pattern_times and max_nb_pattern_times <= nsteps):
		bv_n = BitVecVal(max_nb_pattern_times, nsteps)
		## Build "lines" of pattern matrix as Bit-Vectors
		lines = [BitVec("line_" + str(n), nsteps) for n in range(npatt)]
		cond_ls += [cond_value_vec(lines[n], patterning[i], i, n) for i in range(nsteps) for n in range(npatt)]
		## Build constraint on number of consecutive 1's in lines
		for n in range(npatt):
			[s, bq] = reduce_bv(s, lines[n], 'ReducedLine_' + str(n))
			cond_ls.append(ULE(count_1s(bq), bv_n))
	cond = map_and(cond_ls)
	if (debug):
		print(cond)
	s.add(cond)
	return(s)

#' Update of morphogen diffusion in the fields at a given time step
#'
#' @param s Z3 solver object
#' @param step current time step
#' @param nsteps maximum number of time steps with pattern selection
#' @param nfields number of fields
#' @param Q_previous list of previous system states (bit-vectors of size m) of length nfields
#' @param Q_updated list of system states updated with the patterning effects (bit-vectors of size m) of length nfields
#' @param level_idx_lists lists of 
#'             ([gene identifier, updated expression value, corresponding binary variable], 
#'               list of same-gene binary variable indices) for each gene
#' @param m number of binary variables
#' @param npatt number of possible patterning functions
#' @param patterning list of Z3 bit-vectors that account for the selection of pattern functions 
#'        at each time step
#' @param ei_ls list of Z3 bit-vectors of patterning vectors (associated with each patterning function, of size m)
#'            at each time step, in each field (independent from experiments)
#' @param e_ls list of Z3 bit-vectors of aggregated patterning vectors (of size m)
#'            at each time step, in each field (independent from experiments)
#' @param multi_binary_dict dictionary which contains the full mapping from binary
#'      variables to real genes and the way back
#' @param verbose logical for printing messages
#' @param debug logical for printing (hopefully) useful insights to the code
#' @return s updated solver with patterning conditions
def condition_patterning(s, step, nsteps, nfields, Q_previous, Q_updated, level_idx_lists, m, npatt, patterning, ei_ls, e_ls, multi_binary_dict, verbose=False, debug=False):
	if (step < nsteps):
		## Morphogen "diffusion"
		cond_ls = [
			## If the pth patterning function is selected
			If(cond_value(patterning[step], p, true), 
				## Then compute the corresponding patterning vectors of size m (for each field)
				compute_diffusion(ei_ls[p][step], nfields, Q_previous, level_idx_lists[p], verbose=verbose, debug=debug),
				## Else ignore it (patterning vector set to 0)
				map_and([ei_ls[p][step][nf] == buildZERO(m) for nf in range(nfields)])
			)
			 for p in range(npatt)
		]
		## "Aggregation" function (interactions between morphogens)
		## See list of possible aggregation functions in patterns.py
		params = {"npatt": npatt, "step": step, "vect_ls": ei_ls, "multi_binary_dict": multi_binary_dict, "type": "patterning"}
		cond_ls += [e_ls[step][nf] == aggregate_vectors(params, nf) for nf in range(nfields)]
		## Patterning effect on phenotypes
		## See list of possible patterning application functions in patterns.py
		params = {"Q_previous": Q_previous, "vect_ls": e_ls, "step": step, "multi_binary_dict": multi_binary_dict, "type": "patterning"}
		cond_ls += [apply_vectors(params, nf) == Q_updated[nf] for nf in range(nfields)]
	else:
		cond_ls = [Q_previous[nf] == Q_updated[nf] for nf in range(nfields)]
	cond = map_and(cond_ls)
	if (debug):
		print("\n--- Condition for patterning:")
		print(cond)
	s.add(cond)
	return(s)

#_______________________________#
#  CONDITIONS ON diffusion of   #
#   gene products               #
#_______________________________#
#' These are the logical constraints related to diffusion
#' of gene products (cell signalling)

#' Update of gene product diffusion in the fields at a given time step
#'
#' @param s Z3 solver object
#' @param step current time step
#' @param nfields number of fields
#' @param Q_diffusion list of previous system states (bit-vectors of size m) of length nfields
#' @param Q_next list of system states updated with the patterning effects (bit-vectors of size m) of length nfields
#' @param m number of binary variables
#' @param di_ls list of Z3 bit-vectors of diffusion vectors (associated with each binary variable, of size m)
#'            at each time step, in each field (independent from experiments), from each source field
#' @param d_ls list of Z3 bit-vectors of aggregated diffusion vectors (of size m)
#'            at each time step, in each field (independent from experiments)
#' @param multi_binary_dict dictionary which contains the full mapping from binary
#'      variables to real genes and the way back
#' @param diff_idx_lists lists of 
#'             ([gene identifier, updated expression value, corresponding binary variable], 
#'               list of same-gene binary variable indices) for each (possible source) field for each gene
#' @param verbose logical for printing messages
#' @param debug logical for printing (hopefully) useful insights to the code
#' @return s updated solver with patterning conditions
def condition_diffusion(s, step, nfields, Q_diffusion, Q_next, m, di_ls, d_ls, multi_binary_dict, diff_idx_lists, verbose=False, debug=False):
	## If binary variable #g is active at step #step in source field #s
	is_nactive = lambda so, g, bv : cond_value(Q_diffusion[so], g, bv)
	## Then compute the corresponding diffusion vectors of size m 
	## (for each field)
	is_diff = lambda step, so, g : Implies(is_nactive(so, g, true), compute_diffusion(di_ls[so][g][step], nfields, Q_diffusion, diff_idx_lists[g][so], verbose=verbose, debug=debug))
	## Else ignore it (patterning vector set to 0)
	is_not_diff = lambda step, so, g : Implies(is_nactive(so, g, false), map_and([di_ls[so][g][step][nf] == buildZERO(m) for nf in range(nfields)]))
	diff = lambda step, s, g : And(is_diff(step, so, g), is_not_diff(step, so, g))
	## Gene product "diffusion" for a given gene product/binary variable g and putative source field s
	cond_ls = [diff(step, so, g) for so in range(nfields) for g in range(m)]
	if (debug):
		for nf in range(nfields):
			for g in range(m):
				gene = multi_binary_dict['genes'][g]
				idx = multi_binary_dict['map_to_vgenes'][g][0]
				vgene = multi_binary_dict['vgenes'][idx][0]
				print("\n\n-- Diffusion of gene " 
					+ vgene
					+ " (" + gene + ")"
					+ " from source field of index " + str(nf))
				print(str(Q_diffusion[nf]) + "[" + gene + "] == 1")
				print(" <=> ")
				print("Binary var.: " + str(diff_idx_lists[g][nf][1]))
				print("Diff/ " + str(diff_idx_lists[g][nf][0]))
				so = diff_idx_lists[g][nf][0][0][0]
				for x in diff_idx_lists[g][nf][0]:
					gene = multi_binary_dict['genes'][x[-1]]
					print(str(di_ls[so][g][step][x[0]]) + "[" + gene + "] == " 
					+ str(int(x[1]>0)))
				print("Else " + str(di_ls[so][g][step][nf]) + " == 0\n")
				print(cond_ls[m*nf+g])
	## "Aggregation" function (interactions between gene products)
	## See list of possible aggregation functions in patterns.py
	params = {"nfields": nfields, "m": m, "step": step, "vect_ls": di_ls, "multi_binary_dict": multi_binary_dict, "type": "diffusion"}
	cond_ls += [d_ls[step][nf] == aggregate_vectors(params, nf) for nf in range(nfields)]
	## Patterning effect on phenotypes
	## See list of possible application functions in patterns.py
	params = {"Q_previous": Q_diffusion, "vect_ls": d_ls, "step": step, "multi_binary_dict": multi_binary_dict, "type": "diffusion"}
	## When no diffusion (diffusion_rate = 0) is equivalent to Q_next[nf] == Q_diffusion[nf]
	cond_ls += [apply_vectors(params, nf) == Q_next[nf] for nf in range(nfields)]
	cond = map_and(cond_ls)
	if (debug):
		print("\n--- Condition for diffusion:")
		print(cond)
	s.add(cond)
	return(s)

#_______________________________#
#  CONDITIONS ON uniqueness     #
#_______________________________#
#' These are the logical constraints related to model uniqueness
#' In order to enumerate models, we need to describe which feature
#' (parameter value) must be modified to distinguish between models

#' Build condition of uniqueness
#'
#' @param s Z3 solver object
#' @param bv_vars list of bit-vector variables
#' @param model model solution returned by the solver
#' @param varname variable name on which uniqueness condition should be applied
#' @param verbose logical for printing messages
#' @param debug logical for printing (hopefully) useful insights to the code 
#' @return s updated solver with uniqueness condition
def difference_model(s, bv_vars, model, varname, verbose=False, debug=False):
	var = bv_vars[varname]
	cond_ls = [v != model[v] for v in var]
	cond = map_or(cond_ls)
	if (debug):
		print("\n--- Condition for obtaining a (potentially) next unique model:")
		print(cond)
	s.add(cond)
	return(s)
