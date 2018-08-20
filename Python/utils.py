# -*- coding: utf-8 -*-

from z3 import *
import numpy as np

from globals import *

#######################
## SHORTCUTS         ##
#######################

def ifthenelse(ifc, thenc, elsec=None):
	return(thenc if (ifc) else elsec)

#######################
## DEBUG             ##
#######################

def add_arg_values(s, argnames, args):
	## length(argnames)==length(args) ##
	for i in range(len(args)):
		s += " " + argnames[i] + "=" + str(args[i]) + ", "
	return(s)

def noBV(adj="", argnames=None, args=None):
	s = "ERROR: A " + adj + " Bit-Vector could not be created: "
	return(add_arg_values(s, argnames, args))

def warnBV(warning="", argnames=None, args=None):
	s = "WARNING: " + warning + ": "
	return(add_arg_values(s, argnames, args))

def regularError(msg=""):
	return("ERROR: " + msg + ".")

## Full: get list of binary integers    ##
## corresponding to input Bit-Vector r  ##
## Full iff. a size is given            ##
## Not Full: get integer corresponding  ##
## to input Bit-Vector r                ##
def get_binary_dec(r, size=None):
	n = r.size()
	i = int(str(simplify(BV2Int(r))))
        if (size):
		options = '0' + str(int(size)) + 'b'
		return(concat(reversed(list(format(i, options)))))
	else:
		return(str(i))

#######################
## STRING PROCESSING ##
#######################

## Concatenate strings                                    ##
def concat(x, sep=""):
	return(reduce(lambda a,b : a+sep+b, map(str, x)))

## Sanitize strings using starting and ending characters  ##
get_element_between = lambda x, b, e : concat(x.split(b)[1].split(e)[0].split(" "))

## Sanitize strings using starting character              ##
get_element_after = lambda x, b : concat(x.split(b)[1].split(" "))

## Sanitize strings using ending character              ##
get_element_before = lambda x, b : concat(x.split(b)[0].split(" "))

## Implement custom grep                                  ##
grep = lambda x, pattern : len(x.split(pattern)) > 1

## Count number of occurrences                            ##
count = lambda x, pattern : len(x.split(pattern))-1

#######################
## LIST PROCESSING   ##
#######################

filter_nones = lambda x : filter(lambda y : y!=None, x)

################################
## BIT-VECTORS BUILDING       ##
################################

#' Construct Bit-Vector encoding 2^i of given size
#'
#' @param i log2 of value 
#'          (integer)
#' @param s size of the vector
#'          (positive integer)
#' @return b Z3 Bit-Vector encoding 2^i of size s 
def bv(i, s):
	return(BitVecVal(2**i, s))

#' Construct Bit-Vector 0,0,0,0,...,0
#' of size l
#'
#' @param l size of the vector
#'          (positive integer)
#' @return b zero Z3 Bit-Vector of size l 
def buildZERO(l):
	if (not l):
		print(noBV("ZERO", ["l"], [l]))
		return(None)
	return(BitVecVal(0, l))

#' Construct Bit-Vector 1,1,1,1,...,1
#' of size l
#'
#' @param l size of the vector
#'          (integer)
#' @return b full-one Z3 Bit-Vector of size l 
def buildONE(l):
	if (not l):
		print(noBV("ONE", ["l"], [l]))
		return(None)
	return(BitVecVal(2**l-1, l))

#' Construct Bit-Vector with the 
#' coordinate i removed: 
#' b1, ..., bi-1, bi+1, ..., bn
#'
#' @param b Z3 Bit-Vector
#' @param i index of the coordinate to remove
#'          (integer)
#' @return res Z3 Bit-Vector minus the i^th coordinate
def trim1BV(b, i):
	n = b.size()
	if (i == None or i >= n or i < 0 or n==0):
		print(noBV("TRIMMED", ["i", "n"], [i, n]))
		return(None)
	if (i==0 and n==1):
		print(noBV("EMPTY", ["i", "n"], [i, n]))
		return(None)
	if (i==0):
		return(Extract(n-1, 1, b))
	if (i==n-1):
		return(Extract(n-2, 0, b))
	else:
		return(Concat(Extract(n-1, i+1, b), Extract(i-1, 0, b)))

#' Construct Bit-Vector with only the 
#' coordinate i: bi
#'
#' @param b Z3 Bit-Vector
#' @param i index of the coordinate
#'          (integer)
#' @return res Z3 Bit-Vector of size 1 with i^th coordinate
def extract1BV(b, i):
	n = b.size()
	if (i == None or i >= n or i < 0 or n==0):
		print(noBV("1-SIZED", ["i", "n", "b"], [i, n, b]))
		return(None)
	return(Extract(i, i, b))

#' Construct constant vector using the list 
#' of non-negative bit indices
#'
#' @param ls list of non-negative bit indices
#'           (integer list)
#' @param size size of the final Bit-Vector
#' @return res Z3 Bit-Vector of size size
def idxList2BV(ls, size):
	if (not(ls)):
		return(buildZERO(size))
	if (any([l >= size for l in ls]) or any([l==None for l in ls])):
		print(warnBV("Bit-Vector will be trimmed", ["l", "size"], [l, size]))
	return(BitVecVal(sum([2**i for i in ls]), size))

false = buildZERO(1)
true = buildONE(1)

#' Count the number of 1-bits in Bit-Vector
#'
#' @param q Z3 Bit-Vector
#' @return res Z3 Bit-Vector of the same size which value
#'     is equal to the number of 1's in q
def count_1s(q):
	n = q.size()
	if (n == 1):
		return(q)
	bits = [extract1BV(q, i) for i in range(n)]
	bvs = [Concat(buildZERO(n-1), b) for b in bits]
	res = reduce(lambda a, c: a+c, bvs)
	return(res)	

################################
## OPERATIONS ON BOOLEANS     ##
################################

map_and = lambda ls : reduce(lambda x,y: And(x, y), ls)
map_or = lambda ls : reduce(lambda x,y: Or(x, y), ls)
map_bool_and = lambda ls : reduce(lambda x,y: x and y, ls)

################################
## BIT-VECTORS COMPARISON     ##
################################

#' Test q[c] == q1[c1], where q, q1 are Bit-Vectors
#'
#' @param q Z3 Bit-Vector
#' @param q1 Z3 Bit-Vector
#' @param c index of the bit to be tested in q
#'          (integer)
#' @param c1 index of the bit to be tested in q1
#'          (integer)
#' @return cond condition to add to solver  
def cond_value_vec(q, q1, c, c1=None):
	if (c1 == None):
		c1 = c
	qq = extract1BV(q, c)
	qq1 = extract1BV(q1, c1)
	if (qq == None or qq1 == None):
		print(warnBV("One of the vectors to compare is NULL", ["qq", "qq1"], [qq, qq1]))
		return(False)
	return(qq==qq1)

#' Test q[i] == v, where q is a Bit-Vector
#'
#' @param q Z3 Bit-Vector
#' @param i index of the bit to be tested
#'          (integer)
#' @param v value 
#'          (Z3 bit-vector of size 1)
#' @return res condition to add to solver 
def cond_value(q, i, v):
	if (v==None):
		print(noBV("NULL BIT-VECTOR", ["q", "i"], [q, i]))
		return(None)
	if (v.size()!=1):
		print(noBV("WRONG", ["|v|", "q", "i"], [v.size(), q, i]))
		return(None)		
	qq = extract1BV(q, i)
	if (qq == None):
		print(regularError("NULL extracted Bit-Vector"))
		return(None)
	return(qq==v)

#' Test q[i] == ls[i], for i < q.size(), where q is a Bit-Vector
#' and ls a list of binary values
#'
#' @param q Z3 Bit-Vector
#' @param ls list of binary values
#' @return cond condition to add to solver 
def cond_equal_value(q, ls):
	vector = idxList2BV(filter_nones([ifthenelse(int(ls[i]), i) for i in range(len(ls))]), q.size())
	cond_ls = [cond_value(q, i, ifthenelse(int(ls[i]), true, false)) for i in range(len(ls))]
	return(map_and(cond_ls))
	
#' Test q[idx!=i] == q1[idx!=i]
#' where q, q1 are Bit-Vectors
#'
#' @param q Z3 Bit-Vector
#' @param q1 Z3 Bit-Vector
#' @param i index of the bit to be ruled out
#'          (integer)
#' @return res condition to add to solver 
def equal_except1BV(q, q1, i):
	qq1 = trim1BV(q1, i)
	qq = trim1BV(q, i)
	if (qq1==None or qq==None):
		print(warnBV("One of the vectors to compare is NULL", ["qq", "qq1"], [qq, qq1]))
		return(False)
	return(qq1==qq)

############################################
## MULTI-LEVEL TO BINARY LEVEL MAPPING    ##
############################################
#' These are functions related to conversion of multi-level GRN matrices
#' to lists of binary bit-vector GRNs
#' that is, GRN vectors that will be used to update system states
#' A GRN vector is a binary vector associated with a target gene 
#' (if changes are allowed, also with a time step and field)
#' of length #possibilities 
#' (= |{values for gene #1}| x |{values for gene #2}| x ...)
#' = #values for the vector of binary variables associated with each gene
#' If the ith value of this vector is set to 0 (resp. to 1), 
#' it means that GRF_{target gene}(vector #i) = 0 (resp. = 1)
#' The numbering of Boolean vectors is fixed beforehand 
#' (using gray code, see formula_to_truthtable.py) 

#' @param genes list of (gene, #levels for this gene) pairs
#' @return multi_binary_dict dictionary which contains the full mapping from binary
#'      variables to real genes and the way back:
#' "vgenes": list of (real gene, #levels for this gene) pairs
#' "genes": list of binary variables associated with each concentration level of each gene
#' "map_to_vgenes": list of (real gene identifier, concentration level value) of size #genes
#' "mapping": dictionary of (real gene name, unique integer identifier)
#' "map_list": list of lists of same-gene binary variable indices
def multi_to_binary(genes):
	## List of (gene, #levels) pairs
	vgenes = genes
	## For gene G with n levels, we build n variables G_i
	## such as G_i = 1 iff. level_i < level(G) <= level_i+1 (1 <= i <= n, level_n+1 = +inf)
	## List of equivalent binary gene variable names
	genes = [vgenes[idx_gene][0]+"-"+str(value+1) for idx_gene in range(len(vgenes)) for value in range(vgenes[idx_gene][1])]
	## List of (integer gene identifier, level no.) pairs for each binary variable
	## (|map_to_vgenes| == |genes|)
	map_to_vgenes = [[idx_gene, value+1] for idx_gene in range(len(vgenes)) for value in range(vgenes[idx_gene][1])]
	## Build dictionnaries for mapping ##
	## Link gene integer identifier to gene name
	mapping = dict()
	for i in range(len(vgenes)):
		mapping.setdefault(vgenes[i][0], i)
	## Find all indices of concentration levels/binary variables for this gene         ##
	map_list = [filter_nones([ifthenelse(map_to_vgenes[j][0]==i, j) for j in range(len(genes))]) for i in range(len(vgenes))]
	multi_binary_dict = dict()
	multi_binary_dict.setdefault("vgenes", vgenes)
	multi_binary_dict.setdefault("genes", genes)
	multi_binary_dict.setdefault("map_to_vgenes", map_to_vgenes)
	multi_binary_dict.setdefault("mapping", mapping)
	multi_binary_dict.setdefault("map_list", map_list)
	return(multi_binary_dict)

#' @param vectors list of integer lists (enumeration of all possible gene values)
#' @param m number of binary variables
#' @param multi_binary_dict dictionary which contains the full mapping from binary
#'      variables to real genes and the way back
#' @param verbose logical for printing messages
#' @param debug logical for printing (hopefully) useful insights to the code
#' @return bvectors conversion of every vector in vectors into a Boolean vector
#'         It is the equivalent enumeration of all possible binary variable values
def multi_to_boolean(vectors, m, multi_binary_dict, verbose=False, debug=False):
	if (verbose):
		print("\n--- Conversion from multi-level vectors to Boolean vectors")
	bvectors = []
	## For every integer vector
	for v in vectors:
		## Equivalent Boolean vector initialization
		## of length #(binary variables)
		bvector = [0]*m
		for i in range(len(v)):
				gene = multi_binary_dict.get("mapping")[multi_binary_dict.get("vgenes")[i][0]]
				ids = multi_binary_dict.get("map_list")[gene]
				if (debug):
					print("Coordinate #" + str(i))
					print("-- Gene " + multi_binary_dict.get("vgenes")[i][0] 
						+ " = " + str(v[i]))
					print("Values for > conc. levels: " 
						+ str([multi_binary_dict.get("genes")[idx] 
						+ " = 0" for idx in ids[v[i]:]]))
					print("Values for <= conc. levels: " 
						+ str([multi_binary_dict.get("genes")[idx] 
						+ " = 1" for idx in ids[:v[i]]]))
				for idx in ids[:v[i]]:
					bvector[idx] = 1
		bvectors += [bvector]
	## Compare with another less elegant method
	if (debug):
		print("\n--- Testing conversion to Boolean vectors")
		print("Boolean" + " "*(m-len("Boolean")+3+1) + "Multi-Level")
		for i in range(len(vectors)):
			b1 = concat(bvectors[i])
			vu = concat(vectors[i])
			v1 = list(vu)
			ids = [idx for ls in multi_binary_dict.get("map_list") for idx in ls[1:]]
			for j in range(len(v1)):
				if (v1[j] > 1):
					char = ["0"]*len(multi_binary_dict.get("map_list")[j])
					for idx in range(int(v1[j])):
						char[idx] = "1"
					v1[j] = concat(char)
			v1 = concat(v1)
			print(b1 + " == " + v1 + " (" + vu + ") : " + str(b1 == v1))
	return(bvectors)

#' @param map_to_vgenes array of pairs (real gene identifier, concentration level 
#'           associated to the binary variable) of length m = #binary variable
#' @param grn_value value of gene expression associated with the output for binary variable #j
#'                      for input bvector #i_x
#' @param i_x index of the input Boolean vector
#' @param j index of the binary variable
#' @return i_x if the associated binary variable (#j) is set to 1 with Boolean vector #i_x, otherwise None
def to_binary(map_to_vgenes, grn_value, j, i_x):
	## If the value of the output of multi-level GRN for input bvector #i_x for binary variable j
	## is greater or equal to the associated concentration level threshold
	## Then the corresponding output in the boolean GRN is set to 1
	return(ifthenelse(grn_value >= map_to_vgenes[j][1], i_x))

#' @param GRNs dictionary of GRNs (as described in the observations file)
#' @param map_to_vgenes array of pairs (real gene identifier, concentration level 
#'           associated to the binary variable) of length m = #binary variable
#' @param idx_grn index of the considered GRN in GRNs
#' @param j index of the considered binary variable
#' @return column column of the GRN associated with gene related to binary variable (which is a given
#' concentration level of a real gene) returned as a list of integers
def col_grn(GRNs, map_to_vgenes, idx_grn, j): 
	## GRNs.items() returns a list of pairs (key, item) in dictionary
	## The #idx_grn one is selected, and we are only interested in the item
	grn = GRNs.items()[idx_grn][1]
	## Get the gene integer identifier associated with binary variable #j
	gene = map_to_vgenes[j][0]
	## Return column of GRN associated with considered gene (of length #boolean vectors)
	return(map(lambda x: int(float(x)), list(grn[0:, gene])))

#' @param GRNs dictionary of GRNs (as described in the observations file)
#' @param m number of binary variables
#' @param npossibilities number of Boolean vector values for the m binary variables
#' @param multi_binary_dict dictionary which contains the full mapping from binary
#'      variables to real genes and the way back
#' @param verbose logical for printing messages
#' @param debug logical for printing (hopefully) useful insights to the code
#' @return grn_vector build (variable) Z3 Bit-Vectors associated with GRN of every gene
def compute_grn_bvs(GRNs, m, npossibilities, multi_binary_dict, verbose=False, debug=False):
	map_to_vgenes = multi_binary_dict["map_to_vgenes"]
	vgenes = multi_binary_dict["vgenes"]
	genes = multi_binary_dict["genes"]
	## Build lists of positive coefficient position in each Bit-Vector 
	## associated to every GRN and every binary variable
	if (not debug):
		idx_lists = [
			[
				filter_nones([to_binary(map_to_vgenes, col_grn(GRNs, map_to_vgenes, idx_grn, j)[i_x], j, i_x) 
						## Check all indices of bvectors 
						for i_x in range(npossibilities)
					## For every binary variable 
					## (we want to convert #genes GRN vectors into m GRN vectors)
					]) for j in range(m)] 
				## For every GRN matrix
				for idx_grn in range(len(GRNs))
		]
	else:
		idx_lists = []
		print("\n--- Correct mapping between binary and multi-level variables:")
		## Same code, but unfold
		for idx_grn in range(len(GRNs)):
			idx_grn_ls = []
			for j in range(m):
				idx_list = []
				for i_x in range(npossibilities):
					## Extract column of GRN matrix #idx_grn associated with gene 
					## binary variable #j is a concentration level of this gene
					grn = col_grn(GRNs, map_to_vgenes, idx_grn, j)
					print("grn (\'" + GRNs.keys()[idx_grn] 
						+ "\') #" + str(idx_grn+1) + "/" 
						+ str(len(GRNs)) + " gene (" 
						+ genes[j] + " - " 
						+ vgenes[map_to_vgenes[j][0]][0] + " - " 
						+ str(map_to_vgenes[j][1]) + ") #" + str(j+1) 
						+ "/" + str(m) + " vector #" + str(i_x+1) 
						+ "/" + str(npossibilities) + "=" + str(len(grn)) 
						+ " ? " + str(genes[j] == vgenes[map_to_vgenes[j][0]][0]
						+ "-" + str(map_to_vgenes[j][1])) + " == True")
					idx_list += [to_binary(map_to_vgenes, grn[i_x], j, i_x)]
				idx_grn_ls += [filter_nones(idx_list)]
			idx_lists += [idx_grn_ls]
		## Checking results
		print("\n--- Checking positive indices in associated Bit-Vectors:")
		for i in range(len(idx_lists[0])):
			print("** Gene " + genes[i] + " : " + str(idx_lists[0][i]))
			if (i > 0 and map_to_vgenes[i][1] > 1):
				print("Included in previous one? = " + str(all([x in idx_lists[0][i-1] for x in idx_lists[0][i]])))
	if (debug):
		for idx_grn in range(len(GRNs)):
			print("\n" + " "*len("Max. Values for each gene: ") + concat([x[0] for x in vgenes], ", "))
			print("Max. Values for each gene: " + concat([x[1] for x in vgenes], ", "))
			print("#vectors = " + str(npossibilities))
			for j in range(m):
				print(genes[j] + ": " + concat(idx_lists[idx_grn][j], " "))
				bvv = get_binary_dec(idxList2BV(idx_lists[idx_grn][j], 
					npossibilities), npossibilities)
				bvvect = [bvv[i] == ifthenelse(i in idx_lists[idx_grn][j], "1", "0") 
					for i in range(npossibilities)]
				print(str(all(bvvect)) + " == True: #1's = " + str(count(bvv, "1")) 
					+ " == " + str(len(idx_lists[idx_grn][j])))
	## Build the associated Bit-Vectors
	grns_bv = [[idxList2BV(idx_lists[idx_grn][j], npossibilities) for j in range(m)] for idx_grn in range(len(GRNs))]
	return(grns_bv)

#' @param s Z3 solver object
#' @param Observations list of {name, time step, corresponding GRN set, [phenotypes]}
#' @param GRNs dictionary of GRNs (as described in the observations file)
#' @param m number of binary variables
#' @param k maximum length of state trajectory
#' @param nexp number of observations
#' @param nfields number of fields
#' @param multi_binary_dict dictionary which contains the full mapping from binary
#'      variables to real genes and the way back
#' @param exp_names list of character strings (experiment names)
#' @param npossibilities number of Boolean vector values for the m binary variables
#' @param verbose logical for printing messages
#' @param debug logical for printing (hopefully) useful insights to the code
#' @return grn_vector build (variable) Z3 Bit-Vectors associated with GRN of every gene
def convert_grn_bv(s, Observations, GRNs, m, k, nexp, nfields, multi_binary_dict, exp_names, npossibilities, verbose=False, debug=False):
	## grns_bv is a list of length #experiments of lists of length #steps
	## of lists of length #fields of lists of bit-vectors of length #genes
	grn_name = lambda e,t,nf,j : 'GRN^(' + str(e) + ',t=' + str(t) + ')_(' + str(nf) + ',' + str(j) + ')'
	grns_bv = [[[[
			BitVec(grn_name(e, t, nf, j), npossibilities) 
		for j in range(m)] 
		for nf in range(nfields) ] 
		for t in range(k) ]
		for e in range(nexp) ]
	if (verbose):
		print("\n--- Condition on GRNs (variable)")
	grns = compute_grn_bvs(GRNs, m, npossibilities, multi_binary_dict, verbose=verbose, debug=debug)
	## Get constraints on GRNs in observations
	obs = filter(lambda x:not(x["GRN"]==None), Observations)
	cond_ls = [grns_bv[exp_names.index(oo["name"])][oo["step"]][oo["field"]][j] == grns[GRNs.keys().index(oo["GRN"])][j] 
		for oo in obs for j in range(m)]
	if (debug):
		print("")
		for oo_i in range(len(obs)):
			print("* Observation = #" + str(exp_names.index(obs[oo_i]["name"])) 
				+ " step=" + str(obs[oo_i]["step"]) + " field=" + str(obs[oo_i]["field"]))
			for j in range(m):
				print(multi_binary_dict["genes"][j] + ": " + str(cond_ls[oo_i*m+j]))
	cond = map_and(cond_ls)
	s.add(cond)
	return([s, grns_bv])

#' @param GRNs dictionary of GRNs (as described in the observations file)
#' @param m number of binary variables
#' @param multi_binary_dict dictionary which contains the full mapping from binary
#'      variables to real genes and the way back
#' @param npossibilities number of Boolean vector values for the m binary variables
#' @param verbose logical for printing messages
#' @param debug logical for printing (hopefully) useful insights to the code
#' @return grn_vector build (constant) Z3 Bit-Vectors associated with GRN of every gene
def convert_grn(GRNs, m, multi_binary_dict, npossibilities, verbose=False, debug=False):
	if (verbose):
		print("\n--- Condition on GRNs (constants)")
	grns_bv = compute_grn_bvs(GRNs, m, npossibilities, multi_binary_dict, verbose=verbose, debug=debug)
	if (debug):
		print("\n--- Testing vector equality:")
		for i_grn in range(len(grns_bv)):
			for i_col in range(len(grns_bv[i_grn])):
				## From the newly-built Bit-Vector
				col = grns_bv[i_grn][i_col]
				res_col = get_binary_dec(col, col.size())
				## "Ground truth"
				grn = col_grn(GRNs, multi_binary_dict.get("map_to_vgenes"), i_grn, i_col)
				init_col = concat([int(x >= multi_binary_dict.get("map_to_vgenes")[i_col][1]) for x in grn])
				print("** Gene " + multi_binary_dict.get("genes")[i_col] + " in GRN \'" 
					+ GRNs.keys()[i_grn] + "\': " 
					+ str(res_col == init_col) + " == True")
				print("|BV| = " + str(col.size()) + " == " + str(npossibilities) 
					+ " (size of vector space)" + " == |TT| = " + str(len(init_col)))
				print("Bit-Vector: " + res_col)
				print("Truth Table:" + init_col)
	return(grns_bv)

#' Convert binary -single- trajectory matrix into multi-level trajectory matrix
#'
#' @param trajectory binary matrix of size 
#'        (#steps x #fields x 1 (#trajectories)) x #binary variables
#' @param k number of time steps
#' @param multi_binary_dict dictionary which contains the full mapping from binary
#'      variables to real genes and the way back
#' @param verbose logical for printing messages
#' @param debug logical for printing (hopefully) useful insights to the code
#' @return trajectories_multi integer trajectories matrix
def binary_to_multi(trajectory, k, multi_binary_dict, verbose=False, debug=False):
	nrow, ncol = np.shape(trajectory)
	if (verbose):
		print("\n--- Trajectory matrix")
		print("Dimensions: #steps x #nfields = " + str(nrow) + " x #(binary variables) = " + str(ncol))
	vgenes = [g[0] for g in multi_binary_dict.get("vgenes")]
	nfields = int(nrow/(k+1))
	print("\n#fields == " + str(nfields))
	## Initialize multi-level matrix
	trajectory_multi = np.zeros(((k+1)*nfields, len(vgenes)))
	if (verbose):
		nrow1, ncol1 = np.shape(trajectory_multi)
		print("\n--- Multi-level trajectory matrix")
		print("Dimensions: #steps x #nfields = " + str(nrow1) + " x #(gene variables) = " + str(ncol1))
	for t in range(k+1):
		if (verbose):
			print("-- Time step t=" + str(t))
		for n in range(nfields):
			## Get system states for time step = t, field = n, trajectory # = exp
			binary = map(int, trajectory[t*nfields+n, 0:])
			## Convert binary variables into real gene multi-level variables
			multi = [0]*len(vgenes)
			for i in range(len(binary)):
				if (int(binary[i])):
					idx, value = multi_binary_dict.get("map_to_vgenes")[i]
					multi[idx] = max(multi[idx], value)
			if (verbose):
				print("\n--- FIELD " + str(n+1))
				print(concat([g+" "*max(0, 2-len(g[0])) for g in vgenes], " | "))
				print(concat([str(multi[i])+" "*(len(vgenes[i])-1) for i in range(len(multi))], "  | "))
				print("---")
			if (debug):
				print(concat(multi_binary_dict.get("genes"), " |"))
				print(concat([str(binary[i])+" "*len(multi_binary_dict.get("genes")[i]) for i in range(len(binary))], "|"))
				print(concat(multi))
				print("Dimensions - #genes: " + str(len(multi)) + " == " + str(len(vgenes)) 
					+ ", - #steps: " + str(nfields*t+n) + " < " + str(nfields*(k+1)))
			## State corresponding to field #n and time step t
			trajectory_multi[nfields*t+n, 0:] = multi
	if (debug):
		print("\n-- Final Multi-Level Matrix")
		print("| #steps -> #genes")
		print(trajectory_multi)
	return(trajectory_multi)

##################################
## SAVE RESULTS AS CSV TABLES   ##
##################################

#' Save results as CSV tables
#'
#' @param result matrix to save as CSV
#' @param no solution identifier
#' @param modelname name of model
#' @param filetype (in ["pattern_matrix", "trajectories"])
#' @param patterns list of patterning function (described as dictionary of features)
#' @param multi_binary_dict dictionary which contains the full mapping from binary
#'      variables to real genes and the way back
#' @param bvectors enumeration of Boolean vectors for binary variables
#' @return None (converts and writes result into a proper CSV file, that can be loaded 
#'         using 'numpy.loadtxt' function)
def save_as_csv(result, no, modelname, filetype, patterns, multi_binary_dict, bvectors):
	from subprocess import call
	import os
	from numpy import savetxt, shape, reshape
	path = path_to_results + modelname + "/"
	if (not os.path.isdir("./" + path)):
		call("mkdir " + path, shell=True)
	## File name for the resulting CSV file
	filename = path + "result_" + modelname + "_" + str(no+1) + "_" + filetype + ".csv"
	## - If result is a matrix of trajectories
	if (filetype == "trajectories"):
		## Cut the matrix into nexp x k smaller matrices of size nfields x m
		dims = shape(result)
		f_handle = file(filename, "a")
		for i in range(dims[0]):
			for j in range(dims[1]):
				## File object
				savetxt(f_handle, 
					## Resized smaller matrix
					reshape(result[i, j], (dims[2], dims[3])).astype(int), 
					delimiter=",", fmt='%i', newline="\n", 
					header="Model #" + str(no) + " for " + modelname 
						+ "\n\nTrajectory #" + str(i) + ", Time step #" + str(j) + "\n\n" 
						+ concat(multi_binary_dict.get("genes"), ","), 
					footer="\n " + str(dims[2]) + " fields.")
				f_handle.write("\n\n")
		f_handle.close()
	## - If result is a matrix of pattern selection vectors
	elif (filetype == "pattern_matrix"):
		## Save it as it is, provided some metadata
		savetxt(filename, result.astype(int), delimiter=",", fmt='%i', newline="\n", 
			header="Model #" + str(no) + " for " + modelname + "\n\nPattern matrix\n\nPatterns = \n" 
				+ concat(map(lambda i: "Pattern #" + str(i) + ": " 
					+ str(patterns[i]), range(len(patterns))), "\n")
				+ "\n\n"
				+ concat(range(shape(result)[1]), ","))
	elif (filetype == "grns" and not(result==None)):
		## Cut the matrix into nexp x k x nfields smaller matrices of size #bvectors x m
		nexp, k, nfields, lenbvectors, m = shape(result)
		f_handle = file(filename, "a")
		for i in range(nexp):
			for j in range(k):
				for l in range(nfields):
					## File object
					savetxt(f_handle, 
						## Resized smaller matrix
						reshape(result[i, j, l, 0:, 0:], (lenbvectors, m)).astype(int), 
						delimiter=",", fmt='%i', newline="\n", 
						header="Model #" + str(no) + " for " + modelname 
							+ "\n\nTrajectory #" + str(i) + ", Time step #" + str(j) 
							+ ", Field #" + str(l+1) + "\n\n" 
							+ concat(multi_binary_dict.get("genes"), ","), 
						footer="\n " + str(nfields) + " fields.\n\n" 
						+ "Boolean vectors per line:\n"
						+ concat(multi_binary_dict.get("genes"), ",") + "\n"
						+ concat([concat(bvectors[b], ",") for b in range(lenbvectors)], "\n"))
					f_handle.write("\n\n")
		f_handle.close()
	return(None)
