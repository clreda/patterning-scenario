# -*- coding: utf-8 -*-

from igraph import *
from numpy import shape
from subprocess import call

from pattern_solver import pattern_solver
from utils import *
from globals import *
 
##############################################
## Display the graph visualization of a     ##
## model found by the solver                ##
##############################################

#' Write PNG file using GraphViz
#'
#' @param g igraph object
#' @param filename character string for file name
#' @param remove_dot logical for removing the DOT file
#' @return None (creates DOT file, converts it to PNG, removes DOT file)
def write_dot_file(g, filename, remove_dot=True):
	from os import remove
	dotfile = filename + ".dot"
	g.write_dot(dotfile)
	call("dot -Tpng " + dotfile + " > " + filename + ".png", shell=True)
	if (remove_dot):
		remove(dotfile)
 
#' Display the "Patterning Network" associated with the considered solution
#'
#' @param pattern_matrix numpy matrix of size #patterns x #steps for patterning function selection
#' @param model_id solution identifier
#' @param modelname name of model
#' @param patterns list of patterns (described as dictionaries of features)
#' @param verbose logical for printing messages
#' @param debug logical for printing (hopefully) useful insights to the code
#' @return res igraph associated to the model
def model_to_igraph(pattern_matrix, model_id, modelname, patterns, verbose=False, debug=False):
	if (verbose):
		print("\n-- Conversion to igraph:")
	if (len(shape(pattern_matrix)) > 1):
		npatt, nsteps = shape(pattern_matrix)
	else:
		npatt, nsteps = shape(pattern_matrix)[0], 1
		pattern_matrix = pattern_matrix.reshape((npatt, nsteps))
	## Get the list of morphogens that appear in provided patterns
	morphogens = list(set([di.get("morphogen") for di in patterns]))
	colours = colors[:len(morphogens)]
	## CREATING NODES ------------------------------ ##
	## Step #0
	Nodes = [["Input"]]
	## Step #1, ... #(nsteps-1) where nsteps is the maximum number
	## of pattern selection steps
	## Get all features of a given pattern (except for "axis" one)
	keys = filter(lambda x: not(x=="axis"), di.keys())
	## Label for "No selected pattern" node at step t
	no_pattern_node = lambda t : "(t="+str(t)+")\nNo pattern"
	## Label for one of the selected patterns
	build_node_features = lambda k, di : ifthenelse(not(k=="morphogen"), k+":", "")+str(di[k])
	pattern_node = lambda t, di : "(t="+str(t)+")\n"+concat([build_node_features(k, di) for k in keys], "\n") 
	## Write nodes associated with each pattern and the "no selected pattern" node for each time step
	Nodes += [ ([pattern_node(t, di) for di in patterns] + [no_pattern_node(t)]) for t in range(nsteps)]
	## Step #nsteps
	Nodes += [["Output"]]
	## CREATING NAMES ------------------------------ ##
	Names = [it for t in range(nsteps+2) for it in Nodes[t]]
	## CREATING COLORS ----------------------------- ##
	Colors = ["white"] 
	for t in range(nsteps):
		Colors += [colours[morphogens.index(di.get("morphogen"))] for di in patterns]+["white"]
	Colors += ["white"]
	if (debug):
		for t in range(nsteps+1):
			print("Step #" + str(t))
			print(Nodes[t])
		print("Nodes:")
		print(str(len(Nodes)) + " == " + str(nsteps+2) + ": " + str(len(Nodes) == nsteps+2) + " == True")
		print("Names:")
		print(str(len(Names)) + " == " + str(2+nsteps*(npatt+1)) 
			+ ": " + str(len(Names) == 2+nsteps*(npatt+1)) + " == True")
		print("Colors:")
		print(str(len(Colors)) + " == " + str(len(Names)) 
			+ ": " + str(len(Names) == len(Colors)) + " == True")
		print("\nNodes: " + str(Names))
		print("\nColors: " + str(Colors))
	## CREATING EDGES ------------------------------ ##
	Edges = []
	current = [0]
	for t in range(1, nsteps+2):
		if (t == nsteps+1):
			selected = [len(Names)-1]
		else:
			## Look at the positive values of the pattern matrix
			## i.e. when a pattern is selected at a given time step
			selected = filter(lambda i:pattern_matrix[i, t-1]==1, range(npatt))
			if (not selected):
				## Select "No pattern" node
				selected = [Names.index(Nodes[t][-1])]
			else:
				## Select patterning node corresponding to the given time step
				selected = [Names.index(Nodes[t][i]) for i in selected]
		Edges += [[current_i, i] for i in selected for current_i in current]
		current = selected
	if (debug):
		print("\nEdges:")
		for i in range(len(Edges)):
			print("--Edge #" + str(i))
			print(str(Names[Edges[i][0]]) + " -> (" + str(Names[Edges[i][1]]) + ")")
	## INITIALIZING ATTRIBUTES --------------------- ##
	vertex_attrs = {"label": Names, 
            "size" : [60+5*(len(c)/2-1) for c in Names]*len(Names),
            "color": Colors}
	edge_attrs = {}
	if (verbose):
		print("\n-- Create graph")
	## CREATING GRAPH ------------------------------ ##
	g = Graph(n = len(Names), edges = Edges, directed = True, vertex_attrs = vertex_attrs, edge_attrs = edge_attrs)
	## Delete 0-degree vertices                      ##
	g.vs.select(_degree = 0).delete()
	if (verbose):
		print("\n-- Create DOT file")
	## Get a GraphViz version                        ##
	import os
	path = path_to_results + modelname + "/"
	if (not os.path.isdir("./" + path)):
		call("mkdir " + path, shell=True)
	filename = path + modelname + "_model_" + str(model_id+1)
	write_dot_file(g, filename)
	return(g)

##############################################
## Display evolution of gene expression     ##
## depending on time and position           ##
##############################################

## Return the path for the concentration plot figures
def get_fig_name(path, exp, t, modelname, model_id):
	return(path + "/" + "plot_exp=" + str(exp) 
		+ "_t=" + str(t) + "_" + modelname 
		+ "_" + str(model_id+1) + ".png")

#' Plot evolution of gene expression depending on 
#' time and position (1 plot per time step -one per trajectory) associated
#' with the considered solution
#'
#' @param trajectory numpy binary matrix of size (k x #fields) x #binary variables for one trajectory
#' @param k number of time steps
#' @param exp_id trajectory integer identifier
#' @param model_id solution identifier
#' @param modelname name of model
#' @param multi_binary_dict dictionary which contains the full mapping from binary
#'      variables to real genes and the way back
#' @param show_plot logical to show plots
#' @param plot_fields logical to plot fields instead of position on the x axis
#' @param verbose logical for printing messages
#' @param debug logical for printing (hopefully) useful insights to the code
#' @return res igraph associated to the model
def model_to_plot(trajectory, k, exp_id, model_id, modelname, multi_binary_dict, show_plot=False, plot_fields=False, verbose=False, debug=False):
	import matplotlib.pyplot as plt
	import os
	## Create folder for images
	path_list = [path_to_results[:-1], modelname, modelname + "_" + str(model_id+1)]
	path = concat(path_list, "/")
	for i in range(len(path_list)):
		ptest = concat(path_list[:i+1], "/")
		if (not os.path.isdir("./" + ptest)):
			call("mkdir " + ptest, shell=True)
	nrow, m = shape(trajectory)
	nfields = int(nrow/(k+1))
	if (verbose):
		print("\n--- Solution #" + str(model_id) + " for model " + modelname)
		print("\n--- Plots for trajectory #" + str(exp_id))
		print("\n-- Binary matrix")
		ls = ["steps", "fields", "genes"]
		ls_values = [k, nfields, m]
		for i in range(len(ls)):
			print("#" + ls[i] + " = " + str(ls_values[i]))
	## Convert binary trajectory matrix into multi-level trajectory matrix
	trajectory_multi = binary_to_multi(trajectory, k, multi_binary_dict, verbose=verbose, debug=debug)
	nrow, m = shape(trajectory_multi)
	if (verbose):
		print("\n---- Plotting the resulting multi-level matrix")
		print("\n-- Multi-level matrix")
		ls = ["steps", "fields", "genes"]
		ls_values = [k, nfields, m]
		for i in range(len(ls)):
			print("#" + ls[i] + " = " + str(ls_values[i]))
	gene_colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
	## Either plot by field (of length 1 on the plot) or either
	## plot with field of length 0 on the plot
	## (meaning that the value of a given gene in field i in a given plot
	## is the value of the corresponding function read at x = i)
	x = ifthenelse(plot_fields, range(0, nfields+1), range(1, nfields+1))
	vgenes = [g[0] for g in multi_binary_dict.get("vgenes")]
	shift=1e-2
	flatten = lambda y : [v2 for v1 in y for v2 in v1]
	if (len(vgenes) > 8):
		print("\nWARNING: Not enough colours to plot all genes!")
		u, q = map(int, [len(vgenes)/8, len(vgenes)%8])
		gene_colors = gene_colors*u+gene_colors[:q]
	print("")
	for t in range(k):
		## y[n][i] = value of gene #i in field #n at time step #t
		## (+ shift in order to make it readable on the plot)
		y = [[trajectory_multi[nfields*t+n, i]+(i*shift) for i in range(m)] for n in range(nfields)]
		if (debug):
			print("\nTrajectory #" + str(exp_id) + ", Step t=" + str(t) + ":")
			print("Genes  |" + concat(vgenes, " "))
			for i in range(nfields):
				print("FLD#" + str(i+1) + "  |" 
				+ concat([str(int(y[i][g]))+" "*(len(vgenes[g])-1) for g in range(m)], " "))
		if (verbose):
			print("\nt=" + str(t))
			print("FLD#  " + concat(map(str, range(1, nfields+1)), " | "))
		for i in range(m):
			if (verbose):
				print(vgenes[i][:5] + " "*(max(1, 5-len(vgenes[i])+1)) 
					+ concat(map(str, [int(yy[i]) for yy in y]), " | "))
			if (plot_fields):
				x_get = flatten([[x[j], x[j+1]] for j in range(nfields)])
				y_get = flatten([[yy[i], yy[i]] for yy in y])
			else:
				x_get = x
				y_get = [yy[i] for yy in y]
			plt.plot(x_get, y_get, gene_colors[i])
		plt.title("Evolution of concentration levels at t=" + str(t) 
			+ " in trajectory #" + str(exp_id) + " = f(position)")
		plt.axis([min(x), max(x), -10*shift, max(flatten(y))+1])
		plt.legend(vgenes)
		plt.xlabel('Position (' 
			+ concat(["field #" + str(n+1) for n in range(nfields)], ", ") + ')')
		plt.ylabel('Concentration level thresholds')
		## Draw the separation between fields
		for nf in range(nfields):
			plt.axvline(x[nf], color="r", linestyle=":")
		## Draw the different concentration level values
		for vg in range(len(flatten(y))):
			plt.axhline(vg, color="g", linestyle=":")
		if (show_plot):
			plt.show()
		filename = get_fig_name(path, exp_id, t, modelname, model_id)
		plt.savefig(filename)
		plt.clf()
		print("\n-- Saved in " + filename + "!")
	return(None)

## TODO visualization of GRNs: conversion of truth table into logical formula (previous R code)
