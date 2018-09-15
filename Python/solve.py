# -*- coding: utf-8 -*-

import sys
from z3 import *

from globals import *
from pattern_solver import *
from models import *
from utils import save_as_csv, concat, grep
from visualize import model_to_igraph, model_to_plot

######################################################
## CONVERT TRANEJCTORY MATRIX TO PLOTTABLE MATRIX   ##
######################################################

def to_plottable_matrix(trajectories, model_id, modelname, multi_binary_dict, verbose=False):
	from numpy import shape, zeros, reshape
	## Associated trajectory plots
	nexp, k, nf, m = shape(trajectories)
	for exp in range(nexp):
		trajectory_flat = zeros((k*nf, m))
		for gene_id in range(m):
			trajectory_flat[0:, gene_id] = list(reshape(trajectories[exp, 0:, 0:, gene_id], (k*nf, 1)))
		model_to_plot(trajectory_flat, k-1, exp, model_id, modelname, multi_binary_dict, verbose=verbose)
	print("\nMSG: Plot for solution #" + str(model_id+1) + " of model " + modelname)
	return(None)

##############################################
## SOLVE MODEL call                         ##
##############################################

#________________________#
#      RUN               #
#________________________#

## <!> Run 'clean.sh' beforehand!
def call_run(modelname="toy", plot_trajectories=False, verbose=True):
	names = map(lambda x: modelname+"/"+x, [default_model, default_observations])
	params_in = read_files(names[0], names[1])
	try:
		res = pattern_solver(params_in, verbose=verbose)
	except ValueError:
		print("ERROR: Something occurred during execution")
		return(None)
	if (not res[0]):
		return(None)
	## Dump pattern + trajectory matrices as CSV files
	files = ["pattern_matrix", "trajectories", "grns"]
	## Avoid special characters such as whitespaces
	modelname = concat(modelname.split(" "), "-")
	diffusion_rate = params_in[3]["diffusion-rate"]
	for model_id in range(len(res[0])):
		pattern_matrix, trajectories, grns = res[0][model_id]
		patterns = res[1]
		multi_binary_dict = res[2]
		bvectors = res[3]
		for j in range(len(files)):
			print("\nMSG: Saving \'" + files[j] + "\' object from solution #" + str(model_id+1) + " of model " + modelname)
			save_as_csv(res[0][model_id][j], model_id, modelname, files[j], patterns, multi_binary_dict, bvectors)
		## Visualization of the solution
		model_to_igraph(pattern_matrix, model_id, modelname, patterns, verbose=verbose)
		if (plot_trajectories):
			to_plottable_matrix(trajectories, model_id, modelname, multi_binary_dict)
	print("\nMSG: Results saved as csv files!")
	return(None)

#________________________#
#      PREDICT           #
#________________________#

## <!> Run 'clean.sh' beforehand!
## model_id < 0 means testing the all-zero pattern matrix (no patterning function selected)
def call_predict(modelname, model_id, q0, grn, qf=None, solmax=None, plot_trajectories=False, verbose=True):
	print("\n-----------\nMODEL = " + modelname)
	print("Solution " + ifthenelse(model_id<0, "\'no pattern selection\'", "#" + str(model_id+1)) + "\n-----------")
	from numpy import shape, zeros, reshape
	names = map(lambda x: modelname+"/"+x, [default_model, default_observations])
	params_in = read_files(names[0], names[1], return_phenotypes=True)
	## Delete previous observations
	del params_in[-3]
	## Delete previous fix points
	del params_in[-3]
	nfields = len(params_in[4])
	if (not(nfields == len(q0))):
		print("ERROR: Wrong initial condition vector length: " + str(len(q0)) 
			+ " instead of " + str(nfields))
		return(None)
	if (not(nfields == len(grn))):
		print("ERROR: Wrong initial GRN vector length: " + str(len(q0)) 
			+ " instead of " + str(nfields))
		return(None)
	Patterns = params_in[-1]
	patterning_step = params_in[2]["patterning_step"]
	k = params_in[2]["nsteps"]
	nsteps = ifthenelse(patterning_step==0, k, min(patterning_step, k))
	## Add observations for initial states
	Observations = filter_nones([ifthenelse(len(q0[id_field]) > 0,
			{"name": "Trajectory", 
				"step": 0, 
				"field": id_field, 
				## name of GRN instead of GRN itself: GRNs.get(di.get("GRN"))
				"GRN": grn[id_field],
				"phenotype": Patterns.get(q0[id_field])
			}) for id_field in range(nfields)])
	## Add observations for final states
	if (qf):
		Observations += filter_nones([ifthenelse(len(qf[id_field]) > 0,
			{"name": "Trajectory", 
				"step": k, 
				"field": id_field, 
				## name of GRN instead of GRN itself: GRNs.get(di.get("GRN"))
				"GRN": None,
				"phenotype": Patterns.get(qf[id_field])
			}) for id_field in range(nfields)])
	params_in = params_in[:-2] + [Observations, []] + [params_in[-2]]
	diffusion_rate = params_in[3]["diffusion-rate"]
	if (model_id < 0):
		pattern_matrix = zeros((len(params_in[1]), nsteps))
	else:
		filename = path_to_results + modelname + "_rate=" + str(diffusion_rate) + "/"
		filename += "result_" + concat(modelname.split("-")) + "_" + str(model_id+1) + "_pattern_matrix.csv"
		pattern_matrix = np.loadtxt(open(filename, "rb"), delimiter=",", skiprows=1)
	print("\n* Pattern matrix:")
	print("t= " + concat([str(i)+" "*(3-len(str(i))+1) for i in range(nsteps)]))
	print(pattern_matrix)
	print("\n--- Starting predictions!")
	try:
		res = pattern_solver(params_in, pattern_matrix0=pattern_matrix, solmax=solmax, verbose=True)
	except ValueError:
		print("ERROR: Something occurred during execution")
	 	return(None)
	if (not res[0]):
		return(None)
	## Avoid special characters such as whitespaces
	modelname = "prediction_" + concat(modelname.split(" "), "-") + "_solution=" + str(model_id+1)
	for model_id in range(len(res[0])):
		_, trajectories, grns = res[0][model_id]
		patterns = res[1]
		multi_binary_dict = res[2]
		bvectors = res[3]
		## Dump trajectory matrices as CSV files
		print("\nMSG: Saving \'trajectories\' object from solution #" + str(model_id+1) + " of model " + modelname)
		save_as_csv(res[0][model_id][1], model_id, modelname, "trajectories", patterns, multi_binary_dict, bvectors)
		if (plot_trajectories):
			to_plottable_matrix(trajectories, model_id, modelname, multi_binary_dict)
	print("\nMSG: Results saved as csv files!")
	return(None)

#________________________#
#      INTERFACE         #
#________________________#
 
def printRunSyntaxError(c):
    if (not c):
        print("MSG: If you wanted to run the solver, then you did not use the correct syntax.")
        print("MSG: Correct syntax is \'run [modelname] [--plot]\'.")
        return(True)
    return(False)

def printPredictSyntaxError(c):
    if (not c):
        print("MSG: If you wanted to make predictions with the solver, then you did not use the correct syntax.")
        print("MSG: Correct syntax is:")
	print("\'predict [modelname] [model_id] --q0 [initial system state/condition] --grn [initial GRN conditions] --qf [final system state/condition] --plot [0 or 1]\'.")
	print("(modelname, model_id, q0 and grn are mandatory arguments)")
	print("(Initial (and final) state and GRN names should be defined in the associated \'observations.spec\' file)")
	print("(They should be given for each field separated by \'-\' (null string if no initial condition for a given field)")
        return(True)
    return(False)

def getArgument(x, argv, default):
	xx = "--" + x
	if (any([arg == xx for arg in argv])):
		i = argv.index(xx)
		if (len(argv) < i+2):
			return(default)
		return(argv[i+1])
	return(default)

def build_repeat(argument):
	if (grep(argument, "*")):
		arg, n = argument.split("*")
		return([arg]*int(n))
	return(argument.split("-"))
 
if (len(sys.argv) > 1 and sys.argv[1] == "run"):
	## Run default toy model                               ##
	if (len(sys.argv) == 2):
		call_run()
	elif (len(sys.argv) == 3):
		call_run(sys.argv[2])
	elif (len(sys.argv) == 4 and sys.argv[3] == "--plot"):
		call_run(sys.argv[2], plot_trajectories=True)
	else:
		printRunSyntaxError(False)
elif (len(sys.argv) > 1 and sys.argv[1] == "predict"):
	q0 = build_repeat(getArgument("q0", sys.argv, prediction_initial_phenotypes))
	print(q0)
	qf = build_repeat(getArgument("qf", sys.argv, prediction_final_phenotypes))
	grn = build_repeat(getArgument("grn", sys.argv, prediction_grns))
	solmax = getArgument("solmax", sys.argv, None)
	plot_it = bool(getArgument("plot", sys.argv, ""))
	print("\nq0 = \n" + concat([" field #" + str(i) + ": " + q0[i] for i in range(len(q0))], "\n"))
	print("\nGRN = \n" + concat([" field #" + str(i) + ": " + grn[i] for i in range(len(grn))], "\n"))
	print("\nPlot? = " + str(plot_it))
	## Prediction with first solution of default toy model ##
	if (len(sys.argv) == 2):
		call_predict(prediction_model, prediction_model_id, q0, grn, plot_trajectories=True, solmax=1)
	elif (len(sys.argv) == 6):
		call_predict(sys.argv[2], int(sys.argv[3]), q0, grn)
	elif (len(sys.argv) == 8):
		if ("--solmax" in sys.argv):
			call_predict(sys.argv[2], int(sys.argv[3]), q0, grn, solmax=solmax)
		elif ("--plot" in sys.argv):
			call_predict(sys.argv[2], int(sys.argv[3]), q0, grn, plot_trajectories=plot_it)
		elif ("--qf" in sys.argv):
			call_predict(sys.argv[2], int(sys.argv[3]), q0, grn, qf=qf, plot_trajectories=plot_it)
		else:
			printPredictSyntaxError(False)
	elif (len(sys.argv) == 10):
		call_predict(sys.argv[2], int(sys.argv[3]), q0, grn, plot_trajectories=plot_it, solmax=solmax)
	elif (len(sys.argv) == 12):
		call_predict(sys.argv[2], int(sys.argv[3]), q0, grn, qf=qf, plot_trajectories=plot_it, solmax=solmax)
	else:
		printPredictSyntaxError(False)
else:
	printRunSyntaxError(False)
	printPredictSyntaxError(False)
