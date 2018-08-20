# -*- coding: utf-8 -*-

############
## PATHS  ##
############

## Relative to the current solver/ folder
path_to_models="../models/"
path_to_results="results/"

#################
## PARAMETERS  ##
#################

## Type of conditions for unique solutions
#uniqueness = ["patterning", "states", "grns"]

## Diffusion rate for gene products
## 0 <= diffusion_rate <= 1
## diffusion_rate = 0 <=> No diffusion between neighbouring fields
## diffusion_rate = 1 <=> 100% of element concentration 
## (not taking into account the decrease in concentration due to distance to the source)
## is diffused to the neighbouring fields

## change_grn allows GRNs to change depending on time/position
change_grn = False

#---------------------------------

## See dictionary 'diffusion' in patterns.py
# Function to implement morphogen/gene product diffusion
diffusion_function = "version1"

## See dictionary 'aggregation' in patterns.py
# Function to aggregate parallel patternings
aggregation_function = "logical_sum"

## See dictionary 'application' in patterns.py
# Function to apply patterning changes to a state
application_function = "logical_sum"

## List of colors (must be at least of length #putative morphogens)
colors = ["yellow", "green", "blue", "red", "orange", "grey", "cyan", "pink"]

## Parameters for patterning function selection
## If their value should not be restricted, write 0
# Maximum number of selected patterning functions per level (= time step)
max_nb_patterns_per_level = 0
# Maximum (resp. min) number of selected patterning functions among those provided
max_nb_patterns = 0
min_nb_patterns = 0
# Maximum of times any patterning function is selected (in a row)
max_nb_pattern_times = 1

############
## MODELS ##
############

## In the code

default_model="model.net"
default_observations="observations.spec"
toy_model = "toy/" + default_model
toy_observations = "toy/" + default_observations
gapgene_model = "gap gene sanchez/" + default_model
gapgene_observations = "gap gene sanchez/" + default_observations

## For tests

tests_model = "toy/"
prediction_model = "toy"
prediction_model_id = 0
prediction_initial_phenotypes = "InitialPhenotypes-InitialPhenotypes-InitialPhenotypes"
prediction_grns = "InitialGRNs1-InitialGRNs2-InitialGRNs3"
prediction_final_phenotypes = "FinalPhenotypes1-FinalPhenotypes2-FinalPhenotypes3"
