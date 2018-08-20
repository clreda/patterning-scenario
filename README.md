## Modelling And Prediction of Patterning Phenomena On Drosophila Embryo

This GitHub project provides the code for the Patterning Scenario inference described in the chapter *Modelling And Prediction of Patterning Phenomena On Drosophila Embryo* (non published).

## Requirements & Installation

**Python:** Packages **Z3** (SMT solver) and **igraph** (GRN visualization).

For Debian Linux:

`sudo apt install python-pip`

`sudo python -m pip install z3-solver`

`sudo python -m pip install python-igraph`


## Description

- *models* Example models.

- *Python* Contains the code for the Patterning Scenario inference approach.

- *R* Contains code related to the transformation of results from Boolean network inference (see repository **regulomics/expansion-network** for more information) into the syntax expected in the observations file.


## Patterning Scenario Inference Procedure

Please refer to the report for more information about the implementation. Examples of input files can be found in folder *models*. Drosophila gap-gene segmentation model can be found in the following paper: *Sanchez, L., & Thieffry, D. (2001). A logical analysis of the Drosophila gap-gene system. Journal of theoretical Biology, 211(2), 115-141.*

### Usage

#### Test files

To test function named *function* in the code, type in the terminal (in the "Python" folder):

`python tests.py function`

To see the list of avaiable tests, type the following command:

`python tests.py`


#### Network inference from a model and an experiments file

Tree structure of the model file associated with model named *model*:

```bash
- **/**
-- **models/**
--- **model/**
---- model.net
---- observations.spec
```

Experiments file is *observations.spec*, model file is called *model.net*.

`python solve.py run model [--plot]`


- Option *--plot* plots the gene concentration for each field (per plot), for each time step.

**Example:** To test the **toy** model provided, and plot the protein concentrations, type the following command:

`python solve.py run toy --plot`


#### Predict trajectory from a candidate model, an initial state and GRN, and from a Patterning Scenario

If one wants to confront a given Patterning Scenario with a new initial state, and to generate trajectories with this selection of patterning functions, type the following:

`python solve.py predict toy [--q0 InitialCondition] [--GRN InitialGRN] [--plot 0 or 1]`

where phenotype (resp. GRN) associated with *InitialCondition* (resp. *InitialGRN*) is defined in the observations file. Argument for **--plot** is 1 (if concentration plots should be produced), otherwise 0.

**Clémence Réda, 2018**
