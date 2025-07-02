#!/bin/bash
##filenames
data_filename="./YeastData/yeast_new_ints.txt" # x0 tree filename
frechet_params_filename="./YeastData/yeast_new_ints_FM.txt"
output_filename="./YeastData/yeast_new_ints_MCMC0utput_FM_justdisp.txt"

#Run the MCMC:
args=(
	$data_filename #data_filename
	$frechet_params_filename #INITIAL TREE: write "random" for random tree from the modal topology or "closest" for closest data point to Frechet mean or a name of a file containing a tree
	"0" # Initial squ root t_0 - set to zero if you want to use the Frechet variance
	"50" # Num steps
	"1051" # Seed
	"-n" #
	"500000" # Num interations
	"-t"
	"100" # thin
	"-b"
	"100000" #burnin
	"-o"
	$output_filename #output file
	"-pbg"
	"0.08" #parameter for partial bridge proposal
	"-ptr"
	"0.1" # parameter for t0 random walk poposal
)

java -cp "./dist/BridgingInTreeSpace.jar" bridge/InferBrownianParamsMCMC "${args[@]}"


