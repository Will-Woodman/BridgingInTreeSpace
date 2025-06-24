#!/bin/bash
##filenames
data_filename="/home/c1032934/Documents/Netbeans/TopInf20240503/ExperimentalData/EightYeast/yeast_new_ints.txt" # x0 tree filename
frechet_params_filename="/home/c1032934/Documents/Netbeans/TopInf20240503/ExperimentalData/EightYeast/yeast_new_ints_FM_0624test.txt"
output_filename="/home/c1032934/Documents/Netbeans/TopInf20240503/ExperimentalData/EightYeast/yeast_new_ints_MCMC0_FM_justdisp_624test.txt"

#Run the MCMC:
args=(
	$data_filename #data_filename
	$frechet_params_filename #INITIAL TREE: write "random" for random tree from the modal topology or "closest" for closest data point to Frechet mean or a name of a file containing a tree
	"0" # Initial squ root t_0 - set to zero if you want to use the Frechet variance
	"50" # Num steps
	"1051" # Seed
	"-n" #
	"5000" #"500000" # Num interations
	"-t"
	"1" # "100" # thin
	"-b"
	"100" #"500000" #burnin
	"-o"
	$output_filename #output file
	"-pbg"
	"0.08" #parameter for partial bridge proposal
	"-ptr"
	"0.1" # parameter for t0 random walk poposal
)

java -cp "/home/c1032934/Documents/Netbeans/BridgingInTreeSpace/dist/BridgingInTreeSpace.jar" bridge/InferBrownianParamsMCMC "${args[@]}"


