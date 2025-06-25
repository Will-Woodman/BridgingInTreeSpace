#!/bin/bash
##filenames
data_filename="/home/c1032934/Documents/Netbeans/TopInf20240503/ExperimentalData/EightYeast/yeast_new_ints.txt" # x0 tree filename
frechet_params_filename="/home/c1032934/Documents/Netbeans/TopInf20240503/ExperimentalData/EightYeast/yeast_new_ints_FM_0624test.txt"
frechet_distances_filename="/home/c1032934/Documents/Netbeans/TopInf20240503/ExperimentalData/EightYeast/yeast_new_ints_FM_distances_0624test.txt"
output_filename="/home/c1032934/Documents/Netbeans/TopInf20240503/ExperimentalData/EightYeast/yeast_new_ints_MCMC0624test.txt"
topologies_filename="/home/c1032934/Documents/Netbeans/TopInf20240503/ExperimentalData/EightYeast/yeast_new_ints_MCMC0624test_tops.txt"
modal_top_edges_filename="/home/c1032934/Documents/Netbeans/TopInf20240503/ExperimentalData/EightYeast/yeast_new_ints_MCMC0624test_edges.txt"

#Initially estimate the Frechet mean:
args=(
	$data_filename
	$frechet_params_filename
	$frechet_distances_filename
)

java -cp "./dist/BridgingInTreeSpace.jar" simulateTops/FrechetMeanMainClass "${args[@]}"

#Run the MCMC:
args=(
	$data_filename #data_filename
	"closest" #INITIAL TREE: write "random" for random tree from the modal topology or "closest" for closest data point to Frechet mean
	"0" # Initial squ root t_0 - set to zero if you want to use the Frechet variance
	"50" # Num steps
	"1051" # Seed
	"-n" #
	"500" #"5000000" # Num interations
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
	"-pxg"
	"0.002" #parameters for x0 proposal resampling beginning of bridges
	"0.9"    
)

java -cp "./dist/BridgingInTreeSpace.jar" bridge/InferBrownianParamsMCMC "${args[@]}"

##Count the topologies in the posterior
args=(
	$output_filename #MCMCFilename
	$topologies_filename #TopologiesFilename
        )
        
java -cp "./dist/BridgingInTreeSpace.jar" topologies/countTopologiesAtPoints "${args[@]}"

#Get the edge lengths in the modal topology in the posterior:
args=(
        $output_filename #MCMCFilename
	$modal_top_edges_filename #EdgesFilename
        )
        
java -cp "./dist/BridgingInTreeSpace.jar" topologies/edgeLengthsModalTop "${args[@]}"




