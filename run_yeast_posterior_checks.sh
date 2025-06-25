#!/bin/bash
source_tree_filename="/home/c1032934/Documents/Netbeans/TopInf20240503/ExperimentalData/EightYeast/yeast_new_ints_MCMCOutv2_splitModes_Tree.txt" # source tree filename
output_filename="/home/c1032934/Documents/Netbeans/TopInf20240503/ExperimentalData/EightYeast/yeast_new_ints_forward_sim_trees_test.txt" # data filename
distance_filename="/home/c1032934/Documents/Netbeans/TopInf20240503/ExperimentalData/EightYeast/yeast_new_ints_forward_sim_dists_test.txt" # data filename
seed="1324" #seed for running the random simulations
m=50 #number of steps
disp=0.0169 #dispersion to use for the random walks
numSamples=5000 #number of particules to forward simulate

args=(
        $source_tree_filename # data file name
        $output_filename # file name for outputting the sampled trees
        $seed
        $disp
        $m
        $numSamples
        )
 
java -cp "./dist/BridgingInTreeSpace.jar" simulateTops/GGFRWYeast "${args[@]}"

#now get the distances between the source tree and the simulated particles
args=(
        $output_filename # reread in the simulated trees
        $source_tree_filename # data file name
        $distance_filename  #file for outputting the distances
        )

java -cp "./dist/BridgingInTreeSpace.jar" simulateTops/RobinsonFouldsDataSet "${args[@]}"
