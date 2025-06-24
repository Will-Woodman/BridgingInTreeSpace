#!/bin/bash
 
 data_filename="/home/c1032934/Documents/Netbeans/TopInf20240503/ExperimentalData/EightYeast/yeast_new_ints.txt" # x0 tree filename
 source_tree_filename="/home/c1032934/Documents/Netbeans/TopInf20240503/ExperimentalData/EightYeast/yeast_new_ints_FM_OP.txt" # data filename
 posterior_filename="/data/ww24/ExperimentalData/EightYeast/yeast_new_ints_FM_OP_ChibJeliPostv2_test0624.txt" # output filename for the samples from the posterior

#run the MCMC and independence proposals
args=(
        $data_filename
        $source_tree_filename
        "0.13" # squ root t_0
        "50" # Num steps
        "1260" # Seed
        "-n" #
        "1000" # "1000000" # Num MCMC interations - before thin
        "-t" #
        "1" # "100" # thin
        "-b" #
        "100" # "100000" # burn-in
        "-o" # 
        $posterior_filename # output file for the posterior
        "-pbg" # 
        "0.05" # partial bridge proposal parameter
        "-numProps" # 
        "500" # "50000" # Num independence proposals to run
        $props_filename # output file for the proposals
        )
        
java -cp "/home/c1032934/Documents/Netbeans/BridgingInTreeSpace/dist/BridgingInTreeSpace.jar" MarginalLikelihoods/ChibSampler "${args[@]}"

#calculate the estimates
args=(
	$props_filename
	$posterior_filename
)

java -cp "/home/c1032934/Documents/Netbeans/BridgingInTreeSpace/dist/BridgingInTreeSpace.jar" MarginalLikelihoodCalculations/ChibJeliEstimate "${args[@]}"  >> /home/c1032934/Documents/Netbeans/TopInf20240503/ExperimentalData/EightYeast/ChibJeliEstTest.txt #replace with file name for storing the Chib estimate

java -cp "/home/c1032934/Documents/Netbeans/BridgingInTreeSpace/dist/BridgingInTreeSpace.jar" MarginalLikelihoodCalculations/BridgeSamplingEstimate "${args[@]}"  >> /home/c1032934/Documents/Netbeans/TopInf20240503/ExperimentalData/EightYeast/TunnelEstTest.txt #replace with file name for storing the tunnel estimate





