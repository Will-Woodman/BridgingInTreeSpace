#!/bin/bash
 
 data_filename="/home/c1032934/Documents/Netbeans/TopInf20240503/ExperimentalData/EightYeast/yeast_new_ints.txt" # x0 tree filename
 source_tree_filename="/home/c1032934/Documents/Netbeans/TopInf20240503/ExperimentalData/EightYeast/yeast_new_ints_FM_OP.txt" # data filename
 posterior_filename="/data/ww24/ExperimentalData/EightYeast/yeast_new_ints_FM_OP_StepStone_test0624" # output file for the posterior


args=(
        $data_filename
        $source_tree_filename
        "0.13" # Squ root t_0
        "50" # Num steps
        "1802" # Seed
        "-n" # 
        "50" # "10000" # Num interations - before thin
        "-t" # 
        "1" # "20" # thin
        "-b" # 
        "100" # "10000" # burn-in
        "-o" # 
        $posterior_filename #posterior output files 
        "-pbg" # 
        "0.05" # geometric length bridge prop
        "200" # "50000" # number of Proposals to run
        "0.01" # the first non zero value of beta_k
        ) 
        

java -cp "/home/c1032934/Documents/Netbeans/BridgingInTreeSpace/dist/BridgingInTreeSpace.jar" MarginalLikelihoods/SteppingStoneSampler "${args[@]}"

args=(
 	$posterior_filename
 	)
        
java -cp "/home/c1032934/Documents/Netbeans/BridgingInTreeSpace/dist/BridgingInTreeSpace.jar" MarginalLikelihoodCalculations/StepStoneEstimate "${args[@]}" >> /home/c1032934/Documents/Netbeans/TopInf20240503/ExperimentalData/EightYeast/StepStoneEstTest.txt #replace with file name for storing the Chib estimate
