#!/bin/bash
 
 data_filename="./YeastData/yeast_new_ints.txt" # x0 tree filename
 source_tree_filename="./YeastData/yeast_new_ints_MCMCOutput_ModeTree.txt" # data filename
 posterior_filename="./YeastData/yeast_new_ints_MCMCOutput_ModeTree_StepStone" # output file for the posterior


args=(
        $data_filename
        $source_tree_filename
        "0.13" # Squ root t_0
        "50" # Num steps
        "1802" # Seed
        "-n" # 
        "10000" # Num interations - before thin
        "-t" # 
        "20" # thin
        "-b" # 
        "10000" # burn-in
        "-o" # 
        $posterior_filename #posterior output files 
        "-pbg" # 
        "0.05" # geometric length bridge prop
        "200" # "50000" # number of Proposals to run
        "0.01" # the first non zero value of beta_k
        ) 
        

java -cp "./dist/BridgingInTreeSpace.jar" MarginalLikelihoods/SteppingStoneSampler "${args[@]}"

args=(
 	$posterior_filename
 	)
        
java -cp "./dist/BridgingInTreeSpace.jar" MarginalLikelihoodCalculations/StepStoneEstimate "${args[@]}" >> ./YeastData/StepStoneEstTest.txt #replace with file name for storing the Chib estimate
