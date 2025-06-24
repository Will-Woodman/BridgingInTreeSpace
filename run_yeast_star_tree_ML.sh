#!/bin/bash

data_filename="/home/c1032934/Documents/Netbeans/TopInf20240503/ExperimentalData/EightYeast/yeast_new_ints.txt" # x0 tree filename
output_filename="/home/c1032934/Documents/Netbeans/TopInf20240503/ExperimentalData/EightYeast/yeast_new_ints_star_tree_ML_test.txt" #where to put the output
 
args=(
        $data_filename # data file name
        )
 
java -cp "/home/c1032934/Documents/Netbeans/BridgingInTreeSpace/dist/BridgingInTreeSpace.jar" MarginalLikelihoods/ExactMLStarTree "${args[@]}" >> $output_filename
