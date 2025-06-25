#!/bin/bash

data_filename="./YeastData/yeast_new_ints.txt" # x0 tree filename
output_filename="./YeastData/yeast_new_ints_star_tree_ML_test.txt" #where to put the output
 
args=(
        $data_filename # data file name
        )
 
java -cp "./dist/BridgingInTreeSpace.jar" MarginalLikelihoods/ExactMLStarTree "${args[@]}" >> $output_filename
