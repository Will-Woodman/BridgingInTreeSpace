#!/bin/bash
splits_filename="./YeastData/yeast_new_ints_MCMCOutv2_splitModes.txt"
modal_tree_filename="./YeastData/yeast_new_ints_MCMCOutv2_splitModes_Tree_test0624.txt"
frechet_mean_filename="./YeastData/yeast_new_ints_FM_OP.txt"

#build and save the modal tree:
args=(
        $splits_filename # modal splits
        $modal_tree_filename # modal tree output file
        )
        
java -cp "./dist/BridgingInTreeSpace.jar" simulateTops/MakeTreeFromSplits "${args[@]}"

#now get the sum of the internal edges on the two trees for comparison
args=(
        $modal_tree_filename # x0 tree filename
        )
        
java -cp "./dist/BridgingInTreeSpace.jar" simulateTops/getSumInternalEdges "${args[@]}"

args=(
        $frechet_mean_filename # x0 tree filename
        )
        
java -cp "./dist/BridgingInTreeSpace.jar" simulateTops/getSumInternalEdges "${args[@]}"
