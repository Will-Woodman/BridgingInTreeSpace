#!/bin/bash
splits_filename="./YeastData/yeast_new_ints_MCMCOutput_splitModes.txt"
modal_tree_filename="./YeastData/yeast_new_ints_MCMCOutput_ModeTree.txt"
frechet_mean_filename="./YeastData/yeast_new_ints_FM.txt"

#build and save the modal tree:
args=(
        $splits_filename # modal splits
        $modal_tree_filename # modal tree output file
        )
        
java -cp "./dist/BridgingInTreeSpace.jar" topologies.BuildModeTree "${args[@]}"

#now get the sum of the internal edges on the two trees for comparison
args=(
        $modal_tree_filename # x0 tree filename
        )

echo 'The sum of internal edge lengths in the modal tree is'

java -cp "./dist/BridgingInTreeSpace.jar" BMPaperFiles/GetSumInternalEdges "${args[@]}"

args=(
        $frechet_mean_filename # x0 tree filename
        )
        
echo 'The sum of internal edge lengths in the Frechet mean is'

java -cp "./dist/BridgingInTreeSpace.jar" BMPaperFiles/GetSumInternalEdges "${args[@]}"
