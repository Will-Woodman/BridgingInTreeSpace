args = (
        "/home/c1032934/Documents/Netbeans/TopInf20240503/ExperimentalData/EightYeast/yeast_new_ints_MCMCOutv2_splitModes_Tree.txt" # x0 tree filename
        )
        
java -cp "/home/c1032934/Documents/Netbeans/BridgingInTreeSpace/dist/BridgingInTreeSpace.jar" simulateTops/getSumInternalEdges "${args[@]}"

args = (
        "/home/c1032934/Documents/Netbeans/TopInf20240503/ExperimentalData/EightYeast/yeast_new_ints_FM_OP.txt" # x0 tree filename
        )
        
java -cp "/home/c1032934/Documents/Netbeans/BridgingInTreeSpace/dist/BridgingInTreeSpace.jar" simulateTops/getSumInternalEdges "${args[@]}"
