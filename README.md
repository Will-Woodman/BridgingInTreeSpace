This Github repository includes the code required to perform the statistical inference procedures specified in the paper [Brownian motion, bridges and Bayesian inference in phylogenetic tree space](https://arxiv.org/abs/2506.22135) and my upcoming thesis.

An html file and corresponding RMarkdown file for reproducing the analysis on yeast gene trees detailed in the paper called ExperimentalDataRScript are provided as part of this project.

The code can be downloaded from this Github repository and run from the command line.

Given a set of unrooted phylogenetic trees in Newick string format in a file called Example_trees.txt, withina folder called Example_folder, posterior inference can be run using a .sh script in the following form, where the parameters followed by a comment can be modifed as required:

<pre> ```
#!/bin/bash
##filenames
file_prefix="Example_trees"
folder_name="./Example_folder/"
data_filename="${folder_name}${file_prefix}.txt" # x0 tree filename
frechet_params_filename="${folder_name}${file_prefix}_FM_params.txt"
frechet_distances_filename="${folder_name}${file_prefix}_FM_distances.txt"
output_filename="${folder_name}${file_prefix}_MCMCOutput.txt"
topologies_filename="${folder_name}${file_prefix}_MCMCOutput_tops.txt"
modal_top_edges_filename="${folder_name}${file_prefix}_MCMCOutput_edges.txt"

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
	"5000000" # Num interations
	"-t"
	"100" # thin
	"-b"
	"500000" #burnin
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

``` </pre>

To run the noisy MCMC algorithm for inferring the species tree when gene trees are modelled as draws from a Brownian motion kernel (approxiimated by a random walk kernel, run a .sh script of the following form:



