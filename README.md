# Bridging in tree space

## Table of Contents
- [Installation](#installation)
- [Running the inference procedures](#running-the-inference-procedures)
  - [Posterior inference using bridges](#posterior-inference-using-bridges)
  - [Noisy MCMC for topology inference](#noisy-mcmc-for-topology-inference)
  - [Marginal likelihood for fixed dispersion - Chib and Tunnel](#marginal-likelihood-for-fixed-dispersion---chib-and-tunnel)
  - [Marginal likelihood for fixed dispersion - Stepping Stone](#marginal-likelihood-for-fixed-dispersion---stepping-stone)
  - [Marginal likelihood for unknown dispersion - Chib one block and Tunnel](#marginal-likelihood-for-unknown-dispersion---chib-one-block-and-tunnel)
  - [Marginal likelihood for unknown dispersion - Stepping Stone](#marginal-likelihood-for-unknown-dispersion---stepping-stone)
  - [Marginal likelihood for unknown dispersion - Chib Two Block](#marginal-likelihood-for-unknown-dispersion---chib-two-block)
## Installation


This Github repository includes the code required to perform the statistical inference procedures specified in the paper [Brownian motion, bridges and Bayesian inference in phylogenetic tree space](https://arxiv.org/abs/2506.22135) and my upcoming thesis.

An html file and corresponding RMarkdown file for reproducing the analysis on yeast gene trees detailed in the paper called ExperimentalDataRScript are provided as part of this project. <a href="https://will-woodman.github.io/BridgingInTreeSpace/ExperimentalDataRScript">ExperimentalDataRScript</a>
An RMarkdown file of the same name is given in the source code to reproduce the results.</li>.

The code can be downloaded from this Github repository and run from the command line. The easiest way to do this is with .sh files, where we can specify the various parameters required for the MCMC algorithms. The rest of this readme gives example shell scripts and parameters for running the different MCMC algorithms.

## Running the inference procedures

### Posterior inference using bridges
<details>
<summary>⭐ show/hide</summary>
Given a set of unrooted phylogenetic trees in Newick string format in a file called Example_trees.txt, withina folder called Example_folder, posterior inference can be run using a .sh script in the following form, where the parameters followed by a comment can be modifed as required:

<pre> ```
#!/bin/bash
##filenames
file_prefix="Example_trees"
folder_name="./Example_folder/"
data_filename="${folder_name}${file_prefix}.txt" # data filename
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

</details>

### Marginal likelihood for fixed dispersion - Chib and tunnel

<details>
<summary>⭐ show/hide</summary>
Given a set of unrooted phylogenetic trees in Newick string format in a file called Example_trees.txt, within a folder called Example_folder, to estimate the marginal likelihood using the Chib and tunnel (bridge) methods for a fixed value of dispersion and a fixed source tree given in the file Example_source.txt, a .sh script of the form specified below can be run. Parameters followed by a comment can be modifed as required.

<pre> ```
#!/bin/bash
 
 file_prefix="Example_trees"
 folder_name="./Example_folder/"
 source_tree_filename = "Example_source.txt"
 data_filename="${folder_name}${file_prefix}.txt" # data filename
 source_tree_filename="${folder_name}${source_tree_filename}" # data filename
 posterior_filename="${folder_name}${file_prefix}_Chib_Post.txt" # output filename for the samples from the posterior
 props_filename="${folder_name}${file_prefix}_Chib_Props.txt" # output filename for the samples from the proposals

#run the MCMC and independence proposals
args=(
        $data_filename
        $source_tree_filename
        "0.13" # squ root t_0
        "50" # Num steps
        "1260" # Seed
        "-n" #
        "1000000" # Num MCMC interations - before thin
        "-t" #
        "100" # thin
        "-b" #
        "100000" # burn-in
        "-o" # 
        $posterior_filename # output file for the posterior
        "-pbg" # 
        "0.05" # partial bridge proposal parameter
        "-numProps" # 
        "50000" # Num independence proposals to run
        $props_filename # output file for the proposals
        )
        
java -cp "./dist/BridgingInTreeSpace.jar" MarginalLikelihoods/ChibSampler "${args[@]}"

#calculate the estimates
args=(
	$props_filename
	$posterior_filename
)

java -cp "./dist/BridgingInTreeSpace.jar" MarginalLikelihoodCalculations/ChibJeliEstimate "${args[@]}"  >> ${folder_name}${file_prefix}_ChibEst.txt #replace with file name for storing the Chib estimate

java -cp "./dist/BridgingInTreeSpace.jar" MarginalLikelihoodCalculations/BridgeSamplingEstimate "${args[@]}"  >> ${folder_name}${file_prefix}TunnelEst.txt #replace with file name for storing the tunnel estimate
``` </pre>
</details>

### Marginal likelihood for fixed dispersion - Stepping Stone

<details>
<summary>⭐ show/hide</summary>
Given a set of unrooted phylogenetic trees in Newick string format in a file called Example_trees.txt, within a folder called Example_folder, to estimate the marginal likelihood using the stepping stone method for a fixed value of dispersion and a fixed source tree given in the file Example_source.txt, a .sh script of the form specified below can be run. Parameters followed by a comment can be modifed as required.

<pre> ```
#!/bin/bash
 
 file_prefix="Example_trees"
 folder_name="./Example_folder/"
 source_tree_filename = "Example_source.txt"
 data_filename="${folder_name}${file_prefix}.txt" # data filename
 source_tree_filename="${folder_name}${source_tree_filename}" # data filename
 posterior_filename="${folder_name}${file_prefix}_Chib_Post.txt" # output filename for the samples from the posterior


args=(
        $data_filename
        $source_tree_filename
        "0.13" # Square root of dispersion
        "50" # Num steps in the bridges
        "1802" # Set the seed for the random number generator
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
        
java -cp "./dist/BridgingInTreeSpace.jar" MarginalLikelihoodCalculations/StepStoneEstimate "${args[@]}" >> ${folder_name}${file_prefix}StepStoneEst.txt #file name for storing the Chib estimate

``` </pre>
</details>

### Marginal likelihood for unknown dispersion - Chib one block and tunnel

<details>
<summary>⭐ show/hide</summary>
Given a set of unrooted phylogenetic trees in Newick string format in a file called Example_trees.txt, within a folder called Example_folder, to estimate the marginal likelihood using the Chib one block and tunnel (bridge) methods where the dispersion is considered an unknown nuisance parameter and a fixed source tree is given in the file Example_source.txt, a .sh script of the form specified below can be run. Parameters followed by a comment can be modifed as required.

<pre> ```
#!/bin/bash
 
 file_prefix="Example_trees"
 folder_name="./Example_folder/"
 source_tree_filename = "Example_source.txt"
 data_filename="${folder_name}${file_prefix}.txt" # data filename
 source_tree_filename="${folder_name}${source_tree_filename}" # data filename
 posterior_filename="${folder_name}${file_prefix}_Chib_Post.txt" # output filename for the samples from the posterior
 props_filename="${folder_name}${file_prefix}_Chib_Props.txt" # output filename for the samples from the proposals

#run the MCMC and independence proposals
args=(
        $data_filename
        $source_tree_filename
        "0.13" # squ root t_0
        "50" # Num steps
        "1260" # Seed
        "-n" #
        "1000000" # Num MCMC interations - before thin
        "-t" #
        "100" # thin
        "-b" #
        "100000" # burn-in
        "-o" # 
        $posterior_filename # output file for the posterior
        "-pbg" # 
        "0.05" # partial bridge proposal parameter
        "-numProps" # 
        "50000" # Num independence proposals to run
        $props_filename # output file for the proposals
        )
        
java -cp "./dist/BridgingInTreeSpace.jar" MarginalLikelihoodsWithDisp/TunnelSamplingDisp "${args[@]}"

#calculate the estimates

args=(
	$prop_simple_filename
	$posterior_filename
	$props_filename
	$additional_dispersion_file
)

java -cp "./dist/BridgingInTreeSpace.jar" MarginalLikelihoodCalculationsWithDisp/ChibJeliEstimateDisp "${args[@]}"  >> ${folder_name}${file_prefix}_ChibEst.txt #replace with file name for storing the Chib estimate

java -cp "./dist/BridgingInTreeSpace.jar" MarginalLikelihoodCalculationsWithDisp/TunnelSamplingEstimateDisp "${args[@]}"  >> ${folder_name}${file_prefix}TunnelEst.txt #replace with file name for storing the tunnel estimate
``` </pre>
</details>



### Marginal likelihood for unknown dispersion - Stepping Stone

<details>
<summary>⭐ show/hide</summary>
Given a set of unrooted phylogenetic trees in Newick string format in a file called Example_trees.txt, within a folder called Example_folder, to estimate the marginal likelihood using the stepping stone method where the dispesion is consider an unknown nuisance parameter and a fixed source tree is given in the file Example_source.txt, a .sh script of the form specified below can be run. Parameters followed by a comment can be modifed as required.

<pre> ```
#!/bin/bash
 
 file_prefix="Example_trees"
 folder_name="./Example_folder/"
 source_tree_filename = "Example_source.txt"
 data_filename="${folder_name}${file_prefix}.txt" # data filename
 source_tree_filename="${folder_name}${source_tree_filename}" # data filename
 posterior_filename="${folder_name}${file_prefix}_Chib_Post.txt" # output filename for the samples from the posterior


args=(
        $data_filename
        $source_tree_filename
        "0.13" # Square root of dispersion
        "50" # Num steps in the bridges
        "1802" # Set the seed for the random number generator
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
        

java -cp "./dist/BridgingInTreeSpace.jar" MarginalLikelihoodsWithDisp/SteppingStoneSamplingDisp "${args[@]}"

args=(
	$prop_simple_filename
 	$posterior_filename
 	)
        
java -cp "./dist/BridgingInTreeSpace.jar" MarginalLikelihoodCalculationsWithDisp/StepStoneEstimateDisp "${args[@]}" >> ${folder_name}${file_prefix}StepStoneEst.txt #file name for storing the Chib estimate

``` </pre>
</details>


### Marginal likelihood for unknown dispersion - Chib two block

<details>
<summary>⭐ show/hide</summary>
Given a set of unrooted phylogenetic trees in Newick string format in a file called Example_trees.txt, within a folder called Example_folder, to estimate the marginal likelihood using the Chib two block method where the dispersion is considered an unknown nuisance parameter and a fixed source tree given in the file Example_source.txt, a .sh script of the form specified below can be run. Parameters followed by a comment can be modifed as required.

<pre> ```
#!/bin/bash
 
 file_prefix="Example_trees"
 folder_name="./Example_folder/"
 source_tree_filename = "Example_source.txt"
 data_filename="${folder_name}${file_prefix}.txt" # data filename
 source_tree_filename="${folder_name}${source_tree_filename}" # data filename
 posterior_filename="${folder_name}${file_prefix}_Chib_Post.txt" # output filename for the samples from the posterior
 props_filename="${folder_name}${file_prefix}_Chib_Props.txt" # output filename for the samples from the proposals

#run the MCMC and independence proposals
args=(
        $data_filename
        $source_tree_filename
        "0.13" # squ root t_0
        "50" # Num steps
        "1260" # Seed
        "-n" #
        "1000000" # Num MCMC interations - before thin
        "-t" #
        "100" # thin
        "-b" #
        "100000" # burn-in
        "-o" # 
        $posterior_filename # output file for the posterior
        "-pbg" # 
        "0.05" # partial bridge proposal parameter
        "-numProps" # 
        "50000" # Num independence proposals to run
        $props_filename # output file for the proposals
        )
        
java -cp "./dist/BridgingInTreeSpace.jar" MarginalLikelihoodsWithDisp/ChibTwoBlock "${args[@]}"


#calculate the estimates
args=(
	$fixed_disp_posterior_filename
	$fixed_disp_props_filename
	$posterior_filename
)

java -cp "./dist/BridgingInTreeSpace.jar" MarginalLikelihoodCalculationsWithDisp/ChibTwoBlockEstimate "${args[@]}"  >> ${folder_name}${file_prefix}_ChibEst.txt #replace with file name for storing the Chib estimate

``` </pre>
</details>


### Noisy MCMC for topology inference

<details>
<summary>⭐ show/hide</summary>

Suppose we are given a set of unrooted phylogenetic trees in Newick string format in a file called Example_trees.txt, within a folder called Example_folder. To run the noisy MCMC algorithm for inferring the species tree when gene trees are modelled as draws from a Brownian motion kernel (approximated by a random walk kernel), run a .sh script of the following form (note that in this analysis we manually input an initial tree):

<pre> ```
#!/bin/bash
##filenames
file_prefix="Example_trees"
folder_name="./Example_folder/"
source_tree_filename = "Example_source.txt"
data_filename="${folder_name}${file_prefix}.txt" # data filename
output_filename="${folder_name}${file_prefix}_MCMCOutput.txt"
topologies_filename="${folder_name}${file_prefix}_MCMCOutput_tops.txt"

#Run the MCMC:
args=(
	$data_filename # data filename
	$initial_tree_filename # initial tree filename
	$output_filename # output file for the MCMC
	$topologies_filename # a file to count the topologies in the MCMC output into
	"1" # Initial squ root t_0 
	"50" # Num steps in the random walks
    "10000" # Number of particles in the approx likelihood calcs
	"1" # Number of cores to use in forward simulating particles
	"15000" # Num interations
	"0" #burnin
	"100" # thin
	"0.1" # parameter for the log random walk edge proposal
	"0.1" # parameter for the random walk edge proposal
)

java -cp "./dist/BridgingInTreeSpace.jar" topologies/InferParamsViaApproxLike "${args[@]}"

``` </pre>


</details>

