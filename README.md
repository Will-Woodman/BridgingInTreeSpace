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
 ref_dist_filename="${folder_name}${file_prefix}_Chib_Props.txt" # output filename for the samples from the reference distribution 
 ref_dist_parameters_filename = "${folder_name}${file_prefix}_Chib_Ref_Dist_Params.txt"
 additional_ref_dist_filename="${folder_name}${file_prefix}_Chib_Additional_Disp.txt" # output filename for the additional file for evaluating the density of the reference distribution for dispersion
 prop_simple_filename="${folder_name}${file_prefix}_Chib_Prop_Simple.txt" # output filename for estimating the normalising constant of the reference distribution

String[] args1 = new String[21];
        //arguments for the posterior MCMC sims -- may want to use different parameters for the posterior and reference dists
        args1[0] = args[0]; //data filename
        args1[1] = args[1]; //source tree filename
        args1[2] = "0.0"; // // Initial squ root t_0 - delete this option eventually
        args1[3] = "20"; // Num steps
        args1[4] = "344"; // Seed
        args1[5] = "-n";
        args1[6] = "250000"; // Num interations - before thin
        args1[7] = "-t";
        args1[8] = "20"; // thin
        args1[9] = "-b";
        args1[10] = "10000"; //burn-in
        args1[11] = "-o";
        args1[12] = "/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/MarginalLikelihoods/MultiplePointTest/Post.txt"; //file_name to store posterior sims
        args1[13] = "-pbg";
        args1[14] = "0.01"; //geometric length partial bridge proposal parameter
        args1[15] = "-ptr";
        args1[16] = "0.5"; // dispersion log random walk proposal parameter
        args1[17] = "-prd"; //whether to parametrise the lognormal reference distribution for t0: 1 for yes 0 for no -- advisable to do so
        args1[18] = "1";
        args1[19] =""; //file to store the parameters of the t0 reference distribution in
        args1[20] =""; //file to store the new t0 densities (for the reference dist) in 

        
       
        String[] args2 = new String[19];
        //arguments for the ref dist MCMC sims -- may want to use different parameters for the posterior and reference dists
        args2[0] = args1[0];
        args2[1] = args1[1];
        args2[2] = "0.0"; //Initial squ root t_0 - delete this option eventually
        args2[3] = args1[3]; // Num steps
        args2[4] = "0.0"; //Placeholder for reference dist params
        args2[5] = "0.0"; //Placeholder for reference dist params
        args2[6] = "897"; // Seed
        args2[7] = "-n";
        args2[8] = "250000"; // Num interations - before thin
        args2[9] = "-t";
        args2[10] = "20"; // thin
        args2[11] = "-b";
        args2[12] = "10000";//burnin
        args2[13] = "-o";
        args2[14] = "";//output filename for simple indep props via MCMC
        args2[15] = "-pbg";
        args2[16] = "0.01"; //geometric length partial bridge proposal parameter
        args2[17] = "-ptr"; //dispersion log random walk proposal parameter
        args2[18] = "0.5";
        
        String[] args3 = new String[8];
        //arguments for calculating the proportion of simple proposals
        args3[0] = args1[0];
        args3[1] = args1[1];
        args3[2] = ""; // OutFileName
        args3[3] = "0.0"; //Placeholder for reference dist params
        args3[4] = "0.0"; //Placeholder for reference dist params
        args3[5] = "200"; //numDisps
        args3[6] = "500"; //numProps
        args3[7] = args1[3]; //num steps
	
#run the MCMC and independence proposals
args=(
        $data_filename
        $source_tree_filename
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
		"-ptr"
		"0.5" # dispersion log random walk proposal parameter
        "-prd" 
        "1" # whether to parametrise the lognormal reference distribution for t0: 1 for yes 0 for no -- advisable to do so
        $ref_dist_parameters_filename # file to store the parameters of the t0 reference distribution in
        $additional_ref_dist_filename # file to store the new t0 densities (for the reference dist) in 
		"897" # Seed - ref dist MCMC sims
        "-n"
        "250000" # Num interations - before thin - ref dist MCMC sims
        "-t"
        "20" # thin - ref dist MCMC sims
        "-b"
        "10000" # burnin - ref dist MCMC sims
        "-o"
        $ref_dist_filename # output filename for simple indep props via MCMC
        "-pbg"
        "0.01" # geometric length partial bridge proposal parameter
        "-ptr" 
        "0.5" # dispersion log random walk proposal parameter
	    $prop_simple_filename # Output file name for estimating normalizing constant of the reference distribution
        "200" # number of values of dispersion to use in the numerical integration
        "500" # number of proposals per data point per value of dispersion
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

