# Bridging in tree space

## Table of Contents
- [Installation](#installation)
- [Running the inference procedures](#running-the-inference-procedures)
  - [Simulating a data set under the model](#simulating-a-data-set-under-the-model)
  - [Posterior inference using bridges](#posterior-inference-using-bridges)
  - [Marginal likelihood for fixed dispersion - Chib and Tunnel](#marginal-likelihood-for-fixed-dispersion---chib-and-tunnel)
  - [Marginal likelihood for fixed dispersion - Stepping Stone](#marginal-likelihood-for-fixed-dispersion---stepping-stone)
  - [Marginal likelihood for unknown dispersion - Chib one block and Tunnel](#marginal-likelihood-for-unknown-dispersion---chib-one-block-and-tunnel)
  - [Marginal likelihood for unknown dispersion - Stepping Stone](#marginal-likelihood-for-unknown-dispersion---stepping-stone)
  - [Marginal likelihood for unknown dispersion - Chib Two Block](#marginal-likelihood-for-unknown-dispersion---chib-two-block)
  - [Simulating from the Gaussian kernel distribution using MCMC](#simulating-from-the-gaussian-kernel-distribution-using-mcmc)
  - [Noisy MCMC for topology inference](#noisy-mcmc-for-topology-inference)
    
## Installation
Simulating from the Gaussian kernel distribution using MCMC

This Github repository includes the code required to perform the statistical inference procedures specified in the paper [Brownian motion, bridges and Bayesian inference in phylogenetic tree space](https://arxiv.org/abs/2506.22135) and my upcoming thesis.

An html file and corresponding RMarkdown file for reproducing the analysis on yeast gene trees detailed in the paper called ExperimentalDataRScript are provided as part of this project. <a href="https://will-woodman.github.io/BridgingInTreeSpace/ExperimentalDataRScript">ExperimentalDataRScript</a>
An RMarkdown file of the same name is given in the source code to reproduce the results.</li>.

The code can be downloaded from this Github repository and run from the command line. The easiest way to do this is with .sh files, where we can specify the various parameters required for the MCMC algorithms. The rest of this readme gives example shell scripts and parameters for running the different MCMC algorithms.

We alse give minimal examples of R scripts used to produce the analysis of MCMC runs, although it is likely that any analysis of MCMC runs will be bespoke.

## Running the inference procedures

### Simulating a data set under the model
<details>
<summary>⭐ show/hide</summary>
This file will simulate a random source tree called Example_source in a folder called Example_folder with a specifed number N of taxa from a coalescent model (with edges resampled from a Gamma distributionwith shape-scale parametrisation). It will then produce a data set of n trees by simulating random walks with m steps, started at the source tree and a specified value of dispersion. It also outputs the topologies of the source tree and data as well as geodesic and Robinson Foulds distances from the source tree to the data points.

```bash
#!/bin/bash
##filenames
file_prefix="Example_trees"
source_tree_prefix="Example_source" # data filename
folder_name="./Example_folder/"
source_tree_filename="${folder_name}${source_tree_prefix}.txt" # source tree filename
data_filename="${folder_name}${file_prefix}.txt" # data filename
topologies_filename="${folder_name}${file_prefix}_tops.txt" #file to store the topologies of the simulated data
distances_filename="${folder_name}${file_prefix}_distances_to_source.txt"
source_tree_topology_filename="${folder_name}${file_prefix}_tops.txt" #file to store the topology of the source tree

#Simulate the trees:
args=(
	$source_tree_filename # output file for simulated source trees
	$data_filename # output file for the simulated trees
	$topologies_filename # output file for the topologies of the simulated trees
	$distances_filename # output file for the geodesic distances of the simulated trees to the source        
	$source_tree_topology_filename #file to store the topology of the source tree
	"5" # number of taxa
	"5" # number  of taxa in the trees
	"947" # seed for the random engine
	"100" # number of steps in the random walk
	"2.0" # shape of the gamma distribution for simulating source tree edges
	"0.05" # scale of the gamma distribution for simulating source tree edges
	"0.25" # value of t0 in the random walks
)

java -cp "./dist/BridgingInTreeSpace.jar" simulateTops/simulateSourceTree "${args[@]}"
```

</details>


### Posterior inference using bridges
<details>
<summary>⭐ show/hide</summary>
Given a set of unrooted phylogenetic trees in Newick string format in a file called Example_trees.txt, withina folder called Example_folder, posterior inference can be run using a .sh script in the following form, where the parameters followed by a comment can be modifed as required:

```bash
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
	"5000" # Num interations
	"-t"
	"10" # thin
	"-b"
	"5000" #burnin
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


```

If you have built the mode tree using R you will have outputted a list of splits with their lengths. We can turn this into a phylogenetic tree in Newick string form using the following shell script:

```bash
#!/bin/bash
file_prefix="Example_trees"
folder_name="./Example_folder/"
splits_filename="${folder_name}${file_prefix}_MCMCOutput_mode_tree_edges.txt" # file containing list of splits and lengths
modal_tree_filename="${folder_name}${file_prefix}_MCMCOutput_mode_tree.txt"   # output file to store the mode tree

# build and save the modal tree:
args=(
    $splits_filename     # modal splits
    $modal_tree_filename  # modal tree output file
)

java -cp "./dist/BridgingInTreeSpace.jar" topologies.BuildModeTree "${args[@]}"
```

</details>

### Marginal likelihood for fixed dispersion - Chib and tunnel

<details>
<summary>⭐ show/hide</summary>
Given a set of unrooted phylogenetic trees in Newick string format in a file called Example_trees.txt, within a folder called Example_folder, to estimate the marginal likelihood using the Chib and tunnel (bridge) methods for a fixed value of dispersion and a fixed source tree given in the file Example_source.txt, a .sh script of the form specified below can be run. Parameters followed by a comment can be modifed as required.

```bash
file_prefix="Example_trees"
folder_name="./Example_folder/"
source_tree_filename="Example_source.txt"
data_filename="${folder_name}${file_prefix}.txt" # data filename
source_tree_filename="${folder_name}${source_tree_filename}" # data filename
posterior_filename="${folder_name}${file_prefix}_Chib_Fixed_Disp_Post.txt" # output filename for the samples from the posterior
props_filename="${folder_name}${file_prefix}_Chib_Fixed_Disp_Props.txt" # output filename for the samples from the proposals

#run the MCMC and independence proposals
args=(
        $data_filename
        $source_tree_filename
        "0.13" # squ root t_0
        "50" # Num steps
        "1260" # Seed
        "-n" #
        "10000" # Num MCMC interations - before thin
        "-t" #
        "100" # thin
        "-b" #
        "1000" # burn-in
        "-o" # 
        $posterior_filename # output file for the posterior
        "-pbg" # 
        "0.05" # partial bridge proposal parameter
        "-numProps" # 
        "500" # Num independence proposals to run
        $props_filename # output file for the proposals
        )
        
java -cp "./dist/BridgingInTreeSpace.jar" MarginalLikelihoods/ChibSampler "${args[@]}"

#calculate the estimates
args=(
	$props_filename
	$posterior_filename
)

java -cp "./dist/BridgingInTreeSpace.jar" MarginalLikelihoodCalculations/ChibJeliEstimate "${args[@]}"  >> ${folder_name}${file_prefix}_FixedChibEst.txt #replace with file name for storing the Chib estimate

java -cp "./dist/BridgingInTreeSpace.jar" MarginalLikelihoodCalculations/BridgeSamplingEstimate "${args[@]}"  >> ${folder_name}${file_prefix}_FixedTunnelEst.txt #replace with file name for storing the tunnel estimate

```
</details>

### Marginal likelihood for fixed dispersion - Stepping Stone

<details>
<summary>⭐ show/hide</summary>
Given a set of unrooted phylogenetic trees in Newick string format in a file called Example_trees.txt, within a folder called Example_folder, to estimate the marginal likelihood using the stepping stone method for a fixed value of dispersion and a fixed source tree given in the file Example_source.txt, a .sh script of the form specified below can be run. Parameters followed by a comment can be modifed as required.

```bash
#!/bin/bash
 
 file_prefix="Example_trees"
 folder_name="./Example_folder/"
 source_tree_filename="Example_source.txt"
 data_filename="${folder_name}${file_prefix}.txt" # data filename
 source_tree_filename="${folder_name}${source_tree_filename}" # data filename
 posterior_filename="${folder_name}${file_prefix}_Step_Stone" # output filename for the samples from the posterior

args=(
        $data_filename
        $source_tree_filename
        "0.13" # Square root of dispersion
        "50" # Num steps in the bridges
        "1802" # Set the seed for the random number generator
        "-n" # 
        "1000" # Num interations - before thin
        "-t" # 
        "20" # thin
        "-b" # 
        "1000" # burn-in
        "-o" # 
        $posterior_filename #posterior output files 
        "-pbg" # 
        "0.05" # geometric length bridge prop
        "200" # number of Proposals to run
        "0.01" # the first non zero value of beta_k
        ) 
        

java -cp "./dist/BridgingInTreeSpace.jar" MarginalLikelihoods/SteppingStoneSampler "${args[@]}"

args=(
 	$posterior_filename
 	)
        
java -cp "./dist/BridgingInTreeSpace.jar" MarginalLikelihoodCalculations/StepStoneEstimate "${args[@]}" >> ${folder_name}${file_prefix}_Fixed_Disp_StepStoneEst.txt #file name for storing the stepping stone estimate

```
</details>

### Marginal likelihood for unknown dispersion - Chib one block and tunnel

<details>
<summary>⭐ show/hide</summary>
Given a set of unrooted phylogenetic trees in Newick string format in a file called Example_trees.txt, within a folder called Example_folder, to estimate the marginal likelihood using the Chib one block and tunnel (bridge) methods where the dispersion is considered an unknown nuisance parameter and a fixed source tree is given in the file Example_source.txt, a .sh script of the form specified below can be run. Parameters followed by a comment can be modifed as required.

```bash
#!/bin/bash
 
 file_prefix="Example_trees"
 folder_name="./Example_folder/"
 source_tree_filename="Example_source.txt"
 data_filename="${folder_name}${file_prefix}.txt" # data filename
 source_tree_filename="${folder_name}${source_tree_filename}" # x0 filename
 posterior_filename="${folder_name}${file_prefix}_Chib_Post.txt" # output filename for the samples from the posterior
 ref_dist_filename="${folder_name}${file_prefix}_Chib_Props.txt" # output filename for the samples from the reference distribution 
 ref_dist_parameters_filename="${folder_name}${file_prefix}_Chib_Ref_Dist_Params.txt"
 additional_ref_dist_filename="${folder_name}${file_prefix}_Chib_Additional_Disp.txt" # output filename for the additional file for evaluating the density of the reference distribution for dispersion
 prop_simple_filename="${folder_name}${file_prefix}_Chib_Prop_Simple.txt" # output filename for estimating the normalising constant of the reference distribution
echo $source_tree_filename	

#run the MCMC and independence proposals
args=(
        $data_filename
        $source_tree_filename
        "50" # Num steps
        "1260" # Seed
        "-n" #
        "10000" # Num MCMC interations - before thin
        "-t" #
        "10" # thin
        "-b" #
        "1000" # burn-in
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
        "10000" # Num interations - before thin - ref dist MCMC sims
        "-t"
        "20" # thin - ref dist MCMC sims
        "-b"
        "1000" # burnin - ref dist MCMC sims
        "-o"
        $ref_dist_filename # output filename for simple indep props via MCMC
        "-pbg"
        "0.01" # geometric length partial bridge proposal parameter
        "-ptr" 
        "0.5" # dispersion log random walk proposal parameter
	$prop_simple_filename # Output file name for estimating normalizing constant of the reference distribution
        "20" # number of values of dispersion to use in the numerical integration
        "50" # number of proposals per data point per value of dispersion
        )
        
java -cp "./dist/BridgingInTreeSpace.jar" MarginalLikelihoodsWithDisp/TunnelSamplingDisp "${args[@]}"

#calculate the estimates

args=(
	$prop_simple_filename
	$posterior_filename
	$ref_dist_filename
	$additional_ref_dist_filename
	)

java -cp "./dist/BridgingInTreeSpace.jar" MarginaLikelihoodCalculationsWithDisp/ChibJeliEstimateDisp "${args[@]}"  >> ${folder_name}${file_prefix}_ChibEst.txt #replace with file name for storing the Chib estimate

java -cp "./dist/BridgingInTreeSpace.jar" MarginaLikelihoodCalculationsWithDisp/TunnelSamplingEstimateDisp "${args[@]}"  >> ${folder_name}${file_prefix}TunnelEst.txt #replace with file name for storing the tunnel estimate

```
</details>



### Marginal likelihood for unknown dispersion - Stepping Stone

<details>
<summary>⭐ show/hide</summary>
Given a set of unrooted phylogenetic trees in Newick string format in a file called Example_trees.txt, within a folder called Example_folder, to estimate the marginal likelihood using the stepping stone method where the dispesion is consider an unknown nuisance parameter and a fixed source tree is given in the file Example_source.txt, a .sh script of the form specified below can be run. Parameters followed by a comment can be modifed as required.

```bash
#!/bin/bash

 file_prefix="Example_trees"
 folder_name="./Example_folder/"
 source_tree_filename="Example_source.txt"
 data_filename="${folder_name}${file_prefix}.txt" # data filename
 source_tree_filename="${folder_name}${source_tree_filename}" # data filename
 posterior_filename="${folder_name}${file_prefix}_step_stone_posterior.txt"
 ref_dist_filename="${folder_name}${file_prefix}_step_stone_reference.txt"
 stepping_dist_filename="${folder_name}${file_prefix}_step_stone_disp" # output filename for the samples from the posterior

args=(
        $data_filename
        $source_tree_filename
        "0" # Initial squ root t_0 - set to zero to use the Frechet variance -- recommended to use FV - short initial posterior run
        "20" # Num steps in the bridges
        "1802" # Set the seed for the random number generator - short initial posterior run
        "-n" # 
        "1000" # Num interations - before thin - short initial posterior run
        "-t" # 
        "10" # thin - short initial posterior run
        "-b" # 
        "1000" # burn-in - short initial posterior run
        "-o" # 
        $posterior_filename #posterior output files  - short initial posterior run
        "-pbg" # 
        "0.05" # geometric length bridge prop - short initial posterior run
	"-ptr" #
        "0.5" # log random walk proposal parameter - short initial posterior run
	$ref_dist_filename # Output file for calculatating normalising constant of the reference distribution
	"20" # num of values of dispersion to use in numerical integration for reference dist
        "50" # num of proopsals to calculate the normalising constant of the reference distribution
	"FV" # write FV to use the Frechet variance about x0 or otherwise specidy a number -- recommended to use FV - stepping stone run
	"490" # Seed - stepping stone run
	"-n" #  
	"100" # Num interations - before thin - stepping stone run
	"-t" # 
	"2" # thin - stepping stone run
	"-b" # 
	"1000" # burn-in - stepping stone run
	"-o" # # output file name for the MCMC outputs from each step on the path
	$stepping_dist_filename # 
	"-pbg" # 
	"0.01" # geometric length bridge prop - stepping stone run
	"-ptr" # 
	"0.5" # t0 random walk proposal - stepping stone run
        ) 
        	
 
java -cp "./dist/BridgingInTreeSpace.jar" MarginalLikelihoodsWithDisp/SteppingStoneSamplingDisp "${args[@]}"

args=(
	$ref_dist_filename
 	$stepping_dist_filename
 	)
        
java -cp "./dist/BridgingInTreeSpace.jar" MarginaLikelihoodCalculationsWithDisp/StepStoneEstimateDisp "${args[@]}" >> ${folder_name}${file_prefix}StepStoneEst.txt #file name for storing the Chib estimate

```
</details>


### Marginal likelihood for unknown dispersion - Chib two block

<details>
<summary>⭐ show/hide</summary>
Given a set of unrooted phylogenetic trees in Newick string format in a file called Example_trees.txt, within a folder called Example_folder, to estimate the marginal likelihood using the Chib two block method where the dispersion is considered an unknown nuisance parameter and a fixed source tree given in the file Example_source.txt, a .sh script of the form specified below can be run. Parameters followed by a comment can be modifed as required.

```bash
#!/bin/bash
 
 file_prefix="Example_trees"
 folder_name="./Example_folder/"
 source_tree_filename="Example_source.txt"
 data_filename="${folder_name}${file_prefix}.txt" # data filename
 source_tree_filename="${folder_name}${source_tree_filename}" # data filename
 fixed_disp_posterior_filename="${folder_name}${file_prefix}_Chib_TwoBlock_Post_Fixed.txt" # output filename for the samples from the posterior with fixed dispersion
 props_filename="${folder_name}${file_prefix}_Chib_Two_Block_Props.txt" # output filename for the samples from the proposals
 full_posterior_filename="${folder_name}${file_prefix}_Chib_Two_Block_Post_Full.txt" # output filename for the samples from the proposals

#run the MCMC and independence proposals
args=(
        $data_filename
        $source_tree_filename
        "0.2" # t_0 proposal param in the middle distribution
        "20" # Num steps
        "1260" # Seed
        "-n" #
        "10000" # Num MCMC interations - before thin
        "-t" #
        "100" # thin
        "-b" #
        "1000" # burn-in
        "-o" # 
        $fixed_disp_posterior_filename # output file for the posterior
        "-pbg" # 
        "0.05" # partial bridge proposal parameter
        "-numProps" # 
        "50000" # Num independence proposals to run
        $props_filename # output file for the proposals
	    "897" # Seed - full posterior sims
        "-n" # 
        "25000" # Num interations - before thin - full posterior sims
        "-t" # 
        "20" # thin - full posterior sims
        "-b" # 
        "5000" # burnin - full posterior sims
        "-o" #
        $full_posterior_filename # full posterior file name
        "-pbg" # 
        "0.01" # partial bridge proposal with geometric length - full posterior sims
        "-ptr" # 
        "0.5" # t0 log random walk proposal - full posterior sims
        )
        
java -cp "./dist/BridgingInTreeSpace.jar" MarginalLikelihoodsWithDisp/ChibTwoBlock "${args[@]}"


#calculate the estimates
args=(
	$fixed_disp_posterior_filename
	$props_filename
	$full_posterior_filename
)

java -cp "./dist/BridgingInTreeSpace.jar" MarginaLikelihoodCalculationsWithDisp/ChibTwoBlockEstimate "${args[@]}"  >> ${folder_name}${file_prefix}_ChibTwoBlockEst.txt #replace with file name for storing the Chib estimate

```
</details>


### Noisy MCMC for topology inference

<details>
<summary>⭐ show/hide</summary>

Suppose we are given a set of unrooted phylogenetic trees in Newick string format in a file called Example_trees.txt, within a folder called Example_folder. To run the noisy MCMC algorithm for inferring the species tree when gene trees are modelled as draws from a Brownian motion kernel (approximated by a random walk kernel), run a .sh script of the following form (note that in this analysis we manually input an initial tree):

```bash
#!/bin/bash
##filenames
file_prefix="Example_trees"
folder_name="./Example_folder/"
source_tree_filename="Example_source.txt"
data_filename="${folder_name}${file_prefix}.txt" # data filename
output_filename="${folder_name}${file_prefix}_MCMCOutput.txt"
initial_tree_filename="${folder_name}${source_tree_filename}" #wouldn't usually be the source tree
topologies_filename="${folder_name}${file_prefix}_MCMCOutput_tops.txt"

#Run the MCMC:
args=(
	$data_filename # data filename
	$initial_tree_filename # initial tree filename
	$output_filename # output file for the MCMC
	$topologies_filename # a file to count the topologies in the MCMC output into
	"1" # Initial sql root t_0 
	"50" # Num steps in the random walks
        "1000" # Number of particles in the approx likelihood calcs
	"1" # Number of cores to use in forward simulating particles
	"1500" # Num interations
	"0" #burnin
	"100" # thin
	"0.1" # parameter for the log random walk edge proposal
	"0.1" # parameter for the random walk edge proposal
)

java -cp "./dist/BridgingInTreeSpace.jar" topologies/InferParamsViaApproxLike "${args[@]}"
	
```


</details>

### Simulating from the Gaussian kernel distribution using MCMC
<details>
<summary>⭐ show/hide</summary>
Given an unrooted phylogenetic tree in Newick string format in a file called Example_tree.txt, within a folder called Example_folder, we can simulate from a Gaussian kernel, as in the paper [kdetrees: non-parametric estimation of phylogenetic tree distributions](https://academic.oup.com/bioinformatics/article/30/16/2280/2748204) a .sh script in the following form, where the parameters followed by a comment can be modifed as required:

```bash
#!/bin/bash
##filenames
file_prefix="Example_tree"
folder_name="./Example_folder/"
tree_filename="${folder_name}${file_prefix}.txt" # data filename
output_filename="${folder_name}${file_prefix}_MCMCOutput.txt"
topologies_filename="${folder_name}${file_prefix}_MCMCOutput_tops.txt"
geodesic_distances_filename="${folder_name}${file_prefix}_MCMCOutput_distnaces.txt"

#Run the MCMC:
args=(
	$data_filename # Tree x0 the 'centre' of the distribution
	$output_filename # output file for the simulated trees
	$topologies_filename # output file for the topologies of the simulated trees
	$geodesic_distances_filename # output file for the geodesic distances of the simulated trees to the source        
	"947" # seed for the random engine
	"1.0" # value of t0 in the kernel
	"0.09" # value for the proposal distribution which is a random walk proposal
	"10000" # number of iterations in the MCMC  
	"1000" # number of burn-in iterations in the MCMC    
	"10" # number of thin iterations in the MCMC  
)

java -cp "./dist/BridgingInTreeSpace.jar" bridge/GDKsimulationMain "${args[@]}"

```

A minimal example R script for reproducing the plots in the thesis that show the spread of topologies in the simulated trees and the distribution of distances to the source is given in the following:

```r
Nprime<-2
source_tree_prefix<-"Example_source"
dists<-read.delim(paste0(folder,"/", source_tree_prefix,"_GKD_Sims_distances.txt"),header=FALSE)
chiComp<-function(x){
  return(dchisq(x,Nprime))
}

ggplot(dists)+geom_density(aes(x=as.numeric(V1)^2,colour='GKD'))+geom_function(aes(x=V1,colour='chiSqu'),fun = chiComp)+
  xlab("distance")+ylab('density')+xlim(0,10)

#plot the topologies in the ouput:
topsData<-read.delim(paste0(folder,"/", source_tree_prefix,"_GKD_Sims_tops.txt"),header=FALSE,sep=" ")
colnames(topsData)<-c("Topology","Count")
topsData$Count = as.numeric(topsData$Count)
topsData<-topsData %>% mutate(Props = topsData$Count/sum(topsData$Count))
topsData<-topsData[order(topsData$Count,decreasing=TRUE),]

ggplot(topsData)+geom_col(aes(x=Topology,y=Props))+scale_x_discrete(limits=topsData$Topology) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8))+
  xlab('Topology')+ylab('Proportion')

```

</details>


