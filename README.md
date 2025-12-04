# Bridging in tree space

## Table of Contents
- [Installation](#installation)
- [Summary](#summary)
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
  - [Assessing the number of random walk steps by forward simulation](#assessing-the-number-of-random-walk-steps-by-forward-simulation)
  - [Tuning the bridge proposal](#tuning-the-bridge-proposal)
    

## Summary

This Github repository includes the code required to perform the statistical inference procedures specified in the paper [Brownian motion, bridges and Bayesian inference in phylogenetic tree space](https://arxiv.org/abs/2506.22135) and my upcoming thesis. Particularly we perform inference in [Billera-Holmes-Vogtmann (BHV) tree space](https://www.sciencedirect.com/science/article/pii/S0196885801907596).

We perform inference for Brownian motion kernels in BHV tree space, as an analogue to the Gaussian distribution in Euclidean space. Much of the methodology relies on simulating Brownian bridges, which are noisy walks conditioned on a source tree and a target tree, in BHV tree space.

A html file called ExperimentalDataRScript and a corresponding RMarkdown file for reproducing the analysis on yeast gene trees detailed in the paper [Brownian motion, bridges and Bayesian inference in phylogenetic tree space](https://arxiv.org/abs/2506.22135) are provided as part of this project. <a href="https://will-woodman.github.io/BridgingInTreeSpace/ExperimentalDataRScript">ExperimentalDataRScript</a>. The html file details how to run the analysis and the RMarkdown file can be used to reproduce the plots seen in the paper. The data required to run the yeast gene tree analysis is provided in the folder YeastData.

The simulated data sets used to perform simulation studies in the thesis are given in this repository in the folder SimulatedDataForThesis.

## Installation

The code can be downloaded from this Github repository and run from the command line. In particular the dist folder, which contains all relevant .jar files should be downloaded from this repository. The code has been testing in Java version 17. The easiest way to run the code is with shell script files, where we can specify the various parameters required for the MCMC algorithms. The rest of this readme gives example shell scripts and parameters for running the different MCMC algorithms, where modifiable parameters are clearly commented.

We alse give minimal examples of R scripts used to produce the analysis of MCMC runs. However, it is likely that any future analysis of MCMC runs will be bespoke so these scripts are non-exhaustive.

## Running the inference procedures

### Simulating a data set under the model
<details>
<summary>⭐ show/hide</summary>
This file will simulate a random source tree called Example_source in a folder called Example_folder with a specifed number N of taxa from a coalescent model (with edges resampled from a Gamma distributionwith shape-scale parametrisation). It will then produce a data set of n trees by simulating random walks with m steps, started at the source tree and a specified value of dispersion. It also outputs the topologies of the source tree and data as well as geodesic and Robinson Foulds distances from the source tree to the data points.

```bash
#!/bin/bash
##define the filenames:
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
	"50" # number  of trees to simulate in the data set
	"947" # seed for the random engine
	"100" # number of steps in the random walk
	"2.0" # shape of the gamma distribution for simulating source tree edges
	"0.05" # scale of the gamma distribution for simulating source tree edges
	"0.25" # value of t0 in the random walks
)

java -cp "./dist/BridgingInTreeSpace.jar" simulateTops/simulateSourceTree "${args[@]}"
```


A minimal R script to perform initial analysis of the data set is given in the following. This includes plotting the Robinson Foulds and geodesic distances from the source tree to the data points to understand the dispersion in the data set.

```R

##initial analysis of the data set:
print("spread of the data over the topologies")
folder <- "./Example_folder"
data_file_prefix <- "Example_Trees"
data_filename<-paste0(folder,"/",data_file_prefix,".txt",sep="")
tops_filename<-paste0(folder,"/", data_file_prefix,"_tops.txt")
distances_filename<-paste0(folder,"/", data_file_prefix,"_distances_to_source.txt")
data_filename

Ntaxa<-5

topsData<-read.delim(tops_filename,sep=" ",header=FALSE)
colnames(topsData)<-c("Topology")
topsData<-topsData %>%count(Topology)
colnames(topsData)<-c("Topology","Count")

##plot counts of each topology in the data set
ggplot(topsData)+geom_col(aes(x=Topology,y=Count))+scale_x_discrete(limits=topsData$Topology) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

print("Distribution of Robinson-Foulds distances from data to source tree")
RFDData<-read.delim(distances_filename,sep=" ",header=TRUE)
RFPlotPt<-ggplot(RFDData)+geom_bar(aes(x=RobinsonFouldsDist)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle("Robinson Foulds distances")+
  xlab("Robinson Foulds Distance")+ylab("Count")+theme(legend.position="none",plot.title = element_text(size=11,hjust=0.5),
                                                       axis.title=element_text(size=8.5),axis.text=element_text(size=7.5))+scale_x_continuous(breaks = seq(0, (2*(Ntaxa-3)), by = 2))
RFPlotPt

print("Distribution of Geodesic distances from data to source tree")
GDData<-read.delim(distances_filename,sep=" ",header=TRUE)
GDPlotPt<-ggplot(RFDData)+geom_bar(aes(x=GeodesicDist)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle("Geodesic distances")+
  xlab("Geodesic Distance")+ylab("Count")+theme(legend.position="none",plot.title = element_text(size=11,hjust=0.5),
                                                       axis.title=element_text(size=8.5),axis.text=element_text(size=7.5))

GDPlotPt

```

</details>


### Posterior inference using bridges
<details>
<summary>⭐ show/hide</summary>
Given a set of unrooted phylogenetic trees in Newick string format in a file called Example_trees.txt, within a folder called Example_folder, posterior inference can be run using a .sh script in the following form, where the parameters followed by a comment can be modifed as required. The script will count the topologies in the posterior at various points in order to plot the cumulative proportions of the different topologies. It will also output the edge split lengths in the modal topology in the posterior which are used to calculate modal trees and plot kernel density estimates to understand the uncertainty in the posterior for the source tree.

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
	"0" # Initial squ root t_0 - set to zero if you want to use the Frechet variance (recommended)
	"50" # Num steps
	"1051" # Seed
	"-n" #
	"500000" # Num interations
	"-t"
	"10" # thin
	"-b"
	"50000" #burnin
	"-o"
	$output_filename #output file
	"-pbg"
	"0.1" #parameter for partial bridge proposal
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

Methods to plot, summarise and assess the output of the MCMC runs are given in the following R script. This includes printing the various acceptance rates for MCMC proposals, trace plots and a plot of cumulative proportions of topologies in the posterior to assess convergence, and kernel density estimates of edge lengths in the modal topology in the posterior.

<details>
<summary>⭐ show/hide</summary>

```R
##plot counts of each topology in the data set
ggplot(topsData)+geom_col(aes(x=Topology,y=Count))+scale_x_discrete(limits=topsData$Topology) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

print("Distribution of Robinson-Foulds distances from data to source tree")
RFDData<-read.delim(distances_filename,sep=" ",header=TRUE)
RFPlotPt<-ggplot(RFDData)+geom_bar(aes(x=RobinsonFouldsDist)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle("Robinson Foulds distances")+
  xlab("Robinson Foulds Distance")+ylab("Count")+theme(legend.position="none",plot.title = element_text(size=11,hjust=0.5),
                                                       axis.title=element_text(size=8.5),axis.text=element_text(size=7.5))+scale_x_continuous(breaks = seq(0, (2*(Ntaxa-3)), by = 2))
RFPlotPt

print("Distribution of Geodesic distances from data to source tree")
GDData<-read.delim(distances_filename,sep=" ",header=TRUE)
GDPlotPt<-ggplot(RFDData)+geom_bar(aes(x=GeodesicDist)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle("Geodesic distances")+
  xlab("Geodesic Distance")+ylab("Count")+theme(legend.position="none",plot.title = element_text(size=11,hjust=0.5),
                                                       axis.title=element_text(size=8.5),axis.text=element_text(size=7.5))

GDPlotPt

##Process output of inference:
MCMCoutput_filename<-paste0(folder,"/", data_file_prefix,"_MCMCOutput.txt")
trueT0<-0.25 #obviously only use this if this is simulated data where you know the truth

#Read in the data and get the acceptance rates:
x0Data<-read.delim(MCMCoutput_filename,sep=" ")
x0AcceptanceRate<-x0Data[length(x0Data$source),2]
bridgeAcceptanceRates<-x0Data[!is.na(as.numeric(x0Data$source)),]
Its<-400 ##in this instance after thinning
thin<-1
x0Data<-x0Data[1:Its,]
x0Data$dispersion<-as.numeric(x0Data$dispersion)


##t0 plots:
t0TracePlotPt01<-ggplot(x0Data)+geom_line(aes(y=dispersion,x=(c(1:Its)*thin), size=0.5)+
  xlab("Iteration")+ylab("Dispersion")+theme(legend.position="none",plot.title = element_text(size=22,hjust=0.5),
                                             axis.title=element_text(size=21),axis.text=element_text(size=16.5),axis.title.y = element_text(margin = margin(r=10)))+xlim(c(0,(Its*100+300000))))+
  ggtitle("Traceplot for dispersion")

t0KDEPlotPt01<-ggplot(x0Data)+geom_density(aes(dispersion))+geom_vline(aes(xintercept=trueT0Pt01))+ggtitle("Kernel density estimate for dispersion")+
  xlab("Dispersion")+ylab("Density")+theme(legend.position="none",plot.title = element_text(size=11,hjust=0.5),
                                           axis.title=element_text(size=8.5),axis.text=element_text(size=7.5))##t0KDE

t0KDEPlotPt01

t0AcceptanceRate<-bridgeAcceptanceRates$source[length(bridgeAcceptanceRates$source)]
averageBridgeAcceptanceRate<-mean(as.numeric(bridgeAcceptanceRates$source[1:(length(bridgeAcceptanceRates$source)-1)]))

print(paste("x0 acceptance rate ",x0AcceptanceRate))
print(paste("t0 acceptance rate ",t0AcceptanceRate))
print(paste("average bridge acceptance rate ",averageBridgeAcceptanceRate))

##make a table of the acceptance rates
acceptanceRates<-data.frame(c("x0","t0","bridges"),c(x0AcceptanceRate,t0AcceptanceRate,averageBridgeAcceptanceRate))
colnames(acceptanceRates)<-c("Parameter","Rate")
acceptanceRates


x0Data<-NULL

##now we plot cumulative proportion of each topology in the posterior in order to assess quality of the MCMC ouput:

props_filename<-paste0(folder,"/", data_file_prefix,"_MCMCOutput_tops.txt")
theProps<-read.delim(props_filename,sep=" ")
theProps<-subset(theProps,select=-c(X))

its<-Its/100 ##topologies are counted at every 100 points in the posterior sample

theCols<-colnames(theProps)
theCols<-theCols[-(1)]

for(i in 1:length(theCols))
{
  aColname<- theCols[i]
  theProps[[aColname]]=theProps[[aColname]]/sum(theProps[[aColname]])
}

theProps<-theProps[,c(1,(its+2)-seq(0:(its-1)))]

##change this so it finds the name of the last column itself
theProps<-theProps[order(theProps[[paste0("X",Its*thin,sep="")]],decreasing=TRUE),]
PropsToPlot<-transpose(theProps)
tmydf = setNames(data.frame(t(theProps[,-1])), theProps[,1])
tmydf<-melt(data.table(tmydf))
colnames(tmydf)<-c("topology","proportion")
tmydf[, Iteration := rep(c(1:its), length.out = .N)]

print("Plot the proportion of the top topologies in the posterior:")
TopPropPlotPt<-ggplot(tmydf)+geom_line(aes(x=Iteration*10^4,y=proportion,group=topology))+ggtitle("$t_0=0.01$")+
  xlab("Iteration")+ylab("Proportion")+theme(legend.position="none",plot.title = element_text(size=11,hjust=0.5),
                                             axis.title=element_text(size=8.5),axis.text=element_text(size=7.5))


##now plot kernel density estimates for the edge lengths in the modal topology in the MCMC output
# I have excluded the code to plot the true edge lengths because these example scripts will likely be run on experimental 
# rather than simulated data
edge_filename<-paste0(folder,"/", data_file_prefix,"_MCMCOutput_edges.txt")

theDensity1<-read.delim(edge_filename,sep=" ")

#extract the different splits from the header:
theSplits<-head(read.delim(edgeFilename,header=FALSE,sep="]"),1)
theSplits[1]<-substr(theSplits[1],3,str_length(theSplits[1]))
theSplits<-paste0(theSplits,"]")
theSplits<-head(theSplits,length(theSplits)-1)

nPrime=Ntaxa-3
theDensity1<-theDensity1[]
colnames(theDensity1)<-theSplits##need a better way of doing this...
theDensity1<-theDensity1[,1:nPrime]
theDensity1<-theDensity1[apply(theDensity1,1,min)!=0,]
theDensity <- melt(data.table(theDensity1))

colnames(theDensity)<-c("Split","Length")
print("Plot kernel density estimates of the edge lengths in the posterior")
ggplot(theDensity)+geom_density(aes(Length,color=Split))+ggtitle("KDEs of split lengths in the modal topology in the posterior")

##make a table showing the topologies in the posterior sample, highlighting the true source tree topology (modify if there is no "true" source)
print("Top topologies in the posterior, * denotes source tree topology")
source_tree_prefix<-"Example_tree"
source_tree_top<-read.delim(paste0(folder,"/", source_tree_prefix,"_top.txt"),header=FALSE)[1,1]


topsInPosterior<-data.table(theProps[,1],theProps[,length(colnames(theProps))])
colnames(topsInPosterior)<-c("Topology","Proportion")
topsInPosterior$Topology<-paste(topsInPosterior$Topology,ifelse(topsInPosterior$Topology==sourceTreeTop,"*",""))


topsInPosterior$Proportion<-paste(round(100*as.numeric(topsInPosterior$Proportion), 1), "%", sep="")
print(topsInPosterior)
```

</details>

If you have built the mode tree using R you will have outputted a list of splits with their lengths. We can turn this into a phylogenetic tree in Newick string form using the following shell script. The tree can then be used for example as an alternative to the Frechet mean.

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
        "100000" # Num MCMC interations - before thin
        "-t" #
        "100" # thin
        "-b" #
        "10000" # burn-in
        "-o" # 
        $posterior_filename # output file for the posterior
        "-pbg" # 
        "0.05" # partial bridge proposal parameter
        "-numProps" # 
        "5000" # Num independence proposals to run
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
        "10000" # Num interations - before thin
        "-t" # 
        "20" # thin
        "-b" # 
        "10000" # burn-in
        "-o" # 
        $posterior_filename #posterior output files 
        "-pbg" # 
        "0.05" # geometric length bridge prop
        "2000" # number of Proposals to run
        "0.01" # the first non zero value of beta_k (keep this fixed at 0.01)
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

To obtain the samples for a Chib one block or tunnel estimate, we first simulate from the posterior distribution (on bridges and t0). We use the marginal posterior sample for t0 to parametrise a reference distribution. We then estimate the normalising constant of the reference distribution by directly simulating independence proposals and finally simulate from the reference distribution using MCMC. The following shell script allows different parameters for the two different MCMC runs in case the acceptance rates differ significantly between the two distributions.

```bash
#!/bin/bash
 
 file_prefix="Example_trees"
 folder_name="./Example_folder/"
 source_tree_filename="Example_source.txt"
 data_filename="${folder_name}${file_prefix}.txt" # data filename
 source_tree_filename="${folder_name}${source_tree_filename}" # x0 filename
 posterior_filename="${folder_name}${file_prefix}_Chib_Post.txt" # output filename for the samples from the posterior
 ref_dist_filename="${folder_name}${file_prefix}_Chib_Props.txt" # output filename for the samples from the reference distribution 
 ref_dist_parameters_filename="${folder_name}${file_prefix}_Chib_Ref_Dist_Params.txt" # file to store the parameters of the lognormal reference distribution on t0
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
        "100000" # Num MCMC interations - before thin
        "-t" #
        "10" # thin
        "-b" #
        "10000" # burn-in
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
        "100000" # Num interations - before thin - ref dist MCMC sims
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
        "1000" # number of proposals per data point per value of dispersion
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
Given a set of unrooted phylogenetic trees in Newick string format in a file called Example_trees.txt, within a folder called Example_folder, to estimate the marginal likelihood using the stepping stone method where the dispersion is considered an unknown nuisance parameter and a fixed source tree is given in the file Example_source.txt, a .sh script of the form specified below can be run. Parameters followed by a comment can be modifed as required. 

In order to parametrise a reference distribution for t0, we first simulate from the posterior distribution on bridges and dispersion. We then loop through different distributions on the geometric path between the reference distribution and the posterior.

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

Suppose we are given a set of unrooted phylogenetic trees in Newick string format in a file called Example_trees.txt, within a folder called Example_folder. To run the noisy MCMC algorithm for inferring the species tree when gene trees are modelled as draws from a Brownian motion kernel (approximated by a random walk kernel), run a .sh script of the following form (note that in this analysis we manually input an initial tree). The MCMC is noisy because the likelihood is approximated by forward simulating a number of random walks. The accuracy in the approximation of the likelihood can be modified by increasing the number of forward simulated particles. This methodology is highly parallelisable and can be done parallelised by changing the number of cores used to forward simulate particles.

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
Given an unrooted phylogenetic tree in Newick string format in a file called Example_tree.txt, within a folder called Example_folder, we can simulate from the Gaussian kernel from the paper [kdetrees: non-parametric estimation of phylogenetic tree distributions](https://academic.oup.com/bioinformatics/article/30/16/2280/2748204) using Metropolis-Hastings MCMC. The methodolgy is specified in the thesis. We can do this using
a .sh script in the following form, where the parameters followed by a comment can be modifed as required:

<details> 
<summary>⭐ show/hide</summary>	

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


### Assessing the number of random walk steps by forward simulation
<details>
<summary>⭐ show/hide</summary>
In this section we give details on how to reproduce the analysis in the thesis on how many random walk steps are required to approximate a Brownian motion kernel by forward simulating random walks with different numbers of steps, using both the standard random walk innovations and the innovations defined using the exact distribution on the three spider.

```bash
#!/bin/bash
folder_name="./Example_folder/"
tree_filename="${folder_name}_simulated_source_tree.txt" # file containing list of splits and lengths
exact_dist_counts_filename="${folder_name}topology_counts_exact.txt" # file containing list of splits and lengths
standard_dist_counts_filename="${folder_name}topology_counts_standard.txt" # file containing list of splits and lengths

#build and save the modal tree:
args=(
        $exact_dist_counts_filename # modal tree output file
        $tree_filename # modal splits
        "10" #number of taxa for the source tree
	"10 50 100" #number of steps to use in the random walks
	"1.0" # value of dispersion to use in the random walks
        "10000" # number of random walks to simulate for each number of steps
        "2.0" # shape parameter in gamma dist for source tree edge lengths
        "1.4" # scale parameter in gamma dist for source tree edge lengths
        "true" # whether to sample the steps using the exact dist on the three spider
        "104" # seed for the random number generator 
        )
        
java -cp "./dist/BridgingInTreeSpace.jar" simulateMultipleTopologies.investigateNoOfWalkStepsLoop "${args[@]}"

args=(
        $standard_dist_counts_filename # modal tree output file
        $tree_filename # modal splits
        "10" #number of taxa for the source tree
	"10 50 100" #number of steps to use in the random walks
	"1.0" # value of dispersion to use in the random walks
        "10000" # number of random walks to simulate for each number of steps
        "2.0" # shape parameter in gamma dist for source tree edge lengths
        "1.4" # scale parameter in gamma dist for source tree edge lengths
        "False" # whether to sample the steps using the exact dist on the three spider
        "104" # seed for the random number generator -- leave the same to use the same tree
        )
        

java -cp "./dist/BridgingInTreeSpace.jar" simulateMultipleTopologies.investigateNoOfWalkStepsLoop "${args[@]}"

```

A minimal example R script for reproducing the plots in the thesis that show the spread of topologies in the simulated trees and the distribution of distances to the source is given in the following:

```r
##read in and plot the count data for the standard innovations
countData<-read.delim("Example_folder/topology_counts_standard.txt",header=TRUE,sep=" ")
ggplot(countData)+geom_col(aes(x=Steps,y=Topologies)) + ggtitle("Number of topologies by steps standard distribution")

##read in and plot the count data for the innovations that use the exact distribution on the three spider
countData<-read.delim("Example_folder/topology_counts_exact.txt",header=TRUE,sep=" ")
ggplot(countData)+geom_col(aes(x=Steps,y=Topologies)) + ggtitle("Number of topologies by steps exact distribution")

```

</details>

### Tuning the bridge proposal
<details>
<summary>⭐ show/hide</summary>
In this section we give details on how to simulate a fixed number of iterations of MCMC for a bridge between two fixed endpoints in tree space. However, outputing the topology at each step of the bridge and changing the parameters used in the bridge proposal requires directly modifying the code and so anyone who wanted to repeat the full analysis would need to create a project in an IDE such as Netbeans including the code provided in this repository and make manual changes to the code. Anyone who wants to do this would be advised to contact the author.

```bash
#!/bin/bash
##filenames
file_prefix="Example_trees"
folder_name="./Example_folder/"
source_filename="${folder_name}$SourceTree1.txt" # data filename
target_filename="${folder_name}$TargetTree1.txt" # data filename
output_filename="${folder_name}$Trees1_NoPen_MCMCOutput.txt"

#Run the MCMC:
args=(
	$target_filename #target tree filename
	$source_filename #source tree filename
	"0.5" # Initial squ root t_0
	"50" # Num steps
	"1051" # Seed
	"-n" #
	"400000" # Num interations
	"-t"
	"100" # thin
	"-b"
	"10000" #burnin
	"-o"
	$output_filename #output file
	"-pbg"
	"0.01" #parameter for partial bridge proposal
)

java -cp "./dist/BridgingInTreeSpace.jar" bridge/InferBrownianParamsMCMC "${args[@]}"


```

An example R script for reproducing the plots in the thesis that show the number of topologies in the simulated trees on each step of the bridges is given in the following. This includes printing out the acceptance rates for each MCMC run, producing trace plots for the log likelihood to check issues with convergence and producing the plots showing the different numbers of topologies at each bridge step.

```r
###R script for processing the output of the independence proposal sims
#function for getting the number of unique topologies at each step:
getTheCounts<-function(data,m,name){
  counts<-c(3:(m+2))
  for(i in 3:(m+2)){
    counts[i-2]=length(unique(data[,i]))
    
  }
  counts<-data.frame(c(1:50),counts)
  colnames(counts)<-c("steps","counts")
  counts$name<-name
  return(counts)
  
}

#/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/IndepPropTuning/10Taxa20240910/Trees1_Variable_MCMCoutv1.txt

m<-50
noOfBridges<-4000

##read in the different data sets and print out the acceptance rates in each MCMC run:
theData1<-read.delim("/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/IndepPropTuning/10Taxa20240910/Trees1_1TimesPenalty_MCMCoutv1.txt",sep=" ",header=FALSE)
print(paste0("Acceptance rate for 1 times penalty: ",theData1[4002,10]))

theData2<-read.delim("/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/IndepPropTuning/10Taxa20240910/Trees1_Variable_MCMCoutv1.txt",sep=" ",header=FALSE)
print(paste0("Acceptance rate for 2 times penalty: ",theData2[4002,10]))

theData3<-read.delim("/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/IndepPropTuning/10Taxa20240910/Trees1_3TimesPenalty_MCMCoutv1.txt",sep=" ",header=FALSE)
print(paste0("Acceptance rate for 1 times penalty: ",theData3[4002,10]))

theData4<-read.delim("/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/IndepPropTuning/10Taxa20240910/Trees1_NoPen_MCMCOutput.txt",sep=" ",header=FALSE)
print(paste0("Acceptance rate for no penalty: ",theData4[4002,10]))

##plot traceplots of the log likelihood for each chain to check that there are no issues with convergence
theLikelihood<-data.frame(cbind(c(1:noOfBridges),theData1[2:(noOfBridges+1),53]))
colnames(theLikelihood)<-c("Iteration","Likelihood")

print("Plot of the log likelihood for 1 times penalty")

likePlot1<- ggplot(theLikelihood)+geom_line(aes(x=Iteration,y=Likelihood))+ggtitle("1 times penalty")+
  theme(plot.title = element_text(hjust = 0.5))
likePlot1

theLikelihood<-data.frame(cbind(c(1:noOfBridges),theData2[2:(noOfBridges+1),53]))
colnames(theLikelihood)<-c("Iteration","Likelihood")

print("Plot of the log likelihood for 2 times penalty")

likePlot2<- ggplot(theLikelihood)+geom_line(aes(x=Iteration,y=Likelihood))+ggtitle("2 times penalty")+
  theme(plot.title = element_text(hjust = 0.5))
likePlot2

theLikelihood<-data.frame(cbind(c(1:noOfBridges),theData3[2:(noOfBridges+1),53]))
colnames(theLikelihood)<-c("Iteration","Likelihood")

print("Plot of the likelihood for 3 times penalty")

likePlot3<- ggplot(theLikelihood)+geom_line(aes(x=Iteration,y=Likelihood))+ggtitle("3 times penalty")+
  theme(plot.title = element_text(hjust = 0.5))
likePlot3

theLikelihood<-data.frame(cbind(c(1:noOfBridges),theData4[2:(noOfBridges+1),53]))
colnames(theLikelihood)<-c("Iteration","Likelihood")

print("Plot of the log likelihood for 4 times penalty")

likePlot4<- ggplot(theLikelihood)+geom_line(aes(x=Iteration,y=Likelihood))+ggtitle("4 times penalty")+
  theme(plot.title = element_text(hjust = 0.5))
likePlot4

##count the number of distinct topologies at each step on the bridge for each MCMC run:
theCounts1<-getTheCounts(theData1,m,"1 times")
theCounts2<-getTheCounts(theData2,m,"2 times")
theCounts3<-getTheCounts(theData3,m,"3 times")
theCounts4<-getTheCounts(theData4,m,"4 times")

theCounts<-rbind(theCounts1,theCounts2,theCounts3,theCounts4)

print("Counts of the distinct topologies in the posterior sample by step on the bridge (we use 50 steps in the bridges)")

#final plot:
plotout<-ggplot(theCounts)+geom_line(aes(y=counts,x=steps,linetype=name),position="dodge",size=0.6)+xlab("Step")+ylab("Count")
plotout<-plotout+theme(axis.title=element_text(size=13),axis.text=element_text(size=11),legend.key.size = unit(1.2, 'cm'),legend.text = element_text(size=11),legend.title=element_blank() )
plotout <- plotout + labs(linetype="Tree pair")
plotout

```

</details>
