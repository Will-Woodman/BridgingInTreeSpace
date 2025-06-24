/*
 * InferParamsViaApproxLike.java

    Copyright (C) 2018  Tom M. W. Nye

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact the author at:  <tom.nye@ncl.ac.uk>
                            <http://www.mas.ncl.ac.uk/~ntmwn/>
 */

package topologies;

import MCMC.Chain;
import MCMC.DensityCalculator;
import MCMC.GlobalState;
import MCMC.Kernel;
import MCMC.KernelWrapper;
import MCMC.MetropolisHastingsKernelWrapper;
import MCMC.PositiveParameter;
import MCMC.PosteriorAnalysis;
import static MCMC.PosteriorAnalysis.countTopologies;
import MCMC.State;
import MCMC.SweepKernelWrapper;
import diffbase.TreeState;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import treebase.AlgorithmException;
import treebase.Tree;
import treebase.TreeAsSplits;
import topologies.TopologicalData;

/**
 * Infer params for a set of topologies -- this file loops through different parameters to tune the proposals
 */
public class InferParamsViaApproxLike {
    
    public static void main(String[] args) throws IOException {

        /* All inputs required for run */

        // Data set
//        String dataFilename = "/home/tom/research/Diffusion/bridge/sims/sample5_0.05.txt";
        // Initial tree and sigma
//        String initialTreeFilename = "/home/tom/research/Diffusion/bridge/sims/five_taxon_tree.txt";
        
        double sqrtt0=1.0;
        double fixedVariance = sqrtt0*sqrtt0;

        // Data set
        String dataFilename = "/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/5 Taxon Gene Trees/coal_trees_0.6_0.01_0.01.txt";
        // Initial tree and sigma
        String initialTreeFilename = "/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/5 Taxon Gene Trees/SpeciesTreeTopInitialGuess.txt";
        
        
        /* Read in the data */
        File inFile = new File(dataFilename);
        TopologicalData theData = null;
        try {
            theData = new TopologicalData(inFile);
        }
        catch (java.io.IOException anError) {
            System.out.println("Bad input file. "+anError.getMessage());
            System.exit(1);
        }

        /* Read in an initial tree */
        Tree initialTree = null;
        try {
            initialTree = new Tree(new File(initialTreeFilename));
            initialTree.removeDegreeTwoVertices();
        }
        catch (AlgorithmException anError) {
            System.out.println("Bad initial tree. "+anError.getMessage());
            System.exit(1);
        }  
                
        // Important: number of steps per edge
        int numStepsPerEdge = 50;

        // Prior distribution on x0. No prior on t0
        DensityCalculator thePrior = new BrownianStateForTopologies.ExponentialEdgePrior(0.5, 1.0, false); 
        //DensityCalculator thePrior = new BrownianStateForTopologies.ExponentialEdgePrior(0.2, 1.0, false); 
        /*
        // Likelihood function
        int numParticles = 1000;
        int numProcs = 1;

        // Chain params
        int nits = 4000, burnits = 0, thin = 1;
        */
        
        // Likelihood function
        int numParticles = 1000;
        int numProcs = 1;

        // Chain params
        int nits = 10000, burnits = 0, thin = 1;
        
        String outputFilename = "/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/5 Taxon Gene Trees/coal_trees_0.6_0.01_0.01_AGT1_MCMCoutv3_2024.txt";
        String analysisFilename = "/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/5 Taxon Gene Trees/coal_trees_0.6_0.01_0.01_AGT1_MCMCanalv3_2024.txt";

        /* End inputs ------------------------------------------------------- */

        /* Make the likelihood function */
        LikelihoodForTopologies theLikelihood = new LikelihoodForTopologies(theData, numStepsPerEdge, numParticles, numProcs);
        /* Make an initial state */
        GlobalState[] initialState = new GlobalState[2];
        BrownianStateForTopologies theBrownianState = null;
        for(int i=0; i<2; i++) {
            BrownianStateForTopologies state = new BrownianStateForTopologies(new TreeAsSplits(initialTree), fixedVariance);
            if (i==0) theBrownianState = state;
            java.util.HashSet<State> set = new java.util.HashSet<State>();
            set.add(state);
            initialState[i] = new GlobalState(set);
            try {
                initialState[i].setLogLikelihood(theLikelihood.logDensity(initialState[i], "brownian_motion_parameters"));
            } catch (AlgorithmException ex) {
                System.out.println("Error evaluating initial likelihood.");
            }
       }   

        /* Make the proposals */

        /* Make wrappers */
        ArrayList<KernelWrapper> kernelWrapperList = new ArrayList<KernelWrapper>();        
//        Kernel x0EdgeProp = new TreeState.EdgeLengthProposal(0.1);
//        Kernel edgeLengthProposal = new BrownianStateForTopologies.X0Proposal(x0EdgeProp);
        Kernel edgeLengthProposal = new BrownianStateForTopologies.LogNormalEdgeProposal(0.1);
        KernelWrapper edgePartialWrapper = new MetropolisHastingsKernelWrapper(edgeLengthProposal, thePrior, theLikelihood, theBrownianState.getName());
        kernelWrapperList.add(edgePartialWrapper);
        
        Kernel randomWalkProposal = new BrownianStateForTopologies.RWEdgeProposal(0.1);
        KernelWrapper randomWalkPartialWrapper = new MetropolisHastingsKernelWrapper(randomWalkProposal, thePrior, theLikelihood, theBrownianState.getName());
        kernelWrapperList.add(randomWalkPartialWrapper);
       
        KernelWrapper myKernelWrapper = new SweepKernelWrapper("Overall proposal", kernelWrapperList);


        /* Run the chain */

        File outputFile = new File(outputFilename);
        Chain myChain = new Chain(initialState[0], initialState[1], myKernelWrapper);
        myChain.run(nits, burnits, thin, outputFile);
        
        //output
        
        System.out.println("Topologies in data set:-");
        PosteriorAnalysis.countTopologiesToScreen(new File(dataFilename));
        System.out.println("");
        System.out.println("Topologies in posterior:-");
        PosteriorAnalysis.countTopologiesToScreen(new File(outputFilename));
        PosteriorAnalysis.outputNumericalParamsForModalTopology(new File(analysisFilename), new File(outputFilename));
        System.out.println();
        
        /*print topology outputs to file
    PrintWriter out;
             try {
                 FileWriter out1= new FileWriter(outputFile, true);
                 BufferedWriter out2 =new BufferedWriter(out1);
                 out = new PrintWriter(out2);

             }
             catch (java.io.IOException anErr) {
                 System.out.println("Warning: output to file "+outputFile.getName()+" failed: "
                         +"writing output to console instead.");
                 out = new PrintWriter(System.out);
        
        String outputTopologyCounts ="/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/4 Taxon Gene Trees/coal_trees_0.10_0.2_TopOuts_2Steps.txt";
        File outputTopologyCountsFile = new File(outputTopologyCounts);
                  
        out.println("Topologies in data set:-");
        PosteriorAnalysis.countTopologiesToFile(new File(dataFilename), outputTopologyCountsFile);
        System.out.println("");
        System.out.println("Topologies in posterior:-");
        PosteriorAnalysis.countTopologiesToFile(new File(outputFilename),outputTopologyCountsFile);
        System.out.println();
        
        out.close();
*/
        
    }
    }
    

