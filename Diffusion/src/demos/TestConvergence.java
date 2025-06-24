/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package demos;

import MCMC.PosteriorAnalysis;
import static diffbase.Simulator.sampleRandomWalk;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import simulation.Random;
import treebase.AlgorithmException;
import treebase.Tree;
import treebase.TreeAsSplits;

/**
 *
 * @author ntmwn
 */
public class TestConvergence {
    
    public static void main(String[] args) throws AlgorithmException, IOException {

        Random.setEngine((int)System.currentTimeMillis());

        int numSamples=10000;
        int[] numStepsPerEdge = {25, 50, 100, 200, 500, 1000, 2000};
        double[] sig = {0.025, 0.05, 0.1, 0.15, 200};
        
        double[][] pr = new double[numStepsPerEdge.length][sig.length];
        
       treebase.Tree x0 = new treebase.Tree("((B:1,C:1):1,A:1,(D:1,E:1):1);");
       
       /*treebase.Tree x0 = new treebase.Tree("(1:0.5,2:0.5,3:0.5);");
       */
        x0.removeDegreeTwoVertices();
        
        for (int i=0; i<numStepsPerEdge.length; i++) {
            for (int j=0; j<sig.length; j++) {
                Tree[] trees = sampleRandomWalk(x0, sig[j]*sig[j], numStepsPerEdge[i], false, true, numSamples);
                ArrayList<TreeAsSplits> theTrees = new ArrayList();
                for (int k=0; k<numSamples; k++) {
                    TreeAsSplits t = new TreeAsSplits(trees[k]);
                    theTrees.add(t);
                }
                HashMap<String,Integer> results = PosteriorAnalysis.countTopologies(theTrees);
                String modalTop = results.keySet().iterator().next();
                pr[i][j] = ((double)results.get(modalTop).intValue())/((double)numSamples);
                
                System.out.print(String.format("%7.7f ", pr[i][j]));
            }
            System.out.println("");
        }
        
        
        
    }
        
        


    
}
