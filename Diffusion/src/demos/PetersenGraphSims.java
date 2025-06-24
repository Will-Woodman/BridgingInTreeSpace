/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package demos;

import diffbase.Simulator;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import simulation.Random;
import treebase.AlgorithmException;
import treebase.Split;
import treebase.Tree;
import treebase.TreeAsSplits;

/**
 *
 * @author tom
 */
public class PetersenGraphSims {
    
    public static TreeAsSplits[] topologies;
    
    private static void makeTopologies() throws AlgorithmException {
        
        String[] t = new String[15];
        t[0] = "((B:1,C:1):1,A:1,(D:1,E:1):1);";
        t[1] = "((B:1,D:1):1,A:1,(C:1,E:1):1);";
        t[2] = "((B:1,E:1):1,A:1,(C:1,D:1):1);";
        t[3] = "((A:1,C:1):1,B:1,(D:1,E:1):1);";
        t[4] = "((A:1,D:1):1,B:1,(C:1,E:1):1);";
        t[5] = "((A:1,E:1):1,B:1,(C:1,D:1):1);";
        t[6] = "((A:1,B:1):1,C:1,(D:1,E:1):1);";
        t[7] = "((A:1,D:1):1,C:1,(B:1,E:1):1);";
        t[8] = "((A:1,E:1):1,C:1,(B:1,D:1):1);";
        t[9] = "((A:1,B:1):1,D:1,(C:1,E:1):1);";
        t[10] = "((A:1,C:1):1,D:1,(B:1,E:1):1);";
        t[11] = "((A:1,E:1):1,D:1,(B:1,C:1):1);";
        t[12] = "((A:1,B:1):1,E:1,(C:1,D:1):1);";
        t[13] = "((A:1,C:1):1,E:1,(B:1,D:1):1);";
        t[14] = "((A:1,D:1):1,E:1,(B:1,C:1):1);";
        
        topologies = new TreeAsSplits[15];
        for (int i=0; i<15; i++) {
            topologies[i] = new TreeAsSplits(new Tree(t[i]));
        }

    }
    
    
    public static double[] getTopologyFrequencies(TreeAsSplits[] theTrees) {
        double[] f = new double[15];
        
        for (int i=0; i<theTrees.length; i++) {
            int ind = -1;
            for (int j=0; j<15; j++) {
                if (theTrees[i].matchesTopology(topologies[j])) {
                    ind = j;
                    break;
                }
            }
            f[ind] += 1.0;
        }
        
        /* Normalize -- divide by number of trees */
        for (int j=0; j<15; j++) {
            f[j] = f[j]/theTrees.length; 
        }
        
        return f;
    }
    
    
    public static void main(String[] args) throws IOException, AlgorithmException  {

        Random.setEngine((int) System.currentTimeMillis());
        
        int numBaseTrees = 1;
        int numDispersionValues = 1;
        double minRootDispersion = 10.0;
        double maxRootDispersion = 10.0;
        int numParticles = 10000;
        int numSteps = 40000;
        
        String filename = "/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/petersen.txt";
        
        /* Make base topologies */
        makeTopologies();
        
        /* Make the file */
        java.io.File theFile = new java.io.File(filename);
        java.io.PrintWriter out = new java.io.PrintWriter(new java.io.BufferedWriter(new java.io.FileWriter(theFile)));
        
        /* Make base tree and splits */
        TreeAsSplits baseTree = new TreeAsSplits(new Tree("((B:1,C:1):1,A:1,(D:1,E:1):1);"));
        Split p = new Split("B,C","A,D,E");
        Split q = new Split("A,B,C","D,E");
        
        /* Main loop */
        for (int i=0; i<numBaseTrees; i++) {
            double theta = Math.PI*0.25*((double)(2*i+1))/((double)numBaseTrees);
            TreeAsSplits x0 = baseTree.clone();
            x0.setSplitLength(p, Math.cos(theta));
            x0.setSplitLength(q, Math.sin(theta));
            //System.out.println(x0.toString());
            
            for (int j=0; j<numDispersionValues; j++) {
                //double sd = minRootDispersion + ((double)j)*(maxRootDispersion-minRootDispersion)/((double)(numDispersionValues-1));
                double t0 = 10;
                
                /* Simulate particles */
                TreeAsSplits[] theTrees = Simulator.sampleRandomWalk(x0, t0, numSteps, numParticles, false, false);
                
//                for (int k=0; k<numParticles; k++) {
//                    System.out.println(theTrees[k].toString(false));
//                }
                
                /* Get frequencies */
                double[] topologyProfile = getTopologyFrequencies(theTrees);
                
                /* Output to file */
                out.print(String.format("%7.7f ", Math.cos(theta)));
                out.print(String.format("%7.7f ", Math.sin(theta)));
                out.print(String.format("%7.7f ", t0));
                for (int k=0; k<15; k++) {
                    out.print(String.format("%7.7f ", topologyProfile[k]));
                }
                out.println();
            }
            
        }
        
        /* Finish off */
        out.close();

    }

    
}
