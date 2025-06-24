/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package demos;

import geodesicgraphics.PathFromArrayOfTrees;
import geodesicgraphics.TreeSpacePathViewer;
import geodesics.Geodesic;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import randomwalks.RandomWalkSimulator;
import simulation.Random;
import treebase.TreeAsSplits;

/**
 *
 * @author ntmwn
 */
public class ViewRandomWalk {
    
        public static void main(String[] args) throws treebase.AlgorithmException, java.io.IOException {

        String str = "((A:1,B:1):1,(C:1,D:1):1,(E:1,F:1):2);";
        double t0 = 4.0;

        if (args.length==1) {
            t0 =  Double.parseDouble(args[0]);
            if (t0<=0.0) throw new NumberFormatException();
        }
        else if (args.length==0) {
            System.out.println("\nPlease provide a dispersion variable and Newick string at the command line. \nDisplaying default.\n");
        }
        else {
            str = args[1];
            t0 =  Double.parseDouble(args[0]);
            if (t0<=0.0) throw new NumberFormatException();
        }
        
        treebase.Tree x = new treebase.Tree(str);
        TreeAsSplits x0 = new TreeAsSplits(x);

        final int k=50;
        TreeAsSplits[] theTrees = new TreeAsSplits[k+1];

        
        // Create an array of trees via random walk
        Random.setEngine((int)System.currentTimeMillis());
        theTrees[0] = x0;
        for (int i=1; i<=k; i++) {
            theTrees[i] = RandomWalkSimulator.sampleRandomWalk(theTrees[i-1], t0/k, 1, false);
        }     
        
        PathFromArrayOfTrees thePath = new PathFromArrayOfTrees(theTrees);
        TreeSpacePathViewer theViewer = new TreeSpacePathViewer(thePath,(k+1));

        // Create and set up window.
        javax.swing.JFrame frame = new javax.swing.JFrame("Tree Space Path Viewer");
        frame.setDefaultCloseOperation(javax.swing.JFrame.EXIT_ON_CLOSE);

 
        frame.setContentPane(theViewer);

        // Set position and size
        frame.setLocation(100,100);

        //Display the window.
        frame.pack();
        frame.setVisible(true);

    }

    
}
