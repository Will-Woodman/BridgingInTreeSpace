/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package demos;

import bridge.BridgeWithApproxMVNLike;
import bridge.BridgeState;
import bridge.ForwardStepBridge;
import bridge.StarTree;
import geodesicgraphics.PathFromArrayOfTrees;
import geodesicgraphics.TreeSpacePathViewer;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import randomwalks.RandomWalkSimulator;
import simulation.Random;
import treebase.AlgorithmException;
import treebase.TreeAsSplits;

/**
 *
 * @author ntmwn
 */
public class ViewBrownianBridge {
    
    public static void main(String[] args) throws treebase.AlgorithmException, java.io.IOException {

        String strA = "((A:1,B:1):1,(C:1,D:1):1,(E:1,F:1):2);";
        String strB = "((B:1,C:1):2,(A:1,E:1):1,(D:1,F:1):1);";
        double t0 = 2;

        if (args.length==1) {
            t0 =  Double.parseDouble(args[0]);
            if (t0<=0.0) throw new NumberFormatException();
        }
        else if (args.length==0) {
            System.out.println("\nPlease provide a dispersion variable and file containing two Newick strings at the command line. \nDisplaying default.\n");
        }
        else {
            t0 =  Double.parseDouble(args[0]);
            if (t0<=0.0) throw new NumberFormatException();
            File theFile = new File(args[1]);
            BufferedReader br = new BufferedReader(new FileReader(theFile));
            strA = br.readLine();
            strA = strA.trim();
            strB = br.readLine();
            strB = strB.trim();
            br.close();
        }
        
        treebase.Tree tA = new treebase.Tree(strA);
        treebase.Tree tB = new treebase.Tree(strB);
        TreeAsSplits treeA = new TreeAsSplits(tA);
        TreeAsSplits treeB = new TreeAsSplits(tB);

        final int k=100;
        StarTree.getInstance().setTree(treeA);
        
        // Create an array of trees via bridge
        Random.setEngine((int)System.currentTimeMillis());
        
        ForwardStepBridge template = new BridgeWithApproxMVNLike(treeA, null, k);
        ForwardStepBridge b = null;
        int count = 1;
        while (count>0) {
            try {
                b = template.makeFreshBridge(treeA, treeB);
                b.independenceProposal(t0);
                System.out.println("Bridged in "+count+" attempts.");
                break;
            }
            catch (AlgorithmException ex) {
                count++;
            }
        }
        
        TreeAsSplits[] theTrees = b.getTreePathArray();
        
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
