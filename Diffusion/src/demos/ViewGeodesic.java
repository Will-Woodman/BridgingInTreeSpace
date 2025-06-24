/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package demos;

import geodesicgraphics.TreeSpacePathViewer;
import geodesics.Geodesic;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import treebase.TreeAsSplits;

/**
 *
 * @author ntmwn
 */
public class ViewGeodesic {
    
    public static void main(String[] args) throws treebase.AlgorithmException, java.io.IOException {

        String strA = "((A:1,B:1):1,(C:1,D:1):1,(E:1,F:1):2);";
        String strB = "((B:1,C:1):2,(A:1,E:1):1,(D:1,F:1):1);";

            File theFile = new File("/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/IndepPropTuning/10Taxa/Endpoints2.txt");
            BufferedReader br = new BufferedReader(new FileReader(theFile));
            strA = br.readLine();
            strA = strA.trim();
            strB = br.readLine();
            strB = strB.trim();
            br.close();
       
        
        treebase.Tree tA = new treebase.Tree(strA);
        treebase.Tree tB = new treebase.Tree(strB);
        TreeAsSplits treeA = new TreeAsSplits(tA);
        TreeAsSplits treeB = new TreeAsSplits(tB);
//        System.out.println(treeA.toString());
//        System.out.println(treeB.toString());
//        System.out.println();

        Geodesic h = new Geodesic(treeA, treeB);
        //System.out.println(h.summary());

        TreeSpacePathViewer theViewer = new TreeSpacePathViewer(h, 100); // 50 intermediate trees

        //Create and set up window.
        javax.swing.JFrame frame = new javax.swing.JFrame("Geodesic Viewer");
        frame.setDefaultCloseOperation(javax.swing.JFrame.EXIT_ON_CLOSE);

 
        frame.setContentPane(theViewer);

        // Set position and size
        frame.setLocation(100,100);
        frame.setSize(100,100);

        //Display the window.
        frame.pack();
        frame.setVisible(true);

    }
    
}
