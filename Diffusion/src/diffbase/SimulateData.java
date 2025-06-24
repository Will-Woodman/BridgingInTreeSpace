/*
 * SimulateData.java

    Copyright (C) 2012  Tom M. W. Nye

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

package diffbase;

/**
 * Take a geodesic, evenly distribute points along, and then diffuse out.
 *
 * Comment: it might be better to space out along the geodesic using a normal distribution.
 * This class was basically created to make test datasets for PCA. It probably needs to be deleted!
 */


import geodesics.Geodesic;
import java.io.File;
import java.io.PrintWriter;
import java.io.BufferedWriter;
import java.io.FileWriter;
import treebase.TreeAsSplits;

public class SimulateData {

    public static void simulate(String tA, String tB, int numTrees, int numSteps, double sigma, File theFile) throws treebase.AlgorithmException {
            PrintWriter out;
            try {
                 out = new PrintWriter(new BufferedWriter(new FileWriter(theFile)));
            }
            catch (java.io.IOException anErr) {
                System.out.println("Warning: output to file "+theFile.getName()+" failed: "
                         +"writing output to console instead.");
                out = new PrintWriter(System.out);
            }

            treebase.Tree treeA = new treebase.Tree(tA);
            treeA.removeDegreeTwoVertices();
            TreeAsSplits mapA = new TreeAsSplits(treeA);

            treebase.Tree treeB = new treebase.Tree(tB);
            treeB.removeDegreeTwoVertices();
            TreeAsSplits mapB = new TreeAsSplits(treeB);

            /* Make the geodesic */
            Geodesic g = new Geodesic(mapA, mapB);

            /* Get trees along the geodesic */
            treebase.Tree[] theTrees = g.getCollectionOfTrees(null, numTrees);

            /* Diffuse */
            double var = sigma*sigma;
            for (int i=0; i<numTrees; i++) {
                theTrees[i] = new treebase.TreeWithTopologicalOperations(theTrees[i]);
                Simulator.randomWalkTree((treebase.TreeWithTopologicalOperations) theTrees[i], var, numSteps, false, false);
                out.println(theTrees[i].toString());
            }

            // Flush output
            if (theFile==null) out.flush();
            else out.close();
    }

    public static void main(String[] args) throws treebase.AlgorithmException {

        File theFile = new File("/home/ntmwn/tmp/test_lsqpca.txt");
        int numTrees = 20;
        int numSteps = 50;
        double sigma = 0.001;
        String tA = "((((A:0.2224926152,(B:0.05691056401,(C:0.06023135548,(D:0.09706525234,E:0.03373801815):0.1568416033):0.01468726151):0.03038149411):0.08072250609,(F:0.0006829467136,(G:0.2030891668,H:0.09998635733):0.05055037169):0.04202195657):0.05450168317,I:0.1162593447):0.02166286197,((J:0.2482447313,((((K:0.1051662112,L:0.08885102942):0.06894339356,M:0.03514861643):0.1152782056,N:0.2041272949):0.04201974324,O:0.3101073369):0.008157620067):0.382608338,P:0.006551137986):0.09238393375);";
        String tB = "((((K:0.2342931117,O:0.01629824606):0.07880589533,D:0.2527523843):0.01728966744,(C:0.06558794472,G:0.1292616009):0.0760187108):0.1305966432,(B:0.007208846277,((N:0.04271725053,P:0.0322844578):0.04504326736,((J:0.07568798037,H:0.2558657938):0.1352219474,((E:0.02784406233,(M:0.1049234187,(I:0.07549512256,(A:0.4617796461,L:0.1298683358):0.3370453514):0.2897319592):0.142386648):0.1277288732,F:0.1763035567):0.1660405077):0.06794932405):0.05789672351):0.1901513772);";

        simulate(tA, tB, numTrees, numSteps, sigma, theFile);
    }


}
