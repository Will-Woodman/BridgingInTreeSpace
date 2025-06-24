/*
Class to calculate the exact value of the marginal likelihood for the Brownian motion 
kernel when x_0 is the star tree.

 */
package MarginalLikelihoods;

import java.io.File;
import treebase.AlgorithmException;
import treebase.TreeAsSplits;
import treedatasets.TreeAsSplitsDataSet;

public class ExactMLStarTree {
    
   public static void main(String[] args) throws AlgorithmException {
       String dataFileName=args[0];//file containing the data
       //read in the data
       File inFile = new File(dataFileName);
        TreeAsSplitsDataSet theData = null;
        try {
            theData = new TreeAsSplitsDataSet(inFile);
        }
        catch (java.io.IOException anError) {
            System.out.println("Bad input file. "+anError.getMessage());
            System.exit(1);
        }

        Double squEdge;
        Double ML=0.0;
        Double indML=0.0;
        Integer NPrime=theData.theTrees.get(0).getNumTaxa()-3;

        Double correctionFactor=1.0;
        for(int i=0;i<=2*NPrime+1;i++){
            if(i%2==1) correctionFactor = correctionFactor*i;
        }
        
        //t0Est for star tree:
        Double t0=0.0;
        for(TreeAsSplits dataPoint:theData.theTrees){
            t0+=dataPoint.sumSquaredLengths(true);
        }
        t0= t0/theData.theTrees.size()/NPrime;
        
        //calculate the correction factor for each data point
        correctionFactor=NPrime*Math.log(2)-Math.log(correctionFactor);
        
        //sum the log likelihood over each data point:
        for(TreeAsSplits dataPoint:theData.theTrees ){
        indML=correctionFactor;
        indML-= 0.5*NPrime*Math.log(2*Math.PI*t0);
        squEdge=dataPoint.sumSquaredLengths(true);
        indML-= 0.5*1/t0*squEdge;
        ML+=indML;
   }
    
        System.out.println("the estimated log likelihood for the star tree is" + ML);
        
}
   
}
