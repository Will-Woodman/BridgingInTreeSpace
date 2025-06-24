package MarginalLikelihoodCalculations;



import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
 */

/**
 *
 * @author will
 */
public class ChibJeliEstimateStabilisationTest {

    /**
     * @param args the command line arguments
     */
    
    //Estimate is very slow to run when using all the bridges as a value of Y^*. We try using different numbers of Y^*s to see
    // which gives a good balance of reducing the variance and running quickly
    public static void main(String[] args) throws IOException {
        System.out.println("Chib");

        //the intervals to use: e.g. interval of 2000 means estimate at every 2000th value of the bridge
        Integer[] YStarIntervals = {2000,1000};
        ArrayList<Double[]> allTheEsts = new ArrayList();
        //data file names
        String PosteriorFileName="/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/MarginalLikelihoods/10Taxa/ChibFWM/PostOut";
        String PropsFileName="/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/MarginalLikelihoods/10Taxa/ChibFWM/PropData";
        String outputFilename="/Users/will/Documents/Will PhD/Netbeans Projects/TopInf/MarginalLikelihoods/10Taxa/ChibFWM/StabilizationTest2.txt";
        
         File outputFile = new File(outputFilename);
    PrintWriter out;
             try {
                 FileWriter out1= new FileWriter(outputFile, false);
                 BufferedWriter out2 =new BufferedWriter(out1);
                 out = new PrintWriter(out2);

             }
             catch (java.io.IOException anErr) {
                 System.out.println("Warning: output to file "+outputFile.getName()+" failed: "
                         +"writing output to console instead.");
                 out = new PrintWriter(System.out);
             }
        //loop through the different intervals
        for(Integer i:YStarIntervals){
        //pass to different class to loop through the files for that interval
        Double[] theEsts = ChibEst(PosteriorFileName, PropsFileName, 100, i);
        allTheEsts.add(theEsts);
        System.out.println(theEsts);
        }
        
        //print the estimates to a file
        int numIntervals =YStarIntervals.length;
        int numEsts = allTheEsts.get(0).length;
        for(int i=0;i<numIntervals-1;i++){
            out.print(YStarIntervals[i]+" ");
            }
            out.print(YStarIntervals[numIntervals-1]);
            out.println();
        
        for(int j=0;j<numEsts;j++){
            for(int i=0;i<numIntervals-1;i++){
            out.print(allTheEsts.get(i)[j]+" ");
            }
            out.print(allTheEsts.get(numIntervals-1)[j]);
            out.println();
        }
       
        out.close();          

            
             
    }
    
       public static Double[] ChibEst(String PosteriorFileName,String PropsFileName, Integer numFiles,Integer propYStars) throws IOException {
        Double[] theEsts = new Double[numFiles];
        String ActualPosteriorFileName;
        String ActualPropsFileName;
        for(int j=89;j<numFiles;j++){
        theEsts[j]=0.0;
        ActualPosteriorFileName=PosteriorFileName+j+".txt";
        ActualPropsFileName=PropsFileName+j+".txt";
        File PostFile = new File(ActualPosteriorFileName);
        File PropsFile = new File(ActualPropsFileName);
        theEsts[j]=ChibJeliEstimate.overallEst(ActualPosteriorFileName,ActualPropsFileName,1,propYStars.intValue());
        
  
        
    }
    return(theEsts);
    
    
}
       
}
