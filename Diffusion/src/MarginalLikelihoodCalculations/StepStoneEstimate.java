package MarginalLikelihoodCalculations;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
 */

/**
 *
 * @author will
 * Stepping stone estimate where the first distribution (the independence proposal) has been simulated from
 * directly
 */
public class StepStoneEstimate {

    /**
     Main class to calculate repeated stepping stone estimates of marginal likelihoods:
     */
    public static void main(String[] args) throws IOException {
        System.out.println("Step Stone estimate");
        //work out the number of data points
        String PosteriorFileName=args[0];
        String firstFileName = PosteriorFileName+"0.txt";
        
        File firstFile = new File(firstFileName);
        int n = getNumDataPoints(firstFile);

        int numOfFiles=100;//number of beta_k's sampled from

        //calculate and print the estimate
        double out = overallEst(PosteriorFileName,n, numOfFiles);
        System.out.println(out);
    }
    
    /*
    Calculate the stepping stone estimate of marginal likelihood after running the stepping stone sampler
    given a filename (without the number i and extension .txt),
    the number of data points and number of files outputted by the sampler, i.e. the number of points 
    beta_k chosen by the sampler
    */
    private static double overallEst(String posteriorFilename,int n,int numFiles) throws IOException{
        File PostFile;
        Double[] theEsts = new Double[n];//store individual estimates for each data point
        Double theEst=0.0;//store overall estimate
        //now loop over the data points
        for( int l=0; l<n;l++){
        theEsts[l]=0.0;
        Double[] Est= new Double[numFiles];//store estimated values for each beta_k
        double est=0;//sum to get overall estimate for each data point
        for(int j=0;j<numFiles;j++){
            PostFile = new File(posteriorFilename+j+".txt");
            //read in the relevant column from the posterior file for the current data point
            ArrayList<ArrayList<Double>> PostData= readTheFile(PostFile,l);
            ArrayList<Double> PostLogLike = PostData.get(0);
            Est[j]=0.0;
            //Take care of numerical stability by finding largest value
            double eta=Double.NEGATIVE_INFINITY;
            int numSimples=0;
            for(int i=0;i<PostLogLike.size();i++){
            if((!PostLogLike.get(i).isInfinite())&&(PostLogLike.get(i)>eta)){
                eta= PostLogLike.get(i);

            }

            }
            //sum up the values in the average
            for(int i=0;i<PostLogLike.size();i++){    
            if(!PostLogLike.get(i).isInfinite()){
                Est[j]=Est[j]+Math.exp(PostLogLike.get(i)-eta);
                numSimples++;
            }
            }

            Est[j]=Est[j]/PostLogLike.size();
            //Est[j]=Est[j]/numSimples;
            //add back in the largest value
            est=eta+Math.log(Est[j])+est;
            //System.out.println(Math.exp(eta+Math.log(Est[j])));
        }
        theEsts[l]=est;
        //System.out.println(Math.exp(est));
        theEst+=est;
        }
        return theEst;
    }
    
    //read in the information from the files -- only read for each data point at a time because
    //files could be very large for a high number of data points
    private static ArrayList<ArrayList<Double>> readTheFile(File theFile,int k) throws FileNotFoundException, IOException{
        ArrayList<ArrayList<Double>> out = new ArrayList();
        
       
        ArrayList<Double> logLike = new ArrayList();
        
        BufferedReader br = new BufferedReader(new FileReader(theFile));
        String s;
        int i=0;
        while ((s = br.readLine())!=null) {
            if(i ==0){
                i=i+1;
               
            }
            else{
            String[] ss = s.split(" ");
            if(ss.length==1){//alow for case of one data point
            try{

            logLike.add(Double.parseDouble(ss[k]));
            }
            catch(NumberFormatException anErr){
                logLike.add(Double.NEGATIVE_INFINITY);
            }
             
            }
            else if(!(ss[1].equals("rate"))){//don't read in the lines with acceptance rates
            try{

            logLike.add(Double.parseDouble(ss[k]));
            }
            catch(NumberFormatException anErr){
                logLike.add(Double.NEGATIVE_INFINITY);
            }
             
            }
            }
            
            }
        
        
        out.add(logLike);
        
        return(out);
       
    }
    
   /*auxilliary method to get the number of data points from the direct proposals file -- just counts
    the number of columns in the header
    */
   private static int getNumDataPoints(File theFile) throws FileNotFoundException, IOException{
        
        BufferedReader br = new BufferedReader(new FileReader(theFile));
        String s;
        int i=0;
        s = br.readLine();
            String[] ss = s.split(" ");
           
        Double out =(double) ss.length;
        
        return(out.intValue());
        //return(1);
       
    } 
   
   /*
   //This is to loop through estimates when repeatedly estimating MLs to find the variance 
   public static void main(String[] args) throws IOException {
        System.out.println("Step Stone");
        //work out the number of data points
        String firstFileName = "/data/ww24/MarginalLikelihoods/10Taxa20240529/StepStone/Trees1repeat0r0.txt";
        
        File firstFile = new File(firstFileName);
        int n = getNumDataPoints(firstFile);
        //loop over the separate stepping stone runs to see the variance of the sampler
        for(int k=20;k<21;k++){
        int numOfFiles=100;//number of beta_k's sampled
        String PosteriorFileName="/data/ww24/MarginalLikelihoods/10Taxa20240529/StepStone/Trees1repeat"+k+"r";
        //calculate and print the estimate
        double out = overallEst(PosteriorFileName,n, numOfFiles);
        System.out.println(Math.exp(out));
        }
    }
    */
   
    
}
