package MarginalLikelihoodCalculations;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
 */

/**
 *
 * @author will
 */
public class BridgeSamplingEstimateNumStabMultPoints3Spider {

    /**
     Main class to calculate the tunnel sampling estimate of the marginal likelihood after running the 
     * sampler, in the fixed dispersion case
     */

    public static void main(String[] args) throws IOException {
        double initialVal =0.1;//initial value for the iterative estimate
        int numIts  =100;//number of iterations        
        //work out the number of data points
        String firstFileName = "/data/ww24/MarginalLikelihoods/4TaxaOhne/x"+0+"PropData"+0+".txt";
        File firstFile = new File(firstFileName);
        int n = getNumDataPoints(firstFile);
        for(int i=0;i<7;i++){       
        String OutFileName = "/data/ww24/MarginalLikelihoods/4TaxaOhne/x"+i+"results_tunnel.txt";
        File OutFile = new File (OutFileName);
         PrintWriter out;
             try {
                 FileWriter out1= new FileWriter(OutFile, false);
                 BufferedWriter out2 =new BufferedWriter(out1);
                 out = new PrintWriter(out2);

             }
             catch (java.io.IOException anErr) {
                 System.out.println("Warning: output to file "+OutFile.getName()+" failed: "
                         +"writing output to console instead.");
                 out = new PrintWriter(System.out);
            
}
        for(int j=0;j<100;j++){
        //loop over the different files to estimate the variance of the estimator
        String PosteriorFileName="/data/ww24/MarginalLikelihoods/4TaxaOhne/x"+i+"PostOut"+j+".txt";
        String PropsFileName="/data/ww24/MarginalLikelihoods/4TaxaOhne/x"+i+"PropData"+j+".txt";
        //get the estimate
        double theEst = overallEst(PosteriorFileName,PropsFileName,n,initialVal,numIts);
        out.println(Math.exp(theEst));
        }
        out.close();
    }

    }
    
    private static double overallEst(String posteriorFilename,String propsFilename,int n,double initial,int its) throws IOException{
        double theEst =0;//to hold the overall estimate
        Double[] theEsts = new Double[n];// to hold the individual estimates for each data point
        File PostFile = new File(posteriorFilename);
        File PropsFile = new File(propsFilename);
        for( int k=0; k<n;k++){
        theEsts[k]=0.0;
        //read in the values for the current data point
        ArrayList<ArrayList<Double>> PropData= readTheFile(PropsFile,k);
        ArrayList<ArrayList<Double>> PostData= readTheFile(PostFile,k);
        ArrayList<Double> PropsPropDens = PropData.get(1);
        ArrayList<Double> PropsLogLike = PropData.get(0);
        ArrayList<Double> PostPropDens = PostData.get(1);
        ArrayList<Double> PostLogLike = PostData.get(0);
        
        //combine the prop density and log likelihood for the estimate
        ArrayList<Double> PropsValues = normaliseValues(PropsPropDens,PropsLogLike);
        ArrayList<Double> PostValues = normaliseValues(PostPropDens,PostLogLike);
        
        /*take care of numerical stability by subtracting the median from the posterior
        */
        Double lStar = returnMedian(PostValues);
        
        for(int i =0; i<PropsValues.size();i++){
            Double temp = PropsValues.get(i)-lStar;
            PropsValues.set(i, temp);
        }
        
        for(int i =0; i<PostValues.size();i++){
            Double temp = PostValues.get(i)-lStar;
            PostValues.set(i, temp);
        }

        
        
        double out=initial;
        for(int i=0 ; i<its;i++){
        out=bridgeSamplingIt(PropsValues,PostValues,out);//iterate the bridge sampling estimator
        }
        theEsts[k]=Math.log(out)+lStar;
        //System.out.println(Math.exp(theEsts[k]));
        theEst+=theEsts[k];//add the individual estimate
        }
        return theEst;
    }
    
    //calculate one iteration of the bridgeSamplingEst
     public static double bridgeSamplingIt(ArrayList<Double> ProposalVals,ArrayList<Double> PosteriorVals, double current){
        double M1 = PosteriorVals.size();
        double M2 = ProposalVals.size();
        double s1 = M1 / (M1+ M2);
        double s2 = M2/(M1+M2);
        double denom=0;
        double numNans1=0;
        //calculate the denominator from the posterior values
        for(int i=0; i<M1;i++){
        double Div=1/(s1*Math.exp(PosteriorVals.get(i))+s2*current);
        if(!(Double.isNaN(Div))){
        denom+=Div;
        }
        else{
            numNans1+=1;
        }
        }   
        
        //calculate the numerator from the proposal values
        double numerator=0;
        for(int i=0; i<M2;i++){
        double Div=Math.exp(ProposalVals.get(i))/(s1*Math.exp(ProposalVals.get(i))+s2*current);
        if(!(Double.isNaN(Div))){
        numerator+=Div;
        }
        else{ 
            System.out.println("Shouldn't have NaN value here");
        }
        
    }
        return (double)(M1)/((double)(M2))*Math.exp(Math.log(numerator)-Math.log(denom));
    }

     //divide the log likelihood by the prop density:
    public static ArrayList<Double> normaliseValues(ArrayList<Double> PropDensity, ArrayList<Double> LogLike){
        ArrayList<Double> normalised = new ArrayList();
        int numSamples = PropDensity.size();
        for(int i=0; i<numSamples;i++ ){
            normalised.add(LogLike.get(i)-PropDensity.get(i));
            
        }
        
        return normalised;
        
    } 
    
    //find the median of an array list of Doubles
    public static Double returnMedian(ArrayList<Double> vals){
        int numSamples = vals.size();
        Double median=0.0;
        Collections.sort(vals);
        if(numSamples%2==0){
            int Ind1 = numSamples / 2 - 1;
            int Ind2 = numSamples / 2;
            return (vals.get(Ind1) + vals.get(Ind2)) / 2.0;
        }
        else{
            median = vals.get(numSamples/2);
        }
        
        
        return median;
        
    } 
    
        //read in the information from the files -- only read for each data point at a time because
    //files could be very large for a high number of data points
    public static ArrayList<ArrayList<Double>> readTheFile(File theFile,int k) throws FileNotFoundException, IOException{
        ArrayList<ArrayList<Double>> out = new ArrayList();
        
        
        ArrayList<Double> propDens = new ArrayList();
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
            if(!(ss[1].equals("rate"))){//stop once we get to the part of the file with acceptance rates

            propDens.add(Double.parseDouble(ss[2*k]));
            logLike.add(Double.parseDouble(ss[2*k+1]));
             
            }
           
            
            }
        }
        
        out.add(logLike);
        out.add(propDens);
        
        return(out);
       
    }
    
   public static int getNumDataPoints(File theFile) throws FileNotFoundException, IOException{
        
        BufferedReader br = new BufferedReader(new FileReader(theFile));
        String s;
        //introduce t so we don't end up reading last line with the acceptance rates
        int i=0;
        s = br.readLine();
            String[] ss = s.split(" ");
           
        Double out = ss.length/2.0;
        
        return(out.intValue());
       
    } 
    
}
