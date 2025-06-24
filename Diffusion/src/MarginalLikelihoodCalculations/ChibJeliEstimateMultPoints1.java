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
import java.util.Map;

/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
 */

/**
 *
 * @author will
 */
public class ChibJeliEstimateMultPoints1 {

    /**
     Main class to calculate the Chib one block estimate of the marginal likelihood after running the 
     * sampler, in the fixed dispersion case
     */
    public static void main(String[] args) throws IOException {
        System.out.println("Chib");
        //work out the number of data points
        String firstFileName = args[0]+"PropData"+0+".txt";
        Integer StartPoint = Integer.parseInt(args[1]);
        Integer EndPoint = Integer.parseInt(args[2]);
        File firstFile = new File(firstFileName);
        int n = getNumDataPoints(firstFile);
        File outputFile=new File(args[3]);
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
         
        
        
        for(int j=StartPoint;j<EndPoint;j++){
        //loop over the different files to estimate the variance of the estimator
        String PosteriorFileName=args[0]+"PostOut"+j+".txt";
        String PropsFileName=args[0]+"PropData"+j+".txt";
        
        
        
        int propYstars=500;//estimate the log marginal likelihood at every propYstars-th bridge in the posterior
        
        //calculate the estimate
        double theEst = overallEst(PosteriorFileName,PropsFileName,n,propYstars);
        out.println(theEst);
        }
        out.close();  
    }
    
    /*
    Calculate the Chib one block estimate of the marginal likelihood after running sampler. Input the
    name of the sample posterior and proposals file as well as the number of data points and the interval
    for the number of Y^*'s to use in the estimate (see paper)
    */
    public static double overallEst(String posteriorFilename,String PropsFilename,int n,int prop) throws IOException{
        double theEst =0;//to store the overall estimate
        Double[] theEsts = new Double[n];//to store the individual estimates for each data point
        
        File PostFile = new File(posteriorFilename);
        File PropsFile = new File(PropsFilename);
        for( int k=0; k<n;k++){//loop over the data points
        theEsts[k]=0.0;
        ArrayList<ArrayList<Double>> PropData= readTheFile(PropsFile,k);
        ArrayList<ArrayList<Double>> PostData= readTheFile(PostFile,k);
        
        ArrayList<Double> PropsPropDens = PropData.get(1);//proposal density of bridges simmed under the prop
        ArrayList<Double> PropsLogLike = PropData.get(0);//log likelihood of bridges simmed under the prop
        ArrayList<Double> PostPropDens = PostData.get(1);//proposal density of bridges simmed under the posterior
        ArrayList<Double> PostLogLike = PostData.get(0);//log likelihood of bridges simmed under the posterior
        
        //get an estimate of the posterior density for each of the posterior Y's
        double postSize =PostPropDens.size();
        double numYStars=0.0;//how many bridges at which we estimate the posterior density
        int numNans =0;//need to check if this still occurs -- shouldn't do, in which case delete
        for(int i=0 ; i<postSize;i++){
           if(i%prop==0){ //get the estimate for the ith value of Y^*
           numYStars++;
        double PostEst=ChibJeliIndEst(i,PropsPropDens, PropsLogLike, PostPropDens, PostLogLike );
        if(Double.isNaN(PostEst)){
            numNans+=1;//need to check if this still occurs -- shouldn't do, in which case delete
        }
        else{
        theEsts[k]+=Math.exp(PostLogLike.get(i)-PostEst);
        }
        }
        }
        theEsts[k]=theEsts[k]/(numYStars);//average over the estimates
        theEst+=Math.log(theEsts[k]);//add on the estimated log marginal likelihood for the current data point
        //System.out.println(theEsts[k]);
        }
        return theEst;
    }
    
    public static double propRatio( double YPropDens,double YLogLike,double YStarPropDens, double YStarLogLike){
        //acceptance ratio of moving to YStar given Y
        double alpha = YStarLogLike+YPropDens-YLogLike-YStarPropDens;
        //System.out.println(alpha);
        if(YStarPropDens==Double.NEGATIVE_INFINITY){
            return(Double.NEGATIVE_INFINITY);
        }
        return(Math.min(0.0,alpha));
    }
    
    //estimate of the posterior density at a given Y^*, th ith bridge
    public static double ChibJeliIndEst(int i,ArrayList<Double> PropPropDen,ArrayList<Double> PropLogLik,ArrayList<Double> PostPropDen,ArrayList<Double> PostLogLik ){
        double qYStar = PostPropDen.get(i);
        double numerator = 0;
        double M1 = PostPropDen.size();
        for(int j=0; j<M1;j++){
            numerator+=Math.exp(propRatio(PostPropDen.get(j),PostLogLik.get(j),PostPropDen.get(i),PostLogLik.get(i)));
        }
        
        double M2 = PropPropDen.size();
        double denom = 0;
        for(int j=0; j<M2;j++){
            if(PropLogLik.get(j)==Double.NEGATIVE_INFINITY){
            }
            else{
            denom+=Math.exp(propRatio(PostPropDen.get(i),PostLogLik.get(i),PropPropDen.get(j),PropLogLik.get(j)));
            }
        }
        
        return(qYStar-Math.log(M1)+Math.log(M2)+Math.log(numerator)-Math.log(denom));
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
            if(!(ss[1].equals("rate"))){
            propDens.add(Double.parseDouble(ss[2*k]));//add in the IP density
            logLike.add(Double.parseDouble(ss[2*k+1]));//add in the log likelihood
            }
           
            
            }
        }
        
        out.add(logLike);
        out.add(propDens);
        
        return(out);
       
    }
    
      /*auxilliary method to get the number of data points from the direct proposals file -- just counts
    the number of columns in the header and divides by two
    */ 
   private static int getNumDataPoints(File theFile) throws FileNotFoundException, IOException{
        
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
