/*
    AuxilliaryMethods
    Copyright (C) 2025  William M Woodman

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

    Contact the author at:  <w.m.woodman2@ncl.ac.uk>
                           
 */
package MarginalLikelihoods;

import geodesics.Geodesic;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import treebase.AlgorithmError;
import treebase.TreeAsSplits;

/**
 *
 * @author will
 */
public class AuxilliaryMethods {
    
    /*
    Class to get the parameters of the marginal reference distribution on t0 from a posterior file
    Parameters are maximum likelihood estimators of the parameters of a lognormal distn from the t0 sample
    */
    //from a filename
     public static Double[] getRefDistParameters(String posteriorFileName) throws IOException{
    
    File posteriorFile = new File(posteriorFileName);
    return getRefDistParameters(posteriorFile);
}
    //from a file
    public static Double[] getRefDistParameters(File posteriorFile) throws IOException{
    //get the t0 Sample
    ArrayList<Double> t0Posterior= getT0Posterior(posteriorFile);
    //calculate mu
    int n = t0Posterior.size();
    Double mu =0.0;
    
    for(int i=0; i<n;i++){
       mu+= Math.log(t0Posterior.get(i)); 
    }
    mu/=(double) n;
    
    Double sigma =0.0;
    for(int i=0; i<n;i++){
       sigma+= (Math.log(t0Posterior.get(i))-mu)*(Math.log(t0Posterior.get(i))-mu); 
    }
    
    sigma=Math.sqrt(sigma/(double) n);
    Double[] out = new Double[2];
    out[0]=mu;
    out[1]=sigma;
    
   return out;
}
    //read in entries from the 'dispersion' column of a file (if it has one)
   public static ArrayList<Double> getT0Posterior(File posteriorFile) throws IOException{
       //index of the dispersion column
       int Ind =0;
       ArrayList<Double> posterior = new ArrayList();
       //read in the information from the file
        BufferedReader br = new BufferedReader(new FileReader(posteriorFile));
        String s;
        int i=0;
        while ((s = br.readLine())!=null) {
            if(i ==0){
                String[] ss = s.split(" ");
                int j=0;
                boolean Found=false;
                while(!Found){
                    if(ss[j].equals("Dispersion")||ss[j].equals("dispersion")){
                        //System.out.println(ss[j]);
                        Ind = j;//jth column of the file contains t0 values
                        Found = true;
                    }

                    else{
                        //System.out.println(ss[j]);
                        j=j+1;
                    }
                }
                i=i+1;
               
            }
            else{
            String[] ss = s.split(" ");
            if(!(ss[1].equals("rate"))){//don't add when line is an acceptance rate

            posterior.add(Double.parseDouble(ss[Ind]));
             
            }
           
            
            }
        }
        
     
       
       
       return(posterior);
   }
   
   //classes for calculating the Frechet variance of a data set about some specified tree mu
    public static Double getFrechetVariance (TreeAsSplits mu, ArrayList<TreeAsSplits> theTrees) throws AlgorithmError{
         int numTrees= theTrees.size();
         double totalSquDist=0;
         for(int i=0;i<numTrees;i++){
            TreePair theTreePair = new TreePair(theTrees.get(i),mu);
         totalSquDist+= theTreePair.squDist;
         //System.out.println(theTreePair.squDist); - information about the data set but usually no need to print this
        }
         return(totalSquDist/((double)((mu.getNumTaxa()-3)*numTrees)));
         
     }
    
      private static class TreePair{
      TreeAsSplits treeA;
      TreeAsSplits treeB;
      double squDist;
      
      private TreePair(TreeAsSplits tA, TreeAsSplits tB) throws AlgorithmError{
          treeA=tA;
          treeB=tB;
          Geodesic g = new Geodesic(tA,tB);
          double dist = g.getInternalLength();
          squDist=dist*dist;

             
      }
     
  }
}
