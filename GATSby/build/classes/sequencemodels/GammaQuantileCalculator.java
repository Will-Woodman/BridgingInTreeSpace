/**
    GammaQuantileCalculator
    Likelihood and simulation methods for classic gamma mixture model of rate heterogeneity

    Copyright (C) 2012  Tom M. W. Nye & Sarah E. Heaps

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

package sequencemodels;

/**
 * Calculate quantiles of the gamma distribution.
 * Singleton, immutable class.
 */

import java.util.*;

public final class GammaQuantileCalculator implements java.io.Serializable {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;

    /* Instance variables */
    private final ArrayList<double[][]> lookUpTables;
    private final double[] UPPER_LIMITS = {Double.NaN,0.10,0.13,0.15,0.17,0.18,0.19,0.20};
    /* These are the alpha values below which the quantile function breaks for numCategories=1,...,8.*/

    /**
     * GammaQuantileCalculatorSingletonHolder is loaded on the first execution of
     * GammaQuantileCalculator.getInstance() or the first access to 
     * GammaQuantileCalculatorSingletonHolder.INSTANCE, not before.
     */
    private static class GammaQuantileCalculatorSingletonHolder {

        public static final GammaQuantileCalculator INSTANCE = new GammaQuantileCalculator();
    }

    public static GammaQuantileCalculator getInstance() {
        return GammaQuantileCalculatorSingletonHolder.INSTANCE;
    }

    private GammaQuantileCalculator() {
        /* Set up look-up tables -- read in from data folders */
        lookUpTables = new ArrayList<double[][]>();
        //System.out.println(UPPER_LIMITS.length);
        for(int i=0;i<UPPER_LIMITS.length;i++){
            try{
                lookUpTables.add(this.readLookUpTable(i+1));
            }
            catch(java.io.IOException anEx) {
                System.out.println("Problem reading look-up-table for "+(i+1)+" categories. Exiting.");
                System.exit(1);
            }
        }
    }
    
    // Guarantees singleton-ness of object upon serialization and deserialization
    public Object readResolve() {
        return GammaQuantileCalculatorSingletonHolder.INSTANCE;
    }

    private double[][] readLookUpTable(int numCategories) throws java.io.IOException {
        if(numCategories==1){
            return null;
        }
        String resourcesPath = "data/qGammaLookUpTableFor"+numCategories+"Cats.dat";
        java.io.InputStream stream = GammaQuantileCalculator.class.getResourceAsStream(resourcesPath);
        //java.io.File theFile = new java.io.File("src/data/qGammaLookUpTableFor"+numCategories+"Cats.dat");
        String str;
        java.util.StringTokenizer st;
        double[][] tempTable;

        ArrayList<String> data = new ArrayList<String>();
        //java.io.BufferedReader br = new java.io.BufferedReader(new java.io.FileReader(theFile));
        java.io.BufferedReader br = new java.io.BufferedReader(new java.io.InputStreamReader(stream));
        // Use header to establish number of columns
        str = br.readLine();
        st = new StringTokenizer(str);
        int numCols=st.countTokens();
        // Read data into ArrayList
        str = br.readLine();
        while (str!=null){
            st = new StringTokenizer(str);
            while(st.hasMoreTokens()) data.add(st.nextToken());
            str = br.readLine();
        }
        br.close();
        int numRows=data.size()/numCols;

        //Put the data (currently as strings) into the look-up table
        //and convert to doubles
        tempTable = new double[numRows][numCols];
        for(int i=0;i<numRows;i++){
            for(int j=0;j<numCols;j++){
                tempTable[i][j]=Double.parseDouble(data.remove(0));
            }
        }
        return tempTable;
    }

    public void setRates(double[] rates, double alpha) {
        if(alpha < UPPER_LIMITS[rates.length-1]) {
            // Linearly interpolate
            double increment = lookUpTables.get(rates.length-1)[lookUpTables.get(rates.length-1).length-1][0]/(lookUpTables.get(rates.length-1).length-1.0);
            int line = (int)Math.ceil(alpha/increment);
            double prop = (alpha - lookUpTables.get(rates.length-1)[line-1][0])/increment;
            for(int i=0;i<rates.length;i++) {
                rates[i]= lookUpTables.get(rates.length-1)[line-1][i+1]*(1.0-prop) + lookUpTables.get(rates.length-1)[line][i+1]*prop;
            }
            //System.out.println("Warning: had to use look-up table for alpha = "+alpha);
        }
        else {
            // our original method
            double p, delta=1.0/rates.length;
            double theta=1.0/alpha;
            for (int i=0; i<rates.length; i++) {
                p = (i+0.5)*delta;
                rates[i] = simulation.GammaDistribution.quantile(p, alpha, theta);
            }
            // mrBayes weird method
            /*int K = rates.length;
            for(int i=0; i<K-1; i++) rates[i]=simulation.GammaDistribution.quantile((i+1.0)/K, alpha, 1.0/alpha);
            for(int i=0; i<K-1; i++) rates[i]=simulation.GammaDistribution.cdf(rates[i]*alpha, alpha+1.0, 1.0);
            rates[K-1] = 1.0;
            for(int i=K-1; i>0; i--) {
                rates[i]-= rates[i-1];
                rates[i] *= K;
            }
            rates[0] *= K;*/
        }

    }


    
}
