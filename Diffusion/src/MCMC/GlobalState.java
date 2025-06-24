/*
    GlobalState
    Copyright (C) 2013 Sarah E. Heaps

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

    Contact the author at:  <sarah.heaps@ncl.ac.uk>
 */

package MCMC;

/**
 * Represents the collection of all States in an MCMC chain
 */

import java.util.ArrayList;
import java.util.Iterator;
import java.util.Set;

public class GlobalState extends State {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;
    
    private double logLikelihood;
    
    public GlobalState(Set<State> initialStates) {
        super("Global state");
        logLikelihood = Double.NaN; // Initially set the log-likelihood to to NaN
        Iterator<State> it = initialStates.iterator();
        while(it.hasNext()) {
            State state = it.next();
            subStates.put(state.getName(), state);
        }
    }
    
    /* Record the likelihood when possible */
    public void setLogLikelihood(double l){logLikelihood = l;}
    public double getLogLikelihood(){return logLikelihood;}
    
    @Override
    public String getHeader() {
        String str = "";
        java.util.Iterator<java.util.Map.Entry<String,State>> it = subStates.entrySet().iterator();
        while (it.hasNext()) {
            java.util.Map.Entry<String,State> pairs = it.next();
            str += pairs.getValue().getHeader()+" ";
        }
        return str.trim();
    }
    
    //delete this in a min
    public String getSubstates() {
        String str = "";
        java.util.Iterator<java.util.Map.Entry<String,State>> it = subStates.entrySet().iterator();
        while (it.hasNext()) {
            java.util.Map.Entry<String,State> pairs = it.next();
            str+=pairs.getKey();
        }
        return str.trim();
    }
    
    @Override
    public String getValueString() {
        String str = "";
        java.util.Iterator<java.util.Map.Entry<String,State>> it = subStates.entrySet().iterator();
        while (it.hasNext()) {
            java.util.Map.Entry<String,State> pairs = it.next();
            str += pairs.getValue().getValueString()+" ";
        }
        return str.trim();
    }
    
    /* For recreating states from printed MCMC output*/
    @Override
    public void setValueFromString(String string) throws treebase.AlgorithmException, OutputStringProcessing.InsufficientInformationException {
        String str = new String(string);
        ArrayList<OutputStringProcessing.RecoverStringInformation> recover = new ArrayList<OutputStringProcessing.RecoverStringInformation>();
        ArrayList<String> recoverString = new ArrayList<String>();
        ArrayList<String> subStateNamesList = new ArrayList<String>();
        if(!str.endsWith(" ")) str += " ";
        int startIndex = 0, endIndex = 0;
        java.util.Iterator<java.util.Map.Entry<String,State>> it = subStates.entrySet().iterator();
        while (it.hasNext()) {
            java.util.Map.Entry<String,State> pairs = it.next();
            int tempIndex;
            for(int j=0; j<pairs.getValue().getLengthStringRepr(); j++) {
                tempIndex = str.indexOf(" ",endIndex);
                endIndex = tempIndex+1;
            }
            String subString = str.substring(startIndex, endIndex-1);
            try {
                pairs.getValue().setValueFromString(subString);
            }
            catch(OutputStringProcessing.InsufficientInformationException anEx) {
                recover.addAll(anEx.recover);
                recoverString.addAll(anEx.partString);
                subStateNamesList.addAll(anEx.subStateNames);
            }
            startIndex = endIndex;
        }
        if(recover.size()>0) throw new OutputStringProcessing.InsufficientInformationException(recover, recoverString, subStateNamesList);
    }
    
    @Override
    public int getLengthStringRepr() {
        int sum = 0;
        java.util.Iterator<java.util.Map.Entry<String,State>> it = subStates.entrySet().iterator();
        while (it.hasNext()) {
            java.util.Map.Entry<String,State> pairs = it.next();
            sum += pairs.getValue().getLengthStringRepr();
        }
        return sum;
    }
    
    /* Compute running mean */
    @Override
    public void updateRunningMean(State runningMean, int runningCount) throws treebase.AlgorithmException {
        java.util.Iterator<java.util.Map.Entry<String,State>> it = subStates.entrySet().iterator();
        while (it.hasNext()) {
            java.util.Map.Entry<String,State> pairs = it.next();
            pairs.getValue().updateRunningMean(runningMean.getSubState(pairs.getKey()), runningCount);
        }
    }
    
}
