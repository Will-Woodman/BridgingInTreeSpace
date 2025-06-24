/**
    Alphabet
    Abstract class representing an alphabet of symbols eg. DNA or amino acids

    Copyright (C) 2011  Tom M. W. Nye

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
 * Abstract class representing an alphabet of symbols eg. DNA or amino acids
 * Provides translation between characters and integers.
 *
 * translateLettersToIndex method will need overriding in codon alphabet because
 * it uses triplets
 *
 * Gap is coded as -1, missing as -2
 */

import java.util.ArrayList;
import treebase.AlgorithmException;

public abstract class Alphabet implements java.io.Serializable {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;

    protected ArrayList letters; /** Ensure letters are UPPER CASE */
    public final String gap = "-";
    public final String missing = "?";
    public final static int MISSING_CODE=-2, GAP_CODE=-1;

    public int size() {
        return letters.size();
    }

    public int translateLetterToIndex(String s) throws AlgorithmException {
        if (s.equals(gap)) return GAP_CODE;
        if (s.equals(missing)) return MISSING_CODE;
        int k = letters.indexOf(s.toUpperCase());
        if (k<0) {
            throw new AlgorithmException("Bad translation call to alphabet: "+s);
        }
        return k;
    }

    public String translateIndexToLetter(int k) throws AlgorithmException {
        if (k==GAP_CODE) return gap;
        if (k==MISSING_CODE) return missing;
        if ((k<0)||(k>=letters.size())) {
            throw new AlgorithmException("Bad translation call to alphabet");
        }
        return (String) letters.get(k);
    }

    public int[] translateLettersToIndex(String s) throws AlgorithmException {
        int[] res = new int[s.length()];
        for (int i=0; i<s.length(); i++) {
           res[i] = translateLetterToIndex(s.substring(i, i+1));
        }
        return res;
    }

    public String translateIndicesToLetters(int ind[]) throws AlgorithmException {
        String s = new String();
        for (int i=0; i<ind.length; i++) {
            s +=translateIndexToLetter(ind[i]);
        }
        return s;
    }

    public String getAllLetters() {
        String res = new String();
        for (int i=0; i<letters.size(); i++) {
            res += (String) letters.get(i);
        }
        return res;
    }

}
