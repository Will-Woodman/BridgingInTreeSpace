/**
    DayhoffSixAlphabet

    Copyright (C) 2013  Sarah E. Heaps

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

package sequencemodels;

/**
 * Static Dayhoff-six alphabet -- singleton class
 * 
 * The Dayhoff-six recoding defines the following six groups of amino acids corresponding to the PAM matrix: 
 * 1, cysteine (C)
 * 2, alanine (A), serine (S), threonine (T), proline (P), glycine (G)
 * 3, asparagine (N), aspartic acid (D), glutamic acid (E), glutamine (Q) 
 * 4, histidine (H), arginine (R), lysine (K)
 * 5, methionine (M), isoleucine (I), leucine (L), valine (V)
 * 6, phenylalanine (F), tyrosine (Y), tryptophan (W)
 * 
 */


import java.util.ArrayList;


public class DayhoffSixAlphabet extends Alphabet {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;
    
    private static final String[] data = {"1", "2", "3", "4", "5", "6"};
    
    private static final String[] aminoAcids = {"A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", 
                                                "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"};
            
    private static final String[] dayhoffTargets = {"2", "4", "3", "3", "1", "3", "3", "2", "4", "5",
                                                   "5", "4", "5", "6", "2", "2", "2", "6", "6", "5"};
    
    private static final java.util.HashMap translationTable;
    static {
        translationTable = new java.util.HashMap();
        for (int i=0; i<20; i++) {
            translationTable.put(aminoAcids[i], dayhoffTargets[i]);
        }
    }
    
    protected DayhoffSixAlphabet() {
        // Exists only to defeat instantiation.

        // Set up data
        letters = new ArrayList();
        for (int i=0; i<6; i++) {
            letters.add(data[i]);
        }
    }

    /**
     * DayhoffSixAlphabetSingletonHolder is loaded on the first execution of
     * DayhoffSixAlphabet.getInstance() or the first access to 
     * DayhoffSixAlphabetSingletonHolder.INSTANCE, not before.
     */
    private static class DayhoffSixAlphabetSingletonHolder {

        public static final DayhoffSixAlphabet INSTANCE = new DayhoffSixAlphabet();
    }

    public static DayhoffSixAlphabet getInstance() {
        return DayhoffSixAlphabetSingletonHolder.INSTANCE;
    }
    
    // Guarantees singleton-ness of object upon serialization and deserialization
    public Object readResolve() {
        return DayhoffSixAlphabetSingletonHolder.INSTANCE;
    }
    
    public static String translateAminoAcids(String s) throws treebase.AlgorithmException{
        String res = new String();
        for (int i=0; i<s.length(); i++) {
            Object o = translationTable.get(s.substring(i, i+1));
            if (o==null) throw new treebase.AlgorithmException("Bad amino acid when translating amino acid to Dayhoff-six recoding");
            res += (String) o;
        }
        return res;
    }
    
}
