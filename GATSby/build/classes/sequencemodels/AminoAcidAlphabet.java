/**
    AminoAcidAlphabet
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
 * Static AA alphabet -- singleton class
 */

import java.util.ArrayList;
import java.util.HashMap;

public class AminoAcidAlphabet extends Alphabet {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;

    private static final String[] data = {"A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"};

    private static final String[] codons = {"TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG",
                                            "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG",
                                            "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG",
                                            "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG",
                                            "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG",
                                            "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG",
                                            "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG",
                                            "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG"};
    private static final String[] codonTargets = {"F", "F", "L", "L", "S", "S", "S", "S",
                                                  "Y", "Y", "-", "-", "C", "C", "-", "W",
                                                  "L", "L", "L", "L", "P", "P", "P", "P",
                                                  "H", "H", "Q", "Q", "R", "R", "R", "R",
                                                  "I", "I", "I", "M", "T", "T", "T", "T",
                                                  "N", "N", "K", "K", "S", "S", "R", "R",
                                                  "V", "V", "V", "V", "A", "A", "A", "A",
                                                  "D", "D", "E", "E", "G", "G", "G", "G",}; // "-" used for stop!
    
    private static final java.util.HashMap translationTable;
    static {
        translationTable = new java.util.HashMap();
        for (int i=0; i<64; i++) {
            translationTable.put(codons[i], codonTargets[i]);
        }
    }
    
    private AminoAcidAlphabet() {
        // Exists only to defeat instantiation.
        
        // Set up data
        letters = new ArrayList();
        for (int i=0; i<20; i++) {
            letters.add(data[i]);
        }

    }
    
    /**
     * AminoAcidAlphabetSingletonHolder is loaded on the first execution of
     * AminoAcidAlphabet.getInstance() or the first access to 
     * AminoAcidAlphabetSingletonHolder.INSTANCE, not before.
     */
    private static class AminoAcidAlphabetSingletonHolder {

        public static final AminoAcidAlphabet INSTANCE = new AminoAcidAlphabet();
    }

    public static AminoAcidAlphabet getInstance() {
        return AminoAcidAlphabetSingletonHolder.INSTANCE;
    }


    public static String translateDNA(String s) throws treebase.AlgorithmException{
        String res = new String();
        for (int i=0; i<s.length()-2; i=i+3) {
            Object o = translationTable.get(s.substring(i, i+3));
            if (o==null) throw new treebase.AlgorithmException("Bad codon when translating DNA to amino acid");
            res += (String) o;
        }
        return res;
    }
    
    // Guarantees singleton-ness of object upon serialization and deserialization
    public Object readResolve() {
        return AminoAcidAlphabetSingletonHolder.INSTANCE;
    }


}

