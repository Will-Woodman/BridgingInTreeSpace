/**
    DNAAlphabet

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
 * Static DNA alphabet -- singleton class
 */

import java.util.ArrayList;

public class DNAAlphabet extends Alphabet {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;
    
    private static final String[] data = {"A", "G", "C", "T"};
    
    protected DNAAlphabet() {
        // Exists only to defeat instantiation.

        // Set up data
        letters = new ArrayList();
        for (int i=0; i<4; i++) {
            letters.add(data[i]);
        }
    }

    /**
     * DNAAlphabetSingletonHolder is loaded on the first execution of
     * DNAAlphabet.getInstance() or the first access to 
     * DNAAlphabetSingletonHolder.INSTANCE, not before.
     */
    private static class DNAAlphabetSingletonHolder {

        public static final DNAAlphabet INSTANCE = new DNAAlphabet();
    }

    public static DNAAlphabet getInstance() {
        return DNAAlphabetSingletonHolder.INSTANCE;
    }
    
    // Guarantees singleton-ness of object upon serialization and deserialization
    public Object readResolve() {
        return DNAAlphabetSingletonHolder.INSTANCE;
    }


}


