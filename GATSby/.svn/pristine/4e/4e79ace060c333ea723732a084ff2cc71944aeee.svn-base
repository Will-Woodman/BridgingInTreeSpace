/*
 * Random
    Copyright (C) 2012  Tom M. W. Nye

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

package simulation;

/**
 * Parallel Random number generation. Provides one global random number generator 
 * per thread with sensibly placed seeds.
 */

import cern.jet.random.tdouble.engine.DoubleMersenneTwister;
import cern.jet.random.tdouble.engine.RandomSeedGenerator;

public class Random {
    
    private static final RandomSeedGenerator seedGenerator = new RandomSeedGenerator();
    
    private static synchronized int getNextSeed() {
        int seed = seedGenerator.nextSeed();
        return seed;
    }
    
    private static final ThreadLocal<DoubleMersenneTwister> randomEngine = new ThreadLocal <DoubleMersenneTwister>() {
        @Override 
        protected DoubleMersenneTwister initialValue() {
            return new DoubleMersenneTwister(getNextSeed());
        }
     };

    /** Returns the value in the current thread's copy of the random engine. If the 
        variable has no value for the current thread, it is first initialized 
        using the initialValue method. */
    public static DoubleMersenneTwister getEngine() {
        return randomEngine.get();
    }
    
    /** Sets the current thread's copy of the random engine to a DoubleMersenneTwister
        with the specified seed. Otherwise this is set on the first invocation of 
        getEngine by the initialValue method.*/
    public static void setEngine(int k) {
        randomEngine.set(new DoubleMersenneTwister(k));
    }
    
    /** Sets the current thread's copy of the random engine to the specified 
        DoubleMersenneTwister. Otherwise this is set on the first invocation of 
        getEngine by the initialValue method.*/
    public static void setEngine(DoubleMersenneTwister eng) {
        randomEngine.set(eng);
    }

}

