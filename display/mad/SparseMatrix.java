
/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Adam Roberts.
 *  Code for calculating and storing the sum and max probabilities
 *  was written by Robert Bradley.
 */


package mad;

//import it.unimi.dsi.fastutil.ints.*;
import java.util.*;
import java.lang.Math.*;

public class SparseMatrix{
	
    private int length, width;
    private Vector<Double> values;
    private Map<Integer, Integer> indexMap;

    private Double[] sumProbI;
    private Double[][] maxProbI;

    private Double[] sumProbJ;
    private Double[][] maxProbJ;

    public SparseMatrix(int length, int width){
	this.length = length; // length == length of first sequence + 1
	this.width = width;   // width == length of second sequence + 1

	indexMap = new HashMap<Integer, Integer>();
	values = new Vector<Double>(100, 10);

	sumProbI = new Double[length-1];
	maxProbI = new Double[length-1][2];

	sumProbJ = new Double[width-1];
	maxProbJ = new Double[width-1][2];

    }
    
    public void put(int i, int j, double e){
	int index = j*length+i;

	if (indexMap.containsKey(index)){
	    System.err.println("SANITY CHECK FAILED:  Matrix index repeated (i = " + i + ", j = " + j + ", index = " + index + ")");
	    System.exit(-1);
	}
	
	// check that the probability is valid
	e = SparseMatrix.saneProb (e);

	indexMap.put(index, values.size());
	values.addElement(e);

	// Note that max is over all possible positions in the other sequence, including gapped.
	// The below check over single indices ensures that the possibility of gaps
	// in the other sequence are properly taken into account when finding the max (and alternate max).

	// first sequence
	if (i < length-1) { // i == (length - 1) corresponds to first sequence gapped
                            //   (NB sequence position indices are 0-based)

	    // if max prob undefined, store new value
	    if (maxProbI[i][0] == null) {
		maxProbI[i][0] = (Double) e;
	    }

	    // if new max, then store new value and update alternate max
	    else if (maxProbI[i][0].doubleValue() < e) {
		maxProbI[i][1] = maxProbI[i][0];
		maxProbI[i][0] = (Double) e;
	    }

	    // if new alternate max b/c alternate max is undefined
	    else if (maxProbI[i][1] == null) {
		maxProbI[i][1] = (Double) e;
	    }

	    // if new alternate max
	    else if (maxProbI[i][1].doubleValue() < e) {
		maxProbI[i][1] = (Double) e;
	    }

	}

	// second sequence
	if (j < width-1) { // j == (width - 1) corresponds to second sequence gapped

	    // if max prob undefined, store new value
	    // and initialize alternate max
	    if (maxProbJ[j][0] == null) {
		maxProbJ[j][0] = (Double) e;
	    }

	    // if new max, then store new value and update alternate max
	    else if (maxProbJ[j][0].doubleValue() < e) {
		maxProbJ[j][1] = maxProbJ[j][0];
		maxProbJ[j][0] = (Double) e;
	    }

	    // if new alternate max b/c alternate max is undefined
	    else if (maxProbJ[j][1] == null) {
		maxProbJ[j][1] = (Double) e;
	    }

	    // if new alternate max
	    else if (maxProbJ[j][1].doubleValue() < e) {
		maxProbJ[j][1] = (Double) e;
	    }

	}

	// sum is only over possible aligned character pairs
	// (if gaps are included then all sums be 1. b/c they're normalized probability distributions)
	if ((i < length-1) && (j < width-1)) {

	    if (sumProbI[i] == null)
		sumProbI[i] = (Double) e;
	    else
		sumProbI[i] += (Double) e;

	    if (sumProbJ[j] == null)
		sumProbJ[j] = (Double) e;
	    else
		sumProbJ[j] += (Double) e;

	}

    }
    
    public double get(int i, int j){
	int index = j*length+i;
	if (!indexMap.containsKey(index))
	    return 0.0;
	return SparseMatrix.saneProb (values.get(indexMap.get(index)));
    }
    
    /// Get sum_pos2 P((seq1, pos1) ~ (seq2, pos2)).
    /*
     * Cover the case of no entry for a match probability > 0.01.
     */
    public double getSumProbI(int pos){
	if (sumProbI[pos] == null) // cover case of no input data (no match prob above FSA's threshold)
	    return 0.01;
	return SparseMatrix.saneProb (sumProbI[pos]);
    }

    /// Get sum_pos1 P((seq1, pos1) ~ (seq2, pos2)).
    /*
     * Cover the case of no entry for a match probability > 0.01.
     */
    public double getSumProbJ(int pos){
	if (sumProbJ[pos] == null) // cover case of no input data (no match prob above FSA's threshold)
	    return 0.01;
	return SparseMatrix.saneProb (sumProbJ[pos]);
    }

    /// Get max_pos2 P((seq1, pos1) ~ (seq2, pos2)).
    public double getMaxProbI(int pos){
	// sanity check
	if (maxProbI[pos][0] == null) {
	    System.err.println ("A max prob is null; this should never happen!");
	    System.exit(1);
	}
	return SparseMatrix.saneProb (maxProbI[pos][0]);
    }

    /// Get altmax_pos2 P((seq1, pos1) ~ (seq2, pos2)).
    /*
     * \return second-best probability, or 0.01 if undefined
     * Minimum value of 0.01 is due to the inherent coarse-graining
     * of the probability reported by FSA.
     */
    public double getAltMaxProbI(int pos1, int pos2){

	// sanity check
	if (maxProbI[pos1][0] == null) {
	    System.err.println ("A max prob is null; this should never happen!");
	    System.exit(1);
	}

	// get probability to determine best alternate
	double prob = get(pos1,pos2);

	// if prob is the max, then return the second-best alternate
	if (Math.abs (prob - maxProbI[pos1][0].doubleValue()) < 0.0001) {
	    if (maxProbI[pos1][1] != null)
		return SparseMatrix.saneProb (maxProbI[pos1][1]);
	    else
		return 0.01;
	}
	// else return the max
	else
	    return SparseMatrix.saneProb (maxProbI[pos1][0]);

    }

    /// Get altmax_pos1 P((seq1, pos1) ~ (seq2, pos2)).
    /*
     * \return second-best probability, or 0.01 if undefined
     * Minimum value of 0.01 is due to the inherent coarse-graining
     * of the probability reported by FSA.
     */
    public double getAltMaxProbJ(int pos1, int pos2){

	// sanity check
	if (maxProbJ[pos2][0] == null) {
	    System.err.println ("A max prob is null; this should never happen!");
	    System.exit(1);
	}

	// get probability to determine best alternate
	double prob = get(pos1,pos2);

	// if prob is the max, then return the second-best alternate
	if (Math.abs (prob - maxProbJ[pos2][0].doubleValue()) < 0.0001) {
	    if (maxProbJ[pos2][1] != null)
		return SparseMatrix.saneProb (maxProbJ[pos2][1]);
	    else
		return 0.01;
	}
	// else return the max
	else
	    return SparseMatrix.saneProb (maxProbJ[pos2][0]);

    }

    /// Get max_pos1 P((seq1, pos1) ~ (seq2, pos2)).
    public double getMaxProbJ(int pos){

	// sanity check
	if (maxProbJ[pos][0] == null) {
	    System.err.println ("A max prob is null; this should never happen!");
	    System.exit(1);
	}

	return SparseMatrix.saneProb (maxProbJ[pos][0]);
    }

    /// Bounds-checking for a valid probability.
    static double saneProb (double p) {
	if (p < -0.0001 || p > 1.0001) {
	    System.err.println("Error: probability " + p + " is out of bounds.");
	    System.exit(1);
	}
	p = (p < -0.0001) ? 0. : p;
	p = (p > 1.0001) ? 1.0 : p;
	return p;	
    }

    static double saneProb (Double p) {
	return saneProb (p.doubleValue());
    }

}