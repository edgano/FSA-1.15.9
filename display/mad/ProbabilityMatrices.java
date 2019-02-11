
/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Adam Roberts.
 *  Code for calculating and storing the sum and max probabilities
 *  was written by Robert Bradley.
 */


package mad;

import java.util.*;

public class ProbabilityMatrices{
	int numSeqs;
	int[] seqLengths;
	int[][] index;
	SparseMatrix[] matrices; 

	public ProbabilityMatrices(int[] seqLengths){
		this.seqLengths = seqLengths;
		this.numSeqs = seqLengths.length;
		this.matrices = new SparseMatrix[numSeqs*(numSeqs-1)/2];
		buildIndex();
	}
	
	private void buildIndex(){
		int x = 0;
		this.index = new int[numSeqs][numSeqs];
		for (int i = 0; i < numSeqs; i++)
			for (int j = i+1; j < numSeqs; j++){
				this.index[i][j] = x;
				this.index[j][i] = x++;
			}
	}
	
	public void addElement(int seq1, int seq2, int i1, int i2, double prob){
		
		//Handle gaps and index change from input file
		if (i1 < 0)
			i1 = seqLengths[seq1];
		if (i2 < 0)
			i2 = seqLengths[seq2];
			
		if (seq1 > seq2){
			int temp = i1;
			i1 = i2;
			i2 = temp;
		}
				
		if (!matrixExists(seq1, seq2))
			addMatrix(seq1, seq2);
		
		matrices[getIndex(seq1,seq2)].put(i1, i2, prob);
	}
	
	public double getElement(int seq1, int seq2, int i1, int i2){
		
		if (!matrixExists(seq1, seq2))
			return 0.0;
			
		// Handle gaps
		if (i1 < 0)
			i1 = seqLengths[seq1];
		if (i2 < 0)
			i2 = seqLengths[seq2];
		
		if (seq1 > seq2){
			int temp = i1;
			i1 = i2;
			i2 = temp;
		}	
		return matrices[getIndex(seq1,seq2)].get(i1, i2);
	}

    /// Get sum_pos2 P((seq1, pos1) ~ (seq2, pos2)).
    public double getSumProb(int seq1, int seq2, int pos1){
	
	if (!matrixExists(seq1, seq2)) {
	    System.err.println ("Probability matrix for sequences " + seq1 + " and " + seq2 + " does not exist.");
	    System.exit(1);
	}
	
	// no gaps
	if (pos1 < 0) {
	    System.err.println ("Sum probability is undefined for gaps; offending sequences are " + seq1 + " and " + seq2 + ".");
	    System.exit(1);
	}

	if (seq1 < seq2)
	    return matrices[getIndex(seq1,seq2)].getSumProbI(pos1);
	else
	    return matrices[getIndex(seq2,seq1)].getSumProbJ(pos1);

    }

    /// Get max_pos2 P((seq1, pos1) ~ (seq2, pos2)).
    public double getMaxProb(int seq1, int seq2, int pos1){
	
	if (!matrixExists(seq1, seq2)) {
	    System.err.println ("Probability matrix for sequences " + seq1 + " and " + seq2 + " does not exist.");
	    System.exit(1);
	}
	
	// no gaps
	if (pos1 < 0) {
	    System.err.println ("Max probability is undefined for gaps; offending sequences are " + seq1 + " and " + seq2 + ".");
	    System.exit(1);
	}

	if (seq1 < seq2)
	    return matrices[getIndex(seq1,seq2)].getMaxProbI(pos1);
	else
	    return matrices[getIndex(seq2,seq1)].getMaxProbJ(pos1);

    }

    /// Get altmax_pos2 P((seq1, pos1) ~ (seq2, pos2)).
    public double getAltMaxProb(int seq1, int seq2, int pos1, int pos2){
	
	if (!matrixExists(seq1, seq2)) {
	    System.err.println ("Probability matrix for sequences " + seq1 + " and " + seq2 + " does not exist.");
	    System.exit(1);
	}
	
	// handle gaps
	if (pos1 < 0) {
	    System.err.println ("Alt max probability is undefined for gaps in first sequence; offending sequences are " + seq1 + " and " + seq2 + ".");
	    System.exit(1);
	}
	else if (pos2 < 0) {
	    pos2 = seqLengths[seq2]; // convert to indexing used by SparseMatrix
	}

	if (seq1 < seq2) {
	    return matrices[getIndex(seq1,seq2)].getAltMaxProbI(pos1,pos2);
	}
	else {
	    return matrices[getIndex(seq2,seq1)].getAltMaxProbJ(pos2,pos1);
	}

    }


	
	private int getIndex(int seq1, int seq2){
		return index[seq1][seq2];
	}
	
	private void addMatrix(int seq1, int seq2){
	
		if (seq1 > seq2){
			int temp = seq1;
			seq1 = seq2;
			seq2 = temp;
		}
		matrices[getIndex(seq1, seq2)] = new SparseMatrix(seqLengths[seq1]+1, seqLengths[seq2]+1);
	}
	
	
	public boolean matrixExists(int seq1, int seq2){
		if (seq1 == seq2)
			return false;
		return matrices[getIndex(seq1, seq2)] != null;
	}
	

	
}