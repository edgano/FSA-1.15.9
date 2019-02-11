
/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Adam Roberts.
 *  Code for different accuracy measures was written by Robert Bradley.
 */


package mad;

import java.util.*;
//import it.unimi.dsi.fastutil.objects.*;

public class AlignDAG{
    
    private Node[] V;  // Vertices in our order.
    private Node[] oV; // Vertices in AMAP's order.
    private Node[] startV; // Vertices that begin each sequence
    private HashMap<Node, Node> mergeDict;

    private double acc;
    private double sps;
    private double ppv;
    private double cert;
    private double cons;

    private double accDenom;
    private double spsDenom;
    private double ppvDenom;
    private double certDenom;
    private double consDenom;

    private double deltaAccNorm;
    private double deltaSPSNorm;
    private double deltaPPVNorm;
    private double deltaCertNorm;
    private double deltaConsNorm;

    private double gapFactorCurrent;  /// the current gap factor

    private int colorScheme;
    private int mergeOrder;

    private Set<Node> activeNodes;
    private int numSeqs;

    public AlignDAG (List<String> keys, String[] sequences, HashMap<String, Integer> initDAG, ProbabilityMatrices pm) {
        this.mergeDict = new HashMap<Node, Node>();
        this.numSeqs = sequences.length;

	this.colorScheme = 0;
	this.mergeOrder = 0;

        buildDagFromInit(keys, sequences, initDAG, pm);  
    }
    
    public AlignDAG (List<String> keys, String[] sequences, String[] alignedSeqs, ProbabilityMatrices pm){
    	this.mergeDict = new HashMap<Node, Node>();
    	this.numSeqs = sequences.length;

	this.colorScheme = 0;
	this.mergeOrder = 0;

    	buildDagFromAlign(keys, sequences, alignedSeqs, pm);
    }
    
    /// Update color scheme for all active nodes.
    /*
     * Forces re-coloring of all columns according to the new scheme.
     */
    public void setColorScheme(int scheme) {
	this.colorScheme = scheme;
	Iterator iter = activeNodes.iterator();
	while (iter.hasNext()) {
	    ((Node)iter.next()).setColorScheme(scheme);
	}
    }

    private void buildDagFromAlign(List<String> keys, String[] sequences, String[] alignedSeqs, ProbabilityMatrices pm){

	acc = 0;
	sps = 0;
	ppv = 0;
	cert = 0;
	cons = 0;

	accDenom = 0;
	spsDenom = 0;
	ppvDenom = 0;
	certDenom = 0;
	consDenom = 0;

	deltaAccNorm = 0;
	deltaSPSNorm = 0;
	deltaPPVNorm = 0;
	deltaCertNorm = 0;
	deltaConsNorm = 0;

	gapFactorCurrent = -1;

    	//Find total length of all sequences and initialize nextIndex array
    	int[] nextIndex = new int[numSeqs];

       	for(int i = 0; i < numSeqs; i++)
    		nextIndex[i] = 0;
    	
    	Node lastV = null;
    	
    	int alignLength = alignedSeqs[0].length();
        this.V = new Node[alignLength];
        this.activeNodes = new HashSet<Node>();
        
        //Ensure rows are flush.
        for(int i = 0; i < numSeqs; i++){
	    if (alignedSeqs[i].length() != alignLength){
		System.err.println("ERROR: Alternate alignment rows are not flush.");
		System.exit(-1);
	    }
        }
        
        //Ensure sequences match.
    	for (int i = 0; i < alignLength; i++){
	    int[] seqIndices = new int[numSeqs];
	    int presentSeqs = 0;
	    int seqI = -1;
    		
	    for(int j = 0; j < numSeqs; j++){
		char c = alignedSeqs[j].charAt(i);
		if (c == '-'){
		    seqIndices[j] = -1;
		}
		else if (sequences[j].length() > nextIndex[j] && c == sequences[j].charAt(nextIndex[j])){
		    seqIndices[j] = nextIndex[j]++;
		    presentSeqs++;
		    seqI = j;
		}
		else{
		    System.err.println("ERROR: Sequences in alternate alignment do not match; offending sequence is '" + keys.get(j) + "'.");
		    System.exit(-1);
		}
	    }
    		
	    Node v = new Node(seqIndices, seqI, pm);
	    v.isMerge = (presentSeqs > 1);
	    v.setColorScheme (colorScheme);
	    this.V[i] = v;
	    this.mergeDict.put(v,v);

	    this.acc += v.getAcc();
	    this.sps += v.getSPS();
	    this.ppv += v.getPPV();
	    this.cert += v.getCert();
	    this.cons += v.getCons();

	    this.accDenom += v.getAccDenom();
	    this.spsDenom += v.getSPSDenom();
	    this.ppvDenom += v.getPPVDenom();
	    this.certDenom += v.getCertDenom();
	    this.consDenom += v.getConsDenom();

	    this.activeNodes.add(v);
    		
	    if (lastV != null)
            	lastV.fwdEdge = v;
            lastV = v;
	    
    	}
    	
    	//Ensure we've reached the end of the original sequence.
    	for (int i = 0; i < numSeqs; i++){
	    if (nextIndex[i] != sequences[i].length()){
		System.err.println("ERROR: Sequences in alternate alignment do not match; offending sequence is '" + keys.get(i) + "'.");
		System.exit(-1);	
	    }
    	}
    	
    	this.startV = V;

    }

    private void buildDagFromInit(List<String> keys, String[] sequences, HashMap<String, Integer> initDAG, ProbabilityMatrices pm){
        Node lastV;
        String seq;

	acc = 0;
	sps = 0;
	ppv = 0;
	cert = 0;
	cons = 0;

	accDenom = 0;
	spsDenom = 0;
	ppvDenom = 0;
	certDenom = 0;
	consDenom = 0;

	deltaAccNorm = 0;
	deltaSPSNorm = 0;
	deltaPPVNorm = 0;
	deltaCertNorm = 0;
	deltaConsNorm = 0;

	gapFactorCurrent = -1;

       	//Find total length of all sequences
        int totalLength = 0;
        for (int i = 0; i < this.numSeqs; i++){
        	totalLength += sequences[i].length();
        }
        
        this.V = new Node[totalLength];
        this.oV = new Node[totalLength];
        this.startV = new Node[numSeqs];
        this.activeNodes = new HashSet<Node>();
        
       	int k = 0;
        for (int i = 0; i < this.numSeqs; i++){ // iterate through array of sequences
            lastV = null;
            seq = sequences[i];
            
            for (int j = 0; j < seq.length(); j++){ // iterate through sequence
        		
        		// Initialize empty column and indices arrays
                Node v = new Node(i, j, pm);
		v.setColorScheme (colorScheme);
                this.V[k] = v;
		//                this.oV[initDAG.get(i + "." + (j+1))] = v; //j+1 because index starts at 1
                this.oV[initDAG.get(i + "." + j)] = v; // j because index starts at 0
                this.mergeDict.put(v,v);

		this.acc += v.getAcc();
		this.sps += v.getSPS();
		this.ppv += v.getPPV();
		this.cert += v.getCert();
		this.cons += v.getCons();

		this.accDenom += v.getAccDenom();
		this.spsDenom += v.getSPSDenom();
		this.ppvDenom += v.getPPVDenom();
		this.certDenom += v.getCertDenom();
		this.consDenom += v.getConsDenom();

                this.activeNodes.add(v);
                if (lastV != null)
                    lastV.fwdEdge = v;
                else
                	startV[i] = v;
                lastV = v;
           		k++;
            }
        }
    }
    
    /// Get normalized accuracy for alignment.
    /*
     * \return accuracy, or -1 if undefined
     * Performs bounds-checking.
     */
    public double getAlignAccNorm(){
	if (accDenom == 0.)
	    return -1;
    	return SparseMatrix.saneProb (acc / accDenom);
    }

    /// Get normalized SPS for alignment.
    /*
     * \return SPS, or -1 if undefined
     */
    public double getAlignSPSNorm(){
	if (spsDenom == 0.)
	    return -1;
    	return SparseMatrix.saneProb (sps / spsDenom);
    }

    /// Get normalized PPV for alignment.
    /*
     * \return PPV, or -1 if undefined
     */
    public double getAlignPPVNorm(){
	if (ppvDenom == 0.)
	    return -1;
    	return SparseMatrix.saneProb (ppv / ppvDenom);
    }

    /// Get normalized certainty for alignment.
    /*
     * \return certainty, or -1 if undefined
     */
    public double getAlignCertNorm(){
	if (certDenom == 0.)
	    return -1;
    	return SparseMatrix.saneProb (Node.normalizedCertainty (cert / certDenom));
    }

    /// Get normalized consistency for alignment.
    /*
     * \return consistency, or -1 if undefined
     */
    public double getAlignConsNorm(){
	if (consDenom == 0.)
	    return -1;
    	return SparseMatrix.saneProb (cons / consDenom);
    }

    /// Get the change in normalized accuracy.
    public double getDeltaAccNorm(){
	return deltaAccNorm;
    }

    /// Get the change in normalized SPS.
    public double getDeltaSPSNorm(){
	return deltaSPSNorm;
    }

    /// Get the change in normalized PPV.
    public double getDeltaPPVNorm(){
	return deltaPPVNorm;
    }

    /// Get the change in normalized certainty.
    public double getDeltaCertNorm(){
	return deltaCertNorm;
    }

    /// Get the change in normalized consistency.
    public double getDeltaConsNorm(){
	return deltaConsNorm;
    }

    /// Get the implicit current gap factor.
    /*
     * -1 if Infinity, implicit current gap factor otherwise
     */
    public double getGapFactor(){
	return gapFactorCurrent;
    }

    /// Merge two nodes.
    /*
     * @param i index of source node in oV
     * @param j index of dest node in oV
     */
    public void merge(int i, int j){
	    
	Node n1 = this.mergeDict.get(this.oV[i]);
        Node n2 = this.mergeDict.get(this.oV[j]);
        
        Node m = new Node(n1, n2, mergeOrder++);
	m.setColorScheme (colorScheme);

	// now update accuracies after the merge
	// note that the calculation must be done on the unnormalized accuracies
	deltaAccNorm = m.getAcc() - (n1.getAcc() + n2.getAcc());
	deltaAccNorm /= (1.0 * m.getAccDenom());

	deltaSPSNorm = m.getSPS() - (n1.getSPS() + n2.getSPS());
	deltaSPSNorm /= (1.0 * m.getSPSDenom());

	deltaPPVNorm = m.getPPV() - (n1.getPPV() + n2.getPPV());
	deltaPPVNorm /= (1.0 * m.getPPVDenom());

	deltaCertNorm = m.getCert() - (n1.getCert() + n2.getCert());
	deltaCertNorm /= (1.0 * m.getCertDenom());

	deltaConsNorm = m.getCons() - (n1.getCons() + n2.getCons());
	deltaConsNorm /= (1.0 * m.getConsDenom());

	gapFactorCurrent = m.getWeightTgf();

//	System.err.println ("m.getSPS() = " + m.getSPS()
//			    + "m.getSPSDenom() = " + m.getSPSDenom()
//			    + "; n1.getSPS() = " + n1.getSPS() + "; n2.getSPS() = " + n2.getSPS()
//			    + "; n1.getSPSDenom() = " + n1.getSPSDenom() + "; n2.getSPSDenom() = " + n2.getSPSDenom());
	
	acc += m.getAcc() - (n1.getAcc() + n2.getAcc());
	sps += m.getSPS() - (n1.getSPS() + n2.getSPS());
	ppv += m.getPPV() - (n1.getPPV() + n2.getPPV());
	cert += m.getCert() - (n1.getCert() + n2.getCert());
	cons += m.getCons() - (n1.getCons() + n2.getCons());

	accDenom += m.getAccDenom() - (n1.getAccDenom() + n2.getAccDenom());
	spsDenom += m.getSPSDenom() - (n1.getSPSDenom() + n2.getSPSDenom());
	ppvDenom += m.getPPVDenom() - (n1.getPPVDenom() + n2.getPPVDenom());
	certDenom += m.getCertDenom() - (n1.getCertDenom() + n2.getCertDenom());
	consDenom += m.getConsDenom() - (n1.getConsDenom() + n2.getConsDenom());
        
        activeNodes.remove(n1);
        activeNodes.remove(n2);
        activeNodes.add(m);
        updateMergeDict(m, m);
    }
    
    /// Demerge a nodes.
    /*
     * @param i index of node to demerge in oV
     * @param k index of node in oV
     *
     * The index k is necessary in order to update the gap factor
     * after demerging a node.
     */
    public void demerge(int i, int k){
        
	--mergeOrder;

        Node m = this.mergeDict.get(this.oV[i]);
        Node n1 = m.n1;
        Node n2 = m.n2;

	n1.setColorScheme (colorScheme);
	n2.setColorScheme (colorScheme);
        
	// this calculation must be done on the unnormalized accuracies
        deltaAccNorm = (n1.getAcc() + n2.getAcc()) - m.getAcc();
	deltaAccNorm /= (1.0 * m.getAccDenom());

        deltaSPSNorm = (n1.getSPS() + n2.getSPS()) - m.getSPS();
	deltaSPSNorm /= (1.0 * m.getSPSDenom());

        deltaPPVNorm = (n1.getPPV() + n2.getPPV()) - m.getPPV();
	deltaPPVNorm /= (1.0 * m.getPPVDenom());

        deltaCertNorm = (n1.getCert() + n2.getCert()) - m.getCert();
	deltaCertNorm /= (1.0 * m.getCertDenom());

        deltaConsNorm = (n1.getCons() + n2.getCons()) - m.getCons();
	deltaConsNorm /= (1.0 * m.getConsDenom());

	// decrement the gap factor: look up the gap factor for the previous merge
	if (k >= 0)
	    gapFactorCurrent = (this.mergeDict.get(this.oV[k])).getWeightTgf();
	else // if no previous merge, set to dummy value
	    gapFactorCurrent = -1;

	acc += (n1.getAcc() + n2.getAcc()) - m.getAcc();
	sps += (n1.getSPS() + n2.getSPS()) - m.getSPS();
	ppv += (n1.getPPV() + n2.getPPV()) - m.getPPV();
	cert += (n1.getCert() + n2.getCert()) - m.getCert();
	cons += (n1.getCons() + n2.getCons()) - m.getCons();

	accDenom += (n1.getAccDenom() + n2.getAccDenom()) - m.getAccDenom();
	spsDenom += (n1.getSPSDenom() + n2.getSPSDenom()) - m.getSPSDenom();
	ppvDenom += (n1.getPPVDenom() + n2.getPPVDenom()) - m.getPPVDenom();
	certDenom += (n1.getCertDenom() + n2.getCertDenom()) - m.getCertDenom();
	consDenom += (n1.getConsDenom() + n2.getConsDenom()) - m.getConsDenom();
        
	if (mergeDict.containsKey(m)){
	    System.out.println("found it!");
	}
	
	activeNodes.remove(m);
	activeNodes.add(n1);
	activeNodes.add(n2);
        updateMergeDict(n1, n1);
        updateMergeDict(n2, n2);
    }
    
    private void updateMergeDict(Node n, Node m){
        if (n.isMerge){
            Node n1 = n.n1;
            Node n2 = n.n2;
            
            updateMergeDict(n1, m);
            updateMergeDict(n2, m);
        }
        else 
            this.mergeDict.put(n, m);
    }
    
    public Node[] topoSort(){
    	Stack<Node> toVisit = new Stack<Node>();
    	Node u, v, w;
        Comparator<Node> byRevFinishTime = new NodeFinishTimeComparator();
    	int time = 0;
    	
        Iterator<Node> iter = activeNodes.iterator();
        while (iter.hasNext()) {
        	v = iter.next();
        	v.visited = false;
        	v.finishTime = -1;
    	}
    	
    	for (int i = 0; i < numSeqs; i++){
			u = mergeDict.get(startV[i]);

//    		System.out.print ("outer loop: ");
//			u.output();
			
			if (!u.visited){
				toVisit.push(u);
//	    		System.out.print ("  pushed: ");
//	    		u.output();
			}
			
			while (!toVisit.empty()){
				v = toVisit.pop();
//				System.out.print (" popped: ");
//				v.output();
				if (v.visited){
					if (v.finishTime < 0) // if not finished
						v.finishTime = time++;
				}
				else{
					v.visited = true;
					toVisit.push(v);
//		    		System.out.print ("  pushed: ");
//		    		v.output();
		    		
		    		Node[] edges = processEdges(v.getFwdEdges());
					for(int j = 0; j < edges.length; j++) {
						w  = edges[j];
						
						if (!w.visited){
							toVisit.push(w);			
//							System.out.print ("   pushed ");
//							w.output();
						}
					}
				}
			}
        }

		iter = activeNodes.iterator();
		
		Node[] L = activeNodes.toArray(new Node[0]);
		Arrays.sort(L, byRevFinishTime);
        return L;
    }    
    
    // Maps edges to their merged version and sorts array
    private Node[] processEdges(Collection<Node> edges){
    	Node[] edgeArray = new Node[edges.size()];
    	Iterator<Node> iter = edges.iterator();
    	
    	int i = 0;
    	while (iter.hasNext())
    		edgeArray[i++] = mergeDict.get(iter.next());
    	
    	if (edgeArray.length > 1){
    		Comparator<Node> byNodeSequence = new NodeSequenceComparator();
    		Arrays.sort(edgeArray, byNodeSequence);
    		MySorter sorter = new MySorter();
			edgeArray = sorter.sort(edgeArray);
    	}
    	return edgeArray;	
    }
    

    
    private void checkSort(mad.Node[] nArray){
			int[] maxArray = new int[numSeqs];
			Node v;
			
			for(int i = 0; i < numSeqs; i++){
				maxArray[i] = -1;
			}
			for(int i = 0; i < nArray.length; i++){
				v = nArray[i];
				int[] seqIndices = v.getSeqIndices();
				for(int j = 0; j < numSeqs; j++){
					if (seqIndices[j] != -1 && maxArray[j] > seqIndices[j]){
						System.err.println("Sanity check failed!");
						System.exit(-1);
					}
					else if (seqIndices[j] != -1)
						maxArray[j] = seqIndices[j];
				}
			}
    }
}

class MySorter{
	
	
	public MySorter(){}
	
	public Node[] sort(Node[] nArray){
		boolean isSorted = false;
		int[] maxArray;
		int numSeqs = nArray[0].pm.numSeqs;
		Node v;
		
		while (!isSorted){
			isSorted = true;
			maxArray = new int[numSeqs];
			for(int i = 0; i < numSeqs; i++){
				maxArray[i] = -1;
			}
			for(int i = 0; i < nArray.length; i++){
				v = nArray[i];
				boolean switched = false;
				int[] seqIndices = v.getSeqIndices();
				for(int j = 0; j < numSeqs; j++){
					if (seqIndices[j] != -1 && maxArray[j] > seqIndices[j]){
						if (!switched){
							switchNodes(nArray, i-1, i);
							switched = true;
							isSorted = false;
						}
					}
					else if (seqIndices[j] != -1)
						maxArray[j] = seqIndices[j];
				}
			}
		}	
		return nArray;
	}
	
	
	private Node[] switchNodes(Node[] nArray, int i, int j){
		Node temp = nArray[i];
		nArray[i] = nArray[j];
		nArray[j] = temp;
		return nArray;
	}
}

class NodeFinishTimeComparator implements Comparator<Node> {
	
	public int compare(Node n1, Node n2){
		if (n1.finishTime > n2.finishTime)
			return -1;
		return 1;
	}
}

class NodeSequenceComparator implements Comparator<Node>{
	
	public int compare(Node n1, Node n2){
		int[] seq1Indices = n1.getSeqIndices();
		int[] seq2Indices = n2.getSeqIndices();
		for (int i = 0; i < seq1Indices.length; i++){

			if (seq1Indices[i] != -1 && seq2Indices[i] != -1)
				if (seq1Indices[i] < seq2Indices[i])
					return -1;
				else
					return 1;	
		}
		return 0;
	}
}
