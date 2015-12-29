/**
 * A class representing DNA-Protein complexes.  
 * 
 * This class is used to represent all the DNA protein complexes in the cell.  It contains pointers to the protein
 * and the dna position where it is located.  The DPComplex can be used to represent both polymerases and other 
 * protein-DNA complexes such as repressors.
 *
 * @author Sriram Kosuri
 * @version 1.0
 */

class DPComplex implements Comparable{
    /**
	 * This ID is automatically generated upon generation of the DPComplex and is important in visualization.
	 */
    String ID;
    
	/**
	 * The protein that is in the complex.  Often is a polymerase.
	 */
	Protein myProtein;
	
    /**
	 * The startsite of the protein.  Usually most upstream site position.
	 */
	int start;
	
	/**
	 * The stopsite of the protein.  Usually the most downstream site position.
	 */
    int stop;
    
	/**
	 * The time that this complex is next going to undergo a state change or reaction.  Usually an elongation step
	 */
    double time;
    
	/**
	 * The rate constant of normal elongation of this complex
	 */
	double a;
	
	/**
	 * The pointer to the nascentRNA that is created upon binding.
	 * @see NascentRNA
	 */
    NascentRNA m;
	
	/**
	 * A vector that contains RRComplexes that have finished elongating, but are waiting for the polymerase to finish elongating as well.
	 * @see RRComplex
	 */
    java.util.Vector stalledRibosomes;
	
	/**
	 * A variable tracking state during initiation from a promoter.  0 indicates that the complex is in the process of choosing whether 
	 * to fall off or form an open complex.  1 indicates the complex is deciding between abortive initiation and elongation.  3 indicates
	 * the polymerase is elongating normally.
	 */
    int initState;
	
	/**
	 * A variable that tracks which reaction RNAPElongationEvent should execute when the DPComplex is released from the queue.  0 indicates
	 * that the polymerase has already initialized and is elongating.  2 indicates
	 * that the complex is immediately to fall off.  3 indicates if the polymerase is forming an open complex.  4 indicates the polymerase 
	 * is recycling to the open complex.  5 indicates successfully initiated polymerase
	 *
	 * @see RNAPElongationEvent#move()
	 */
    int whichReaction;
	
	/**
	 * The promoter that this complex initially formed on.  All DNA protein complexes must bind at a promoter, even if it is acting as 
	 * just a protein-DNA interaction.
	 *
	 * @see myPromoter
	 */
    Promoter myPromoter;

	/**
	 * The phage or piece of DNA that this complex resides on.
	 *
	 * @see Phage
	 */
    Phage myPhage;
	
	/**
	 * The reaction that is called when the protein releases the DNA.
	 *
	 */
    DPComplexReleaseRxn releaseRxn;
    
	/**
	 * This variable is true when an upstream complex causes the protein to release teh DNA because of a collision.
	 */
	boolean released; 
	
	/**
	 * The only constructor to create a DPComplex.
	 * @param id The ID of the DPComplex.
	 * @param prot The protein that is binding the DNA.
	 * @param myPromoter The promoter at which the protein is binding.
	 * @param start The most downstream site of the protein on the DNA
	 * @param stop The most upstream site of the protein on the DNA
	 * @param a The elongation rate of this complex
	 * @param rand A random number to be able to calculate when the complex will react next 
	 * @param time The current time the Complex is created
	 * @param myPhage the Phage/DNA that the promoter and protein are on.
	 */
    public DPComplex(String id, Protein prot, Promoter myPromoter, int start,int stop, double a, double rand, double time, Phage myPhage){
        ID = id;
        myProtein = prot;
        this.start = start;
        this.stop = stop;
        this.a = a;
        stalledRibosomes = new java.util.Vector();
        initState = 0;
        this.myPromoter = myPromoter;
        this.myPhage = myPhage;
        whichReaction = 0;
        calculateTime(rand, time);
        for(int i=stop;i<=start;i++){
            myPhage.myDNAComplexes[i] = this;
        }
        releaseRxn = myPhage.releaseRxn[myProtein.getID()-1];
        released = false;
    }
	
	/**
	 * Compares two DPComplexes' times to see which should execute first.
	 * @param r The DPComplex to compare to
	 * @return -1 if this complex should execute first, 1 if r should execute first, 0 if they are equal
	 */
    public int compareTo(Object r){
        DPComplex rhs = (DPComplex) r;
        if((this.time - rhs.time)<0)
            return -1;
        else if((this.time - rhs.time)>0)
            return 1;
        else
            return 0;
    }
    
	/**
	 * An method to access the DPComplexReleaseRxn
	 * @return the DPComplexReleaseRxn that should be called upon the release of this DPComplex
	 * @see DPComplexReleaseRxn
	 */
    public DPComplexReleaseRxn getReleaseRxn(){
        return releaseRxn;
    }
    
    /**
	 * The method that encodes some of the logic of initiation of a polymerase from a promoter, and eventually 
	 * elongation.
	 */
    public void calculateTime(double random, double t){
        
        if(initState==0){
            //calculate both on and off times.. whichever is least is what is going to react
            double aoff = myPromoter.aoff[myProtein.getID()-1];
            double ainiton = myPromoter.ainiton[myProtein.getID()-1];
            double time1 = (1/aoff)*Math.log(1/random) + t; // fall off
            double time2 = (1/ainiton)*Math.log(1/myPhage.myCell.mySim.gint.myIPQ.getRandomDouble()) + t; //initiation
            double runoffPercent = myPromoter.runoffPercent[myProtein.getID()-1];
            
			//decide if going to fall off or inititalize            
            if(time1<time2){
				//fall off
                time = time1;
                whichReaction = 2;
            }else{
				//initialize
				time = time2;
				whichReaction = 3;
            }
        }else if(initState==1){
			//decide if it is going to abortive initiation
			double aelong = myPromoter.aelongstart[myProtein.getID()-1];
			double arecyc = myPromoter.arecyc[myProtein.getID()-1];
            double time1 = (1/(aelong))*Math.log(1/myPhage.myCell.mySim.gint.myIPQ.getRandomDouble()) + t; // elongate
			double time2 = (1/(arecyc))*Math.log(1/myPhage.myCell.mySim.gint.myIPQ.getRandomDouble()) + t; // recycle

			double runoffPercent = myPromoter.runoffPercent[myProtein.getID()-1];

			if (myPhage.myCell.mySim.gint.myIPQ.getRandomDouble()<=(double)(runoffPercent/100)){
				//elongate
				time = time1;
				whichReaction = 5;
			}else{
				//recycle
				time = time2;
				whichReaction = 4;
			}
        }else{
            //called when normal elongation is taking place initState = 3
            time = (1/a)*Math.log(1/random) + t;
        }
    }
    
}

/**
 * This class is used to used to create instances of ribosomes elongating on nascentRNA.
 *
 * @author Sriram Kosuri
 * @see NascentRNA
 * @version 1.0
 */

class RRComplex implements Comparable{
    
	
	/**
	 * The number of bases that the ribosome takes to clear the ribosome binding site.
	 */
    int initstepsize;

    /**
	 * The number of bases that the ribosome must traverse to go from clearance of the RBS to the termination site.
	 */
	int elongstepsize;
	
	/**
	 * The time that this complex is next going to undergo a state change or reaction.  
	 */
	double time;
    
	/**
	 * The rate constant representing the process of the ribosome binding, assembling, and initiating on the RBS
	 */
	double a;
    
	/**
	 * Indicates if this Complex is in the initiation state (still occupying an RBS).
	 */
	boolean init;
	
	/**
	 * The NascentRBS that the ribosome initially bound to.
	 */
    NascentRBS rbs;
	
	/**
	 * An indicator of whether the RNA polymerase that is transcribing the nascentRNA has transcribed the entire coding domain being translated
	 */
    boolean totallyentered;
    
	/**
	 * The NascentRNA that the ribosome is translating along.
	 */
	NascentRNA nrna;
	/**
	 * Indicates if this complex is still in the pre-initiation phase (ribosome assembling on the RNA).
	 */
	boolean pre_init;
	
	/**
	 * The only constructor to create an RRComplex.
	 * 
	 * @param rbs
	 * @param a
	 * @param initstepsize
	 * @param elongstepsize
	 */
    public RRComplex(NascentRBS rbs,double a, int initstepsize, int elongstepsize){
        this.a = a;
        this.init = true;
        this.initstepsize = initstepsize;
        this.elongstepsize = elongstepsize;
        this.rbs = rbs;
        totallyentered = false;
        nrna = null;
		pre_init = true;
    }
	
	/**
	 * Compares two RRComplexes' times to see which should execute first.
	 * @param r The RRComplex to compare to
	 * @return -1 if this complex should execute first, 1 if r should execute first, 0 if they are equal
	 */
    public int compareTo(Object r){
        RRComplex rhs = (RRComplex) r;
        if((this.time - rhs.time)<0)
            return -1;
        else if((this.time - rhs.time)>0)
            return 1;
        else
            return 0;
    }
}
