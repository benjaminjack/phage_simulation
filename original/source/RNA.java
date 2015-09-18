import java.util.Vector;
/**
 * General class for RNA molecules.
 *
 * @author Sriram Kosuri
 * @version 1.0
 * @see TabascoSimulator
 */
public class RNA implements Molecule{
    /**
	 * The most upstream position of the RNA in the DNA coordinates from where it is transcribed.
	 */
	int start;
	/**
	 * The most downstream position of the RNA in the DNA coordinates from where it is transcribed.
	 */	
	int stop;
	/**
	 * The ID of the RNA molecule.
	 */
	int ID;
	/**
	 * The number of copies of the RNA molecule.
	 */
	private int copyNumber;
    
	/**
	 * A constructor for an RNA molecule
	 * @param start The most upstream position of the RNA in the DNA coordinates from where it is transcribed.
	 * @param stop The most downstream position of the RNA in the DNA coordinates from where it is transcribed.
	 */
    public RNA(int start, int stop){
        this.start = start;
        this.stop = stop;
    }

	/**
		* A constructor for an RNA molecule
	 * @param start The most upstream position of the RNA in the DNA coordinates from where it is transcribed.
	 * @param stop The most downstream position of the RNA in the DNA coordinates from where it is transcribed.
	 * @param ID the ID of the RNA molecule
	 */	
    public RNA(int start, int stop, int ID){
        this(start, stop);
        this.ID = ID;
    }

	/**
		* A constructor for an RNA molecule
	 * @param ID the ID of the RNA molecule
	 */	
    public RNA(int ID){
        this.ID = ID;
    }
    
	/**
	 * Returns the most upstream position of the RNA in the DNA coordinates from where it is transcribed.
	 * @return the most upstream position of the RNA in the DNA coordinates from where it is transcribed.
	 */
    public int getStart(){
        return start;
    }
	
	/**
	 * Increments the number of copies of RNA by inc.
	 * @param inc The number of copies of RNA to increment the current copy number by.
	 */
    public void incrementCopyNumber(int inc){
        copyNumber += inc;
    }

	/**
	 * Returns the current copy number of the RNA
	 * @return the current copy number of the RNA.
	 */
    public int getCopyNumber(){
        return copyNumber;
    }
    
	/**
	 * Returns the most dowsntream position of the RNA in the DNA coordinates from where it is transcribed.
	 * @return the most downstream position of the RNA in the DNA coordinates from where it is transcribed.
	 */
    public int getStop(){
        return stop;
    }

	/**
	 * Returns the current length of the RNA.  Simply returns (stop-start+1).
	 */
    public int getLength(){
        return stop-start+1;
    }
    
	/**
	 * Returns the ID of the RNA.
	 * @return the ID of the RNA.
	 */
    public int getID(){
        return ID;
    }

    
}

/**
 * A class used for tracking RNA that is still attached to a polymerase.  
 * @author Sriram  Kosuri 
 * @version 1.0
 * @see RNA
 * @see DPComplex
 */
class NascentRNA extends RNA {
	/**
	 * The transcribing polymerase that this RNA is attached to.
	 */
    DPComplex dp;
    /**
     * The constructor for a Nascent RNA. 
	 * @param start The most upstream position of the RNA in the DNA coordinates from where it is transcribed.
	 * @param stop The most downstream position of the RNA in the DNA coordinates from where it is transcribed.
	 * @param dp The transcribing polymerase that this RNA is attached to.
	 */
    public NascentRNA(int start, int stop, DPComplex dp){
        super(start,stop);
        this.dp = dp;
    }
    
}

/**
 * A class that represents ribosome binding sites that are in a cell and no longer bound to the DNA and polymerase that is being transcribed.  A NascentRBS is converted
 * to a MatureRBS upon the RNA polymerase reaching the stop site of translation.
 * @author Sriram  Kosuri 
 * @version 1.0
 * @see RNA
 * @see NascentRBS
 */
class MatureRBS extends RNA implements Molecule{
    /**
	 * The position on the DNA where translation begins.
	 */
	int startsite;
	/**
     * The reaction to be called upon the creation of a ribosome binding site usually from a ribosome that is releasing it from a bound state.
	 */
    Reaction entryReaction;
	/**
	 * The reaction to be called upon conversion of a nascentRBS to a matureRBS.
	 */	
    Reaction conversionReaction;
	/**
	 * The translation stop that corresponds to this ribosome binding site.
	 */
	TranslationStop ts;
	/**
	 * The number of bases that are needed for the ribosome to elongate past to clear the ribosome binding site.
	 */
	int initstepsize;
	/**
	 * The number of bases that are needed for the ribosome to go from the initial jump to the translation stop.
	 */
	int elongstepsize;
	/**
	 * The nascentRBS that corresponds to this matureRBS.  The nascentRBS
	 */
	NascentRBS nrbs;
	/**
	 * The second order binding constant of the ribosome to this ribosome binding site.
	 */
	double strength;
	/**
	 * The degradation rate of the RBS.  Degradation of mRNA occurs by degrading ribosome binding sites.  
	 */
	double mrnadeg;
	/**
     * Upon binding of the mRBS, the complex is converted to an initrna.
	 */
	GammaJumpRNA initrna;
	/**
	 * The reaction that is called upon binding of the mRBS by the ribosome.
	 */
	Reaction mrbsBindingRxn;
	/**
	 * The first-order rate constant determining the time it takes for the ribosome to form a stable elongation complex.
	 */
	double initRateConstant;
    
	/**
	 * The constructor for a MatureRBS.   Only called once at the start of simulation for each coding domain.
	 * @param start The most upstream position of the mRBS on the DNA.
	 * @param stop The most downstream position of the mRBS on the DNA
	 * @param startsite The startsite of translation
	 * @param ID The ID of the RBS
	 * @param ts The translation stop site that corresponds to this mRBS
	 * @param in Defines the initstepsize; the number of bases that are needed for the ribosome to go from the initial jump to the translation stop.
	 * @param el Defines the elongstepsize; the number of bases that are needed for the ribosome to go from the initial jump to the translation stop.
	 * @param strength The second order binding constant of the ribosome to this ribosome binding site.
	 * @param mrnadeg The degradation rate of the RBS.  Degradation of mRNA occurs by degrading ribosome binding sites.
	 * @param initRateConstant The first-order rate constant determining the time it takes for the ribosome to form a stable elongation complex.
	 */
    public MatureRBS(int start, int stop, int startsite, int ID,TranslationStop ts,int in, int el,double strength,double mrnadeg,double initRateConstant){
        super(start,stop,ID);
        this.startsite = startsite;
        this.ts = ts;
        this.initstepsize = in;
        this.elongstepsize = el;
        this.strength = strength;
        this.mrnadeg = mrnadeg;
		this.initRateConstant = initRateConstant;
    }

	/**
	 * Defines the Entry reaction, which causes an matureRBS to increment in quantity
	 * @param r The reaction to be set as the EntryReaction.
	 */
    public void setEntryReaction(Reaction r){
        entryReaction = r;
    }

	/**
	 * Defines the Converstion reaction, which causes a nascentRBS to be converted to a matureRBS.
	 * @param r The reaction to be set as the converstion reaction.
	 */
    public void setConversionReaction(Reaction r){
        conversionReaction = r;
    }

}
/**
 * A helper class to be used as entries into the NascentRBS#myRBSInstances.  The purpose of this class is to track instantiations of nascentRBS and the DPComplex 
 * that is currently transcribing the mRNA that the nascentRBS is on.  In addition, if the nascentRBS is bound, the RRComplex associated with it is also tracked.
 * @version 1.0
 * @see RNA
 * @see DPComplex
 * @see RRComplex
 * @see NascentRBS
 */
class RBSInstance {
    /**
	 * The DPComplex that is currently transcribing the RNA that this RBSInstance is encoded on.
	 */
	DPComplex dp;
    /**
	 * If this RBSInstance is bound by a ribosome, rr points to the RRComplex that is associated with it.
	 */
    RRComplex rr;
    /**
	 * Tracks whether the RBSInstance is bound or not by a ribosome.
	 */
    boolean isBound;
    
	/**
	 * The constructor to create an RBSInstance.
	 * @param dp The DPComplex that is currently tarnscribing the RNA that this RBSInstance is encoded on.
	 */
    public RBSInstance(DPComplex dp){
        this.dp = dp;
        rr = null;
        isBound = false;
    }

    
}

/**
 * A class that represents ribosome binding sites that are on RNA that is still being transcribed.  Each instance of an nascentRBS is tracked individually so that if
 * the open reading frame being transcribed is not completed, translation from the nascentRBS can be stopped.
 * @version 1.0
 * @see RNA
 * @see DPComplex
 * @see RRComplex
 */
class NascentRBS implements Molecule{
	/**
	 * The matureRBS that this nascentRBS corresponds to.  Upon hitting the corresponding translation stop, nascentRBSs are converted to matureRBSs.
	 */
    MatureRBS mrbs;
	/**
	 * The number of copies of the nascentRBS that are not bound.
	 */
    int copyNumber;
	/**
	 * The reaction to call upon creation of a nascentRBS
	 */
    NascentRBSEntryRxn entryReaction;
	/**
	 * A vector tracking all the existing corresponding nascentRBSs.  This vector contains a series of RBSInstance instances.  The total of this vector is the number of bound
	 * and unbound nascentRBSs.
	 */
    Vector myRBSInstances;
	/**
	 * The reaction to call upon loss of a nascentRBS.
	 */
    NascentRBSLossRxn lossReaction;
	/**
	 * The ID of this nascentRBS.
	 */
    int ID;
    
	/**
	 * The constructor of an NascentRBS.  Only called once at the start of simulation for each coding domain.
	 */
    public NascentRBS(int ID, MatureRBS r){
        this.mrbs = r;
        copyNumber = 0;
        myRBSInstances = new Vector();
    }
	
	/**
	 * Returns the ID of the NascentRBS.
	 */
    public int getID(){
        return ID;
    }
    
	/**
	 * Increments the copy number of the nascentRBS.  Be careful to first create an RBSInstance before calling this method.
	 * @param inc The number of nascentRBS to increment by.  This number can be negative to decrement.
	 */
	
    public void incrementCopyNumber(int inc){
        copyNumber +=inc;
    }

	/**
	 * Returns the copy number of unbound nascentRBSs.
	 * @return the copy number of unbound nascentRBSs
	 */
    public int getCopyNumber(){
        return copyNumber;
    }

}

/**
 * A site that tracks the translation stop site for an open reading frame.  Each ribosome binding site has a corresonding translation stop site.
 * @version 1.0
 * @see RNA
 * @see NascentRBS
 * @see MatureRBS
 */
class TranslationStop extends RNA {
    /**
	 * The stopsite of the translation stop in DNA coordinates.
	 */
	int stopsite;
	/**
     * The protein product that is created upon translation of the open reading frame.
	 */
    Protein product;
	/**
	 * The matureRBS that corresponds to open reading frame here.  This link allows overlapping gene products (rather than relying on the previously transcribing matureRBS)
	 */
    MatureRBS mrbs;

	/**
	 * The constructor for the translation stop.  Called once at the start of simulation for each coding domain.
	 * @param stopsite the site on the DNA where translation stops.
	 * @param prod The protein product encoded by coding domain ending in this translation stop
	 * @param ID The ID of the stopsite
	 */
    public TranslationStop(int stopsite, Protein prod, int ID){
        super(stopsite,stopsite,ID);
        this.stopsite = stopsite;
        this.product = prod;
    }

}

/**
 * A class representing RNA complexes with ribosomes that are elongated in jumps simulated by the Gamma Distribution.  Tracks initiation complexes and elongation complexes of 
 * RNA and ribosome.
 * @version 1.0
 * @see RNA
 */
class GammaJumpRNA extends RNA{
    //must set this before start
    /**
	 * The reaction to call upon initiation or elongation completing.
	 */
	GammaJumpRxn rxnToUpdate;
    //initiation or elongation complex (true or false)
    /**
	 * Indicates whether the complex is an initiation (true) or elongation (false) complex.
	 */
	boolean initiation;
    /**
	 * The Mature Ribosome binding site that is bound if this is an initiation complex.  Also, it is used to determine the protein product if elongation occurs.
	 */
    MatureRBS mrbs;
    /**
	 * A pointer to the simulation.
	 */
    TabascoSimulator mySim;
	
	/**
	 * The constructor for a GammaJumpRNA.  
	 * @param ID the ID of the GammaJumpRNA
	 * @param mrbs The mature RBS associated with this GammaJumpRNA
	 * @param init Determines whether this is an initiation or elongation complex.
	 * @param mySim A pointer to the simulation.
	 */
    public GammaJumpRNA(int ID, MatureRBS mrbs, boolean init, TabascoSimulator mySim){
        super(ID);
        this.mrbs = mrbs;
        this.initiation = init;
		this.mySim = mySim;
	}
    
	/**
	 * Increments the copynumber.
	 * @inc the Number to increment the copy number by (can be negative to decrement).
	 */
    public void incrementCopyNumber(int inc){
        super.incrementCopyNumber(inc);
        rxnToUpdate.incremented(inc);
		Reaction r = (Reaction) rxnToUpdate;
		for(int i=0;i<(r.updateRxn).size();i++){
			int reacid = ((Integer) r.updateRxn.elementAt(i)).intValue();
			mySim.gint.myIPQ.updateOther(reacid, mySim.gint.t);
		}
    }

}
