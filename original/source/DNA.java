/**
 * Generic class that describes features on a piece of DNA (or Phage).
 *
 * @author Sriram Kosuri
 * @version 1.0
 * @see Phage
 */
public class DNA implements Molecule{
	/**
	 * The startsite of the DNA feature.  Usually most upstream site position.
	 */
	int start;
	/**
	 * The stopsite of the DNA feature.  Usually the most downstream site position.
	 */
	int stop;
	/**
	 * The length of the DNA feature.
	 */
	int length;
	/**
	 * The ID of the feature
	 */
	int ID;
	/**
	 * The number of copies of this DNA feature.  This is zero if the DNA feature has not entered the cell yet,
	 * or is blocked by other a DPComplex
	 * @see DPComplex
	 */
	int copyNumber;
	
	/**
	 * The reaction that needs to be called when the DNA feature enters the cell or is unblocked by a DPComplex
	 * @see DPComplex
	 */
	private Reaction entryReaction;
	
	/**
	 * The reaction that needs to be called when the DNA feature is blocked by a DPComplex
	 */
	private Reaction lossReaction;
    
	/**
	 * If the DNA feature is being blocked by more than one DPComplex, the second blocking will increment
	 * this value by 1.  If this value is greater than zero, a DPComplex that unblocks this feature will
	 * decrement this number.
	 */ 
	int stillBlocked;
	
	/**
	 * An indicator of whether this DNA feature has entered the cell.
	 */
    boolean hasEntered;
    
	/**
	 * The default constructure of a DNA feature.
	 * @param start the most upstream position of the DNA feature
	 * @param stop the most downstream position of the DNA feature.
	 * @param ID the ID of the DNA feature.
	 */
    public DNA(int start, int stop, int ID){
        this.start = start;
        this.stop = stop;
        this.length = stop-start+1;
        this.ID = ID;
        copyNumber = 0;
        stillBlocked = 0;
        hasEntered = false;
    }

	/**
	 * Retrieves the start position of the DNA feature.
	 * @return the start position of the DNA feature.
	 */
    public int getStart(){
        return start;
    }
	/**
	 * Retrieves the stop position of the DNA feature.
	 * @return the start position of the DNA feature
	 */
    public int getStop(){
        return stop;
    }
	/**
	 * Retrieves the copy number of the DNA feature.
	 * @return the copy number of the DNA feature
	 */
    public int getCopyNumber(){
        return copyNumber;
    }
	
	/**
	 * Increments the copy number of the DNA feature.
	 * @param inc the value to increment the current copyNumber by
	 */
    public void incrementCopyNumber(int inc){
        copyNumber += inc;
        if(copyNumber>1||copyNumber<0){
            System.out.println("ARG: In increment copy Number");
            System.exit(9);
        }
    }
	/**
	 * Retrieves the ID of this DNA feature.
	 * @return the ID of this DNA feature
	 */
    public int getID(){
        return ID;
    }
	/**
	 * Retrieves the reaction to execute upon blocking by a DPComplex.
	 * @return the reaction to execute upon blocking by a DPComplex
	 */
    public Reaction getLossReaction(){
        return lossReaction;
    }
	/**
	 * Sets the lossReaction as the reaction that is input.
	 * @param rec the reaction that is to be set as this DNA feature's lossReaction
	 */
    public void setLossReaction(Reaction rec){
        lossReaction = rec;
    }
	/**
	 * Retrieves the reaction to execute upon unblocking or entry of the DNA feature.
	 * @return the reaction to execute upon unblocking or entry of the DNA feature.
	 */
    public Reaction getEntryReaction(){
        return entryReaction;
    }
	
	/**
	 * Sets the entryReaction as teh reaction that is input.
	 * @param rec the reaction that is to be set as this DNA feature's entryReaction
	 */
    public void setEntryReaction(Reaction rec){
        entryReaction = rec;
    }
    
}

/**
 * Abstract class that describes a Promoter.  This is used to describe promoters that initiate transcription as well
 * as just protein binding sites that do not.  All promoters are parameterized with the constants aon, aoff, ainiton
 * aelongstart, runoffPercent, arecyc for each possible polymerase.
 *
 * @author Sriram Kosuri
 * @version 1.0
 */

abstract class Promoter extends DNA{
    /**
	 * The name of the promoter.  For output purposes.
	 */
	private String name;
    /**
	 * The startsite of transcription for this promoter
	 */
	private int startsite;
	/**
	 * The second-order binding constant describing the binding of protein with this promoters to form an 
	 * initiation complex.  The vector describes values for each RNA polymerase encoded.
	 */
    double[] aon;
    /**
	 * The first-order constant describing the unbinding of the initiation complex to free the promoter. 
	 * The vector describes values for each RNA polymerase encoded.
	 */
	double[] aoff;
	/**
	 * The first-order constant for transition from an intiation complex to an elongation complex. 
	 * The vector describes values for each RNA polymerase encoded.
	 */
    double[] ainiton;
	
	/**
	 * The first-order constant for the elongation complex to transcribe the first 8 bases and maintain a stable
	 * elongation complex.  The vector describes values for each RNA polymerase encoded.
	 */
    double[] aelongstart;
	/**
	 * The percentage of initiation complexes that go on to transcribe.  The other polymerases will recycle back to
	 * the initiation complex, taking arecyc as the time it takes for them to recyle.  The vector describes values for
	 * each RNA polymerase encoded.
	 */
    double[] runoffPercent;
	/**
	 * The first-order constant describing the process of how long the polymerase takes to recycle back to the initiation
	 * complex once it has decided to recycle instead of forming a stable elongation complex.  The vector describes values
	 * for each RNA polymerase encoded.
	 */
	double[] arecyc;
	/**
	 * A variable that tracks the number of initiated polymerases from this promoter.  Used for output purposes.
	 */
    private int initiated;
	
	/**
	 * The only Promoter constructor.
	 * @param start The most upstream position of the promoter.
	 * @param stop The most downstream position of the promoter.
	 * @param ID The promoter's ID.
	 * @param startsite The initiation site for transcription.
	 * @param name The name of the promoter for output purposes.
	 * @param aon Vector containing values for aon for all polymerases.
	 * @param aoff Vector containing values for aoff for all polymerases.
	 * @param ainiton Vector containing values for ainiton for all polymerases.
	 * @param runoffpercent Vector containing values for runoffpercent for all polymerases.
	 * @param aelongstart Vector containing values for aelongstart for all polymerases.
	 * @param arecyc Vector containing values for arecyc for all polymerases.
	 */
    public Promoter(int start, int stop, int ID ,int startsite, String name,double[] aon,double[] aoff,double[] ainiton,double[] runoffpercent,double[] aelongstart,double[] arecyc){
        super(start,stop, ID);
        this.name = name;
        this.startsite = startsite;
        this.aon = aon;
        this.aoff = aoff;
        this.ainiton = ainiton;
        this.runoffPercent = runoffpercent;
        this.aelongstart = aelongstart;
		this.arecyc = arecyc;
		this.initiated = 0;
    }
	
	/**
	 * Resets the initiated variable upon output
	 * @return the current value of initiated.  The number of polymerase initiating transcription since the last time resetInitiated() was invoked.
	 */
	public int resetInitiated(){
		int ret = this.initiated;
		initiated = 0;
		return ret;
	}
	/**
	 * Increments the initiated variable by 1.
	 */
	public void incrementInitiated(){
		initiated += 1;
	}
	
	/**
	 * Returns the name of the promoter
	 * @return the name of the promoter
	 */
    public String getName(){
        return name;
    };

	/**
	 * Returns the Position of the start site of transcription.
	 * @return the position of the start site of transcription.
	 */
    public int getStartSite(){
        return startsite;
    }
	
	/**
	 * Sets the position of the of the start site of transcription.
	 * @param startsite the position of the start site of transcription
	 */
    public void setStartSite(int startsite){
        this.startsite = startsite;
    }

	/**
	 * Overrides the string representation of this object to describe the promoter.  This was used for debugging purposes.
	 * @return A string describing the promoter.
	 */
    public String toString(){
        return ("PROMOTER\tID: "+getID()+"\tname: " + getName() + "\tstart: " + getStart() + "\tstop: " + getStop() + "\tstartsite: " + getStartSite());
    }
}

/**
 * A class for T7 promoters.  This class, while still used, is legacy.  There is no need to distinguish a T7 promoter from any other promoter, since now all 
 * promoters need to be parameterized for all polymerases.
 *
 * @author Sriram Kosuri
 * @version 1.0
 */
class T7Promoter extends Promoter{

    public T7Promoter(int start, int stop, int ID,int startsite, String name, double[] aon, double[] aoff, double[] ainiton, double[] runoffPercent, double[] aelongstart, double[] arecyc){
        super(start,stop, ID, startsite,name,aon,aoff,ainiton,runoffPercent,aelongstart,arecyc);
    }

    public String toString(){
        return ("T7_PROMOTER\tID: "+getID()+"\tname: " + getName() + "\tstart: " + getStart() + "\tstop: " + getStop() + "\tstartsite: " + getStartSite());
    }
        
}

/**
 * A class for E. coli promoters.  This class, while still used, is legacy.  There is no need to distinguish a E. coli promoter from any other promoter, since now all 
 * promoters need to be parameterized for all polymerases.
 *
 * @author Sriram Kosuri
 * @version 1.0
 */
class EcoliPromoter extends Promoter{

    public EcoliPromoter(int start, int stop, int ID, int startsite, String name, double[] aon,double[] aoff,double[] ainiton,double[] runoffPercent,double[] aelongstart,double[] arecyc){
        super(start,stop,ID,startsite,name,aon,aoff,ainiton,runoffPercent,aelongstart,arecyc);
    }

    
    public String toString(){
        return ("ECOLI_PROMOTER\tID: "+getID()+"\tname: " + getName() + "\tstart: " + getStart() + "\tstop: " + getStop() + "\tstartsite: " + getStartSite());
    }
    
}

/**
 * A class that describes the DNA encoding the function of a transcriptional terminator. 
 *
 * @author Sriram Kosuri
 * @version 1.0
 */
class Terminator extends DNA{
	/**
	 * The name of the terminator
	 */
    private String name;
    /**
	 * The terminator efficiency for the polymerases in polID
	 */
	private double[] efficiency;
	/**
	 * The polymerases that are affected by this terminator
	 */
    int[] polID;
	/**
	 * The position on the DNA of transcription termination.
	 */	
    private int stopsite;
    /**
	 * The only constructor for the Terminator class.
	 *
	 * @param start The most upstream position of the terminator.
	 * @param stop The most downstream positino of the terminator
	 * @param ID the terminator's ID
	 * @param name The name of the terminator for output purposes
	 * @param polID A vector that contains integers refering to the ID of polymerases affected by this terminator
	 * @param efficiency A vector that contains an number from 0 (no termination) to 1 (100% termination) describing how efficient the terminator is.  The indices for polID and efficiency are synchronized
	 *        such that the first entry of efficiency refers to the polymerase specified by the first entry of polID.    
	 */
    public Terminator(int start, int stop, int ID,int stopsite, String name, int[] polID, double[] efficiency){
        super(start,stop,ID);
        this.name = name;
        this.efficiency = efficiency;
        this.polID = polID;
        this.stopsite = stopsite;
    }

	/**
	 * returns the Vector of efficiencies for polymerases corresponding to the vector polID.
	 * @return the vector of efficiencies for polymerases corresponding to the vector polID.
	 */
    public double[] getEfficiency(){
        return efficiency;
    }

	/**
	 * Returns the name of the promoter.
	 * @return the name of the promoter.
	 */
    public String getName(){
        return name;
    }
    
	/**
	 * Returns the vector containing the polymeraseID's that are affected by this terminator
	 * @return the vector containing the polymeraseID's that are affected by this terminator
	 */
    public int[] getPolID(){
        return polID;
    }

	/**
	 * Sets the position on the DNA where the polymerase stops transcription.
	 * @param stopsite the position on the DNA where the polymerase stops transcription.
	 */
    public void setStopSite(int stopsite){
        this.stopsite = stopsite;
    }

	/** 
	 * Returns the position on the DNA where the polymerase stops transcription.
	 * @return the position on the DNA where the polymerase stops transcription.
     */
    public int getStopSite(){
        return stopsite;
    }
    
	/**
     * Overrides the string representation of this object to describe the terminator.  This was used for debugging purposes.
	 * @return A string describing the terminator.
	 */
    public String toString(){
        return ("TERMINATOR\tID: "+getID()+"\tname: " + getName() + "\tstart: " + getStart() + "\tstop: " + getStop() + "\tpolID: " + getPolID() + "\tStopSite: " + getStopSite());
    }
    
}

