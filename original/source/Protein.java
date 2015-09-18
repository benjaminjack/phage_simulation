/**
 * A class for proteins that are contained in the cell.  
 * @author Sriram Kosuri
 * @version 1.0
 * @see Molecule
 */
public class Protein implements Molecule{
	/**
	 * The name of the protein
	 */
    String name;
	/**
	 * The number of copies of the protein
	 */
    protected int copyNumber;
    /**
	 * The ID of the protein
	 */
	private int ID;
    /**
	 * Depracated.  Used for visualization purposes. 
	 */
	private int organism;
    /**
	 * The reaction that is executed every time this protein is produced.  
	 */
	Reaction gainReaction;
    /**
	 * Only Contructor for protein
	 * @param name The name of the protein; used for output
	 * @copyNumber The number of copies of this protein initially.
	 * @organism A depracated input.  Used for labeling purposes only.
	 */
	
	/**
	 * The length in amino acids for this protein.  Used for energy tracking purposes.
	 */
	int length;
	
    public Protein(String name, int copyNumber, int ID, int organism){
        this.name = name;
        this.copyNumber = copyNumber;
        this.organism = organism;
        this.ID = ID;
		this.length = 0;
    }
	/**
	 * Returns the name of the protein.
	 * @return the name of the protein.
	 */
    public String getName(){
        return name;
    }
	/**
	 * Returns the number of copies of the protein
	 * @return the number of copies of the protein
	 */
    public int getCopyNumber(){
        return copyNumber;
    }
	
	/**
	 * Incremements the copy number of the protein by inc.
	 * @param inc Increments the copy number of protein by inc.  This can be negative or positive.
	 */
    public void incrementCopyNumber(int inc){
        this.copyNumber += inc;
    }

	/**
	 * Returns the ID of the protein.
	 * @return the ID of the protein.  This ID is indexed to the position in the myMolecules vector minus 1.
	 */
    public int getID(){
        return ID;
    }
	
	/**
	 * Returns the organism of the protein. Depracated.
	 * @return the organism ID of the protein.
	 */
    public int getOrganism(){
        return organism;
    }
}

/**
 * A class for proteins that are contained in the cell.  
 * @author Sriram Kosuri
 * @version 1.0
 * @see Protein
 */
class Polymerase extends Protein{

    /**
	 * The reaction that is called when the Polymerase is lost. 
	 */
	private FirstOrderRxn0 lossReaction;
    /**
	 * The rate at which the polymerase traverses the DNA.  Can be set to zero (e.g., a repressor)
	 */
	double speed;
	/**
	 * The number of bases that the polymerase occupies on the DNA.
	 */
	int footprint;
    
	
	/**
	 * The only constructor for the polymerase.  
	 * @param name The name of the polymerase
	 * @param copyNumber The initial copyNumber of the Polymerase.
	 * @param ID The ID of the polymerase.  Also corresponds to the myCell#myMolecules vector.
	 * @param organism The organism of the polymerase.  Used for visualization.
	 * @param speed The rate at which the polymerase traverses the DNA.
	 * @param footprint The number of bases on the DNA that the polymerase covers.  
	 */
    public Polymerase(String name, int copyNumber, int ID, int organism, double speed, int footprint){
        super(name,copyNumber,ID,organism);
        this.speed = speed;
		this.footprint = footprint;
		this.length = 0;
    }
    /**
	 * Set the reaction that is called when a polymerase is decremented.
	 * @param r the Reaction that is called when the polymerase is decremented.
	 */
    public void setLossReaction(FirstOrderRxn0 r){
        lossReaction = r;
    }
	
	/**
     * Returns the reaction that is called when the polymerase is decremented.
	 * @param r the Reaction that is called when the polymerase is decrememented.
	 */
    public Reaction getLossReaction(){
        return lossReaction;
    }
}

/**
 * A class for proteins that are contained in the cell.  
 * @author Sriram Kosuri
 * @version 1.0
 * @see Protein
 */
class Ribosome extends Protein{
    /**
	 * The rate at which the ribosome traverses the mRNA.  bp/s.
	 */
	double speed;
	/**
	 * The number of bases of mRNA that the ribosome covers.
	 */
    int footprint;

	/**
     * The only constructor for the ribosome.  
	 * @param name The name of the ribosome
	 * @param copyNumber The initial copyNumber of the ribosome.
	 * @param ID The ID of the ribosome.  Also corresponds to the myCell#myMolecules vector.
	 * @param organism The organism of the ribosome. Depracated.
	 * @param speed The rate at which the ribosome traverses the RNA (bp/s).
	 * @param footprint The number of bases on the RNA that the ribosome covers.  
	 */
    public Ribosome(String name, int copyNumber,int ID, int organism, double speed,int footprint){
        super(name,copyNumber,ID,organism);
        this.speed = speed;
        this.footprint = footprint;
    }

	/**
		* Incremements the copy number of the ribosome by inc.
	 * @param inc Increments the copy number of ribosome by inc.  This can be negative or positive.
	 */
    public void incrementCopyNumber(int inc){
        this.copyNumber += inc;
    }
}

