import java.util.*;
import java.io.*;

/**
 * The class that is a container for the DNA being simulated.  Each instance of this class within a cell has its own 
 * set of genetic elements and mapping to position, set of complexes on the DNA, et cetera.
 *
 * @author Sriram Kosuri
 * @version 1.0
 * @see Cell
 */
public class Phage {
    /** 
	 * A pointer to the cell that this DNA/Phage is contained within.
	 */
	Cell myCell;
	/**
	 * The ID of the Phage.  Not currently used.
	 */
    int ID;
	/**
	 * A 2-D array containing the genetic elements important to transcription and translation. The 0 index of the first dimension contains the promoter elements, while the 1 index of 
	 * same dimension contains the terminators. 
	 */
    DNA[][] myDnaElements;
    /**
	 * Contains alternating elements of ribosome binding sites and termination codons.  The odd indeces are ribosome binding sites , and the even indices are translation stops, starting with index 0.
	 */
	RNA[] myRnaElements;
	/**
	 * Tracks the entry of the genome into the cell.  The number refers to the highest DNA number that has entered the cell.
	 */
    int enteredDNAindex;

	/**
	 * A vector, of length of the DNA molecule, that contains the ID for DNA features that are encoded upon the DNA.  DNA elements thus cannot overlap.
	 */
    short[] myDNAFeatures;
    
	/**
	 * A vector, of length of the DNA molecule, that contains the address to complexes that are encoded upon the DNA.  Complexes on the DNA cannot overlap.
	 */
	DPComplex[] myDNAComplexes;
	/**
	 * A vector, of length of the DNA molecule, that contains the ID for the RNA features that are encoded upon the RNA.  RNA elements thus cannot overlap.
	 */
    short[] myRNAFeatures;
	/**
	 * The length of the DNA molecule that this instance of the Phage class represents.
	 */
    public int GENOME_LENGTH;
	/**
	 * A priority queue containing the Nascent Ribosome binding sites that are associated with this DNA molecules.  
	 */
    PriorityQueue nrbs_PQ;
	/**
	 * The rate at which the DNA enters the cell up to base denoted by the autoEntryOffSite variable.
	 */
    double autoEntryRateConstant;
	/**
	 * The site at which entry of the DNA based upon the 
	 */
    int autoEntryOffSite;
	/**
	 * The number of ribosomes that have elongated beyond the point that the polymerases they are transcribing from have elongated.
	 */
    int stalled_rib;

	/**
	 * The reaction that is called upon release of a polymerase.  The index of the polymerase is co-indexed with the releaseRxn. 
	 */
	DPComplexReleaseRxn[] releaseRxn; 
   
   /**
	* The only constructor for a Phage instance.  
	* @param mycell The Cell that this Phage/DNA is contained within
	* @param ID The ID of the Phage/DNA
	*
	*/
   public Phage(Cell mycell, int ID){
        try{
            this.myCell = mycell;
            enteredDNAindex=0;
            this.ID = ID;
        }catch(Exception e){
            e.printStackTrace();
            System.exit(9);
        }
        stalled_rib = 0;
    }
    
   /**
	* Takes a position on the DNA and returns what DNA feature is encoded upon it, if any.
	* @param featurePosition the position on the DNA 
	* @returns the DNA feature that is encoded.  Returns null if there is no element at that position.
	*/
    public DNA convertFeatureToElement(int featurePosition){
        // retieve the ID of the Feature at given position
        int rid = (int) this.myDNAFeatures[featurePosition];
        int i = ((int) Math.floor(rid/1000)); //-1 , this is a relic from when the it was the 1st dimension in the 2darray.
        int j = (rid%1000) - 1;
        DNA d = null;
        if (rid<=0){
            System.out.println("Should not get here if checking for zeros outside of this method.");
            return null;
           
        }
        
        //there are now possible 6 "types" of promoters(i.e. Ecoli, T7, etc).
        if (i<7){
            d = (DNA) this.myDnaElements[0][j];
        }else{
            //Id 7000 are terminators
            if (i==7){
                d = (DNA) this.myDnaElements[1][j];
            }else{
                //ID = 8000 are RnaseII sites.
                if (i==8){
                    d = (DNA) this.myDnaElements[2][j];
                }
            }
        }
        
        return d;
        
    }
   
    /**
	 * A method that is run at the creation of the variables that track the DNA at single base resolution.
	 */
    protected void createTrackers(){
        myDNAFeatures = new short[GENOME_LENGTH+1];
        myDNAComplexes = new DPComplex[GENOME_LENGTH+1];
        //create myDNAFeatures
        for(int i=0;i<(myDnaElements.length);i++){
            for(int j=0;j<myDnaElements[i].length;j++){
                DNA d = (DNA) myDnaElements[i][j];
                if(d.start>0){
                    for(int k=d.start;k<=d.stop;k++){
                        myDNAFeatures[k] = (short) d.ID;
                    }
                }else{
                    for(int k=d.stop*-1;k<=(d.start*-1);k++){
                        myDNAFeatures[k] = (short) d.ID;
                    }
                }
            }
        }
        //myDNAComplexes is initially empty;
        //myRNAFeatures
        myRNAFeatures = new short[GENOME_LENGTH+1];
        for(int i=0;i<myRnaElements.length;i++){
            RNA r = (RNA) myRnaElements[i];
            if(r.start>0){
                for(int j=r.start;j<=r.stop;j++){
                    myRNAFeatures[j] = (short) r.ID;
                }
            }else{
                for(int j=r.stop*-1;j<=(r.start*-1);j++){
                    myRNAFeatures[j] = (short) r.ID;
                }
            }
        }
    }
    
    
}


