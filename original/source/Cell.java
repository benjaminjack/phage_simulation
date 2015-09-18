import java.io.*;

/**
 * A class representing a cell.
 *
 * The cell class contains the arrays for the Molecules in the cell, the cell volume, and 
 * how the transcribing polymerases interact on the DNA.  Also includes the pointers to
 * the DNA/phage that are entering or currently in the cell.
 *
 * @author Sriram Kosuri
 * @version 1.0
 * @see Molecule
 * @see Phage
 */

public class Cell {
	
	/**
	 * The vector containing all the molecules.  Indexed according to the moleculeID-1.
	 * @see Molecule
	 */
    Molecule[] myMolecules;
	
	/**
	 * A vector that tracks energy usage for the myMolecules vector.  For now, tracks only the amino acids 
	 * used by each protein
	 */
	int moleculeEnergy[];
	
	/**
	 * The volume of the cell.  Used to calculate mesoscopic rate constants from macroscopic rate constants
	 */
	double VOLUME;
    
	/**
	 * The number of phage infecting the cell.
	 * @see Phage
	 */
	int NUM_OF_PHAGE;

	/**
	 * A vector containing all the phage or DNA that are infecting or contained in the cell.
	 * @see Phage
	 */
    Phage myPhages[];
    
	/**
	 * A pointer to the simulation this cell is contained in.
	 */
	TabascoSimulator mySim;
    
	/**
	 * A matrix of integers that describes how each polymerase interacts with one another.  The matrix is N by N
	 * matrix where N is the number of polymerases.  The first index refers to the upstream polymerase index and
	 * the second refers to the downstream polymerase index.  The integer at that value can be either 1, indicating
	 * downstream polymerase will fall off, -1, indicating the upstream polymerase will fall off, or 0 indicating
	 * the downstream polymerase will just follow at the upstream polymerase speed.  
	 */
	int[][] polymeraseInteractionMatrix;
	
	/**
	 * A vector containing the number of transcription events that have occured indexed to each polymerase since
	 * the last time resetTranscriptionEvents() was called.
	 * @see resetTranscriptionEvents()
	 */
    long[] transcriptionEvents;
	
	/**
	 * Creates an empty Cell structure.  In order to fill the contents of this 
	 * class, the TabascoXML.fillCell() method needs to be invoked.
	 *
	 * @param mysim The pointer to the simulator that this cell resides in.
	 * @see TabascoXML#fillCell()
	 */
    public Cell(TabascoSimulator mysim) throws IOException{
        mySim = mysim;
    }
    
	/**
	 * Resets the counter on the number of transcription events that have taken place for each polymerase. 
	 * 
	 * @return a vector containing the number of transcription elongation events that have occured since the method was last invoked 
	 */
	public long[] resetTranscriptionEvents(){
		long[] t = new long[transcriptionEvents.length];
		for(int i=0;i<transcriptionEvents.length;i++){
			t[i] = transcriptionEvents[i];
			transcriptionEvents[i] = 0;
		}
		return t;
	}
	
}
