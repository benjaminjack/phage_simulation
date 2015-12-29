
/**
 * An interface that is used in order to track various molecules such as proteins, rna species, and ribosomes.  This is
 * to allow them a common interface and to be all stored in the Cell#myMolecules vector.  Essentially, anything that
 * is tracked at the species level is treated as a Molecule.
 * @author Sriram Kosuri
 * @version 1.0
 */
interface Molecule {
	/**
	 * Increments the copy number of the molecule by inc.
	 * @param inc The amount to increment the molecule by.  This number can be negative.
	 */
    void incrementCopyNumber(int inc);
	/**
	 * 
	 * @return the number of copies of the molecule.
	 */
    int getCopyNumber();
	/**
	 * Returns the ID of the molecule.
	 * @return the ID of the molecule
	 */
    int getID();
}

/**
 * The interface is used for reactions that follow Gamma distributed kinetics rather than exponential kinetics 
 */
interface GammaJumpRxn {
	/**
	 *  Increments the copy number of the molecule by inc.
	 *  @param inc The number of copies to increment the current copy number by.
	 */
    public void incremented(int inc);
}