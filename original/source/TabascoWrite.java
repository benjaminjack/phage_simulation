import java.io.*;
import java.util.StringTokenizer;

/**
 *  TabascoWrite is the class that is used to write DNA output files of simulations.  
 *
 * @author Jason Kelly
 * @author Sriram Kosuri
 * @version 1.0
 */

public class TabascoWrite{
    
    //Variables
    private String fileName;
    private BufferedWriter bwDNA;
    
    //CONSTRUCTOR
    public TabascoWrite(String filename) throws IOException{
        fileName = filename;
        bwDNA = new BufferedWriter(new StringFileWriter(fileName));
        
    }
    
    //METHODS
    
    /* Writes the first line which may contain general information which is not changed at
     * each new time point.
     */
    
    public void writeFirstLine(Phage thephage) throws IOException{
        int timestep = thephage.myCell.mySim.timeStep;
        int timeend = thephage.myCell.mySim.timeEnd;
        short[] dnafeatures = thephage.myDNAFeatures;
        int numofcoding = thephage.myRnaElements.length/2;
       
        // ADD write timeSTEp and then finsihs fixing read
        bwDNA.write("%GeneralInfo\t");
        writeNumberOfTimepoints(timeend/timestep+1); //add one because first time point is at ~zero seconds
        writeGenomeSize(dnafeatures.length-1);
        writeTimeStep(timestep);
        writeTimeEnd(timeend);
        writeNumberOfCodingRgns(numofcoding);
        writeDNAFeatures(dnafeatures);
        
        bwDNA.write("\n");
        
    }
    
    /* Writes one main line.  Writes the results of one timestep to the output file.
     *
     */
    
    public void writeMainLine(double time, Phage thephage, int whichphage) throws IOException{
        
        int entry = thephage.enteredDNAindex;
        DPComplex[] complexes = thephage.myDNAComplexes;
        Molecule[] molecules = thephage.myCell.myMolecules;
        int lastprotID = thephage.myCell.mySim.LASTPROTEINID;
        int numofcoding = thephage.myRnaElements.length/2;
        
        
        writeTime(time);
        writeEntry(entry);
        writeComplexes(complexes);
        writeMRNA(molecules, thephage, lastprotID, whichphage, numofcoding);
        bwDNA.write("\n");
    }
    
    /* Individual writers, each of these methods writes one variable or set of variables
     * to the current bwDNA.
     */
    public void writeNumberOfCodingRgns(int numOfCodingRgn) throws IOException{
        bwDNA.write("%NumberOfCodingRgns\t" + numOfCodingRgn +  "\t");
    }
    
    public void writeTimeStep(int timeStep) throws IOException{
        bwDNA.write("%TimeStep\t" + timeStep +  "\t");
    }
    
    public void writeTimeEnd(int timeEnd) throws IOException{
        bwDNA.write("%TimeEnd\t" + timeEnd +  "\t");
    }
    
    
    public void writeGenomeSize(int GenomeSize) throws IOException{
        bwDNA.write("%GenomeSize\t" + GenomeSize +  "\t");
    }
    
    
    public void writeNumberOfTimepoints(int numOfTimePoints) throws IOException{
        bwDNA.write("%NumberOfTimePoints\t" + numOfTimePoints +  "\t");
    }
    
    public void writeTime(double time) throws IOException{
        bwDNA.write("%Time\t" + time + "\t");
    }
    
    public void writeEntry(int entry) throws IOException{
        bwDNA.write("%Entry\t" + entry + "\t");
    }
    
    public void writeComplexes(DPComplex[] DNAComplexArray) throws IOException{
        
        bwDNA.write("%Complexes\t");
        for (int i=0;i<(DNAComplexArray.length-1); i++){
            String writeout;
            if(DNAComplexArray[i]==null){
                writeout = "0-0-0";
            }else{
                writeout = DNAComplexArray[i].ID;
            }
            bwDNA.write(writeout + "\t");
        }
    }
    
    public void writeDNAFeatures(short[] DNAFeatures) throws IOException{
        
        bwDNA.write("%DNAFeatures\t");
        for (int i=0;i<(DNAFeatures.length-1); i++){
            bwDNA.write(DNAFeatures[i] + "\t");
        }
        
    }
    
    public void writeMRNA(Molecule[] MoleculesArray, Phage thephage , int lastProtID, int whichPhage, int numOfCodingRgn) throws IOException{
        
        // Right now there is an intitiation and an elongation gammaJumpRNA,
        // with a different binding schema this would be different
        int NUM_OF_GJRNA_FOR_EACH_CODING_RGN = 2;
        
        //LastProtID shows the start of the GammaJumpRNA's
        int firstGJRNAid = lastProtID+1+(whichPhage-1)*NUM_OF_GJRNA_FOR_EACH_CODING_RGN ;
        
        bwDNA.write("%mRNA\t");
        
        for (int i=0;i<numOfCodingRgn;i++){

            MatureRBS tempRBS = (MatureRBS) thephage.myRnaElements[i*2+1];
			int NascentRBScopyNumber = tempRBS.nrbs.myRBSInstances.size();

            bwDNA.write(tempRBS.startsite + "\t");
            bwDNA.write(tempRBS.elongstepsize + "\t");
            bwDNA.write(NascentRBScopyNumber + "\t");
        
		}
        
    }
    
    public void cleanUp() throws IOException{
        bwDNA.close();
    }
    
}