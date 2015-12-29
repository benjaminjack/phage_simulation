
import java.io.*;
import java.util.StringTokenizer;
import javax.swing.*;
import java.awt.*;
import java.awt.geom.*;

/**
 *  TabascoRead is the class that is used to read DNA output files of simulations.  
 *
 * @author Jason Kelly
 * @author Sriram Kosuri
 * @version 1.0
 */

public class TabascoRead{
    int[] dnaArray;
    int timeStep;
    int timeEnd;
    int numberOfTimePoints;
    int numberOfCodingRgns;
    int genomeSize;
    String brFileLoc;
    
    double currentTimeF;
    int enteredDNAF;
    int[] mRNAStartF;
    int[] mRNALengthF;
    int[] mRNACopyNumberF; 
    String[] complexArrayF;
    
    double currentTime;
    BufferedReader brDNA;
    
    
    public TabascoRead(String filename) throws IOException{
        
        brFileLoc = filename;
        brDNA = new BufferedReader(new FileReader(filename));
        
        String firstline = brDNA.readLine(); // Read in the first line.
        StringTokenizer firstline_st = new StringTokenizer(firstline,"\t\n\r\f");
        String currentToken;
        //String currentToken = firstline_st.nextToken(); 
        
        while (firstline_st.hasMoreTokens()){
            
            currentToken = firstline_st.nextToken();
            
            if (currentToken.startsWith("%NumberOfTimePoints")){
                numberOfTimePoints = Integer.parseInt(firstline_st.nextToken());
            }else{
                if (currentToken.startsWith("%GenomeSize")){
                    genomeSize = Integer.parseInt(firstline_st.nextToken());
                }else{
                    if (currentToken.startsWith("%TimeStep")){
                        timeStep = Integer.parseInt(firstline_st.nextToken());
                    }else{
                        if(currentToken.startsWith("%TimeEnd")){
                            timeEnd = Integer.parseInt(firstline_st.nextToken());
                        }else{
                            if(currentToken.startsWith("%NumberOfCodingRgns")){
                                numberOfCodingRgns = Integer.parseInt(firstline_st.nextToken());
                            }else{
                                
                                if(currentToken.startsWith("%DNAFeatures")){
                                    dnaArray = new int[genomeSize];
                                    for (int iDNA=0;iDNA<genomeSize;iDNA++){
                                        dnaArray[iDNA] = Integer.parseInt(firstline_st.nextToken());
                                    }
                                }
                            }
                        }
                    }
                }
            }
            
            
            
        } //while loop
        
        //complexArray = new int[genomeSize];  // initialize the complex array.
        complexArrayF = new String[genomeSize]; // initialize the complex array.

        mRNAStartF = new int[numberOfCodingRgns];
        mRNALengthF = new int[numberOfCodingRgns];
        mRNACopyNumberF = new int[numberOfCodingRgns];
    }
    
    // FAST Methods
    
    public void loadNextDNAArrayFast(){
        try{
            
            String currentToken;
            String fastLine = brDNA.readLine();
            StringTokenizer fastLine_st = new StringTokenizer(fastLine,"\t\n\r\f");
            
            while (fastLine_st.hasMoreTokens()){
                
                currentToken = fastLine_st.nextToken();
                if (currentToken.startsWith("%Time")){
                    currentTimeF = Double.parseDouble(fastLine_st.nextToken());
                }else{
                    if (currentToken.startsWith("%Entry")){
                        enteredDNAF = Integer.parseInt(fastLine_st.nextToken());
                    }else{
                        if (currentToken.startsWith("%Complexes")){
                            for (int iDNA=0; iDNA<genomeSize;iDNA++){
//                                complexArrayF[iDNA] = Integer.parseInt(fastLine_st.nextToken());
                                complexArrayF[iDNA] = fastLine_st.nextToken();
                            }
                        }else{
                            
                            if (currentToken.startsWith("%mRNA")){
                                for (int iRNA=0; iRNA<numberOfCodingRgns;iRNA++){
                                    //System.out.println(iRNA);
                                   // System.out.println(fastLine_st.nextToken());
                                    mRNAStartF[iRNA] = Integer.parseInt(fastLine_st.nextToken());
                                    mRNALengthF[iRNA] = Integer.parseInt(fastLine_st.nextToken());
                                    mRNACopyNumberF[iRNA] = Integer.parseInt(fastLine_st.nextToken());
                                    // you will have to pass the codind rgn # to the get methods for the RNA arrays.
                                }
                                
                            }
                        }
                    }
                    
                }
                  
            }
        }catch(Exception e){
            e.printStackTrace();
            // ADD a return for if the tokenizer comes up with an error.
            System.exit(9);
        }
    }
	
	//clear mem Fast method
	
	public void clearMemFast(){
		int[] mRNAStartF = null;
		int[] mRNALengthF = null;
		int[] mRNACopyNumberF = null; 
		String[] complexArrayF = null;
		int[] dnaArray = null;
	}
        // get methods for FAST
        public double getTimeFast(){
            return currentTimeF;
        }
        
        public int getEnteredDNAFast(){
            return enteredDNAF;
        }
       
        public String[] getComplexArrayFast(){
            return complexArrayF;
        }
            
        public int[] getMRNAStartFast(){
            return mRNAStartF;
        }
        public int[] getMRNALengthFast(){
            return mRNALengthF;
        }
        public int[] getMRNACopyNumberFast(){
            return mRNACopyNumberF;
        }
        
        
        public double getTime(int stepNumber) throws IOException{
            
            StringTokenizer mainLine_st = returnTokenRow(stepNumber);
            String currentToken = mainLine_st.nextToken();
            
            
            while (mainLine_st.hasMoreTokens()){
                
                currentToken = mainLine_st.nextToken();
                if (currentToken.startsWith("%Time")){
                    currentTime = Double.parseDouble(mainLine_st.nextToken());
                    return currentTime;
                }
                
            }
            return 9999999; //this is equivalent of no time being listed in the output file.
        }
        
        
        public int getEnteredDNA(int stepNumber) throws IOException{
            int enteredDNA;
            
            StringTokenizer mainLine_st = returnTokenRow(stepNumber);
            String currentToken = mainLine_st.nextToken();
            
            
            while (mainLine_st.hasMoreTokens()){
                currentToken = mainLine_st.nextToken();
                if (currentToken.startsWith("%Entry")){
                    enteredDNA = Integer.parseInt(mainLine_st.nextToken());
                    return enteredDNA;
                }
                
            }
            return 999999;
        }
        
        
        public int getMRNAStart(int stepNumber, int codingRgn) throws IOException{
            
            StringTokenizer mainLine_st = returnTokenRow(stepNumber);
            String currentToken = mainLine_st.nextToken();
            
            int[] mRNAStart = new int[numberOfCodingRgns];
            
            
            while (mainLine_st.hasMoreTokens()){
                currentToken = mainLine_st.nextToken();
                if (currentToken.startsWith("%mRNA")){
                    for (int iRNA=0; iRNA<numberOfCodingRgns;iRNA++){
                        mRNAStart[iRNA] = Integer.parseInt(mainLine_st.nextToken());
                        mainLine_st.nextToken();
                        mainLine_st.nextToken();
                    }
                    return mRNAStart[codingRgn-1]; //return the mRNA start position from the coding region you care about.
                }
                
            }
            return 999999;
        }
        
        public int getMRNALength(int stepNumber, int codingRgn) throws IOException{
            
            StringTokenizer mainLine_st = returnTokenRow(stepNumber);
            String currentToken = mainLine_st.nextToken();
            
            int[] mRNALength = new int[numberOfCodingRgns];
            
            
            while (mainLine_st.hasMoreTokens()){
                currentToken = mainLine_st.nextToken();
                if (currentToken.startsWith("%mRNA")){
                    for (int iRNA=0; iRNA<numberOfCodingRgns;iRNA++){
                        mainLine_st.nextToken();
                        mRNALength[iRNA] = Integer.parseInt(mainLine_st.nextToken());
                        mainLine_st.nextToken();
                    }
                    return mRNALength[codingRgn-1]; //return the mRNA start position from the coding region you care about.
                }
                
            }
            return 9999999;
        }
        
        public int getMRNACopyNumber(int stepNumber, int codingRgn) throws IOException{
            
            StringTokenizer mainLine_st = returnTokenRow(stepNumber);
            String currentToken = mainLine_st.nextToken();
            
            int[] mRNACopyNumber = new int[numberOfCodingRgns];

            while (mainLine_st.hasMoreTokens()){
                currentToken = mainLine_st.nextToken();
                if (currentToken.startsWith("%mRNA")){
                    for (int iRNA=0; iRNA<numberOfCodingRgns;iRNA++){
                        mainLine_st.nextToken();
                        mainLine_st.nextToken();
                        mRNACopyNumber[iRNA] = Integer.parseInt(mainLine_st.nextToken());
                    }
                    return mRNACopyNumber[codingRgn-1]; //return the mRNA start position from the coding region you care about.
                }
                
            }
            return 999999;
        }
        
    /*    
        public int[] getComplexArray(int stepNumber) throws IOException{
            
			StringTokenizer mainLine_st = returnTokenRow(stepNumber);
            String currentToken = mainLine_st.nextToken();
            
            while (mainLine_st.hasMoreTokens()){
                currentToken = mainLine_st.nextToken();
                if (currentToken.startsWith("%Complexes")){
                    for (int iDNA=0; iDNA<genomeSize;iDNA++){
                        complexArray[iDNA] = Integer.parseInt(mainLine_st.nextToken());
                    }
                    return complexArray;
                }
                
                
            }
            return null;
        }
      */
	
        // GET METHODS
        
        public int getTimeStep(){
            return timeStep;
        }
        
        public int getNumberOfCodingRgns(){
            return numberOfCodingRgns;
        }
        
        public int getTimeEnd(){
            return timeEnd;
        }
        
        public int getTimePoints(){
            return numberOfTimePoints;
        }
        
        public int getGenomeSize(){
            return genomeSize;
        }
        
        public int[] getDNAArray(){
            return dnaArray;
        }
        
        //UTILITIES
        
        // find the given row and return it as a string tokenizer.
        public StringTokenizer returnTokenRow(int rowNum) throws IOException{
            //BufferedReader brDNA = originalBr;
            
            brDNA = resetBr(brDNA,brFileLoc);
            
            //move to the line corresponding to the given step number.
            for(int i=0;i<rowNum;i++){
                brDNA.readLine();
            }
            
            String mainLine = brDNA.readLine();
            StringTokenizer newTokenizer = new StringTokenizer(mainLine,"\t\n\r\f");
            
            return newTokenizer;
        }
        
        public BufferedReader resetBr(BufferedReader brToBeReset, String filename)throws IOException{
            
            brToBeReset.close();
            BufferedReader newBr = new BufferedReader(new FileReader(filename));
            return newBr;
        }
        
        public void cleanUp(){
            try{
                brDNA.close();
            }catch(Exception e){
                e.printStackTrace();
                System.exit(9);
            }
            
        }
        
        public boolean IsThereNextComplexArray(){
            //Is there a  use for this method?  Maybe if you didn't want to deal with making a loop based on the time points could just
            // go to end of complex arrays instead...
            return true;
        }
        
    }
    
    
    
    
    
    
    
    
