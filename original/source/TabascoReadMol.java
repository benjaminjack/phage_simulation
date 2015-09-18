
import java.io.*;
import java.util.StringTokenizer;
import javax.swing.*;
import java.awt.*;
import java.awt.geom.*;
import java.util.*;

/**
 *  TabascoReadMol is the class that is used to read molecule output files of simulations.  
 *
 * @author Jason Kelly
 * @author Sriram Kosuri
 * @version 1.0
 */


public class TabascoReadMol{
    
    double time;
    double abortPercent;
    double runOffPercent;
    int[] polymeraseCopyNumbers; //need an array for polymerases since now will have differnt types
    int ribosomeCopyNumber;
    int[] customProteinCopyNumbers;
	int[] rnaCopyNumbers;
    // these will need to be set in constructor eventaully
    String[] polymeraseNames;
    String[] customProtNames;
	String[] rnaNames;
    int numberOfCustomProt;
    int numberOfPolymerase;
    int numberOfColiPolymerase;
    int numberOfT7Polymerase;
	int num_of_time_points;
	int num_of_proteins;
	int num_of_coding;
	
    int RIBOSOME_TYPE_NUM = 1; //this should always be this
    // Should have variables for the total number of lines in the file (#of timesteps)
    
    BufferedReader brMol;
    
    public TabascoReadMol(String MolFileLoc){
        try{

            
            brMol = new BufferedReader(new FileReader(MolFileLoc));
            
            String firstline = brMol.readLine(); // Read in the first line.
            StringTokenizer firstline_st = new StringTokenizer(firstline,"\t\n\r\f");
            
            String currentToken = firstline_st.nextToken(); //read in time
            currentToken = firstline_st.nextToken(); // get rid of ribosome
            
            Vector tempNames = new Vector();
            int polyCount = 0;
           // String currentToken = "";
            
            while (firstline_st.hasMoreTokens()){
                currentToken = firstline_st.nextToken();
                if (currentToken.startsWith("NumColi")){
                    numberOfColiPolymerase = Integer.parseInt(firstline_st.nextToken());
                }else if (currentToken.startsWith("NumT7")){
					numberOfT7Polymerase = Integer.parseInt(firstline_st.nextToken());
				}else if(currentToken.startsWith("NumOfTimePoints")){
					num_of_time_points = Integer.parseInt(firstline_st.nextToken());
				}else if(currentToken.startsWith("NumOfProteins")){
					num_of_proteins = Integer.parseInt(firstline_st.nextToken());
				}else if(currentToken.startsWith("NumOfCoding")){
					num_of_coding = Integer.parseInt(firstline_st.nextToken());
				}
				else{
					tempNames.addElement(currentToken);
				}
            }
            Iterator tempNamesIter = tempNames.iterator();
            numberOfPolymerase = numberOfColiPolymerase + numberOfT7Polymerase;
            numberOfCustomProt = num_of_proteins - numberOfPolymerase - 1; //for ribosomes
            
            polymeraseCopyNumbers = new int[numberOfPolymerase];
            customProteinCopyNumbers = new int[numberOfCustomProt];
            polymeraseNames = new String[numberOfPolymerase];
            customProtNames = new String[numberOfCustomProt];
            rnaNames = new String[num_of_coding];
			rnaCopyNumbers = new int[num_of_coding];
			/*
			System.out.println("numberOfPolymerase: " + numberOfPolymerase);
			System.out.println("numberOfCustomProt: " + numberOfCustomProt);
			System.out.println("numberOfCoding: " + num_of_coding);
			System.out.println("tempNames.size: " + tempNames.size());
			*/
            for (int i=0;i<(numberOfPolymerase);i++){
                polymeraseNames[i] = (String) tempNamesIter.next();
            }
            for (int i=0;i<numberOfCustomProt;i++){
                customProtNames[i] = (String) tempNamesIter.next();
            }
            for (int i=0;i<num_of_coding;i++){
				rnaNames[i] = (String) tempNamesIter.next();
			}
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(9);
        }
    }
    
	public void close() throws IOException{
		brMol.close();
	}
	
    // This method retrieves the next complex array.
    public void loadNextMoleculeArrayFast(){
        try{
            
            // this method loads the variables which can then be read out with the get methods.
            String mainLine = brMol.readLine();
            StringTokenizer mainLine_st = new StringTokenizer(mainLine,"\t\n\r\f");
			int counter = 0;
            time = Double.parseDouble(mainLine_st.nextToken());
            ribosomeCopyNumber = Math.round(Float.parseFloat(mainLine_st.nextToken()));
			
            for (int i=0;i<numberOfPolymerase;i++){
                polymeraseCopyNumbers[i] = Math.round(Float.parseFloat(mainLine_st.nextToken()));
				counter++;
			}
            
            for (int i=0;i<numberOfCustomProt;i++){
                customProteinCopyNumbers[i] = Math.round(Float.parseFloat(mainLine_st.nextToken()));
            }
            
			for (int i=0;i<num_of_coding;i++){
				rnaCopyNumbers[i] = Math.round(Float.parseFloat(mainLine_st.nextToken()));
			}
			
        }
        catch(Exception e){
            e.printStackTrace();
            // ADD a return for if the tokenizer comes up with an error.
            System.exit(9);
        }
        
        
    }
    
    // GET METHODS
    
    public double getTime(){
        return time;
    }
    /*
    public double getAbortPercent(){
        return abortPercent;
    }
    
    public double getRunOffPercent(){
        return runOffPercent;
    }
    */
    public int[] getPolymeraseCopyNumbers(){
        return polymeraseCopyNumbers;
    }
    
    public int getRibosomeCopyNumber(){
        return ribosomeCopyNumber;
    }
    
    public int[] getCustomProteinCopyNumbers(){
        return customProteinCopyNumbers;
    }
    
    public int getNumberOfColiPolymerase(){
        return numberOfColiPolymerase;
    }
    
    public int getNumberOfT7Polymerase(){
        return numberOfT7Polymerase;
    }
    
	public int getNumberOfRNA(){
		return num_of_coding;
	}
    
	public int[] getRNACopyNumbers(){
		return rnaCopyNumbers;
	}
	
	public String[] getRNANames(){
		return rnaNames;
	}
	
    public String[] getPolymeraseNames(){
        return polymeraseNames;
    }
    
    public String[] getCustomProtNames(){
        return customProtNames;
    }
	
	public int getNumOfTimePoints(){
		return num_of_time_points;
	}
    
    
    // This main is just for testing.
    public static void main(String args[]){
        String OUTPUT_FILE_NAME = "output_sim1_Mol";
        //String INPUT_FILE_LOC = "C:\\Documents and Settings\\Jason\\My Documents\\UROP\\Endy\\Summer 2003\\molecular crowding\\" + INPUT_FILE_NAME;
        String OUTPUT_FILE_LOC = "C:\\Documents and Settings\\Jason\\Desktop\\TabascoCurrent\\output\\" + OUTPUT_FILE_NAME;
        try{
            // TabascoAnim aFrame = new TabascoAnim(myPhage.myDNAFeatures,Complex2DArray);
            
            //BufferedReader bDNA = new BufferedReader(new FileReader(OUTPUT_FILE_LOC+".txt"));
            String fileSpot = OUTPUT_FILE_LOC+".txt";
            
            TabascoReadMol tRead = new TabascoReadMol(fileSpot);
            for (int i=0;i<10;i++){
                tRead.loadNextMoleculeArrayFast();
                System.out.println(tRead.getTime());
                System.out.println(tRead.getPolymeraseCopyNumbers()[1]);
            }
        }
        
        catch(Exception e){
            e.printStackTrace();
            System.exit(9);
        }
        
    }
    
}







