//
//  t7v3.java
//

//import java.awt.*;
import java.io.*;
import javax.swing.*;
import java.awt.*;
import java.awt.geom.*;
import java.util.*;
import com.sun.image.codec.jpeg.*;
import java.util.zip.Deflater;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

/**
 *  TabascoSimulator is the class used to run simulations and produce output files.  
 * @version 1.0
 * @author Sriram Kosuri
 */
public class TabascoSimulator {
    /**
     * The seed used to initiate the random number generator.  
     */
    long RANDOMSEED;
    /**
     * Pointer to Gillespie integrator instance in this simulation
     */
    GIntegrator gint;
    /**
     * Pointer to the Cell class for this simulation
     */
    Cell myCell;
    /**
     * XML Object to read individual phage XML files
     */
    Vector phageXMOs;
    /**
     * The prefix for the location of the output files.
     */
    String OUTPUT_FILE_LOC;
    /**
     * The name of this simulation.  Currently, this is not output in a meaningful manner.  
     */
    String SIM_NAME;
    /**
     * The ID of the first E. coli polymerase in the molecule array.  For imaging purposes.
     */
    int COLIPOLYMERASEID = -1;
    /**
     * The ID of the first T7 polymerase in the molecule array.  For imaging purposes
     */
    int T7POLYMERASEID = -1;
    /**
     * The ID of the Ribosome in the molecule array
     */
    int RIBOSOMEID = -1;
    /**
     * The ID of the first Protein in the molecule array.
     */
    int FIRSTPROTEINID = -1; // need to figure out how to read these from myCell.myMolecules
    /**
     * The ID of the last Protein in the molecule array.
     */
    int LASTPROTEINID = -1;
    /**
     * The ID of the first polymerase in the molecule array
     */
    int FIRSTPOLYMERASEID = -1;
    /**
     * The ID of the last polymerase in the molecule array 
     */
    int LASTPOLYMERASEID = -1;
    /**
     * The number of E. coli polymerases.  For imaging purposes
     */
    int NUM_OF_COLIPOL = 0;
    /**
     * The number of T7 RNA polymerases. For imaging purposes
     */
    int NUM_OF_T7POL = 0;
    /**
     * The ID of the first Gamma Jump RNA in the molecular array.
     */
    int FIRSTGJRNAID = -1;
    /**
     * For initial protein/RNA production, the number of molecules before the value is recorded in the output file.
     */
    int MIN_NUM_FOR_INIT_OUTPUT = 1;
    /**
     * The time step for outputting the simulation.  In 1/100 sec.
     */
    int timeStep = 200; // in 1/100 sec
    /**
     * The time step for outputting the simulation.  In 1/100 sec. 
     */
    int timeEnd = 30000; // in 1/100 sec
    /**
     * The prefix ID for T7 RNA polymerase promoters.  For imaging purposes. 
     */
    int T7_PRO_ID = 1000;
    /**
     * The prefix ID for E. coli RNA polymerase promoters.  For imaging purposes.
     */
    int COLI_PRO_ID = 2000;
    /**
     * A flag on whether or not to write the DNA position file.
     */
    boolean writeDNA = false;
    /**
     * A flag on whether or not to write the RNA initiation file
     */
    boolean writeInitRNA = false;
    /**
     * A flag on whether or not to write the energy file.
     */
    boolean writeEnergy = false;
    /**
     * A flag on whether or not to include bound Ribosome binding sites in deciding to degrade RNA.
     */
    boolean compDeg = false;
	
    /**
     * The method run from the command line to run Tabasco Simulations.  
     * @param args The arguments passed to the simulator from the command line.  The common usage is the input file, the output file prefix, and an optional random seed.
     */
    public static void main(String args[]) {
		
        int NUM_OF_SIMS;
        TabascoSimulator mySim;
        Reaction[] cellReactions;
        Vector phageReactions;
        int reactionid = 1;
        int moleculeID = -1;
        int totalphage = -1;
        XMLObject simxmo;
	String filenameDNA[];
	TabascoWrite[] tWrite;
		
        try{
            System.out.println("Initializing...");
	    mySim = new TabascoSimulator();
			
	    if(args.length==2||args.length==3){
		simxmo = new XMLObject(new File(args[0]));
		mySim.OUTPUT_FILE_LOC = args[1];
	    }else{
		throw new IOException("usage: java TabascoSimulator input_file output_file_prefix [RANDOM#]");
	    }
	    NUM_OF_SIMS = TabascoXML.fillSim(mySim,simxmo);
            int simCounter=0;
			
	    do{
		mySim = new TabascoSimulator();
		simxmo = new XMLObject(new File(args[0]));
		mySim.OUTPUT_FILE_LOC = args[1];
				
		NUM_OF_SIMS = TabascoXML.fillSim(mySim,simxmo);
		if(args.length==3){
		    mySim.RANDOMSEED = Long.parseLong(args[2]);
		}
		reactionid = 1;
		moleculeID = -1;
		totalphage = -1;
				
                mySim.myCell = new Cell(mySim);
                
                cellReactions = TabascoXML.fillCell(mySim.myCell,simxmo, reactionid);
                reactionid = cellReactions.length+1;
                Vector thephages = new Vector();
                phageReactions = new Vector();
                totalphage = 0;

                //make instances of phage and fill reactions
                for(int i=0;i<mySim.phageXMOs.size();i++) {
                    XMLObject phagexmo = (XMLObject) mySim.phageXMOs.elementAt(i);
                    int tempPhage = Integer.parseInt(phagexmo.getAttribute("copy_number"));
                    for(int j=0;j<tempPhage;j++){
			moleculeID = mySim.myCell.myMolecules.length+1;
                        totalphage++;
                        Phage p = new Phage(mySim.myCell,((i+1)*1000+j));
                        thephages.add(p);
                        TabascoXML.fillPhage(mySim.myCell,p,phagexmo,moleculeID);
                        p.createTrackers();
                        Reaction[] tempreact = TabascoXML.fillPhageReactions(p, reactionid,moleculeID);
                        phageReactions.add(tempreact);
                        reactionid += tempreact.length;
                    }
                }

                //copy phage vector to mycell array
                mySim.myCell.myPhages = new Phage[thephages.size()];
                for(int i=0;i<thephages.size();i++){
                    mySim.myCell.myPhages[i] = (Phage) thephages.elementAt(i);
                }

                //combine reactions into uber-reaction array
                Vector uberReactionVector = new Vector();
                for(int i=0;i<cellReactions.length;i++){
                    uberReactionVector.add(cellReactions[i]);
                }
                for(int i=0;i<phageReactions.size();i++){
                    Reaction[] phagetempreactions = (Reaction[]) phageReactions.elementAt(i);
                    for(int j=0;j<phagetempreactions.length;j++){
                        uberReactionVector.add(phagetempreactions[j]);
                    }
                }

                //turn uberreactionvector to array
                Reaction[] allReactions = new Reaction[uberReactionVector.size()];
                for(int i=0;i<uberReactionVector.size();i++){
                    allReactions[i] = (Reaction) uberReactionVector.elementAt(i);
                }

                //Initialize a's
                for(int i=0;i<allReactions.length;i++){
                    allReactions[i].calculateA();
                }

                //make gintegrator
                mySim.gint = new GIntegrator(allReactions,mySim.RANDOMSEED+simCounter);

                double timeout=0;
				
		tWrite = new TabascoWrite[totalphage];
		
                String filenameMol = mySim.OUTPUT_FILE_LOC + "Mol_sim" + (simCounter+1) + ".txt"; 
		BufferedWriter bwMol = new BufferedWriter(new StringFileWriter(filenameMol));
		
		int numofcoding = mySim.myCell.myPhages[0].myRnaElements.length/2;
		
		if(mySim.writeDNA){
		    filenameDNA = new String[totalphage];
		    for(int i=0;i<totalphage;i++){
			filenameDNA[i] = mySim.OUTPUT_FILE_LOC + "DNA_phage" + (i+1) +"_sim" + (simCounter+1) + ".txt";		      
			tWrite[i] = new TabascoWrite(filenameDNA[i]);
		    }		
		    for(int i=0;i<totalphage;i++){
			tWrite[i].writeFirstLine(mySim.myCell.myPhages[i]);
		    }
		}
		
                bwMol.write(">time\t");

                for(int i=0;i<mySim.LASTPROTEINID;i++){
                    bwMol.write(((Protein) mySim.myCell.myMolecules[i]).name+"\t");
                }
		for(int i=0;i<(numofcoding);i++){
		    MatureRBS tempRBS = (MatureRBS) mySim.myCell.myPhages[0].myRnaElements[i*2+1];
		    bwMol.write("RNA: "+tempRBS.ts.product.getName() + "\t");
		}
		
                bwMol.write("NumColiPoly:\t" + (mySim.NUM_OF_COLIPOL) + "\t" + "NumT7Poly\t" +
                            (mySim.NUM_OF_T7POL) + "\t");
		int timestep = mySim.timeStep;
		int timeend = mySim.timeEnd;
		int num_of_time_points = (timeend/timestep+1); //add one because first time point is at ~zero seconds
		bwMol.write("NumOfTimePoints:\t" + num_of_time_points + "\t");
		bwMol.write("NumOfProteins:\t" + mySim.LASTPROTEINID + "\tNumOfCoding\t" + numofcoding);
                bwMol.write("\n");
                int countLoop = 0;
		
		BufferedWriter bwInitRNA;
		BufferedWriter bwEnergy;
		
		if(mySim.writeInitRNA){
		    String filenameInitRNA = mySim.OUTPUT_FILE_LOC + "Init_Rna_sim"+(simCounter+1)+".txt";
		    bwInitRNA = new BufferedWriter(new StringFileWriter(filenameInitRNA));
		    
		    bwInitRNA.write("time");
		    for(int i=0;i<mySim.myCell.myPhages[0].myDnaElements[0].length;i++){
			bwInitRNA.write("\t"+((Promoter) mySim.myCell.myPhages[0].myDnaElements[0][i]).getName());
		    }
		    for(int i=mySim.FIRSTPOLYMERASEID;i<=mySim.LASTPOLYMERASEID;i++){
			bwInitRNA.write("\t" + ((Polymerase) mySim.myCell.myMolecules[i-1]).getName());
		    }
		    bwInitRNA.write("\n");
		}else{
		    bwInitRNA = null;
		}
                
		if(mySim.writeEnergy){
		    String filenameMolEnergy = mySim.OUTPUT_FILE_LOC + "MolEnergy_sim" + (simCounter+1) + ".txt";
		    bwEnergy = new BufferedWriter(new StringFileWriter(filenameMolEnergy));
		    
		    bwEnergy.write(">time");
		    for(int i=0;i<mySim.LASTPROTEINID;i++){
			bwEnergy.write("\t"+((Protein) mySim.myCell.myMolecules[i]).name);
					    
		    }
		    bwEnergy.write("\n");
		}else{
		    bwEnergy = null;
		}
				
				
		while((mySim.gint.t*100)<mySim.timeEnd){
					
		    mySim.gint.runLoop();
                    countLoop++;
		    if((mySim.gint.t*100)>timeout){
			System.out.println("time = " + mySim.gint.t);
						
			//if it skipped more than 1 timestep b/c no reactions were occurring...
                        while ((mySim.gint.t*100-timeout)>mySim.timeStep) {
			    if(mySim.writeDNA){    
				for(int i=0;i<totalphage;i++){
				    tWrite[i].writeMainLine(mySim.gint.t, mySim.myCell.myPhages[i],i+1);
				}
			    }
                            bwMol.write(mySim.gint.t+"");
                            for(int i=0;i<mySim.LASTPROTEINID;i++){
                                bwMol.write("\t"+mySim.myCell.myMolecules[i].getCopyNumber());
                            }
							
			    for(int i=0;i<numofcoding;i++){
				int MatureRBSCopyNumber = 0;
				int GammaJumpRNACopyNumber = 0;
				for(int j=0;j<mySim.myCell.myPhages.length;j++){
				    MatureRBSCopyNumber += mySim.myCell.myPhages[j].myRnaElements[i*2+1].getCopyNumber();
				    GammaJumpRNACopyNumber += mySim.myCell.myMolecules[(mySim.LASTPROTEINID)+i*2+j*(numofcoding*2)].getCopyNumber();
				}
				int totalRNA = MatureRBSCopyNumber + GammaJumpRNACopyNumber;
				bwMol.write("\t" + totalRNA);
			    }							
							
                            bwMol.write("\n");
                            timeout+=mySim.timeStep;
                        }



                        //write DNA
			if(mySim.writeDNA){
			    for(int i=0;i<totalphage;i++){
				tWrite[i].writeMainLine(mySim.gint.t, mySim.myCell.myPhages[i],i+1);
			    }
			}
						
                        //write molecules
                        bwMol.write(mySim.gint.t+"");
                        for(int i=0;i<mySim.LASTPROTEINID;i++){
                            bwMol.write("\t"+mySim.myCell.myMolecules[i].getCopyNumber());
                        }

                        countLoop=0;
                        timeout+=mySim.timeStep;
											
			//write mrna's for mol file
			for(int i=0;i<numofcoding;i++){
			    int MatureRBSCopyNumber = 0;
			    int GammaJumpRNACopyNumber = 0;
			    for(int j=0;j<mySim.myCell.myPhages.length;j++){
				MatureRBSCopyNumber += mySim.myCell.myPhages[j].myRnaElements[i*2+1].getCopyNumber();
				GammaJumpRNACopyNumber += mySim.myCell.myMolecules[(mySim.LASTPROTEINID)+i*2+j*(numofcoding*2)].getCopyNumber();
			    }
			    int totalRNA = MatureRBSCopyNumber + GammaJumpRNACopyNumber;
							
			    bwMol.write("\t" + totalRNA);
			}
			bwMol.write("\n");
			if(mySim.writeInitRNA){
			    bwInitRNA.write(mySim.gint.t+"");
			    int tempPromInit = 0;
			    for(int i=0;i<mySim.myCell.myPhages[0].myDnaElements[0].length;i++){
				tempPromInit = 0;
				for(int j=0;j<mySim.myCell.myPhages.length;j++){
				    tempPromInit += ((Promoter) (mySim.myCell.myPhages[j].myDnaElements[0][i])).resetInitiated();
				}
				bwInitRNA.write("\t" + tempPromInit);
			    }
			    long[] transEvents = mySim.myCell.resetTranscriptionEvents();
			    for(int i=0;i<transEvents.length;i++){
				bwInitRNA.write("\t" + transEvents[i]);
			    }						
			    bwInitRNA.write("\n");
			}
			if(mySim.writeEnergy){
			    bwEnergy.write(mySim.gint.t+"");
			    for(int i=0;i<mySim.LASTPROTEINID;i++){
				bwEnergy.write("\t"+mySim.myCell.moleculeEnergy[i]);
			    }
			    bwEnergy.write("\n");
			}
						
                    }
                }
                timeout=0;
                System.out.println("Done iteration "+(simCounter+1));
                bwMol.close();
		if(mySim.writeInitRNA)
		    bwInitRNA.close();
				
		if(mySim.writeEnergy)
		    bwEnergy.close();
				
		if(mySim.writeDNA){
		    for(int i=0;i<totalphage;i++){
			tWrite[i].cleanUp();
		    }
		}
		simCounter++;

	    }while(simCounter<NUM_OF_SIMS);
			
	    System.out.println("finished simulation");

	    System.exit(0);
			
        }catch(Exception e){
            e.printStackTrace();
            System.exit(9);
        }
    }
    /**
     * Averages the energy output files from individual simulations. 
     * @param iterations The number of simulations to average
     * @param OUTPUT_FILE_LOC The prefix for the output file location
     * @param STD_FILE_LOC The standard deviation output file location
     * @param AVG_FILE_LOC THe average output file location 
     * @param num_of_time_points The number of timepoints in each simulation
     */
    public static void combineEnergyData(int iterations, String OUTPUT_FILE_LOC, String STD_FILE_LOC, String AVG_FILE_LOC, int num_of_time_points) throws IOException{
	BufferedWriter bw = new BufferedWriter(new StringFileWriter(AVG_FILE_LOC));
        BufferedWriter bw2 = new BufferedWriter(new StringFileWriter(STD_FILE_LOC));
        BufferedReader brmol;
	String filenames[] = new String[iterations];
	for(int p=0;p<iterations;p++){
	    filenames[p] = OUTPUT_FILE_LOC + "MolEnergy_sim" + (p+1) + ".txt";
	}
	//read first line
	brmol = new BufferedReader(new FileReader(filenames[0]));
	String firstline = brmol.readLine();
	StringTokenizer st = new StringTokenizer(firstline,"\t");
	int columns = st.countTokens();
	brmol.close();
		
	//write first line
	bw.write(firstline + "\n");
	bw2.write(firstline + "\n");
	double data [][][] = new double[iterations][num_of_time_points][columns];
	double avg[][] = new double[num_of_time_points][columns];
	double std_dev[][] = new double[num_of_time_points][columns];
		
	for(int i=0;i<iterations;i++){
	    brmol = new BufferedReader(new FileReader(filenames[i]));
	    String currentLine = brmol.readLine();
	    for(int j=0;j<num_of_time_points;j++){
		currentLine = brmol.readLine();
		st = new StringTokenizer(currentLine);
		for(int k=0;k<columns;k++){
		    data[i][j][k] = Double.parseDouble(st.nextToken());
		}
	    }
	    brmol.close();
	}
	//take sum
	for(int i=0;i<iterations;i++){
	    for(int j=0;j<num_of_time_points;j++){
		for(int k=0;k<columns;k++){
		    avg[j][k] += (data[i][j][k]);
		}
	    }
	}
		
	//divide by iterations and write
	for(int j=0;j<num_of_time_points;j++){
	    for(int k=0;k<columns;k++){
		avg[j][k] /= iterations;
		bw.write(avg[j][k]+"\t");
	    }
	    bw.write("\n");
	}
	bw.close();
		
	//take sum of differences
	for(int i=0;i<iterations;i++){
	    for(int j=0;j<num_of_time_points;j++){
		for(int k=0;k<columns;k++){
                    std_dev[j][k] += Math.pow(((data[i][j][k])-avg[j][k]),2);
                }
            }
	}
		
	//divide by iterations, sqrt and write
	for(int j=0;j<num_of_time_points;j++){
	    for(int k=0;k<columns;k++){
		std_dev[j][k] /= iterations;
		std_dev[j][k] = Math.sqrt(std_dev[j][k]);
		bw2.write(std_dev[j][k]+"\t");
	    }
	    bw2.write("\n");
	}
		
	bw2.close();
		
    }
    /**
     * Averages the RNA initiation output files from individual simulations. 
     * @param iterations The number of simulations to average
     * @param OUTPUT_FILE_LOC The prefix for the output file location
     * @param STD_FILE_LOC The standard deviation output file location
     * @param AVG_FILE_LOC THe average output file location 
     * @param num_of_time_points The number of timepoints in each simulation
     */
    public static void combineInitRnaData(int iterations, String OUTPUT_FILE_LOC, String STD_FILE_LOC, String AVG_FILE_LOC,int num_of_time_points) throws IOException{
	BufferedWriter bw = new BufferedWriter(new StringFileWriter(AVG_FILE_LOC));
        BufferedWriter bw2 = new BufferedWriter(new StringFileWriter(STD_FILE_LOC));
        BufferedReader brmol;
	String filenames[] = new String[iterations];
	for(int p=0;p<iterations;p++){
	    filenames[p] = OUTPUT_FILE_LOC + "Init_Rna_sim" + (p+1) + ".txt";
	}

	brmol = new BufferedReader(new FileReader(filenames[0]));
	//read first line
	String firstline = brmol.readLine();
	StringTokenizer st = new StringTokenizer(firstline,"\t");
	int columns = st.countTokens();
	brmol.close();
		
	//write first line
	bw.write(firstline + "\n");
	bw2.write(firstline + "\n");
	double data [][][] = new double[iterations][num_of_time_points][columns];
	double avg[][] = new double[num_of_time_points][columns];
	double std_dev[][] = new double[num_of_time_points][columns];
		
	for(int i=0;i<iterations;i++){
	    brmol = new BufferedReader(new FileReader(filenames[i]));
	    String currentLine = brmol.readLine();
	    for(int j=0;j<num_of_time_points;j++){
		currentLine = brmol.readLine();
		st = new StringTokenizer(currentLine);
		for(int k=0;k<columns;k++){
		    data[i][j][k] = Double.parseDouble(st.nextToken());
		}
	    }
	    brmol.close();
	}
		
	//take sum
	for(int i=0;i<iterations;i++){
	    for(int j=0;j<num_of_time_points;j++){
		for(int k=0;k<columns;k++){
		    avg[j][k] += (data[i][j][k]);
		}
	    }
	}
		
	//divide by iterations and write
	for(int j=0;j<num_of_time_points;j++){
	    for(int k=0;k<columns;k++){
		avg[j][k] /= iterations;
		bw.write(avg[j][k]+"\t");
	    }
	    bw.write("\n");
	}
	bw.close();
		
	//take sum of differences
	for(int i=0;i<iterations;i++){
	    for(int j=0;j<num_of_time_points;j++){
		for(int k=0;k<columns;k++){
                    std_dev[j][k] += Math.pow(((data[i][j][k])-avg[j][k]),2);
                }
            }
	}
		
	//divide by iterations, sqrt and write
	for(int j=0;j<num_of_time_points;j++){
	    for(int k=0;k<columns;k++){
		std_dev[j][k] /= iterations;
		std_dev[j][k] = Math.sqrt(std_dev[j][k]);
		bw2.write(std_dev[j][k]+"\t");
	    }
	    bw2.write("\n");
	}
		
	bw2.close();
		
    }
    /**
     * Averages the molecule file data across multiple simulations
     * @param iterations The number of simulations to average
     * @param OUTPUT_FILE_LOC The prefix for the output file loation
     * @param STD_FILE_LOC The location of the standard deviation output file
     * @param AVG_FILE_LOC The location of the average output file
     */
    public static int combineData(int iterations, String OUTPUT_FILE_LOC, String STD_FILE_LOC, String AVG_FILE_LOC) throws IOException{
	BufferedWriter bw = new BufferedWriter(new StringFileWriter(AVG_FILE_LOC));
        BufferedWriter bw2 = new BufferedWriter(new StringFileWriter(STD_FILE_LOC));
        TabascoReadMol[] brmol = new TabascoReadMol[iterations];

	for(int p=0;p<iterations;p++){
            brmol[p] = new TabascoReadMol(OUTPUT_FILE_LOC + "Mol_sim" + (p+1) + ".txt");
        }
		
	String polnames[] = brmol[0].getPolymeraseNames();
	String custprotnames[] = brmol[0].getCustomProtNames();
	String rnanames[] = brmol[0].getRNANames();
	//need to write RNA's
		
	int num_of_time_points = brmol[0].getNumOfTimePoints();
	int num_of_pol = polnames.length;
	int num_of_custprot = custprotnames.length;
	int num_of_prot = num_of_pol + num_of_custprot + 1; //ribosomes
	int num_of_rna = rnanames.length;
		
	bw.write("time\tRibosomes\t");
	bw2.write("time\tRibosomes\t");
		
	for(int i=0;i<num_of_pol;i++){
	    bw.write(polnames[i] + "\t");
	    bw2.write(polnames[i] + "\t");
	}
		
	for(int i=0;i<num_of_custprot;i++){
	    bw.write(custprotnames[i] + "\t");
	    bw2.write(custprotnames[i] + "\t");
	}
		
	for(int i=0;i<num_of_rna;i++){
	    bw.write(rnanames[i] + "\t");
	    bw2.write(rnanames[i] + "\t");
	}

	bw.write("NumColiPoly:\t" + (brmol[0].getNumberOfColiPolymerase()) + "\t" + "NumT7Poly\t" + (brmol[0].getNumberOfT7Polymerase()) + "\t");
	bw.write("NumOfTimePoints:\t" + num_of_time_points + "\t");
	bw.write("NumOfProteins:\t" + num_of_prot + "\tNumOfCoding\t" + num_of_rna);

	bw2.write("NumColiPoly:\t" + (brmol[0].getNumberOfColiPolymerase()) + "\t" + "NumT7Poly\t" + (brmol[0].getNumberOfT7Polymerase()) + "\t");
	bw2.write("NumOfTimePoints:\t" + num_of_time_points + "\t");
	bw2.write("NumOfProteins:\t" + num_of_prot + "\tNumOfCoding\t" + num_of_rna);
		
		
		
	bw.write("\n");
	bw2.write("\n");
		
	int columns = custprotnames.length + polnames.length + num_of_rna + 2; //the two is for time & ribosomes
	for(int iterate_time=0;iterate_time<num_of_time_points;iterate_time++){
	    double data[][] = new double[iterations][columns];
	    double avg[] = new double[columns];
	    double std_dev[] = new double[columns];			
	    for(int i=0;i<iterations;i++){
		brmol[i].loadNextMoleculeArrayFast();
	    }
	    for(int i=0;i<iterations;i++){
		int species_count = 0;
		data[i][species_count] = brmol[i].getTime();species_count++;
		data[i][species_count] = brmol[i].getRibosomeCopyNumber();species_count++;
		int[] polcopies = brmol[i].getPolymeraseCopyNumbers();
		int[] custprotcopies = brmol[i].getCustomProteinCopyNumbers();
		int[] rnacopies = brmol[i].getRNACopyNumbers();
				
		for(int p=0;p<num_of_pol;p++){
		    data[i][species_count] = (double) polcopies[p];species_count++;
		}
		for(int p=0;p<num_of_custprot;p++){
		    data[i][species_count] = (double) custprotcopies[p];species_count++;
		}
		for(int p=0;p<num_of_rna;p++){
		    data[i][species_count] = (double) rnacopies[p];species_count++;
		}
	    }
			
	    //take sum
            for(int j=0;j<iterations;j++){
                for(int k=0;k<columns;k++){
                    avg[k] += (data[j][k]);
                }
            }
            
            //divide by iterations and write
            for(int m=0;m<columns;m++){
                avg[m] /= iterations;
                bw.write(avg[m]+"\t");
            }
         
	    //take sum of differences
            for(int j=0;j<iterations;j++){
                for(int k=0;k<columns;k++){
                    std_dev[k] += Math.pow(((data[j][k])-avg[k]),2);
                }
            }
            //divide by iterations, sqrt and write
            for(int m=0;m<columns;m++){
                std_dev[m] /= iterations;
                std_dev[m] = Math.sqrt(std_dev[m]);
                bw2.write(std_dev[m]+"\t");
            }
			
	    bw.write("\n");
            bw2.write("\n");

	}
	for(int i=0;i<iterations;i++){
	    brmol[i].close();
	}
		
		
        bw.close();
        bw2.close();
	return num_of_time_points;
    }

	
    /**
     * A class that was used for testing of the gamma distribution.
     */
    private static void testGamma(){
        IndexedPriorityQueue ipq = new IndexedPriorityQueue(1,15084735);
        for(int i=0;i<5000;i++){
            System.out.println(ipq.getGammaTime(50,50));
	    double random = ipq.getRandomDouble();
	    double time = (1)*Math.log(1/random);
	    System.out.print("\t" + time + "\n");
        }
    }
}


/**
 *  An extension of FileWriter that allows writing strings. 
 * @version 1.0
 * @author Sriram Kosuri
 */

class StringFileWriter extends FileWriter{
    /**
     * Only constructor for StringFileWriter
     * @param fileName name of file to write.  
     */
    public StringFileWriter(String fileName) throws IOException{
        super(fileName);
    }
    /**
     * Writes a string to the file
     * @param s String to write to file
     */
    public void write(String s) throws IOException{
        super.write(s,0,s.length());
    }
    
}

