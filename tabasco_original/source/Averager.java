import java.io.*;

/**
 * Averages existing output files
 *
 * @author Sriram Kosuri
 * @version 1.0
 * @see TabascoSimulator
 */

public class Averager {
/**
 * Run from command line to average a set of output files.  Outputs to text files in the directory the program is run.
 * @param args The command line input specifying the number of iterations, and the output file prefix.  
 * Assumes after the prefix are Mol_sim#.txt or RNA_sim#.txt.  
 */
    public static void main(String args[]) {
		try{
			if(args.length!=2){
				System.out.println("Usage: java Averager #_of_iterations output_file_prefix");
				System.exit(0);
			}
			String STD_FILE_LOC = "output_STDDEV.txt";
			String AVG_FILE_LOC = "output_AVG.txt";
			String STD_RNAINITFILE_LOC = "rna_init_STDDEV.txt";
			String AVG_RNAINITFILE_LOC = "rna_init_AVG.txt";
			String AVG_ENERGY_LOC = "prot_energy_AVG.txt";
			String STD_ENERGY_LOC = "prot_energy_STDDEV.txt";
			int iterations = Integer.parseInt(args[0]);
			String OUTPUT_FILE_LOC = args[1];
			int num_of_time_points = TabascoSimulator.combineData(iterations, OUTPUT_FILE_LOC, STD_FILE_LOC, AVG_FILE_LOC);
			System.out.println("Averaged Molecule File");
			try{
			    TabascoSimulator.combineInitRnaData(iterations, OUTPUT_FILE_LOC, STD_RNAINITFILE_LOC, AVG_RNAINITFILE_LOC, num_of_time_points);
			    System.out.println("Averaged RNA Initiation Files");
			}catch(FileNotFoundException e){
			    System.out.println("No RNA Initiation Files found");
			}
			try{
			    TabascoSimulator.combineEnergyData(iterations, OUTPUT_FILE_LOC, STD_ENERGY_LOC, AVG_ENERGY_LOC, num_of_time_points);
			    System.out.println("Averaged Energy Files");
			}catch(FileNotFoundException e){
			    System.out.println("No Energy Files Found");
			}
			System.out.println("Done!");
		}catch(Exception e){
			e.printStackTrace();
			System.exit(9);
		}
	}
    
}
