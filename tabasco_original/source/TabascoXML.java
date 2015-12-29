import java.util.*;
import java.io.*;
/**
 * TabascoXML parses the XML input file that specifies simulation parameters
 * @version 1.0
 * @author Sriram Kosuri
 */
public class TabascoXML {
	/**
	 * Parses parameters to initialize the Simulator class
	 * @param mySim A pointer to the simulation class to initialize
	 * @param simxmo The XML object from which to parse information
	 */
    public static int fillSim(TabascoSimulator mySim, XMLObject simxmo) throws IOException{
        XMLObject simparamxmo = simxmo.getChild("EXECUTE-SIMULATION");
        int NUM_OF_SIMS = Integer.parseInt(simparamxmo.getAttribute("runs"));
        mySim.RANDOMSEED = Long.parseLong(simparamxmo.getAttribute("random_seed"));
        mySim.SIM_NAME = simparamxmo.getAttribute("simulation_name");
        mySim.timeStep = Integer.parseInt(simparamxmo.getAttribute("time_step"))*100;
        mySim.timeEnd = Integer.parseInt(simparamxmo.getAttribute("time_end"))*100;
		String wd = simparamxmo.getAttribute("write_DNA");
		String wir = simparamxmo.getAttribute("write_InitRNA");
		String we = simparamxmo.getAttribute("write_Energy");
		String compDeg = simparamxmo.getAttribute("compDeg");
				
		if(wd.toLowerCase().equals("yes")){
			mySim.writeDNA = true;
		}
		if(wir.toLowerCase().equals("yes")){
			mySim.writeInitRNA = true;
		}
		if(we.toLowerCase().equals("yes")){
			mySim.writeEnergy = true;
		}
		if(compDeg.toLowerCase().equals("yes")){
			mySim.compDeg = true;
		}
		
		return NUM_OF_SIMS;
    }

	/** 
 	 * Parses parameters to initialize the Cell class
	 * @param myCell A pointer to the cell class to initialize
	 * @param simxmo The XML object from which to parse information
	 * @param reactionid The reactionID to begin adding reactions
	 * @return A vector containing reactions to add to the simulation class.
	 */	
    public static Reaction[] fillCell(Cell myCell, XMLObject simxmo, int reactionid) throws IOException{
        Vector reactions = new Vector();
        Vector molecule = new Vector();
        int moleculeID = 1;
        final double AVOGADRO = 6.022e23;
        int T7_ORGANISM_ID = 1;

        XMLObject simparamxmo = simxmo.getChild("EXECUTE-SIMULATION");
        XMLObject cellxmo =(XMLObject) (simparamxmo.getNamedChildren("CELL")).firstElement();

        // Set the volume and growth rate of the cell
        myCell.VOLUME = Double.parseDouble(cellxmo.getAttribute("volume"));

        //Create ribosomes
        try{
            XMLObject ribosomexmo =(XMLObject) (cellxmo.getNamedChildren("RIBOSOME")).firstElement();
            int init = Integer.parseInt(ribosomexmo.getAttribute("initial_number"));
            double speed = Double.parseDouble(ribosomexmo.getAttribute("speed"));
            int footprint = Integer.parseInt(ribosomexmo.getAttribute("footprint"));
            Ribosome t = new Ribosome("Ribosome",init,moleculeID,2,speed,footprint);moleculeID++;
            molecule.add(t);
            myCell.mySim.RIBOSOMEID = moleculeID-1;
        }catch(Exception e){
            e.printStackTrace();
            System.out.println("error in creating ribosomes: check input file");
            System.exit(0);
        }

        myCell.mySim.FIRSTPOLYMERASEID = moleculeID;
        //create polymerases

        for(Iterator i = cellxmo.getNamedChildren("POLYMERASE").iterator(); i.hasNext();) {
            XMLObject polymerasexmo = (XMLObject) i.next();
            int init = Integer.parseInt(polymerasexmo.getAttribute("n"));
            double speed = Double.parseDouble(polymerasexmo.getAttribute("speed"));
            int ID = Integer.parseInt(polymerasexmo.getAttribute("id"));
            int organism_id = Integer.parseInt(polymerasexmo.getAttribute("organism_id"));
			int footprint = Integer.parseInt(polymerasexmo.getAttribute("footprint"));
            if(organism_id==1){
                myCell.mySim.NUM_OF_T7POL++;
            }else{
                myCell.mySim.NUM_OF_COLIPOL++;
            }
            String name = polymerasexmo.getAttribute("name");
            if(ID!=moleculeID){
                throw new IOException("Error in polymerase ID scheme: must start with ID 2 and progress incrementally");
            }
            
            Polymerase p = new Polymerase(name,init,moleculeID,organism_id,speed,footprint);moleculeID++;
            molecule.add(p);

        }
        myCell.mySim.LASTPOLYMERASEID = moleculeID-1;
		myCell.transcriptionEvents = new long[myCell.mySim.LASTPOLYMERASEID - myCell.mySim.FIRSTPOLYMERASEID + 1];
        //get phage
        myCell.mySim.phageXMOs = cellxmo.getNamedChildren("DNA_SYSTEM");

        //get polymerase interaction model
        String interaction_type = cellxmo.getAttribute("polymerase_interaction_model");
        int[][] polintermatrix;
        
        if(interaction_type.equals("DOWNSTREAM_FALL_OFF")){
            int matrixdim = myCell.mySim.LASTPOLYMERASEID - myCell.mySim.FIRSTPOLYMERASEID + 1;
            polintermatrix = new int[matrixdim][matrixdim];
            for(int i=0;i<polintermatrix.length;i++){
                for(int j=0;j<polintermatrix[0].length;j++){
                    polintermatrix[i][j] = -1;
                }
            }
        }else if(interaction_type.equals("UPSTREAM_FALL_OFF")){
            int matrixdim = myCell.mySim.LASTPOLYMERASEID - myCell.mySim.FIRSTPOLYMERASEID + 1;
            polintermatrix = new int[matrixdim][matrixdim];
            for(int i=0;i<polintermatrix.length;i++){
                for(int j=0;j<polintermatrix[0].length;j++){
                    polintermatrix[i][j] = 1;
                }
            }
        }else if(interaction_type.equals("TRAFFIC_JAM")){
            int matrixdim = myCell.mySim.LASTPOLYMERASEID - myCell.mySim.FIRSTPOLYMERASEID + 1;
            polintermatrix = new int[matrixdim][matrixdim];
            for(int i=0;i<polintermatrix.length;i++){
                for(int j=0;j<polintermatrix[0].length;j++){
                    polintermatrix[i][j] = 0;
                }
            }
        }else if(interaction_type.equals("CUSTOM")){
			polintermatrix = new int[6][6];
			
			for(int i=0;i<polintermatrix.length;i++){
                for(int j=0;j<polintermatrix[0].length;j++){
                    polintermatrix[i][j] = 0;
                }
            }
			
            for(int i=4;i<polintermatrix.length;i++){
                for(int j=0;j<4;j++){
                    polintermatrix[i][j] = -1;
                }
            }
		}else{
            throw new IOException("Polymerase interaction model input incorrectly");
        }
        myCell.polymeraseInteractionMatrix = polintermatrix;

        //CUSTOM SPECIES INPUT
        myCell.mySim.FIRSTPROTEINID = moleculeID;
        for(Iterator i = cellxmo.getNamedChildren("SPECIES").iterator(); i.hasNext();) {
            XMLObject speciesxmo = (XMLObject) i.next();
            String name = speciesxmo.getAttribute("name");
            int init = Integer.parseInt(speciesxmo.getAttribute("n"));
            int ID = Integer.parseInt(speciesxmo.getAttribute("ID"));
            int organism = Integer.parseInt(speciesxmo.getAttribute("organism"));
            Protein p = new Protein(name,init,ID,organism);
            molecule.add(p);moleculeID++;
        }
        myCell.mySim.LASTPROTEINID = moleculeID-1;

        //create myMolecules array
        myCell.myMolecules = new Molecule[molecule.size()];
		//create moleculeEnergy array
		myCell.moleculeEnergy = new int[molecule.size()];
        for(int i=0;i<molecule.size();i++){
            myCell.myMolecules[i] = (Molecule) molecule.elementAt(i);
        }
    
        //Custom Reactions
        XMLObject reactionsXMO = cellxmo.getChild("REACTIONS");
        for(Iterator i = reactionsXMO.getNamedChildren("REACTION").iterator(); i.hasNext();) {
            XMLObject reactXMO = (XMLObject) i.next();
            String reaction_desc = reactXMO.getContents().trim();
            StringTokenizer st = new StringTokenizer(reaction_desc," >+;\t\r\n\f",true);
            //take in first number
            int reacid1 = Integer.parseInt(st.nextToken());
            int reacid2 = 0;
            String n = st.nextToken();
            if(n.equals(">")){
            }else{
                if(!n.equals("+")){
                    throw new IOException("First number not followed by + or > ");
                }
                n = st.nextToken();
                reacid2 = Integer.parseInt(n);
                n = st.nextToken();
                if(!n.equals(">")){
                    throw new IOException("Too many reactants");
                }
            }
            n = st.nextToken();
            int prodid1 = Integer.parseInt(n);
            int prodid2 = 0;
            n = st.nextToken();
            if(n.equals(";")){
            }else{
                if(!n.equals("+")){
                    throw new IOException("First number in reactants not followed by + or ; ");
                }
                n = st.nextToken();
                prodid2 = Integer.parseInt(n);
                n = st.nextToken();
                if(!n.equals(";")){
                    throw new IOException("Too many products");
                }
            }
            double macroRateConstant = Double.parseDouble(st.nextToken());
            double rateConstant = macroRateConstant;        //This is for 1st order reaction
            double secondOrderRateConstant = (rateConstant/(AVOGADRO*myCell.VOLUME));   // for 2nd order rections
            Reaction r;
            if(reacid1==0){
                //no reactants.. must be only one product
                if(prodid1==0||prodid2!=0){
                    throw new IOException("Zero Order Reactions can only have one product");
                }
                r = new ZeroOrderRxn(reactionid,myCell.myMolecules[prodid1-1],rateConstant);
            }else{
                if(reacid2==0){
                    //only 1 reactant
                    if(prodid1==0){
                        //no products
                        //System.out.println("myMolecules size = " + myCell.myMolecules.length);
						r = new FirstOrderRxn0(reactionid,myCell.myMolecules[reacid1-1], rateConstant);
						
                    }else{
                        if(prodid2==0){
                            //one product
                            r = new FirstOrderRxn1(reactionid,myCell.myMolecules[reacid1-1],myCell.myMolecules[prodid1-1],rateConstant);
                        }else{
                            //two products
                            r = new FirstOrderRxn2(reactionid,myCell.myMolecules[reacid1-1],myCell.myMolecules[prodid1-1],myCell.myMolecules[prodid2-1],rateConstant);
                        }
                    }
                }else{
                    //two reactants
                    if(prodid1==0){
                        //no products -- doesn't happen
                        throw new IOException("Does not support second order reactions with no product");
                    }else{
                        if(prodid2==0){
                            //one product
                            r = new SecondOrderRxn(reactionid,myCell.myMolecules[reacid1-1],myCell.myMolecules[reacid2-1],myCell.myMolecules[prodid1-1],secondOrderRateConstant);
                        }else{
                            //two products
                            r = new SecondOrderRxn2(reactionid,myCell.myMolecules[reacid1-1],myCell.myMolecules[reacid2-1],myCell.myMolecules[prodid1-1],myCell.myMolecules[prodid2-1],secondOrderRateConstant);
                        }
                    }
                }
            }
            if(r==null){
                throw new IOException("Something went wrong");
            }else{
                reactions.add(r);reactionid++;
            }
        }

        // Added by Jason Kelly July 15, 2003
        // Polymerase Gain Reactions
        for (int i = myCell.mySim.FIRSTPOLYMERASEID-1;i<myCell.mySim.LASTPOLYMERASEID;i++){
            reactions.add(new ZeroOrderRxn(reactionid,myCell.myMolecules[i],0));reactionid++;
            ((Protein) myCell.myMolecules[i]).gainReaction = (Reaction) reactions.elementAt(reactionid-2);
        }
        // End Addition by Jason Kelly July 15, 2003

        //Polymerase Loss Reactions
        for (int i = myCell.mySim.FIRSTPOLYMERASEID-1;i<myCell.mySim.LASTPOLYMERASEID;i++){
            reactions.add(new FirstOrderRxn0(reactionid,myCell.myMolecules[i],0));reactionid++;
            ((Polymerase) myCell.myMolecules[i]).setLossReaction((FirstOrderRxn0) reactions.elementAt(reactionid-2));
        }
        
        //Ribosome Gain Reactions
        reactions.add(new ZeroOrderRxn(reactionid,myCell.myMolecules[myCell.mySim.RIBOSOMEID-1],0));reactionid++;
        ((Protein) myCell.myMolecules[myCell.mySim.RIBOSOMEID-1]).gainReaction = (Reaction) reactions.elementAt(reactionid-2);

        //Protein Entry Reactions
        for(int i=myCell.mySim.FIRSTPROTEINID;i<=myCell.mySim.LASTPROTEINID;i++){
            Protein prot = (Protein) myCell.myMolecules[i-1];
            reactions.add(new ZeroOrderRxn(reactionid,prot,0));reactionid++;
            prot.gainReaction = (Reaction) reactions.elementAt(reactionid-2);
        }
        Object[] tempvect = reactions.toArray();
        Reaction[] reactvect = new Reaction[tempvect.length];
        for(int i=0;i<reactvect.length;i++){
            reactvect[i] = (Reaction) tempvect[i];
        }
        return reactvect;
    }


	
	/** 
	* Parses parameters to initialize the Phage class
	* @param myCell A pointer to the cell class the phage is contained within
	* @param myPhage A pointer to the phage class to initialize
	* @param phagexmo The XML object from which to parse information
	* @param moleculeID The moleculeID to begin adding molecules (Cell class already added molecules)
	*/	
    public static void fillPhage(Cell myCell, Phage myPhage, XMLObject phagexmo, int moleculeID) throws IOException{
        Vector[] tempDnaElements = new Vector[2];
        tempDnaElements[0] = new Vector(); // Promoters
        tempDnaElements[1] = new Vector(); // Terminators
        //tempDnaElements[2] = new Vector(); // RnaseIII Sites
                                           //tempDnaElements[3] = new Vector();
                                           //int T7_PRO_ID = 1000;
                                           //int COLI_PRO_ID = 2000;


        Vector molecule = new Vector(); //this is only for the GammaJMPRxn molecules, all other (free-species) molecules were defined in myCell

        final double AVOGADRO = 6.022e23;
        int proCount = 0; //this is a counter so can keep track of the ID for both types of promoters
        int NUM_OF_RIB = 1; //will always be one.

        //get genome length
        myPhage.GENOME_LENGTH = Integer.parseInt(phagexmo.getAttribute("genome_length"));

        //promoters
		int ID = myPhage.myCell.mySim.T7_PRO_ID;
        for(Iterator it = phagexmo.getNamedChildren("PROMOTER").iterator(); it.hasNext();) {
            XMLObject promoterxmo = (XMLObject) it.next();
            int polyIndex = 0;
            int numOfPol = myPhage.myCell.mySim.LASTPOLYMERASEID - myPhage.myCell.mySim.FIRSTPOLYMERASEID + 1 + NUM_OF_RIB;
            int start = Integer.parseInt(promoterxmo.getAttribute("start"));
            int stop = Integer.parseInt(promoterxmo.getAttribute("stop"));
            int startsite = Integer.parseInt(promoterxmo.getAttribute("startsite"));
            int organism_id = Integer.parseInt(promoterxmo.getAttribute("organism_id"));
            String name = promoterxmo.getAttribute("name");
            double[] aon = new double[numOfPol];
            double[] aoff = new double[numOfPol];
            double[] ainiton = new double[numOfPol];
            double[] runoffPercent = new double[numOfPol];
            double[] aelong = new double[numOfPol];
			double[] arecyc = new double[numOfPol];
            //Values are set to -1 to initialize the matrix with something that will throw an error
            // if ever a call is made to a polymerase that shouldnt ever bind to this promoter.
            for (int i=0;i<numOfPol;i++){
                aon[i] = -1;
                aoff[i] = -1;
                ainiton[i] = -1;
                runoffPercent[i] = -1;
                aelong[i] = -1;
				arecyc[i] = -1;
            }
            for(Iterator iter = promoterxmo.getNamedChildren("POLYMERASE").iterator(); iter.hasNext();){
                XMLObject prompolxmo = (XMLObject) iter.next();
                int polID = Integer.parseInt(prompolxmo.getAttribute("polymeraseID"));
                double macro_aon =  Double.parseDouble(prompolxmo.getAttribute("aon"));
                aon[polID-1] = (macro_aon/(AVOGADRO*myCell.VOLUME));
                aoff[polID-1] = Double.parseDouble(prompolxmo.getAttribute("aoff"));
                ainiton[polID-1] = Double.parseDouble(prompolxmo.getAttribute("ainiton"));
                runoffPercent[polID-1] = Double.parseDouble(prompolxmo.getAttribute("runoff_percent"));
                aelong[polID-1] = Double.parseDouble(prompolxmo.getAttribute("aelong"));
				arecyc[polID-1] = Double.parseDouble(prompolxmo.getAttribute("arecyc"));
            }
            Promoter t;
            if(organism_id==1){
                t = new T7Promoter(start,stop,++ID,startsite,name,aon,aoff, ainiton, runoffPercent, aelong, arecyc);
            }else{
                t = new EcoliPromoter(start,stop,++ID,startsite,name,aon,aoff, ainiton, runoffPercent, aelong, arecyc);
            }
            ((Vector) tempDnaElements[0]).add(t);
            proCount++;
        }

        //terminators
		int termID = 7000;
        for(Iterator it = phagexmo.getNamedChildren("TERMINATOR").iterator(); it.hasNext();) {
            XMLObject termxmo = (XMLObject) it.next();
            int start = Integer.parseInt(termxmo.getAttribute("start"));
            int stop = Integer.parseInt(termxmo.getAttribute("stop"));
            int stopsite = Integer.parseInt(termxmo.getAttribute("stopsite"));
            String name = termxmo.getAttribute("name");
            Vector polIDVect = new Vector();
            Vector efficiencyVect = new Vector();
            for(Iterator iter = termxmo.getNamedChildren("POLYMERASE").iterator(); iter.hasNext();) {
                XMLObject termpolxmo = (XMLObject) iter.next();
                int polID = Integer.parseInt(termpolxmo.getAttribute("ID"));
                polIDVect.add(new Integer(polID));
                double efficiency = Double.parseDouble(termpolxmo.getAttribute("efficiency"));
                efficiencyVect.add(new Double(efficiency));
            }
            int[] polID_arry = new int[polIDVect.size()];
            double[] eff_arry = new double[efficiencyVect.size()];
            for(int i=0;i<polID_arry.length;i++){
                polID_arry[i] = ((Integer) polIDVect.elementAt(i)).intValue();
                eff_arry[i] = ((Double) efficiencyVect.elementAt(i)).doubleValue();
            }
            Terminator t = new Terminator(start,stop,++termID,stopsite,name,polID_arry,eff_arry);
            ((Vector) tempDnaElements[1]).add(t);

        }

        //RBS/mRNA/Translation Stop
        int index = 1;
        Vector mre = new Vector();
        for(Iterator it = phagexmo.getNamedChildren("RBS").iterator(); it.hasNext();) {
            XMLObject rbsxmo = (XMLObject) it.next();
            int start = Integer.parseInt(rbsxmo.getAttribute("start"));
            int stop = Integer.parseInt(rbsxmo.getAttribute("stop"));
            int startsite = Integer.parseInt(rbsxmo.getAttribute("startsite"));
            int initstepsize = Integer.parseInt(rbsxmo.getAttribute("initstepsize"));
			double initRateConstant = Double.parseDouble(rbsxmo.getAttribute("initRateConstant"));
            double macroStrength = Double.parseDouble(rbsxmo.getAttribute("strength"));
            double strength = macroStrength/(AVOGADRO*myCell.VOLUME);
            int stopsite = Integer.parseInt(rbsxmo.getAttribute("stopsite"));
            int protID = Integer.parseInt(rbsxmo.getAttribute("protid"));
            double mrnadeg = Double.parseDouble(rbsxmo.getAttribute("mrnadeg"));
            int elongstepsize = stopsite - startsite-initstepsize;
			TranslationStop ts = new TranslationStop(stopsite,(Protein)myCell.myMolecules[protID-1],(index)*-1);
            mre.add(ts);
            MatureRBS mrbs = new MatureRBS(start,stop,startsite,index,(TranslationStop) mre.elementAt((2*index)-2),initstepsize,elongstepsize,strength,mrnadeg,initRateConstant);
            mre.add(mrbs);
            index++;
            ts.mrbs = mrbs;
            mrbs.nrbs = new NascentRBS(index-1, mrbs);
            molecule.add(new GammaJumpRNA(moleculeID,mrbs,true,myCell.mySim));moleculeID++;
            molecule.add(new GammaJumpRNA(moleculeID,mrbs,false,myCell.mySim));moleculeID++;
            ((Protein) myCell.myMolecules[protID-1]).length = (initstepsize+elongstepsize)/3;
        }

        //make myRNA Elements
        myPhage.myRnaElements = new RNA[mre.size()];
        for(int i=0;i<mre.size();i++){
            myPhage.myRnaElements[i] = (RNA) mre.elementAt(i);
        }

        //dna entry
        myPhage.autoEntryRateConstant = Double.parseDouble(phagexmo.getAttribute("entry_rate_constant"));
        myPhage.autoEntryOffSite = Integer.parseInt(phagexmo.getAttribute("entry_offsite"));
        
        // Fix the myPhage.myCell.myMolecules array by adding the GammaJMPRxns in the molecule vector
        Molecule[] myMoleculeTemp = myPhage.myCell.myMolecules;
        myPhage.myCell.myMolecules = new Molecule[molecule.size()+myMoleculeTemp.length];
        int molCount = 0;
        for(molCount=0;molCount<myMoleculeTemp.length;molCount++){
            myPhage.myCell.myMolecules[molCount] = (Molecule) myMoleculeTemp[molCount];
        }
        for(int i=molCount;i<(molCount + molecule.size());i++){
            myPhage.myCell.myMolecules[i] = (Molecule) molecule.elementAt(i-molCount);
        }
        
        //Convert myDNAELements from vector to array.
        myPhage.myDnaElements = new DNA[tempDnaElements.length][];
        
        for (int i=0;i<myPhage.myDnaElements.length;i++){
            myPhage.myDnaElements[i] = new DNA[tempDnaElements[i].size()];
            for (int j=0;j<tempDnaElements[i].size();j++){
                myPhage.myDnaElements[i][j] =  (DNA) tempDnaElements[i].elementAt(j);
            }
        }
    }
    
	/** 
	* Parses parameters to fill in reactions from phage class
	* @param myPhage A pointer to the phage that is being initialized
	* @param reactionid The reactionID to add reactions from.
	* @param GJRNAINDEX The index in the molecules vector that the GammaJumpRNA start.
	* @return A vector containing reactions to add to the simulation class.
	*/		
    public static Reaction[] fillPhageReactions(Phage myPhage, int reactionid,  int GJRNAINDEX){
        //rnapEVRxn acts as a pointer to the elongation event rxn to be used for passing to other rxns.
        RNAPElongationEvent rnapEVRxn;
        int NUM_OF_RIBO = 1;
        Vector reactions = new Vector();
        //GP16Entry
        reactions.add(new GP16EntryEvent(reactionid,myPhage.autoEntryRateConstant, myPhage.autoEntryOffSite,myPhage));reactionid++;
        //System.out.println("Input part 11");
        
        //promoter and Terminator entry events
        for(int j=0;j<2;j++){
            for(int i=0;i<myPhage.myDnaElements[j].length;i++){
                reactions.add(new ZeroOrderRxn(reactionid,(DNA) myPhage.myDnaElements[j][i],0));
                ((DNA) myPhage.myDnaElements[j][i]).setEntryReaction((Reaction) reactions.lastElement());reactionid++;
                reactions.add(new FirstOrderRxn0(reactionid,(DNA) myPhage.myDnaElements[j][i],0));
                ((DNA) myPhage.myDnaElements[j][i]).setLossReaction((Reaction) reactions.lastElement());reactionid++;
            }
        }
        
        //rnap elongation
        rnapEVRxn = new RNAPElongationEvent(reactionid,myPhage,myPhage.myCell);
        reactions.add(rnapEVRxn);reactionid++;
        
        //promoter binding reactions
        
        // for loop goes through each promoter
        for(int i=0;i<(myPhage.myDnaElements[0].length);i++){
            //for loop to cycle through the various polymerases that act with that promoter
            for (int j = myPhage.myCell.mySim.FIRSTPOLYMERASEID-1;j<myPhage.myCell.mySim.LASTPOLYMERASEID;j++){
                //if statement to check if the rates are = -1 (i.e. not intitialized) otherwise make a reaction.
                if (((Promoter) myPhage.myDnaElements[0][i]).aon[j] != -1){
                    reactions.add(new PromoterBindingRxn(reactionid,myPhage.myCell.myMolecules[j],(DNA) myPhage.myDnaElements[0][i], rnapEVRxn,myPhage));reactionid++;
                }
                
            }
            
        }
        //E. coli promoters
        /*
        for(int i=0;i<(myPhage.myDnaElements[1].size());i++){
            reactions.add(new PromoterBindingRxn(reactionid,myPhage.myCell.myMolecules[myPhage.myCell.mySim.COLIPOLYMERASEID-1],(DNA) myPhage.myDnaElements[1].elementAt(i), rnapEVRxn,myPhage));reactionid++;
        }*/
        
        //termination reactions
        //int countRelease = 0;
        myPhage.releaseRxn = new DPComplexReleaseRxn[myPhage.myCell.mySim.LASTPOLYMERASEID];
        //for (int i = myPhage.myCell.mySim.FIRSTPOLYMERASEID-1;i<myPhage.myCell.mySim.LASTPOLYMERASEID;i++){
        //The first poistion is empty ebcause prot.ID-1 is used to index this, and ribosome is in the first slot.
        for(int i=NUM_OF_RIBO;i<myPhage.myCell.mySim.LASTPOLYMERASEID;i++){
            Polymerase pol = (Polymerase) myPhage.myCell.myMolecules[i];
            DPComplexReleaseRxn termrelease = new DPComplexReleaseRxn(reactionid,pol,0,rnapEVRxn);
            myPhage.releaseRxn[i] = termrelease;
            reactions.add(termrelease);reactionid++;
            //countRelease++;
        }
        /*
        Polymerase polt7 = (Polymerase) myPhage.myCell.myMolecules[myPhage.myCell.mySim.T7POLYMERASEID-1];
        DPComplexReleaseRxn termrelease2 = new DPComplexReleaseRxn(reactionid,polt7,0,rnapEVRxn);
        myPhage.t7ReleaseRxn = termrelease2;
        reactions.add(termrelease2);reactionid++;
         */
        // What the heck is htis for?? (jrk)
        /*
        for(int i=0;i<(myPhage.myDnaElements[2].size());i++){
            Terminator term = (Terminator) myPhage.myDnaElements[2].elementAt(i);
            Polymerase pol = null;
            if(term.organism==1){
                pol = (Polymerase) myPhage.myCell.myMolecules[myPhage.myCell.mySim.T7POLYMERASEID-1];
            }else if(term.organism==2){
                pol = (Polymerase) myPhage.myCell.myMolecules[myPhage.myCell.mySim.COLIPOLYMERASEID-1];
            }
        }*/
        
        //MatureRBS entry rxns, conversion rxns, and NascentRBSEntryRxns
        for(int i=1;i<myPhage.myRnaElements.length;i+=2){
            //mature rbs conversion rxn
            reactions.add(new RBSConversionRxn(reactionid,((MatureRBS) myPhage.myRnaElements[i]).nrbs,myPhage.myRnaElements[i], 0));reactionid++;
            ((MatureRBS) myPhage.myRnaElements[i]).setConversionReaction((Reaction) reactions.lastElement());
            //mature rbs entry rxn
            reactions.add(new ZeroOrderRxn(reactionid,myPhage.myRnaElements[i],0));reactionid++;
            ((MatureRBS) myPhage.myRnaElements[i]).setEntryReaction((Reaction) reactions.lastElement());
            //nrbs entry reaction
            reactions.add(new NascentRBSEntryRxn(reactionid,((MatureRBS) myPhage.myRnaElements[i]).nrbs,0));reactionid++;
            ((MatureRBS) myPhage.myRnaElements[i]).nrbs.entryReaction = (NascentRBSEntryRxn) reactions.lastElement();
            reactions.add(new NascentRBSLossRxn(reactionid,((MatureRBS) myPhage.myRnaElements[i]).nrbs,0));reactionid++;
            ((MatureRBS) myPhage.myRnaElements[i]).nrbs.lossReaction = (NascentRBSLossRxn) reactions.lastElement();
        }
        
        //MrnaTrackerReactions
        MrnaTrackerRxn mtrTemp = new MrnaTrackerRxn(reactionid,0,(Ribosome) myPhage.myCell.myMolecules[myPhage.myCell.mySim.RIBOSOMEID-1], myPhage, myPhage.myCell);
        rnapEVRxn.mtr = mtrTemp;
        reactions.add(mtrTemp);
        reactionid++;
        //System.out.println(myPhage.myCell.myMolecules[myPhage.myCell.mySim.RIBOSOMEID-1]);
        
        //Ribosome Binding Reaction
        Ribosome ribTemp = (Ribosome) myPhage.myCell.myMolecules[myPhage.myCell.mySim.RIBOSOMEID-1];
        
        int gjrnaindex = GJRNAINDEX;
        for(int i=1;i<myPhage.myRnaElements.length;i+=2){
            MatureRBS mrbs = (MatureRBS) myPhage.myRnaElements[i];
            GammaJumpRNA initrna = (GammaJumpRNA) myPhage.myCell.myMolecules[gjrnaindex-1];gjrnaindex+=2;
            reactions.add(new MatureRBSBindingRxn(reactionid,ribTemp,mrbs,initrna,mrbs.strength));reactionid++;
            mrbs.mrbsBindingRxn = (Reaction) reactions.lastElement();
            reactions.add(new NascentRBSBindingRxn(reactionid,ribTemp,mrbs.nrbs,mrbs.strength,mtrTemp));reactionid++;
        }
        
        //Mrna Degradation Reactions
        for(int i=1;i<myPhage.myRnaElements.length;i+=2){
            MatureRBS mrbs = (MatureRBS) myPhage.myRnaElements[i];
            reactions.add(new FirstOrderRxn0(reactionid,mrbs,mrbs.mrnadeg));reactionid++;
            //System.out.println("The mrna degradation rate is: " +mrbs.mrnadeg);
        }
        
        gjrnaindex = GJRNAINDEX;
        //Mature Mrna Initation and Elongation Reactions and initMRNA degradation reactions
        ribTemp = (Ribosome) myPhage.myCell.myMolecules[myPhage.myCell.mySim.RIBOSOMEID-1];
        for(int i=1;i<myPhage.myRnaElements.length;i+=2){
            MatureRBS mrbs = (MatureRBS) myPhage.myRnaElements[i];
            GammaJumpRNA initrna = (GammaJumpRNA) myPhage.myCell.myMolecules[gjrnaindex-1];gjrnaindex++;
            mrbs.initrna = initrna;
            GammaJumpRNA elongrna = (GammaJumpRNA) myPhage.myCell.myMolecules[gjrnaindex-1];gjrnaindex++;
			InitiationMRBSRxn imr = new InitiationMRBSRxn(reactionid,initrna,elongrna,mrbs,ribTemp.speed,mrbs.initstepsize, myPhage.myCell,mrbs.initRateConstant);
            reactions.add(imr);reactionid++;
            initrna.rxnToUpdate = (GammaJumpRxn) reactions.lastElement();
			if(!myPhage.myCell.mySim.compDeg){
				reactions.add(new InitiationMRBSDegRxn(reactionid,initrna,ribTemp,mrbs.mrnadeg,imr));reactionid++;
			}
			reactions.add(new ElongationMRBSRxn(reactionid,elongrna,mrbs.ts.product,ribTemp,ribTemp.speed,mrbs.elongstepsize,myPhage.myCell));reactionid++;
            elongrna.rxnToUpdate = (GammaJumpRxn) reactions.lastElement();
        }
        
        //return the reaction vectors
        Reaction[] retreac = new Reaction[reactions.size()];
        for(int i=0;i<reactions.size();i++){
            retreac[i] = (Reaction) reactions.elementAt(i);
        }
        return retreac;
    }
}
