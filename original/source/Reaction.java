import java.util.Vector;
import java.io.IOException;

/**
 * A general class for Reactions.  Most instantiations of Reaction actually instantiate from subclasses.  This class should perhaps become abstract.
 * @version 1.0
 * @author Sriram Kosuri
 */
public class Reaction {
	/**
	 * The reactants of the reaction whose copy numbers will be decremented upon reaction execution.
	 */
    Molecule[] Reactants;
	/**
     * The products of the reaction whose copy numbers will be incremented upon reaction execution.
	 */
    Molecule[] Products;
	/**
	 * The reactions ID.  Relates to the indices in the IndexedPriorityQueue.
	 */
    protected int ID;
	/**
	 * The value of the mesoscopic rate constant of this reaction.
	 */
    double rateConstant;
	/**
	 * The next time that this reaction is set to execute.  
	 */
    double time;
	/**
	 * A list of reactions that should be updated upon this reaction's execution.
	 */
    Vector updateRxn;
	/**
	 * The propensity of the reaction.  The propensity is usually defined as the mesoscopic rate constant times the number of molecules of each reactant.
	 */
    double a;
    
	/**
	 * A default empty constructor.  Assumes very little about the number of reactants and products.  
	 * @param ID the reaction ID
	 * @param rc the rate constant of the reaction.
	 */
    protected Reaction(int ID, double rc){
		this.ID = ID;
		rateConstant = rc;
		updateRxn = new Vector();
    }
    
	/**
	 * A constructor that takes in reactant and product molecules as well.
	 * @param ID the reaction ID
	 * @param reac The reactants to be decremented upon execution of this reaction.
	 * @param prod The products to be incremented upon execution of this reaction.
	 * @param rc The mesoscopic rate constant of the reaction.
	 */
    public Reaction(int ID, Molecule[] reac, Molecule[] prod, double rc){
        this.ID = ID;
        this.Reactants = reac;
        this.Products = prod;
        this.updateRxn = new Vector();
        this.rateConstant = rc;
    }
    
	/**
	 * Returns the ID of the reaction.  This is linked to the address in Gintegrator.
	 * @return the ID of the reaction.
	 */
    public int getID(){
        return ID;
    }
	
	/**
	 * Returns the time that this reaction is next expected to execute.  
	 * @return the time that this reaction is next expected to execute.  
	 */
    public double getTime(){
		return time;
    }

	/**
	 * Sets the time that this reaction is next expected to execute.
	 * @param t The time that this reaction is next expected to execute.
	 */
    public void setTime(double t){
        time = t;
    }

	/**
	 * Returns the current propensity of the reaction.  The propensity, a,  is usually defined as the mesoscopic rate constant times the number of molecules of each reactant.
	 * @return the current propensity of the reaction.
	 */
    public double getA(){
        return a;
    }

	/**
	 * Compares this reaction time to execute to another reaction.  It subtracts rhs.time from the current time of this reaction.  This is useful for ordering in the IndexedPriorityQueue.
	 * @param rhs The reaction to compare time to next execution with.
	 */
    public double compareTo(Reaction rhs){
        return (this.time - rhs.time);
    }
    
	/**
	 * Calculate the propensity of the reaction.  This is accomplished here as the mesoscopic rate constant times the number of molecules of each reactant.
	 */
    public void calculateA(){
        a = rateConstant;
        for(int i=0;i<Reactants.length;i++){
            a *= Reactants[i].getCopyNumber();
        }
    }
	
	/**
	 * Calculate the time the reaction will next execute.  Assumes an exponential distribution, witht he propensity being the average.  
	 * @param random A randomly generated double value between [0,1) that will be used to choose a time based on the distribution.
	 * @param t The current time of the simulation.
	 */
    public void calculateTime(double random, double t){
        time = (1/a)*Math.log(1/random) + t;
    }
    
	/**
	 * Executes the reaction.  Increments the copy numbers of the reactants and decrements the copy numbers of the products.
	 */
    public void execute(){
        //update reactants
        for(int i=0;i<Reactants.length;i++){
            Reactants[i].incrementCopyNumber(-1);
        }
        //update products
        for(int i=0;i<Products.length;i++){
            Products[i].incrementCopyNumber(1);
        }
    }

	/**
	 * Returns the set of reactants that this reaction depends on.  In this case, just the reactants.
	 * @return an array of reactant molecules that this reaction depends on.
	 */ 
    public Molecule[] dependsOn(){
        return Reactants;
    }

	/**
	 * Returns all the molecules that are affected upon the execution of this reaction.  In this case, all the reactants and products.
	 * @return an array of reactant molecules that this reaction affects.
	 */
    public Molecule[] affects(){
        int length = Reactants.length+Products.length;
        Molecule[] aff = new Molecule[length];
        int i;
        for(i=0;i<Reactants.length;i++){
            aff[i] = Reactants[i];
        }
        for(int j=0;j<Products.length;j++){
            aff[i+j] = Products[j];
        }
        return aff;
    }

	/**
	 * A method that is used upon initialization to populate the updateRxn vector.  Checks if this reaction affects a particular reaction R (in which case this the particular 
	 * reaction will be added to the reaction's updateRxn vector) and if that particular reaction R affects this reaction (in which case the particular reactions updateRxn vector will be appended).  
	 */
    public void findAffects(Reaction r){
        //First see if this object affects r
        Molecule[] aff = this.affects();
        Molecule[] dep = r.dependsOn();
        for(int i=0;i<dep.length;i++){
            for(int j=0;j<aff.length;j++){
                if(dep[i].equals(aff[j])){
                    boolean already_added = false;
                    Integer newint = new Integer(r.getID());
                    for(int k=0;k<this.updateRxn.size();k++){
                        if(this.updateRxn.elementAt(k).equals(newint))
                            already_added=true;
                    }
                    if(!already_added)
                        this.updateRxn.add(newint);
                }

            }
        }
        
        
        //Now see if r affects this object
        if(this.ID!=r.ID){
            aff = r.affects();
            dep = this.dependsOn();
            for(int i=0;i<dep.length;i++){
                for(int j=0;j<aff.length;j++){
                    if(dep[i].equals(aff[j])){
                        boolean already_added = false;
                        Integer newint = new Integer(this.getID());
                        for(int k=0;k<r.updateRxn.size();k++){
                            if(r.updateRxn.elementAt(k).equals(newint))
                                already_added=true;
                        }
                        if(!already_added)
                            r.updateRxn.add(newint);
                    }
                }
            }
        }
    }
    
    
}

/**
 * This reaction is a specialized reaction that tracks the movements and positions of polymerases and proteins on the DNA.  
 * @version 1.0
 * @author Sriram Kosuri
 */
class RNAPElongationEvent extends Reaction{
	/**
	 * A placeholder array of size zero.  This reaction does not update other reactions through the normal method of affects and dependencies.
	 */
    protected Molecule[] affdep;
	/**
	 * The priorityqueue that stores the complexes that this reaction is tracking on the DNA
	 */	
    protected PriorityQueue myPQ;
	/**
	 * A pointer to the DNA molecule that this reaction tracks complexes on.
	 */
    protected Phage myPhage;
	/**
	 * A pointer to the cell #myPhage is contained within.
	 */	
    protected Cell myCell;
	/**
     * A pointer reaction to the corresponding tracker reaction of the complexes on the RNA.
	 */	
    MrnaTrackerRxn mtr;
    
	/**
	 * A general constructor that should be labeled private in the future.  
	 * @param ID the reaction ID
	 * @param myPhage the DNA molecule that this reactions tracks complexes on.
	 * @param myCell the cell that #myPhage is contained within.
	 */
    public RNAPElongationEvent(int ID, Phage myPhage, Cell myCell){
        super(ID,null,null,0);
        affdep = new Molecule[0];
        myPQ = new PriorityQueue();
        this.myPhage = myPhage;
        updateRxn.add(new Integer(ID));
        a=0;

        this.myCell = myCell;

    }

	/**
	 * The main constructor of this class.  
	 * @param ID the reaction ID
	 * @param myPhage the DNA molecule that this reactions tracks complexes on.
	 * @param myCell the cell that #myPhage is contained within.
	 * @param mtr The reaction to the corresponding tracker reaction of the complexes on the RNA generated on #myPhage.
	 */
    public RNAPElongationEvent(int ID, Phage myPhage, Cell myCell, MrnaTrackerRxn mtr){
        this(ID,myPhage,myCell);
        this.mtr = mtr;
    }

	/**
	 * The method that is called when a new complex forms on the DNA.
	 * @param dp The complex that is to be added to #myPQ.
	 */
    public void addComplex(DPComplex dp){
        myPQ.insert(dp);
        try{
            time = ((DPComplex) myPQ.findMin()).time;
            myCell.mySim.gint.myIPQ.justUpdate(ID);
        }catch(Exception e){
            e.printStackTrace();
            System.exit(9);
        }
    }
	/**
     * Overrides the inherited method so that this method does nothing.
	 */
    public void calculateA(){
    }

	/**
	 * Overrides the inherited method so that this method returns an empty array.
	 * @return An empty array.
	 */	
    public Molecule[] affects() {
        return affdep;
    }

	/**
	 * Overrides the inherited method so that this method returns an empty array.
	 * @return An empty array.
	 */	
    public Molecule[] dependsOn() {
        return affdep;
    }
	
	/**
	 * Overrides thte inherited method.  This method deletes the complex that just reacted, recalculates it's reaction time, reinserts in into the PriorityQueue, and then takes the new complex that will react
	 * next as this reaction's new time to react.
	 * @param random The random number that will be used to update the executing DPComplex's next reaction time.
	 * @param t The current time of the simulation.
	 */
    public void calculateTime(double random, double t){

		try{
            DPComplex dp = (DPComplex) myPQ.deleteMin();
            dp.calculateTime(random,time);
            myPQ.insert(dp);
            time = ((DPComplex) myPQ.findMin()).time;
        }catch(java.io.IOException e){
            time = Double.POSITIVE_INFINITY;
        }
    }

	/**
	 * Overrides the inherited method.  This method takes the complex at the top of the PriorityQueue, and executes it using the #move() method.
	 */
    public void execute(){
        try{
        DPComplex dp = (DPComplex) myPQ.findMin();
        move(dp);
        
        }catch(Exception e){
            e.printStackTrace();
            System.exit(9);
        }
    }
	
	//called to terminate an mRNA
	/**
	 * The method is called when a particular polymerase falls off the DNA.  In particular it is coded to remove all ribosomes translating incomplete mRNA.  
	 * @param dp The DPComplex to be terminated.
	 */
    public void terminateDPComplex(DPComplex dp) throws IOException{
		//this polymerase falls off ... taken from termination of polymerase
		//adds to the polymerase count and removes polymerase from the pq that it belongs to.
		if(dp.whichReaction!=1){
			DPComplexReleaseRxn r = dp.getReleaseRxn();
            r.execute(dp);
            for(int i=0;i<(r.updateRxn).size();i++){
                int reacid = ((Integer) r.updateRxn.elementAt(i)).intValue();
                if(reacid!=r.getID())
                    myCell.mySim.gint.myIPQ.updateOther(reacid, myCell.mySim.gint.t);
            }
            //clear myDNAComplexes vector
            for(int k=dp.stop;k<=dp.start;k++){
                myPhage.myDNAComplexes[k] = null;
            }
                        
            //read promoter complex
            if(dp.myPromoter.stillBlocked>0){
				dp.myPromoter.stillBlocked--;
				
			}else{
				if(dp.myPromoter.copyNumber==0){
					
					Reaction er = dp.myPromoter.getEntryReaction();
					er.execute();
					for(int i=0;i<(er.updateRxn).size();i++){
						int reacid = ((Integer) er.updateRxn.elementAt(i)).intValue();
						if(reacid!=er.getID())
							myCell.mySim.gint.myIPQ.updateOther(reacid, myCell.mySim.gint.t);
					}
				}else{
					System.out.println("ARG: Reaction uncovering dna... copy # = " + dp.myPromoter.getCopyNumber());
					System.out.println(dp.myPromoter);
					System.exit(9);
				}
            }
		}else{
			//remove stalled ribosomes
			Reaction ribgainreac = mtr.rib.gainReaction;
			for(int i=0;i<dp.stalledRibosomes.size();i++){			
				//return a ribosome
				ribgainreac.execute();
				myPhage.stalled_rib--;
			}

			dp.stalledRibosomes.removeAllElements();
			for(int j=0;j<(ribgainreac.updateRxn).size();j++){
				int reacid = ((Integer) ribgainreac.updateRxn.elementAt(j)).intValue();
				if(reacid!=ribgainreac.getID())
					myCell.mySim.gint.myIPQ.updateOther(reacid, myCell.mySim.gint.t);
			}
			
			//remove any nascent ribosomes that may be associated with
			//this transcribing polymerase
			for(int i=1;i<=mtr.myPQ.size();i++){
				RRComplex rr = (RRComplex) mtr.myPQ.elementAt(i);
				if(rr.nrna.dp==dp){
					//!totallyentered_init
					if(!rr.totallyentered&&rr.init){
						//return a nascent rbs such that it can be removed later
						//Nascent rbs entry reaction
						Reaction r = rr.rbs.entryReaction;
						r.execute();
						for(int j=0;j<(r.updateRxn).size();j++){
							int reacid = ((Integer) r.updateRxn.elementAt(j)).intValue();
							if(reacid!=r.getID())
								myCell.mySim.gint.myIPQ.updateOther(reacid, myCell.mySim.gint.t);
						}
						for(int k=0;k<rr.rbs.myRBSInstances.size();k++){
							RBSInstance rbsi = (RBSInstance) rr.rbs.myRBSInstances.elementAt(k);
							if(dp==rbsi.dp){
								if(rbsi.isBound){
									rbsi.isBound = false;
									rbsi.rr=null;
									break;
								}else{
									System.out.println("Bound/Unbound nrbs mismatch: exiting");
									System.exit(9);
								}
							}
						}
						
						
						//return a ribosome
						//return ribosome for both init and elong complexes only if cds has not entered.
						ribgainreac.execute();
						for(int j=0;j<(ribgainreac.updateRxn).size();j++){
							int reacid = ((Integer) ribgainreac.updateRxn.elementAt(j)).intValue();
							if(reacid!=ribgainreac.getID())
								myCell.mySim.gint.myIPQ.updateOther(reacid, myCell.mySim.gint.t);
						}
						
						//delete rr from mrna_tracker_rxn
						//remove the RRComplex from the priority queue only if cds is not entered.
						if(i!=1){
							mtr.myPQ.delete(i);
						}else
							mtr.myPQ.deleteMin();
						
						i=0;
					}
					//!totallyentered_elong
					else if(!rr.totallyentered&&!rr.init){
						//return a ribosome
						//return ribosome for both init and elong complexes only if cds has not entered.
						ribgainreac.execute();
						for(int j=0;j<(ribgainreac.updateRxn).size();j++){
							int reacid = ((Integer) ribgainreac.updateRxn.elementAt(j)).intValue();
							if(reacid!=ribgainreac.getID())
								myCell.mySim.gint.myIPQ.updateOther(reacid, myCell.mySim.gint.t);
						}
						
						//delete rr from mrna_tracker_rxn
						//remove the RRComplex from the priority queue only if cds is not entered.
						if(i!=1){
							mtr.myPQ.delete(i);
						}else
							mtr.myPQ.deleteMin();
						
						i=0;
					}
					//totallyentered_init
					else if(rr.totallyentered&&rr.init){
						//do nothing... will take care of itself.
					}
					//totallyentered_elong
					else if(rr.totallyentered&&!rr.init){
						//do nothing... will take care of itself.
					}
					
					
				}//dp = dp
			}//end for loop
			
			//remove from nrbs... if it is in the middle of a coding region
			for(int i=0;i<myPhage.myRnaElements.length;i+=2){
				TranslationStop ts = (TranslationStop) myPhage.myRnaElements[i];
				MatureRBS mrbs = (MatureRBS) myPhage.myRnaElements[i+1];
				int startsite_of_protein = mrbs.startsite;
				int endsite_of_protein = ts.stopsite;
				//if myRBSInstance is initiated but not completed
				if(dp.m.start<=startsite_of_protein && dp.m.stop<endsite_of_protein){
					for(int j=0;j<mrbs.nrbs.myRBSInstances.size();j++){
						RBSInstance rbsi = (RBSInstance) mrbs.nrbs.myRBSInstances.elementAt(j);
						if(dp==rbsi.dp){
							if(rbsi.rr==null||!rbsi.rr.totallyentered){
								//remove a ribosome binding site through reaction
								Reaction r = mrbs.nrbs.lossReaction;
								r.execute();
								for(int k=0;k<(r.updateRxn).size();k++){
									int reacid = ((Integer) r.updateRxn.elementAt(k)).intValue();
								if(reacid!=r.getID())
									myCell.mySim.gint.myIPQ.updateOther(reacid, myCell.mySim.gint.t);
								}
							
								//remove from vector
								mrbs.nrbs.myRBSInstances.remove(j);j--;
							}
						}
					}
				}
			}
						
			//remove from mrnatrackerrxn
			DPComplexReleaseRxn dpreleasereac = dp.getReleaseRxn();
			dpreleasereac.execute(dp);			
			for(int n=0;n<(dpreleasereac.updateRxn).size();n++){
				int reacid = ((Integer) dpreleasereac.updateRxn.elementAt(n)).intValue();
				if(reacid!=dpreleasereac.getID())
					myCell.mySim.gint.myIPQ.updateOther(reacid, myCell.mySim.gint.t);
			}
			
			
			//reset time on mtr
			myCell.mySim.gint.myIPQ.updateOther(mtr.getID(), myCell.mySim.gint.t);
			
			//check uncovering of any genetic elements. be careful of boundaries.
			int last_updated = -1;
			for(int k=dp.stop;k<=dp.start;k++){
				//need some form of memory to make sure i don't continually 'uncover' same element
				if(myPhage.myDNAFeatures[k]!=0&&myPhage.myDNAFeatures[k]!=last_updated){
					DNA d = myPhage.convertFeatureToElement(k);
					last_updated = myPhage.myDNAFeatures[k];
					if(d.stillBlocked>0){
						d.stillBlocked--;
					}else{
						if(d.getCopyNumber()==0){
							Reaction dnaentryreac = d.getEntryReaction();
							dnaentryreac.execute();
							//for entry events only
							for(int q=0;q<(dnaentryreac.updateRxn).size();q++){
								int reacid = ((Integer) dnaentryreac.updateRxn.elementAt(q)).intValue();
								if(reacid!=dnaentryreac.getID())
									myCell.mySim.gint.myIPQ.updateOther(reacid, myCell.mySim.gint.t);
							}
						}else{
							System.out.println("ARG: Reaction uncovering dna... copy # = " + d.getCopyNumber());
							System.out.println(d);
							System.exit(9);
						}
					}
					
				}
				
			}
			
			//set to released state
			dp.released = true;
			//remove polymerase from myDNAComplexes array.
			for(int k=dp.stop;k<=dp.start;k++){
				myPhage.myDNAComplexes[k] = null;
			}
		}
		
	}
    
	//called to terminate an mRNA
	/**
	 * The method to update the state of a DPComplex.  This method is where most of the logic of transcription is hard-coded into the simulation  
	 * @param dp The DPComplex that is to be updated.
	 */	
    void move(DPComplex dp) throws IOException{

		if(dp.whichReaction==2){
            //this is if the complex is going to immediately fall off
            //remove complex and add polymerase
            DPComplexReleaseRxn r = dp.getReleaseRxn();
            r.execute(dp);
            for(int i=0;i<(r.updateRxn).size();i++){
                int reacid = ((Integer) r.updateRxn.elementAt(i)).intValue();
                if(reacid!=r.getID())
                    myCell.mySim.gint.myIPQ.updateOther(reacid, myCell.mySim.gint.t);
            }
            //clear myDNAComplexes vector
            for(int k=dp.stop;k<=dp.start;k++){
                myPhage.myDNAComplexes[k] = null;
            }
            
           
            
            //read promoter complex
			if(dp.myPromoter.stillBlocked>0){
				dp.myPromoter.stillBlocked--;
				
			}else{
				if(dp.myPromoter.copyNumber==0){
					Reaction er = dp.myPromoter.getEntryReaction();
					er.execute();
					for(int i=0;i<(er.updateRxn).size();i++){
						int reacid = ((Integer) er.updateRxn.elementAt(i)).intValue();
						if(reacid!=er.getID())
							myCell.mySim.gint.myIPQ.updateOther(reacid, myCell.mySim.gint.t);
					}
				}else{
					System.out.println("ARG: Reaction uncovering dna... copy # = " + dp.myPromoter.getCopyNumber());
					System.out.println(dp.myPromoter);
					System.exit(9);
				}
			}
            return;
        }else if(dp.whichReaction==3){
            //this occurs if the complex is going towards initiation phase
             //update dp.whichReaction
            dp.initState=1;
            
            
            return;
        }else if(dp.whichReaction==4){
            //this occurs when the complex is falls off and recycles to open complex
            
            
            dp.initState = 0;
            return;
        }else if(dp.whichReaction==5){
            //this occurs if the compex has successfully initiated and is about to continue to make mRNA
            
            dp.myPromoter.incrementInitiated();
			
            dp.whichReaction = 1;
            dp.initState = 3;
            for(int i=0;i<13;i++){
                if(!dp.released)
                    this.move(dp);
            }
            return;
        }
        else{
            //check if it's blocked by another complex
            if(myPhage.myDNAComplexes[dp.start+1]!=null){
                DPComplex dpDownstream = myPhage.myDNAComplexes[dp.start+1];
                int downstreamProteinIndex =  dpDownstream.myProtein.getID() - myPhage.myCell.mySim.FIRSTPOLYMERASEID;
                int currentProteinIndex = dp.myProtein.getID() - myPhage.myCell.mySim.FIRSTPOLYMERASEID;
                int interactionindex = myPhage.myCell.polymeraseInteractionMatrix[currentProteinIndex][downstreamProteinIndex];
                if(interactionindex==0){
                    return;
                } 
				if(interactionindex==-1){
					//signal downstream polymerase to fall off
					DPComplex downstream_dp = myPhage.myDNAComplexes[dp.start+1];
					terminateDPComplex(downstream_dp);
				}
				if(interactionindex==1){
                    //terminate mrna and dp complex
					terminateDPComplex(dp);
					//exit move method
                    return;
                }
            }//end: check if polymerase is blocked

            //actually move
            myPhage.myDNAComplexes[dp.start+1] = myPhage.myDNAComplexes[dp.start];
            dp.start++;
            myPhage.myDNAComplexes[dp.stop]=null;
            dp.stop++;
			myPhage.myCell.transcriptionEvents[dp.myProtein.getID() - myPhage.myCell.mySim.FIRSTPOLYMERASEID]++;
            
            //check if we've caused entry
            if(myPhage.enteredDNAindex<dp.start){
                myPhage.enteredDNAindex = dp.start;
            }

            //check Covering of Elements
            if(myPhage.myDNAFeatures[dp.start]!=0&&myPhage.myDNAFeatures[dp.start-1]!=myPhage.myDNAFeatures[dp.start]){

                DNA d = myPhage.convertFeatureToElement(dp.start);
                
                if(d.getCopyNumber()>0){
                    Reaction r = d.getLossReaction();
                    r.execute();
                    //for entry events only
                    for(int i=0;i<(r.updateRxn).size();i++){
                        int reacid = ((Integer) r.updateRxn.elementAt(i)).intValue();
                        if(reacid!=r.getID())
                            myCell.mySim.gint.myIPQ.updateOther(reacid, myCell.mySim.gint.t);
                    }
                }else{
                    if(d.hasEntered){
                        d.stillBlocked++;
                    }else{
                        d.hasEntered = true;
                    }
                }
            }

            //check Uncovering of Elements
            if(myPhage.myDNAFeatures[dp.stop-1]!=0&&myPhage.myDNAFeatures[dp.stop]!=myPhage.myDNAFeatures[dp.stop-1]){
                DNA d = myPhage.convertFeatureToElement(dp.stop-1);
                
                if(d.stillBlocked>0){
                    d.stillBlocked--;
                }else{
                    if(d.getCopyNumber()==0){
                        Reaction r = d.getEntryReaction();
                        r.execute();
                        //for entry events only
                        for(int i=0;i<(r.updateRxn).size();i++){
                            int reacid = ((Integer) r.updateRxn.elementAt(i)).intValue();
                            if(reacid!=r.getID())
                                myCell.mySim.gint.myIPQ.updateOther(reacid, myCell.mySim.gint.t);
                        }
                    }else{
                        System.out.println("ARG: Reaction uncovering dna... copy # = " + d.getCopyNumber());
                        System.out.println(d);
                        System.exit(9);
                    }
                }
            }

            //check Termination
            if(myPhage.myDNAFeatures[dp.start]>7000&&myPhage.myDNAFeatures[dp.start]<8000){

                Terminator term = (Terminator) myPhage.convertFeatureToElement(dp.start);
                
                if(term.getStopSite()==dp.start){
                    boolean samePolID = false;
                    int[] polIDarry = term.getPolID();
                    double termEff = 0;
                    for(int i=0;i<polIDarry.length;i++){
                        if(polIDarry[i]==dp.myProtein.getID()){
                            samePolID = true;
                            double[] eff = term.getEfficiency();
                            termEff = eff[i];
                            break;
                        }
                    }
                    if(samePolID){
                        if(myCell.mySim.gint.myIPQ.getRandomDouble()<termEff){
							terminateDPComplex(dp);							
                        }	
                    }
                }
            }


        

            //check falloff
            //no need to check for nascent RNA's
            if(dp.start>=myPhage.GENOME_LENGTH){
				terminateDPComplex(dp);
            }

            //update mrna's
            ++dp.m.stop;

            //now check for mRNA elements entry
            int eri = dp.m.stop;
            if(myPhage.myRNAFeatures[eri]!=0){
                if(myPhage.myRNAFeatures[eri+1]!=myPhage.myRNAFeatures[eri]){
                    if(myPhage.myRNAFeatures[eri]>0){
                        //hit a translation start
                        MatureRBS rbs = (MatureRBS) myPhage.myRnaElements[(myPhage.myRNAFeatures[eri]*2)-1];
                        NascentRBSEntryRxn r = rbs.nrbs.entryReaction;
                        r.execute(dp);
                        for(int i=0;i<(r.updateRxn).size();i++){
                            int reacid = ((Integer) r.updateRxn.elementAt(i)).intValue();
                            if(reacid!=r.getID())
                                myCell.mySim.gint.myIPQ.updateOther(reacid, myCell.mySim.gint.t);
                        }
                    }else{
                        //hit a translation stop
                        //make sure this translation stop is for the particular translation start
                        //remove stalled ribosomes
                        TranslationStop ts = (TranslationStop) myPhage.myRnaElements[(myPhage.myRNAFeatures[eri]*-2)-2];
                        for(int i=0;i<dp.stalledRibosomes.size();i++){
                            RRComplex rr = (RRComplex) dp.stalledRibosomes.elementAt(i);
                            if(ts==rr.rbs.mrbs.ts){
                                mtr.rrexecute(rr);
                                myPhage.stalled_rib--;
                                dp.stalledRibosomes.remove(i);
                                i--;
                            }
                        }
						
                        //update NascentRBS and make it a matureRBS
                        //we must see if the rbs is free or not
                        //go through the nrbs's rbsinstances and see if where the dp complexes match... if they don't we have a problem
                        boolean anyRBS = false;
                        boolean rbsIsFree = true;
                        RBSInstance rbsi = null;
                        for(int i=0;i<ts.mrbs.nrbs.myRBSInstances.size();i++){
                            rbsi = (RBSInstance) ts.mrbs.nrbs.myRBSInstances.elementAt(i);
                            if(dp==rbsi.dp){
                                anyRBS = true;
                                if(rbsi.isBound==true){
                                    rbsIsFree = false;
                                }
                                break;
                            }
                        }					
						
                        if(anyRBS){
                            if(rbsIsFree){ // The RBS is free...
                                Reaction r = ts.mrbs.conversionReaction;
                                r.execute();
                                for(int i=0;i<(r.updateRxn).size();i++){
                                    int reacid = ((Integer) r.updateRxn.elementAt(i)).intValue();
                                    if(reacid!=r.getID())
                                        myCell.mySim.gint.myIPQ.updateOther(reacid, myCell.mySim.gint.t);
                                }
                                //if i get an null pointer.. i also got a problem
                                ts.mrbs.nrbs.myRBSInstances.remove(rbsi);
                            }else{  //The RBS is bound.
                                rbsi.rr.totallyentered = true;
                            }
                        }
						
						for(int i=1;i<=mtr.myPQ.size();i++){
							RRComplex rr_test = (RRComplex) mtr.myPQ.elementAt(i);
							if(rr_test.rbs==ts.mrbs.nrbs){
								rr_test.totallyentered = true;
							}
						}
						
                    }
                }
            }
        }
    }
}

/**
 * A class that is used to cause polymerase independent entry of the genome into the cell.
 * @version 1.0
 * @author Sriram Kosuri
 */
class GP16EntryEvent extends Reaction {
	/**
	 * The DNA molecule that this reaction acts on.
	 */
    protected Phage myPhage;
	/**
	 * A placeholder array that is empty.
	 */
    protected Molecule[] affdep;
	/**
     *  The position on the DNA that this reaction turns off.  This constant is not the length of the entire DNA molecule if only a portion is of the DNA enters through this kind of mechanism.
	 */
    int offindex;
    
	/**
	 * The only constructor for this reaction.  
	 * @param ID The ID of this reaction.
	 * @param rc The rate constant that describes the rate of entry of the DNA into the cell
	 * @param offindex The position of the DNA at which this entry mechanisms stops
	 * @param myPhage The DNA that this reaction causes to enter the cell
	 */
    public GP16EntryEvent(int ID,double rc,int offindex,Phage myPhage){
        super(ID,null,null,rc);
        this.myPhage = myPhage;
        updateRxn.add(new Integer(ID));
        affdep = new Molecule[0];
        this.offindex = offindex;
    }

    
	/**
	 * The method to update the entered DNA index upon execution of this reaction.
	 */
    public void execute(){
        int edi = myPhage.enteredDNAindex++;
        if(myPhage.myDNAFeatures[edi]>0){
            if(myPhage.myDNAFeatures[edi+1]==0){
                DNA d = myPhage.convertFeatureToElement(edi);
                
                if(d.hasEntered==true){
                }else{
                    d.hasEntered=true;
                    Reaction r = d.getEntryReaction();
                    r.execute();
                    //for entry events only
                    for(int i=0;i<(r.updateRxn).size();i++){
                        int rid = ((Integer) r.updateRxn.elementAt(i)).intValue();
                        if(rid!=r.getID())
                            myPhage.myCell.mySim.gint.myIPQ.updateOther(rid, myPhage.myCell.mySim.gint.t);
                    }
                }
            }
        }
    }
	/**
	 * Sets the rate constant to either the defined #rateConstant or to zero if the DNA has entered beyond the @offIndex
	 */ 
    public void calculateA(){
        if(myPhage.enteredDNAindex<offindex)
            a = rateConstant;
        else
            a=0;
    }
	/**
     * Overrides the inherited method so that this method returns an empty array.
	 * @return An empty array.
	 */		
    public Molecule[] affects() {
        return affdep;
    }

	/**
	 * Overrides the inherited method so that this method returns an empty array.
	 * @return An empty array.
	 */		
    public Molecule[] dependsOn() {
        return affdep;
    }
    
}
/**
 * A class used to instantiate reactions that have no reactants, i.e., reactions that generate something from nothing.  This reaction is usually called within another reaction to create some molecule, e.g., a protein upon translation.
 * @version 1.0
 * @author Sriram Kosuri
 */
class ZeroOrderRxn extends Reaction {
	/**
	 * A constructor that takes in reactant and product molecules as well.
	 * @param ID the reaction ID
	 * @param prod The product to be incremented upon execution of this reaction.
	 * @param rc The mesoscopic rate constant of the reaction.
	 */
    public ZeroOrderRxn(int ID, Molecule prod, double rc){
        super(ID,rc);
        Reactants = new Molecule[0];
        Products = new Molecule[1];
        Products[0] = prod;
    }
	
	/**
	 * Overrides the parent method to only return the generated product
	 * @return an array of size one containing the generated product.
	 */
    public Molecule[] affects(){
        return Products;
    }
}

class FirstOrderRxn0 extends Reaction {
    
    public FirstOrderRxn0(int ID, Molecule reac, double rc){
        super(ID,rc);
        Reactants = new Molecule[] {reac};
        Products = new Molecule[0];
    }

    public Molecule[] affects(){
        return Reactants;
    }

}

class FirstOrderRxn1 extends Reaction {

    public FirstOrderRxn1(int ID, Molecule reac, Molecule prod, double rc){
        super(ID,rc);
        Reactants = new Molecule[] {reac};
        Products = new Molecule[] {prod};
    }

}

class FirstOrderRxn2 extends Reaction {
    public FirstOrderRxn2(int ID, Molecule reac, Molecule prod1, Molecule prod2, double rc){
        super(ID,rc);
        Reactants = new Molecule[] {reac};
        Products = new Molecule[] {prod1, prod2};
    }
}

class SecondOrderRxn2 extends Reaction {
    public SecondOrderRxn2(int ID, Molecule reac1, Molecule reac2, Molecule prod1, Molecule prod2, double rc){
        super(ID,rc);
        Reactants = new Molecule[] {reac1,reac2};
        Products = new Molecule[] {prod1, prod2};
    }
	
	public void calculateA(){
        a = rateConstant;
		
        if ((Reactants[0].getCopyNumber()!=0)&&(Reactants[0]==Reactants[1])){    
			a = a*Reactants[0].getCopyNumber()*(Reactants[1].getCopyNumber()-1);
        }else{
			a = a*Reactants[0].getCopyNumber()*(Reactants[1].getCopyNumber());
        }    
        
	}
	
}


class SecondOrderRxn extends Reaction {
    public SecondOrderRxn(int ID, Molecule reac1, Molecule reac2, Molecule prod, double rc){
        super(ID,rc);
        Reactants = new Molecule[] {reac1,reac2};
        Products = new Molecule[] {prod};
    }
    
      public void calculateA(){
        a = rateConstant;
       
        if ((Reactants[0].getCopyNumber()!=0)&&(Reactants[0]==Reactants[1])){    
        a = a*Reactants[0].getCopyNumber()*(Reactants[1].getCopyNumber()-1);
        }else{
        a = a*Reactants[0].getCopyNumber()*(Reactants[1].getCopyNumber());
        }    
        
      }
}

class RBSConversionRxn extends Reaction {
    public RBSConversionRxn(int ID, Molecule reac, Molecule prod, double rc){
        super(ID,rc);
        Reactants = new Molecule[] {reac};
        Products = new Molecule[] {prod};
        a = 0;
        rc = 0;
    }

    public void calculateA(){};
}

class ElongationMRBSRxn extends Reaction implements GammaJumpRxn{
    Cell myCell;
    int elonglength;
    PriorityQueue myPQ;
    boolean initialtime;
	Protein prot;
    

    public ElongationMRBSRxn(int ID,Molecule mrnaelong,Protein prod, Ribosome rib,double rc,int elonglength,Cell myCell){
        super(ID,rc);
        this.elonglength = elonglength;
        this.myCell = myCell;
        a = 0;
		prot = prod;
        Reactants = new Molecule[] {mrnaelong};
        Products = new Molecule[] {prod,rib};
        time = Double.POSITIVE_INFINITY;
		myPQ = new PriorityQueue();
    }

    public void calculateA(){

    }

	
    public void incremented(int inc){
        try{
            if(inc==1){
                double newtime = myCell.mySim.gint.myIPQ.getGammaTime(rateConstant,elonglength)+myCell.mySim.gint.t;
                myPQ.insert(new Double(newtime));
                time = ((Double) myPQ.findMin()).doubleValue();
            }
        }catch(Exception e){
            e.printStackTrace();
            System.exit(9);
        }
    }
    
    
    public void calculateTime(double random, double t){

    }

    public void execute(){
        super.execute();
		myCell.moleculeEnergy[prot.getID()-1] += prot.length;
        try{
            myPQ.deleteMin();
            try{
                time = ((Double) myPQ.findMin()).doubleValue();
            }catch(java.io.IOException e){
                time = Double.POSITIVE_INFINITY;
            }
        }catch(Exception e){
            e.printStackTrace();
            System.exit(9);
        }
    }
    
}

class InitiationMRBSDegRxn extends FirstOrderRxn1{
	InitiationMRBSRxn initmrbsrxn;
	
	public InitiationMRBSDegRxn(int ID, Molecule mrnainit, Molecule rib, double rc, InitiationMRBSRxn initmrbsrxn){
		super(ID,mrnainit,rib,rc);
		this.initmrbsrxn = initmrbsrxn;
	}
	
	public void execute(){
		super.execute();
		initmrbsrxn.removeOne();
	}
}

class InitiationMRBSRxn extends Reaction implements GammaJumpRxn{
    Cell myCell;
    int initlength;
    PriorityQueue myPQ;
    GammaJumpRNA mrnainit;
	GammaJumpRNA mrnaelong; 
	double initRateConstant;
	
    public InitiationMRBSRxn(int ID,Molecule mrnainit, Molecule mrnaelong,MatureRBS mrbs,double rc, int initlength,Cell myCell, double initRateConstant){
        super(ID,rc);
        this.mrnainit = (GammaJumpRNA) mrnainit;
		this.initlength = initlength;
        this.myCell = myCell;
		this.mrnaelong = (GammaJumpRNA) mrnaelong;
        a=0;
        Reactants = new Molecule[] {mrnainit};
        Products = new Molecule[] {mrnaelong, mrbs};
        time = Double.POSITIVE_INFINITY;
		this.initRateConstant = initRateConstant;
		myPQ = new PriorityQueue();
    }

    public void calculateA(){

    }
	
	public void removeOne(){
		try{
			if(myPQ.size()==0){
				System.out.println("Deleting1 mrnainit " +mrnainit.getID() + " that is currently at " + mrnainit.getCopyNumber() + " at time " + myCell.mySim.gint.t);
				System.out.println("Trying to remove an initRBS that doesn't exist");
				System.out.println(myPQ.size() + " =? " + mrnainit.getCopyNumber());
				System.exit(9);
			}else if(myPQ.size()==1){
				myPQ.deleteMin();
				time = Double.POSITIVE_INFINITY;
			}else{
				int index = (int) Math.floor(myPQ.size()*myCell.mySim.gint.myIPQ.getRandomDouble()) + 1;
				if(index==1){
					myPQ.deleteMin();
					time = ((InitiationEntry) myPQ.findMin()).time;
				}else{
					myPQ.delete(index);
					time = ((InitiationEntry) myPQ.findMin()).time;
				}
			}
		}catch(Exception e){
            e.printStackTrace();
            System.exit(9);
        }			
	}
	
    public void incremented(int inc){
        try{
            if(inc==1){
				double newtime = (1/initRateConstant)*Math.log(1/myCell.mySim.gint.myIPQ.getRandomDouble()) + myCell.mySim.gint.t;;
                myPQ.insert(new InitiationEntry(newtime,true));
                time = ((InitiationEntry) myPQ.findMin()).time;
            }   
		}catch(Exception e){
            e.printStackTrace();
            System.exit(9);
        }
    }
    
    public void calculateTime(double random, double t){
		try{
			if(myPQ.size()>0)
				time = ((InitiationEntry) myPQ.findMin()).time;
			else
				time = Double.POSITIVE_INFINITY;
		}catch(IOException e){
			e.printStackTrace();
			System.exit(9);
		}
	}

    public void execute(){
        try{
            InitiationEntry ie = (InitiationEntry) myPQ.deleteMin();
			if(ie.starting_off){
				ie.starting_off = false;
                ie.time = myCell.mySim.gint.myIPQ.getGammaTime(rateConstant,initlength)+myCell.mySim.gint.t;
				myPQ.insert(ie);
			}else{
				super.execute();
			}
            if(myPQ.size()>0){
                time = ((InitiationEntry) myPQ.findMin()).time;
            }else{
                time = Double.POSITIVE_INFINITY;
            }
        }catch(Exception e){
            e.printStackTrace();
            System.exit(9);
        }
    }
}

class InitiationEntry implements Comparable{
	double time;
	boolean starting_off; // starting_off = true if initiating,   = false when clearing while elongating 
	
	public InitiationEntry(double time, boolean starting_off){
		this.time = time;
		this.starting_off = starting_off;
	}
	
	public int compareTo(Object r){
        InitiationEntry rhs = (InitiationEntry) r;
        if((this.time - rhs.time)<0)
            return -1;
        else if((this.time - rhs.time)>0)
            return 1;
        else
            return 0;
    }
	
}

class MatureRBSBindingRxn extends Reaction {
    public MatureRBSBindingRxn(int ID,Ribosome rib,MatureRBS mrbs,Molecule mrnainit, double rc){
        super(ID,rc);
        Reactants = new Molecule[] {rib,mrbs};
        Products = new Molecule[] {mrnainit};
    }

    public void execute(){
        super.execute();
    }

    public void calculateTime(double random, double t){
        super.calculateTime(random,t);
    }
}

class NascentRBSEntryRxn extends ZeroOrderRxn{
    public NascentRBSEntryRxn(int ID, NascentRBS nrbs, double rc){
        super(ID,nrbs,rc);
    }

    public void execute(DPComplex dp){
        super.execute();
        ((NascentRBS) Products[0]).myRBSInstances.add(new RBSInstance(dp));
    }
}

class NascentRBSLossRxn extends FirstOrderRxn0{
    public NascentRBSLossRxn(int ID, NascentRBS nrbs, double rc){
        super(ID,nrbs,rc);
    }

    public void execute(){
        super.execute();
    }
}

//binding reaction for nascent mrna's
class NascentRBSBindingRxn extends Reaction {

    MrnaTrackerRxn mtr;

    public NascentRBSBindingRxn(int ID, Ribosome rib, NascentRBS nrbs, double rc, MrnaTrackerRxn mtr){
        super(ID, rc);
        Reactants = new Molecule[] {rib,nrbs};
        Products = new Molecule[0];
        this.mtr = mtr;
    }
    
    public void execute(){
        //update Reactants
        //this order matters... b/c of adding rbsinstance.. must pick before nrbs quantity is updated
        mtr.rbsBinding((NascentRBS) Reactants[1]);
        for(int i=0;i<Reactants.length;i++){
            Reactants[i].incrementCopyNumber(-1);
        }
    }

    public Molecule[] affects(){
        return Reactants;
    }
    
}

class PromoterBindingRxn extends Reaction {
    RNAPElongationEvent myElongRxn;
    Phage myPhage;
    double speed;
    // This was added by Jason Kelly 7/18/2003 as a check for determining the Ainit,on and off rates 
    // from abortive inititation %
    static int complexID;
    // end addition
    
    private PromoterBindingRxn(int ID, Molecule pol, Molecule dna, double rc, RNAPElongationEvent e, Phage myPhage, double speed){
        super(ID,rc);
        Reactants = new Molecule[] {pol, dna};
        Products = new Molecule[0];
        myElongRxn = e;
        updateRxn.add(new Integer(e.getID()));
        this.myPhage = myPhage;
        this.speed = speed;
        complexID = 1;
    }

    private PromoterBindingRxn(int ID, Molecule pol, Molecule dna, RNAPElongationEvent e, Phage myPhage, double speed){
        this(ID,pol,dna,((Promoter) dna).aon[pol.getID()-1],e,myPhage,speed);
    }

    public PromoterBindingRxn(int ID, Molecule pol, Molecule dna, RNAPElongationEvent e, Phage myPhage){
        this(ID,pol,dna,((Promoter) dna).aon[pol.getID()-1],e,myPhage,((Polymerase) pol).speed);        
    }
    
    public void execute(){
        //update reactants
        for(int i=0;i<Reactants.length;i++){
            Reactants[i].incrementCopyNumber(-1);
        }
       //Added by Jason 7/18/2003
        Polymerase pol = (Polymerase) Reactants[0];
        String dpid = pol.getOrganism() + "-" + pol.getID() + "-" + complexID%10000;
        
        
        DPComplex dp = new DPComplex(dpid, pol,(Promoter) Reactants[1],((DNA) Reactants[1]).stop, ((DNA) Reactants[1]).stop-pol.footprint, speed, myPhage.myCell.mySim.gint.myIPQ.getRandomDouble(), myPhage.myCell.mySim.gint.t,myPhage);
        myElongRxn.addComplex(dp);
        complexID++;
        
        //Create mRNA
        int startmrna = ((Promoter) Reactants[1]).getStartSite();
        dp.m = new NascentRNA(startmrna, startmrna, dp);
        //createMrna()
    }

        public Molecule[] affects(){
        return Reactants;
    }
    
}

class DPComplexReleaseRxn extends Reaction {
    RNAPElongationEvent myElongRxn;
    
    public DPComplexReleaseRxn(int ID, Molecule pol, double rc, RNAPElongationEvent e){
        super(ID,rc);
        myElongRxn = e;
        Reactants = new Molecule[0];
        Products = new Molecule[] {pol};
    }

    public void execute(){
        System.out.println("Should never get to this line in DPComplexReleaseRxn's execute() function");
        System.exit(9);
    }
    
    //This method should be sufficient
    public void execute(DPComplex dp){
        try{
            //update products
            for(int i=0;i<Products.length;i++){
                Products[i].incrementCopyNumber(1);
            }
            //update RNAP elongation event queue
			for(int i=1;i<=myElongRxn.myPQ.size();i++){
				DPComplex dptest = (DPComplex) myElongRxn.myPQ.elementAt(i);
				if(dptest==dp){
					if(i==1){
						myElongRxn.myPQ.deleteMin();
					}else{
						myElongRxn.myPQ.delete(i);
					}
					return;
				}
			}
			
			//myElongRxn.myCell.mySim.gint.updateSelf()
			System.out.println("DIDN'T REMOVE DP COMPLEX... OOPS..");
			System.exit(9);
			
        }catch(Exception e){
            e.printStackTrace();
            System.exit(9);
        }
    }

    public Molecule[] affects(){
        return Products;
    }
    
}

class MrnaTrackerRxn extends Reaction {

    Phage myPhage;
    //Vector freeRBS;
    PriorityQueue myPQ;
    Ribosome rib;
    Reaction nascentRBSEntryRxn;
    Cell myCell;
    
    public MrnaTrackerRxn(int ID, double rc, Ribosome rib, Phage myPhage, Cell myCell){
        super(ID, rc);
        this.myPhage = myPhage;
        myPQ= new PriorityQueue();
        this.rib = rib;
        Reactants = new Molecule[0];
        Products = new Molecule[0];
        updateRxn.add(new Integer(ID));
        myPhage.nrbs_PQ = this.myPQ;
        this.myCell = myCell;
    }

    public void calculateA(){
        
    }

    
    public void rbsBinding(NascentRBS nrbs){
        RRComplex rric = new RRComplex(nrbs,rib.speed,nrbs.mrbs.initstepsize,nrbs.mrbs.elongstepsize);
        rric.time = (1/nrbs.mrbs.initRateConstant)*Math.log(1/myCell.mySim.gint.myIPQ.getRandomDouble()) + myCell.mySim.gint.t;

        //select a particular rbs binding --
        int index = (int) Math.floor(nrbs.copyNumber*myCell.mySim.gint.myIPQ.getRandomDouble());
        int k = 0;
        RBSInstance rbsi = null;
        //this should work no matter the size of nrbs... if it throws an exception here.. something is wrong
        for(int i=0;k<=index;i++){
            rbsi = ((RBSInstance) nrbs.myRBSInstances.elementAt(i));
            if(!rbsi.isBound){
                k++;
            }
        }
        //now add
        rbsi.isBound = true;
        rbsi.rr = rric;
        //now give rr an mrna
        rbsi.rr.nrna = rbsi.dp.m;
        myPQ.insert(rric);
		myCell.mySim.gint.myIPQ.updateSelf(this,myCell.mySim.gint.t);
    }

    public void calculateTime(double random, double t){
        try{
            //find minimum
            RRComplex rr = (RRComplex) myPQ.findMin();
            time = rr.time;
        }catch(java.io.IOException e){
            time = Double.POSITIVE_INFINITY;
        }
    }

    public void execute(){
        try{
            //find minimum
            RRComplex rr = (RRComplex) myPQ.deleteMin();
			
            //check if it is pre-init
			if(rr.pre_init){
				rr.pre_init = false;
				rr.time = myCell.mySim.gint.myIPQ.getGammaTime(rr.a,rr.initstepsize)+myCell.mySim.gint.t;
				myPQ.insert(rr);
				return;
			}
			
			//execute reaction
            //if it's initiation.. delete from queue and then add elongation to queue, return NascentRBS
            if(rr.init){
                if(!rr.totallyentered){
                    //change boolean markers
                    rr.init=false;
                    //free rbs

                    Reaction r = rr.rbs.entryReaction;
                    r.execute();
                    for(int i=0;i<(r.updateRxn).size();i++){
                        int reacid = ((Integer) r.updateRxn.elementAt(i)).intValue();
                        if(reacid!=r.getID())
                            myCell.mySim.gint.myIPQ.updateOther(reacid, myCell.mySim.gint.t);
                    }

                    //update rbsi
                    for(int i=0;i<rr.rbs.myRBSInstances.size();i++){
                        RBSInstance rbsi = (RBSInstance) rr.rbs.myRBSInstances.elementAt(i);
                        if(rbsi.rr!=null&&rbsi.rr==rr){
                            rbsi.rr = null;
                            rbsi.isBound = false;
                            break;
                        }
                    }
                    
                    //recalculate time
                    rr.time = myCell.mySim.gint.myIPQ.getGammaTime(rr.a,rr.elongstepsize)+myCell.mySim.gint.t;
                    //insert into priority queue
                    myPQ.insert(rr);
                }else{
                    //totally entered
                    rr.init = false;

                    //remove rbsi
                    RBSInstance rbsi = null;
                    for(int i=0;i<rr.rbs.myRBSInstances.size();i++){
                        rbsi = (RBSInstance) rr.rbs.myRBSInstances.elementAt(i);
                        if(rbsi.rr!=null&&rbsi.rr==rr){
                            break;
                        }
                    }

                    //if i get an null pointer.. i've got a problem
                    rr.rbs.myRBSInstances.remove(rbsi);

                    //free Mature RB
                    MatureRBS mrbs = rr.rbs.mrbs;
                    Reaction r = mrbs.entryReaction;
                    r.execute();
                    for(int i=0;i<(r.updateRxn).size();i++){
                        int reacid = ((Integer) r.updateRxn.elementAt(i)).intValue();
                        if(reacid!=r.getID())
                            myCell.mySim.gint.myIPQ.updateOther(reacid, myCell.mySim.gint.t);
                    }
					//just add an mrna elong to react
					((InitiationMRBSRxn) mrbs.initrna.rxnToUpdate).mrnaelong.incrementCopyNumber(1);
                }
            }else{
				//!rr.init
                //check if polymerase has elongated
                if(rr.totallyentered){
                    rrexecute(rr);
                }else{
                    //this means polymerase hasn't gotten to where it wants to go yet...
                    //so this part should put the complex into a queue waiting to be popped out.
                    int ribspacing = (rr.nrna.dp.stalledRibosomes.size()+1)*this.rib.footprint;
                    int ribspace = rr.nrna.dp.start - rr.rbs.mrbs.stop;
                    if(ribspacing < ribspace){
                        //add to queue
                        rr.nrna.dp.stalledRibosomes.add(rr);
                        myPhage.stalled_rib++;
                    }else{
                        //free ribosome
                        Reaction r = rib.gainReaction;
                        r.execute();
                        for(int i=0;i<(r.updateRxn).size();i++){
                            int reacid = ((Integer) r.updateRxn.elementAt(i)).intValue();
                            if(reacid!=r.getID())
                                myCell.mySim.gint.myIPQ.updateOther(reacid, myCell.mySim.gint.t);
                        }
                    }
                }
            }
			
			
			
        }catch(Exception e){
            e.printStackTrace();
            System.exit(9);
        }
    }

    public void rrexecute(RRComplex rr){
        //free ribosome
        Reaction r = rib.gainReaction;
        r.execute();
        for(int i=0;i<(r.updateRxn).size();i++){
            int reacid = ((Integer) r.updateRxn.elementAt(i)).intValue();
            if(reacid!=r.getID())
                myCell.mySim.gint.myIPQ.updateOther(reacid, myCell.mySim.gint.t);
        }
        //add protein
        Reaction r2 = rr.rbs.mrbs.ts.product.gainReaction;
        r2.execute();
        for(int i=0;i<(r2.updateRxn).size();i++){
            int reacid = ((Integer) r2.updateRxn.elementAt(i)).intValue();
            if(reacid!=r2.getID())
                myCell.mySim.gint.myIPQ.updateOther(reacid, myCell.mySim.gint.t);
        }
		//add to energy counter
		myCell.moleculeEnergy[rr.rbs.mrbs.ts.product.getID()-1] += rr.rbs.mrbs.ts.product.length;
    }

}	
