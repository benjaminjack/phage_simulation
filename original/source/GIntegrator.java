
import java.io.IOException;
import java.util.Random;

/**
 * Class that orchestrates the Gillespie-based integrator and stores the priority queue.
 *
 * @author Sriram Kosuri
 * @version 1.0
 * @see IndexedPriorityQueue
 */

public class GIntegrator {
	/**
     * The IndexedPriorityQueue where all the reactions, including specialized DNA entries are stored.  This queue sorts
	 * these reactions such that the reaction to be executed next always is in the first position.
	 */
    IndexedPriorityQueue myIPQ;
    /**
	 * The current time of the simulation.
	 */
	double t;
	/**
	 * A tracker to account for the number of reactions that have occurred.
	 */
    int counter;

	/**
	 * The only constructor of GIntegrator.
	 * @param thereacs A vector of reactions that are defined at the start of simulation using TabascoXML.
	 * @param key A key to create the random number generator.
	 */
    public GIntegrator(Reaction[] thereacs, long key){
        myIPQ = new IndexedPriorityQueue(thereacs.length, key);
        t = 0;
        try{
            for(int i=0;i<thereacs.length;i++){
                myIPQ.toss(thereacs[i]);
            }
            myIPQ.buildHeap();
        }catch(IOException e){
            e.printStackTrace();
            System.exit(9);
        }
        counter =0;
    }
    /**
	 * The main loop that is run find the next reaction to execute, update the time, execute the reaction, and reinsert into the priority queue.
	 */
    public void runLoop(){
        //get reaction with minimum tau
        Reaction r = myIPQ.findMin();
		//update time
		if(r.time<t){
			System.out.println("Reaction is trying to go back in time");
			System.exit(9);
		}		
		
        t = r.time;
        //update reactants
        r.execute();
        //Find reactions that this affects and update their times
        for(int i=0;i<(r.updateRxn).size();i++){
            int rid = ((Integer) r.updateRxn.elementAt(i)).intValue();
            if(rid!=r.getID())
                myIPQ.updateOther(rid, t);
            else
                myIPQ.updateSelf(r,t);
        }
        counter++;

    }

}

/**
 * An indexed priority queue as specified by Michael Gibson for the Next Reaction Method.
 *
 * @author Mark Allen Weiss
 * @author Sriram Kosuri
 * @version 1.0
 * @see IndexedPriorityQueue
 */
class IndexedPriorityQueue {
	/**
	 * Inaccessible vector containing the reactions to in a priority queue structure.  Element 1 is the reaction that is to execute next.  In
	 * the future, it may not be necessary to use this class, because the dependency graph uses pointers to the actual reaction, not reaction ID's.
	 * There is no need to index the priority queue.
	 */
    private Reaction[] pq;
	/**
	 * A vector where the position corresponds to the reaction index, and the value at the position corresponds to the 
	 * position in the pq vector.
	 */
    private int[] index;
	/**
	 * Number of reaction sin the pq vector
	 */ 
    private int currentSize;      // Number of elements in heap
    /**
	 * The random number generator used for this simulation.
	 */
	private Random randomNG;

    
    /**
	 * Construct the binary heap.
     * @param capacity the capacity of the binary heap.
     * @param key The key to the random number generator
     */
    public IndexedPriorityQueue( int capacity, long key )
    {
        currentSize = 0;
        pq = new Reaction[ capacity + 1 ];
        index = new int[ capacity + 1 ];
        randomNG = new Random(key);
    }
	
	/**
	 * Returns a pseudo-random uniform double value ranging from [0,1).
     * @return the pseudo-random double ranging from [0,1).
	 */
    public double getRandomDouble(){
        return randomNG.nextDouble();
    }

	/*
	 * Returns a pseudo-random number distributed according to the Gamma distribution.  Uses the rejection method as outlined
	 * in Numerical Recipes in C (Chapter 7, section 3).  
	 * @param ratec The rate constant for an individual exponential reaction event (e.g., one translation elongation event)
	 * @param numofevents The number of exponential events that happen in a row. 
	 * @return a pseudo-random number that represents the amount of time it takes for (numofevents) reactions of rate (ratec)
	 */
    public double getGammaTime(double ratec, int numofevents){
        int j;
        double am,e,s,v1,v2,x,y;
        double mean = 1/ratec;

        if(numofevents < 6){
            x=1.0;
            for (j=1;j<=numofevents;j++)
                x *= getRandomDouble();
            x = -Math.log(x);
        }else {
            do {
                do {
                    do {
                        v1 = getRandomDouble();
                        v2 = 2*getRandomDouble()-1;
                    }while (v1*v1+v2*v2 > 1);
                    y = v2/v1;
                    am=numofevents-1;
                    s=Math.sqrt(2*am+1);
                    x=s*y+am;
                }while (x <=0);
                e=(1+y*y)*Math.exp(am*Math.log(x/am)-s*y);
            }while (getRandomDouble() > e);
        }
        return mean*x;
    }

	
	/**
	 * Calculates the the time for the reaction to occur.  Should be deprecated, as it is only used at first instantiation of the 
	 * reactions.  After such time, you must add the current time of the simulation
	 * @see Reaction#calculateTime(Reaction, double)
	 */ 
	private void calculateTime(Reaction r){
		r.time = (1/r.a)*Math.log(1/getRandomDouble());
    }
    

    /**
     * Insert into the priority queue, maintaining heap order.
     * Duplicates are allowed.  Not used in Tabasco.
     * @param x the item to insert.
     * @exception IOException if container is full.
     */
    private void insert( Reaction x ) throws IOException
    {
        calculateTime(x);
        if( isFull( ) )
            throw new IOException( );
        // Percolate up
        int hole = ++currentSize;
        for( ; hole > 1 && x.compareTo( pq[ hole / 2 ] ) < 0; hole /= 2 ){
            pq[ hole ] = pq[ hole / 2 ];
            index[pq[hole].getID()]=hole;
        }
        pq[ hole ] = x;
        index[pq[hole].getID()]=hole;
        //Update dependency graph
        for(int i=1;i<=currentSize;i++){
            x.findAffects(pq[i]);
        }
    }

    public void toss( Reaction x ) throws IOException
    {
        calculateTime(x);
        if(isFull()){
            throw new IOException();
        }
        int hole = ++currentSize;
        pq[hole] = x;
        index[pq[hole].getID()] = hole;
        //Update dependency graph
        for(int i=1;i<=currentSize;i++){
            x.findAffects(pq[i]);
        }
    }

    /**
        * Find the smallest item in the priority queue.
     * @return the smallest item, or null, if empty.
     */
    public Reaction findMin( )
    {
        if( isEmpty( ) )
            return null;
        return pq[ 1 ];
    }

    /**
        * Remove the smallest item from the priority queue.
     * @return the smallest item, or null, if empty.
     */
    public Reaction deleteMin( )
    {
        if( isEmpty( ) )
            return null;
        Reaction minItem = findMin( );
        pq[ 1 ] = pq[ currentSize-- ];
        index[minItem.getID()] = -1;
        percolateDown( 1 );
        return minItem;
    }

    /**
        * Establish heap order property from an arbitrary
     * arrangement of items. Runs in linear time.
     */
    public void buildHeap( )
    {
        for( int i = currentSize / 2; i > 0; i-- )
            percolateDown( i );
    }

    /**
        * Test if the priority queue is logically empty.
     * @return true if empty, false otherwise.
     */
    public boolean isEmpty( )
    {
        return currentSize == 0;
    }

    /**
        * Test if the priority queue is logically full.
     * @return true if full, false otherwise.
     */
    public boolean isFull( )
    {
        return currentSize == pq.length - 1;
    }
    /**
        * Make the priority queue logically empty.
     */
    public void makeEmpty( )
    {
        currentSize = 0;
    }


    /**
        *Swaps tree nodes
     * @param i swaps tree node with index i with j
     * @param j swaps tree node with index j with i
     */
    public void swap(int i, int j){
        Reaction tmp = pq[i];
        //update position i
        pq[i] = pq[j];
        index[pq[i].getID()] = i;
        //update position j
        pq[j] = tmp;
        index[pq[j].getID()] = j;
    }

	/*
	 * This method is used for updating a reaction (with ID=ID), upon execution of another reaction, and reordering the pq.
	 * @param ID ID of reaction to be updated
	 * @param t the current time of the simulation
	 */
	
    public void updateOther(int ID, double t){
        int n = index[ID];
        double aold = pq[n].a;
        if(aold==0){
            Reaction r = pq[n];
            r.calculateA();
            r.calculateTime(getRandomDouble(),t);
        }else{
            pq[n].calculateA();
            double anew = pq[n].a;
            if(anew==0){
                pq[n].time = Double.POSITIVE_INFINITY;
            }else{
                double newtime = (aold/anew)*(pq[n].time - t) + t;
                pq[n].time = newtime;
            }
        }
        update_aux(n);
    }
	/*
	 * This method is used for updating a reaction (r) upon its own execution, and reordering the pq.
	 * @param r the reaction to be updated.
	 * @param t the current time of the simulation.
	 */
	
    public void updateSelf(Reaction r, double t){
        r.calculateA();
        r.calculateTime(getRandomDouble(),t);
        update_aux(index[r.getID()]);
    }

	/*
	 * A helper method to check what the minimum of two children are from any particular node on the binary tree. 
	 * @param n The index of the parent in pq.
	 * @return the index in pq of the smallest child of entry n.  Returns -1 if there are no children.
	 */
    protected int getMinimumChildIndex(int n){
        if(n*2<=currentSize){
            if((n*2+1)<=currentSize){
                if((pq[n*2].time)<(pq[n*2+1].time)){
                    return n*2;
                }else{
                    return (n*2+1);
                }
            }else{
                return n*2;
            }
        }else{
            return -1;
        }
    }
	/*
	 * A helper method to update the position of a particular reaction within pq without updating its current time.
	 * @param ID the ID of the reaction to update the position of.
	 */
    public void justUpdate(int ID){
        update_aux(index[ID]);
    }
    
	/*
	 * A recursive helper method to update the position of an entry n within the pq.  Will check the nodes children and parents.
	 * @param n the reaction to update.
	 */
    private void update_aux(int n){
        if((n/2)>0&&(pq[n].time < pq[n/2].time)){
            swap(n, n/2);
            update_aux(n/2);
            return;
        }
        int min = getMinimumChildIndex(n);
        if(min!=-1&&(pq[n].time>pq[min].time)){
            swap(n, min);
            update_aux(min);
        }else{
            return;
        }
    }

    /**
        * Internal method to percolate down in the heap.
     * @param hole the index at which the percolate begins.
     */
    private void percolateDown( int hole )
    {
        int child;
        Reaction tmp = pq[ hole ];
        for( ; hole * 2 <= currentSize; hole = child )
        {
            child = hole * 2;
            if( child != currentSize && pq[ child + 1 ].compareTo( pq[ child ] ) < 0 )
                child++;
            if( pq[ child ].compareTo( tmp ) < 0 ){
                pq[ hole ] = pq[ child ];
                index[pq[hole].getID()] = hole;
            }else
                break;
        }
        pq[ hole ] = tmp;
        index[pq[hole].getID()] = hole;
    }
	
	/**
	 * A test method that overides the toString method.  Outputs the reactions that are in the queue and their order.
	 * returns an string describing the contents of the pq. 
	 */
    public String toString(){
        String s1 = ("This is what is contained in the IndexedPriorityQueue\n");
        String s2 = ("The number of Reactions are: " + currentSize+"\n");
        String arry = "";
        String arry2 = "";
        String arrylabel = ("The contents of pq are: \n");
        String arry2label = ("The contents of index are: \n");
        for(int i=1;i<=currentSize;i++){
            arry = arry + "Position "+i+" is: " + pq[i] + "\n";
        }
        for(int i=1;i<=currentSize;i++){
            arry2 = arry2 + "Position "+i+" is: " + index[i] +"\n";
        }
        return s1+s2+arrylabel+arry+arry2label+arry2;
    }

	/**
	 * An internal method to test what reactions depend on what. 
	 */
    private void testDependencies(){
        for(int i=1;i<=currentSize;i++){
            System.out.println("UpdateRxn for Reaction " + i);
            for(int j=0;j<pq[i].updateRxn.size();j++){
                System.out.println(pq[i].updateRxn.elementAt(j));
            }
        }
    }

}

