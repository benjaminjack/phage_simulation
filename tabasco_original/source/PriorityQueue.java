// PriorityQueue class
//
// CONSTRUCTION: empty or with initial array.
//
// ******************PUBLIC OPERATIONS*********************
// void insert( x )       --> Insert x
// Comparable deleteMin( )--> Return and remove smallest item
// Comparable findMin( )  --> Return smallest item
// boolean isEmpty( )     --> Return true if empty; else false
// void makeEmpty( )      --> Remove all items
// ******************ERRORS********************************
// Throws UnderflowException for findMin and deleteMin when empty
/**
 * Implements a binary heap.
 * Note that all "matching" is based on the compareTo method.
 * @author Mark Allen Weiss
 * @author Sriram Kosuri
 */
import java.io.IOException;
public class PriorityQueue {
	/**
	 * The intial capacity of the queue, before it needs to automatically double.
	 */
	private static final int DEFAULT_CAPACITY = 100;
	/**
	 * The size of the priority queue.
	 */
    private int currentSize;      // Number of elements in heap
	/**
	 * The array that contains the elements in the queue.  
	 */
    private Comparable [ ] array; // The heap array

    /**
        * Construct the binary heap.
     */
    public PriorityQueue(){
        currentSize = 0;
        array = new Comparable[ DEFAULT_CAPACITY + 1 ];
    }

    /**
        * Insert into the priority queue.
     * Duplicates are allowed.
     * @param x the item to insert.
     * @return null, signifying that decreaseKey cannot be used.
     */
    public void insert( Comparable x )
    {
        if( currentSize + 1 == array.length )
            doubleArray( );
        // Percolate up
        int hole = ++currentSize;
        array[ 0 ] = x;
        for( ; x.compareTo( array[ hole / 2 ] ) < 0; hole /= 2 ){
            array[ hole ] = array[ hole / 2 ];
		}
        array[ hole ] = x;
    }

    /**
        * Retrieve a particular element in the priority queue
     * @param i the index 
     * @return the Comparable element requested
     * @throws IOException if index is outside valid indices.
     */

    public Comparable elementAt(int i) throws java.io.IOException
    {
        if(i<1||i>currentSize)
            throw new IOException("Index outside array size");
        else
            return array[i];
    }

    
    /**
        * Find the smallest item in the priority queue.
     * @return the smallest item.
     * @throws UnderflowException if empty.
     */
    public Comparable findMin( ) throws java.io.IOException
    {
        if( isEmpty( ) )
            throw new java.io.IOException( "Empty binary heap" );
        return array[ 1 ];
    }

    /**
        * Remove the smallest item from the priority queue.
     * @return the smallest item.
     * @throws UnderflowException if empty.
     */
    public Comparable deleteMin( ) throws java.io.IOException
    {
        Comparable minItem = findMin( );
        array[ 1 ] = array[ currentSize-- ];
        percolateDown( 1 );
        return minItem;
    }

    /**
        * Removes an object from the priority queue
     * @param index the index of the parameter to be removed
     * @return the deleted item
     * @throws IOException if object is the first object
     */
    public Comparable delete(int hole) throws java.io.IOException
    {
        //if it is the first index.. throw an IOException
        if(hole==1){
            throw new java.io.IOException("Trying to remove minimum index: use deleteMin() instead");
        }
        //swap last element into items position and then percolate down
        Comparable rem = array[hole];
		promote(hole);
		deleteMin();
        return rem;
    }
    
    /**
        * Establish heap order property from an arbitrary
     * arrangement of items. Runs in linear time.
     */
    private void buildHeap( )
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
        * Returns size.
     * @return current size.
     */
    public int size( )
    {
        return currentSize;
    }

    /**
        * Make the priority queue logically empty.
     */
    public void makeEmpty( )
    {
        currentSize = 0;
    }

    
    /**
        * Internal method to percolate down in the heap.
     * @param hole the index at which the percolate begins.
     */
    private void percolateDown( int hole )
    {
        int child;
        Comparable tmp = array[ hole ];
        for( ; hole * 2 <= currentSize; hole = child ){
            child = hole * 2;
            if( (child < currentSize) && (array[child + 1].compareTo(array[child]) < 0) )
                child++;
			
            if( array[ child ].compareTo( tmp ) < 0 )
				array[ hole ] = array[ child ];
            else
                break;
		}
        array[hole] = tmp;
    }
    /**
        * Internal method to extend array.
     */
    private void doubleArray( )
    {
        Comparable [ ] newArray;
        newArray = new Comparable[ array.length * 2 ];
        for( int i = 0; i < array.length; i++ )
            newArray[ i ] = array[ i ];
        array = newArray;
    }
	
	/**
		* Private method to help with arbritary delete. percolate all the way up essentially
	 */
	private void promote(int prom)
	{
		while( (prom > 1) ){
			Comparable temp = array[prom];
			array[prom] = array[prom/2];
			array[prom/2] = temp;
			prom = prom/2;
		}
		
	}

}