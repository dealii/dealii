/*----------------------------   sparsity_pattern.h     ---------------------------*/
/*      $Id$                 */
#ifndef __sparsity_pattern_H
#define __sparsity_pattern_H
/*----------------------------   sparsity_pattern.h     ---------------------------*/



#include <base/exceptions.h>
#include <base/subscriptor.h>
#include <lac/forward_declarations.h>

#include <vector>




/**
 * Structure representing the sparsity pattern of a sparse matrix.
 * 
 * The following picture will illustrate the relation between the
 * #SparsityPattern# an the #SparseMatrix#.
 *
 * \begin{verbatim}
 *  SparsityPattern:                               \
 *                                                 |
 *              _________________________          |
 *  rowstart   |0 | 4| 8|11|13|....                |
 *             |__|__|__|__|__|__________          | 
 *              |   \  \                           |
 *              |    \  \__                        | 
 *              |     \    \                       |
 *              |      \    \__                    |
 *             \ /      \      \                   |
 *              |        \      \__                |     
 *              |         \        \               |
 *              |          \        \__            |   
 *              |           \	       \	   |            
 *              0___________4___________8____       \ Position         
 *  colnums    |3 | 2| 9|17| 1| 4| 6| 8| 4|..       /         
 *             /__|/_|__|__|__|__|__|__|__|__      |         
 *            /   /                                |        
 *           / \______  _____/ \_____  _____/      |         
 *          /         \/             \/            |                 
 *         /   /  row = 0        row = 1           |    
 *        /   /                                    |
 *       /   /                                     |
 *      /   /                                      | 
 *     /   /___colnums[1]                          |  
 *    /                                            |
 *   /_________colnums[0]                          |
 *                                                 |                    
 *                                                /                    
 * \end{verbatim}
 *
 * \begin{verbatim}
 * For row = 0
 *   
 * it exists: (0| 3) = colnums[0]
 *            (0| 2) = colnums[1]
 *            (0| 9) = colnums[2]
 *            (0|17) = colnums[3]
 *
 * For row = 1
 *   
 * it exists: (1| 1) = colnums[4]
 *            (1| 4) = colnums[5]
 * ....
 *
 * \end{verbatim}
 *
 * \begin{verbatim}
 * SparseMatrix:                                  \
 *                                                 |
 *              _____________________________      |
 *  val        |  |  |  |  |  |  |  |  | 3|..       \ Value
 *             |__|__|__|__|__|__|__|__|__|__       /
 *                                                 |
 *                                                 |
 *                                                /
 * \end{verbatim}
 *
 * If you want to get the #3# you need to get its position in the
 * table above and its value by returning the value of the element on
 * which the pointer shows, using #*val#.  For example #val[8]=3#. Its
 * position is #colnums[8]# so #row=2#. In other words, if you want to get
 * the element #a_{24}# you know that #row=2#. To get the element, a
 * search of #4# form #colnums[rowstart[2]]# to #colnums[rowstart[3]]# is
 * needed. Then #a_{24}=val[number of the found element] = 3#.
 *
 *
 * @author Wolfgang Bangerth and others
 */
class SparsityPattern : public Subscriptor
{
  public:
				     /**
				      * Define a value which is used
				      * to indicate that a certain
				      * value in the #colnums# array
				      * is unused, i.e. does not
				      * represent a certain column
				      * number index.
				      *
				      * Indices with this invalid
				      * value are used to insert new
				      * entries to the sparsity
				      * pattern using the #add# member
				      * function, and are removed when
				      * calling #compress#.
				      *
				      * You should not assume that the
				      * variable declared here has a
				      * certain value. The
				      * initialization is given here
				      * only to enable the compiler to
				      * perform some optimizations,
				      * but the actual value of the
				      * variable may change over time.
				      */
    static const unsigned int invalid_entry = static_cast<unsigned int>(-1);
    
				     /**
				      * Initialize the matrix empty,
				      * that is with no memory
				      * allocated. This is useful if
				      * you want such objects as
				      * member variables in other
				      * classes. You can make the
				      * structure usable by calling
				      * the #reinit# function.
				      */
    SparsityPattern ();
    
				     /**
				      * Copy constructor. This constructor is
				      * only allowed to be called if the matrix
				      * structure to be copied is empty. This is
				      * so in order to prevent involuntary
				      * copies of objects for temporaries, which
				      * can use large amounts of computing time.
				      * However, copy constructors are needed
				      * if yo want to use the STL data types
				      * on classes like this, e.g. to write
				      * such statements like
				      * #v.push_back (SparsityPattern());#,
				      * with #v# a vector of #SparsityPattern#
				      * objects.
				      *
				      * Usually, it is sufficient to use the
				      * explicit keyword to disallow unwanted
				      * temporaries, but for the STL vectors,
				      * this does not work. Since copying a
				      * structure like this is not useful
				      * anyway because multiple matrices can
				      * use the same sparsity structure, copies
				      * are only allowed for empty objects, as
				      * described above.
				      */
    SparsityPattern (const SparsityPattern &);

				     /**
				      * Initialize a rectangular matrix with
				      * #m# rows and #n# columns.
				      * The matrix may contain at most #max_per_row#
				      * nonzero entries per row.
				      */
    SparsityPattern (const unsigned int m,
		     const unsigned int n,
		     const unsigned int max_per_row);

				     /**
				      * Initialize a rectangular
				      * matrix with #m# rows and #n#
				      * columns.  The maximal number
				      * of nonzero entries for each
				      * row is given by the
				      * #row_lengths# array.
				      */
    SparsityPattern (const unsigned int          m,
		     const unsigned int          n,
		     const vector<unsigned int> &row_lengths);
    
				     /**
				      * Initialize a square matrix of dimension
				      * #n# with at most #max_per_row#
				      * nonzero entries per row.
				      */
    SparsityPattern (const unsigned int n,
		     const unsigned int max_per_row);

				     /**
				      * Initialize a square
				      * matrix with #m# rows and #m#
				      * columns.  The maximal number
				      * of nonzero entries for each
				      * row is given by the
				      * #row_lengths# array.
				      */
    SparsityPattern (const unsigned int          m,
		     const vector<unsigned int> &row_lengths);

				     /**
				      * Copy operator. For this the same holds
				      * as for the copy constructor: it is
				      * declared, defined and fine to be called,
				      * but the latter only for empty objects.
				      */
    SparsityPattern & operator = (const SparsityPattern &);

				     /**
				      * Make a copy with extra off-diagonals.
				      *
				      * This constructs objects intended for
				      * the application of the ILU(n)-method
				      * or other incomplete decompositions.
				      * Therefore, additional to the original
				      * entry structure, space for
				      * #extra_off_diagonals#
				      * side-diagonals is provided on both
				      * sides of the main diagonal.
				      *
				      * #max_per_row# is the maximum number of
				      * nonzero elements per row which this
				      * structure is to hold. It is assumed
				      * that this number is sufficiently large
				      * to accomodate both the elements in
				      * #original# as well as the new
				      * off-diagonal elements created by this
				      * constructor. You will usually want to
				      * give the same number as you gave for
				      * #original# plus the number of side
				      * diagonals times two. You may however
				      * give a larger value if you wish to add
				      * further nonzero entries for the
				      * decomposition based on other criteria
				      * than their being on side-diagonals.
				      *
				      * This function requires that #original#
				      * refer to a square matrix structure.
				      * It shall be compressed. The matrix 
				      * structure is not compressed
				      * after this function finishes.
				      */
    SparsityPattern (const SparsityPattern& original,
		     const unsigned int        max_per_row,
		     const unsigned int        extra_off_diagonals);
    
				     /**
				      * Destructor.
				      */
    ~SparsityPattern ();
    
				     /**
				      * Reallocate memory and set up data
				      * structures for a new matrix with
				      * #m# rows and #n# columns,
				      * with at most #max_per_row#
				      * nonzero entries per row.
				      *
				      * This function simply maps its
				      * operations to the other
				      * #reinit# function.
				      */
    void reinit (const unsigned int m,
		 const unsigned int n,
		 const unsigned int max_per_row);

				     /**
				      * Reallocate memory for a matrix
				      * of size #m \times n#. The
				      * number of entries for each row
				      * is taken from the array
				      * #row_lengths# which has to
				      * give this number of each row
				      * #i=1...m#.
				      *
				      * If #m*n==0# all memory is freed,
				      * resulting in a total reinitialization
				      * of the object. If it is nonzero, new
				      * memory is only allocated if the new
				      * size extends the old one. This is done
				      * to save time and to avoid fragmentation
				      * of the heap.
				      */
    void reinit (const unsigned int          m,
		 const unsigned int          n,
		 const vector<unsigned int> &row_lengths);
    
				     /**
				      * This function compresses the sparsity
				      * structure that this object represents.
				      * It does so by eliminating unused
				      * entries and sorting the remaining
				      * ones to allow faster access by usage
				      * of binary search algorithms. A special
				      * sorting scheme is used for the diagonal
				      * entry of square matrices, which is
				      * always the first entry of each row.
				      *
				      * The memory which is no more
				      * needed is released.
				      *
				      * #SparseMatrix# objects require the
				      * #SparsityPattern# objects they are
				      * initialized with to be compressed, to
				      * reduce memory requirements.
				      */
    void compress ();

				     /**
				      * Return whether the object is empty. It
				      * is empty if no memory is allocated,
				      * which is the same as that both
				      * dimensions are zero.
				      */
    bool empty () const;

				     /**
				      * Return the maximum number of entries per
				      * row. Before compression, this equals the
				      * number given to the constructor, while
				      * after compression, it equals the maximum
				      * number of entries actually allocated by
				      * the user.
				      */
    unsigned int max_entries_per_row () const;

				     /**
				      * Return the index of the matrix
				      * element with row number #i# and
				      * column number #j#. If the matrix
				      * element is not a nonzero one,
				      * return #SparsityPattern::invalid_entry#.
				      *
				      * This function is usually called
				      * by the #operator()# of the
				      * #SparseMatrix#. It shall only be
				      * called for compressed sparsity
				      * patterns, since in this case
				      * searching whether the entry
				      * exists can be done quite fast
				      * with a binary sort algorithm
				      * because the column numbers are
				      * sorted.
				      */
    unsigned int operator() (const unsigned int i, const unsigned int j) const;

				     /**
				      * Add a nonzero entry to the matrix.
				      * This function may only be called
				      * for non-compressed sparsity patterns.
				      *
				      * If the entry already exists, nothing
				      * bad happens.
				      */
    void add (const unsigned int i, const unsigned int j);
    
				     /**
				      * Print the sparsity of the matrix
				      * in a format that #gnuplot# understands
				      * and which can be used to plot the
				      * sparsity pattern in a graphical
				      * way. The format consists of pairs
				      * #i j# of nonzero elements, each
				      * representing one entry of this
				      * matrix, one per line of the output
				      * file. Indices are counted from
				      * zero on, as usual. Since sparsity
				      * patterns are printed in the same
				      * way as matrices are displayed, we
				      * print the negative of the column
				      * index, which means that the
				      * #(0,0)# element is in the top left
				      * rather than in the bottom left
				      * corner.
				      *
				      * Print the sparsity pattern in
				      * gnuplot by setting the data style
				      * to dots or points and use the
				      * #plot# command.
				      */
    void print_gnuplot (ostream &out) const;

				     /**
				      * Return number of rows of this
				      * matrix, which equals the dimension
				      * of the image space.
				      */
    unsigned int n_rows () const;

				     /**
				      * Return number of columns of this
				      * matrix, which equals the dimension
				      * of the range space.
				      */
    unsigned int n_cols () const;

				     /**
				      * Number of entries in a specific row.
				      */
    unsigned int row_length (const unsigned int row) const;

				     /**
				      * Access to column number field.
				      * Return the column number of
				      * the #index#th entry in #row#.
				      */
    unsigned int column_number (const unsigned int row,
				const unsigned int index) const;

				     /**
				      * Compute the bandwidth of the matrix
				      * represented by this structure. The
				      * bandwidth is the maximum of
				      * $|i-j|$ for which the index pair
				      * $(i,j)$ represents a nonzero entry
				      * of the matrix.
				      */
    unsigned int bandwidth () const;

				     /**
				      * Return the number of nonzero elements of
				      * this matrix. Actually, it returns the
				      * number of entries in the sparsity
				      * pattern; if any of the entries should
				      * happen to be zero, it is counted
				      * anyway.
				      *
				      * This function may only be called if the
				      * matrix struct is compressed. It does not
				      * make too much sense otherwise anyway.
				      */
    unsigned int n_nonzero_elements () const;

				     /**
				      * Return whether the structure is 
				      * compressed or not.
				      */
    bool is_compressed () const;
    
				     /**
				      * This is kind of an expert mode. Get 
				      * access to the rowstart array, but
				      * read-only.
				      *
				      * Use of this function is highly
				      * deprecated. Use #row_length#
				      * and #column_number# instead.
				      *
				      * Though the return value is declared
				      * #const#, you should be aware that it
				      * may change if you call any nonconstant
				      * function of objects which operate on
				      * it.
				      *
				      * You should use this interface very
				      * carefully and only if you are absolutely
				      * sure to know what you do. You should
				      * also note that the structure of these
				      * arrays may change over time.
				      * If you change the layout yourself, you
				      * should also rename this function to
				      * avoid programs relying on outdated
				      * information!  */
    const unsigned int * get_rowstart_indices () const;

				     /**
				      * This is kind of an expert mode: get
				      * access to the colnums array, but
				      * readonly.
				      *
				      * Use of this function is highly
				      * deprecated. Use #row_length#
				      * and #column_number# instead.
				      *
				      * Though the return value is declared
				      * #const#, you should be aware that it
				      * may change if you call any nonconstant
				      * function of objects which operate on
				      * it.
				      *
				      * You should use this interface very
				      * carefully and only if you are absolutely
				      * sure to know what you do. You should
				      * also note that the structure of these
				      * arrays may change over time.
				      * If you change the layout yourself, you
				      * should also rename this function to
				      * avoid programs relying on outdated
				      * information!
				      */
    const unsigned int * get_column_numbers () const;
    
    
				     /**
				      * Exception
				      */
    DeclException1 (ExcInvalidNumber,
		    int,
		    << "The provided number is invalid here: " << arg1);
    				     /**
				      * Exception
				      */
    DeclException2 (ExcInvalidIndex,
		    int, int,
		    << "The given index " << arg1
		    << " should be less than " << arg2 << ".");
				     /**
				      * Exception
				      */
    DeclException2 (ExcNotEnoughSpace,
		    int, int,
		    << "Upon entering a new entry to row " << arg1
		    << ": there was no free entry any more. " << endl
		    << "(Maximum number of entries for this row: "
		    << arg2 << "; maybe the matrix is already compressed?)");
				     /**
				      * Exception
				      */
    DeclException0 (ExcNotCompressed);
				     /**
				      * Exception
				      */
    DeclException0 (ExcMatrixIsCompressed);
				     /**
				      * Exception
				      */
    DeclException0 (ExcEmptyObject);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInternalError);
				     /**
				      * Exception
				      */
    DeclException0 (ExcIO);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidConstructorCall);
				     /**
				      * Exception
				      */
    DeclException0 (ExcNotSquare);
    
  private:
				     /**
				      * Maximum number of rows that can
				      * be stored in the #row_start# array.
				      * Since reallocation of that array
				      * only happens if the present one is
				      * too small, but never when the size
				      * of this matrix structure shrinks,
				      * #max_dim# might be larger than
				      * #rows# and in this case #row_start#
				      * has more elements than are used.
				      */
    unsigned int max_dim;

				     /**
				      * Number of rows that this sparsity
				      * structure shall represent.
				      */
    unsigned int rows;

				     /**
				      * Number of columns that this sparsity
				      * structure shall represent.
				      */
    unsigned int cols;

				     /**
				      * Size of the actually allocated array
				      * #colnums#. Here, the same applies as
				      * for the #rowstart# array, i.e. it
				      * may be larger than the actually used
				      * part of the array.
				      */
    unsigned int max_vec_len;

				     /**
				      * Maximum number of elements per
				      * row. This is set to the value
				      * given to the #reinit# function
				      * (or to the constructor), or to
				      * the maximum row length
				      * computed from the vectors in
				      * case the more flexible
				      * constructors or reinit
				      * versions are called. Its value
				      * is more or less meaningsless
				      * after #compress()# has been
				      * called.
				      */
    unsigned int max_row_length;

				     /**
				      * Array which hold for each row which
				      * is the first element in #colnums#
				      * belonging to that row. Note that
				      * the size of the array is one larger
				      * than the number of rows, because
				      * the last element is used for
				      * #row=rows#, i.e. the row past the
				      * last used one. The value of
				      * #rowstart[rows]# equals the index
				      * of the element past the end in
				      * #colnums#; this way, we are able to
				      * write loops like
				      * #for (i=rowstart[k]; i<rowstart[k+1]; ++i)#
				      * also for the last row.
				      *
				      * Note that the actual size of the
				      * allocated memory may be larger than
				      * the region that is used. The actual
				      * number of elements that was allocated
				      * is stored in #max_dim#.
				      */
    unsigned int *rowstart;

				     /**
				      * Array of column numbers. In this array,
				      * we store for each non-zero element its
				      * column number. The column numbers for
				      * the elements in row #r# are stored
				      * within the index range
				      * #rowstart[r]...rowstart[r+1]#. Therefore
				      * to find out whether a given element
				      * #(r,c)# exists, we have to check
				      * whether the column number #c# exists in
				      * the abovementioned range within this
				      * array. If it exists, say at position
				      * #p# within this array, the value of
				      * the respective element in the sparse
				      * matrix will also be at position #p#
				      * of the values array of that class.
				      *
				      * At the beginning, all elements of
				      * this array are set to #-1# indicating
				      * invalid (unused) column numbers
				      * (however, note that if this object
				      * refers to a square matrix, the diagonal
				      * elements are preset, see below). Now, if
				      * nonzero elements are added, one column
				      * number in the row's respective range
				      * after the other is set to the column
				      * number of the added element. When
				      * compress is called, unused elements
				      * (indicated by column numbers #-1#)
				      * are eliminated by copying the column
				      * number of subsequent rows and the
				      * column numbers within each row (with
				      * the exception of the diagonal element)
				      * are sorted, such that finding whether
				      * an element exists and determining its
				      * position can be done by a binary search.
				      *
				      * If this object represents a square
				      * matrix, the first element in each
				      * row always denotes the diagonal
				      * element, i.e.
				      * #colnums[rowstart[r]]==r#.
				      */
    unsigned int *colnums;

				     /**
				      * Store whether the #compress# function
				      * was called for this object.
				      */
    bool compressed;

    
    template <typename number> friend class SparseMatrix;
};




/*---------------------- Inline functions -----------------------------------*/



inline
unsigned int
SparsityPattern::n_rows () const
{
  return rows;
};



inline
unsigned int
SparsityPattern::n_cols () const
{
  return cols;
};



inline
bool
SparsityPattern::is_compressed () const
{
  return compressed;
};



inline
const unsigned int *
SparsityPattern::get_rowstart_indices () const
{
  return rowstart;
};



inline
const unsigned int *
SparsityPattern::get_column_numbers () const
{
  return colnums;
};



inline
unsigned int
SparsityPattern::row_length (const unsigned int row) const
{
  Assert(row<rows, ExcIndexRange(row,0,rows));
  return rowstart[row+1]-rowstart[row];
}



inline
unsigned int
SparsityPattern::column_number (const unsigned int row,
				const unsigned int index) const
{
  Assert(row<rows, ExcIndexRange(row,0,rows));
  Assert(index<row_length(row), ExcIndexRange(index,0,row_length(row)));

  return colnums[rowstart[row]+index];
}



inline
unsigned int
SparsityPattern::n_nonzero_elements () const
{
  Assert ((rowstart!=0) && (colnums!=0), ExcEmptyObject());  
  Assert (compressed, ExcNotCompressed());
  return rowstart[rows]-rowstart[0];
};



/*----------------------------   sparsity_pattern.h     ---------------------------*/
/* end of #ifndef __sparsity_pattern_H */
#endif
/*----------------------------   sparsity_pattern.h     ---------------------------*/
