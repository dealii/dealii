//----------------------------  sparsity_pattern.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparsity_pattern.h  ---------------------------
#ifndef __deal2__sparsity_pattern_h
#define __deal2__sparsity_pattern_h


#include <base/exceptions.h>
#include <base/subscriptor.h>

template <typename number> class SparseMatrix;
class CompressedSparsityPattern;

#include <vector>
#include <iterator>


/**
 * Structure representing the sparsity pattern of a sparse matrix.
 * 
 * The following picture will illustrate the relation between the
 * @p{SparsityPattern} an the @p{SparseMatrix}.
 *
 * @begin{verbatim}
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
 * @end{verbatim}
 *
 * @begin{verbatim}
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
 * @end{verbatim}
 *
 * @begin{verbatim}
 * SparseMatrix:                                  \
 *                                                 |
 *              _____________________________      |
 *  val        |  |  |  |  |  |  |  |  | 3|..       \ Value
 *             |__|__|__|__|__|__|__|__|__|__       /
 *                                                 |
 *                                                 |
 *                                                /
 * @end{verbatim}
 *
 * If you want to get the @p{3} you need to get its position in the
 * table above and its value by returning the value of the element on
 * which the pointer shows, using @p{*val}.  For example @p{val[8]=3}. Its
 * position is @p{colnums[8]} so @p{row=2}. In other words, if you want to get
 * the element @p{a_{24}} you know that @p{row=2}. To get the element, a
 * search of @p{4} form @p{colnums[rowstart[2]]} to @p{colnums[rowstart[3]]} is
 * needed. Then @p{a_{24}=val[number of the found element] = 3}.
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
				      * value in the @p{colnums} array
				      * is unused, i.e. does not
				      * represent a certain column
				      * number index.
				      *
				      * Indices with this invalid
				      * value are used to insert new
				      * entries to the sparsity
				      * pattern using the @p{add} member
				      * function, and are removed when
				      * calling @p{compress}.
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
				      * the @p{reinit} function.
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
				      * @p{v.push_back (SparsityPattern());},
				      * with @p{v} a vector of @p{SparsityPattern}
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
				      * @p{m} rows and @p{n} columns.
				      * The matrix may contain at most @p{max_per_row}
				      * nonzero entries per row.
				      */
    SparsityPattern (const unsigned int m,
		     const unsigned int n,
		     const unsigned int max_per_row);

				     /**
				      * Initialize a rectangular
				      * matrix with @p{m} rows and @p{n}
				      * columns.  The maximal number
				      * of nonzero entries for each
				      * row is given by the
				      * @p{row_lengths} array.
				      */
    SparsityPattern (const unsigned int               m,
		     const unsigned int               n,
		     const std::vector<unsigned int> &row_lengths);
    
				     /**
				      * Initialize a square matrix of dimension
				      * @p{n} with at most @p{max_per_row}
				      * nonzero entries per row.
				      */
    SparsityPattern (const unsigned int n,
		     const unsigned int max_per_row);

				     /**
				      * Initialize a square
				      * matrix with @p{m} rows and @p{m}
				      * columns.  The maximal number
				      * of nonzero entries for each
				      * row is given by the
				      * @p{row_lengths} array.
				      */
    SparsityPattern (const unsigned int               m,
		     const std::vector<unsigned int> &row_lengths);

				     /**
				      * Make a copy with extra off-diagonals.
				      *
				      * This constructs objects intended for
				      * the application of the ILU(n)-method
				      * or other incomplete decompositions.
				      * Therefore, additional to the original
				      * entry structure, space for
				      * @p{extra_off_diagonals}
				      * side-diagonals is provided on both
				      * sides of the main diagonal.
				      *
				      * @p{max_per_row} is the maximum number of
				      * nonzero elements per row which this
				      * structure is to hold. It is assumed
				      * that this number is sufficiently large
				      * to accomodate both the elements in
				      * @p{original} as well as the new
				      * off-diagonal elements created by this
				      * constructor. You will usually want to
				      * give the same number as you gave for
				      * @p{original} plus the number of side
				      * diagonals times two. You may however
				      * give a larger value if you wish to add
				      * further nonzero entries for the
				      * decomposition based on other criteria
				      * than their being on side-diagonals.
				      *
				      * This function requires that @p{original}
				      * refer to a square matrix structure.
				      * It shall be compressed. The matrix 
				      * structure is not compressed
				      * after this function finishes.
				      */
    SparsityPattern (const SparsityPattern  &original,
		     const unsigned int      max_per_row,
		     const unsigned int      extra_off_diagonals);
    
				     /**
				      * Destructor.
				      */
    ~SparsityPattern ();

				     /**
				      * Copy operator. For this the same holds
				      * as for the copy constructor: it is
				      * declared, defined and fine to be called,
				      * but the latter only for empty objects.
				      */
    SparsityPattern & operator = (const SparsityPattern &);
    
				     /**
				      * Reallocate memory and set up data
				      * structures for a new matrix with
				      * @p{m} rows and @p{n} columns,
				      * with at most @p{max_per_row}
				      * nonzero entries per row.
				      *
				      * This function simply maps its
				      * operations to the other
				      * @p{reinit} function.
				      */
    void reinit (const unsigned int m,
		 const unsigned int n,
		 const unsigned int max_per_row);

				     /**
				      * Reallocate memory for a matrix
				      * of size @p{m \times n}. The
				      * number of entries for each row
				      * is taken from the array
				      * @p{row_lengths} which has to
				      * give this number of each row
				      * @p{i=1...m}.
				      *
				      * If @p{m*n==0} all memory is freed,
				      * resulting in a total reinitialization
				      * of the object. If it is nonzero, new
				      * memory is only allocated if the new
				      * size extends the old one. This is done
				      * to save time and to avoid fragmentation
				      * of the heap.
				      */
    void reinit (const unsigned int               m,
		 const unsigned int               n,
		 const std::vector<unsigned int> &row_lengths);
    
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
				      * @p{SparseMatrix} objects require the
				      * @p{SparsityPattern} objects they are
				      * initialized with to be compressed, to
				      * reduce memory requirements.
				      */
    void compress ();

				     /**
				      * This function can be used as a
				      * replacement for @ref{reinit},
				      * subsequent calls to @ref{add}
				      * and a final call to
				      * @ref{close} if you know
				      * exactly in advance the entries
				      * that will form the matrix
				      * sparsity pattern.
				      *
				      * The first two parameters
				      * determine the size of the
				      * matrix. For the two last ones,
				      * note that a sparse matrix can
				      * be described by a sequence of
				      * rows, each of which is
				      * represented by a sequence of
				      * pairs of column indices and
				      * values. In the present
				      * context, the @ref{begin} and
				      * @ref{end} parameters designate
				      * iterators (of forward iterator
				      * type) into a container, one
				      * representing one row. The
				      * distance between @ref{begin}
				      * and @ref{end} should therefore
				      * be equal to
				      * @ref{n_rows}. These iterators
				      * may be iterators of
				      * @p{std::vector},
				      * @p{std::list}, pointers into a
				      * C-style array, or any other
				      * iterator satisfying the
				      * requirements of a forward
				      * iterator. The objects pointed
				      * to by these iterators
				      * (i.e. what we get after
				      * applying @p{operator*} or
				      * @p{operator->} to one of these
				      * iterators) must be a container
				      * itself that provides functions
				      * @p{begin} and @p{end}
				      * designating a range of
				      * iterators that describe the
				      * contents of one
				      * line. Dereferencing these
				      * inner iterators must either
				      * yield a pair of an unsigned
				      * integer as column index and a
				      * value of arbitrary type (such
				      * a type would be used if we
				      * wanted to describe a sparse
				      * matrix with one such object),
				      * or simply an unsigned integer
				      * (of we only wanted to describe
				      * a sparsity pattern). The
				      * function is able to determine
				      * itself whether an unsigned
				      * integer or a pair is what we
				      * get after dereferencing the
				      * inner iterators, through some
				      * template magic.
				      *
				      * While the order of the outer
				      * iterators denotes the
				      * different rows of the matrix,
				      * the order of the inner
				      * iterator denoting the columns
				      * does not matter, as they are
				      * sorted internal to this
				      * function anyway.
				      *
				      * Since that all sounds very
				      * complicated, consider the
				      * following example code, which
				      * may be used to fill a sparsity
				      * pattern:
				      * @begin{verbatim}
				      * std::vector<std::vector<unsigned int> > column_indices (n_rows);
				      * for (unsigned int row=0; row<n_rows; ++row)
				      *         // generate necessary columns in this row
				      *   fill_row (column_indices[row]);
				      *
				      * sparsity.copy_from (n_rows, n_cols,
				      *                     column_indices.begin(),
				      *                     column_indices.end());
				      * @end{verbatim}
				      *
				      * Note that this example works
				      * since the iterators
				      * dereferenced yield containers
				      * with functions @p{begin} and
				      * @p{end} (namely
				      * @p{std::vector}s), and the
				      * inner iterators dereferenced
				      * yield unsigned integers as
				      * column indices. Note that we
				      * could have replaced each of
				      * the two @p{std::vector}
				      * occurrences by @p{std::list},
				      * and the inner one by
				      * @p{std::set} as well.
				      *
				      * Another example would be as
				      * follows, where we initialize a
				      * whole matrix, not only a
				      * sparsity pattern:
				      * @begin{verbatim}
				      * std::vector<std::map<unsigned int,double> > entries (n_rows);
				      * for (unsigned int row=0; row<n_rows; ++row)
				      *         // generate necessary pairs of columns
				      *         // and corresponding values in this row
				      *   fill_row (entries[row]);
				      *
				      * sparsity.copy_from (n_rows, n_cols,
				      *                     column_indices.begin(),
				      *                     column_indices.end());
				      * matrix.reinit (sparsity);
				      * matrix.copy_from (column_indices.begin(),
				      *                   column_indices.end());
				      * @end{verbatim}
				      *
				      * This example works because
				      * dereferencing iterators of the
				      * inner type yields a pair of
				      * unsigned integers and a value,
				      * the first of which we take as
				      * column index. As previously,
				      * the outer @p{std::vector}
				      * could be replaced by
				      * @p{std::list}, and the inner
				      * @p{std::map<unsigned int,double>}
				      * could be replaced by
				      * @p{std::vector<std::pair<unsigned int,double> >},
				      * or a list or set of such
				      * pairs, as they all return
				      * iterators that point to such
				      * pairs.
				      */
    template <typename ForwardIterator>
    void copy_from (const unsigned int    n_rows,
		    const unsigned int    n_cols,
		    const ForwardIterator begin,
		    const ForwardIterator end);

				     /**
				      * Copy data from an object of
				      * type
				      * @ref{CompressedSparsityPattern}. Previous
				      * content of this object is
				      * lost.
				      */
    void copy_from (const CompressedSparsityPattern &csp);
    
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
				      * element with row number @p{i} and
				      * column number @p{j}. If the matrix
				      * element is not a nonzero one,
				      * return @p{SparsityPattern::invalid_entry}.
				      *
				      * This function is usually called
				      * by the @p{operator()} of the
				      * @p{SparseMatrix}. It shall only be
				      * called for compressed sparsity
				      * patterns, since in this case
				      * searching whether the entry
				      * exists can be done quite fast
				      * with a binary sort algorithm
				      * because the column numbers are
				      * sorted.
				      */
    unsigned int operator() (const unsigned int i, 
			     const unsigned int j) const;

				     /**
				      * Add a nonzero entry to the matrix.
				      * This function may only be called
				      * for non-compressed sparsity patterns.
				      *
				      * If the entry already exists, nothing
				      * bad happens.
				      */
    void add (const unsigned int i, 
	      const unsigned int j);
    
				     /**
				      * Make the sparsity pattern
				      * symmetric by adding the
				      * sparsity pattern of the
				      * transpose object.
				      *
				      * This function throws an
				      * exception if the sparsity
				      * pattern does not represent a
				      * square matrix.
				      */
    void symmetrize ();
    
				     /**
				      * Print the sparsity of the matrix
				      * in a format that @p{gnuplot} understands
				      * and which can be used to plot the
				      * sparsity pattern in a graphical
				      * way. The format consists of pairs
				      * @p{i j} of nonzero elements, each
				      * representing one entry of this
				      * matrix, one per line of the output
				      * file. Indices are counted from
				      * zero on, as usual. Since sparsity
				      * patterns are printed in the same
				      * way as matrices are displayed, we
				      * print the negative of the column
				      * index, which means that the
				      * @p{(0,0)} element is in the top left
				      * rather than in the bottom left
				      * corner.
				      *
				      * Print the sparsity pattern in
				      * gnuplot by setting the data style
				      * to dots or points and use the
				      * @p{plot} command.
				      */
    void print_gnuplot (std::ostream &out) const;

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
				      * the @p{index}th entry in @p{row}.
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
				      * deprecated. Use @p{row_length}
				      * and @p{column_number} instead.
				      *
				      * Though the return value is declared
				      * @p{const}, you should be aware that it
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
				      * deprecated. Use @p{row_length}
				      * and @p{column_number} instead.
				      *
				      * Though the return value is declared
				      * @p{const}, you should be aware that it
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
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      */
    unsigned int memory_consumption () const;

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
		    << ": there was no free entry any more. " << std::endl
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
    DeclException0 (ExcInvalidConstructorCall);
				     /**
				      * Exception
				      */
    DeclException0 (ExcNotSquare);
				     /**
				      * Exception
				      */
    DeclException2 (ExcIteratorRange,
		    int, int,
		    << "The iterators denote a range of " << arg1
		    << " elements, but the given number of rows was " << arg2);
    
  private:
				     /**
				      * Maximum number of rows that can
				      * be stored in the @p{row_start} array.
				      * Since reallocation of that array
				      * only happens if the present one is
				      * too small, but never when the size
				      * of this matrix structure shrinks,
				      * @p{max_dim} might be larger than
				      * @p{rows} and in this case @p{row_start}
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
				      * @p{colnums}. Here, the same applies as
				      * for the @p{rowstart} array, i.e. it
				      * may be larger than the actually used
				      * part of the array.
				      */
    unsigned int max_vec_len;

				     /**
				      * Maximum number of elements per
				      * row. This is set to the value
				      * given to the @p{reinit} function
				      * (or to the constructor), or to
				      * the maximum row length
				      * computed from the vectors in
				      * case the more flexible
				      * constructors or reinit
				      * versions are called. Its value
				      * is more or less meaningsless
				      * after @p{compress()} has been
				      * called.
				      */
    unsigned int max_row_length;

				     /**
				      * Array which hold for each row which
				      * is the first element in @p{colnums}
				      * belonging to that row. Note that
				      * the size of the array is one larger
				      * than the number of rows, because
				      * the last element is used for
				      * @p{row=rows}, i.e. the row past the
				      * last used one. The value of
				      * @p{rowstart[rows]} equals the index
				      * of the element past the end in
				      * @p{colnums}; this way, we are able to
				      * write loops like
				      * @p{for (i=rowstart[k]; i<rowstart[k+1]; ++i)}
				      * also for the last row.
				      *
				      * Note that the actual size of the
				      * allocated memory may be larger than
				      * the region that is used. The actual
				      * number of elements that was allocated
				      * is stored in @p{max_dim}.
				      */
    unsigned int *rowstart;

				     /**
				      * Array of column numbers. In this array,
				      * we store for each non-zero element its
				      * column number. The column numbers for
				      * the elements in row @p{r} are stored
				      * within the index range
				      * @p{rowstart[r]...rowstart[r+1]}. Therefore
				      * to find out whether a given element
				      * @p{(r,c)} exists, we have to check
				      * whether the column number @p{c} exists in
				      * the abovementioned range within this
				      * array. If it exists, say at position
				      * @p{p} within this array, the value of
				      * the respective element in the sparse
				      * matrix will also be at position @p{p}
				      * of the values array of that class.
				      *
				      * At the beginning, all elements of
				      * this array are set to @p{-1} indicating
				      * invalid (unused) column numbers
				      * (however, note that if this object
				      * refers to a square matrix, the diagonal
				      * elements are preset, see below). Now, if
				      * nonzero elements are added, one column
				      * number in the row's respective range
				      * after the other is set to the column
				      * number of the added element. When
				      * compress is called, unused elements
				      * (indicated by column numbers @p{-1})
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
				      * @p{colnums[rowstart[r]]==r}.
				      */
    unsigned int *colnums;

				     /**
				      * Store whether the @p{compress} function
				      * was called for this object.
				      */
    bool compressed;

				     /**
				      * Optimized replacement for
				      * @p{std::lower_bound} for
				      * searching within the range of
				      * column indices. Slashes
				      * execution time by
				      * approximately one half for the
				      * present application, partly
				      * because we have eliminated
				      * templates and the compiler
				      * seems to be able to optimize
				      * better, and partly because the
				      * binary search is replaced by a
				      * linear search for small loop
				      * lengths.
				      */
    static
    const unsigned int * const
    optimized_lower_bound (const unsigned int *first,
			   const unsigned int *last,
			   const unsigned int &val);

				     /**
				      * Helper function to get the
				      * column index from a
				      * dereferenced iterator in the
				      * @ref{copy_from} function, if
				      * the inner iterator type points
				      * to plain unsigned integers.
				      */
    static
    unsigned int
    get_column_index_from_iterator (const unsigned int i);
    
				     /**
				      * Helper function to get the
				      * column index from a
				      * dereferenced iterator in the
				      * @ref{copy_from} function, if
				      * the inner iterator type points
				      * to pairs of unsigned integers
				      * and some other value.
				      */
    template <typename value>
    unsigned int
    get_column_index_from_iterator (const std::pair<unsigned int, value> &i);

				     /**
				      * Likewise, but sometimes needed
				      * for certain types of
				      * containers that make the first
				      * element of the pair constant
				      * (such as @p{std::map}).
				      */
    template <typename value>
    unsigned int
    get_column_index_from_iterator (const std::pair<const unsigned int, value> &i);
    
				     /**
				      * Make all sparse matrices
				      * friends of this class.
				      */
    template <typename number> friend class SparseMatrix;
};


/*---------------------- Inline functions -----------------------------------*/


inline
const unsigned int * const
SparsityPattern::optimized_lower_bound (const unsigned int *first,
					const unsigned int *last,
					const unsigned int &val)
{
				   // this function is mostly copied
				   // over from the STL __lower_bound
				   // function, but with template args
				   // replaced by the actual data
				   // types needed here, and above all
				   // with a rolled out search on the
				   // last four elements
  unsigned int len = last-first;

  if (len==0)
    return first;
  
  while (true)
    {
				       // if length equals 8 or less,
				       // then do a rolled out
				       // search. use a switch without
				       // breaks for that and roll-out
				       // the loop somehow
      if (len < 8)
	{
	  switch (len)
	    {
	      case 7:
		    if (*first >= val)
		      return first;
		    ++first;
	      case 6:
		    if (*first >= val)
		      return first;
		    ++first;
	      case 5:
		    if (*first >= val)
		      return first;
		    ++first;
	      case 4:
		    if (*first >= val)
		      return first;
		    ++first;
	      case 3:
		    if (*first >= val)
		      return first;
		    ++first;
	      case 2:
		    if (*first >= val)
		      return first;
		    ++first;
	      case 1:
		    if (*first >= val)
		      return first;
		    return first+1;
	      default:
						     // indices seem
						     // to not be
						     // sorted
						     // correctly!? or
						     // did len
						     // become==0
						     // somehow? that
						     // shouln't have
						     // happened
		    Assert (false, ExcInternalError());
	    };
	};
      
      
      
      const unsigned int         half   = len >> 1;
      const unsigned int * const middle = first + half;
      
				       // if the value is larger than
				       // that pointed to by the
				       // middle pointer, then the
				       // insertion point must be
				       // right of it
      if (*middle < val)
	{
	  first = middle + 1;
	  len  -= half + 1;
	}
      else
	len = half;
    }
}



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



inline
unsigned int
SparsityPattern::get_column_index_from_iterator (const unsigned int i)
{
  return i;
};



template <typename value>
inline
unsigned int
SparsityPattern::get_column_index_from_iterator (const std::pair<unsigned int, value> &i)
{
  return i.first;
};



template <typename value>
inline
unsigned int
SparsityPattern::get_column_index_from_iterator (const std::pair<const unsigned int, value> &i)
{
  return i.first;
};



template <typename ForwardIterator>
void
SparsityPattern::copy_from (const unsigned int    n_rows,
			    const unsigned int    n_cols,
			    const ForwardIterator begin,
			    const ForwardIterator end)
{
  Assert (static_cast<unsigned int>(std::distance (begin, end)) == n_rows,
	  ExcIteratorRange (std::distance (begin, end), n_rows));
  
				   // first determine row lengths for
				   // each row. if the matrix is
				   // square, then we might have to
				   // add an additional entry for the
				   // diagonal, if that is not yet
				   // present. as we have to call
				   // compress anyway later on, don't
				   // bother to check whether that
				   // diagonal entry is in a certain
				   // row or not
  const bool is_square = (n_rows == n_cols);
  std::vector<unsigned int> row_lengths;
  row_lengths.reserve(n_rows);
  for (ForwardIterator i=begin; i!=end; ++i)
    row_lengths.push_back (std::distance (i->begin(), i->end())
			   +
			   (is_square ? 1 : 0));
  reinit (n_rows, n_cols, row_lengths);

				   // now enter all the elements into
				   // the matrix. note that if the
				   // matrix is square, then we
				   // already have the diagonal
				   // element preallocated
				   //
				   // for use in the inner loop, we
				   // define a typedef to the type of
				   // the inner iterators
  unsigned int row = 0;
  typedef typename std::iterator_traits<ForwardIterator>::value_type::const_iterator inner_iterator;
  for (ForwardIterator i=begin; i!=end; ++i, ++row)
    {
      unsigned int *cols = &colnums[rowstart[row]] + (is_square ? 1 : 0);
      const inner_iterator end_of_row = i->end();
      for (inner_iterator j=i->begin(); j!=end_of_row; ++j)
	{
	  const unsigned int col = get_column_index_from_iterator(*j);
	  Assert (col < n_cols, ExcInvalidIndex(col,n_cols));
	  
	  if ((col!=row) || !is_square)
	    *cols++ = col;
	};
    };

				   // finally compress
				   // everything. this also sorts the
				   // entries within each row
  compress ();
};


#endif
