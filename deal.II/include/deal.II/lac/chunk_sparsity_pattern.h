//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__chunk_sparsity_pattern_h
#define __deal2__chunk_sparsity_pattern_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/vector_slice.h>

#include <deal.II/lac/sparsity_pattern.h>

#include <vector>
#include <iostream>

DEAL_II_NAMESPACE_OPEN


template <typename> class ChunkSparseMatrix;


/*! @addtogroup Sparsity
 *@{
 */


/**
 * Structure representing the sparsity pattern of a sparse matrix.
 *
 * This class is an example of the "static" type of @ref Sparsity.
 *
 * It uses the compressed row storage (CSR) format to store data.
 *
 * @author Wolfgang Bangerth, 2008
 */
class ChunkSparsityPattern : public Subscriptor
{
  public:
    
				     /**
				      * Define a value which is used
				      * to indicate that a certain
				      * value in the colnums array
				      * is unused, i.e. does not
				      * represent a certain column
				      * number index.
				      *
				      * Indices with this invalid
				      * value are used to insert new
				      * entries to the sparsity
				      * pattern using the add() member
				      * function, and are removed when
				      * calling compress().
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
    static const unsigned int invalid_entry = SparsityPattern::invalid_entry;

				     /**
				      * Initialize the matrix empty,
				      * that is with no memory
				      * allocated. This is useful if
				      * you want such objects as
				      * member variables in other
				      * classes. You can make the
				      * structure usable by calling
				      * the reinit() function.
				      */
    ChunkSparsityPattern ();
    
				     /**
				      * Copy constructor. This
				      * constructor is only allowed to
				      * be called if the matrix
				      * structure to be copied is
				      * empty. This is so in order to
				      * prevent involuntary copies of
				      * objects for temporaries, which
				      * can use large amounts of
				      * computing time.  However, copy
				      * constructors are needed if yo
				      * want to use the STL data types
				      * on classes like this, e.g. to
				      * write such statements like
				      * <tt>v.push_back
				      * (ChunkSparsityPattern());</tt>,
				      * with <tt>v</tt> a vector of
				      * ChunkSparsityPattern objects.
				      *
				      * Usually, it is sufficient to
				      * use the explicit keyword to
				      * disallow unwanted temporaries,
				      * but for the STL vectors, this
				      * does not work. Since copying a
				      * structure like this is not
				      * useful anyway because multiple
				      * matrices can use the same
				      * sparsity structure, copies are
				      * only allowed for empty
				      * objects, as described above.
				      */
    ChunkSparsityPattern (const ChunkSparsityPattern &);

				     /**
				      * Initialize a rectangular
				      * matrix.
				      *
				      * @arg m number of rows
				      * @arg n number of columns
				      * @arg max_per_row maximum
				      * number of nonzero entries per row
				      *
				      * @arg optimize_diagonal store
				      * diagonal entries first in row;
				      * see optimize_diagonal(). This
				      * takes effect for quadratic
				      * matrices only.
				      */
    ChunkSparsityPattern (const unsigned int m,
			  const unsigned int n,
			  const unsigned int max_chunks_per_row,
			  const unsigned int chunk_size,
			  const bool optimize_diagonal = true);

				     /**
				      * Initialize a rectangular
				      * matrix.
				      *
				      * @arg m number of rows
				      * @arg n number of columns
				      *
				      * @arg row_lengths possible
				      * number of nonzero entries for
				      * each row.  This vector must
				      * have one entry for each row.
				      *
				      * @arg optimize_diagonal store
				      * diagonal entries first in row;
				      * see optimize_diagonal(). This
				      * takes effect for quadratic
				      * matrices only.
				      */
    ChunkSparsityPattern (const unsigned int               m,
			  const unsigned int               n,
			  const std::vector<unsigned int>& row_lengths,
			  const unsigned int chunk_size,
			  const bool optimize_diagonal = true);
    
				     /**
				      * Initialize a quadratic matrix
				      * of dimension <tt>n</tt> with
				      * at most <tt>max_per_row</tt>
				      * nonzero entries per row.
				      *
				      * This constructor automatically
				      * enables optimized storage of
				      * diagonal elements. To avoid
				      * this, use the constructor
				      * taking row and column numbers
				      * separately.
				      */
    ChunkSparsityPattern (const unsigned int n,
			  const unsigned int max_per_row,
			  const unsigned int chunk_size);

				     /**
				      * Initialize a quadratic matrix.
				      *
				      * @arg m number of rows and columns
				      *
				      * @arg row_lengths possible
				      * number of nonzero entries for
				      * each row.  This vector must
				      * have one entry for each row.
				      *
				      * @arg optimize_diagonal store
				      * diagonal entries first in row;
				      * see optimize_diagonal().
				      */
    ChunkSparsityPattern (const unsigned int               m,
			  const std::vector<unsigned int>& row_lengths,
			  const unsigned int               chunk_size,
			  const bool optimize_diagonal = true);
    
				     /**
				      * Destructor.
				      */
    ~ChunkSparsityPattern ();

				     /**
				      * Copy operator. For this the
				      * same holds as for the copy
				      * constructor: it is declared,
				      * defined and fine to be called,
				      * but the latter only for empty
				      * objects.
				      */
    ChunkSparsityPattern & operator = (const ChunkSparsityPattern &);
    
				     /**
				      * Reallocate memory and set up data
				      * structures for a new matrix with
				      * <tt>m </tt>rows and <tt>n</tt> columns,
				      * with at most <tt>max_per_row</tt>
				      * nonzero entries per row.
				      *
				      * This function simply maps its
				      * operations to the other
				      * <tt>reinit</tt> function.
				      */
    void reinit (const unsigned int m,
		 const unsigned int n,
		 const unsigned int max_per_row,
		 const unsigned int chunk_size,
		 const bool optimize_diagonal = true);

				     /**
				      * Reallocate memory for a matrix
				      * of size <tt>m x n</tt>. The
				      * number of entries for each row
				      * is taken from the array
				      * <tt>row_lengths</tt> which has to
				      * give this number of each row
				      * <tt>i=1...m</tt>.
				      *
				      * If <tt>m*n==0</tt> all memory is freed,
				      * resulting in a total reinitialization
				      * of the object. If it is nonzero, new
				      * memory is only allocated if the new
				      * size extends the old one. This is done
				      * to save time and to avoid fragmentation
				      * of the heap.
				      *
				      * If the number of rows equals
				      * the number of columns and the
				      * last parameter is true,
				      * diagonal elements are stored
				      * first in each row to allow
				      * optimized access in relaxation
				      * methods of SparseMatrix.
				      */
    void reinit (const unsigned int               m,
		 const unsigned int               n,
		 const std::vector<unsigned int> &row_lengths,
		 const unsigned int chunk_size,
		 const bool optimize_diagonal = true);

				     /**
				      * Same as above, but with a
				      * VectorSlice argument instead.
				      */
    void reinit (const unsigned int               m,
		 const unsigned int               n,
		 const VectorSlice<const std::vector<unsigned int> > &row_lengths,
		 const unsigned int chunk_size,
		 const bool optimize_diagonal = true);
    
				     /**
				      * This function compresses the sparsity
				      * structure that this object represents.
				      * It does so by eliminating unused
				      * entries and sorting the remaining ones
				      * to allow faster access by usage of
				      * binary search algorithms. A special
				      * sorting scheme is used for the
				      * diagonal entry of quadratic matrices,
				      * which is always the first entry of
				      * each row.
				      *
				      * The memory which is no more
				      * needed is released.
				      *
				      * SparseMatrix objects require the
				      * ChunkSparsityPattern objects they are
				      * initialized with to be compressed, to
				      * reduce memory requirements.
				      */
    void compress ();
    
    				     /**
				      * This function can be used as a
				      * replacement for reinit(),
				      * subsequent calls to add() and
				      * a final call to close() if you
				      * know exactly in advance the
				      * entries that will form the
				      * matrix sparsity pattern.
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
				      * context, the begin() and
				      * end() parameters designate
				      * iterators (of forward iterator
				      * type) into a container, one
				      * representing one row. The
				      * distance between begin()
				      * and end() should therefore
				      * be equal to
				      * n_rows(). These iterators
				      * may be iterators of
				      * <tt>std::vector</tt>,
				      * <tt>std::list</tt>, pointers into a
				      * C-style array, or any other
				      * iterator satisfying the
				      * requirements of a forward
				      * iterator. The objects pointed
				      * to by these iterators
				      * (i.e. what we get after
				      * applying <tt>operator*</tt> or
				      * <tt>operator-></tt> to one of these
				      * iterators) must be a container
				      * itself that provides functions
				      * <tt>begin</tt> and <tt>end</tt>
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
				      * @code
				      * std::vector<std::vector<unsigned int> > column_indices (n_rows);
				      * for (unsigned int row=0; row<n_rows; ++row)
				      *         // generate necessary columns in this row
				      *   fill_row (column_indices[row]);
				      *
				      * sparsity.copy_from (n_rows, n_cols,
				      *                     column_indices.begin(),
				      *                     column_indices.end());
				      * @endcode
				      *
				      * Note that this example works
				      * since the iterators
				      * dereferenced yield containers
				      * with functions <tt>begin</tt> and
				      * <tt>end</tt> (namely
				      * <tt>std::vector</tt>s), and the
				      * inner iterators dereferenced
				      * yield unsigned integers as
				      * column indices. Note that we
				      * could have replaced each of
				      * the two <tt>std::vector</tt>
				      * occurrences by <tt>std::list</tt>,
				      * and the inner one by
				      * <tt>std::set</tt> as well.
				      *
				      * Another example would be as
				      * follows, where we initialize a
				      * whole matrix, not only a
				      * sparsity pattern:
				      * @code
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
				      * @endcode
				      *
				      * This example works because
				      * dereferencing iterators of the
				      * inner type yields a pair of
				      * unsigned integers and a value,
				      * the first of which we take as
				      * column index. As previously,
				      * the outer <tt>std::vector</tt>
				      * could be replaced by
				      * <tt>std::list</tt>, and the inner
				      * <tt>std::map<unsigned int,double></tt>
				      * could be replaced by
				      * <tt>std::vector<std::pair<unsigned int,double> ></tt>,
				      * or a list or set of such
				      * pairs, as they all return
				      * iterators that point to such
				      * pairs.
				      */
    template <typename ForwardIterator>
    void copy_from (const unsigned int    n_rows,
		    const unsigned int    n_cols,
		    const ForwardIterator begin,
		    const ForwardIterator end,
		    const unsigned int chunk_size,
		    const bool optimize_diagonal = true);

				     /**
				      * Copy data from an object of type
				      * CompressedSparsityPattern,
				      * CompressedSetSparsityPattern or
				      * CompressedSimpleSparsityPattern.
				      * Previous content of this object is
				      * lost, and the sparsity pattern is in
				      * compressed mode afterwards.
				      */
    template <typename SparsityType>
    void copy_from (const SparsityType &csp,
		    const unsigned int  chunk_size,
		    const bool          optimize_diagonal = true);

				     /**
				      * Take a full matrix and use its
				      * nonzero entries to generate a
				      * sparse matrix entry pattern
				      * for this object.
				      *
				      * Previous content of this
				      * object is lost, and the
				      * sparsity pattern is in
				      * compressed mode afterwards.
				      */
    template <typename number>
    void copy_from (const FullMatrix<number> &matrix,
		    const unsigned int chunk_size,
		    const bool optimize_diagonal = true);
    
				     /**
				      * Return whether the object is empty. It
				      * is empty if no memory is allocated,
				      * which is the same as that both
				      * dimensions are zero.
				      */
    bool empty () const;

				     /**
				      * Return the chunk size given as
				      * argument when constructing this
				      * object.
				      */
    unsigned int get_chunk_size () const;
    
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
				      * quadratic matrix.
				      */
    void symmetrize ();
    
				     /**
				      * Return number of rows of this
				      * matrix, which equals the dimension
				      * of the image space.
				      */
    inline unsigned int n_rows () const;

				     /**
				      * Return number of columns of this
				      * matrix, which equals the dimension
				      * of the range space.
				      */
    inline unsigned int n_cols () const;

    				     /**
				      * Check if a value at a certain
				      * position may be non-zero.
				      */
    bool exists (const unsigned int i,
                 const unsigned int j) const;

				     /**
				      * Number of entries in a specific row.
				      */
    unsigned int row_length (const unsigned int row) const;

				     /**
				      * Compute the bandwidth of the matrix
				      * represented by this structure. The
				      * bandwidth is the maximum of $|i-j|$
				      * for which the index pair $(i,j)$
				      * represents a nonzero entry of the
				      * matrix. Consequently, the maximum
				      * bandwidth a $n\times m$ matrix can
				      * have is $\max\{n-1,m-1\}$.
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
				      * Determine whether the matrix
				      * uses special convention for
				      * quadratic matrices.
				      *
				      * A return value <tt>true</tt> means
				      * that diagonal elements are stored
				      * first in each row. A number of
				      * functions in this class and the
				      * library in general, for example
				      * relaxation methods like Jacobi() and
				      * SOR(), require this to make their
				      * operations more efficient, since they
				      * need to quickly access the diagonal
				      * elements and do not have to search for
				      * them if they are the first element of
				      * each row. A side effect of this scheme
				      * is that each row contains at least one
				      * element, even if the row is empty
				      * (i.e. the diagonal element exists, but
				      * has value zero).
				      *
				      * A return value <tt>false</tt> means
				      * that diagonal elements are stored
				      * anywhere in the row, or not at all. In
				      * particular, a row or even the whole
				      * matrix may be empty. This can be used
				      * if you have block matrices where the
				      * off-diagonal blocks are quadratic but
				      * are never used for operations like the
				      * ones mentioned above. In this case,
				      * some memory can be saved by not using
				      * the diagonal storage optimization.
                                      */
    bool optimize_diagonal () const;

				     /**
				      * Return whether this object stores only
				      * those entries that have been added
				      * explicitly, or if the sparsity pattern
				      * contains elements that have been added
				      * through other means (implicitly) while
				      * building it. For the current class,
				      * the result is false because we store
				      * entire chunks, not individual
				      * elements, and adding one entry to the
				      * sparsity pattern requires also adding
				      * all the other elements of a chunk. The
				      * only exception is if
				      * <code>chunk_size==1</code>, the
				      * sparsity pattern is nonsymmetric or
				      * optimize_diag has been set to false.
				      *
				      * This function mainly serves the
				      * purpose of describing the current
				      * class in cases where several kinds of
				      * sparsity patterns can be passed as
				      * template arguments.
				      */
    bool stores_only_added_elements () const;
    
				     /**
				      * Write the data of this object
				      * en bloc to a file. This is
				      * done in a binary mode, so the
				      * output is neither readable by
				      * humans nor (probably) by other
				      * computers using a different
				      * operating system of number
				      * format.
				      *
				      * The purpose of this function
				      * is that you can swap out
				      * matrices and sparsity pattern
				      * if you are short of memory,
				      * want to communicate between
				      * different programs, or allow
				      * objects to be persistent
				      * across different runs of the
				      * program.
				      */
    void block_write (std::ostream &out) const;

				     /**
				      * Read data that has previously
				      * been written by block_write()
				      * from a file. This is done
				      * using the inverse operations
				      * to the above function, so it
				      * is reasonably fast because the
				      * bitstream is not interpreted
				      * except for a few numbers up
				      * front.
				      *
				      * The object is resized on this
				      * operation, and all previous
				      * contents are lost.
				      *
				      * A primitive form of error
				      * checking is performed which
				      * will recognize the bluntest
				      * attempts to interpret some
				      * data as a vector stored
				      * bitwise to a file, but not
				      * more.
				      */
    void block_read (std::istream &in);

				     /**
				      * Print the sparsity of the
				      * matrix. The output consists of
				      * one line per row of the format
				      * <tt>[i,j1,j2,j3,...]</tt>. <i>i</i>
				      * is the row number and
				      * <i>jn</i> are the allocated
				      * columns in this row.
				      */
    void print (std::ostream &out) const;

				     /**
				      * Print the sparsity of the matrix
				      * in a format that <tt>gnuplot</tt> understands
				      * and which can be used to plot the
				      * sparsity pattern in a graphical
				      * way. The format consists of pairs
				      * <tt>i j</tt> of nonzero elements, each
				      * representing one entry of this
				      * matrix, one per line of the output
				      * file. Indices are counted from
				      * zero on, as usual. Since sparsity
				      * patterns are printed in the same
				      * way as matrices are displayed, we
				      * print the negative of the column
				      * index, which means that the
				      * <tt>(0,0)</tt> element is in the top left
				      * rather than in the bottom left
				      * corner.
				      *
				      * Print the sparsity pattern in
				      * gnuplot by setting the data style
				      * to dots or points and use the
				      * <tt>plot</tt> command.
				      */
    void print_gnuplot (std::ostream &out) const;

				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object. See
				      * MemoryConsumption.
				      */
    std::size_t memory_consumption () const;

				     /** @addtogroup Exceptions
				      * @{ */
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
				      * This exception is thrown if
				      * the matrix does not follow the
				      * convention of storing diagonal
				      * elements first in row. Refer
				      * to
				      * SparityPattern::optimize_diagonal()
				      * for more information.
				      */
    DeclException0 (ExcDiagonalNotOptimized);
				     /**
				      * Exception
				      */
    DeclException2 (ExcIteratorRange,
		    int, int,
		    << "The iterators denote a range of " << arg1
		    << " elements, but the given number of rows was " << arg2);
                                     /**
                                      * Exception
                                      */
    DeclException0 (ExcMETISNotInstalled);
                                     /**
                                      * Exception
                                      */
    DeclException1 (ExcInvalidNumberOfPartitions,
                    int,
                    << "The number of partitions you gave is " << arg1
                    << ", but must be greater than zero.");
                                     /**
                                      * Exception
                                      */
    DeclException2 (ExcInvalidArraySize,
                    int, int,
                    << "The array has size " << arg1 << " but should have size "
                    << arg2);
				     //@}
  private:
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
				      * The size of chunks.
				      */
    unsigned int chunk_size;

				     /**
				      * The reduced sparsity pattern. We store
				      * only which chunks exist, with each
				      * chunk a block in the matrix of size
				      * chunk_size by chunk_size.
				      */
    SparsityPattern sparsity_pattern;

				     /**
				      * Make all the chunk sparse matrix kinds
				      * friends.
				      */
    template <typename> friend class ChunkSparseMatrix;
};


/*@}*/
/*---------------------- Inline functions -----------------------------------*/

#ifndef DOXYGEN


inline
unsigned int
ChunkSparsityPattern::n_rows () const
{
  return rows;
}


inline
unsigned int
ChunkSparsityPattern::n_cols () const
{
  return cols;
}



inline
unsigned int
ChunkSparsityPattern::get_chunk_size () const
{
  return chunk_size;
}



inline
bool
ChunkSparsityPattern::is_compressed () const
{
  return sparsity_pattern.compressed;
}


inline
bool
ChunkSparsityPattern::optimize_diagonal () const
{
  return sparsity_pattern.diagonal_optimized;
}



template <typename ForwardIterator>
void
ChunkSparsityPattern::copy_from (const unsigned int    n_rows,
				 const unsigned int    n_cols,
				 const ForwardIterator begin,
				 const ForwardIterator end,
				 const unsigned int    chunk_size,
				 const bool optimize_diag)
{
  Assert (static_cast<unsigned int>(std::distance (begin, end)) == n_rows,
	  ExcIteratorRange (std::distance (begin, end), n_rows));

				   // first determine row lengths for
				   // each row. if the matrix is
				   // quadratic, then we might have to
				   // add an additional entry for the
				   // diagonal, if that is not yet
				   // present. as we have to call
				   // compress anyway later on, don't
				   // bother to check whether that
				   // diagonal entry is in a certain
				   // row or not
  const bool is_square = optimize_diag && (n_rows == n_cols);
  std::vector<unsigned int> row_lengths;
  row_lengths.reserve(n_rows);
  for (ForwardIterator i=begin; i!=end; ++i)
    row_lengths.push_back (std::distance (i->begin(), i->end())
			   +
			   (is_square ? 1 : 0));
  reinit (n_rows, n_cols, row_lengths, chunk_size, is_square);

				   // now enter all the elements into
				   // the matrix
  unsigned int row = 0;
  typedef typename std::iterator_traits<ForwardIterator>::value_type::const_iterator inner_iterator;
  for (ForwardIterator i=begin; i!=end; ++i, ++row)
    {
      const inner_iterator end_of_row = i->end();
      for (inner_iterator j=i->begin(); j!=end_of_row; ++j)
	{
	  const unsigned int col
	    = internal::SparsityPatternTools::get_column_index_from_iterator(*j);
	  Assert (col < n_cols, ExcInvalidIndex(col,n_cols));

	  add (row, col);
	}
    }

				   // finally compress
				   // everything. this also sorts the
				   // entries within each row
  compress ();
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
