//----------------------------  compressed_sparsity_pattern.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  compressed_sparsity_pattern.h  ---------------------------
#ifndef __deal2__compressed_sparsity_pattern_h
#define __deal2__compressed_sparsity_pattern_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/subscriptor.h>

template <typename number> class SparseMatrix;

#include <vector>
#include <set>


/**
 * This class acts as an intermediate form of the
 * @ref{SparsityPattern} class. From the interface it mostly
 * represents a @ref{SparsityPattern} object that is kept compressed
 * at all times. However, since the final sparsity pattern is not
 * known while constructing it, keeping the pattern compressed at all
 * times can only be achieved at the expense of either increased
 * memory or run time consumption upon use. The main purpose of this
 * class is to avoid some memory bottlenecks, so we chose to implement
 * it memory conservative, but the chosen data format is too unsuited
 * to be used for actual matrices. It is therefore necessary to first
 * copy the data of this object over to an object of type
 * @ref{SparsityPattern} before using it in actual matrices.
 *
 * Another viewpoint is that this class does not need up front
 * allocation of a certain amount of memory, but grows as necessary.
 *
 *
 * @sect3{Rationale}
 *
 * When constructing the sparsity pattern of a matrix, you usually
 * first have to provide an empty sparsity pattern object with a fixed
 * maximal number of entries per row. To find out about this maximal
 * row length, one usually calls the function
 * @ref{DoFHandler}@p{::max_couplings_per_dof} which returns an
 * estimate for that quantity. While this estimate is usually quite
 * good in 2d and exact in 1d, it is often significantly too large in
 * 3d and especially for higher order elements. Furthermore, normally
 * only a small fraction of the rows of a matrix will end up having
 * the maximal number of nonzero entries per row (usually those nodes
 * adjacent to hanging nodes), most have much less. In effect, the
 * empty @ref{SparsityPattern} object has allocated much too much
 * memory. Although this unnecessarily allocated memory is later freed
 * when @ref{SparsityPattern}@p{::compress} is called, this
 * overallocation has, with higher order elements and in 3d, sometimes
 * been so large that the program aborted due to lack of memory.
 *
 * This class therefore provides an alternative representation of a
 * sparsity pattern: we don't specify a maximal row length initially,
 * but store a set of column indices indicating possible nonzero
 * entries in the sparsity pattern for each row. This is very much
 * like the final "compressed" format used in the
 * @ref{SparsityPattern} object after compression, but uses a less
 * compact memory storage format, since the exact number of entries
 * per row is only known a posteriori and since it may change (for the
 * @ref{SparsityPattern} class, no more changes are allowed after
 * compressing it). We can therefore not store all the column indices
 * in a big array, but have to use a vector of sets. This can later be
 * used to actually initialize a @ref{SparsityPattern} object with the
 * then final set of necessary indices.
 *
 *
 * @sect3{Interface}
 *
 * Since this class is intended as an intermediate replacement of the
 * @ref{SparsityPattern} class, it has mostly the same interface, with
 * small changes where necessary. In particular, the @ref{add}
 * function, and the functions inquiring properties of the sparsity
 * pattern are the same.
 *
 *
 * @sect3{Usage}
 *
 * Use this class as follows:
 * @begin{verbatim}
 * CompressedSparsityPattern compressed_pattern (dof_handler.n_dofs());
 * DoFTools::make_sparsity_pattern (dof_handler,
 *                                  compressed_pattern);
 * constraints.condense (compressed_pattern);
 *
 * SparsityPattern sp;
 * sp.copy_from (compressed_pattern);
 * @end{verbatim}
 *
 *
 * @author Wolfgang Bangerth, 2001
 */
class CompressedSparsityPattern : public Subscriptor
{
  public:
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
    CompressedSparsityPattern ();
    
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
				      * @p{v.push_back (CompressedSparsityPattern());},
				      * with @p{v} a vector of @p{CompressedSparsityPattern}
				      * objects.
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
    CompressedSparsityPattern (const CompressedSparsityPattern &);

				     /**
				      * Initialize a rectangular
				      * matrix with @p{m} rows and
				      * @p{n} columns.
				      */
    CompressedSparsityPattern (const unsigned int m,
			       const unsigned int n);
    
				     /**
				      * Initialize a square matrix of
				      * dimension @p{n}.
				      */
    CompressedSparsityPattern (const unsigned int n);

				     /**
				      * Copy operator. For this the
				      * same holds as for the copy
				      * constructor: it is declared,
				      * defined and fine to be called,
				      * but the latter only for empty
				      * objects.
				      */
    CompressedSparsityPattern & operator = (const CompressedSparsityPattern &);
    
				     /**
				      * Reallocate memory and set up
				      * data structures for a new
				      * matrix with @p{m} rows and
				      * @p{n} columns, with at most
				      * @p{max_per_row} nonzero
				      * entries per row.
				      */
    void reinit (const unsigned int m,
		 const unsigned int n);
    
				     /**
				      * Since this object is kept
				      * compressed at all times anway,
				      * this function does nothing,
				      * but is declared to make the
				      * interface of this class as
				      * much alike as that of the
				      * @ref{SparsityPattern} class.
				      */
    void compress ();
    
				     /**
				      * Return whether the object is
				      * empty. It is empty if no
				      * memory is allocated, which is
				      * the same as that both
				      * dimensions are zero.
				      */
    bool empty () const;

				     /**
				      * Return the maximum number of
				      * entries per row. Note that
				      * this number may change as
				      * entries are added.
				      */
    unsigned int max_entries_per_row () const;

				     /**
				      * Add a nonzero entry to the
				      * matrix. If the entry already
				      * exists, nothing bad happens.
				      */
    void add (const unsigned int i, 
	      const unsigned int j);
    
				     /**
				      * Check if a value at a certain
				      * position may be non-zero.
				      */
    bool exists (const unsigned int i, const unsigned int j) const;
    
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
				      * Exception
				      */
    DeclException2 (ExcInvalidIndex,
		    int, int,
		    << "The given index " << arg1
		    << " should be less than " << arg2 << ".");
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
				      * Actual data: store for each
				      * row the set of nonzero
				      * entries.
				      */
    std::vector<std::set<unsigned int> > column_indices;

    friend class SparsityPattern;
};


/*---------------------- Inline functions -----------------------------------*/


inline
unsigned int
CompressedSparsityPattern::n_rows () const
{
  return rows;
};


inline
unsigned int
CompressedSparsityPattern::n_cols () const
{
  return cols;
};

#endif
