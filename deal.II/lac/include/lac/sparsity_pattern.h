//----------------------------  sparsity_pattern.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparsity_pattern.h  ---------------------------
#ifndef __deal2__sparsity_pattern_h
#define __deal2__sparsity_pattern_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/subscriptor.h>

#include <vector>
#include <iostream>


class SparsityPattern;
template <typename number> class FullMatrix;
template <typename number> class SparseMatrix;
class CompressedSparsityPattern;



/*! @addtogroup Matrix1
 *@{
 */

namespace internals
{
  namespace SparsityPatternIterators
  {
                                     // forward declaration
    class Iterator;
    
                                     /**
                                      * Accessor class for iterators into
                                      * sparsity patterns. This class is also
                                      * the base class for both const and
                                      * non-const accessor classes into sparse
                                      * matrices.
                                      *
                                      * Note that this class only allow read
                                      * access to elements, providing their
                                      * row and column number. It does not
                                      * allow to modify the sparsity pattern
                                      * itself.
                                      */
    class Accessor
    {
      public:
                                         /**
                                          * Constructor.
                                          */
        Accessor (const SparsityPattern *matrix,
                  const unsigned int     row,
                  const unsigned int     index);

                                         /**
                                          * Constructor. Construct the end
                                          * accessor for the given sparsity
                                          * pattern.
                                          */
        Accessor (const SparsityPattern *matrix);

                                         /**
                                          * Row number of the element
                                          * represented by this object. This
                                          * function can only be called for
                                          * entries for which is_valid_entry()
                                          * is true.
                                          */
        unsigned int row () const;

                                         /**
                                          * Index in row of the element
                                          * represented by this object. This
                                          * function can only be called for
                                          * entries for which is_valid_entry()
                                          * is true.
                                          */
        unsigned int index () const;

                                         /**
                                          * Column number of the element
                                          * represented by this object. This
                                          * function can only be called for
                                          * entries for which is_valid_entry() is
                                          * true.
                                          */
        unsigned int column () const;

                                         /**
                                          * Return whether the sparsity
                                          * pattern entry pointed to by this
                                          * iterator is valid or not. Note
                                          * that after compressing the
                                          * sparsity pattern, all entries are
                                          * valid. However, before
                                          * compression, the sparsity pattern
                                          * allocated some memory to be used
                                          * while still adding new nonzero
                                          * entries; if you create iterators
                                          * in this phase of the sparsity
                                          * pattern's lifetime, you will
                                          * iterate over elements that are not
                                          * valid. If this is so, then this
                                          * function will return false.
                                          */
        inline bool is_valid_entry () const;


                                         /**
                                          * Comparison. True, if
                                          * both iterators point to
                                          * the same matrix
                                          * position.
                                          */
	bool operator == (const Accessor &) const;


                                         /**
                                          * Comparison
                                          * operator. Result is true
                                          * if either the first row
                                          * number is smaller or if
                                          * the row numbers are
                                          * equal and the first
                                          * index is smaller.
                                          *
                                          * This function is only valid if
                                          * both iterators point into the same
                                          * sparsity pattern.
                                          */
	bool operator < (const Accessor &) const;

                                         /**
                                          * Move the accessor to the next
                                          * nonzero entry in the matrix. This
                                          * function should actually only be
                                          * called from iterator classes, but
                                          * due to various problems with
                                          * friends and templates in gcc 2.95,
                                          * we can't make it protected. Don't
                                          * use it in your own programs.
                                          */
        void advance ();

                                         /**
                                          * Exception
                                          */
        DeclException0 (ExcInvalidIterator);
        
      protected:
                                         /**
                                          * The sparsity pattern we operate on
                                          * accessed.
                                          */
        const SparsityPattern * sparsity_pattern;

                                         /**
                                          * Current row number.
                                          */
        unsigned int a_row;

                                         /**
                                          * Current index in row.
                                          */
        unsigned int a_index;

                                         /**
                                          * Grant access to iterator class.
                                          */
        friend class Iterator;
    };



                                     /**
				      * STL conforming iterator walking over
				      * the elements of a sparsity pattern.
				      */
    class Iterator
    {
      public:
                                         /**
                                          * Constructor. Create an iterator
                                          * into the sparsity pattern @p sp for the
                                          * given row and the index within it.
                                          */ 
	Iterator (const SparsityPattern *sp,
                  const unsigned int     row,
                  const unsigned int     index);
        
                                         /**
                                          * Prefix increment.
                                          */
	Iterator& operator++ ();

                                         /**
                                          * Postfix increment.
                                          */
	Iterator operator++ (int);

                                         /**
                                          * Dereferencing operator.
                                          */
	const Accessor & operator* () const;

                                         /**
                                          * Dereferencing operator.
                                          */
	const Accessor * operator-> () const;

                                         /**
                                          * Comparison. True, if
                                          * both iterators point to
                                          * the same matrix
                                          * position.
                                          */
	bool operator == (const Iterator&) const;

                                         /**
                                          * Inverse of <tt>==</tt>.
                                          */
	bool operator != (const Iterator&) const;

                                         /**
                                          * Comparison
                                          * operator. Result is true
                                          * if either the first row
                                          * number is smaller or if
                                          * the row numbers are
                                          * equal and the first
                                          * index is smaller.
                                          *
                                          * This function is only valid if
                                          * both iterators point into the same
                                          * matrix.
                                          */
	bool operator < (const Iterator&) const;

      private:
                                         /**
                                          * Store an object of the
                                          * accessor class.
                                          */
        Accessor accessor;
    };
    
  }
}



/**
 * Structure representing the sparsity pattern of a sparse matrix.
 * 
 * The following picture will illustrate the relation between the
 * SparsityPattern an the SparseMatrix.
 *
 * @verbatim
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
 * @endverbatim
 *
 * @verbatim
 * For row = 0
 *   
 * there are: (0| 3) = colnums[0]
 *            (0| 2) = colnums[1]
 *            (0| 9) = colnums[2]
 *            (0|17) = colnums[3]
 *
 * For row = 1
 *   
 * there are: (1| 1) = colnums[4]
 *            (1| 4) = colnums[5]
 * ....
 *
 * @endverbatim
 *
 * @verbatim
 * SparseMatrix:                                  \
 *                                                 |
 *              _____________________________      |
 *  val        |  |  |  |  |  |  |  |  | 3|..       \ Value
 *             |__|__|__|__|__|__|__|__|__|__       /
 *                                                 |
 *                                                 |
 *                                                /
 * @endverbatim
 *
 * If you want to get the <tt>3</tt> you need to get its position in
 * the table above and its value by returning the value of the element
 * on which the pointer shows, using <tt>*val</tt>. For example
 * <tt>val[8]=3</tt>. Its position is <tt>colnums[8]</tt> so
 * <tt>row=2</tt>. In other words, if you want to get the element
 * <i>a<sub>24</sub></i> you know that <tt>row=2</tt>. To get the
 * element, a search of <tt>4</tt> form <tt>colnums[rowstart[2]]</tt>
 * to <tt>colnums[rowstart[3]]</tt> is needed. Then
 * <i>a<sub>24</sub></i>=<tt>val[number of the found element] =
 * 3</tt>.
 *
 *
 * @author Wolfgang Bangerth, Guido Kanschat and others
 */
class SparsityPattern : public Subscriptor
{
  public:
                                     /**
                                      * Typedef an iterator class that allows
                                      * to walk over all nonzero elements of a
                                      * sparsity pattern.
                                      */
    typedef
    internals::SparsityPatternIterators::Iterator
    const_iterator;
    
                                     /**
                                      * Typedef an iterator class that allows
                                      * to walk over all nonzero elements of a
                                      * sparsity pattern.
                                      *
                                      * Since the iterator does not allow to
                                      * modify the sparsity pattern, this type
                                      * is the same as that for @p
                                      * const_iterator.
                                      */
    typedef
    internals::SparsityPatternIterators::Iterator
    iterator;
    
    
				     /**
				      * Define a value which is used
				      * to indicate that a certain
				      * value in the #colnums array
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
    static const unsigned int invalid_entry = deal_II_numbers::invalid_unsigned_int;
    
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
    SparsityPattern ();
    
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
				      * (SparsityPattern());</tt>,
				      * with <tt>v</tt> a vector of
				      * SparsityPattern objects.
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
    SparsityPattern (const SparsityPattern &);

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
    SparsityPattern (const unsigned int m,
		     const unsigned int n,
		     const unsigned int max_per_row,
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
    SparsityPattern (const unsigned int               m,
		     const unsigned int               n,
		     const std::vector<unsigned int> &row_lengths,
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
    SparsityPattern (const unsigned int n,
		     const unsigned int max_per_row);

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
    SparsityPattern (const unsigned int               m,
		     const std::vector<unsigned int> &row_lengths,
		     const bool optimize_diagonal = true);

				     /**
				      * Make a copy with extra off-diagonals.
				      *
				      * This constructs objects intended for
				      * the application of the ILU(n)-method
				      * or other incomplete decompositions.
				      * Therefore, additional to the original
				      * entry structure, space for
				      * <tt>extra_off_diagonals</tt>
				      * side-diagonals is provided on both
				      * sides of the main diagonal.
				      *
				      * <tt>max_per_row</tt> is the
				      * maximum number of nonzero
				      * elements per row which this
				      * structure is to hold. It is
				      * assumed that this number is
				      * sufficiently large to
				      * accomodate both the elements
				      * in <tt>original</tt> as well
				      * as the new off-diagonal
				      * elements created by this
				      * constructor. You will usually
				      * want to give the same number
				      * as you gave for
				      * <tt>original</tt> plus the
				      * number of side diagonals times
				      * two. You may however give a
				      * larger value if you wish to
				      * add further nonzero entries
				      * for the decomposition based on
				      * other criteria than their
				      * being on side-diagonals.
				      *
				      * This function requires that
				      * <tt>original</tt> refers to a
				      * quadratic matrix structure.
				      * It must be compressed. The
				      * matrix structure is not
				      * compressed after this function
				      * finishes.
				      */
    SparsityPattern (const SparsityPattern  &original,
		     const unsigned int      max_per_row,
		     const unsigned int      extra_off_diagonals);
    
				     /**
				      * Destructor.
				      */
    ~SparsityPattern ();

				     /**
				      * Copy operator. For this the
				      * same holds as for the copy
				      * constructor: it is declared,
				      * defined and fine to be called,
				      * but the latter only for empty
				      * objects.
				      */
    SparsityPattern & operator = (const SparsityPattern &);
    
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
				      * IF the number of rows equals
				      * the number of columns and the
				      * last parameter is true,
				      * diagonals elements are stored
				      * first in each row to allow
				      * optimized access in relaxation
				      * methods of SparseMatrix.
				      */
    void reinit (const unsigned int               m,
		 const unsigned int               n,
		 const std::vector<unsigned int> &row_lengths,
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
				      * SparsityPattern objects they are
				      * initialized with to be compressed, to
				      * reduce memory requirements.
				      */
    void compress ();

				     /**
				      * STL-like iterator with the first entry
				      * of the matrix. The resulting iterator
				      * can be used to walk over all nonzero
				      * entries of the sparsity pattern.
				      */
    inline iterator begin () const;

				     /**
				      * Final iterator.
				      */
    inline iterator end () const;

				     /**
				      * STL-like iterator with the first entry
				      * of row <tt>r</tt>.
				      *
				      * Note that if the given row is empty,
				      * i.e. does not contain any nonzero
				      * entries, then the iterator returned by
				      * this function equals
				      * <tt>end(r)</tt>. Note also that the
				      * iterator may not be dereferencable in
				      * that case.
				      */
    inline iterator begin (const unsigned int r) const;

				     /**
				      * Final iterator of row <tt>r</tt>. It
				      * points to the first element past the
				      * end of line @p r, or past the end of
				      * the entire sparsity pattern.
				      *
				      * Note that the end iterator is not
				      * necessarily dereferencable. This is in
				      * particular the case if it is the end
				      * iterator for the last row of a matrix.
				      */
    inline iterator end (const unsigned int r) const;
    
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
		    const bool optimize_diagonal = true);

				     /**
				      * Copy data from an object of
				      * type
				      * CompressedSparsityPattern.
				      * Previous content of this
				      * object is lost, and the
				      * sparsity pattern is in
				      * compressed mode afterwards.
				      */
    void copy_from (const CompressedSparsityPattern &csp,
		    const bool optimize_diagonal = true);

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
		    const bool optimize_diagonal = true);
    
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
				      * element with row number <tt>i</tt>
				      * and column number <tt>j</tt>. If
				      * the matrix element is not a
				      * nonzero one, return
				      * SparsityPattern::invalid_entry.
				      *
				      * This function is usually
				      * called by the
				      * SparseMatrix::operator()(). It
				      * may only be called for
				      * compressed sparsity patterns,
				      * since in this case searching
				      * whether the entry exists can
				      * be done quite fast with a
				      * binary sort algorithm because
				      * the column numbers are sorted.
				      *
				      * If <tt>m</tt> is the number of
				      * entries in <tt>row</tt>, then the
				      * complexity of this function is
				      * <i>log(m)</i> if the sparsity
				      * pattern is compressed.
				      *
				      * @deprecated Use
				      * SparseMatrix::const_iterator
				      */
    unsigned int operator() (const unsigned int i, 
			     const unsigned int j) const;

				     /**
				      * This is the inverse operation
				      * to operator()(): given a
				      * global index, find out row and
				      * column of the matrix entry to
				      * which it belongs. The returned
				      * value is the pair composed of
				      * row and column index.
				      *
				      * This function may only be
				      * called if the sparsity pattern
				      * is closed. The global index
				      * must then be between zero and
				      * n_nonzero_elements().
				      *
				      * If <tt>N</tt> is the number of
				      * rows of this matrix, then the
				      * complexity of this function is
				      * <i>log(N)</i>.
				      */
    std::pair<unsigned int, unsigned int>
    matrix_position (const unsigned int global_index) const;
    
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
    bool exists (const unsigned int i, const unsigned int j) const;
    

				     /**
				      * Number of entries in a specific row.
				      */
    inline unsigned int row_length (const unsigned int row) const;

				     /**
				      * Access to column number field.
				      * Return the column number of
				      * the <tt>index</tt>th entry in
				      * <tt>row</tt>. Note that if
				      * diagonal elements are
				      * optimized, the first element
				      * in each row is the diagonal
				      * element,
				      * i.e. <tt>column_number(row,0)==row</tt>.
				      *
				      * If the sparsity pattern is
				      * already compressed, then
				      * (except for the diagonal
				      * element), the entries are
				      * sorted by columns,
				      * i.e. <tt>column_number(row,i)</tt>
				      * <tt><</tt> <tt>column_number(row,i+1)</tt>.
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
                                      * Use the METIS partitioner to generate
                                      * a partitioning of the degrees of
                                      * freedom represented by this sparsity
                                      * pattern. In effect, we view this
                                      * sparsity pattern as a graph of
                                      * connections between various degrees of
                                      * freedom, where each nonzero entry in
                                      * the sparsity pattern corresponds to an
                                      * edge between two nodes in the
                                      * connection graph. The goal is then to
                                      * decompose this graph into groups of
                                      * nodes so that a minimal number of
                                      * edges are cut by the boundaries
                                      * between node groups. This partitioning
                                      * is done by METIS. Note that METIS can
                                      * only partition symmetric sparsity
                                      * patterns, and that of course the
                                      * sparsity pattern has to be square. We
                                      * do not check for symmetry of the
                                      * sparsity pattern, since this is an
                                      * expensive operation, but rather leave
                                      * this as the responsibility of caller
                                      * of this function.
                                      *
                                      * After calling this function, the
                                      * output array will have values between
                                      * zero and @p n_partitions-1 for each
                                      * node (i.e. row or column of the
                                      * matrix).
                                      *
                                      * This function will generate an error
                                      * if METIS is not installed unless
                                      * @p n_partitions is one. I.e., you can
                                      * write a program so that it runs in the
                                      * single-processor single-partition case
                                      * without METIS installed, and only
                                      * requires METIS when multiple
                                      * partitions are required.
                                      *
                                      * Note that the sparsity pattern itself
                                      * is not changed by calling this
                                      * function. However, you will likely use
                                      * the information generated by calling
                                      * this function to renumber degrees of
                                      * freedom, after which you will of
                                      * course have to regenerate the sparsity
                                      * pattern.
                                      *
                                      * This function will rarely be called
                                      * separately, since in finite element
                                      * methods you will want to partition the
                                      * mesh, not the matrix. This can be done
                                      * by calling
                                      * @p GridTools::partition_triangulation.
                                      */
    void partition (const unsigned int         n_partitions,
                    std::vector<unsigned int> &partition_indices) const;
    
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
    unsigned int memory_consumption () const;

				     /**
				      * This is kind of an expert mode. Get 
				      * access to the rowstart array, but
				      * read-only.
				      *
				      * Use of this function is highly
				      * deprecated. Use @p row_length and
				      * @p column_number instead. Also, using
				      * iterators may get you most of the
				      * information you may want.
				      *
				      * Though the return value is declared
				      * <tt>const</tt>, you should be aware that it
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
    inline const unsigned int * get_rowstart_indices () const;

				     /**
				      * @deprecated. Use @p row_length and
				      * @p column_number instead. Also, using
				      * iterators may get you most of the
				      * information you may want.
				      *
				      * This is kind of an expert mode: get
				      * access to the colnums array, but
				      * readonly.
				      *
				      * Though the return value is declared
				      * <tt>const</tt>, you should be aware that it
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
    inline const unsigned int * get_column_numbers () const;

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
    DeclException0 (ExcNotQuadratic);
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
  private:
				     /**
				      * Maximum number of rows that can
				      * be stored in the #rowstart array.
				      * Since reallocation of that array
				      * only happens if the present one is
				      * too small, but never when the size
				      * of this matrix structure shrinks,
				      * #max_dim might be larger than
				      * #rows and in this case #rowstart
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
				      * #colnums. Here, the same applies as
				      * for the #rowstart array, i.e. it
				      * may be larger than the actually used
				      * part of the array.
				      */
    unsigned int max_vec_len;

				     /**
				      * Maximum number of elements per
				      * row. This is set to the value
				      * given to the reinit() function
				      * (or to the constructor), or to
				      * the maximum row length
				      * computed from the vectors in
				      * case the more flexible
				      * constructors or reinit
				      * versions are called. Its value
				      * is more or less meaningsless
				      * after compress() has been
				      * called.
				      */
    unsigned int max_row_length;

				     /**
				      * Array which hold for each row
				      * which is the first element in
				      * #colnums belonging to that
				      * row. Note that the size of the
				      * array is one larger than the
				      * number of rows, because the
				      * last element is used for
				      * <tt>row</tt>=#rows, i.e. the
				      * row past the last used
				      * one. The value of
				      * #rowstart[#rows]} equals the
				      * index of the element past the
				      * end in #colnums; this way, we
				      * are able to write loops like
				      * <tt>for (i=rowstart[k];
				      * i<rowstart[k+1]; ++i)</tt>
				      * also for the last row.
				      *
				      * Note that the actual size of the
				      * allocated memory may be larger than
				      * the region that is used. The actual
				      * number of elements that was allocated
				      * is stored in #max_dim.
				      */
    unsigned int *rowstart;

				     /**
				      * Array of column numbers. In
				      * this array, we store for each
				      * non-zero element its column
				      * number. The column numbers for
				      * the elements in row <i>r</i>
				      * are stored within the index
				      * range
				      * #rowstart[<i>r</i>]...#rowstart[<i>r+1</i>]. Therefore
				      * to find out whether a given
				      * element (<i>r,c</i>) exists,
				      * we have to check whether the
				      * column number <i>c</i> exists
				      * in the abovementioned range
				      * within this array. If it
				      * exists, say at position
				      * <i>p</i> within this array,
				      * the value of the respective
				      * element in the sparse matrix
				      * will also be at position
				      * <i>p</i> of the values array
				      * of that class.
				      *
				      * At the beginning, all elements
				      * of this array are set to
				      * @p -1 indicating invalid
				      * (unused) column numbers
				      * (diagonal elements are preset
				      * if optimized storage is
				      * requested, though). Now, if
				      * nonzero elements are added,
				      * one column number in the row's
				      * respective range after the
				      * other is set to the column
				      * number of the added
				      * element. When compress is
				      * called, unused elements
				      * (indicated by column numbers
				      * @p -1) are eliminated by
				      * copying the column number of
				      * subsequent rows and the column
				      * numbers within each row (with
				      * possible exception of the
				      * diagonal element) are sorted,
				      * such that finding whether an
				      * element exists and determining
				      * its position can be done by a
				      * binary search.
				      */
    unsigned int *colnums;

				     /**
				      * Store whether the compress()
				      * function was called for this
				      * object.
				      */
    bool compressed;

				     /**
				      * Is special treatment of
				      * diagonals enabled?
				      */
    bool diagonal_optimized;
    
				     /**
				      * Optimized replacement for
				      * <tt>std::lower_bound</tt> for
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
    const unsigned int *
    optimized_lower_bound (const unsigned int *first,
			   const unsigned int *last,
			   const unsigned int &val);

				     /**
				      * Helper function to get the
				      * column index from a
				      * dereferenced iterator in the
				      * copy_from() function, if
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
				      * copy_from() function, if
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
				      * (such as <tt>std::map</tt>).
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


/*@}*/
/*---------------------- Inline functions -----------------------------------*/

/// @if NoDoc


namespace internals
{
  namespace SparsityPatternIterators
  {    
    inline
    Accessor::
    Accessor (const SparsityPattern *sparsity_pattern,
              const unsigned int     r,
              const unsigned int     i)
                    :
                    sparsity_pattern(sparsity_pattern),
                    a_row(r),
                    a_index(i)
    {}


    inline
    Accessor::
    Accessor (const SparsityPattern *sparsity_pattern)
                    :
                    sparsity_pattern(sparsity_pattern),
                    a_row(sparsity_pattern->n_rows()),
                    a_index(0)
    {}


    inline
    unsigned int
    Accessor::row() const
    {
      Assert (is_valid_entry() == true, ExcInvalidIterator());

      return a_row;
    }


    inline
    unsigned int
    Accessor::column() const
    {
      Assert (is_valid_entry() == true, ExcInvalidIterator());
      
      return (sparsity_pattern
              ->get_column_numbers()[sparsity_pattern
                                     ->get_rowstart_indices()[a_row]+a_index]);
    }


    inline
    unsigned int
    Accessor::index() const
    {
      Assert (is_valid_entry() == true, ExcInvalidIterator());

      return a_index;
    }



    inline
    bool
    Accessor::is_valid_entry () const
    {
      return (sparsity_pattern
              ->get_column_numbers()[sparsity_pattern
                                     ->get_rowstart_indices()[a_row]+a_index]
              != SparsityPattern::invalid_entry);
    }



    inline
    bool
    Accessor::operator == (const Accessor& other) const
    {      
      return (sparsity_pattern  == other.sparsity_pattern &&
              a_row   == other.a_row &&
              a_index == other.a_index);
    }



    inline
    bool
    Accessor::operator < (const Accessor& other) const
    {      
      Assert (sparsity_pattern == other.sparsity_pattern,
              ExcInternalError());
      
      return (a_row < other.a_row ||
              (a_row == other.a_row &&
               a_index < other.a_index));
    }
    

    inline
    void
    Accessor::advance ()
    {      
      Assert (a_row < sparsity_pattern->n_rows(), ExcIteratorPastEnd());
  
      ++a_index;

                                       // if at end of line: cycle until we
                                       // find a row with a nonzero number of
                                       // entries
      while (a_index >= sparsity_pattern->row_length(a_row))
        {
          a_index = 0;
          ++a_row;

                                           // if we happened to find the end
                                           // of the matrix, then stop here
          if (a_row == sparsity_pattern->n_rows())
            break;
        }
    }
    


    inline
    Iterator::Iterator (const SparsityPattern *sparsity_pattern,
                        const unsigned int     r,
                        const unsigned int     i)
                    :
                    accessor(sparsity_pattern, r, i)
    {}



    inline
    Iterator &
    Iterator::operator++ ()
    {
      accessor.advance ();
      return *this;
    }



    inline
    Iterator
    Iterator::operator++ (int)
    {
      const Iterator iter = *this;
      accessor.advance ();
      return iter;
    }



    inline
    const Accessor &
    Iterator::operator* () const
    {
      return accessor;
    }



    inline
    const Accessor *
    Iterator::operator-> () const
    {
      return &accessor;
    }


    inline
    bool
    Iterator::operator == (const Iterator& other) const
    {
      return (accessor == other.accessor);
    }



    inline
    bool
    Iterator::operator != (const Iterator& other) const
    {
      return ! (*this == other);
    }


    inline
    bool
    Iterator::operator < (const Iterator& other) const
    {
      return accessor < other.accessor;
    }
    
  }
}




inline
SparsityPattern::iterator
SparsityPattern::begin () const
{
                                   // search for the first line with a nonzero
                                   // number of entries
  for (unsigned int r=0; r<n_rows(); ++r)
    if (row_length(r) > 0)
      return iterator(this, r, 0);

                                   // alright, this matrix is completely
                                   // empty. that's strange but ok. simply
                                   // return the end() iterator
  return end();
}


inline
SparsityPattern::iterator
SparsityPattern::end () const
{
  return iterator(this, n_rows(), 0);
}



inline
SparsityPattern::iterator
SparsityPattern::begin (const unsigned int r) const
{
  Assert (r<n_rows(), ExcIndexRange(r,0,n_rows()));

  if (row_length(r) > 0)
    return iterator(this, r, 0);
  else
    return end (r);
}



inline
SparsityPattern::iterator
SparsityPattern::end (const unsigned int r) const
{
  Assert (r<n_rows(), ExcIndexRange(r,0,n_rows()));

                                   // place the iterator on the first entry
                                   // past this line, or at the end of the
                                   // matrix
  for (unsigned int i=r+1; i<n_rows(); ++i)
    if (row_length(i) > 0)
      return iterator(this, i, 0);

                                   // if there is no such line, then take the
                                   // end iterator of the matrix
  return end();
}



inline
const unsigned int *
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
}


inline
unsigned int
SparsityPattern::n_cols () const
{
  return cols;
}


inline
bool
SparsityPattern::is_compressed () const
{
  return compressed;
}


inline
bool
SparsityPattern::optimize_diagonal () const
{
  return diagonal_optimized;
}


inline
const unsigned int *
SparsityPattern::get_rowstart_indices () const
{
  return rowstart;
}


inline
const unsigned int *
SparsityPattern::get_column_numbers () const
{
  return colnums;
}


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
}



inline
unsigned int
SparsityPattern::
get_column_index_from_iterator (const unsigned int i)
{
  return i;
}



template <typename value>
inline
unsigned int
SparsityPattern::
get_column_index_from_iterator (const std::pair<unsigned int, value> &i)
{
  return i.first;
}



template <typename value>
inline
unsigned int
SparsityPattern::
get_column_index_from_iterator (const std::pair<const unsigned int, value> &i)
{
  return i.first;
}



template <typename ForwardIterator>
void
SparsityPattern::copy_from (const unsigned int    n_rows,
			    const unsigned int    n_cols,
			    const ForwardIterator begin,
			    const ForwardIterator end,
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
  reinit (n_rows, n_cols, row_lengths, is_square);

				   // now enter all the elements into
				   // the matrix. note that if the
				   // matrix is quadratic, then we
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
}


/// @endif

#endif
