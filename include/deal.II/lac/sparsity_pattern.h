// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii__sparsity_pattern_h
#define dealii__sparsity_pattern_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>
#include <boost/serialization/array.hpp>
#include <boost/serialization/split_member.hpp>

#include <vector>
#include <iostream>

DEAL_II_NAMESPACE_OPEN

class SparsityPattern;
class ChunkSparsityPattern;
template <typename number> class FullMatrix;
template <typename number> class SparseMatrix;
template <typename number> class SparseLUDecomposition;
template <typename number> class SparseILU;
template <typename VectorType> class VectorSlice;

namespace ChunkSparsityPatternIterators
{
  class Accessor;
}


/*! @addtogroup Sparsity
 *@{
 */

namespace internals
{
  namespace SparsityPatternTools
  {
    /**
     * Declare type for container size.
     */
    typedef types::global_dof_index size_type;

    /**
     * Helper function to get the column index from a dereferenced iterator in
     * the copy_from() function, if the inner iterator type points to plain
     * unsigned integers.
     */
    size_type
    get_column_index_from_iterator (const size_type i);

    /**
     * Helper function to get the column index from a dereferenced iterator in
     * the copy_from() function, if the inner iterator type points to pairs of
     * unsigned integers and some other value.
     */
    template <typename value>
    size_type
    get_column_index_from_iterator (const std::pair<size_type, value> &i);

    /**
     * Likewise, but sometimes needed for certain types of containers that
     * make the first element of the pair constant (such as
     * <tt>std::map</tt>).
     */
    template <typename value>
    size_type
    get_column_index_from_iterator (const std::pair<const size_type, value> &i);

  }
}


/**
 * Iterators on objects of type SparsityPattern.
 */
namespace SparsityPatternIterators
{
  // forward declaration
  class Iterator;

  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Accessor class for iterators into sparsity patterns. This class is also
   * the base class for both const and non-const accessor classes into sparse
   * matrices.
   *
   * Note that this class only allows read access to elements, providing their
   * row and column number (or alternatively the index within the complete
   * sparsity pattern). It does not allow modifying the sparsity pattern
   * itself.
   *
   * @author Wolfgang Bangerth
   * @date 2004
   */
  class Accessor
  {
  public:
    /**
     * Constructor.
     */
    Accessor (const SparsityPattern *matrix,
              const std::size_t      index_within_sparsity);

    /**
     * Constructor. Construct the end accessor for the given sparsity pattern.
     */
    Accessor (const SparsityPattern *matrix);

    /**
     * Row number of the element represented by this object. This function can
     * only be called for entries for which is_valid_entry() is true.
     */
    size_type row () const;

    /**
     * Index within the current row of the element represented by this object.
     * This function can only be called for entries for which is_valid_entry()
     * is true.
     */
    size_type index () const;

    /**
     * Column number of the element represented by this object. This function
     * can only be called for entries for which is_valid_entry() is true.
     */
    size_type column () const;

    /**
     * Return whether the sparsity pattern entry pointed to by this iterator
     * is valid or not. Note that after compressing the sparsity pattern, all
     * entries are valid. However, before compression, the sparsity pattern
     * allocated some memory to be used while still adding new nonzero
     * entries; if you create iterators in this phase of the sparsity
     * pattern's lifetime, you will iterate over elements that are not valid.
     * If this is so, then this function will return false.
     */
    bool is_valid_entry () const;

    /**
     * Comparison. True, if both iterators point to the same matrix position.
     */
    bool operator == (const Accessor &) const;

    /**
     * Comparison operator. Result is true if either the first row number is
     * smaller or if the row numbers are equal and the first index is smaller.
     *
     * This function is only valid if both iterators point into the same
     * sparsity pattern.
     */
    bool operator < (const Accessor &) const;

  protected:
    /**
     * The sparsity pattern we operate on accessed.
     */
    const SparsityPattern *sparsity_pattern;

    /**
     * Index in global sparsity pattern. This index represents the location
     * the iterator/accessor points to within the array of the SparsityPattern
     * class that stores the column numbers. It is also the index within the
     * values array of a sparse matrix that stores the corresponding value of
     * this site.
     */
    std::size_t index_within_sparsity;

    /**
     * Move the accessor to the next nonzero entry in the matrix.
     */
    void advance ();

    /**
     * Grant access to iterator class.
     */
    friend class Iterator;

    /**
     * Grant access to accessor class of ChunkSparsityPattern.
     */
    friend class ChunkSparsityPatternIterators::Accessor;
  };



  /**
   * An iterator class for walking over the elements of a sparsity pattern.
   *
   * The typical use for these iterators is to iterate over the elements of a
   * sparsity pattern (or, since they also serve as the basis for iterating
   * over the elements of an associated matrix, over the elements of a sparse
   * matrix), or over the elements of individual rows. There is no guarantee
   * that the elements of a row are actually traversed in an order in which
   * column numbers monotonically increase. See the documentation of the
   * SparsityPattern class for more information.
   *
   * @note This class operates directly on the internal data structures of the
   * SparsityPattern class. As a consequence, some operations are cheap and
   * some are not. In particular, it is cheap to access the column index of
   * the sparsity pattern entry pointed to. On the other hand, it is expensive
   * to access the row index (this requires $O(\log(N))$ operations for a
   * matrix with $N$ row). As a consequence, when you design algorithms that
   * use these iterators, it is common practice to not loop over <i>all</i>
   * elements of a sparsity pattern at once, but to have an outer loop over
   * all rows and within this loop iterate over the elements of this row. This
   * way, you only ever need to dereference the iterator to obtain the column
   * indices whereas the (expensive) lookup of the row index can be avoided by
   * using the loop index instead.
   */
  class Iterator
  {
  public:
    /**
     * Constructor. Create an iterator into the sparsity pattern @p sp for the
     * given global index (i.e., the index of the given element counting from
     * the zeroth row).
     */
    Iterator (const SparsityPattern *sp,
              const std::size_t      index_within_sparsity);

    /**
     * Prefix increment.
     */
    Iterator &operator++ ();

    /**
     * Postfix increment.
     */
    Iterator operator++ (int);

    /**
     * Dereferencing operator.
     */
    const Accessor &operator* () const;

    /**
     * Dereferencing operator.
     */
    const Accessor *operator-> () const;

    /**
     * Comparison. True, if both iterators point to the same matrix position.
     */
    bool operator == (const Iterator &) const;

    /**
     * Inverse of <tt>==</tt>.
     */
    bool operator != (const Iterator &) const;

    /**
     * Comparison operator. Result is true if either the first row number is
     * smaller or if the row numbers are equal and the first index is smaller.
     *
     * This function is only valid if both iterators point into the same
     * matrix.
     */
    bool operator < (const Iterator &) const;

    /**
     * Return the distance between the current iterator and the argument. The
     * distance is given by how many times one has to apply operator++ to the
     * current iterator to get the argument (for a positive return value), or
     * operator-- (for a negative return value).
     */
    int operator - (const Iterator &p) const;

  private:
    /**
     * Store an object of the accessor class.
     */
    Accessor accessor;
  };
}



/**
 * A class that can store which elements of a matrix are nonzero (or, in fact,
 * <i>may</i> be nonzero) and for which we have to allocate memory to store
 * their values. This class is an example of the "static" type of sparsity
 * patters (see
 * @ref Sparsity).
 * It uses the <a
 * href="https://en.wikipedia.org/wiki/Sparse_matrix">compressed row storage
 * (CSR)</a> format to store data, and is used as the basis for the
 * SparseMatrix class.
 *
 * The elements of a SparsityPattern, corresponding to the places where
 * SparseMatrix objects can store nonzero entries, are stored row-by-row.
 * Within each row, elements are generally stored left-to-right in increasing
 * column index order; the exception to this rule is that if the matrix is
 * square (n_rows() == n_columns()), then the diagonal entry is stored as the
 * first element in each row to make operations like applying a Jacobi or SSOR
 * preconditioner faster. As a consequence, if you traverse the elements of a
 * row of a SparsityPattern with the help of iterators into this object (using
 * SparsityPattern::begin and SparsityPattern::end) you will find that the
 * elements are not sorted by column index within each row whenever the matrix
 * is square (the first item will be the diagonal, followed by the other
 * entries sorted by column index).
 *
 * @note While this class forms the basis upon which SparseMatrix objects base
 * their storage format, and thus plays a central role in setting up linear
 * systems, it is rarely set up directly due to the way it stores its
 * information. Rather, one typically goes through an intermediate format
 * first, see for example the step-2 tutorial program as well as the
 * documentation module
 * @ref Sparsity.
 *
 * @author Wolfgang Bangerth, Guido Kanschat and others
 */
class SparsityPattern : public Subscriptor
{
public:
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Typedef an iterator class that allows to walk over all nonzero elements
   * of a sparsity pattern.
   */
  typedef
  SparsityPatternIterators::Iterator
  const_iterator;

  /**
   * Typedef an iterator class that allows to walk over all nonzero elements
   * of a sparsity pattern.
   *
   * Since the iterator does not allow to modify the sparsity pattern, this
   * type is the same as that for @p const_iterator.
   */
  typedef
  SparsityPatternIterators::Iterator
  iterator;


  /**
   * Define a value which is used to indicate that a certain value in the
   * #colnums array is unused, i.e. does not represent a certain column number
   * index.
   *
   * Indices with this invalid value are used to insert new entries to the
   * sparsity pattern using the add() member function, and are removed when
   * calling compress().
   *
   * You should not assume that the variable declared here has a certain
   * value. The initialization is given here only to enable the compiler to
   * perform some optimizations, but the actual value of the variable may
   * change over time.
   */
  static const size_type invalid_entry = numbers::invalid_size_type;

  /**
   * @name Construction and setup Constructors, destructor; functions
   * initializing, copying and filling an object.
   */
// @{
  /**
   * Initialize the matrix empty, that is with no memory allocated. This is
   * useful if you want such objects as member variables in other classes. You
   * can make the structure usable by calling the reinit() function.
   */
  SparsityPattern ();

  /**
   * Copy constructor. This constructor is only allowed to be called if the
   * matrix structure to be copied is empty. This is so in order to prevent
   * involuntary copies of objects for temporaries, which can use large
   * amounts of computing time. However, copy constructors are needed if one
   * wants to place a SparsityPattern in a container, e.g., to write such
   * statements like <tt>v.push_back (SparsityPattern());</tt>, with
   * <tt>v</tt> a vector of SparsityPattern objects.
   *
   * Usually, it is sufficient to use the explicit keyword to disallow
   * unwanted temporaries, but this does not work for <tt>std::vector</tt>s.
   * Since copying a structure like this is not useful anyway because multiple
   * matrices can use the same sparsity structure, copies are only allowed for
   * empty objects, as described above.
   */
  SparsityPattern (const SparsityPattern &);

  /**
   * Initialize a rectangular pattern of size <tt>m x n</tt>.
   *
   * @param[in] m The number of rows.
   * @param[in] n The number of columns.
   * @param[in] max_per_row Maximum number of nonzero entries per row.
   */
  SparsityPattern (const size_type m,
                   const size_type n,
                   const unsigned int max_per_row);


  /**
   * Initialize a rectangular pattern of size <tt>m x n</tt>.
   *
   * @param[in] m The number of rows.
   * @param[in] n The number of columns.
   * @param[in] row_lengths Possible number of nonzero entries for each row.
   * This vector must have one entry for each row.
   */
  SparsityPattern (const size_type               m,
                   const size_type               n,
                   const std::vector<unsigned int> &row_lengths);

  /**
   * Initialize a quadratic pattern of dimension <tt>m</tt> with at most
   * <tt>max_per_row</tt> nonzero entries per row.
   *
   * This constructor automatically enables optimized storage of diagonal
   * elements. To avoid this, use the constructor taking row and column
   * numbers separately.
   */
  SparsityPattern (const size_type m,
                   const unsigned int max_per_row);

  /**
   * Initialize a quadratic pattern of size <tt>m x m</tt>.
   *
   * @param[in] m The number of rows and columns.
   * @param[in] row_lengths Maximum number of nonzero entries for each row.
   * This vector must have one entry for each row.
   */
  SparsityPattern (const size_type m,
                   const std::vector<unsigned int> &row_lengths);

  /**
   * Make a copy with extra off-diagonals.
   *
   * This constructs objects intended for the application of the ILU(n)-method
   * or other incomplete decompositions.  Therefore, additional to the
   * original entry structure, space for <tt>extra_off_diagonals</tt> side-
   * diagonals is provided on both sides of the main diagonal.
   *
   * <tt>max_per_row</tt> is the maximum number of nonzero elements per row
   * which this structure is to hold. It is assumed that this number is
   * sufficiently large to accommodate both the elements in <tt>original</tt>
   * as well as the new off-diagonal elements created by this constructor. You
   * will usually want to give the same number as you gave for
   * <tt>original</tt> plus the number of side diagonals times two. You may
   * however give a larger value if you wish to add further nonzero entries
   * for the decomposition based on other criteria than their being on side-
   * diagonals.
   *
   * This function requires that <tt>original</tt> refers to a quadratic
   * matrix structure.  It must be compressed. The matrix structure is not
   * compressed after this function finishes.
   */
  SparsityPattern (const SparsityPattern  &original,
                   const unsigned int        max_per_row,
                   const size_type        extra_off_diagonals);

  /**
   * Destructor.
   */
  ~SparsityPattern ();

  /**
   * Copy operator. For this the same holds as for the copy constructor: it is
   * declared, defined and fine to be called, but the latter only for empty
   * objects.
   */
  SparsityPattern &operator = (const SparsityPattern &);

  /**
   * Reallocate memory and set up data structures for a new matrix with <tt>m
   * </tt>rows and <tt>n</tt> columns, with at most <tt>max_per_row</tt>
   * nonzero entries per row.
   *
   * This function simply maps its operations to the other <tt>reinit</tt>
   * function.
   */
  void reinit (const size_type m,
               const size_type n,
               const unsigned int max_per_row);


  /**
   * Reallocate memory for a matrix of size <tt>m x n</tt>. The number of
   * entries for each row is taken from the array <tt>row_lengths</tt> which
   * has to give this number of each row <tt>i=1...m</tt>.
   *
   * If <tt>m*n==0</tt> all memory is freed, resulting in a total
   * reinitialization of the object. If it is nonzero, new memory is only
   * allocated if the new size extends the old one. This is done to save time
   * and to avoid fragmentation of the heap.
   *
   * If the number of rows equals the number of columns and the last parameter
   * is true, diagonal elements are stored first in each row to allow
   * optimized access in relaxation methods of SparseMatrix.
   */
  void reinit (const size_type               m,
               const size_type               n,
               const std::vector<unsigned int> &row_lengths);


  /**
   * Same as above, but with a VectorSlice argument instead.
   */
  void reinit (const size_type                                   m,
               const size_type                                   n,
               const VectorSlice<const std::vector<unsigned int> > &row_lengths);

  /**
   * This function compresses the sparsity structure that this object
   * represents.  It does so by eliminating unused entries and sorting the
   * remaining ones to allow faster access by usage of binary search
   * algorithms. A special sorting scheme is used for the diagonal entry of
   * quadratic matrices, which is always the first entry of each row.
   *
   * The memory which is no more needed is released.
   *
   * SparseMatrix objects require the SparsityPattern objects they are
   * initialized with to be compressed, to reduce memory requirements.
   */
  void compress ();


  /**
   * This function can be used as a replacement for reinit(), subsequent calls
   * to add() and a final call to close() if you know exactly in advance the
   * entries that will form the matrix sparsity pattern.
   *
   * The first two parameters determine the size of the matrix. For the two
   * last ones, note that a sparse matrix can be described by a sequence of
   * rows, each of which is represented by a sequence of pairs of column
   * indices and values. In the present context, the begin() and end()
   * parameters designate iterators (of forward iterator type) into a
   * container, one representing one row. The distance between begin() and
   * end() should therefore be equal to n_rows(). These iterators may be
   * iterators of <tt>std::vector</tt>, <tt>std::list</tt>, pointers into a
   * C-style array, or any other iterator satisfying the requirements of a
   * forward iterator. The objects pointed to by these iterators (i.e. what we
   * get after applying <tt>operator*</tt> or <tt>operator-></tt> to one of
   * these iterators) must be a container itself that provides functions
   * <tt>begin</tt> and <tt>end</tt> designating a range of iterators that
   * describe the contents of one line. Dereferencing these inner iterators
   * must either yield a pair of an unsigned integer as column index and a
   * value of arbitrary type (such a type would be used if we wanted to
   * describe a sparse matrix with one such object), or simply an unsigned
   * integer (of we only wanted to describe a sparsity pattern). The function
   * is able to determine itself whether an unsigned integer or a pair is what
   * we get after dereferencing the inner iterators, through some template
   * magic.
   *
   * While the order of the outer iterators denotes the different rows of the
   * matrix, the order of the inner iterator denoting the columns does not
   * matter, as they are sorted internal to this function anyway.
   *
   * Since that all sounds very complicated, consider the following example
   * code, which may be used to fill a sparsity pattern:
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
   * Note that this example works since the iterators dereferenced yield
   * containers with functions <tt>begin</tt> and <tt>end</tt> (namely
   * <tt>std::vector</tt>s), and the inner iterators dereferenced yield
   * unsigned integers as column indices. Note that we could have replaced
   * each of the two <tt>std::vector</tt> occurrences by <tt>std::list</tt>,
   * and the inner one by <tt>std::set</tt> as well.
   *
   * Another example would be as follows, where we initialize a whole matrix,
   * not only a sparsity pattern:
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
   * This example works because dereferencing iterators of the inner type
   * yields a pair of unsigned integers and a value, the first of which we
   * take as column index. As previously, the outer <tt>std::vector</tt> could
   * be replaced by <tt>std::list</tt>, and the inner <tt>std::map<unsigned
   * int,double></tt> could be replaced by <tt>std::vector<std::pair<unsigned
   * int,double> ></tt>, or a list or set of such pairs, as they all return
   * iterators that point to such pairs.
   */
  template <typename ForwardIterator>
  void copy_from (const size_type n_rows,
                  const size_type n_cols,
                  const ForwardIterator begin,
                  const ForwardIterator end);

  /**
   * Copy data from an object of type DynamicSparsityPattern. Although not a
   * compressed sparsity pattern, this function is also instantiated if the
   * argument is of type SparsityPattern (i.e., the current class). Previous
   * content of this object is lost, and the sparsity pattern is in compressed
   * mode afterwards.
   */
  template <typename SparsityPatternType>
  void copy_from (const SparsityPatternType &dsp);


  /**
   * Take a full matrix and use its nonzero entries to generate a sparse
   * matrix entry pattern for this object.
   *
   * Previous content of this object is lost, and the sparsity pattern is in
   * compressed mode afterwards.
   */
  template <typename number>
  void copy_from (const FullMatrix<number> &matrix);

  /**
   * Make the sparsity pattern symmetric by adding the sparsity pattern of the
   * transpose object.
   *
   * This function throws an exception if the sparsity pattern does not
   * represent a quadratic matrix.
   */
  void symmetrize ();

  /**
   * Add a nonzero entry to the matrix.  This function may only be called for
   * non-compressed sparsity patterns.
   *
   * If the entry already exists, nothing bad happens.
   */
  void add (const size_type i,
            const size_type j);

  /**
   * Add several nonzero entries to the specified matrix row.  This function
   * may only be called for non-compressed sparsity patterns.
   *
   * If some of the entries already exist, nothing bad happens.
   */
  template <typename ForwardIterator>
  void add_entries (const size_type row,
                    ForwardIterator begin,
                    ForwardIterator end,
                    const bool      indices_are_sorted = false);

// @}




  /**
   * @name Iterators
   */
// @{

  /**
   * Iterator starting at the first entry of the matrix. The resulting
   * iterator can be used to walk over all nonzero entries of the sparsity
   * pattern.
   *
   * Note the discussion in the general documentation of this class about the
   * order in which elements are accessed.
   */
  iterator begin () const;

  /**
   * Final iterator.
   */
  iterator end () const;

  /**
   * Iterator starting at the first entry of row <tt>r</tt>.
   *
   * Note that if the given row is empty, i.e. does not contain any nonzero
   * entries, then the iterator returned by this function equals
   * <tt>end(r)</tt>. Note also that the iterator may not be dereferencable in
   * that case.
   *
   * Note also the discussion in the general documentation of this class about
   * the order in which elements are accessed.
   */
  iterator begin (const size_type r) const;

  /**
   * Final iterator of row <tt>r</tt>. It points to the first element past the
   * end of line @p r, or past the end of the entire sparsity pattern.
   *
   * Note that the end iterator is not necessarily dereferencable. This is in
   * particular the case if it is the end iterator for the last row of a
   * matrix.
   */
  iterator end (const size_type r) const;


// @}
  /**
   * @name Querying information
   */
// @{
  /**
   * Test for equality of two SparsityPatterns.
   */
  bool operator == (const SparsityPattern &)  const;

  /**
   * Return whether the object is empty. It is empty if no memory is
   * allocated, which is the same as that both dimensions are zero.
   */
  bool empty () const;

  /**
   * Return the maximum number of entries per row. Before compression, this
   * equals the number given to the constructor, while after compression, it
   * equals the maximum number of entries actually allocated by the user.
   */
  size_type max_entries_per_row () const;

  /**
   * Compute the bandwidth of the matrix represented by this structure. The
   * bandwidth is the maximum of $|i-j|$ for which the index pair $(i,j)$
   * represents a nonzero entry of the matrix. Consequently, the maximum
   * bandwidth a $n\times m$ matrix can have is $\max\{n-1,m-1\}$, a diagonal
   * matrix has bandwidth 0, and there are at most $2*q+1$ entries per row if
   * the bandwidth is $q$. The returned quantity is sometimes called "half
   * bandwidth" in the literature.
   */
  size_type bandwidth () const;

  /**
   * Return the number of nonzero elements of this matrix. Actually, it
   * returns the number of entries in the sparsity pattern; if any of the
   * entries should happen to be zero, it is counted anyway.
   *
   * This function may only be called if the matrix struct is compressed. It
   * does not make too much sense otherwise anyway.
   */
  std::size_t n_nonzero_elements () const;

  /**
   * Return whether the structure is compressed or not.
   */
  bool is_compressed () const;

  /**
   * Return number of rows of this matrix, which equals the dimension of the
   * image space.
   */
  size_type n_rows () const;

  /**
   * Return number of columns of this matrix, which equals the dimension of
   * the range space.
   */
  size_type n_cols () const;

  /**
   * Number of entries in a specific row.
   */
  unsigned int row_length (const size_type row) const;

  /**
   * Return whether this object stores only those entries that have been added
   * explicitly, or if the sparsity pattern contains elements that have been
   * added through other means (implicitly) while building it. For the current
   * class, the result is false if and only if it is square because it then
   * unconditionally stores the diagonal entries whether they have been added
   * explicitly or not.
   *
   * This function mainly serves the purpose of describing the current class
   * in cases where several kinds of sparsity patterns can be passed as
   * template arguments.
   */
  bool stores_only_added_elements () const;

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object. See MemoryConsumption.
   */
  std::size_t memory_consumption () const;

// @}
  /**
   * @name Accessing entries
   */
// @{
  /**
   * Return the index of the matrix element with row number <tt>i</tt> and
   * column number <tt>j</tt>. If the matrix element is not a nonzero one,
   * return SparsityPattern::invalid_entry.
   *
   * This function is usually called by the SparseMatrix::operator()(). It may
   * only be called for compressed sparsity patterns, since in this case
   * searching whether the entry exists can be done quite fast with a binary
   * sort algorithm because the column numbers are sorted.
   *
   * If <tt>m</tt> is the number of entries in <tt>row</tt>, then the
   * complexity of this function is <i>log(m)</i> if the sparsity pattern is
   * compressed.
   *
   * @note This function is not cheap since it has to search through all of
   * the elements of the given row <tt>i</tt> to find whether index <tt>j</tt>
   * exists. Thus, it is more expensive than necessary in cases where you want
   * to loop over all of the nonzero elements of this sparsity pattern (or of
   * a sparse matrix associated with it) or of a single row. In such cases, it
   * is more efficient to use iterators over the elements of the sparsity
   * pattern or of the sparse matrix.
   */
  size_type operator() (const size_type i,
                        const size_type j) const;

  /**
   * This is the inverse operation to operator()(): given a global index, find
   * out row and column of the matrix entry to which it belongs. The returned
   * value is the pair composed of row and column index.
   *
   * This function may only be called if the sparsity pattern is closed. The
   * global index must then be between zero and n_nonzero_elements().
   *
   * If <tt>N</tt> is the number of rows of this matrix, then the complexity
   * of this function is <i>log(N)</i>.
   */
  std::pair<size_type, size_type>
  matrix_position (const std::size_t global_index) const;

  /**
   * Check if a value at a certain position may be non-zero.
   */
  bool exists (const size_type i,
               const size_type j) const;

  /**
   * The index of a global matrix entry in its row.
   *
   * This function is analogous to operator(), but it computes the index not
   * with respect to the total field, but only with respect to the row
   * <tt>j</tt>.
   */
  size_type row_position(const size_type i,
                         const size_type j) const;

  /**
   * Access to column number field.  Return the column number of the
   * <tt>index</tt>th entry in <tt>row</tt>. Note that if diagonal elements
   * are optimized, the first element in each row is the diagonal element,
   * i.e. <tt>column_number(row,0)==row</tt>.
   *
   * If the sparsity pattern is already compressed, then (except for the
   * diagonal element), the entries are sorted by columns, i.e.
   * <tt>column_number(row,i)</tt> <tt><</tt> <tt>column_number(row,i+1)</tt>.
   */
  size_type column_number (const size_type row,
                           const unsigned int index) const;


// @}
  /**
   * @name Input/Output
   */
// @{
  /**
   * Write the data of this object en bloc to a file. This is done in a binary
   * mode, so the output is neither readable by humans nor (probably) by other
   * computers using a different operating system or number format.
   *
   * The purpose of this function is that you can swap out matrices and
   * sparsity pattern if you are short of memory, want to communicate between
   * different programs, or allow objects to be persistent across different
   * runs of the program.
   */
  void block_write (std::ostream &out) const;

  /**
   * Read data that has previously been written by block_write() from a file.
   * This is done using the inverse operations to the above function, so it is
   * reasonably fast because the bitstream is not interpreted except for a few
   * numbers up front.
   *
   * The object is resized on this operation, and all previous contents are
   * lost.
   *
   * A primitive form of error checking is performed which will recognize the
   * bluntest attempts to interpret some data as a vector stored bitwise to a
   * file, but not more.
   */
  void block_read (std::istream &in);

  /**
   * Print the sparsity of the matrix. The output consists of one line per row
   * of the format <tt>[i,j1,j2,j3,...]</tt>. <i>i</i> is the row number and
   * <i>jn</i> are the allocated columns in this row.
   */
  void print (std::ostream &out) const;

  /**
   * Print the sparsity of the matrix in a format that <tt>gnuplot</tt>
   * understands and which can be used to plot the sparsity pattern in a
   * graphical way. The format consists of pairs <tt>i j</tt> of nonzero
   * elements, each representing one entry of this matrix, one per line of the
   * output file. Indices are counted from zero on, as usual. Since sparsity
   * patterns are printed in the same way as matrices are displayed, we print
   * the negative of the column index, which means that the <tt>(0,0)</tt>
   * element is in the top left rather than in the bottom left corner.
   *
   * Print the sparsity pattern in gnuplot by setting the data style to dots
   * or points and use the <tt>plot</tt> command.
   */
  void print_gnuplot (std::ostream &out) const;

  /**
   * Prints the sparsity of the matrix in a .svg file which can be opened in a
   * web browser. The .svg file contains squares which correspond to the
   * entries in the matrix. An entry in the matrix which contains a non-zero
   * value corresponds with a red square while a zero-valued entry in the
   * matrix correspond with a white square.
   */
  void print_svg (std::ostream &out) const;


  /**
   * Write the data of this object to a stream for the purpose of
   * serialization
   */
  template <class Archive>
  void save (Archive &ar, const unsigned int version) const;

  /**
   * Read the data of this object from a stream for the purpose of
   * serialization
   */
  template <class Archive>
  void load (Archive &ar, const unsigned int version);

  BOOST_SERIALIZATION_SPLIT_MEMBER()

// @}

  /**
   * @addtogroup Exceptions
   * @{
   */
  /**
   * You tried to add an element to a row, but there was no space left.
   */
  DeclException2 (ExcNotEnoughSpace,
                  int, int,
                  << "Upon entering a new entry to row " << arg1
                  << ": there was no free entry any more. " << std::endl
                  << "(Maximum number of entries for this row: "
                  << arg2 << "; maybe the matrix is already compressed?)");
  /**
   * The operation is only allowed after the SparsityPattern has been set up
   * and compress() was called.
   */
  DeclException0 (ExcNotCompressed);
  /**
   * This operation changes the structure of the SparsityPattern and is not
   * possible after compress() has been called.
   */
  DeclException0 (ExcMatrixIsCompressed);
  /**
   * Exception
   */
  DeclException0 (ExcInvalidConstructorCall);
  /**
   * This exception is thrown if the matrix does not follow the convention of
   * storing diagonal elements first in row. Refer to
   * SparityPattern::optimize_diagonal() for more information.
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
  DeclException1 (ExcInvalidNumberOfPartitions,
                  int,
                  << "The number of partitions you gave is " << arg1
                  << ", but must be greater than zero.");
  //@}
private:
  /**
   * Maximum number of rows that can be stored in the #rowstart array.  Since
   * reallocation of that array only happens if the present one is too small,
   * but never when the size of this matrix structure shrinks, #max_dim might
   * be larger than #rows and in this case #rowstart has more elements than
   * are used.
   */
  size_type max_dim;

  /**
   * Number of rows that this sparsity structure shall represent.
   */
  size_type rows;

  /**
   * Number of columns that this sparsity structure shall represent.
   */
  size_type cols;

  /**
   * Size of the actually allocated array #colnums. Here, the same applies as
   * for the #rowstart array, i.e. it may be larger than the actually used
   * part of the array.
   */
  std::size_t max_vec_len;

  /**
   * Maximum number of elements per row. This is set to the value given to the
   * reinit() function (or to the constructor), or to the maximum row length
   * computed from the vectors in case the more flexible constructors or
   * reinit versions are called. Its value is more or less meaningless after
   * compress() has been called.
   */
  unsigned int max_row_length;

  /**
   * Array which hold for each row which is the first element in #colnums
   * belonging to that row. Note that the size of the array is one larger than
   * the number of rows, because the last element is used for
   * <tt>row</tt>=#rows, i.e. the row past the last used one. The value of
   * #rowstart[#rows]} equals the index of the element past the end in
   * #colnums; this way, we are able to write loops like <tt>for
   * (i=rowstart[k]; i<rowstart[k+1]; ++i)</tt> also for the last row.
   *
   * Note that the actual size of the allocated memory may be larger than the
   * region that is used. The actual number of elements that was allocated is
   * stored in #max_dim.
   */
  std::size_t *rowstart;

  /**
   * Array of column numbers. In this array, we store for each non-zero
   * element its column number. The column numbers for the elements in row
   * <i>r</i> are stored within the index range
   * #rowstart[<i>r</i>]...#rowstart[<i>r+1</i>]. Therefore to find out
   * whether a given element (<i>r,c</i>) exists, we have to check whether the
   * column number <i>c</i> exists in the above-mentioned range within this
   * array. If it exists, say at position <i>p</i> within this array, the
   * value of the respective element in the sparse matrix will also be at
   * position <i>p</i> of the values array of that class.
   *
   * At the beginning, all elements of this array are set to @p -1 indicating
   * invalid (unused) column numbers (diagonal elements are preset if
   * optimized storage is requested, though). Now, if nonzero elements are
   * added, one column number in the row's respective range after the other is
   * set to the column number of the added element. When compress is called,
   * unused elements (indicated by column numbers @p -1) are eliminated by
   * copying the column number of subsequent rows and the column numbers
   * within each row (with possible exception of the diagonal element) are
   * sorted, such that finding whether an element exists and determining its
   * position can be done by a binary search.
   */
  size_type *colnums;

  /**
   * Store whether the compress() function was called for this object.
   */
  bool compressed;

  /**
   * Is special treatment of diagonals enabled?
   */
  bool store_diagonal_first_in_row;

  /**
   * Make all sparse matrices friends of this class.
   */
  template <typename number> friend class SparseMatrix;
  template <typename number> friend class SparseLUDecomposition;
  template <typename number> friend class SparseILU;
  template <typename number> friend class ChunkSparseMatrix;

  friend class ChunkSparsityPattern;

  /**
   * Also give access to internal details to the iterator/accessor classes.
   */
  friend class SparsityPatternIterators::Iterator;
  friend class SparsityPatternIterators::Accessor;
  friend class ChunkSparsityPatternIterators::Accessor;
};


/*@}*/
/*---------------------- Inline functions -----------------------------------*/

#ifndef DOXYGEN


namespace SparsityPatternIterators
{
  inline
  Accessor::
  Accessor (const SparsityPattern *sparsity_pattern,
            const std::size_t      i)
    :
    sparsity_pattern(sparsity_pattern),
    index_within_sparsity(i)
  {}


  inline
  Accessor::
  Accessor (const SparsityPattern *sparsity_pattern)
    :
    sparsity_pattern(sparsity_pattern),
    index_within_sparsity(sparsity_pattern->rowstart[sparsity_pattern->rows])
  {}


  inline
  bool
  Accessor::is_valid_entry () const
  {
    return (index_within_sparsity < sparsity_pattern->rowstart[sparsity_pattern->rows]
            &&
            sparsity_pattern->colnums[index_within_sparsity]
            != SparsityPattern::invalid_entry);
  }


  inline
  size_type
  Accessor::row() const
  {
    Assert (is_valid_entry() == true, ExcInvalidIterator());

    const std::size_t *insert_point =
      std::upper_bound(sparsity_pattern->rowstart,
                       sparsity_pattern->rowstart + sparsity_pattern->rows + 1,
                       index_within_sparsity);
    return insert_point - sparsity_pattern->rowstart - 1;
  }


  inline
  size_type
  Accessor::column() const
  {
    Assert (is_valid_entry() == true, ExcInvalidIterator());

    return (sparsity_pattern->colnums[index_within_sparsity]);
  }


  inline
  size_type
  Accessor::index() const
  {
    Assert (is_valid_entry() == true, ExcInvalidIterator());

    return index_within_sparsity - sparsity_pattern->rowstart[row()];
  }




  inline
  bool
  Accessor::operator == (const Accessor &other) const
  {
    return (sparsity_pattern == other.sparsity_pattern &&
            index_within_sparsity == other.index_within_sparsity);
  }



  inline
  bool
  Accessor::operator < (const Accessor &other) const
  {
    Assert (sparsity_pattern == other.sparsity_pattern,
            ExcInternalError());

    return index_within_sparsity < other.index_within_sparsity;
  }


  inline
  void
  Accessor::advance ()
  {
    Assert (index_within_sparsity < sparsity_pattern->rowstart[sparsity_pattern->rows],
            ExcIteratorPastEnd());
    ++index_within_sparsity;
  }



  inline
  Iterator::Iterator (const SparsityPattern *sparsity_pattern,
                      const std::size_t      i)
    :
    accessor(sparsity_pattern, i)
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
  Iterator::operator == (const Iterator &other) const
  {
    return (accessor == other.accessor);
  }



  inline
  bool
  Iterator::operator != (const Iterator &other) const
  {
    return ! (*this == other);
  }


  inline
  bool
  Iterator::operator < (const Iterator &other) const
  {
    return accessor < other.accessor;
  }


  inline
  int
  Iterator::operator - (const Iterator &other) const
  {
    Assert (accessor.sparsity_pattern == other.accessor.sparsity_pattern,
            ExcInternalError());

    return (*this)->index_within_sparsity - other->index_within_sparsity;
  }
}



inline
SparsityPattern::iterator
SparsityPattern::begin () const
{
  return iterator(this, rowstart[0]);
}


inline
SparsityPattern::iterator
SparsityPattern::end () const
{
  return iterator(this, rowstart[rows]);
}



inline
SparsityPattern::iterator
SparsityPattern::begin (const size_type r) const
{
  Assert (r<n_rows(), ExcIndexRangeType<size_type>(r,0,n_rows()));

  return iterator(this, rowstart[r]);
}



inline
SparsityPattern::iterator
SparsityPattern::end (const size_type r) const
{
  Assert (r<n_rows(), ExcIndexRangeType<size_type>(r,0,n_rows()));

  return iterator(this, rowstart[r+1]);
}



inline
SparsityPattern::size_type
SparsityPattern::n_rows () const
{
  return rows;
}


inline
SparsityPattern::size_type
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
SparsityPattern::stores_only_added_elements () const
{
  return (store_diagonal_first_in_row == false);
}



inline
unsigned int
SparsityPattern::row_length (const size_type row) const
{
  Assert(row<rows, ExcIndexRangeType<size_type>(row,0,rows));
  return rowstart[row+1]-rowstart[row];
}



inline
SparsityPattern::size_type
SparsityPattern::column_number (const size_type row,
                                const unsigned int index) const
{
  Assert(row<rows, ExcIndexRangeType<size_type>(row,0,rows));
  Assert(index<row_length(row), ExcIndexRange(index,0,row_length(row)));

  return colnums[rowstart[row]+index];
}


inline
std::size_t
SparsityPattern::n_nonzero_elements () const
{
  Assert ((rowstart!=0) && (colnums!=0), ExcEmptyObject());
  Assert (compressed, ExcNotCompressed());
  return rowstart[rows]-rowstart[0];
}



template <class Archive>
inline
void
SparsityPattern::save (Archive &ar, const unsigned int) const
{
  // forward to serialization function in the base class.
  ar   &static_cast<const Subscriptor &>(*this);

  ar &max_dim &rows &cols &max_vec_len &max_row_length &compressed &store_diagonal_first_in_row;

  ar &boost::serialization::make_array(rowstart, max_dim + 1);
  ar &boost::serialization::make_array(colnums, max_vec_len);
}



template <class Archive>
inline
void
SparsityPattern::load (Archive &ar, const unsigned int)
{
  // forward to serialization function in the base class.
  ar   &static_cast<Subscriptor &>(*this);

  ar &max_dim &rows &cols &max_vec_len &max_row_length &compressed &store_diagonal_first_in_row;

  rowstart = new std::size_t [max_dim + 1];
  colnums = new size_type [max_vec_len];

  ar &boost::serialization::make_array(rowstart, max_dim + 1);
  ar &boost::serialization::make_array(colnums, max_vec_len);
}



inline
bool
SparsityPattern::operator == (const SparsityPattern &sp2)  const
{
  // it isn't quite necessary to compare *all* member variables. by only
  // comparing the essential ones, we can say that two sparsity patterns are
  // equal even if one is compressed and the other is not (in which case some
  // of the member variables are not yet set correctly)
  if (rows != sp2.rows ||
      cols != sp2.cols ||
      compressed != sp2.compressed ||
      store_diagonal_first_in_row != sp2.store_diagonal_first_in_row)
    return false;

  for (size_type i = 0; i < rows+1; ++i)
    if (rowstart[i] != sp2.rowstart[i])
      return false;

  for (size_type i = 0; i < rowstart[rows]; ++i)
    if (colnums[i] != sp2.colnums[i])
      return false;

  return true;
}



namespace internal
{
  namespace SparsityPatternTools
  {
    /**
     * Declare type for container size.
     */
    typedef types::global_dof_index size_type;

    inline
    size_type
    get_column_index_from_iterator (const size_type i)
    {
      return i;
    }



    template <typename value>
    inline
    size_type
    get_column_index_from_iterator (const std::pair<size_type, value> &i)
    {
      return i.first;
    }



    template <typename value>
    inline
    size_type
    get_column_index_from_iterator (const std::pair<const size_type, value> &i)
    {
      return i.first;
    }
  }
}



template <typename ForwardIterator>
void
SparsityPattern::copy_from (const size_type       n_rows,
                            const size_type       n_cols,
                            const ForwardIterator begin,
                            const ForwardIterator end)
{
  Assert (static_cast<size_type>(std::distance (begin, end)) == n_rows,
          ExcIteratorRange (std::distance (begin, end), n_rows));

  // first determine row lengths for each row. if the matrix is quadratic,
  // then we might have to add an additional entry for the diagonal, if that
  // is not yet present. as we have to call compress anyway later on, don't
  // bother to check whether that diagonal entry is in a certain row or not
  const bool is_square = (n_rows == n_cols);
  std::vector<unsigned int> row_lengths;
  row_lengths.reserve(n_rows);
  for (ForwardIterator i=begin; i!=end; ++i)
    row_lengths.push_back (std::distance (i->begin(), i->end())
                           +
                           (is_square ? 1 : 0));
  reinit (n_rows, n_cols, row_lengths);

  // now enter all the elements into the matrix. note that if the matrix is
  // quadratic, then we already have the diagonal element preallocated
  //
  // for use in the inner loop, we define a typedef to the type of the inner
  // iterators
  size_type row = 0;
  typedef typename std::iterator_traits<ForwardIterator>::value_type::const_iterator inner_iterator;
  for (ForwardIterator i=begin; i!=end; ++i, ++row)
    {
      size_type *cols = &colnums[rowstart[row]] + (is_square ? 1 : 0);
      const inner_iterator end_of_row = i->end();
      for (inner_iterator j=i->begin(); j!=end_of_row; ++j)
        {
          const size_type col
            = internal::SparsityPatternTools::get_column_index_from_iterator(*j);
          Assert (col < n_cols, ExcIndexRange(col,0,n_cols));

          if ((col!=row) || !is_square)
            *cols++ = col;
        }
    }

  // finally compress everything. this also sorts the entries within each row
  compress ();
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
