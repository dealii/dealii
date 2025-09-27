// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2000 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_sparsity_pattern_h
#define dealii_sparsity_pattern_h


#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/linear_index_iterator.h>

#include <deal.II/lac/sparsity_pattern_base.h>

// boost::serialization::make_array used to be in array.hpp, but was
// moved to a different file in BOOST 1.64
#include <boost/version.hpp>
#if BOOST_VERSION >= 106400
#  include <boost/serialization/array_wrapper.hpp>
#else
#  include <boost/serialization/array.hpp>
#endif
#include <boost/serialization/split_member.hpp>

#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
class SparsityPattern;
class DynamicSparsityPattern;
class ChunkSparsityPattern;
template <typename number>
class FullMatrix;
template <typename number>
class SparseMatrix;
template <typename number>
class SparseLUDecomposition;
template <typename number>
class SparseILU;

namespace ChunkSparsityPatternIterators
{
  class Accessor;
}
#endif

/**
 * @addtogroup Sparsity
 * @{
 */

namespace internals
{
  namespace SparsityPatternTools
  {
    /**
     * Declare type for container size.
     */
    using size_type = types::global_dof_index;

    /**
     * Helper function to get the column index from a dereferenced iterator in
     * the copy_from() function, if the inner iterator type points to plain
     * unsigned integers.
     */
    size_type
    get_column_index_from_iterator(const size_type i);

    /**
     * Helper function to get the column index from a dereferenced iterator in
     * the copy_from() function, if the inner iterator type points to pairs of
     * unsigned integers and some other value.
     */
    template <typename value>
    size_type
    get_column_index_from_iterator(const std::pair<size_type, value> &i);

    /**
     * Likewise, but sometimes needed for certain types of containers that
     * make the first element of the pair constant (such as
     * <tt>std::map</tt>).
     */
    template <typename value>
    size_type
    get_column_index_from_iterator(const std::pair<const size_type, value> &i);

  } // namespace SparsityPatternTools
} // namespace internals


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
  using size_type = types::global_dof_index;

  /**
   * Accessor class for iterators into sparsity patterns. This class is also
   * the base class for both const and non-const accessor classes into sparse
   * matrices.
   *
   * Note that this class only allows read access to elements, providing their
   * row and column number (or alternatively the index within the complete
   * sparsity pattern). It does not allow modifying the sparsity pattern
   * itself.
   */
  class Accessor
  {
  public:
    /**
     * Size type of SparsityPattern.
     */
    using size_type = SparsityPatternIterators::size_type;

    /**
     * Constructor.
     */
    Accessor(const SparsityPattern *matrix, const std::size_t linear_index);

    /**
     * Constructor. Construct the end accessor for the given sparsity pattern.
     */
    Accessor(const SparsityPattern *matrix);

    /**
     * Default constructor creating a dummy accessor. This constructor is here
     * only to be able to store accessors in STL containers such as
     * `std::vector`.
     */
    Accessor();

    /**
     * Row number of the element represented by this object. This function can
     * only be called for entries for which is_valid_entry() is true.
     */
    size_type
    row() const;

    /**
     * Index within the current row of the element represented by this object.
     * This function can only be called for entries for which is_valid_entry()
     * is true.
     */
    size_type
    index() const;

    /**
     * This function returns the how-many'th entry within the entire sparsity
     * pattern the current iterator points to. While the order of entries in
     * a sparsity pattern is generally not important, this function allows
     * indexing entries of the sparsity pattern using a linear index.
     *
     * This function can only be called for entries for which is_valid_entry()
     * is true.
     */
    size_type
    global_index() const;

    /**
     * Column number of the element represented by this object. This function
     * can only be called for entries for which is_valid_entry() is true.
     */
    size_type
    column() const;

    /**
     * Return whether the sparsity pattern entry pointed to by this iterator
     * is valid or not. Note that after compressing the sparsity pattern, all
     * entries are valid. However, before compression, the sparsity pattern
     * allocated some memory to be used while still adding new nonzero
     * entries; if you create iterators in this phase of the sparsity
     * pattern's lifetime, you will iterate over elements that are not valid.
     * If this is so, then this function will return false.
     */
    bool
    is_valid_entry() const;

    /**
     * Comparison. True, if both iterators point to the same matrix position.
     */
    bool
    operator==(const Accessor &) const;

    /**
     * Comparison operator. Result is true if either the first row number is
     * smaller or if the row numbers are equal and the first index is smaller.
     *
     * This function is only valid if both iterators point into the same
     * sparsity pattern.
     */
    bool
    operator<(const Accessor &) const;

  protected:
    DeclExceptionMsg(DummyAccessor,
                     "The instance of this class was initialized"
                     " without SparsityPattern object, which"
                     " means that it is a dummy accessor that can"
                     " not do any operations.");

    /**
     * The sparsity pattern we operate on accessed.
     */
    const SparsityPattern *container;

    /**
     * Index in global sparsity pattern. This index represents the location
     * the iterator/accessor points to within the array of the SparsityPattern
     * class that stores the column numbers. It is also the index within the
     * values array of a sparse matrix that stores the corresponding value of
     * this site.
     */
    std::size_t linear_index;

    /**
     * Move the accessor to the next nonzero entry in the matrix.
     */
    void
    advance();

    // Grant access to iterator class.
    friend class LinearIndexIterator<Iterator, Accessor>;

    // Grant access to accessor class of ChunkSparsityPattern.
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
  class Iterator : public LinearIndexIterator<Iterator, Accessor>
  {
  public:
    /**
     * Size type.
     */
    using size_type = types::global_dof_index;

    /**
     * Type of the stored pointer.
     */
    using container_pointer_type = SparsityPattern *;

    /**
     * Constructor. Create an iterator into the sparsity pattern @p sp for the
     * given global index (i.e., the index of the given element counting from
     * the zeroth row).
     */
    Iterator(const SparsityPattern *sp, const std::size_t linear_index);

    /**
     * Constructor. Create an iterator into the sparsity pattern @p sp for
     * a given accessor.
     */
    Iterator(const Accessor &accessor);
  };
} // namespace SparsityPatternIterators

/**
 * This class stores a sparsity pattern in
 * the <a
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
 * While this class forms the basis upon which SparseMatrix objects base
 * their storage format, and thus plays a central role in setting up linear
 * systems, it is rarely set up directly due to the way it stores its
 * information. Rather, one typically goes through an intermediate format
 * first, see for example the step-2 tutorial program as well as the
 * documentation topic
 * @ref Sparsity.
 *
 * You can iterate over entries in the pattern using begin(), end(),
 * begin(row), and end(row). These functions return an iterator of type
 * SparsityPatternIterators::Iterator. When dereferencing an iterator @p it,
 * you have access to the member functions in
 * SparsityPatternIterators::Accessor, like <tt>it->column()</tt> and
 * <tt>it->row()</tt>.
 */
class SparsityPattern : public SparsityPatternBase
{
public:
  /**
   * Declare type for container size.
   */
  using size_type = types::global_dof_index;

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
  static constexpr size_type invalid_entry = numbers::invalid_size_type;

  /**
   * Typedef an iterator class that allows to walk over all nonzero elements
   * of a sparsity pattern.
   */
  using const_iterator = SparsityPatternIterators::Iterator;

  /**
   * Typedef an iterator class that allows to walk over all nonzero elements
   * of a sparsity pattern.
   *
   * Since the iterator does not allow to modify the sparsity pattern, this
   * type is the same as that for @p const_iterator.
   */
  using iterator = SparsityPatternIterators::Iterator;

  /**
   * @name Iterators
   *
   * @{
   */

  /**
   * Iterator starting at the first entry of the matrix. The resulting
   * iterator can be used to walk over all nonzero entries of the sparsity
   * pattern.
   *
   * The order in which elements are accessed depends on the storage scheme
   * implemented by derived classes.
   */
  iterator
  begin() const;

  /**
   * Final iterator.
   */
  iterator
  end() const;

  /**
   * Iterator starting at the first entry of row <tt>r</tt>.
   *
   * Note that if the given row is empty, i.e. does not contain any nonzero
   * entries, then the iterator returned by this function equals
   * <tt>end(r)</tt>. Note also that the iterator may not be dereferenceable in
   * that case.
   *
   * The order in which elements are accessed depends on the storage scheme
   * implemented by derived classes.
   */
  iterator
  begin(const size_type r) const;

  /**
   * Final iterator of row <tt>r</tt>. It points to the first element past the
   * end of line @p r, or past the end of the entire sparsity pattern.
   *
   * Note that the end iterator is not necessarily dereferenceable. This is in
   * particular the case if it is the end iterator for the last row of a
   * matrix.
   */
  iterator
  end(const size_type r) const;

  /**
   * @}
   */

  /**
   * @name Construction and setup
   *
   * Constructors, destructor, functions initializing, copying and filling an
   * object.
   */
  /**
   * @{
   */

  /**
   * Initialize the matrix empty, that is with no memory allocated. This is
   * useful if you want such objects as member variables in other classes. You
   * can make the structure usable by calling the reinit() function.
   */
  SparsityPattern();

  /**
   * Copy constructor. This constructor is only allowed to be called if the
   * matrix structure to be copied is empty. This is so in order to prevent
   * involuntary copies of objects for temporaries, which can use large
   * amounts of computing time. However, copy constructors are needed if one
   * wants to place a SparsityPattern in a container, e.g., to write such
   * statements like <tt>v.push_back (SparsityPattern());</tt>, with
   * <tt>v</tt> a std::vector of SparsityPattern objects.
   *
   * Usually, it is sufficient to use the explicit keyword to disallow
   * unwanted temporaries, but this does not work for <tt>std::vector</tt>s.
   * Since copying a structure like this is not useful anyway because multiple
   * matrices can use the same sparsity structure, copies are only allowed for
   * empty objects, as described above.
   */
  SparsityPattern(const SparsityPattern &);

  /**
   * Initialize a rectangular pattern of size <tt>m x n</tt>.
   *
   * @param[in] m The number of rows.
   * @param[in] n The number of columns.
   * @param[in] max_per_row Maximum number of nonzero entries per row.
   */
  SparsityPattern(const size_type    m,
                  const size_type    n,
                  const unsigned int max_per_row);


  /**
   * Initialize a rectangular pattern of size <tt>m x n</tt>.
   *
   * @param[in] m The number of rows.
   * @param[in] n The number of columns.
   * @param[in] row_lengths Possible number of nonzero entries for each row.
   * This vector must have one entry for each row.
   */
  SparsityPattern(const size_type                  m,
                  const size_type                  n,
                  const std::vector<unsigned int> &row_lengths);

  /**
   * Initialize a quadratic pattern of dimension <tt>m</tt> with at most
   * <tt>max_per_row</tt> nonzero entries per row.
   *
   * This constructor automatically enables optimized storage of diagonal
   * elements. To avoid this, use the constructor taking row and column
   * numbers separately.
   */
  SparsityPattern(const size_type m, const unsigned int max_per_row);

  /**
   * Initialize a quadratic pattern of size <tt>m x m</tt>.
   *
   * @param[in] m The number of rows and columns.
   * @param[in] row_lengths Maximum number of nonzero entries for each row.
   * This vector must have one entry for each row.
   */
  SparsityPattern(const size_type                  m,
                  const std::vector<unsigned int> &row_lengths);

  /**
   * Make a copy with extra off-diagonals.
   *
   * This constructs objects intended for the application of the ILU(n)-method
   * or other incomplete decompositions.  Therefore, additional to the
   * original entry structure, space for <tt>extra_off_diagonals</tt>
   * side-diagonals is provided on both sides of the main diagonal.
   *
   * <tt>max_per_row</tt> is the maximum number of nonzero elements per row
   * which this structure is to hold. It is assumed that this number is
   * sufficiently large to accommodate both the elements in <tt>original</tt>
   * as well as the new off-diagonal elements created by this constructor. You
   * will usually want to give the same number as you gave for
   * <tt>original</tt> plus the number of side diagonals times two. You may
   * however give a larger value if you wish to add further nonzero entries
   * for the decomposition based on other criteria than their being on
   * side-diagonals.
   *
   * This function requires that <tt>original</tt> refers to a quadratic
   * matrix structure.  It must be compressed. The matrix structure is not
   * compressed after this function finishes.
   */
  SparsityPattern(const SparsityPattern &original,
                  const unsigned int     max_per_row,
                  const size_type        extra_off_diagonals);

  /**
   * Copy operator. For this the same holds as for the copy constructor: it is
   * declared, defined and fine to be called, but the latter only for empty
   * objects.
   */
  SparsityPattern &
  operator=(const SparsityPattern &);

  /**
   * Reallocate memory for a matrix of size @p m times @p n. The number of
   * entries for each row is taken from the ArrayView @p row_lengths which
   * has to give this number of each row $i=0\ldots m-1$.
   *
   * If <tt>m*n==0</tt> all memory is freed, resulting in a total
   * reinitialization of the object. If it is nonzero, new memory is only
   * allocated if the new size extends the old one. This is done to save time
   * and to avoid fragmentation of the heap.
   */
  void
  reinit(const size_type                      m,
         const size_type                      n,
         const ArrayView<const unsigned int> &row_lengths);

  /**
   * Same as the other reinit(), but uses a std::vector instead of an ArrayView.
   */
  void
  reinit(const size_type                  m,
         const size_type                  n,
         const std::vector<unsigned int> &row_lengths);

  /**
   * Reallocate memory and set up data structures for a new matrix with @p m
   * rows and @p n columns, with at most @p max_per_row
   * nonzero entries per row.
   *
   * This function simply maps its operations to the other reinit()
   * function.
   */
  void
  reinit(const size_type m, const size_type n, const unsigned int max_per_row);

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
  void
  compress();


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
  void
  copy_from(const size_type       n_rows,
            const size_type       n_cols,
            const ForwardIterator begin,
            const ForwardIterator end);

  /**
   * Copy data from a DynamicSparsityPattern. Previous content of this object
   * is lost, and the sparsity pattern is in compressed mode afterwards.
   */
  void
  copy_from(const DynamicSparsityPattern &dsp);

  /**
   * Copy data from a SparsityPattern. Previous content of this object is
   * lost, and the sparsity pattern is in compressed mode afterwards.
   */
  void
  copy_from(const SparsityPattern &sp);

  /**
   * Take a full matrix and use its nonzero entries to generate a sparse
   * matrix entry pattern for this object.
   *
   * Previous content of this object is lost, and the sparsity pattern is in
   * compressed mode afterwards.
   *
   * Once you have built a sparsity pattern with this function, you
   * probably want to attach a SparseMatrix object to it. The original
   * `matrix` object can then be copied into this SparseMatrix object
   * using the version of SparseMatrix::copy_from() that takes a
   * FullMatrix object as argument. Through this procedure, you can
   * convert a FullMatrix into a SparseMatrix.
   */
  template <typename number>
  void
  copy_from(const FullMatrix<number> &matrix);

  /**
   * Add a nonzero entry to the matrix. This function may only be called for
   * non-compressed sparsity patterns.
   *
   * If the entry already exists, nothing bad happens.
   */
  void
  add(const size_type i, const size_type j);

  /**
   * Add several nonzero entries to the specified matrix row.  This function
   * may only be called for non-compressed sparsity patterns.
   *
   * If some of the entries already exist, nothing bad happens.
   */
  template <typename ForwardIterator>
  void
  add_entries(const size_type row,
              ForwardIterator begin,
              ForwardIterator end,
              const bool      indices_are_sorted = false);

  virtual void
  add_row_entries(const size_type                  &row,
                  const ArrayView<const size_type> &columns,
                  const bool indices_are_sorted = false) override;

  using SparsityPatternBase::add_entries;

  /**
   * Make the sparsity pattern symmetric by adding the sparsity pattern of the
   * transpose object.
   *
   * This function throws an exception if the sparsity pattern does not
   * represent a quadratic matrix.
   */
  void
  symmetrize();

  /**
   * @}
   */


  /**
   * @name Querying information
   */
  /**
   * @{
   */

  /**
   * Return the number of nonzero elements of this matrix. Actually, it
   * returns the number of entries in the sparsity pattern; if any of the
   * entries should happen to be zero, it is counted anyway.
   *
   * This function may only be called if the matrix struct is compressed. It
   * does not make too much sense otherwise anyway.
   */
  std::size_t
  n_nonzero_elements() const;

  /**
   * Return whether the object is empty. It is empty if no memory is
   * allocated, which is the same as that both dimensions are zero.
   */
  bool
  empty() const;

  /**
   * Compute the bandwidth of the matrix represented by this structure. The
   * bandwidth is the maximum of $|i-j|$ for which the index pair $(i,j)$
   * represents a nonzero entry of the matrix. Consequently, the maximum
   * bandwidth a $n\times m$ matrix can have is $\max\{n-1,m-1\}$, a diagonal
   * matrix has bandwidth 0, and there are at most $2*q+1$ entries per row if
   * the bandwidth is $q$. The returned quantity is sometimes called "half
   * bandwidth" in the literature.
   */
  size_type
  bandwidth() const;

  /**
   * Return whether the structure is compressed or not.
   */
  bool
  is_compressed() const;

  /**
   * Return the maximum number of entries per row. Before compression, this
   * equals the number given to the constructor, while after compression, it
   * equals the maximum number of entries actually allocated by the user.
   */
  size_type
  max_entries_per_row() const;

  /**
   * Test for equality of two SparsityPatterns.
   */
  bool
  operator==(const SparsityPattern &) const;

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
  bool
  stores_only_added_elements() const;

  /**
   * Number of entries in a specific row.
   */
  unsigned int
  row_length(const size_type row) const;

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object. See MemoryConsumption.
   */
  std::size_t
  memory_consumption() const;

  /**
   * @}
   */

  /**
   * @name Accessing entries
   *
   * @{
   */

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
  size_type
  operator()(const size_type i, const size_type j) const;

  /**
   * Check if a value at a certain position may be non-zero.
   */
  bool
  exists(const size_type i, const size_type j) const;

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
  size_type
  column_number(const size_type row, const unsigned int index) const;

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
  matrix_position(const std::size_t global_index) const;

  /**
   * The index of a global matrix entry in its row.
   *
   * This function is analogous to operator(), but it computes the index not
   * with respect to the total field, but only with respect to the row
   * <tt>j</tt>.
   */
  size_type
  row_position(const size_type i, const size_type j) const;

  /**
   * @}
   */

  /**
   * @name Input/Output
   */
  /**
   * @{
   */

  /**
   * Print the sparsity of the matrix. The output consists of one line per row
   * of the format <tt>[i,j1,j2,j3,...]</tt>. <i>i</i> is the row number and
   * <i>jn</i> are the allocated columns in this row.
   */
  void
  print(std::ostream &out) const;

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
  void
  print_gnuplot(std::ostream &out) const;

  /**
   * Prints the sparsity of the matrix in a .svg file which can be opened in a
   * web browser. The .svg file contains squares which correspond to the
   * entries in the matrix. An entry in the matrix which contains a non-zero
   * value corresponds with a red square while a zero-valued entry in the
   * matrix correspond with a white square.
   */
  void
  print_svg(std::ostream &out) const;

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
  void
  block_write(std::ostream &out) const;

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
  void
  block_read(std::istream &in);

  /**
   * Write the data of this object to a stream for the purpose of
   * serialization
   */
  template <class Archive>
  void
  save(Archive &ar, const unsigned int version) const;

  /**
   * Read the data of this object from a stream for the purpose of
   * serialization
   */
  template <class Archive>
  void
  load(Archive &ar, const unsigned int version);

#ifdef DOXYGEN
  /**
   * Write and read the data of this object from a stream for the purpose
   * of serialization.
   */
  template <class Archive>
  void
  serialize(Archive &archive, const unsigned int version);
#else
  // This macro defines the serialize() method that is compatible with
  // the templated save() and load() method that have been implemented.
  BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif

  /**
   * @}
   */

  /**
   * @addtogroup Exceptions
   *
   * @{
   */
  /**
   * Exception
   */
  DeclException2(ExcIteratorRange,
                 int,
                 int,
                 << "The iterators denote a range of " << arg1
                 << " elements, but the given number of rows was " << arg2);
  /**
   * Exception
   */
  DeclException1(ExcInvalidNumberOfPartitions,
                 int,
                 << "The number of partitions you gave is " << arg1
                 << ", but must be greater than zero.");

  /**
   * You tried to add an element to a row, but there was no space left.
   */
  DeclException2(ExcNotEnoughSpace,
                 int,
                 int,
                 << "Upon entering a new entry to row " << arg1
                 << ": there was no free entry any more. " << std::endl
                 << "(Maximum number of entries for this row: " << arg2
                 << "; maybe the matrix is already compressed?)");

  /**
   * This operation changes the structure of the SparsityPattern and is not
   * possible after compress() has been called.
   */
  DeclExceptionMsg(
    ExcMatrixIsCompressed,
    "The operation you attempted changes the structure of the SparsityPattern "
    "and is not possible after compress() has been called.");

  /**
   * The operation is only allowed after the SparsityPattern has been set up
   * and compress() was called.
   */
  DeclExceptionMsg(
    ExcNotCompressed,
    "The operation you attempted is only allowed after the SparsityPattern "
    "has been set up and compress() was called.");

  /**
   * @}
   */
private:
  /**
   * Is special treatment of diagonals enabled?
   */
  bool store_diagonal_first_in_row;

  /**
   * Maximum number of rows that can be stored in the #rowstart array.  Since
   * reallocation of that array only happens if the present one is too small,
   * but never when the size of this matrix structure shrinks, #max_dim might
   * be larger than #rows and in this case #rowstart has more elements than
   * are used.
   */
  size_type max_dim;

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
  std::unique_ptr<std::size_t[]> rowstart;

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
  std::unique_ptr<size_type[]> colnums;

  /**
   * Store whether the compress() function was called for this object.
   */
  bool compressed;

  // Make all sparse matrices friends of this class.
  template <typename number>
  friend class SparseMatrix;
  template <typename number>
  friend class SparseLUDecomposition;
  template <typename number>
  friend class SparseILU;
  template <typename number>
  friend class ChunkSparseMatrix;

  friend class ChunkSparsityPattern;
  friend class DynamicSparsityPattern;

  // Also give access to internal details to the iterator/accessor classes.
  friend class SparsityPatternIterators::Iterator;
  friend class SparsityPatternIterators::Accessor;
  friend class ChunkSparsityPatternIterators::Accessor;
};


/**
 * @}
 */
/*---------------------- Inline functions -----------------------------------*/

#ifndef DOXYGEN


namespace SparsityPatternIterators
{
  inline Accessor::Accessor(const SparsityPattern *sparsity_pattern,
                            const std::size_t      i)
    : container(sparsity_pattern)
    , linear_index(i)
  {}



  inline Accessor::Accessor(const SparsityPattern *sparsity_pattern)
    : container(sparsity_pattern)
    , linear_index(container->rowstart[container->rows])
  {}



  inline Accessor::Accessor()
    : container(nullptr)
    , linear_index(numbers::invalid_size_type)
  {}



  inline bool
  Accessor::is_valid_entry() const
  {
    Assert(container != nullptr, DummyAccessor());
    return (linear_index < container->rowstart[container->rows] &&
            container->colnums[linear_index] != SparsityPattern::invalid_entry);
  }



  inline size_type
  Accessor::row() const
  {
    Assert(is_valid_entry() == true, ExcInvalidIterator());

    const std::size_t *insert_point =
      std::upper_bound(container->rowstart.get(),
                       container->rowstart.get() + container->rows + 1,
                       linear_index);
    return insert_point - container->rowstart.get() - 1;
  }



  inline size_type
  Accessor::column() const
  {
    Assert(is_valid_entry() == true, ExcInvalidIterator());

    return (container->colnums[linear_index]);
  }



  inline size_type
  Accessor::index() const
  {
    Assert(is_valid_entry() == true, ExcInvalidIterator());

    return linear_index - container->rowstart[row()];
  }



  inline size_type
  Accessor::global_index() const
  {
    Assert(is_valid_entry() == true, ExcInvalidIterator());

    return linear_index;
  }



  inline bool
  Accessor::operator==(const Accessor &other) const
  {
    return (container == other.container && linear_index == other.linear_index);
  }



  inline bool
  Accessor::operator<(const Accessor &other) const
  {
    Assert(container != nullptr, DummyAccessor());
    Assert(other.container != nullptr, DummyAccessor());
    Assert(container == other.container, ExcInternalError());

    return linear_index < other.linear_index;
  }



  inline void
  Accessor::advance()
  {
    Assert(container != nullptr, DummyAccessor());
    Assert(linear_index < container->rowstart[container->rows],
           ExcIteratorPastEnd());
    ++linear_index;
  }


  inline Iterator::Iterator(const SparsityPattern *sp,
                            const std::size_t      linear_index)
    : LinearIndexIterator<Iterator, Accessor>(Accessor(sp, linear_index))
  {}


  inline Iterator::Iterator(const Accessor &accessor)
    : LinearIndexIterator<Iterator, Accessor>(accessor)
  {}


} // namespace SparsityPatternIterators



inline std::size_t
SparsityPattern::n_nonzero_elements() const
{
  Assert(compressed, ExcNotCompressed());

  if ((rowstart != nullptr) && (colnums != nullptr))
    return rowstart[rows] - rowstart[0];
  else
    // the object is empty or has zero size
    return 0;
}



inline bool
SparsityPattern::is_compressed() const
{
  return compressed;
}



inline bool
SparsityPattern::stores_only_added_elements() const
{
  return (store_diagonal_first_in_row == false);
}



inline unsigned int
SparsityPattern::row_length(const size_type row) const
{
  AssertIndexRange(row, rows);
  return rowstart[row + 1] - rowstart[row];
}



inline SparsityPattern::size_type
SparsityPattern::column_number(const size_type    row,
                               const unsigned int index) const
{
  AssertIndexRange(row, rows);
  AssertIndexRange(index, row_length(row));

  return colnums[rowstart[row] + index];
}



inline SparsityPattern::iterator
SparsityPattern::begin() const
{
  if (n_rows() > 0)
    return {this, rowstart[0]};
  else
    return end();
}



inline SparsityPattern::iterator
SparsityPattern::end() const
{
  if (n_rows() > 0)
    return {this, rowstart[rows]};
  else
    return {nullptr, 0};
}



inline SparsityPattern::iterator
SparsityPattern::begin(const size_type r) const
{
  AssertIndexRange(r, n_rows());

  return {this, rowstart[r]};
}



inline SparsityPattern::iterator
SparsityPattern::end(const size_type r) const
{
  AssertIndexRange(r, n_rows());

  return {this, rowstart[r + 1]};
}



namespace internal
{
  namespace SparsityPatternTools
  {
    /**
     * Declare type for container size.
     */
    using size_type = types::global_dof_index;

    inline size_type
    get_column_index_from_iterator(const size_type i)
    {
      return i;
    }



    template <typename value>
    inline size_type
    get_column_index_from_iterator(const std::pair<size_type, value> &i)
    {
      return i.first;
    }



    template <typename value>
    inline size_type
    get_column_index_from_iterator(const std::pair<const size_type, value> &i)
    {
      return i.first;
    }
  } // namespace SparsityPatternTools
} // namespace internal



template <typename ForwardIterator>
void
SparsityPattern::copy_from(const size_type       n_rows,
                           const size_type       n_cols,
                           const ForwardIterator begin,
                           const ForwardIterator end)
{
  Assert(static_cast<size_type>(std::distance(begin, end)) == n_rows,
         ExcIteratorRange(std::distance(begin, end), n_rows));

  // first determine row lengths for each row. if the matrix is quadratic,
  // then we might have to add an additional entry for the diagonal, if that
  // is not yet present. as we have to call compress anyway later on, don't
  // bother to check whether that diagonal entry is in a certain row or not
  const bool                is_square = (n_rows == n_cols);
  std::vector<unsigned int> row_lengths;
  row_lengths.reserve(n_rows);
  for (ForwardIterator i = begin; i != end; ++i)
    row_lengths.push_back(std::distance(i->begin(), i->end()) +
                          (is_square ? 1 : 0));
  reinit(n_rows, n_cols, row_lengths);

  // now enter all the elements into the matrix. note that if the matrix is
  // quadratic, then we already have the diagonal element preallocated
  //
  // for use in the inner loop, we define an alias to the type of the inner
  // iterators
  size_type row = 0;
  using inner_iterator =
    typename std::iterator_traits<ForwardIterator>::value_type::const_iterator;
  for (ForwardIterator i = begin; i != end; ++i, ++row)
    {
      size_type           *cols = &colnums[rowstart[row]] + (is_square ? 1 : 0);
      const inner_iterator end_of_row = i->end();
      for (inner_iterator j = i->begin(); j != end_of_row; ++j)
        {
          const size_type col =
            internal::SparsityPatternTools::get_column_index_from_iterator(*j);
          AssertIndexRange(col, n_cols);

          if ((col != row) || !is_square)
            *cols++ = col;
        }
    }

  // finally compress everything. this also sorts the entries within each row
  compress();
}



template <class Archive>
inline void
SparsityPattern::save(Archive &ar, const unsigned int) const
{
  // forward to serialization function in the base class.
  ar &boost::serialization::base_object<const EnableObserverPointer>(*this);

  ar &max_dim &rows &cols &max_vec_len &max_row_length &compressed;

  if (max_dim != 0)
    ar &boost::serialization::make_array(rowstart.get(), max_dim + 1);
  else
    Assert(rowstart.get() == nullptr, ExcInternalError());

  if (max_vec_len != 0)
    ar &boost::serialization::make_array(colnums.get(), max_vec_len);
  else
    Assert(colnums.get() == nullptr, ExcInternalError());
  ar &store_diagonal_first_in_row;
}



template <class Archive>
inline void
SparsityPattern::load(Archive &ar, const unsigned int)
{
  // forward to serialization function in the base class.
  ar &boost::serialization::base_object<EnableObserverPointer>(*this);

  ar &max_dim &rows &cols &max_vec_len &max_row_length &compressed;

  if (max_dim != 0)
    {
      rowstart = std::make_unique<std::size_t[]>(max_dim + 1);
      ar &boost::serialization::make_array(rowstart.get(), max_dim + 1);
    }
  else
    rowstart.reset();

  if (max_vec_len != 0)
    {
      colnums = std::make_unique<size_type[]>(max_vec_len);
      ar &boost::serialization::make_array(colnums.get(), max_vec_len);
    }
  else
    colnums.reset();
  ar &store_diagonal_first_in_row;
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
