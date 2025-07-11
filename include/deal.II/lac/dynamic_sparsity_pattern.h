// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_dynamic_sparsity_pattern_h
#define dealii_dynamic_sparsity_pattern_h


#include <deal.II/base/config.h>

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/sparsity_pattern_base.h>

#include <algorithm>
#include <iostream>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
class DynamicSparsityPattern;
#endif

/**
 * @addtogroup Sparsity
 * @{
 */


/**
 * Iterators on objects of type DynamicSparsityPattern.
 */
namespace DynamicSparsityPatternIterators
{
  // forward declaration
  class Iterator;

  /**
   * Declare type for container size.
   */
  using size_type = types::global_dof_index;

  /**
   * Accessor class for iterators into objects of type DynamicSparsityPattern.
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
     * Constructor.
     */
    Accessor(const DynamicSparsityPattern *sparsity_pattern,
             const size_type               row,
             const unsigned int            index_within_row);

    /**
     * Constructor. Construct the end accessor for the given sparsity pattern.
     */
    Accessor(const DynamicSparsityPattern *sparsity_pattern);

    /**
     * Default constructor creating a dummy accessor. This constructor is here
     * only to be able to store accessors in STL containers such as
     * `std::vector`.
     */
    Accessor();

    /**
     * Row number of the element represented by this object.
     */
    size_type
    row() const;

    /**
     * Index within the current row of the element represented by this object.
     */
    size_type
    index() const;

    /**
     * Column number of the element represented by this object.
     */
    size_type
    column() const;

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
                     " without DynamicSparsityPattern object, which"
                     " means that it is a dummy accessor that can"
                     " not do any operations.");

    /**
     * The sparsity pattern we operate on accessed.
     */
    const DynamicSparsityPattern *sparsity_pattern;

    /**
     * The row we currently point into.
     */
    size_type current_row;

    /**
     * A pointer to the element within the current row that we currently point
     * to.
     */
    std::vector<size_type>::const_iterator current_entry;

    /**
     * A pointer to the end of the current row. We store this to make
     * comparison against the end of line iterator cheaper as it otherwise
     * needs to do the IndexSet translation from row index to the index within
     * the 'lines' array of DynamicSparsityPattern.
     */
    std::vector<size_type>::const_iterator end_of_row;

    /**
     * Move the accessor to the next nonzero entry in the matrix.
     */
    void
    advance();

    // Grant access to iterator class.
    friend class Iterator;
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
   * DynamicSparsityPattern class. As a consequence, some operations are cheap
   * and some are not. In particular, it is cheap to access the column index
   * of the sparsity pattern entry pointed to. On the other hand, it is
   * expensive to compute the distance between two iterators. As a
   * consequence, when you design algorithms that use these iterators, it is
   * common practice to not loop over <i>all</i> elements of a sparsity
   * pattern at once, but to have an outer loop over all rows and within this
   * loop iterate over the elements of this row. This way, you only ever need
   * to dereference the iterator to obtain the column indices whereas the
   * (expensive) lookup of the row index can be avoided by using the loop
   * index instead.
   */
  class Iterator
  {
  public:
    /**
     * Constructor. Create an iterator into the sparsity pattern @p sp for the
     * given global index (i.e., the index of the given element counting from
     * the zeroth row).
     */
    Iterator(const DynamicSparsityPattern *sp,
             const size_type               row,
             const unsigned int            index_within_row);

    /**
     * Constructor. Create an invalid (end) iterator into the sparsity pattern
     * @p sp.
     */
    Iterator(const DynamicSparsityPattern *sp);

    /**
     * Default constructor creating an invalid iterator. This constructor is
     * here only to be able to store iterators in STL containers such as
     * `std::vector`.
     */
    Iterator() = default;

    /**
     * Prefix increment.
     */
    Iterator &
    operator++();

    /**
     * Postfix increment.
     */
    Iterator
    operator++(int);

    /**
     * Dereferencing operator.
     */
    const Accessor &
    operator*() const;

    /**
     * Dereferencing operator.
     */
    const Accessor *
    operator->() const;

    /**
     * Comparison. True, if both iterators point to the same matrix position.
     */
    bool
    operator==(const Iterator &) const;

    /**
     * Inverse of <tt>==</tt>.
     */
    bool
    operator!=(const Iterator &) const;

    /**
     * Comparison operator. Result is true if either the first row number is
     * smaller or if the row numbers are equal and the first index is smaller.
     *
     * This function is only valid if both iterators point into the same
     * matrix.
     */
    bool
    operator<(const Iterator &) const;

    /**
     * Return the distance between the current iterator and the argument. The
     * distance is given by how many times one has to apply operator++ to the
     * current iterator to get the argument (for a positive return value), or
     * operator-- (for a negative return value).
     */
    int
    operator-(const Iterator &p) const;

  private:
    /**
     * Store an object of the accessor class.
     */
    Accessor accessor;
  };
} // namespace DynamicSparsityPatternIterators


/**
 * This class acts as an intermediate form of the SparsityPattern class. From
 * the interface it mostly represents a SparsityPattern object that is kept
 * compressed at all times. However, since the final sparsity pattern is not
 * known while constructing it, keeping the pattern compressed at all times
 * can only be achieved at the expense of either increased memory or run time
 * consumption upon use. The main purpose of this class is to avoid some
 * memory bottlenecks, so we chose to implement it memory conservative. The
 * chosen data format is too unsuited to be used for actual matrices, though.
 * It is therefore necessary to first copy the data of this object over to an
 * object of type SparsityPattern before using it in actual matrices.
 *
 * Another viewpoint is that this class does not need up front allocation of a
 * certain amount of memory, but grows as necessary.  An extensive description
 * of sparsity patterns can be found in the documentation of the
 * @ref Sparsity
 * topic.
 *
 * This class is an example of the "dynamic" type of
 * @ref Sparsity.
 * It is used in most tutorial programs in one way or another.
 *
 * <h3>Interface</h3>
 *
 * Since this class is intended as an intermediate replacement of the
 * SparsityPattern class, it has mostly the same interface, with small changes
 * where necessary. In particular, the add() function, and the functions
 * inquiring properties of the sparsity pattern are the same.
 *
 *
 * <h3>Usage</h3>
 *
 * Usage of this class is explained in step-2 (without constraints) and step-6
 * (with AffineConstraints) and typically looks as follows:
 * @code
 * DynamicSparsityPattern dynamic_pattern (dof_handler.n_dofs());
 * DoFTools::make_sparsity_pattern (dof_handler,
 *                                  dynamic_pattern,
 *                                  constraints);
 * SparsityPattern sp;
 * sp.copy_from (dynamic_pattern);
 * @endcode
 */
class DynamicSparsityPattern : public SparsityPatternBase
{
public:
  /**
   * Declare the type for container size.
   */
  using size_type = types::global_dof_index;

  /**
   * Typedef an for iterator class that allows to walk over all nonzero
   * elements of a sparsity pattern.
   *
   * Since the iterator does not allow to modify the sparsity pattern, this
   * type is the same as that for @p const_iterator.
   */
  using iterator = DynamicSparsityPatternIterators::Iterator;

  /**
   * Typedef for an iterator class that allows to walk over all nonzero
   * elements of a sparsity pattern.
   */
  using const_iterator = DynamicSparsityPatternIterators::Iterator;

  /**
   * Initialize as an empty object. This is useful if you want such objects as
   * member variables in other classes. You can make the structure usable by
   * calling the reinit() function.
   */
  DynamicSparsityPattern();

  /**
   * Copy constructor. This constructor is only allowed to be called if the
   * sparsity structure to be copied is empty. This is so in order to prevent
   * involuntary copies of objects for temporaries, which can use large
   * amounts of computing time.  However, copy constructors are needed if you
   * want to place a DynamicSparsityPattern in a container, e.g. to write such
   * statements like <tt>v.push_back (DynamicSparsityPattern());</tt>, with @p
   * v a vector of @p DynamicSparsityPattern objects.
   */
  DynamicSparsityPattern(const DynamicSparsityPattern &);

  /**
   * Initialize a rectangular sparsity pattern with @p m rows and @p n
   * columns. The @p rowset restricts the storage to elements in rows of this
   * set.  Adding elements outside of this set has no effect. The default
   * argument keeps all entries.
   */
  DynamicSparsityPattern(const size_type m,
                         const size_type n,
                         const IndexSet &rowset = IndexSet());

  /**
   * Create a square SparsityPattern using the given index set. The total size
   * is given by the size of @p indexset and only rows corresponding to
   * indices in @p indexset are stored on the current processor.
   */
  DynamicSparsityPattern(const IndexSet &indexset);

  /**
   * Initialize a square pattern of dimension @p n.
   */
  DynamicSparsityPattern(const size_type n);

  /**
   * Copy operator. For this the same holds as for the copy constructor: it is
   * declared, defined and fine to be called, but the latter only for empty
   * objects.
   */
  DynamicSparsityPattern &
  operator=(const DynamicSparsityPattern &);

  /**
   * Reallocate memory and set up data structures for a new sparsity pattern
   * with @p m rows and @p n columns. The @p rowset restricts the storage to
   * elements in rows of this set.  Adding elements outside of this set has no
   * effect. The default argument keeps all entries.
   */
  void
  reinit(const size_type m,
         const size_type n,
         const IndexSet &rowset = IndexSet());

  /**
   * Since this object is kept compressed at all times anyway, this function
   * does nothing, but is declared to make the interface of this class as much
   * alike as that of the SparsityPattern class.
   */
  void
  compress();

  /**
   * Return whether the object is empty. It is empty if no memory is
   * allocated, which is the same as that both dimensions are zero.
   */
  bool
  empty() const;

  /**
   * Return the maximum number of entries per row. Note that this number may
   * change as entries are added.
   */
  size_type
  max_entries_per_row() const;

  /**
   * Add a nonzero entry. If the entry already exists, this call does nothing.
   */
  void
  add(const size_type i, const size_type j);

  /**
   * Add several nonzero entries to the specified row. Already existing
   * entries are ignored.
   */
  template <typename ForwardIterator>
  void
  add_entries(const size_type row,
              ForwardIterator begin,
              ForwardIterator end,
              const bool      indices_are_unique_and_sorted = false);

  virtual void
  add_row_entries(const size_type                  &row,
                  const ArrayView<const size_type> &columns,
                  const bool indices_are_sorted = false) override;

  using SparsityPatternBase::add_entries;

  /**
   * Check if a value at a certain position may be non-zero.
   */
  bool
  exists(const size_type i, const size_type j) const;

  /**
   * Return a view of this sparsity pattern.
   * That is, for all rows in @p rows extract non-empty columns.
   * The resulting sparsity pattern will have number of rows equal
   * `rows.n_elements()`.
   */
  DynamicSparsityPattern
  get_view(const IndexSet &rows) const;

  /**
   * Make the sparsity pattern symmetric by adding the sparsity pattern of the
   * transpose object.
   *
   * This function throws an exception if the sparsity pattern does not
   * represent a square matrix.
   */
  void
  symmetrize();

  /**
   * Construct and store in this object the sparsity pattern corresponding to
   * the product of @p left and @p right sparsity pattern.
   */
  template <typename SparsityPatternTypeLeft, typename SparsityPatternTypeRight>
  void
  compute_mmult_pattern(const SparsityPatternTypeLeft  &left,
                        const SparsityPatternTypeRight &right);

  /**
   * Construct and store in this object the sparsity pattern corresponding to
   * the product of transposed @p left and non-transpose @p right sparsity pattern.
   */
  template <typename SparsityPatternTypeLeft, typename SparsityPatternTypeRight>
  void
  compute_Tmmult_pattern(const SparsityPatternTypeLeft  &left,
                         const SparsityPatternTypeRight &right);

  /**
   * Print the sparsity pattern. The output consists of one line per row of
   * the format <tt>[i,j1,j2,j3,...]</tt>. <i>i</i> is the row number and
   * <i>jn</i> are the allocated columns in this row.
   */
  void
  print(std::ostream &out) const;

  /**
   * Print the sparsity pattern in a format that @p gnuplot understands and
   * which can be used to plot the sparsity pattern in a graphical way. The
   * format consists of pairs <tt>i j</tt> of nonzero elements, each
   * representing one entry, one per line of the output file. Indices are
   * counted from zero on, as usual. Since sparsity patterns are printed in
   * the same way as matrices are displayed, we print the negative of the
   * column index, which means that the <tt>(0,0)</tt> element is in the top
   * left rather than in the bottom left corner.
   *
   * Print the sparsity pattern in gnuplot by setting the data style to dots
   * or points and use the @p plot command.
   */
  void
  print_gnuplot(std::ostream &out) const;

  /**
   * Number of entries in a specific row. This function can only be called if
   * the given row is a member of the index set of rows that we want to store.
   */
  size_type
  row_length(const size_type row) const;

  /**
   * Clear all entries stored in a specific row.
   */
  void
  clear_row(const size_type row);

  /**
   * Access to column number field.  Return the column number of the @p
   * indexth entry in @p row.
   */
  size_type
  column_number(const size_type row, const size_type index) const;

  /**
   * Return index of column @p col in row @p row. If the column does not
   * exist in this sparsity pattern, the returned value will be
   * 'numbers::invalid_size_type'.
   */
  size_type
  column_index(const size_type row, const size_type col) const;

  /**
   * @name Iterators
   * @{
   */

  /**
   * Iterator starting at the first entry of the matrix. The resulting
   * iterator can be used to walk over all nonzero entries of the sparsity
   * pattern.
   *
   * Note the discussion in the general documentation of this class about the
   * order in which elements are accessed.
   *
   * @note If the sparsity pattern has been initialized with an IndexSet that
   * denotes which rows to store, then iterators will simply skip over rows
   * that are not stored. In other words, they will look like empty rows, but
   * no exception will be generated when iterating over such rows.
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
   * Note also the discussion in the general documentation of this class about
   * the order in which elements are accessed.
   *
   * @note If the sparsity pattern has been initialized with an IndexSet that
   * denotes which rows to store, then iterators will simply skip over rows
   * that are not stored. In other words, they will look like empty rows, but
   * no exception will be generated when iterating over such rows.
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
   * Compute the bandwidth of the matrix represented by this structure. The
   * bandwidth is the maximum of $|i-j|$ for which the index pair $(i,j)$
   * represents a nonzero entry of the matrix.
   */
  size_type
  bandwidth() const;

  /**
   * Return the number of nonzero elements allocated through this sparsity
   * pattern.
   */
  size_type
  n_nonzero_elements() const;

  /**
   * Return the IndexSet that sets which rows are active on the current
   * processor. It corresponds to the IndexSet given to this class in the
   * constructor or in the reinit function.
   */
  const IndexSet &
  row_index_set() const;

  /**
   * Return the IndexSet that contains entries for all columns in which at least
   * one element exists in this sparsity pattern.
   *
   * @note In a parallel context, this only considers the locally stored rows.
   */
  IndexSet
  nonempty_cols() const;

  /**
   * Return the IndexSet that contains entries for all rows in which at least
   * one element exists in this sparsity pattern.
   *
   * @note In a parallel context, this only considers the locally stored rows.
   */
  IndexSet
  nonempty_rows() const;

  /**
   * return whether this object stores only those entries that have been added
   * explicitly, or if the sparsity pattern contains elements that have been
   * added through other means (implicitly) while building it. For the current
   * class, the result is always true.
   *
   * This function mainly serves the purpose of describing the current class
   * in cases where several kinds of sparsity patterns can be passed as
   * template arguments.
   */
  static bool
  stores_only_added_elements();

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  size_type
  memory_consumption() const;

private:
  /**
   * A flag that stores whether any entries have been added so far.
   */
  bool have_entries;

  /**
   * A set that contains the valid rows.
   */

  IndexSet rowset;


  /**
   * Store some data for each row describing which entries of this row are
   * nonzero. Data is stored sorted in the @p entries std::vector.  The vector
   * per row is dynamically growing upon insertion doubling its memory each
   * time.
   */
  struct Line
  {
  public:
    /**
     * Storage for the column indices of this row. This array is always kept
     * sorted.
     */
    std::vector<size_type> entries;

    /**
     * Add the given column number to this line.
     */
    void
    add(const size_type col_num);

    /**
     * Add the columns specified by the iterator range to this line.
     */
    template <typename ForwardIterator>
    void
    add_entries(ForwardIterator begin,
                ForwardIterator end,
                const bool      indices_are_sorted);

    /**
     * estimates memory consumption.
     */
    size_type
    memory_consumption() const;
  };


  /**
   * Actual data: store for each row the set of nonzero entries.
   */
  std::vector<Line> lines;

  // make the accessor class a friend
  friend class DynamicSparsityPatternIterators::Accessor;
};

/** @} */
/*---------------------- Inline functions -----------------------------------*/


namespace DynamicSparsityPatternIterators
{
  inline Accessor::Accessor(const DynamicSparsityPattern *sparsity_pattern,
                            const size_type               row,
                            const unsigned int            index_within_row)
    : sparsity_pattern(sparsity_pattern)
    , current_row(row)
    , current_entry(
        ((sparsity_pattern->rowset.size() == 0) ?
           sparsity_pattern->lines[current_row].entries.begin() :
           sparsity_pattern
             ->lines[sparsity_pattern->rowset.index_within_set(current_row)]
             .entries.begin()) +
        index_within_row)
    , end_of_row(
        (sparsity_pattern->rowset.size() == 0) ?
          sparsity_pattern->lines[current_row].entries.end() :
          sparsity_pattern
            ->lines[sparsity_pattern->rowset.index_within_set(current_row)]
            .entries.end())
  {
    AssertIndexRange(current_row, sparsity_pattern->n_rows());
    Assert((sparsity_pattern->rowset.size() == 0) ||
             sparsity_pattern->rowset.is_element(current_row),
           ExcMessage("You can't create an iterator into a "
                      "DynamicSparsityPattern's row that is not "
                      "actually stored by that sparsity pattern "
                      "based on the IndexSet argument to it."));
    AssertIndexRange(
      index_within_row,
      ((sparsity_pattern->rowset.size() == 0) ?
         sparsity_pattern->lines[current_row].entries.size() :
         sparsity_pattern
           ->lines[sparsity_pattern->rowset.index_within_set(current_row)]
           .entries.size()));
  }


  inline Accessor::Accessor(const DynamicSparsityPattern *sparsity_pattern)
    : sparsity_pattern(sparsity_pattern)
    , current_row(numbers::invalid_size_type)
    , current_entry()
    , end_of_row()
  {}



  inline Accessor::Accessor()
    : sparsity_pattern(nullptr)
    , current_row(numbers::invalid_size_type)
    , current_entry()
    , end_of_row()
  {}


  inline size_type
  Accessor::row() const
  {
    Assert(sparsity_pattern != nullptr, DummyAccessor());
    Assert(current_row < sparsity_pattern->n_rows(), ExcInternalError());

    return current_row;
  }


  inline size_type
  Accessor::column() const
  {
    Assert(sparsity_pattern != nullptr, DummyAccessor());
    Assert(current_row < sparsity_pattern->n_rows(), ExcInternalError());

    return *current_entry;
  }


  inline size_type
  Accessor::index() const
  {
    Assert(sparsity_pattern != nullptr, DummyAccessor());
    Assert(current_row < sparsity_pattern->n_rows(), ExcInternalError());

    return (current_entry -
            ((sparsity_pattern->rowset.size() == 0) ?
               sparsity_pattern->lines[current_row].entries.begin() :
               sparsity_pattern
                 ->lines[sparsity_pattern->rowset.index_within_set(current_row)]
                 .entries.begin()));
  }



  inline bool
  Accessor::operator==(const Accessor &other) const
  {
    Assert(sparsity_pattern != nullptr, DummyAccessor());
    Assert(other.sparsity_pattern != nullptr, DummyAccessor());
    // compare the sparsity pattern the iterator points into, the
    // current row, and the location within this row. ignore the
    // latter if the row is past-the-end because in that case the
    // current_entry field may not point to a deterministic location
    return (sparsity_pattern == other.sparsity_pattern &&
            current_row == other.current_row &&
            ((current_row == numbers::invalid_size_type) ||
             (current_entry == other.current_entry)));
  }



  inline bool
  Accessor::operator<(const Accessor &other) const
  {
    Assert(sparsity_pattern != nullptr, DummyAccessor());
    Assert(other.sparsity_pattern != nullptr, DummyAccessor());
    Assert(sparsity_pattern == other.sparsity_pattern, ExcInternalError());

    // if *this is past-the-end, then it is less than no one
    if (current_row == numbers::invalid_size_type)
      return (false);
    // now *this should be an valid value
    Assert(current_row < sparsity_pattern->n_rows(), ExcInternalError());

    // if other is past-the-end
    if (other.current_row == numbers::invalid_size_type)
      return (true);
    // now other should be an valid value
    Assert(other.current_row < sparsity_pattern->n_rows(), ExcInternalError());

    // both iterators are not one-past-the-end
    return ((current_row < other.current_row) ||
            ((current_row == other.current_row) &&
             (current_entry < other.current_entry)));
  }


  inline void
  Accessor::advance()
  {
    Assert(sparsity_pattern != nullptr, DummyAccessor());
    Assert(current_row < sparsity_pattern->n_rows(), ExcInternalError());

    // move to the next element in this row
    ++current_entry;

    // if this moves us beyond the end of the row, go to the next row
    // if possible, or set the iterator to an invalid state if not.
    //
    // going to the next row is a bit complicated because we may have
    // to skip over empty rows, and because we also have to avoid rows
    // that aren't listed in a possibly passed IndexSet argument of
    // the sparsity pattern. consequently, rather than trying to
    // duplicate code here, just call the begin() function of the
    // sparsity pattern itself
    if (current_entry == end_of_row)
      {
        if (current_row + 1 < sparsity_pattern->n_rows())
          *this = *sparsity_pattern->begin(current_row + 1);
        else
          *this = Accessor(sparsity_pattern); // invalid object
      }
  }



  inline Iterator::Iterator(const DynamicSparsityPattern *sparsity_pattern,
                            const size_type               row,
                            const unsigned int            index_within_row)
    : accessor(sparsity_pattern, row, index_within_row)
  {}



  inline Iterator::Iterator(const DynamicSparsityPattern *sparsity_pattern)
    : accessor(sparsity_pattern)
  {}



  inline Iterator &
  Iterator::operator++()
  {
    accessor.advance();
    return *this;
  }



  inline Iterator
  Iterator::operator++(int)
  {
    const Iterator iter = *this;
    accessor.advance();
    return iter;
  }



  inline const Accessor &
  Iterator::operator*() const
  {
    return accessor;
  }



  inline const Accessor *
  Iterator::operator->() const
  {
    return &accessor;
  }


  inline bool
  Iterator::operator==(const Iterator &other) const
  {
    return (accessor == other.accessor);
  }



  inline bool
  Iterator::operator!=(const Iterator &other) const
  {
    return !(*this == other);
  }


  inline bool
  Iterator::operator<(const Iterator &other) const
  {
    return accessor < other.accessor;
  }


  inline int
  Iterator::operator-(const Iterator &other) const
  {
    Assert(accessor.sparsity_pattern == other.accessor.sparsity_pattern,
           ExcInternalError());
    DEAL_II_NOT_IMPLEMENTED();

    return 0;
  }
} // namespace DynamicSparsityPatternIterators


inline void
DynamicSparsityPattern::Line::add(const size_type j)
{
  // first check the last element (or if line is still empty)
  if ((entries.empty()) || (entries.back() < j))
    {
      entries.push_back(j);
      return;
    }

  // do a binary search to find the place where to insert:
  std::vector<size_type>::iterator it =
    Utilities::lower_bound(entries.begin(), entries.end(), j);

  // If this entry is a duplicate, exit immediately
  if (*it == j)
    return;

  // Insert at the right place in the vector. Vector grows automatically to
  // fit elements. Always doubles its size.
  entries.insert(it, j);
}



inline void
DynamicSparsityPattern::add(const size_type i, const size_type j)
{
  AssertIndexRange(i, n_rows());
  AssertIndexRange(j, n_cols());

  if (rowset.size() > 0 && !rowset.is_element(i))
    return;

  have_entries = true;

  const size_type rowindex =
    rowset.size() == 0 ? i : rowset.index_within_set(i);
  lines[rowindex].add(j);
}



template <typename ForwardIterator>
inline void
DynamicSparsityPattern::add_entries(const size_type row,
                                    ForwardIterator begin,
                                    ForwardIterator end,
                                    const bool      indices_are_sorted)
{
  AssertIndexRange(row, rows);

  if (rowset.size() > 0 && !rowset.is_element(row))
    return;

  if (!have_entries && begin < end)
    have_entries = true;

  const size_type rowindex =
    rowset.size() == 0 ? row : rowset.index_within_set(row);
  lines[rowindex].add_entries(begin, end, indices_are_sorted);
}



inline types::global_dof_index
DynamicSparsityPattern::row_length(const size_type row) const
{
  AssertIndexRange(row, n_rows());

  if (!have_entries)
    return 0;

  if (rowset.size() > 0 && !rowset.is_element(row))
    return 0;

  const size_type rowindex =
    rowset.size() == 0 ? row : rowset.index_within_set(row);
  return lines[rowindex].entries.size();
}



inline types::global_dof_index
DynamicSparsityPattern::column_number(const size_type row,
                                      const size_type index) const
{
  AssertIndexRange(row, n_rows());
  Assert(rowset.size() == 0 || rowset.is_element(row), ExcInternalError());

  const size_type local_row =
    rowset.size() != 0u ? rowset.index_within_set(row) : row;
  AssertIndexRange(index, lines[local_row].entries.size());
  return lines[local_row].entries[index];
}



inline DynamicSparsityPattern::iterator
DynamicSparsityPattern::begin() const
{
  if (n_rows() > 0)
    return begin(0);
  else
    return end();
}


inline DynamicSparsityPattern::iterator
DynamicSparsityPattern::end() const
{
  return {this};
}



inline DynamicSparsityPattern::iterator
DynamicSparsityPattern::begin(const size_type r) const
{
  AssertIndexRange(r, n_rows());

  if (!have_entries)
    return {this};

  if (rowset.size() > 0)
    {
      // We have an IndexSet that describes the locally owned set. For
      // performance reasons we need to make sure that we don't do a
      // linear search over 0..n_rows(). Instead, find the first entry
      // >= row r in the locally owned set (this is done in log
      // n_ranges time inside at()). From there, we move forward until
      // we find a non-empty row. By iterating over the IndexSet instead
      // of incrementing the row index, we potentially skip over entries
      // not in the rowset.
      IndexSet::ElementIterator it = rowset.at(r);
      if (it == rowset.end())
        return end(); // we don't own any row between r and the end

      // Instead of using row_length(*it)==0 in the while loop below,
      // which involves an expensive index_within_set() call, we
      // look at the lines vector directly. This works, because we are
      // walking over this vector entry by entry anyways.
      size_type rowindex = rowset.index_within_set(*it);

      while (it != rowset.end() && lines[rowindex].entries.empty())
        {
          ++it;
          ++rowindex;
        }

      if (it == rowset.end())
        return end();
      else
        return {this, *it, 0};
    }

  // Without an index set we have to do a linear search starting at
  // row r until we find a non-empty one. We will check the lines vector
  // directly instead of going through the slower row_length() function
  size_type row = r;

  while (row < n_rows() && lines[row].entries.empty())
    {
      ++row;
    }

  if (row == n_rows())
    return {this};
  else
    return {this, row, 0};
}



inline DynamicSparsityPattern::iterator
DynamicSparsityPattern::end(const size_type r) const
{
  AssertIndexRange(r, n_rows());

  const size_type row = r + 1;
  if (row == n_rows())
    return {this};
  else
    return begin(row);
}



inline const IndexSet &
DynamicSparsityPattern::row_index_set() const
{
  return rowset;
}



inline bool
DynamicSparsityPattern::stores_only_added_elements()
{
  return true;
}


DEAL_II_NAMESPACE_CLOSE

#endif
