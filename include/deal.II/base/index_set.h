// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_index_set_h
#define dealii_index_set_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>

#include <boost/serialization/vector.hpp>

#include <algorithm>
#include <iterator>
#include <vector>


#ifdef DEAL_II_WITH_TRILINOS
#  include <Epetra_Map.h>
#  ifdef DEAL_II_TRILINOS_WITH_TPETRA
#    include <Tpetra_Map.hpp>
#  endif
#endif

#if defined(DEAL_II_WITH_MPI) || defined(DEAL_II_WITH_PETSC)
#  include <mpi.h>
#else
using MPI_Comm = int;
#  ifndef MPI_COMM_WORLD
#    define MPI_COMM_WORLD 0
#  endif
#endif

DEAL_II_NAMESPACE_OPEN

/**
 * A class that represents a subset of indices among a larger set. For
 * example, it can be used to denote the set of degrees of freedom within the
 * range $[0,\text{dof\_handler.n\_dofs})$ that belongs to a particular
 * subdomain, or those among all degrees of freedom that are stored on a
 * particular processor in a distributed parallel computation.
 *
 * This class can represent a collection of half-open ranges of indices as
 * well as individual elements. For practical purposes it also stores the
 * overall range these indices can assume. In other words, you need to specify
 * the size of the index space $[0,\text{size})$ of which objects of this
 * class are a subset.
 *
 * There are two ways to iterate over the IndexSets: First, begin() and end()
 * allow iteration over individual indices in the set. Second,
 * begin_interval() and end_interval() allow iteration over the half-open
 * ranges as described above.
 *
 * The data structures used in this class along with a rationale can be found
 * in the
 * @ref distributed_paper "Distributed Computing paper".
 *
 * @author Wolfgang Bangerth, 2009
 */
class IndexSet
{
public:
  // forward declarations:
  class ElementIterator;
  class IntervalIterator;

  /**
   * @p size_type is the type used for storing the size and the individual
   * entries in the IndexSet.
   */
  using size_type = types::global_dof_index;

  /**
   * One can see an IndexSet as a container of size size(), where the elements
   * of the containers are bool values that are either false or true,
   * depending on whether a particular index is an element of the IndexSet or
   * not. In other words, an IndexSet is a bit like a vector in which the
   * elements we store are booleans. In this view, the correct local alias
   * indicating the type of the elements of the vector would then be @p bool.
   *
   * On the other hand, @p bool has the disadvantage that it is not a
   * numerical type that, for example, allows multiplication with a @p double.
   * In other words, one can not easily use a vector of booleans in a place
   * where other vectors are allowed. Consequently, we declare the type of the
   * elements of such a vector as a signed integer. This uses the fact that in
   * the C++ language, booleans are implicitly convertible to integers. In
   * other words, declaring the type of the elements of the vector as a signed
   * integer is only a small lie, but it is a useful one.
   */
  using value_type = signed int;


  /**
   * Default constructor.
   */
  IndexSet();

  /**
   * Constructor that also sets the overall size of the index range.
   */
  explicit IndexSet(const size_type size);

  /**
   * Copy constructor.
   */
  IndexSet(const IndexSet &) = default;

  /**
   * Copy assignment operator.
   */
  IndexSet &
  operator=(const IndexSet &) = default;

  /**
   * Move constructor. Create a new IndexSet by transferring the internal data
   * of the input set.
   */
  IndexSet(IndexSet &&is) noexcept;

  /**
   * Move assignment operator. Transfer the internal data of the input set into
   * the current one.
   */
  IndexSet &
  operator=(IndexSet &&is) noexcept;

#ifdef DEAL_II_WITH_TRILINOS
  /**
   * Constructor from a Trilinos Epetra_BlockMap.
   */
  explicit IndexSet(const Epetra_BlockMap &map);
#endif

  /**
   * Remove all indices from this index set. The index set retains its size,
   * however.
   */
  void
  clear();

  /**
   * Set the maximal size of the indices upon which this object operates.
   *
   * This function can only be called if the index set does not yet contain
   * any elements.  This can be achieved by calling clear(), for example.
   */
  void
  set_size(const size_type size);

  /**
   * Return the size of the index space of which this index set is a subset
   * of.
   *
   * Note that the result is not equal to the number of indices within this
   * set. The latter information is returned by n_elements().
   */
  size_type
  size() const;

  /**
   * Add the half-open range $[\text{begin},\text{end})$ to the set of indices
   * represented by this class.
   * @param[in] begin The first element of the range to be added.
   * @param[in] end The past-the-end element of the range to be added.
   */
  void
  add_range(const size_type begin, const size_type end);

  /**
   * Add an individual index to the set of indices.
   */
  void
  add_index(const size_type index);

  /**
   * Add a whole set of indices described by dereferencing every element of
   * the iterator range <code>[begin,end)</code>.
   *
   * @param[in] begin Iterator to the first element of range of indices to be
   * added
   * @param[in] end The past-the-end iterator for the range of elements to be
   * added. @pre The condition <code>begin@<=end</code> needs to be satisfied.
   */
  template <typename ForwardIterator>
  void
  add_indices(const ForwardIterator &begin, const ForwardIterator &end);

  /**
   * Add the given IndexSet @p other to the current one, constructing the
   * union of *this and @p other.
   *
   * If the @p offset argument is nonzero, then every index in @p other is
   * shifted by @p offset before being added to the current index set. This
   * allows to construct, for example, one index set from several others that
   * are supposed to represent index sets corresponding to different ranges
   * (e.g., when constructing the set of nonzero entries of a block vector
   * from the sets of nonzero elements of the individual blocks of a vector).
   *
   * This function will generate an exception if any of the (possibly shifted)
   * indices of the @p other index set lie outside the range
   * <code>[0,size())</code> represented by the current object.
   */
  void
  add_indices(const IndexSet &other, const unsigned int offset = 0);

  /**
   * Return whether the specified index is an element of the index set.
   */
  bool
  is_element(const size_type index) const;

  /**
   * Return whether the index set stored by this object defines a contiguous
   * range. This is true also if no indices are stored at all.
   */
  bool
  is_contiguous() const;

  /**
   * Return whether the index set stored by this object contains no elements.
   * This is similar, but faster than checking <code>n_elements() == 0</code>.
   */
  bool
  is_empty() const;

  /**
   * Return whether the IndexSets are ascending with respect to MPI process
   * number and 1:1, i.e., each index is contained in exactly one IndexSet
   * (among those stored on the different processes), each process stores
   * contiguous subset of indices, and the index set on process $p+1$ starts
   * at the index one larger than the last one stored on process $p$.
   * In case there is only one MPI process, this just means that the IndexSet
   * is complete.
   */
  bool
  is_ascending_and_one_to_one(const MPI_Comm &communicator) const;

  /**
   * Return the number of elements stored in this index set.
   */
  size_type
  n_elements() const;

  /**
   * Return the global index of the local index with number @p local_index
   * stored in this index set. @p local_index obviously needs to be less than
   * n_elements().
   */
  size_type
  nth_index_in_set(const unsigned int local_index) const;

  /**
   * Return the how-manyth element of this set (counted in ascending order) @p
   * global_index is. @p global_index needs to be less than the size(). This
   * function returns numbers::invalid_dof_index if the index @p global_index is not actually
   * a member of this index set, i.e. if is_element(global_index) is false.
   */
  size_type
  index_within_set(const size_type global_index) const;

  /**
   * Each index set can be represented as the union of a number of contiguous
   * intervals of indices, where if necessary intervals may only consist of
   * individual elements to represent isolated members of the index set.
   *
   * This function returns the minimal number of such intervals that are
   * needed to represent the index set under consideration.
   */
  unsigned int
  n_intervals() const;

  /**
   * This function returns the local index of the beginning of the largest
   * range.
   */
  unsigned int
  largest_range_starting_index() const;

  /**
   * Compress the internal representation by merging individual elements with
   * contiguous ranges, etc. This function does not have any external effect.
   */
  void
  compress() const;

  /**
   * Comparison for equality of index sets. This operation is only allowed if
   * the size of the two sets is the same (though of course they do not have
   * to have the same number of indices).
   */
  bool
  operator==(const IndexSet &is) const;

  /**
   * Comparison for inequality of index sets. This operation is only allowed
   * if the size of the two sets is the same (though of course they do not
   * have to have the same number of indices).
   */
  bool
  operator!=(const IndexSet &is) const;

  /**
   * Return the intersection of the current index set and the argument given,
   * i.e. a set of indices that are elements of both index sets. The two index
   * sets must have the same size (though of course they do not have to have
   * the same number of indices).
   */
  IndexSet operator&(const IndexSet &is) const;

  /**
   * This command takes an interval <tt>[begin, end)</tt> and returns the
   * intersection of the current index set with the interval, shifted to the
   * range <tt>[0, end-begin)</tt>.
   *
   * In other words, the result of this operation is the intersection of the
   * set represented by the current object and the interval <tt>[begin,
   * end)</tt>, as seen <i>within the interval <tt>[begin, end)</tt></i> by
   * shifting the result of the intersection operation to the left by
   * <tt>begin</tt>. This corresponds to the notion of a <i>view</i>: The
   * interval <tt>[begin, end)</tt> is a <i>window</i> through which we see
   * the set represented by the current object.
   */
  IndexSet
  get_view(const size_type begin, const size_type end) const;

  /**
   * Remove all elements contained in @p other from this set. In other words,
   * if $x$ is the current object and $o$ the argument, then we compute $x
   * \leftarrow x \backslash o$.
   */
  void
  subtract_set(const IndexSet &other);

  /**
   * Remove and return the last element of the last range.
   * This function throws an exception if the IndexSet is empty.
   */
  size_type
  pop_back();

  /**
   * Remove and return the first element of the first range.
   * This function throws an exception if the IndexSet is empty.
   */
  size_type
  pop_front();

  /**
   * Fill the given vector with all indices contained in this IndexSet.
   */
  void
  fill_index_vector(std::vector<size_type> &indices) const;

  /**
   * Fill the given vector with either zero or one elements, providing a
   * binary representation of this index set. The given vector is assumed to
   * already have the correct size.
   *
   * The given argument is filled with integer values zero and one, using
   * <code>vector.operator[]</code>. Thus, any object that has such an
   * operator can be used as long as it allows conversion of integers zero and
   * one to elements of the vector. Specifically, this is the case for classes
   * Vector, BlockVector, but also std::vector@<bool@>, std::vector@<int@>,
   * and std::vector@<double@>.
   */
  template <typename VectorType>
  void
  fill_binary_vector(VectorType &vector) const;

  /**
   * Output a text representation of this IndexSet to the given stream. Used
   * for testing.
   */
  template <class StreamType>
  void
  print(StreamType &out) const;

  /**
   * Write the IndexSet into a text based file format, that can be read in
   * again using the read() function.
   */
  void
  write(std::ostream &out) const;

  /**
   * Construct the IndexSet from a text based representation given by the
   * stream @p in written by the write() function.
   */
  void
  read(std::istream &in);

  /**
   * Write the IndexSet into a binary, compact representation, that can be
   * read in again using the block_read() function.
   */
  void
  block_write(std::ostream &out) const;

  /**
   * Construct the IndexSet from a binary representation given by the stream
   * @p in written by the write_block() function.
   */
  void
  block_read(std::istream &in);

#ifdef DEAL_II_WITH_TRILINOS
  /**
   * Given an MPI communicator, create a Trilinos map object that represents a
   * distribution of vector elements or matrix rows in which we will locally
   * store those elements or rows for which we store the index in the current
   * index set, and all the other elements/rows elsewhere on one of the other
   * MPI processes.
   *
   * The last argument only plays a role if the communicator is a parallel
   * one, distributing computations across multiple processors. In that case,
   * if the last argument is false, then it is assumed that the index sets
   * this function is called with on all processors are mutually exclusive but
   * together enumerate each index exactly once. In other words, if you call
   * this function on two processors, then the index sets this function is
   * called with must together have all possible indices from zero to
   * size()-1, and no index must appear in both index sets. This corresponds,
   * for example, to the case where we want to split the elements of vectors
   * into unique subsets to be stored on different processors -- no element
   * should be owned by more than one processor, but each element must be
   * owned by one.
   *
   * On the other hand, if the second argument is true, then the index sets
   * can be overlapping, and they also do not need to span the whole index
   * set. This is a useful operation if we want to create vectors that not
   * only contain the locally owned indices, but for example also the elements
   * that correspond to degrees of freedom located on ghost cells. Another
   * application of this method is to select a subset of the elements of a
   * vector, e.g. for extracting only certain solution components.
   */
  Epetra_Map
  make_trilinos_map(const MPI_Comm &communicator = MPI_COMM_WORLD,
                    const bool      overlapping  = false) const;

#  ifdef DEAL_II_TRILINOS_WITH_TPETRA
  Tpetra::Map<int, types::global_dof_index>
  make_tpetra_map(const MPI_Comm &communicator = MPI_COMM_WORLD,
                  const bool      overlapping  = false) const;
#  endif
#endif


  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  std::size_t
  memory_consumption() const;

  DeclException1(ExcIndexNotPresent,
                 size_type,
                 << "The global index " << arg1
                 << " is not an element of this set.");

  /**
   * Write or read the data of this object to or from a stream for the purpose
   * of serialization
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version);


  /**
   * @name Iterators
   * @{
   */

  /**
   * Dereferencing an IntervalIterator will return a reference to an object of
   * this type. It allows access to a contiguous interval $[a,b[$ (also called
   * a range) of the IndexSet being iterated over.
   */
  class IntervalAccessor
  {
  public:
    /**
     * Construct a valid accessor given an IndexSet and the index @p range_idx
     * of the range to point to.
     */
    IntervalAccessor(const IndexSet *idxset, const size_type range_idx);

    /**
     * Construct an invalid accessor for the IndexSet.
     */
    explicit IntervalAccessor(const IndexSet *idxset);

    /**
     * Number of elements in this interval.
     */
    size_type
    n_elements() const;

    /**
     * If true, we are pointing at a valid interval in the IndexSet.
     */
    bool
    is_valid() const;

    /**
     * Return an iterator pointing at the first index in this interval.
     */
    ElementIterator
    begin() const;

    /**
     * Return an iterator pointing directly after the last index in this
     * interval.
     */
    ElementIterator
    end() const;

    /**
     * Return the index of the last index in this interval.
     */
    size_type
    last() const;

  private:
    /**
     * Private copy constructor.
     */
    IntervalAccessor(const IntervalAccessor &other);
    /**
     * Private copy operator.
     */
    IntervalAccessor &
    operator=(const IntervalAccessor &other);

    /**
     * Test for equality, used by IntervalIterator.
     */
    bool
    operator==(const IntervalAccessor &other) const;
    /**
     * Smaller-than operator, used by IntervalIterator.
     */
    bool
    operator<(const IntervalAccessor &other) const;
    /**
     * Advance this accessor to point to the next interval in the @p
     * index_set.
     */
    void
    advance();
    /**
     * Reference to the IndexSet.
     */
    const IndexSet *index_set;

    /**
     * Index into index_set.ranges[]. Set to numbers::invalid_dof_index if
     * invalid or the end iterator.
     */
    size_type range_idx;

    friend class IntervalIterator;
  };

  /**
   * Class that represents an iterator pointing to a contiguous interval
   * $[a,b[$ as returned by IndexSet::begin_interval().
   */
  class IntervalIterator
  {
  public:
    /**
     * Construct a valid iterator pointing to the interval with index @p
     * range_idx.
     */
    IntervalIterator(const IndexSet *idxset, const size_type range_idx);

    /**
     * Construct an invalid iterator (used as end()).
     */
    explicit IntervalIterator(const IndexSet *idxset);

    /**
     * Construct an empty iterator.
     */
    IntervalIterator();

    /**
     * Copy constructor from @p other iterator.
     */
    IntervalIterator(const IntervalIterator &other) = default;

    /**
     * Assignment of another iterator.
     */
    IntervalIterator &
    operator=(const IntervalIterator &other) = default;

    /**
     * Prefix increment.
     */
    IntervalIterator &
    operator++();

    /**
     * Postfix increment.
     */
    IntervalIterator
    operator++(int);

    /**
     * Dereferencing operator, returns an IntervalAccessor.
     */
    const IntervalAccessor &operator*() const;

    /**
     * Dereferencing operator, returns a pointer to an IntervalAccessor.
     */
    const IntervalAccessor *operator->() const;

    /**
     * Comparison.
     */
    bool
    operator==(const IntervalIterator &) const;

    /**
     * Inverse of <tt>==</tt>.
     */
    bool
    operator!=(const IntervalIterator &) const;

    /**
     * Comparison operator.
     */
    bool
    operator<(const IntervalIterator &) const;

    /**
     * Return the distance between the current iterator and the argument. The
     * distance is given by how many times one has to apply operator++ to the
     * current iterator to get the argument (for a positive return value), or
     * operator-- (for a negative return value).
     */
    int
    operator-(const IntervalIterator &p) const;

    /**
     * Mark the class as forward iterator and declare some alias which are
     * standard for iterators and are used by algorithms to enquire about the
     * specifics of the iterators they work on.
     */
    using iterator_category = std::forward_iterator_tag;
    using value_type        = IntervalAccessor;
    using difference_type   = std::ptrdiff_t;
    using pointer           = IntervalAccessor *;
    using reference         = IntervalAccessor &;

  private:
    /**
     * Accessor that contains what IndexSet and interval we are pointing at.
     */
    IntervalAccessor accessor;
  };

  /**
   * Class that represents an iterator pointing to a single element in the
   * IndexSet as returned by IndexSet::begin().
   */
  class ElementIterator
  {
  public:
    /**
     * Construct an iterator pointing to the global index @p index in the
     * interval @p range_idx
     */
    ElementIterator(const IndexSet *idxset,
                    const size_type range_idx,
                    const size_type index);

    /**
     * Construct an iterator pointing to the end of the IndexSet.
     */
    explicit ElementIterator(const IndexSet *idxset);

    /**
     * Dereferencing operator. The returned value is the index of the element
     * inside the IndexSet.
     */
    size_type operator*() const;

    /**
     * Does this iterator point to an existing element?
     */
    bool
    is_valid() const;

    /**
     * Prefix increment.
     */
    ElementIterator &
    operator++();

    /**
     * Postfix increment.
     */
    ElementIterator
    operator++(int);

    /**
     * Comparison.
     */
    bool
    operator==(const ElementIterator &) const;

    /**
     * Inverse of <tt>==</tt>.
     */
    bool
    operator!=(const ElementIterator &) const;

    /**
     * Comparison operator.
     */
    bool
    operator<(const ElementIterator &) const;

    /**
     * Return the distance between the current iterator and the argument. In
     * the expression <code>it_left-it_right</code> the distance is given by
     * how many times one has to apply operator++ to the right operand @p
     * it_right to get the left operand @p it_left (for a positive return
     * value), or to @p it_left to get the @p it_right (for a negative return
     * value).
     */
    std::ptrdiff_t
    operator-(const ElementIterator &p) const;

    /**
     * Mark the class as forward iterator and declare some alias which are
     * standard for iterators and are used by algorithms to enquire about the
     * specifics of the iterators they work on.
     */
    using iterator_category = std::forward_iterator_tag;
    using value_type        = size_type;
    using difference_type   = std::ptrdiff_t;
    using pointer           = size_type *;
    using reference         = size_type &;

  private:
    /**
     * Advance iterator by one.
     */
    void
    advance();

    /**
     * The parent IndexSet.
     */
    const IndexSet *index_set;
    /**
     * Index into index_set.ranges.
     */
    size_type range_idx;
    /**
     * The global index this iterator is pointing at.
     */
    size_type idx;
  };

  /**
   * Return an iterator that points at the first index that is contained in
   * this IndexSet.
   */
  ElementIterator
  begin() const;

  /**
   * Return an element iterator pointing to the element with global index
   * @p global_index or the next larger element if the index is not in the
   * set. This is equivalent to
   * @code
   * auto p = begin();
   * while (*p<global_index)
   *   ++p;
   * return p;
   * @endcode
   *
   * If there is no element in this IndexSet at or behind @p global_index,
   * this method will return end().
   */
  ElementIterator
  at(const size_type global_index) const;

  /**
   * Return an iterator that points one after the last index that is contained
   * in this IndexSet.
   */
  ElementIterator
  end() const;

  /**
   * Return an Iterator that points at the first interval of this IndexSet.
   */
  IntervalIterator
  begin_intervals() const;

  /**
   * Return an Iterator that points one after the last interval of this
   * IndexSet.
   */
  IntervalIterator
  end_intervals() const;

  /**
   * @}
   */

private:
  /**
   * A type that denotes the half open index range <code>[begin,end)</code>.
   *
   * The nth_index_in_set denotes the how many-th index within this IndexSet
   * the first element of the current range is. This information is only
   * accurate if IndexSet::compress() has been called after the last
   * insertion.
   */
  struct Range
  {
    size_type begin;
    size_type end;

    size_type nth_index_in_set;

    /**
     * Default constructor. Since there is no useful choice for a default
     * constructed interval, this constructor simply creates something that
     * resembles an invalid range. We need this constructor for serialization
     * purposes, but the invalid range should be filled with something read
     * from the archive before it is used, so we should hopefully never get to
     * see an invalid range in the wild.
     */
    Range();

    /**
     * Constructor. Create a half-open interval with the given indices.
     *
     * @param i1 Left end point of the interval.
     * @param i2 First index greater than the last index of the indicated
     * range.
     */
    Range(const size_type i1, const size_type i2);

    friend inline bool
    operator<(const Range &range_1, const Range &range_2)
    {
      return (
        (range_1.begin < range_2.begin) ||
        ((range_1.begin == range_2.begin) && (range_1.end < range_2.end)));
    }

    static bool
    end_compare(const IndexSet::Range &x, const IndexSet::Range &y)
    {
      return x.end < y.end;
    }

    static bool
    nth_index_compare(const IndexSet::Range &x, const IndexSet::Range &y)
    {
      return (x.nth_index_in_set + (x.end - x.begin) <
              y.nth_index_in_set + (y.end - y.begin));
    }

    friend inline bool
    operator==(const Range &range_1, const Range &range_2)
    {
      return ((range_1.begin == range_2.begin) && (range_1.end == range_2.end));
    }

    static std::size_t
    memory_consumption()
    {
      return sizeof(Range);
    }

    /**
     * Write or read the data of this object to or from a stream for the
     * purpose of serialization
     */
    template <class Archive>
    void
    serialize(Archive &ar, const unsigned int version);
  };

  /**
   * A set of contiguous ranges of indices that make up (part of) this index
   * set. This variable is always kept sorted.
   *
   * The variable is marked "mutable" so that it can be changed by compress(),
   * though this of course doesn't change anything about the external
   * representation of this index set.
   */
  mutable std::vector<Range> ranges;

  /**
   * True if compress() has been called after the last change in the set of
   * indices.
   *
   * The variable is marked "mutable" so that it can be changed by compress(),
   * though this of course doesn't change anything about the external
   * representation of this index set.
   */
  mutable bool is_compressed;

  /**
   * The overall size of the index range. Elements of this index set have to
   * have a smaller number than this value.
   */
  size_type index_space_size;

  /**
   * This integer caches the index of the largest range in @p ranges. This
   * gives <tt>O(1)</tt> access to the range with most elements, while general
   * access costs <tt>O(log(n_ranges))</tt>. The largest range is needed for
   * the methods @p is_element(), @p index_within_set(), @p nth_index_in_set.
   * In many applications, the largest range contains most elements (the
   * locally owned range), whereas there are only a few other elements
   * (ghosts).
   */
  mutable size_type largest_range;

  /**
   * A mutex that is used to synchronize operations of the do_compress()
   * function that is called from many 'const' functions via compress().
   */
  mutable Threads::Mutex compress_mutex;

  /**
   * Actually perform the compress() operation.
   */
  void
  do_compress() const;
};


/**
 * Create and return an index set of size $N$ that contains every single index
 * within this range. In essence, this function returns an index set created
 * by
 * @code
 *  IndexSet is (N);
 *  is.add_range(0, N);
 * @endcode
 * This function exists so that one can create and initialize index sets that
 * are complete in one step, or so one can write code like
 * @code
 *   if (my_index_set == complete_index_set(my_index_set.size())
 *     ...
 * @endcode
 *
 * @relatesalso IndexSet
 */
inline IndexSet
complete_index_set(const unsigned int N)
{
  IndexSet is(N);
  is.add_range(0, N);
  is.compress();
  return is;
}

/* ------------------ inline functions ------------------ */


/* IntervalAccessor */

inline IndexSet::IntervalAccessor::IntervalAccessor(
  const IndexSet *          idxset,
  const IndexSet::size_type range_idx)
  : index_set(idxset)
  , range_idx(range_idx)
{
  Assert(range_idx < idxset->n_intervals(),
         ExcInternalError("Invalid range index"));
}



inline IndexSet::IntervalAccessor::IntervalAccessor(const IndexSet *idxset)
  : index_set(idxset)
  , range_idx(numbers::invalid_dof_index)
{}



inline IndexSet::IntervalAccessor::IntervalAccessor(
  const IndexSet::IntervalAccessor &other)
  : index_set(other.index_set)
  , range_idx(other.range_idx)
{
  Assert(range_idx == numbers::invalid_dof_index || is_valid(),
         ExcMessage("invalid iterator"));
}



inline IndexSet::size_type
IndexSet::IntervalAccessor::n_elements() const
{
  Assert(is_valid(), ExcMessage("invalid iterator"));
  return index_set->ranges[range_idx].end - index_set->ranges[range_idx].begin;
}



inline bool
IndexSet::IntervalAccessor::is_valid() const
{
  return index_set != nullptr && range_idx < index_set->n_intervals();
}



inline IndexSet::ElementIterator
IndexSet::IntervalAccessor::begin() const
{
  Assert(is_valid(), ExcMessage("invalid iterator"));
  return {index_set, range_idx, index_set->ranges[range_idx].begin};
}



inline IndexSet::ElementIterator
IndexSet::IntervalAccessor::end() const
{
  Assert(is_valid(), ExcMessage("invalid iterator"));

  // point to first index in next interval unless we are the last interval.
  if (range_idx < index_set->ranges.size() - 1)
    return {index_set, range_idx + 1, index_set->ranges[range_idx + 1].begin};
  else
    return index_set->end();
}



inline IndexSet::size_type
IndexSet::IntervalAccessor::last() const
{
  Assert(is_valid(), ExcMessage("invalid iterator"));

  return index_set->ranges[range_idx].end - 1;
}



inline IndexSet::IntervalAccessor &
IndexSet::IntervalAccessor::operator=(const IndexSet::IntervalAccessor &other)
{
  index_set = other.index_set;
  range_idx = other.range_idx;
  Assert(range_idx == numbers::invalid_dof_index || is_valid(),
         ExcMessage("invalid iterator"));
  return *this;
}



inline bool
IndexSet::IntervalAccessor::
operator==(const IndexSet::IntervalAccessor &other) const
{
  Assert(index_set == other.index_set,
         ExcMessage(
           "Can not compare accessors pointing to different IndexSets"));
  return range_idx == other.range_idx;
}



inline bool
IndexSet::IntervalAccessor::
operator<(const IndexSet::IntervalAccessor &other) const
{
  Assert(index_set == other.index_set,
         ExcMessage(
           "Can not compare accessors pointing to different IndexSets"));
  return range_idx < other.range_idx;
}



inline void
IndexSet::IntervalAccessor::advance()
{
  Assert(
    is_valid(),
    ExcMessage(
      "Impossible to advance an IndexSet::IntervalIterator that is invalid"));
  ++range_idx;

  // set ourselves to invalid if we walk off the end
  if (range_idx >= index_set->ranges.size())
    range_idx = numbers::invalid_dof_index;
}


/* IntervalIterator */

inline IndexSet::IntervalIterator::IntervalIterator(
  const IndexSet *          idxset,
  const IndexSet::size_type range_idx)
  : accessor(idxset, range_idx)
{}



inline IndexSet::IntervalIterator::IntervalIterator()
  : accessor(nullptr)
{}



inline IndexSet::IntervalIterator::IntervalIterator(const IndexSet *idxset)
  : accessor(idxset)
{}



inline IndexSet::IntervalIterator &
IndexSet::IntervalIterator::operator++()
{
  accessor.advance();
  return *this;
}



inline IndexSet::IntervalIterator
IndexSet::IntervalIterator::operator++(int)
{
  const IndexSet::IntervalIterator iter = *this;
  accessor.advance();
  return iter;
}



inline const IndexSet::IntervalAccessor &IndexSet::IntervalIterator::
                                         operator*() const
{
  return accessor;
}



inline const IndexSet::IntervalAccessor *IndexSet::IntervalIterator::
                                         operator->() const
{
  return &accessor;
}



inline bool
IndexSet::IntervalIterator::
operator==(const IndexSet::IntervalIterator &other) const
{
  return accessor == other.accessor;
}



inline bool
IndexSet::IntervalIterator::
operator!=(const IndexSet::IntervalIterator &other) const
{
  return !(*this == other);
}



inline bool
IndexSet::IntervalIterator::
operator<(const IndexSet::IntervalIterator &other) const
{
  return accessor < other.accessor;
}



inline int
IndexSet::IntervalIterator::
operator-(const IndexSet::IntervalIterator &other) const
{
  Assert(accessor.index_set == other.accessor.index_set,
         ExcMessage(
           "Can not compare iterators belonging to different IndexSets"));

  const size_type lhs = (accessor.range_idx == numbers::invalid_dof_index) ?
                          accessor.index_set->ranges.size() :
                          accessor.range_idx;
  const size_type rhs =
    (other.accessor.range_idx == numbers::invalid_dof_index) ?
      accessor.index_set->ranges.size() :
      other.accessor.range_idx;

  if (lhs > rhs)
    return static_cast<int>(lhs - rhs);
  else
    return -static_cast<int>(rhs - lhs);
}



/* ElementIterator */

inline IndexSet::ElementIterator::ElementIterator(
  const IndexSet *          idxset,
  const IndexSet::size_type range_idx,
  const IndexSet::size_type index)
  : index_set(idxset)
  , range_idx(range_idx)
  , idx(index)
{
  Assert(range_idx < index_set->ranges.size(),
         ExcMessage(
           "Invalid range index for IndexSet::ElementIterator constructor."));
  Assert(
    idx >= index_set->ranges[range_idx].begin &&
      idx < index_set->ranges[range_idx].end,
    ExcInternalError(
      "Invalid index argument for IndexSet::ElementIterator constructor."));
}



inline IndexSet::ElementIterator::ElementIterator(const IndexSet *idxset)
  : index_set(idxset)
  , range_idx(numbers::invalid_dof_index)
  , idx(numbers::invalid_dof_index)
{}



inline bool
IndexSet::ElementIterator::is_valid() const
{
  Assert((range_idx == numbers::invalid_dof_index &&
          idx == numbers::invalid_dof_index) ||
           (range_idx < index_set->ranges.size() &&
            idx < index_set->ranges[range_idx].end),
         ExcInternalError("Invalid ElementIterator state."));

  return (range_idx < index_set->ranges.size() &&
          idx < index_set->ranges[range_idx].end);
}



inline IndexSet::size_type IndexSet::ElementIterator::operator*() const
{
  Assert(
    is_valid(),
    ExcMessage(
      "Impossible to dereference an IndexSet::ElementIterator that is invalid"));
  return idx;
}



inline bool
IndexSet::ElementIterator::
operator==(const IndexSet::ElementIterator &other) const
{
  Assert(index_set == other.index_set,
         ExcMessage(
           "Can not compare iterators belonging to different IndexSets"));
  return range_idx == other.range_idx && idx == other.idx;
}



inline void
IndexSet::ElementIterator::advance()
{
  Assert(
    is_valid(),
    ExcMessage(
      "Impossible to advance an IndexSet::ElementIterator that is invalid"));
  if (idx < index_set->ranges[range_idx].end)
    ++idx;
  // end of this range?
  if (idx == index_set->ranges[range_idx].end)
    {
      // point to first element in next interval if possible
      if (range_idx < index_set->ranges.size() - 1)
        {
          ++range_idx;
          idx = index_set->ranges[range_idx].begin;
        }
      else
        {
          // we just fell off the end, set to invalid:
          range_idx = numbers::invalid_dof_index;
          idx       = numbers::invalid_dof_index;
        }
    }
}



inline IndexSet::ElementIterator &
IndexSet::ElementIterator::operator++()
{
  advance();
  return *this;
}



inline IndexSet::ElementIterator
IndexSet::ElementIterator::operator++(int)
{
  const IndexSet::ElementIterator it = *this;
  advance();
  return it;
}



inline bool
IndexSet::ElementIterator::
operator!=(const IndexSet::ElementIterator &other) const
{
  return !(*this == other);
}



inline bool
IndexSet::ElementIterator::
operator<(const IndexSet::ElementIterator &other) const
{
  Assert(index_set == other.index_set,
         ExcMessage(
           "Can not compare iterators belonging to different IndexSets"));
  return range_idx < other.range_idx ||
         (range_idx == other.range_idx && idx < other.idx);
}



inline std::ptrdiff_t
IndexSet::ElementIterator::
operator-(const IndexSet::ElementIterator &other) const
{
  Assert(index_set == other.index_set,
         ExcMessage(
           "Can not compare iterators belonging to different IndexSets"));
  if (*this == other)
    return 0;
  if (!(*this < other))
    return -(other - *this);

  // only other can be equal to end() because of the checks above.
  Assert(is_valid(), ExcInternalError());

  // Note: we now compute how far advance *this in "*this < other" to get other,
  // so we need to return -c at the end.

  // first finish the current range:
  std::ptrdiff_t c = index_set->ranges[range_idx].end - idx;

  // now walk in steps of ranges (need to start one behind our current one):
  for (size_type range = range_idx + 1;
       range < index_set->ranges.size() && range <= other.range_idx;
       ++range)
    c += index_set->ranges[range].end - index_set->ranges[range].begin;

  Assert(
    other.range_idx < index_set->ranges.size() ||
      other.range_idx == numbers::invalid_dof_index,
    ExcMessage(
      "Inconsistent iterator state. Did you invalidate iterators by modifying the IndexSet?"));

  // We might have walked too far because we went until the end of
  // other.range_idx, so walk backwards to other.idx:
  if (other.range_idx != numbers::invalid_dof_index)
    c -= index_set->ranges[other.range_idx].end - other.idx;

  return -c;
}


/* Range */

inline IndexSet::Range::Range()
  : begin(numbers::invalid_dof_index)
  , end(numbers::invalid_dof_index)
  , nth_index_in_set(numbers::invalid_dof_index)
{}



inline IndexSet::Range::Range(const size_type i1, const size_type i2)
  : begin(i1)
  , end(i2)
  , nth_index_in_set(numbers::invalid_dof_index)
{}



/* IndexSet itself */

inline IndexSet::IndexSet()
  : is_compressed(true)
  , index_space_size(0)
  , largest_range(numbers::invalid_unsigned_int)
{}



inline IndexSet::IndexSet(const size_type size)
  : is_compressed(true)
  , index_space_size(size)
  , largest_range(numbers::invalid_unsigned_int)
{}



inline IndexSet::IndexSet(IndexSet &&is) noexcept
  : ranges(std::move(is.ranges))
  , is_compressed(is.is_compressed)
  , index_space_size(is.index_space_size)
  , largest_range(is.largest_range)
{
  is.ranges.clear();
  is.is_compressed    = true;
  is.index_space_size = 0;
  is.largest_range    = numbers::invalid_unsigned_int;

  compress();
}



inline IndexSet &
IndexSet::operator=(IndexSet &&is) noexcept
{
  ranges           = std::move(is.ranges);
  is_compressed    = is.is_compressed;
  index_space_size = is.index_space_size;
  largest_range    = is.largest_range;

  is.ranges.clear();
  is.is_compressed    = true;
  is.index_space_size = 0;
  is.largest_range    = numbers::invalid_unsigned_int;

  compress();

  return *this;
}



inline IndexSet::ElementIterator
IndexSet::begin() const
{
  compress();
  if (ranges.size() > 0)
    return {this, 0, ranges[0].begin};
  else
    return end();
}



inline IndexSet::ElementIterator
IndexSet::at(const size_type global_index) const
{
  compress();
  Assert(global_index < size(),
         ExcIndexRangeType<size_type>(global_index, 0, size()));

  if (ranges.empty())
    return end();

  std::vector<Range>::const_iterator main_range =
    ranges.begin() + largest_range;

  Range r(global_index, global_index + 1);
  // This optimization makes the bounds for lower_bound smaller by checking
  // the largest range first.
  std::vector<Range>::const_iterator range_begin, range_end;
  if (global_index < main_range->begin)
    {
      range_begin = ranges.begin();
      range_end   = main_range;
    }
  else
    {
      range_begin = main_range;
      range_end   = ranges.end();
    }

  // This will give us the first range p=[a,b[ with b>=global_index using
  // a binary search
  const std::vector<Range>::const_iterator p =
    Utilities::lower_bound(range_begin, range_end, r, Range::end_compare);

  // We couldn't find a range, which means we have no range that contains
  // global_index and also no range behind it, meaning we need to return end().
  if (p == ranges.end())
    return end();

  // Finally, we can have two cases: Either global_index is not in [a,b[,
  // which means we need to return an iterator to a because global_index, ...,
  // a-1 is not in the IndexSet (if branch). Alternatively, global_index is in
  // [a,b[ and we will return an iterator pointing directly at global_index
  // (else branch).
  if (global_index < p->begin)
    return {this, static_cast<size_type>(p - ranges.begin()), p->begin};
  else
    return {this, static_cast<size_type>(p - ranges.begin()), global_index};
}



inline IndexSet::ElementIterator
IndexSet::end() const
{
  compress();
  return IndexSet::ElementIterator(this);
}



inline IndexSet::IntervalIterator
IndexSet::begin_intervals() const
{
  compress();
  if (ranges.size() > 0)
    return IndexSet::IntervalIterator(this, 0);
  else
    return end_intervals();
}



inline IndexSet::IntervalIterator
IndexSet::end_intervals() const
{
  compress();
  return IndexSet::IntervalIterator(this);
}



inline void
IndexSet::clear()
{
  // reset so that there are no indices in the set any more; however,
  // as documented, the index set retains its size
  ranges.clear();
  is_compressed = true;
  largest_range = numbers::invalid_unsigned_int;
}



inline void
IndexSet::set_size(const size_type sz)
{
  Assert(ranges.empty(),
         ExcMessage("This function can only be called if the current "
                    "object does not yet contain any elements."));
  index_space_size = sz;
  is_compressed    = true;
}



inline IndexSet::size_type
IndexSet::size() const
{
  return index_space_size;
}



inline void
IndexSet::compress() const
{
  if (is_compressed == true)
    return;

  do_compress();
}



inline void
IndexSet::add_index(const size_type index)
{
  Assert(index < index_space_size,
         ExcIndexRangeType<size_type>(index, 0, index_space_size));

  const Range new_range(index, index + 1);
  if (ranges.size() == 0 || index > ranges.back().end)
    ranges.push_back(new_range);
  else if (index == ranges.back().end)
    ranges.back().end++;
  else
    ranges.insert(Utilities::lower_bound(ranges.begin(),
                                         ranges.end(),
                                         new_range),
                  new_range);
  is_compressed = false;
}



inline void
IndexSet::add_range(const size_type begin, const size_type end)
{
  Assert((begin < index_space_size) ||
           ((begin == index_space_size) && (end == index_space_size)),
         ExcIndexRangeType<size_type>(begin, 0, index_space_size));
  Assert(end <= index_space_size,
         ExcIndexRangeType<size_type>(end, 0, index_space_size + 1));
  Assert(begin <= end, ExcIndexRangeType<size_type>(begin, 0, end));

  if (begin != end)
    {
      const Range new_range(begin, end);

      // the new index might be larger than the last index present in the
      // ranges. Then we can skip the binary search
      if (ranges.size() == 0 || begin > ranges.back().end)
        ranges.push_back(new_range);
      else
        ranges.insert(Utilities::lower_bound(ranges.begin(),
                                             ranges.end(),
                                             new_range),
                      new_range);
      is_compressed = false;
    }
}



template <typename ForwardIterator>
inline void
IndexSet::add_indices(const ForwardIterator &begin, const ForwardIterator &end)
{
  if (begin == end)
    return;

  // identify ranges in the given iterator range by checking whether some
  // indices happen to be consecutive. to avoid quadratic complexity when
  // calling add_range many times (as add_range() going into the middle of an
  // already existing range must shift entries around), we first collect a
  // vector of ranges.
  std::vector<std::pair<size_type, size_type>> tmp_ranges;
  bool                                         ranges_are_sorted = true;
  for (ForwardIterator p = begin; p != end;)
    {
      const size_type begin_index = *p;
      size_type       end_index   = begin_index + 1;
      ForwardIterator q           = p;
      ++q;
      while ((q != end) && (*q == end_index))
        {
          ++end_index;
          ++q;
        }

      tmp_ranges.emplace_back(begin_index, end_index);
      p = q;

      // if the starting index of the next go-around of the for loop is less
      // than the end index of the one just identified, then we will have at
      // least one pair of ranges that are not sorted, and consequently the
      // whole collection of ranges is not sorted.
      if (p != end && *p < end_index)
        ranges_are_sorted = false;
    }

  if (!ranges_are_sorted)
    std::sort(tmp_ranges.begin(), tmp_ranges.end());

  // if we have many ranges, we first construct a temporary index set (where
  // we add ranges in a consecutive way, so fast), otherwise, we work with
  // add_range(). the number 9 is chosen heuristically given the fact that
  // there are typically up to 8 independent ranges when adding the degrees of
  // freedom on a 3D cell or 9 when adding degrees of freedom of faces. if
  // doing cell-by-cell additions, we want to avoid repeated calls to
  // IndexSet::compress() which gets called upon merging two index sets, so we
  // want to be in the other branch then.
  if (tmp_ranges.size() > 9)
    {
      IndexSet tmp_set(size());
      tmp_set.ranges.reserve(tmp_ranges.size());
      for (const auto &i : tmp_ranges)
        tmp_set.add_range(i.first, i.second);
      this->add_indices(tmp_set);
    }
  else
    for (const auto &i : tmp_ranges)
      add_range(i.first, i.second);
}



inline bool
IndexSet::is_element(const size_type index) const
{
  if (ranges.empty() == false)
    {
      compress();

      // fast check whether the index is in the largest range
      Assert(largest_range < ranges.size(), ExcInternalError());
      if (index >= ranges[largest_range].begin &&
          index < ranges[largest_range].end)
        return true;

      // get the element after which we would have to insert a range that
      // consists of all elements from this element to the end of the index
      // range plus one. after this call we know that if p!=end() then
      // p->begin<=index unless there is no such range at all
      //
      // if the searched for element is an element of this range, then we're
      // done. otherwise, the element can't be in one of the following ranges
      // because otherwise p would be a different iterator
      //
      // since we already know the position relative to the largest range (we
      // called compress!), we can perform the binary search on ranges with
      // lower/higher number compared to the largest range
      std::vector<Range>::const_iterator p = std::upper_bound(
        ranges.begin() +
          (index < ranges[largest_range].begin ? 0 : largest_range + 1),
        index < ranges[largest_range].begin ? ranges.begin() + largest_range :
                                              ranges.end(),
        Range(index, size() + 1));

      if (p == ranges.begin())
        return ((index >= p->begin) && (index < p->end));

      Assert((p == ranges.end()) || (p->begin > index), ExcInternalError());

      // now move to that previous range
      --p;
      Assert(p->begin <= index, ExcInternalError());

      return (p->end > index);
    }

  // didn't find this index, so it's not in the set
  return false;
}



inline bool
IndexSet::is_contiguous() const
{
  compress();
  return (ranges.size() <= 1);
}



inline bool
IndexSet::is_empty() const
{
  return ranges.empty();
}



inline IndexSet::size_type
IndexSet::n_elements() const
{
  // make sure we have non-overlapping ranges
  compress();

  size_type v = 0;
  if (!ranges.empty())
    {
      Range &r = ranges.back();
      v        = r.nth_index_in_set + r.end - r.begin;
    }

#ifdef DEBUG
  size_type s = 0;
  for (const auto &range : ranges)
    s += (range.end - range.begin);
  Assert(s == v, ExcInternalError());
#endif

  return v;
}



inline unsigned int
IndexSet::n_intervals() const
{
  compress();
  return ranges.size();
}



inline unsigned int
IndexSet::largest_range_starting_index() const
{
  Assert(ranges.empty() == false, ExcMessage("IndexSet cannot be empty."));

  compress();
  const std::vector<Range>::const_iterator main_range =
    ranges.begin() + largest_range;

  return main_range->nth_index_in_set;
}



inline IndexSet::size_type
IndexSet::nth_index_in_set(const unsigned int n) const
{
  Assert(n < n_elements(), ExcIndexRangeType<size_type>(n, 0, n_elements()));

  compress();

  // first check whether the index is in the largest range
  Assert(largest_range < ranges.size(), ExcInternalError());
  std::vector<Range>::const_iterator main_range =
    ranges.begin() + largest_range;
  if (n >= main_range->nth_index_in_set &&
      n < main_range->nth_index_in_set + (main_range->end - main_range->begin))
    return main_range->begin + (n - main_range->nth_index_in_set);

  // find out which chunk the local index n belongs to by using a binary
  // search. the comparator is based on the end of the ranges. Use the
  // position relative to main_range to subdivide the ranges
  Range r(n, n + 1);
  r.nth_index_in_set = n;
  std::vector<Range>::const_iterator range_begin, range_end;
  if (n < main_range->nth_index_in_set)
    {
      range_begin = ranges.begin();
      range_end   = main_range;
    }
  else
    {
      range_begin = main_range + 1;
      range_end   = ranges.end();
    }

  const std::vector<Range>::const_iterator p =
    Utilities::lower_bound(range_begin, range_end, r, Range::nth_index_compare);

  Assert(p != ranges.end(), ExcInternalError());
  return p->begin + (n - p->nth_index_in_set);
}



inline IndexSet::size_type
IndexSet::index_within_set(const size_type n) const
{
  // to make this call thread-safe, compress() must not be called through this
  // function
  Assert(is_compressed == true, ExcMessage("IndexSet must be compressed."));
  Assert(n < size(), ExcIndexRangeType<size_type>(n, 0, size()));

  // return immediately if the index set is empty
  if (is_empty())
    return numbers::invalid_dof_index;

  // check whether the index is in the largest range. use the result to
  // perform a one-sided binary search afterward
  Assert(largest_range < ranges.size(), ExcInternalError());
  std::vector<Range>::const_iterator main_range =
    ranges.begin() + largest_range;
  if (n >= main_range->begin && n < main_range->end)
    return (n - main_range->begin) + main_range->nth_index_in_set;

  Range                              r(n, n);
  std::vector<Range>::const_iterator range_begin, range_end;
  if (n < main_range->begin)
    {
      range_begin = ranges.begin();
      range_end   = main_range;
    }
  else
    {
      range_begin = main_range + 1;
      range_end   = ranges.end();
    }

  std::vector<Range>::const_iterator p =
    Utilities::lower_bound(range_begin, range_end, r, Range::end_compare);

  // if n is not in this set
  if (p == range_end || p->end == n || p->begin > n)
    return numbers::invalid_dof_index;

  Assert(p != ranges.end(), ExcInternalError());
  Assert(p->begin <= n, ExcInternalError());
  Assert(n < p->end, ExcInternalError());
  return (n - p->begin) + p->nth_index_in_set;
}



inline bool
IndexSet::operator==(const IndexSet &is) const
{
  Assert(size() == is.size(), ExcDimensionMismatch(size(), is.size()));

  compress();
  is.compress();

  return ranges == is.ranges;
}



inline bool
IndexSet::operator!=(const IndexSet &is) const
{
  Assert(size() == is.size(), ExcDimensionMismatch(size(), is.size()));

  compress();
  is.compress();

  return ranges != is.ranges;
}



template <typename Vector>
void
IndexSet::fill_binary_vector(Vector &vector) const
{
  Assert(vector.size() == size(), ExcDimensionMismatch(vector.size(), size()));

  compress();
  // first fill all elements of the vector with zeroes.
  std::fill(vector.begin(), vector.end(), 0);

  // then write ones into the elements whose indices are contained in the
  // index set
  for (const auto &range : ranges)
    for (size_type i = range.begin; i < range.end; ++i)
      vector[i] = 1;
}



template <class StreamType>
inline void
IndexSet::print(StreamType &out) const
{
  compress();
  out << "{";
  std::vector<Range>::const_iterator p;
  for (p = ranges.begin(); p != ranges.end(); ++p)
    {
      if (p->end - p->begin == 1)
        out << p->begin;
      else
        out << "[" << p->begin << "," << p->end - 1 << "]";

      if (p != --ranges.end())
        out << ", ";
    }
  out << "}" << std::endl;
}



template <class Archive>
inline void
IndexSet::Range::serialize(Archive &ar, const unsigned int)
{
  ar &begin &end &nth_index_in_set;
}



template <class Archive>
inline void
IndexSet::serialize(Archive &ar, const unsigned int)
{
  ar &ranges &is_compressed &index_space_size &largest_range;
}

DEAL_II_NAMESPACE_CLOSE

#endif
