// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2014 by the deal.II authors
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

#ifndef __deal2__index_set_h
#define __deal2__index_set_h

#include <deal.II/base/config.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/exceptions.h>
#include <vector>
#include <algorithm>

#ifdef DEAL_II_WITH_TRILINOS
#  include <Epetra_Map.h>
#endif

#if defined(DEAL_II_WITH_MPI) || defined(DEAL_II_WITH_PETSC)
#include <mpi.h>
#else
typedef int MPI_Comm;
#  ifndef MPI_COMM_WORLD
#    define MPI_COMM_WORLD 0
#  endif
#endif

DEAL_II_NAMESPACE_OPEN

/**
 * A class that represents a subset of indices among a larger set. For
 * example, it can be used to denote the set of degrees of freedom
 * within the range $[0,\text{dof\_handler.n\_dofs})$ that belongs to
 * a particular subdomain, or those among all degrees of freedom that
 * are stored on a particular processor in a distributed parallel
 * computation.
 *
 * This class can represent a collection of half-open ranges of
 * indices as well as individual elements. For practical purposes it
 * also stores the overall range these indices can assume. In other
 * words, you need to specify the size of the index space
 * $[0,\text{size})$ of which objects of this class are a subset.
 *
 * The data structures used in this class along with a rationale can be found
 * in the @ref distributed_paper "Distributed Computing paper".
 *
 * @author Wolfgang Bangerth, 2009
 */
class IndexSet
{
public:
  /**
   * Default constructor.
   */
  IndexSet ();

  /**
   * Constructor that also sets the overall size of the index range.
   */
  explicit IndexSet (const types::global_dof_index size);

  /**
   * Remove all indices from this index set. The index set retains its size,
   * however.
   */
  void clear ();

  /**
   * Set the maximal size of the indices upon which this object operates.
   *
   * This function can only be called if the index set does not yet contain
   * any elements.  This can be achieved by calling clear(), for example.
   */
  void set_size (const types::global_dof_index size);

  /**
   * Return the size of the index space of which this index set is a subset
   * of.
   *
   * Note that the result is not equal to the number of indices within this
   * set. The latter information is returned by n_elements().
   */
  types::global_dof_index size () const;

  /**
   * Add the half-open range $[\text{begin},\text{end})$ to the set of indices
   * represented by this class.
   */
  void add_range (const types::global_dof_index begin,
                  const types::global_dof_index end);

  /**
   * Add an individual index to the set of indices.
   */
  void add_index (const types::global_dof_index index);

  /**
   * Add a whole set of indices described by dereferencing every element of
   * the the iterator range <code>[begin,end)</code>.
   */
  template <typename ForwardIterator>
  void add_indices (const ForwardIterator &begin,
                    const ForwardIterator &end);

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
  void add_indices(const IndexSet &other,
                   const unsigned int offset = 0);

  /**
   * Return whether the specified index is an element of the index set.
   */
  bool is_element (const types::global_dof_index index) const;

  /**
   * Return whether the index set stored by this object defines a contiguous
   * range. This is true also if no indices are stored at all.
   */
  bool is_contiguous () const;

  /**
   * Return the number of elements stored in this index set.
   */
  types::global_dof_index n_elements () const;

  /**
   * Return the global index of the local index with number @p local_index
   * stored in this index set. @p local_index obviously needs to be less than
   * n_elements().
   */
  types::global_dof_index nth_index_in_set (const unsigned int local_index) const;

  /**
   * Return the how-manyth element of this set (counted in ascending order) @p
   * global_index is. @p global_index needs to be less than the size(). This
   * function throws an exception if the index @p global_index is not actually
   * a member of this index set, i.e. if is_element(global_index) is false.
   */
  types::global_dof_index index_within_set (const types::global_dof_index global_index) const;

  /**
   * Each index set can be represented as the union of a number of contiguous
   * intervals of indices, where if necessary intervals may only consist of
   * individual elements to represent isolated members of the index set.
   *
   * This function returns the minimal number of such intervals that are
   * needed to represent the index set under consideration.
   */
  unsigned int n_intervals () const;

  /**
   * Compress the internal representation by merging individual elements with
   * contiguous ranges, etc. This function does not have any external effect.
   */
  void compress () const;

  /**
   * Comparison for equality of index sets. This operation is only allowed if
   * the size of the two sets is the same (though of course they do not have
   * to have the same number of indices).
   */
  bool operator == (const IndexSet &is) const;

  /**
   * Comparison for inequality of index sets. This operation is only allowed
   * if the size of the two sets is the same (though of course they do not
   * have to have the same number of indices).
   */
  bool operator != (const IndexSet &is) const;

  /**
   * Return the intersection of the current index set and the argument given,
   * i.e. a set of indices that are elements of both index sets. The two index
   * sets must have the same size (though of course they do not have to have
   * the same number of indices).
   */
  IndexSet operator & (const IndexSet &is) const;

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
  IndexSet get_view (const types::global_dof_index begin,
                     const types::global_dof_index end) const;


  /**
   * Removes all elements contained in @p other from this set. In other words,
   * if $x$ is the current object and $o$ the argument, then we compute $x
   * \leftarrow x \backslash o$.
   */
  void subtract_set (const IndexSet &other);


  /**
   * Fills the given vector with all indices contained in this IndexSet.
   */
  void fill_index_vector(std::vector<types::global_dof_index> &indices) const;

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
  template <typename Vector>
  void fill_binary_vector (Vector &vector) const;

  /**
   * Outputs a text representation of this IndexSet to the given stream. Used
   * for testing.
   */
  template <class STREAM>
  void print(STREAM &out) const;

  /**
   * Writes the IndexSet into a text based file format, that can be read in
   * again using the read() function.
   */
  void write(std::ostream &out) const;

  /**
   * Constructs the IndexSet from a text based representation given by the
   * stream @p in written by the write() function.
   */
  void read(std::istream &in);

  /**
   * Writes the IndexSet into a binary, compact representation, that can be
   * read in again using the block_read() function.
   */
  void block_write(std::ostream &out) const;

  /**
   * Constructs the IndexSet from a binary representation given by the stream
   * @p in written by the write_block() function.
   */
  void block_read(std::istream &in);

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
   * this function is called on on all processors are mutually exclusive but
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
   * can be overlapping, though they still need to contain each index exactly
   * once on all processors taken together. This is a useful operation if we
   * want to create vectors that not only contain the locally owned indices,
   * but for example also the elements that correspond to degrees of freedom
   * located on ghost cells.
   */
  Epetra_Map make_trilinos_map (const MPI_Comm &communicator = MPI_COMM_WORLD,
                                const bool      overlapping  = false) const;
#endif


  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  std::size_t memory_consumption () const;

  DeclException1 (ExcIndexNotPresent, types::global_dof_index,
                  << "The global index " << arg1
                  << " is not an element of this set.");

  /**
   * Write or read the data of this object to or
   * from a stream for the purpose of serialization
   */
  template <class Archive>
  void serialize (Archive &ar, const unsigned int version);

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
    types::global_dof_index begin;
    types::global_dof_index end;

    types::global_dof_index nth_index_in_set;

    /**
     * Default constructor. Since there is no useful choice for
     * a default constructed interval, this constructor simply
     * creates something that resembles an invalid range. We
     * need this constructor for serialization purposes, but the
     * invalid range should be filled with something read from
     * the archive before it is used, so we should hopefully
     * never get to see an invalid range in the wild.
     **/
    Range ();

    /**
     * Constructor. Create a half-open interval with the given indices.
     *
     * @param i1 Left end point of the interval.
     * @param i2 First index greater than the last index of the indicated range.
     **/
    Range (const types::global_dof_index i1,
           const types::global_dof_index i2);

    friend
    inline bool operator< (const Range &range_1,
                           const Range &range_2)
    {
      return ((range_1.begin < range_2.begin)
              ||
              ((range_1.begin == range_2.begin)
               &&
               (range_1.end < range_2.end)));
    }

    static bool end_compare(const IndexSet::Range &x, const IndexSet::Range &y)
    {
      return x.end < y.end;
    }

    static bool nth_index_compare (const IndexSet::Range &x,
                                   const IndexSet::Range &y)
    {
      return (x.nth_index_in_set+(x.end-x.begin) <
              y.nth_index_in_set+(y.end-y.begin));
    }

    friend
    inline bool operator== (const Range &range_1,
                            const Range &range_2)
    {
      return ((range_1.begin == range_2.begin)
              &&
              (range_1.end == range_2.end));
    }

    std::size_t memory_consumption () const
    {
      return sizeof(Range);
    }

    /**
     * Write or read the data of this object to or from a stream for the
     * purpose of serialization
     */
    template <class Archive>
    void serialize (Archive &ar, const unsigned int version);
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
  types::global_dof_index index_space_size;

  /**
   * This integer caches the index of the largest range in @p ranges. This
   * gives <tt>O(1)</tt> access to the range with most elements, while general
   * access costs <tt>O(log(n_ranges))</tt>. The largest range is needed for
   * the methods @p is_element(), @p index_within_set(), @p
   * nth_index_in_set. In many applications, the largest range contains most
   * elements (the locally owned range), whereas there are only a few other
   * elements (ghosts).
   */
  mutable types::global_dof_index largest_range;

  /**
   * Actually perform the compress() operation.
   */
  void do_compress() const;
};


/**
 * Create and return an index set of size $N$ that contains every
 * single index within this range. In essence, this function
 * returns an index set created by
 * @code
 *  IndexSet is (N);
 *  is.add_range(0, N);
 * @endcode
 * This function exists so that one can create and initialize
 * index sets that are complete in one step, or so one can write
 * code like
 * @code
 *   if (my_index_set == complete_index_set(my_index_set.size())
 *     ...
 * @endcode
 *
 * @relates IndexSet
 */
inline
IndexSet complete_index_set (const unsigned int N)
{
  IndexSet is (N);
  is.add_range(0, N);
  return is;
}

/* ------------------ inline functions ------------------ */

inline
IndexSet::Range::Range ()
  :
  begin(numbers::invalid_dof_index),
  end(numbers::invalid_dof_index)
{}


inline
IndexSet::Range::Range (const types::global_dof_index i1,
                        const types::global_dof_index i2)
  :
  begin(i1),
  end(i2)
{}


inline
IndexSet::IndexSet ()
  :
  is_compressed (true),
  index_space_size (0),
  largest_range (deal_II_numbers::invalid_unsigned_int)
{}



inline
IndexSet::IndexSet (const types::global_dof_index size)
  :
  is_compressed (true),
  index_space_size (size),
  largest_range (deal_II_numbers::invalid_unsigned_int)
{}



inline
void
IndexSet::clear ()
{
  ranges.clear ();
  largest_range = 0;
  is_compressed = true;
}


inline
void
IndexSet::set_size (const types::global_dof_index sz)
{
  Assert (ranges.empty(),
          ExcMessage ("This function can only be called if the current "
                      "object does not yet contain any elements."));
  index_space_size = sz;
  is_compressed = true;
}



inline
types::global_dof_index
IndexSet::size () const
{
  return index_space_size;
}



inline
void
IndexSet::compress () const
{
  if (is_compressed == true)
    return;

  do_compress();
}



inline
void
IndexSet::add_index (const types::global_dof_index index)
{
  Assert (index < index_space_size,
          ExcIndexRangeType<types::global_dof_index> (index, 0, index_space_size));

  const Range new_range(index, index+1);
  if (ranges.size() == 0 || index > ranges.back().end)
    ranges.push_back(new_range);
  else if (index == ranges.back().end)
    ranges.back().end++;
  else
    ranges.insert (Utilities::lower_bound (ranges.begin(),
                                           ranges.end(),
                                           new_range),
                   new_range);
  is_compressed = false;
}



template <typename ForwardIterator>
inline
void
IndexSet::add_indices (const ForwardIterator &begin,
                       const ForwardIterator &end)
{
  // insert each element of the range. if some of them happen to be
  // consecutive, merge them to a range
  for (ForwardIterator p=begin; p!=end;)
    {
      const types::global_dof_index begin_index = *p;
      types::global_dof_index       end_index   = begin_index + 1;
      ForwardIterator q = p;
      ++q;
      while ((q != end) && (*q == end_index))
        {
          ++end_index;
          ++q;
        }

      add_range (begin_index, end_index);
      p = q;
    }
}



inline
bool
IndexSet::is_element (const types::global_dof_index index) const
{
  if (ranges.empty() == false)
    {
      compress ();

      // fast check whether the index is in the largest range
      Assert (largest_range < ranges.size(), ExcInternalError());
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
      std::vector<Range>::const_iterator
      p = std::upper_bound (ranges.begin() + (index<ranges[largest_range].begin?
                                              0 : largest_range+1),
                            index<ranges[largest_range].begin ?
                            ranges.begin() + largest_range:
                            ranges.end(),
                            Range (index, size()+1));

      if (p == ranges.begin())
        return ((index >= p->begin) && (index < p->end));

      Assert ((p == ranges.end()) || (p->begin > index),
              ExcInternalError());

      // now move to that previous range
      --p;
      Assert (p->begin <= index, ExcInternalError());

      return (p->end > index);
    }

  // didn't find this index, so it's not in the set
  return false;
}



inline
bool
IndexSet::is_contiguous () const
{
  compress ();
  return (ranges.size() <= 1);
}



inline
types::global_dof_index
IndexSet::n_elements () const
{
  // make sure we have non-overlapping ranges
  compress ();

  types::global_dof_index v = 0;
  if (!ranges.empty())
    {
      Range &r = ranges.back();
      v = r.nth_index_in_set + r.end - r.begin;
    }

#ifdef DEBUG
  types::global_dof_index s = 0;
  for (std::vector<Range>::iterator range = ranges.begin();
       range != ranges.end();
       ++range)
    s += (range->end - range->begin);
  Assert(s==v, ExcInternalError());
#endif

  return v;
}



inline
unsigned int
IndexSet::n_intervals () const
{
  compress ();
  return ranges.size();
}



inline
types::global_dof_index
IndexSet::nth_index_in_set (const unsigned int n) const
{
  // to make this call thread-safe, compress() must not be called through this
  // function
  Assert (is_compressed == true, ExcMessage ("IndexSet must be compressed."));
  Assert (n < n_elements(), ExcIndexRangeType<types::global_dof_index> (n, 0, n_elements()));

  // first check whether the index is in the largest range
  Assert (largest_range < ranges.size(), ExcInternalError());
  std::vector<Range>::const_iterator main_range=ranges.begin()+largest_range;
  if (n>=main_range->nth_index_in_set &&
      n<main_range->nth_index_in_set+(main_range->end-main_range->begin))
    return main_range->begin + (n-main_range->nth_index_in_set);

  // find out which chunk the local index n belongs to by using a binary
  // search. the comparator is based on the end of the ranges. Use the
  // position relative to main_range to subdivide the ranges
  Range r (n,n+1);
  r.nth_index_in_set = n;
  std::vector<Range>::const_iterator range_begin, range_end;
  if (n<main_range->nth_index_in_set)
    {
      range_begin = ranges.begin();
      range_end   = main_range;
    }
  else
    {
      range_begin = main_range + 1;
      range_end   = ranges.end();
    }

  const std::vector<Range>::const_iterator
  p = Utilities::lower_bound(range_begin, range_end, r,
                             Range::nth_index_compare);

  Assert (p != ranges.end(), ExcInternalError());
  return p->begin + (n-p->nth_index_in_set);
}



inline
types::global_dof_index
IndexSet::index_within_set (const types::global_dof_index n) const
{
  // to make this call thread-safe, compress() must not be called through this
  // function
  Assert (is_compressed == true, ExcMessage ("IndexSet must be compressed."));
  Assert (is_element(n) == true, ExcIndexNotPresent (n));
  Assert (n < size(), ExcIndexRangeType<types::global_dof_index> (n, 0, size()));

  // check whether the index is in the largest range. use the result to
  // perform a one-sided binary search afterward
  Assert (largest_range < ranges.size(), ExcInternalError());
  std::vector<Range>::const_iterator main_range=ranges.begin()+largest_range;
  if (n >= main_range->begin && n < main_range->end)
    return (n-main_range->begin) + main_range->nth_index_in_set;

  Range r(n, n);
  std::vector<Range>::const_iterator range_begin, range_end;
  if (n<main_range->begin)
    {
      range_begin = ranges.begin();
      range_end   = main_range;
    }
  else
    {
      range_begin = main_range + 1;
      range_end   = ranges.end();
    }

  std::vector<Range>::const_iterator
  p = Utilities::lower_bound(range_begin, range_end, r,
                             Range::end_compare);

  Assert(p!=ranges.end(), ExcInternalError());
  Assert(p->begin<=n, ExcInternalError());
  Assert(n<p->end, ExcInternalError());
  return (n-p->begin) + p->nth_index_in_set;
}



inline
bool
IndexSet::operator == (const IndexSet &is) const
{
  Assert (size() == is.size(),
          ExcDimensionMismatch (size(), is.size()));

  compress ();
  is.compress ();

  return ranges == is.ranges;
}



inline
bool
IndexSet::operator != (const IndexSet &is) const
{
  Assert (size() == is.size(),
          ExcDimensionMismatch (size(), is.size()));

  compress ();
  is.compress ();

  return ranges != is.ranges;
}



template <typename Vector>
void
IndexSet::fill_binary_vector (Vector &vector) const
{
  Assert (vector.size() == size(),
          ExcDimensionMismatch (vector.size(), size()));

  compress();
  // first fill all elements of the vector with zeroes.
  std::fill (vector.begin(), vector.end(), 0);

  // then write ones into the elements whose indices are contained in the
  // index set
  for (std::vector<Range>::iterator it = ranges.begin();
       it != ranges.end();
       ++it)
    for (types::global_dof_index i=it->begin; i<it->end; ++i)
      vector[i] = 1;
}



template <class STREAM>
inline
void
IndexSet::print (STREAM &out) const
{
  compress();
  out << "{";
  std::vector<Range>::const_iterator p;
  for (p = ranges.begin(); p != ranges.end(); ++p)
    {
      if (p->end-p->begin==1)
        out << p->begin;
      else
        out << "[" << p->begin << "," << p->end-1 << "]";

      if (p !=--ranges.end())
        out << ", ";
    }
  out << "}" << std::endl;
}



template <class Archive>
inline
void
IndexSet::Range::serialize (Archive &ar, const unsigned int)
{
  ar &begin &end &nth_index_in_set;
}



template <class Archive>
inline
void
IndexSet::serialize (Archive &ar, const unsigned int)
{
  ar &ranges &is_compressed &index_space_size &largest_range;
}

DEAL_II_NAMESPACE_CLOSE

#endif
