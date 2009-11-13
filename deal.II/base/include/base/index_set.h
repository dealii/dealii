//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__index_set_h
#define __deal2__index_set_h

#include <base/config.h>
#include <base/exceptions.h>

#include <vector>
#include <algorithm>

#ifdef DEAL_II_USE_TRILINOS
#  include <Epetra_Map.h>
#endif

#if defined(DEAL_II_COMPILER_SUPPORTS_MPI) || defined(DEAL_II_USE_PETSC)
#include <mpi.h>
#else
typedef int MPI_Comm;
#define MPI_COMM_WORLD 0
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
				      * Constructor that also sets the
				      * overall size of the index
				      * range.
				      */
    IndexSet (const unsigned int size);

				     /**
				      * Set the maximal size of the
				      * indices upon which this object
				      * operates.
				      *
				      * This function can only be
				      * called if the index set does
				      * not yet contain any elements.
				      */
    void set_size (const unsigned int size);

				     /**
				      * Return the size of the index
				      * space of which this index set
				      * is a subset of.
				      */
    unsigned int size () const;

				     /**
				      * Add the half-open range
				      * $[\text{begin},\text{end})$ to
				      * the set of indices represented
				      * by this class.
				      */
    void add_range (const unsigned int begin,
		    const unsigned int end);

				     /**
				      * Add an individual index to the
				      * set of indices.
				      */
    void add_index (const unsigned int index);

				     /**
				      * Add a whole set of indices
				      * described by dereferencing
				      * every element of the the
				      * iterator range
				      * <code>[begin,end)</code>.
				      */
    template <typename ForwardIterator>
    void add_indices (const ForwardIterator &begin,
		      const ForwardIterator &end);

				     /**
				      * Return whether the specified
				      * index is an element of the
				      * index set.
				      */
    bool is_element (const unsigned int index) const;

				     /**
				      * Return whether the index set
				      * stored by this object defines
				      * a contiguous range. This is
				      * true also if no indices are
				      * stored at all.
				      */
    bool is_contiguous () const;

				     /**
				      * Return the number of elements
				      * stored in this index set.
				      */
    unsigned int n_elements () const;

				     /**
				      * Return the nth index stored in
				      * this index set. @p n obviously
				      * needs to be less than
				      * n_elements().
				      */
    unsigned int nth_index_in_set (const unsigned int n) const;

				     /**
				      * Return the how-manyth element
				      * of this set (counted in
				      * ascending order) @p n is. @p n
				      * needs to be less than the
				      * size(). This function throws
				      * an exception if the index @p n
				      * is not actually a member of
				      * this index set, i.e. if
				      * is_element(n) is false.
				      */
    unsigned int index_within_set (const unsigned int n) const;

				     /**
				      * Compress the internal
				      * representation by merging
				      * individual elements with
				      * contiguous ranges, etc. This
				      * function does not have any
				      * external effect.
				      */
    void compress () const;

#ifdef DEAL_II_USE_TRILINOS
				     /**
				      * Given an MPI communicator,
				      * create a Trilinos map object
				      * that represents a distribution
				      * of vector elements or matrix
				      * rows in which we will locally
				      * store those elements or rows
				      * for which we store the index
				      * in the current index set, and
				      * all the other elements/rows
				      * elsewhere on one of the other
				      * MPI processes.
				      *
				      * The last argument only plays a
				      * role if the communicator is a
				      * parallel one, distributing
				      * computations across multiple
				      * processors. In that case, if
				      * the last argument is false,
				      * then it is assumed that the
				      * index sets this function is
				      * called on on all processors
				      * are mutually exclusive but
				      * together enumerate each index
				      * exactly once. In other words,
				      * if you call this function on
				      * two processors, then the index
				      * sets this function is called
				      * with must together have all
				      * possible indices from zero to
				      * size()-1, and no index must
				      * appear in both index
				      * sets. This corresponds, for
				      * example, to the case where we
				      * want to split the elements of
				      * vectors into unique subsets to
				      * be stored on different
				      * processors -- no element
				      * should be owned by more than
				      * one processor, but each
				      * element must be owned by one.
				      *
				      * On the other hand, if the
				      * second argument is true, then
				      * the index sets can be
				      * overlapping, though they still
				      * need to contain each index
				      * exactly once on all processors
				      * taken together. This is a
				      * useful operation if we want to
				      * create vectors that not only
				      * contain the locally owned
				      * indices, but for example also
				      * the elements that correspond
				      * to degrees of freedom located
				      * on ghost cells.
				      */
    Epetra_Map make_trilinos_map (const MPI_Comm &communicator = MPI_COMM_WORLD,
				  const bool      overlapping  = false) const;
#endif

  private:
				     /**
				      * A type that denotes the half
				      * open index range
				      * <code>[begin,end)</code>.
				      *
				      * The nth_index_in_set denotes
				      * the how many-th index within
				      * this IndexSet the first
				      * element of the current range
				      * is. This information is only
				      * accurate if
				      * IndexSet::compress() has been
				      * called after the last
				      * insertion.
				      */
    struct Range
    {
	unsigned int begin;
	unsigned int end;

	unsigned int nth_index_in_set;

	Range (const unsigned int i1,
	       const unsigned int i2);

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
    };

				     /**
				      * A set of contiguous ranges of
				      * indices that make up (part of)
				      * this index set. This variable
				      * is always kept sorted.
				      *
				      * The variable is marked
				      * "mutable" so that it can be
				      * changed by compress(), though
				      * this of course doesn't change
				      * anything about the external
				      * representation of this index
				      * set.
				      */
    mutable std::vector<Range> ranges;

				     /**
				      * True if compress() has been
				      * called after the last change
				      * in the set of indices.
				      *
				      * The variable is marked
				      * "mutable" so that it can be
				      * changed by compress(), though
				      * this of course doesn't change
				      * anything about the external
				      * representation of this index
				      * set.
				      */
    mutable bool is_compressed;

				     /**
				      * The overall size of the index
				      * range. Elements of this index
				      * set have to have a smaller
				      * number than this value.
				      */
    unsigned int index_space_size;
};


/* ------------------ inline functions ------------------ */

inline
IndexSet::Range::Range (const unsigned int i1,
			const unsigned int i2)
		:
		begin(i1),
		end(i2)
{}



inline
IndexSet::IndexSet ()
		:
		is_compressed (true),
		index_space_size (0)
{}



inline
IndexSet::IndexSet (const unsigned int size)
		:
		is_compressed (true),
		index_space_size (size)
{}



inline
void
IndexSet::set_size (const unsigned int sz)
{
  Assert (ranges.size() == 0,
	  ExcMessage ("This function can only be called if the current "
		      "object does not yet contain any elements."));
  index_space_size = sz;
  is_compressed = true;
}



inline
unsigned int
IndexSet::size () const
{
  return index_space_size;
}


inline
void
IndexSet::add_range (const unsigned int begin,
		     const unsigned int end)
{
  Assert (begin < index_space_size,
	  ExcIndexRange (begin, 0, index_space_size));
  Assert (end <= index_space_size,
	  ExcIndexRange (end, 0, index_space_size+1));
  Assert (begin <= end,
	  ExcIndexRange (begin, 0, end));

  if (begin != end)
    {
      const Range new_range(begin,end);
      ranges.insert (std::lower_bound (ranges.begin(),
				       ranges.end(),
				       new_range),
		     new_range);
      is_compressed = false;
    }
}



inline
void
IndexSet::add_index (const unsigned int index)
{
  Assert (index < index_space_size,
	  ExcIndexRange (index, 0, index_space_size));

  const Range new_range(index, index+1);
  ranges.insert (std::lower_bound (ranges.begin(),
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
				   // insert each element of the
				   // range. if some of them happen to
				   // be consecutive, merge them to a
				   // range
  for (ForwardIterator p=begin; p!=end;)
    {
      const unsigned int begin_index = *p;
      unsigned int       end_index   = begin_index + 1;
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
IndexSet::is_element (const unsigned int index) const
{
  if (ranges.size() > 0)
    {
      compress ();

				       // get the element after which
				       // we would have to insert a
				       // range that consists of all
				       // elements from this element
				       // to the end of the index
				       // range plus one. after this
				       // call we know that if
				       // p!=end() then
				       // p->begin<=index unless there
				       // is no such range at all
				       //
				       // if the searched for element
				       // is an element of this range,
				       // then we're done. otherwise,
				       // the element can't be in one
				       // of the following ranges
				       // because otherwise p would be
				       // a different iterator
      std::vector<Range>::const_iterator
	p = std::upper_bound (ranges.begin(),
			      ranges.end(),
			      Range (index, size()+1));

      if (p == ranges.begin())
	return ((index >= p->begin) && (index < p->end));

      Assert ((p == ranges.end()) || (p->begin > index),
		ExcInternalError());

				       // now move to that previous
				       // range
      --p;
      Assert (p->begin <= index, ExcInternalError());

      return (p->end > index);
    }

				   // didn't find this index, so it's
				   // not in the set
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
unsigned int
IndexSet::n_elements () const
{
				   // make sure we have
				   // non-overlapping ranges
  compress ();

  unsigned int s = 0;
  for (std::vector<Range>::iterator range = ranges.begin();
       range != ranges.end();
       ++range)
    s += (range->end - range->begin);

  return s;
}



inline
unsigned int
IndexSet::nth_index_in_set (const unsigned int n) const
{
  Assert (n < n_elements(), ExcIndexRange (n, 0, n_elements()));

//TODO: this could be done more efficiently by using a binary search
  for (std::vector<Range>::const_iterator p = ranges.begin();
       p != ranges.end(); ++p)
    {
      Assert (n >= p->nth_index_in_set, ExcInternalError());
      if (n < p->nth_index_in_set + (p->end-p->begin))
	return p->begin + (n-p->nth_index_in_set);
    }

  Assert (false, ExcInternalError());
  return numbers::invalid_unsigned_int;
}



inline
unsigned int
IndexSet::index_within_set (const unsigned int n) const
{
  Assert (is_element(n) == true,
	  ExcMessage ("Given number is not an element of this set."));
  Assert (n < size(), ExcIndexRange (n, 0, size()));

//TODO: this could be done more efficiently by using a binary search
  for (std::vector<Range>::const_iterator p = ranges.begin();
       p != ranges.end(); ++p)
    {
      Assert (n >= p->begin, ExcInternalError());
      if (n < p->end)
	return (n-p->begin) + p->nth_index_in_set;
    }

  Assert (false, ExcInternalError());
  return numbers::invalid_unsigned_int;
}





DEAL_II_NAMESPACE_CLOSE

#endif

