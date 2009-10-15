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

#include <set>


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

  private:
				     /**
				      * A type that denotes an index
				      * range.
				      */
    typedef std::pair<unsigned int, unsigned int> ContiguousRange;

				     /**
				      * A structure that provides an
				      * ordering to ranges.
				      */
    struct RangeComparison
    {
	bool operator() (const ContiguousRange &range_1,
			 const ContiguousRange &range_2) const;
    };

				     /**
				      * A set of contiguous ranges of
				      * indices that make up (part of)
				      * this index set.
				      *
				      * The variable is marked
				      * "mutable" so that it can be
				      * changed by compress(), though
				      * this of course doesn't change
				      * anything about the external
				      * representation of this index
				      * set.
				      */
    mutable std::set<ContiguousRange, RangeComparison> contiguous_ranges;

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
bool
IndexSet::RangeComparison::operator() (const ContiguousRange &range_1,
				       const ContiguousRange &range_2) const
{
  return ((range_1.first < range_2.first)
	  ||
	  ((range_1.first == range_2.first)
	   &&
	   (range_1.second < range_2.second)));
}



inline
IndexSet::IndexSet ()
		:
		index_space_size (0)
{}



inline
IndexSet::IndexSet (const unsigned int size)
		:
		index_space_size (size)
{}



inline
void
IndexSet::set_size (const unsigned int sz)
{
  Assert (contiguous_ranges.size() == 0,
	  ExcMessage ("This function can only be called if the current "
		      "object does not yet contain any elements."));
  index_space_size = sz;
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
				       // if it turns out to be a
				       // single element then add that
				       // separately
      if (end == begin+1)
	add_index (begin);
      else
	contiguous_ranges.insert (ContiguousRange(begin,end));
    }
}



inline
void
IndexSet::add_index (const unsigned int index)
{
  Assert (index < index_space_size,
	  ExcIndexRange (index, 0, index_space_size));

  contiguous_ranges.insert (ContiguousRange(index, index+1));
}



inline
bool
IndexSet::is_element (const unsigned int index) const
{
//TODO: since the list of ranges is sorted, we be faster here by doing a binary search
				   // if not, we need to walk the list
				   // of contiguous ranges:
  if (contiguous_ranges.size() > 0)
    for (std::set<ContiguousRange, RangeComparison>::const_iterator
	   i = contiguous_ranges.begin();
	 i != contiguous_ranges.end();
	 ++i)
      if ((index >= i->first) && (index < i->second))
	return true;

				   // didn't find this index, so it's
				   // not in the set
  return false;
}



inline
bool
IndexSet::is_contiguous () const
{
  compress ();
  return (contiguous_ranges.size() <= 1);
}



inline
unsigned int
IndexSet::n_elements () const
{
				   // make sure we have
				   // non-overlapping ranges and
				   // individual indices
  compress ();

  unsigned int s = 0;
  for (std::set<ContiguousRange, RangeComparison>::iterator
	 range = contiguous_ranges.begin();
       range != contiguous_ranges.end();
       ++range)
    s += (range->second - range->first);

  return s;
}



inline
unsigned int
IndexSet::nth_index_in_set (const unsigned int n) const
{
  Assert (n < n_elements(), ExcIndexRange (n, 0, n_elements()));

  Assert (is_contiguous(), ExcNotImplemented());
  return (n+contiguous_ranges.begin()->first);
}



inline
unsigned int
IndexSet::index_within_set (const unsigned int n) const
{
  Assert (is_element(n) == true,
	  ExcMessage ("Given number is not an element of this set."));
  Assert (n < size(), ExcIndexRange (n, 0, size()));

  Assert (is_contiguous(), ExcNotImplemented());
  return (n-contiguous_ranges.begin()->first);
}





DEAL_II_NAMESPACE_CLOSE

#endif

