//---------------------------------------------------------------------------
//      $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2006, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <base/index_set.h>

DEAL_II_NAMESPACE_OPEN

void
IndexSet::compress () const
{
				   // first see if any of the
				   // individual entries can be merged
				   // to form contiguous ranges
  {
    std::set<unsigned int>::iterator p = individual_indices.begin();
    while (p != individual_indices.end())
      {
	unsigned int n = 0;
	std::set<unsigned int>::iterator q = p;
	++q;
	while ((q != individual_indices.end())
	       &&
	       (*q == *p+n+1))
	  {
	    ++q;
	    ++n;
	  }

	if (n>0)
	  {
	    const unsigned int range_end = *q;

					     // add the new range,
					     // delete the elements of
					     // the range from the set
					     // of individual indices
	    contiguous_ranges.insert (ContiguousRange(*p, *p+n+1));
	    individual_indices.erase (p, q);

					     // reset p to the element
					     // past the last
	    p = individual_indices.lower_bound (range_end);
	  }
	else
	  ++p;
      }
  }


				   // next see if we can roll any of
				   // the individual entries into the
				   // contiguous ranges
  {
    std::set<ContiguousRange, RangeComparison>::iterator
      next_range            = contiguous_ranges.begin();
    std::set<unsigned int>::iterator
      next_individual_index = individual_indices.begin();

    while ((next_range != contiguous_ranges.end())
	   &&
	   (next_individual_index != individual_indices.end()))
      {
					 // if the next individual
					 // index is beyond the upper
					 // bound of the next range,
					 // then we need to inspect
					 // the next range
	if (*next_individual_index > next_range->second)
	  {
	    ++next_range;
	    continue;
	  }

					 // if, on the other hand, the
					 // next individual index is
					 // below the lower bound of
					 // the next range, then we
					 // need to inspect the next
					 // index
	if (*next_individual_index < next_range->first-1)
	  {
	    ++next_individual_index;
	    continue;
	  }

					 // now we know that the
					 // individual index is within
					 // the next range, or
					 // possibly one element
					 // before, or one element
					 // after. so they can be
					 // merged. note that ranges
					 // are half open, so in the
					 // upper bound we have to
					 // consider the next index+1
	ContiguousRange new_range (std::min (*next_individual_index,
					     next_range->first),
				   std::max (*next_individual_index+1,
					     next_range->second));

					 // erase old index and range,
					 // and reset iterators
					 //
					 // note that the following
					 // index can't possibly be
					 // *next_index+1, since we
					 // have merged contiguous
					 // sets of individual indices
					 // above. however, there may
					 // be more individual indices
					 // that lie within the
					 // current range, so we'll
					 // delete all of them
	const unsigned int old_index = *next_individual_index;
	individual_indices.erase (next_individual_index);
	next_individual_index = individual_indices.lower_bound (old_index);

	contiguous_ranges.erase (next_range);
	contiguous_ranges.insert (new_range);

	next_range = contiguous_ranges.lower_bound (new_range);
      }
  }


				   // finally see if any of the
				   // contiguous ranges can be
				   // merged. since they are sorted by
				   // their first index, determining
				   // overlap isn't all that hard
  for (std::set<ContiguousRange, RangeComparison>::iterator
	 i = contiguous_ranges.begin();
       i != contiguous_ranges.end();
       ++i)
    {
      std::set<ContiguousRange, RangeComparison>::const_iterator
	next = i;
      ++next;

      if (next != contiguous_ranges.end())
	if (next->first <= i->second)
	  {
	    const unsigned int first_index = i->first;
	    const unsigned int last_index  = next->second;

					     // delete the two old
					     // ranges and insert the
					     // new range
	    contiguous_ranges.erase (i, ++next);
	    contiguous_ranges.insert (ContiguousRange(first_index,
						      last_index));

					     // then set iterator to
					     // the same position as
					     // before
	    i = contiguous_ranges.lower_bound (ContiguousRange(first_index,
							       last_index));
	  }
    }
}


DEAL_II_NAMESPACE_CLOSE
