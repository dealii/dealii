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
