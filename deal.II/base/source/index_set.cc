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
				   // see if any of the
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

      unsigned int first_index = i->first;
      unsigned int last_index  = i->second;

				       // see if we can merge any of
				       // the following ranges
      bool can_merge = false;
      while (next != contiguous_ranges.end() &&
	     (next->first <= last_index))
	{
	  last_index = next->second;
	  ++next;
	  can_merge = true;
	}

      if (can_merge == true)
	{
					   // delete the old
					   // ranges and insert the
					   // new range
	  contiguous_ranges.erase (i, next);
	  contiguous_ranges.insert (ContiguousRange(first_index,
						    last_index));

					   // then set iterator to
					   // the same position as
					   // before. in the next
					   // step, the loop is then
					   // going to increase 'i'
					   // by one. note that we
					   // don't need to
					   // reconsider the same
					   // iterator for merging
					   // with the next range
					   // because we have made
					   // sure we gobble up as
					   // many of the following
					   // ranges as possible
	  i = contiguous_ranges.lower_bound (ContiguousRange(first_index,
							     last_index));
	}
    }
}


DEAL_II_NAMESPACE_CLOSE
