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
  if (is_compressed == true)
    return;
				   // see if any of the
				   // contiguous ranges can be
				   // merged. since they are sorted by
				   // their first index, determining
				   // overlap isn't all that hard
  for (std::set<Range>::iterator
	 i = ranges.begin();
       i != ranges.end();
       ++i)
    {
      std::set<Range>::const_iterator
	next = i;
      ++next;

      unsigned int first_index = i->begin;
      unsigned int last_index  = i->end;

				       // see if we can merge any of
				       // the following ranges
      bool can_merge = false;
      while (next != ranges.end() &&
	     (next->begin <= last_index))
	{
	  last_index = std::max (last_index, next->end);
	  ++next;
	  can_merge = true;
	}

      if (can_merge == true)
	{
					   // delete the old
					   // ranges and insert the
					   // new range
	  ranges.erase (i, next);
	  ranges.insert (Range(first_index,
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
	  i = ranges.lower_bound (Range(first_index,
					last_index));
	}
    }

  is_compressed = true;
}


DEAL_II_NAMESPACE_CLOSE
