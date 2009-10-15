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
  for (std::vector<Range>::iterator
	 i = ranges.begin();
       i != ranges.end(); )
    {
      std::vector<Range>::iterator
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
					   // delete the old ranges
					   // and insert the new range
					   // in place of the previous
					   // one
	  *i = Range(first_index, last_index);
	  i = ranges.erase (i+1, next);
	}
      else
	++i;
    }


				   // now compute indices within set
  unsigned int next_index = 0;
  for (std::vector<Range>::iterator
	 i = ranges.begin();
       i != ranges.end();
       ++i)
    {
      i->nth_index_in_set = next_index;
      next_index += (i->end - i->begin);
    }
  is_compressed = true;

				   // check that next_index is
				   // correct. needs to be after the
				   // previous statement because we
				   // otherwise will get into an
				   // endless loop
  Assert (next_index == n_elements(), ExcInternalError());
}


DEAL_II_NAMESPACE_CLOSE
