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


				   // we can roll any of
				   // the individual entries into the
				   // contiguous ranges
}


DEAL_II_NAMESPACE_CLOSE
