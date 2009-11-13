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

#ifdef DEAL_II_USE_TRILINOS
#  ifdef DEAL_II_COMPILER_SUPPORTS_MPI
#    include <Epetra_MpiComm.h>
#  endif
#  include <Epetra_SerialComm.h>
#endif

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



#ifdef DEAL_II_USE_TRILINOS

Epetra_Map
IndexSet::make_trilinos_map (const MPI_Comm &communicator,
			     const bool overlapping) const
{
  compress ();

  if ((is_contiguous() == true) && (!overlapping))
    return Epetra_Map (size(),
		       n_elements(),
		       0,
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
		       Epetra_MpiComm(communicator));
#else
		       Epetra_SerialComm());
#endif
  else
    {
      std::vector<int> indices;
      indices.reserve(n_elements());
      for (std::vector<Range>::iterator
	     i = ranges.begin();
	   i != ranges.end();
	   ++i)
	for (unsigned int j=i->begin; j<i->end; ++j)
	  indices.push_back (j);
      Assert (indices.size() == n_elements(), ExcInternalError());
  
      return Epetra_Map (-1,
			 n_elements(),
			 &indices[0],
			 0,
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
			 Epetra_MpiComm(communicator));
#else
			 Epetra_SerialComm());
#endif
    }
}

  
#endif


DEAL_II_NAMESPACE_CLOSE
