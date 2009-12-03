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
#  include <Epetra_Map.h>
#endif

DEAL_II_NAMESPACE_OPEN

void
IndexSet::do_compress () const
{
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



IndexSet
IndexSet::operator & (const IndexSet &is) const
{
  Assert (size() == is.size(),
	  ExcDimensionMismatch (size(), is.size()));

  compress ();
  is.compress ();

  std::vector<Range>::const_iterator r1 = ranges.begin(),
				     r2 = is.ranges.begin();
  IndexSet result (size());

  while ((r1 != ranges.end())
	 &&
	 (r2 != is.ranges.end()))
    {
				       // if r1 and r2 do not overlap
				       // at all, then move the
				       // pointer that sits to the
				       // left of the other up by one
      if (r1->end <= r2->begin)
	++r1;
      else if (r2->end <= r1->begin)
	++r2;
      else
	{
					   // the ranges must overlap
					   // somehow
	  Assert (((r1->begin <= r2->begin) &&
		   (r1->end > r2->begin))
		  ||
		  ((r2->begin <= r1->begin) &&
		   (r2->end > r1->begin)),
		  ExcInternalError());

					   // add the overlapping
					   // range to the result
	  result.add_range (std::max (r1->begin,
				      r2->begin),
			    std::min (r1->end,
				      r2->end));

					   // now move that iterator
					   // that ends earlier one
					   // up. note that it has to
					   // be this one because a
					   // subsequent range may
					   // still have a chance of
					   // overlapping with the
					   // range that ends later
	  if (r1->end <= r2->end)
	    ++r1;
	  else
	    ++r2;
	}
    }

  result.compress ();
  return result;
}



IndexSet
IndexSet::get_view (const unsigned int begin,
		    const unsigned int end) const
{
  Assert (begin <= end,
	  ExcMessage ("End index needs to be larger or equal to begin index!"));
  Assert (end <= size(),
	  ExcMessage ("Given range exceeds index set dimension"));

  IndexSet result (end-begin);
  std::vector<Range>::const_iterator r1 = ranges.begin();

  while (r1 != ranges.end())
    {
      if ((r1->end > begin)
	  &&
	  (r1->begin < end))
	{
	  result.add_range (std::max(r1->begin, begin)-begin,
			    std::min(r1->end, end)-begin);

	}
      else
	if (r1->begin >= end)
	  break;

      ++r1;
    }

  result.compress();
  return result;
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
