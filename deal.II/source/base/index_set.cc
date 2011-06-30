//---------------------------------------------------------------------------
//      $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2006, 2008, 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/index_set.h>
#include <list>

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


				   // now compute indices within set and the
				   // range with most elements
  unsigned int next_index = 0, largest_range_size = 0;
  for (std::vector<Range>::iterator
	 i = ranges.begin();
       i != ranges.end();
       ++i)
    {
      Assert(i->begin < i->end, ExcInternalError());

      i->nth_index_in_set = next_index;
      next_index += (i->end - i->begin);
      if (i->end - i->begin > largest_range_size)
	{
	  largest_range_size = i->end - i->begin;
	  largest_range = i - ranges.begin();
	}
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



unsigned int
IndexSet::n_intervals () const
{
  compress ();
  return ranges.size();
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



void
IndexSet::write(std::ostream & out) const
{
  compress();
  out << size() << " ";
  out << ranges.size() << std::endl;
  std::vector<Range>::const_iterator r = ranges.begin();
  for ( ;r!=ranges.end(); ++r)
    {
      out << r->begin << " " << r->end << std::endl;
    }
}



void
IndexSet::read(std::istream & in)
{
  unsigned int s, numranges, b, e;
  in >> s >> numranges;
  ranges.clear();
  set_size(s);
  for (unsigned int i=0;i<numranges;++i)
    {
      in >> b >> e;
      add_range(b,e);
    }
}


void
IndexSet::block_write(std::ostream & out) const
{
  Assert (out, ExcIO());
  out.write(reinterpret_cast<const char*>(&index_space_size),
	    sizeof(index_space_size));
  size_t n_ranges = ranges.size();
  out.write(reinterpret_cast<const char*>(&n_ranges),
	    sizeof(n_ranges));
  if (ranges.empty() == false)
    out.write (reinterpret_cast<const char*>(&*ranges.begin()),
	       ranges.size() * sizeof(Range));
  Assert (out, ExcIO());
}

void
IndexSet::block_read(std::istream & in)
{
  unsigned int size;
  size_t n_ranges;
  in.read(reinterpret_cast<char*>(&size), sizeof(size));
  in.read(reinterpret_cast<char*>(&n_ranges), sizeof(n_ranges));
				   // we have to clear ranges first
  ranges.clear();
  set_size(size);
  ranges.resize(n_ranges, Range(0,0));
  if (n_ranges)
    in.read(reinterpret_cast<char*>(&*ranges.begin()),
	    ranges.size() * sizeof(Range));
}



void
IndexSet::subtract_set (const IndexSet & other)
{
  compress();
  other.compress();
  is_compressed = false;


				   // we save new ranges to be added to our
				   // IndexSet in an temporary list and add
				   // all of them in one go at the end. This
				   // is necessary because a growing ranges
				   // vector invalidates iterators.
  std::list<Range> temp_list;

  std::vector<Range>::iterator own_it = ranges.begin();
  std::vector<Range>::iterator other_it = other.ranges.begin();

  while (own_it != ranges.end() && other_it != other.ranges.end())
    {
				       //advance own iterator until we get an
				       //overlap
      if (own_it->end <= other_it->begin)
	{
	  ++own_it;
	  continue;
	}
				       //we are done with other_it, so advance
      if (own_it->begin >= other_it->end)
	{
	  ++other_it;
	  continue;
	}

				       //Now own_it and other_it overlap.
				       //First save the part of own_it that is
				       //before other_it (if not empty).
      if (own_it->begin < other_it->begin)
	{
	  Range r(own_it->begin, other_it->begin);
	  r.nth_index_in_set = 0; //fix warning of unused variable
	  temp_list.push_back(r);
	}
				       // change own_it to the sub range
				       // behind other_it. Do not delete
				       // own_it in any case. As removal would
				       // invalidate iterators, we just shrink
				       // the range to an empty one.
      own_it->begin = other_it->end;
      if (own_it->begin > own_it->end)
	{
	  own_it->begin = own_it->end;
	  ++own_it;
	}

				       // continue without advancing
				       // iterators, the right one will be
				       // advanced next.
    }

				   // Now delete all empty ranges we might
				   // have created.
  for (std::vector<Range>::iterator it = ranges.begin();
       it != ranges.end(); )
    {
      if (it->begin >= it->end)
	it = ranges.erase(it);
      else
	++it;
    }

				   // done, now add the temporary ranges
  for (std::list<Range>::iterator it = temp_list.begin();
       it != temp_list.end();
       ++it)
    add_range(it->begin, it->end);

  compress();
}


void IndexSet::fill_index_vector(std::vector<unsigned int> & indices) const
{
  compress();

  indices.clear();
  indices.reserve(n_elements());

  for (std::vector<Range>::iterator it = ranges.begin();
       it != ranges.end();
       ++it)
    for (unsigned int i=it->begin; i<it->end; ++i)
      indices.push_back (i);

  Assert (indices.size() == n_elements(), ExcInternalError());
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
      std::vector<unsigned int> indices;
      fill_index_vector(indices);

      int * indices_ptr = reinterpret_cast<int*>(&indices[0]);

      return Epetra_Map (-1,
			 n_elements(),
			 indices_ptr,
			 0,
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
			 Epetra_MpiComm(communicator));
#else
			 Epetra_SerialComm());
      (void)communicator;
#endif
    }
}


#endif


std::size_t
IndexSet::memory_consumption () const
{
  return MemoryConsumption::memory_consumption (ranges) +
    MemoryConsumption::memory_consumption (is_compressed) +
    MemoryConsumption::memory_consumption (index_space_size);
}



DEAL_II_NAMESPACE_CLOSE
