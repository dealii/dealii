//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/compressed_simple_sparsity_pattern.h>
#include <base/memory_consumption.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <functional>

DEAL_II_NAMESPACE_OPEN



template <typename ForwardIterator>
void
CompressedSimpleSparsityPattern::Line::add_entries (ForwardIterator begin,
						    ForwardIterator end,
						    const bool      indices_are_sorted)
{
  const int n_elements = end - begin;
  if (n_elements <= 0)
    return;

  const unsigned int stop_size = entries.size() + n_elements;

  if (indices_are_sorted == true && n_elements > 3)
    {
				   // in debug mode, check whether the
				   // indices really are sorted.
#ifdef DEBUG
      {
	ForwardIterator test = begin, test1 = begin;
	++test1;
	for ( ; test1 != end; ++test, ++test1)
	  Assert (*test1 > *test, ExcInternalError());
      }
#endif

      if (entries.size() == 0 || entries.back() < *begin)
	{
	  entries.insert(entries.end(), begin, end);
	  return;
	}

				   // find a possible insertion point for
				   // the first entry. check whether the
				   // first entry is a duplicate before
				   // actually doing something.
      ForwardIterator my_it = begin;
      unsigned int col = *my_it;
      std::vector<unsigned int>::iterator it =
	std::lower_bound(entries.begin(), entries.end(), col);
      while (*it == col)
	{
	  ++my_it;
	  if (my_it == end)
	    break;
	  col = *my_it;
				   // check the very next entry in the
				   // current array
	  ++it;
	  if (it == entries.end())
	    break;
	  if (*it > col)
	    break;
	  if (*it == col)
	    continue;
				   // ok, it wasn't the very next one, do a
				   // binary search to find the insert point
	  it = std::lower_bound(it, entries.end(), col);
	  if (it == entries.end())
	    break;
	}
				   // all input entries were duplicates.
      if (my_it == end)
	return;

				   // resize vector by just inserting the
				   // list
      const unsigned int pos1 = it - entries.begin();
      Assert (pos1 <= entries.size(), ExcInternalError());
      entries.insert (it, my_it, end);
      it = entries.begin() + pos1;
      Assert (entries.size() >= (unsigned int)(it-entries.begin()), ExcInternalError());

				   // now merge the two lists.
      std::vector<unsigned int>::iterator it2 = it + (end-my_it);

				   // as long as there are indices both in
				   // the end of the entries list and in the
				   // input list
      while (my_it != end && it2 != entries.end())
	{
	  if (*my_it < *it2)
	    *it++ = *my_it++;
	  else if (*my_it == *it2)
	    {
	      *it++ = *it2++;
	      ++my_it;
	    }
	  else
	    *it++ = *it2++;
	}
				   // in case there are indices left in the
				   // input list
      while (my_it != end)
	*it++ = *my_it++;

				   // in case there are indices left in the
				   // end of entries
      while (it2 != entries.end())
	*it++ = *it2++;

				   // resize and return
      const unsigned int new_size = it - entries.begin();
      Assert (new_size <= stop_size, ExcInternalError());
      entries.resize (new_size);
      return;
    }

				   // unsorted case or case with too few
				   // elements
  ForwardIterator my_it = begin;

				   // If necessary, increase the size of the
				   // array.
  if (stop_size > entries.capacity())
    entries.reserve (stop_size);

  unsigned int col = *my_it;
  std::vector<unsigned int>::iterator it, it2;
				   // insert the first element as for one
				   // entry only first check the last
				   // element (or if line is still empty)
  if ( (entries.size()==0) || ( entries.back() < col) ) {
    entries.push_back(col);
    it = entries.end()-1;
  }
  else {
				   // do a binary search to find the place
				   // where to insert:
    it2 = std::lower_bound(entries.begin(), entries.end(), col);

				   // If this entry is a duplicate, continue
				   // immediately Insert at the right place
				   // in the vector. Vector grows
				   // automatically to fit elements. Always
				   // doubles its size.
    if (*it2 != col)
      it = entries.insert(it2, col);
    else
      it = it2;
  }

  ++my_it;
				   // Now try to be smart and insert with
				   // bias in the direction we are
				   // walking. This has the advantage that
				   // for sorted lists, we always search in
				   // the right direction, what should
				   // decrease the work needed in here.
  for ( ; my_it != end; ++my_it)
    {
      col = *my_it;
				   // need a special insertion command when
				   // we're at the end of the list
      if (col > entries.back()) {
	entries.push_back(col);
	it = entries.end()-1;
      }
				   // search to the right (preferred search
				   // direction)
      else if (col > *it) {
      	it2 = std::lower_bound(it++, entries.end(), col);
	if (*it2 != col)
	  it = entries.insert(it2, col);
      }
				   // search to the left
      else if (col < *it) {
	it2 = std::lower_bound(entries.begin(), it, col);
	if (*it2 != col)
	  it = entries.insert(it2, col);
      }
				   // if we're neither larger nor smaller,
				   // then this was a duplicate and we can
				   // just continue.
    }
}


unsigned int
CompressedSimpleSparsityPattern::Line::memory_consumption () const
{
  return entries.capacity()*sizeof(unsigned int)+sizeof(Line);
}


CompressedSimpleSparsityPattern::CompressedSimpleSparsityPattern ()
                :
		rows(0),
		cols(0),
		rowset(0)
{}



CompressedSimpleSparsityPattern::
CompressedSimpleSparsityPattern (const CompressedSimpleSparsityPattern &s)
                :
		Subscriptor(),
		rows(0),
		cols(0),
		rowset(0)
{
  Assert (s.rows == 0, ExcInvalidConstructorCall());
  Assert (s.cols == 0, ExcInvalidConstructorCall());
}



CompressedSimpleSparsityPattern::CompressedSimpleSparsityPattern (const unsigned int m,
								  const unsigned int n,
								  const IndexSet & rowset_
)
		:
                rows(0),
                cols(0),
		rowset(0)
{
  reinit (m,n, rowset_);
}



CompressedSimpleSparsityPattern::CompressedSimpleSparsityPattern (const unsigned int n)
		:
                rows(0),
                cols(0),
		rowset(0)
{
  reinit (n,n);
}



CompressedSimpleSparsityPattern &
CompressedSimpleSparsityPattern::operator = (const CompressedSimpleSparsityPattern &s)
{
  Assert (s.rows == 0, ExcInvalidConstructorCall());
  Assert (s.cols == 0, ExcInvalidConstructorCall());

  Assert (rows == 0, ExcInvalidConstructorCall());
  Assert (cols == 0, ExcInvalidConstructorCall());

  return *this;
}



void
CompressedSimpleSparsityPattern::reinit (const unsigned int m,
					 const unsigned int n,
					 const IndexSet & rowset_)
{
  rows = m;
  cols = n;
  rowset=rowset_;

  Assert(rowset.size()==0 || rowset.size() == m, ExcInvalidConstructorCall());
  
  std::vector<Line> new_lines (rowset.size()==0 ? rows : rowset.n_elements());
  lines.swap (new_lines);
}



void
CompressedSimpleSparsityPattern::compress ()
{}



bool
CompressedSimpleSparsityPattern::empty () const
{
  return ((rows==0) && (cols==0));
}



unsigned int
CompressedSimpleSparsityPattern::max_entries_per_row () const
{
  unsigned int m = 0;
  for (unsigned int i=0; i<lines.size(); ++i)
    {
      m = std::max (m, static_cast<unsigned int>(lines[i].entries.size()));
    }

  return m;
}



bool
CompressedSimpleSparsityPattern::exists (const unsigned int i,
					 const unsigned int j) const
{
  Assert (i<rows, ExcIndexRange(i, 0, rows));
  Assert (j<cols, ExcIndexRange(j, 0, cols));
  Assert( rowset.size()==0 || rowset.is_element(i), ExcInternalError());

  const unsigned int rowindex =
    rowset.size()==0 ? i : rowset.index_within_set(i);

  return std::binary_search (lines[rowindex].entries.begin(),
                             lines[rowindex].entries.end(),
                             j);
}



void
CompressedSimpleSparsityPattern::symmetrize ()
{
  Assert (rows==cols, ExcNotQuadratic());

				   // loop over all elements presently
				   // in the sparsity pattern and add
				   // the transpose element. note:
				   //
				   // 1. that the sparsity pattern
				   // changes which we work on, but
				   // not the present row
				   //
				   // 2. that the @p{add} function can
				   // be called on elements that
				   // already exist without any harm
  for (unsigned int row=0; row<lines.size(); ++row)
    {
      const unsigned int rowindex =
	rowset.size()==0 ? row : rowset.nth_index_in_set(row);

      for (std::vector<unsigned int>::const_iterator
             j=lines[row].entries.begin();
           j != lines[row].entries.end();
           ++j)
				       // add the transpose entry if
				       // this is not the diagonal
        if (rowindex != *j)
          add (*j, rowindex);
    }
}



void
CompressedSimpleSparsityPattern::print (std::ostream &out) const
{
  for (unsigned int row=0; row<lines.size(); ++row)
    {
      out << '[' << (rowset.size()==0 ? row : rowset.nth_index_in_set(row));

      for (std::vector<unsigned int>::const_iterator
             j=lines[row].entries.begin();
           j != lines[row].entries.end(); ++j)
        out << ',' << *j;

      out << ']' << std::endl;
    }

  AssertThrow (out, ExcIO());
}



void
CompressedSimpleSparsityPattern::print_gnuplot (std::ostream &out) const
{
  for (unsigned int row=0; row<lines.size(); ++row)
    {
      const unsigned int rowindex =
	rowset.size()==0 ? row : rowset.nth_index_in_set(row);

      for (std::vector<unsigned int>::const_iterator
             j=lines[row].entries.begin();
           j != lines[row].entries.end(); ++j)
                                         // while matrix entries are usually
                                         // written (i,j), with i vertical and
                                         // j horizontal, gnuplot output is
                                         // x-y, that is we have to exchange
                                         // the order of output
        out << *j << " "
	    << -static_cast<signed int>(rowindex)
	    << std::endl;
    }


  AssertThrow (out, ExcIO());
}



unsigned int
CompressedSimpleSparsityPattern::bandwidth () const
{
  unsigned int b=0;
  for (unsigned int row=0; row<lines.size(); ++row)
    {
      const unsigned int rowindex =
	rowset.size()==0 ? row : rowset.nth_index_in_set(row);

      for (std::vector<unsigned int>::const_iterator
             j=lines[row].entries.begin();
           j != lines[row].entries.end(); ++j)
	if (static_cast<unsigned int>(std::abs(static_cast<int>(rowindex-*j))) > b)
	  b = std::abs(static_cast<signed int>(rowindex-*j));
    }

  return b;
}



unsigned int
CompressedSimpleSparsityPattern::n_nonzero_elements () const
{
  unsigned int n=0;
  for (unsigned int i=0; i<lines.size(); ++i)
    {
      n += lines[i].entries.size();
    }

  return n;
}


unsigned int
CompressedSimpleSparsityPattern::memory_consumption () const
{
				   //TODO: IndexSet...
  unsigned int mem = sizeof(CompressedSimpleSparsityPattern);
  for (unsigned int i=0; i<lines.size(); ++i)
    mem += MemoryConsumption::memory_consumption (lines[i]);

  return mem;
}


// explicit instantiations
template void CompressedSimpleSparsityPattern::Line::add_entries(unsigned int *,
								 unsigned int *,
								 const bool);
template void CompressedSimpleSparsityPattern::Line::add_entries(const unsigned int *,
								 const unsigned int *,
								 const bool);
#ifndef DEAL_II_VECTOR_ITERATOR_IS_POINTER
template void CompressedSimpleSparsityPattern::Line::
add_entries(std::vector<unsigned int>::iterator,
	    std::vector<unsigned int>::iterator,
	    const bool);
#endif

DEAL_II_NAMESPACE_CLOSE
