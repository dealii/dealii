//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/compressed_simple_sparsity_pattern.h>

#include <iostream>
#include <iomanip>
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

  const unsigned int n_cols = static_cast<unsigned int>(n_elements);
  const unsigned int stop_size = entries.size() + n_cols;

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

				   // resize vector by just inserting the
				   // list
      const unsigned int col = *begin;
      std::vector<unsigned int>::iterator it = 
	std::lower_bound(entries.begin(), entries.end(), col);
      const unsigned int pos1 = it - entries.begin();
      entries.insert (it, begin, end);
      it = entries.begin() + pos1;

				   // now merge the two lists.
      ForwardIterator my_it = begin;
      std::vector<unsigned int>::iterator it2 = it + n_cols;

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



CompressedSimpleSparsityPattern::CompressedSimpleSparsityPattern ()
                :
		rows(0),
		cols(0)
{}



CompressedSimpleSparsityPattern::
CompressedSimpleSparsityPattern (const CompressedSimpleSparsityPattern &s)
                :
		Subscriptor(),
		rows(0),
		cols(0)
{
  Assert (s.rows == 0, ExcInvalidConstructorCall());
  Assert (s.cols == 0, ExcInvalidConstructorCall());
}



CompressedSimpleSparsityPattern::CompressedSimpleSparsityPattern (const unsigned int m,
						      const unsigned int n) 
		:
                rows(0),
                cols(0)
{
  reinit (m,n);
}



CompressedSimpleSparsityPattern::CompressedSimpleSparsityPattern (const unsigned int n)
		:
                rows(0),
                cols(0)
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
				   const unsigned int n)
{
  rows = m;
  cols = n;

  std::vector<Line> new_lines (rows);
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
  for (unsigned int i=0; i<rows; ++i)
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

  return std::binary_search (lines[i].entries.begin(),
                             lines[i].entries.end(),
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
  for (unsigned int row=0; row<rows; ++row)
    {
      for (std::vector<unsigned int>::const_iterator
             j=lines[row].entries.begin();
           j != lines[row].entries.end();
           ++j)
				       // add the transpose entry if
				       // this is not the diagonal
        if (row != *j)
          add (*j, row);
    }
}



void
CompressedSimpleSparsityPattern::print (std::ostream &out) const
{
  for (unsigned int row=0; row<rows; ++row)
    {
      out << '[' << row;

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
  for (unsigned int row=0; row<rows; ++row)
    {
      for (std::vector<unsigned int>::const_iterator
             j=lines[row].entries.begin();
           j != lines[row].entries.end(); ++j)
                                         // while matrix entries are usually
                                         // written (i,j), with i vertical and
                                         // j horizontal, gnuplot output is
                                         // x-y, that is we have to exchange
                                         // the order of output
        out << *j << " " << -static_cast<signed int>(row) << std::endl;
    }
      

  AssertThrow (out, ExcIO());
}



unsigned int
CompressedSimpleSparsityPattern::bandwidth () const
{
  unsigned int b=0;
  for (unsigned int row=0; row<rows; ++row)
    {
      for (std::vector<unsigned int>::const_iterator
             j=lines[row].entries.begin();
           j != lines[row].entries.end(); ++j)
        if (static_cast<unsigned int>(std::abs(static_cast<int>(row-*j))) > b)
          b = std::abs(static_cast<signed int>(row-*j));
    }
  
  return b;
}



unsigned int
CompressedSimpleSparsityPattern::n_nonzero_elements () const
{
  unsigned int n=0;
  for (unsigned int i=0; i<rows; ++i)
    {
      n += lines[i].entries.size();
    }
  
  return n;
}


// explicit instantiations
template void CompressedSimpleSparsityPattern::Line::add_entries(unsigned int *,
								 unsigned int *,
								 const bool);
template void CompressedSimpleSparsityPattern::Line::
add_entries(std::vector<unsigned int>::iterator,
	    std::vector<unsigned int>::iterator,
	    const bool);


DEAL_II_NAMESPACE_CLOSE
