// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#include <deal.II/base/memory_consumption.h>

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

  const size_type stop_size = entries.size() + n_elements;

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
      size_type col = *my_it;
      std::vector<size_type>::iterator it =
        Utilities::lower_bound(entries.begin(), entries.end(), col);
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
          it = Utilities::lower_bound(it, entries.end(), col);
          if (it == entries.end())
            break;
        }
      // all input entries were duplicates.
      if (my_it == end)
        return;

      // resize vector by just inserting the
      // list
      const size_type pos1 = it - entries.begin();
      Assert (pos1 <= entries.size(), ExcInternalError());
      entries.insert (it, my_it, end);
      it = entries.begin() + pos1;
      Assert (entries.size() >= (size_type)(it-entries.begin()), ExcInternalError());

      // now merge the two lists.
      std::vector<size_type>::iterator it2 = it + (end-my_it);

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
      const size_type new_size = it - entries.begin();
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

  size_type col = *my_it;
  std::vector<size_type>::iterator it, it2;
  // insert the first element as for one
  // entry only first check the last
  // element (or if line is still empty)
  if ( (entries.size()==0) || ( entries.back() < col) )
    {
      entries.push_back(col);
      it = entries.end()-1;
    }
  else
    {
      // do a binary search to find the place
      // where to insert:
      it2 = Utilities::lower_bound(entries.begin(), entries.end(), col);

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
      if (col > entries.back())
        {
          entries.push_back(col);
          it = entries.end()-1;
        }
      // search to the right (preferred search
      // direction)
      else if (col > *it)
        {
          it2 = Utilities::lower_bound(it++, entries.end(), col);
          if (*it2 != col)
            it = entries.insert(it2, col);
        }
      // search to the left
      else if (col < *it)
        {
          it2 = Utilities::lower_bound(entries.begin(), it, col);
          if (*it2 != col)
            it = entries.insert(it2, col);
        }
      // if we're neither larger nor smaller,
      // then this was a duplicate and we can
      // just continue.
    }
}


CompressedSimpleSparsityPattern::size_type
CompressedSimpleSparsityPattern::Line::memory_consumption () const
{
  return entries.capacity()*sizeof(size_type)+sizeof(Line);
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



CompressedSimpleSparsityPattern::CompressedSimpleSparsityPattern (const size_type m,
    const size_type n,
    const IndexSet &rowset_
                                                                 )
  :
  rows(0),
  cols(0),
  rowset(0)
{
  reinit (m,n, rowset_);
}


CompressedSimpleSparsityPattern::CompressedSimpleSparsityPattern (const IndexSet &rowset_)
  :
  rows(0),
  cols(0),
  rowset(0)
{
  reinit (rowset_.size(), rowset_.size(), rowset_);
}


CompressedSimpleSparsityPattern::CompressedSimpleSparsityPattern (const size_type n)
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
CompressedSimpleSparsityPattern::reinit (const size_type m,
                                         const size_type n,
                                         const IndexSet &rowset_)
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



CompressedSimpleSparsityPattern::size_type
CompressedSimpleSparsityPattern::max_entries_per_row () const
{
  size_type m = 0;
  for (size_type i=0; i<lines.size(); ++i)
    {
      m = std::max (m, static_cast<size_type>(lines[i].entries.size()));
    }

  return m;
}



bool
CompressedSimpleSparsityPattern::exists (const size_type i,
                                         const size_type j) const
{
  Assert (i<rows, ExcIndexRange(i, 0, rows));
  Assert (j<cols, ExcIndexRange(j, 0, cols));
  Assert( rowset.size()==0 || rowset.is_element(i), ExcInternalError());

  const size_type rowindex =
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
  for (size_type row=0; row<lines.size(); ++row)
    {
      const size_type rowindex =
        rowset.size()==0 ? row : rowset.nth_index_in_set(row);

      for (std::vector<size_type>::const_iterator
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
  for (size_type row=0; row<lines.size(); ++row)
    {
      out << '[' << (rowset.size()==0 ? row : rowset.nth_index_in_set(row));

      for (std::vector<size_type >::const_iterator
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
  for (size_type row=0; row<lines.size(); ++row)
    {
      const size_type rowindex =
        rowset.size()==0 ? row : rowset.nth_index_in_set(row);

      for (std::vector<size_type >::const_iterator
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



CompressedSimpleSparsityPattern::size_type
CompressedSimpleSparsityPattern::bandwidth () const
{
  size_type b=0;
  for (size_type row=0; row<lines.size(); ++row)
    {
      const size_type rowindex =
        rowset.size()==0 ? row : rowset.nth_index_in_set(row);

      for (std::vector<size_type>::const_iterator
           j=lines[row].entries.begin();
           j != lines[row].entries.end(); ++j)
        if (static_cast<size_type>(std::abs(static_cast<int>(rowindex-*j))) > b)
          b = std::abs(static_cast<signed int>(rowindex-*j));
    }

  return b;
}



CompressedSimpleSparsityPattern::size_type
CompressedSimpleSparsityPattern::n_nonzero_elements () const
{
  size_type n=0;
  for (size_type i=0; i<lines.size(); ++i)
    {
      n += lines[i].entries.size();
    }

  return n;
}


CompressedSimpleSparsityPattern::size_type
CompressedSimpleSparsityPattern::memory_consumption () const
{
  //TODO: IndexSet...
  size_type mem = sizeof(CompressedSimpleSparsityPattern);
  for (size_type i=0; i<lines.size(); ++i)
    mem += MemoryConsumption::memory_consumption (lines[i]);

  return mem;
}


// explicit instantiations
template void CompressedSimpleSparsityPattern::Line::add_entries(size_type *,
    size_type *,
    const bool);
template void CompressedSimpleSparsityPattern::Line::add_entries(const size_type *,
    const size_type *,
    const bool);
#ifndef DEAL_II_VECTOR_ITERATOR_IS_POINTER
template void CompressedSimpleSparsityPattern::Line::
add_entries(std::vector<size_type>::iterator,
            std::vector<size_type>::iterator,
            const bool);
#endif

DEAL_II_NAMESPACE_CLOSE
