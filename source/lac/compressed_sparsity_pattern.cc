// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
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

#include <deal.II/lac/compressed_sparsity_pattern.h>

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <functional>

DEAL_II_NAMESPACE_OPEN


#ifdef DEAL_II_MSVC
__declspec(selectany)  // weak external binding because of multiple link error
#endif
const unsigned int CompressedSparsityPattern::Line::cache_size;



// This function was originally inlined, because it was called very
// often without any need, because cache_entries==0. On the other
// hand, it is very long and causes linker warnings on certain
// systems.

// Therefore, we require now, that the caller checks if this function
// is necessary and only calls it, if it is actually used. Since it is
// called less often, it is removed from the inlined section.
void
CompressedSparsityPattern::Line::flush_cache () const
{
  // Make sure the caller checked
  // necessity of this function.
  Assert(cache_entries != 0, ExcInternalError());

  // first sort the entries in the cache, so
  // that it is easier to merge it with the
  // main array. note that due to the way
  // add() inserts elements, there can be no
  // duplicates in the cache
  //
  // do the sorting in a way that is fast for
  // the small cache sizes we have
  // here. basically, use bubble sort
  switch (cache_entries)
    {
    case 1:
    {
      break;
    }

    case 2:
    {
      if (cache[1] < cache[0])
        std::swap (cache[0], cache[1]);
      break;
    }

    case 3:
    {
      if (cache[1] < cache[0])
        std::swap (cache[0], cache[1]);
      if (cache[2] < cache[1])
        std::swap (cache[1], cache[2]);
      if (cache[1] < cache[0])
        std::swap (cache[0], cache[1]);
      break;
    }

    case 4:
    case 5:
    case 6:
    case 7:
    {
      for (unsigned int i=0; i<cache_entries; ++i)
        for (unsigned int j=i+1; j<cache_entries; ++j)
          if (cache[j] < cache[i])
            std::swap (cache[i], cache[j]);
      break;
    }

    default:
    {
      std::sort (&cache[0], &cache[cache_entries]);
      break;
    }
    }

  // TODO: could use the add_entries
  // function of the constraint line for
  // doing this, but that one is
  // non-const. Still need to figure out
  // how to do that.

  // next job is to merge the two
  // arrays. special case the case that the
  // original array is empty.
  if (entries.size() == 0)
    {
      entries.resize (cache_entries);
      for (unsigned int i=0; i<cache_entries; ++i)
        entries[i] = cache[i];
    }
  else
    {
      // first count how many of the cache
      // entries are already in the main
      // array, so that we can efficiently
      // allocate memory
      unsigned int n_new_entries = 0;
      {
        unsigned int cache_position = 0;
        unsigned int entry_position = 0;
        while ((entry_position<entries.size()) &&
               (cache_position<cache_entries))
          {
            ++n_new_entries;
            if (entries[entry_position] < cache[cache_position])
              ++entry_position;
            else if (entries[entry_position] == cache[cache_position])
              {
                ++entry_position;
                ++cache_position;
              }
            else
              ++cache_position;
          }

        // scoop up leftovers in arrays
        n_new_entries += (entries.size() - entry_position) +
                         (cache_entries - cache_position);
      }

      // then allocate new memory and merge
      // arrays, if there are any entries at
      // all that need to be merged
      Assert (n_new_entries >= entries.size(),
              ExcInternalError());
      if (n_new_entries > entries.size())
        {
          std::vector<types::global_dof_index> new_entries;
          new_entries.reserve (n_new_entries);
          unsigned int cache_position = 0;
          unsigned int entry_position = 0;
          while ((entry_position<entries.size()) &&
                 (cache_position<cache_entries))
            if (entries[entry_position] < cache[cache_position])
              {
                new_entries.push_back (entries[entry_position]);
                ++entry_position;
              }
            else if (entries[entry_position] == cache[cache_position])
              {
                new_entries.push_back (entries[entry_position]);
                ++entry_position;
                ++cache_position;
              }
            else
              {
                new_entries.push_back (cache[cache_position]);
                ++cache_position;
              }

          // copy remaining elements from the
          // array that we haven't
          // finished. note that at most one
          // of the following loops will run
          // at all
          for (; entry_position < entries.size(); ++entry_position)
            new_entries.push_back (entries[entry_position]);
          for (; cache_position < cache_entries; ++cache_position)
            new_entries.push_back (cache[cache_position]);

          Assert (new_entries.size() == n_new_entries,
                  ExcInternalError());

          // finally swap old and new array,
          // and set cache size to zero
          new_entries.swap (entries);
        }
    }

  cache_entries = 0;
}



template <typename ForwardIterator>
void
CompressedSparsityPattern::Line::add_entries (ForwardIterator begin,
                                              ForwardIterator end,
                                              const bool indices_are_sorted)
{
  // use the same code as when flushing the
  // cache in case we have many (more than
  // three) entries in a sorted
  // list. Otherwise, go on to the single
  // add() function.
  const int n_elements = end - begin;
  if (n_elements <= 0)
    return;

  const unsigned int n_cols = static_cast<unsigned int>(n_elements);

  if (indices_are_sorted == true)
    {

      // next job is to merge the two
      // arrays. special case the case that the
      // original array is empty.
      if (entries.size() == 0)
        {
          entries.resize (n_cols);
          ForwardIterator my_it = begin;
          for (unsigned int i=0; i<n_cols; ++i)
            entries[i] = *my_it++;
        }
      else
        {
          // first count how many of the cache
          // entries are already in the main
          // array, so that we can efficiently
          // allocate memory
          unsigned int n_new_entries = 0;
          {
            unsigned int entry_position = 0;
            ForwardIterator my_it = begin;
            while ((entry_position<entries.size()) &&
                   (my_it != end))
              {
                ++n_new_entries;
                if (entries[entry_position] < *my_it)
                  ++entry_position;
                else if (entries[entry_position] == *my_it)
                  {
                    ++entry_position;
                    ++my_it;
                  }
                else
                  ++my_it;
              }

            // scoop up leftovers in arrays
            n_new_entries += (entries.size() - entry_position) +
                             (end - my_it);
          }

          // then allocate new memory and merge
          // arrays, if there are any entries at
          // all that need to be merged
          Assert (n_new_entries >= entries.size(),
                  ExcInternalError());
          if (n_new_entries > entries.size())
            {
              std::vector<types::global_dof_index> new_entries;
              new_entries.reserve (n_new_entries);
              ForwardIterator my_it = begin;
              unsigned int entry_position = 0;
              while ((entry_position<entries.size()) &&
                     (my_it != end))
                if (entries[entry_position] < *my_it)
                  {
                    new_entries.push_back (entries[entry_position]);
                    ++entry_position;
                  }
                else if (entries[entry_position] == *my_it)
                  {
                    new_entries.push_back (entries[entry_position]);
                    ++entry_position;
                    ++my_it;
                  }
                else
                  {
                    new_entries.push_back (*my_it);
                    ++my_it;
                  }

              // copy remaining elements from the
              // array that we haven't
              // finished. note that at most one
              // of the following loops will run
              // at all
              for (; entry_position < entries.size(); ++entry_position)
                new_entries.push_back (entries[entry_position]);
              for (; my_it != end; ++my_it)
                new_entries.push_back (*my_it);

              Assert (new_entries.size() == n_new_entries,
                      ExcInternalError());

              // finally swap old and new array,
              // and set cache size to zero
              new_entries.swap (entries);
            }
        }
      return;
    }

  // otherwise, insert the indices one
  // after each other
  for (ForwardIterator it = begin; it != end; ++it)
    add (*it);
}



CompressedSparsityPattern::CompressedSparsityPattern ()
  :
  rows(0),
  cols(0)
{}



CompressedSparsityPattern::
CompressedSparsityPattern (const CompressedSparsityPattern &s)
  :
  Subscriptor(),
  rows(0),
  cols(0)
{
  Assert (s.rows == 0, ExcInvalidConstructorCall());
  Assert (s.cols == 0, ExcInvalidConstructorCall());
}



CompressedSparsityPattern::CompressedSparsityPattern (const size_type m,
                                                      const size_type n)
  :
  rows(0),
  cols(0)
{
  reinit (m,n);
}



CompressedSparsityPattern::CompressedSparsityPattern (const size_type n)
  :
  rows(0),
  cols(0)
{
  reinit (n,n);
}



CompressedSparsityPattern &
CompressedSparsityPattern::operator = (const CompressedSparsityPattern &s)
{
  Assert (s.rows == 0, ExcInvalidConstructorCall());
  Assert (s.cols == 0, ExcInvalidConstructorCall());

  Assert (rows == 0, ExcInvalidConstructorCall());
  Assert (cols == 0, ExcInvalidConstructorCall());

  return *this;
}



void
CompressedSparsityPattern::reinit (const size_type m,
                                   const size_type n)
{
  rows = m;
  cols = n;

  std::vector<Line> new_lines (rows);
  lines.swap (new_lines);
}



void
CompressedSparsityPattern::compress ()
{}



bool
CompressedSparsityPattern::empty () const
{
  return ((rows==0) && (cols==0));
}



CompressedSparsityPattern::size_type
CompressedSparsityPattern::max_entries_per_row () const
{
  size_type m = 0;
  for (size_type i=0; i<rows; ++i)
    {
      if (lines[i].cache_entries != 0)
        lines[i].flush_cache ();
      m = std::max (m, static_cast<size_type>(lines[i].entries.size()));
    }

  return m;
}



bool
CompressedSparsityPattern::exists (const size_type i,
                                   const size_type j) const
{
  Assert (i<rows, ExcIndexRange(i, 0, rows));
  Assert (j<cols, ExcIndexRange(j, 0, cols));

  if (lines[i].cache_entries != 0)
    lines[i].flush_cache();
  return std::binary_search (lines[i].entries.begin(),
                             lines[i].entries.end(),
                             j);
}



void
CompressedSparsityPattern::symmetrize ()
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
  for (size_type row=0; row<rows; ++row)
    {
      if (lines[row].cache_entries != 0)
        lines[row].flush_cache ();
      for (std::vector<size_type>::const_iterator
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
CompressedSparsityPattern::print (std::ostream &out) const
{
  for (size_type row=0; row<rows; ++row)
    {
      if (lines[row].cache_entries != 0)
        lines[row].flush_cache ();

      out << '[' << row;

      for (std::vector<size_type>::const_iterator
           j=lines[row].entries.begin();
           j != lines[row].entries.end(); ++j)
        out << ',' << *j;

      out << ']' << std::endl;
    }

  AssertThrow (out, ExcIO());
}



void
CompressedSparsityPattern::print_gnuplot (std::ostream &out) const
{
  for (size_type row=0; row<rows; ++row)
    {
      if (lines[row].cache_entries != 0)
        lines[row].flush_cache ();
      for (std::vector<size_type>::const_iterator
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



CompressedSparsityPattern::size_type
CompressedSparsityPattern::bandwidth () const
{
  size_type b=0;
  for (size_type row=0; row<rows; ++row)
    {
      if (lines[row].cache_entries != 0)
        lines[row].flush_cache ();

      for (std::vector<size_type>::const_iterator
           j=lines[row].entries.begin();
           j != lines[row].entries.end(); ++j)
        if (static_cast<size_type>(std::abs(static_cast<int>(row-*j))) > b)
          b = std::abs(static_cast<signed int>(row-*j));
    }

  return b;
}



CompressedSparsityPattern::size_type
CompressedSparsityPattern::n_nonzero_elements () const
{
  size_type n=0;
  for (size_type i=0; i<rows; ++i)
    {
      if (lines[i].cache_entries != 0)
        lines[i].flush_cache ();
      n += lines[i].entries.size();
    }

  return n;
}


// explicit instantiations
template void CompressedSparsityPattern::Line::add_entries(size_type *,
                                                           size_type *,
                                                           const bool);
template void CompressedSparsityPattern::Line::add_entries(const size_type *,
                                                           const size_type *,
                                                           const bool);
#ifndef DEAL_II_VECTOR_ITERATOR_IS_POINTER
template void CompressedSparsityPattern::Line::
add_entries(std::vector<size_type>::iterator,
            std::vector<size_type>::iterator,
            const bool);
#endif

DEAL_II_NAMESPACE_CLOSE
