//----------------------------  compressed_sparsity_pattern.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  compressed_sparsity_pattern.cc  ---------------------------


#include <lac/compressed_sparsity_pattern.h>

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <functional>


const unsigned int CompressedSparsityPattern::Line::cache_size;



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



CompressedSparsityPattern::CompressedSparsityPattern (const unsigned int m,
						      const unsigned int n) 
		:
                rows(0),
                cols(0)
{
  reinit (m,n);
}



CompressedSparsityPattern::CompressedSparsityPattern (const unsigned int n)
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
CompressedSparsityPattern::reinit (const unsigned int m,
				   const unsigned int n)
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



unsigned int
CompressedSparsityPattern::max_entries_per_row () const
{
  unsigned int m = 0;
  for (unsigned int i=0; i<rows; ++i)
    {
      lines[i].flush_cache ();
      m = std::max (m, static_cast<unsigned int>(lines[i].entries.size()));
    }
  
  return m;
}



bool 
CompressedSparsityPattern::exists (const unsigned int i,
				   const unsigned int j) const
{
  Assert (i<rows, ExcInvalidIndex(i,rows));
  Assert (j<cols, ExcInvalidIndex(j,cols));

  lines[i].flush_cache();
  return std::binary_search (lines[i].entries.begin(),
                             lines[i].entries.end(),
                             j);
}



void
CompressedSparsityPattern::symmetrize ()
{
  Assert (rows==cols, ExcNotSquare());

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
      lines[row].flush_cache ();
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
CompressedSparsityPattern::print_gnuplot (std::ostream &out) const
{ 
  for (unsigned int row=0; row<rows; ++row)
    {
      lines[row].flush_cache ();
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
CompressedSparsityPattern::bandwidth () const
{
  unsigned int b=0;
  for (unsigned int row=0; row<rows; ++row)
    {
      lines[row].flush_cache ();

      for (std::vector<unsigned int>::const_iterator
             j=lines[row].entries.begin();
           j != lines[row].entries.end(); ++j)
        if (static_cast<unsigned int>(std::abs(static_cast<int>(row-*j))) > b)
          b = std::abs(static_cast<signed int>(row-*j));
    }
  
  return b;
}



unsigned int
CompressedSparsityPattern::n_nonzero_elements () const
{
  unsigned int n=0;
  for (unsigned int i=0; i<rows; ++i)
    {
      lines[i].flush_cache ();
      n += lines[i].entries.size();
    }
  
  return n;
}
