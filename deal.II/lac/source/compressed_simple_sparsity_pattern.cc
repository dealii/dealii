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

DEAL_II_NAMESPACE_CLOSE
