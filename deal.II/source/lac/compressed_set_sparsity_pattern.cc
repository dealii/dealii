//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2005, 2006, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/compressed_set_sparsity_pattern.h>

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <functional>

DEAL_II_NAMESPACE_OPEN



CompressedSetSparsityPattern::CompressedSetSparsityPattern ()
                :
		rows(0),
		cols(0)
{}



CompressedSetSparsityPattern::
CompressedSetSparsityPattern (const CompressedSetSparsityPattern &s)
                :
		Subscriptor(),
		rows(0),
		cols(0)
{
  Assert (s.rows == 0, ExcInvalidConstructorCall());
  Assert (s.cols == 0, ExcInvalidConstructorCall());
}



CompressedSetSparsityPattern::CompressedSetSparsityPattern (const unsigned int m,
						      const unsigned int n) 
		:
                rows(0),
                cols(0)
{
  reinit (m,n);
}



CompressedSetSparsityPattern::CompressedSetSparsityPattern (const unsigned int n)
		:
                rows(0),
                cols(0)
{
  reinit (n,n);
}



CompressedSetSparsityPattern &
CompressedSetSparsityPattern::operator = (const CompressedSetSparsityPattern &s)
{
  Assert (s.rows == 0, ExcInvalidConstructorCall());
  Assert (s.cols == 0, ExcInvalidConstructorCall());

  Assert (rows == 0, ExcInvalidConstructorCall());
  Assert (cols == 0, ExcInvalidConstructorCall());

  return *this;
}



void
CompressedSetSparsityPattern::reinit (const unsigned int m,
				   const unsigned int n)
{
  rows = m;
  cols = n;

  std::vector<Line> new_lines (rows);
  lines.swap (new_lines);
}



void
CompressedSetSparsityPattern::compress ()
{}



bool
CompressedSetSparsityPattern::empty () const
{
  return ((rows==0) && (cols==0));
}



unsigned int
CompressedSetSparsityPattern::max_entries_per_row () const
{
  unsigned int m = 0;
  for (unsigned int i=0; i<rows; ++i)
    {
      m = std::max (m, static_cast<unsigned int>(lines[i].entries.size()));
    }
  
  return m;
}



bool 
CompressedSetSparsityPattern::exists (const unsigned int i,
				   const unsigned int j) const
{
  Assert (i<rows, ExcIndexRange(i, 0, rows));
  Assert (j<cols, ExcIndexRange(j, 0, cols));

  return (lines[i].entries.find (j) != lines[i].entries.end ());
}



void
CompressedSetSparsityPattern::symmetrize ()
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
      for (std::set<unsigned int>::const_iterator
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
CompressedSetSparsityPattern::print (std::ostream &out) const
{ 
  AssertThrow (out, ExcIO());

  for (unsigned int row=0; row<rows; ++row)
    {
      out << '[' << row;
      
      for (std::set<unsigned int>::const_iterator
             j=lines[row].entries.begin();
           j != lines[row].entries.end(); ++j)
        out << ',' << *j;

      out << ']' << std::endl;
    }

  AssertThrow (out, ExcIO());
}



void
CompressedSetSparsityPattern::print_gnuplot (std::ostream &out) const
{ 
  AssertThrow (out, ExcIO());

  for (unsigned int row=0; row<rows; ++row)
    {
      for (std::set<unsigned int>::const_iterator
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
CompressedSetSparsityPattern::bandwidth () const
{
  unsigned int b=0;
  for (unsigned int row=0; row<rows; ++row)
    {
      for (std::set<unsigned int>::const_iterator
             j=lines[row].entries.begin();
           j != lines[row].entries.end(); ++j)
        if (static_cast<unsigned int>(std::abs(static_cast<int>(row-*j))) > b)
          b = std::abs(static_cast<signed int>(row-*j));
    }
  
  return b;
}



unsigned int
CompressedSetSparsityPattern::n_nonzero_elements () const
{
  unsigned int n=0;
  for (unsigned int i=0; i<rows; ++i)
    {
      n += lines[i].entries.size();
    }
  
  return n;
}

DEAL_II_NAMESPACE_CLOSE
