//----------------------------  compressed_sparsity_pattern.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002 by the deal.II authors
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

CompressedSparsityPattern::CompressedSparsityPattern () :
		rows(0),
		cols(0)
{};



CompressedSparsityPattern::CompressedSparsityPattern (const CompressedSparsityPattern &s) :
		Subscriptor(),
		rows(0),
		cols(0)
{
  Assert (s.rows == 0, ExcInvalidConstructorCall());
  Assert (s.cols == 0, ExcInvalidConstructorCall());
};



CompressedSparsityPattern::CompressedSparsityPattern (const unsigned int m,
						      const unsigned int n) 
		: rows(0),
		  cols(0)
{
  reinit (m,n);
};



CompressedSparsityPattern::CompressedSparsityPattern (const unsigned int n)
		: rows(0),
		  cols(0)
{
  reinit (n,n);
};



CompressedSparsityPattern &
CompressedSparsityPattern::operator = (const CompressedSparsityPattern &s)
{
  Assert (s.rows == 0, ExcInvalidConstructorCall());
  Assert (s.cols == 0, ExcInvalidConstructorCall());

  Assert (rows == 0, ExcInvalidConstructorCall());
  Assert (cols == 0, ExcInvalidConstructorCall());

  return *this;
};



void
CompressedSparsityPattern::reinit (const unsigned int m,
				   const unsigned int n)
{
  rows = m;
  cols = n;

  std::vector<std::set<unsigned int> > new_column_indices (rows);
  column_indices.swap (new_column_indices);
}



void
CompressedSparsityPattern::compress ()
{};



bool
CompressedSparsityPattern::empty () const
{
  return ((rows==0) && (cols==0));
};



unsigned int
CompressedSparsityPattern::max_entries_per_row () const 
{
  unsigned int m = 0;
  for (unsigned int i=1; i<rows; ++i)
    m = std::max (m, static_cast<unsigned int>(column_indices[i].size()));

  return m;
};



void
CompressedSparsityPattern::add (const unsigned int i,
				const unsigned int j)
{
  Assert (i<rows, ExcInvalidIndex(i,rows));
  Assert (j<cols, ExcInvalidIndex(j,cols));

				   // the std::set automatically
				   // assures uniqueness and
				   // sortedness of the column indices
  column_indices[i].insert (j);
};


bool 
CompressedSparsityPattern::exists (const unsigned int i,
				   const unsigned int j) const
{
  Assert (i<rows, ExcInvalidIndex(i,rows));
  Assert (j<cols, ExcInvalidIndex(j,cols));
//TODO: Ask Wolfgang on how to use column_indices
  Assert(false, ExcNotImplemented());
  return true;
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
    for (std::set<unsigned int>::const_iterator i=column_indices[row].begin();
	 i!=column_indices[row].end(); ++i)
				       // add the transpose entry if
				       // this is not the diagonal
      if (row != *i)
	add (*i, row);
};



void
CompressedSparsityPattern::print_gnuplot (std::ostream &out) const
{
  for (unsigned int row=0; row<rows; ++row)
    for (std::set<unsigned int>::const_iterator i=column_indices[row].begin();
	 i!=column_indices[row].end(); ++i)
				       // while matrix entries are
				       // usually written (i,j),
				       // with i vertical and j
				       // horizontal, gnuplot output
				       // is x-y, that is we have to
				       // exchange the order of
				       // output
      out << *i << " " << -static_cast<signed int>(row) << std::endl;

  AssertThrow (out, ExcIO());
}



unsigned int
CompressedSparsityPattern::row_length (const unsigned int row) const
{
  return column_indices[row].size();
};



unsigned int
CompressedSparsityPattern::column_number (const unsigned int row,
					  const unsigned int index) const
{
  Assert (index < column_indices[row].size(),
	  ExcIndexRange (index, 0, column_indices[row].size()));
  std::set<unsigned int>::const_iterator p = column_indices[row].begin();
  std::advance (p, index);
  return *p;
};



unsigned int
CompressedSparsityPattern::bandwidth () const
{
  unsigned int b=0;
  for (unsigned int row=0; row<rows; ++row)
    for (std::set<unsigned int>::const_iterator i=column_indices[row].begin();
	 i!=column_indices[row].end(); ++i)
      if (static_cast<unsigned int>(std::abs(static_cast<int>(row-*i))) > b)
	b = std::abs(static_cast<signed int>(row-*i));

  return b;
};



unsigned int
CompressedSparsityPattern::n_nonzero_elements () const
{
  unsigned int n=0;
  for (unsigned int i=0; i<rows; ++i)
    n += column_indices[i].size();
  return n;
};
