//----------------------------  sparsity_pattern.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparsity_pattern.cc  ---------------------------


#include <lac/sparse_matrix.h>

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <numeric>


SparsityPattern::SparsityPattern () :
		max_dim(0),
		max_vec_len(0),
		rowstart(0),
		colnums(0)
{
  reinit (0,0,0);
};



SparsityPattern::SparsityPattern (const SparsityPattern &s) :
		Subscriptor(),
		max_dim(0),
		max_vec_len(0),
		rowstart(0),
		colnums(0)
{
  Assert (s.rowstart == 0, ExcInvalidConstructorCall());
  Assert (s.colnums == 0, ExcInvalidConstructorCall());
  Assert (s.rows == 0, ExcInvalidConstructorCall());
  Assert (s.cols == 0, ExcInvalidConstructorCall());
  
  reinit (0,0,0);
};



SparsityPattern::SparsityPattern (const unsigned int m,
				  const unsigned int n,
				  const unsigned int max_per_row) 
		: max_dim(0),
		  max_vec_len(0),
		  rowstart(0),
		  colnums(0)
{
  reinit (m,n,max_per_row);
};



SparsityPattern::SparsityPattern (const unsigned int          m,
				  const unsigned int          n,
				  const vector<unsigned int> &row_lengths) 
		: max_dim(0),
		  max_vec_len(0),
		  rowstart(0),
		  colnums(0)
{
  reinit (m, n, row_lengths);
};



SparsityPattern::SparsityPattern (const unsigned int n,
				  const unsigned int max_per_row)
		: max_dim(0),
		  max_vec_len(0),
		  rowstart(0),
		  colnums(0)
{
  reinit (n,n,max_per_row);
};



SparsityPattern::SparsityPattern (const unsigned int          m,
				  const vector<unsigned int> &row_lengths) 
		: max_dim(0),
		  max_vec_len(0),
		  rowstart(0),
		  colnums(0)
{
  reinit (m, m, row_lengths);
};



SparsityPattern::SparsityPattern (const SparsityPattern &original,
				  const unsigned int        max_per_row,
				  const unsigned int        extra_off_diagonals)
		: max_dim(0),
		  max_vec_len(0),
		  rowstart(0),
		  colnums(0)
{
  Assert (original.rows==original.cols, ExcNotSquare());
  Assert (original.is_compressed(), ExcNotCompressed());
  
  reinit (original.rows, original.cols, max_per_row);

				   // now copy the entries from
				   // the other object
  for (unsigned int row=0; row<original.rows; ++row)
    {
				       // copy the elements of this row
				       // of the other object 
				       // 
				       // note that the first object actually
				       // is the main-diagonal element,
				       // which we need not copy
				       //
				       // we do the copying in two steps:
				       // first we note that the elements in
				       // #original# are sorted, so we may
				       // first copy all the elements up to
				       // the first side-diagonal one which
				       // is to be filled in. then we insert
				       // the side-diagonals, finally copy
				       // the rest from that element onwards
				       // which is not a side-diagonal any
				       // more.
      const unsigned int * const
	original_row_start = &original.colnums[original.rowstart[row]] + 1;
				       // the following requires that
				       // #original# be compressed since
				       // otherwise there might be invalid_entry's
      const unsigned int * const
	original_row_end   = &original.colnums[original.rowstart[row+1]];

				       // find pointers before and
				       // after extra
				       // off-diagonals. if at top or
				       // bottom of matrix, then set
				       // these pointers such that no
				       // copying is necessary (see
				       // the #copy# commands)
      const unsigned int * const
	original_last_before_side_diagonals
	= (row > extra_off_diagonals ?
	   lower_bound (original_row_start,
			original_row_end,
			row-extra_off_diagonals) :
	   original_row_start);
      
      const unsigned int * const
	original_first_after_side_diagonals
	= (row < rows-extra_off_diagonals-1 ?
	   upper_bound (original_row_start,
			original_row_end,
			row+extra_off_diagonals) :
	   original_row_end);

				       // find first free slot. the
				       // first slot in each row is
				       // the diagonal element
      unsigned int * next_free_slot = &colnums[rowstart[row]] + 1;

				       // copy elements before side-diagonals
      next_free_slot = copy (original_row_start,
			     original_last_before_side_diagonals,
			     next_free_slot);

				       // insert left and right side-diagonals
      for (unsigned int i=1; i<=min(row,extra_off_diagonals);
	   ++i, ++next_free_slot)
	*next_free_slot = row-i;
      for (unsigned int i=1; i<=min(extra_off_diagonals, rows-row-1);
	   ++i, ++next_free_slot)
	*next_free_slot = row+i;

				       // copy rest
      next_free_slot = copy (original_first_after_side_diagonals,
			     original_row_end,
			     next_free_slot);

				       // this error may happen if the
				       // sum of previous elements per row
				       // and those of the new diagonals
				       // exceeds the maximum number of
				       // elements per row given to this
				       // constructor
      Assert (next_free_slot <= &colnums[rowstart[row+1]],
	      ExcNotEnoughSpace (0,rowstart[row+1]-rowstart[row]));
    };
};



SparsityPattern::~SparsityPattern ()
{
  if (rowstart != 0)  delete[] rowstart;
  if (colnums != 0)   delete[] colnums;
}



SparsityPattern &
SparsityPattern::operator = (const SparsityPattern &s)
{
  Assert (s.rowstart == 0, ExcInvalidConstructorCall());
  Assert (s.colnums == 0, ExcInvalidConstructorCall());
  Assert (s.rows == 0, ExcInvalidConstructorCall());
  Assert (s.cols == 0, ExcInvalidConstructorCall());

  Assert (rowstart == 0, ExcInvalidConstructorCall());
  Assert (colnums == 0, ExcInvalidConstructorCall());
  Assert (rows == 0, ExcInvalidConstructorCall());
  Assert (cols == 0, ExcInvalidConstructorCall());

  return *this;
};



void
SparsityPattern::reinit (const unsigned int m,
			 const unsigned int n,
			 const unsigned int max_per_row)
{
				   // simply map this function to the
				   // other #reinit# function
  const vector<unsigned int> row_lengths (m, max_per_row);
  reinit (m, n, row_lengths);
};



void
SparsityPattern::reinit (const unsigned int m,
			 const unsigned int n,
			 const vector<unsigned int> &row_lengths)
{
  Assert (((m==0) && (n==0)) || (*max_element(row_lengths.begin(), row_lengths.end()) > 0),
	  ExcInvalidNumber(*max_element(row_lengths.begin(), row_lengths.end())));
  Assert (row_lengths.size() == m, ExcInvalidNumber (m));
	  
  rows = m;
  cols = n;

				   // delete empty matrices
  if ((m==0) || (n==0))
    {
      if (rowstart)  delete[] rowstart;
      if (colnums)   delete[] colnums;
      rowstart = 0;
      colnums = 0;
      max_vec_len = max_dim = rows = cols = 0;
				       // if dimension is zero: ignore
				       // max_per_row
      max_row_length = 0;
      compressed = false;
      return;
    };

				   // find out how many entries we
				   // need in the #colnums# array. if
				   // this number is larger than
				   // #max_vec_len#, then we will need
				   // to reallocate memory
				   //
				   // note that the number of elements
				   // is bounded by the number of
				   // columns
  unsigned int vec_len = 0;
  for (unsigned int i=0; i<m; ++i)
    vec_len += min(row_lengths[i], n);

  max_row_length = (row_lengths.size() == 0 ?
		    0 :
		    min (*max_element(row_lengths.begin(), row_lengths.end()),
			 n));


				   // allocate memory for the rowstart
				   // values, if necessary
  if (rows > max_dim)
    {
      if (rowstart) delete[] rowstart;
      max_dim = rows;
      rowstart = new unsigned int[max_dim+1];
    };

				   // allocate memory for the column
				   // numbers if necessary
  if (vec_len > max_vec_len)
    {
      if (colnums) delete[] colnums;
      max_vec_len = vec_len;
      colnums = new unsigned int[max_vec_len];
    };

				   // set the rowstart array 
  rowstart[0] = 0;
  for (unsigned int i=1; i<=rows; ++i)
    rowstart[i] = rowstart[i-1]+min(row_lengths[i-1],n);
  Assert (rowstart[rows]==vec_len, ExcInternalError());

				   // preset the column numbers by a
				   // value indicating it is not in
				   // use
  fill_n (&colnums[0], vec_len, invalid_entry);

				   // if the matrix is square: let the
				   // first entry in each row be the
				   // diagonal value
  if (rows == cols)
    for (unsigned int i=0;i<rows;i++)
      colnums[rowstart[i]] = i;

  compressed = false;
}



void
SparsityPattern::compress ()
{
  Assert ((rowstart!=0) && (colnums!=0), ExcEmptyObject());

				   // do nothing if already compressed
  if (compressed)
    return;
  
  unsigned int next_free_entry = 0,
	       next_row_start  = 0,
	       row_length      = 0;

				   // first find out how many non-zero
				   // elements there are, in order to
				   // allocate the right amount of
				   // memory
  const unsigned int
    nonzero_elements = count_if (&colnums[rowstart[0]],
				 &colnums[rowstart[rows]],
				 bind2nd(not_equal_to<unsigned int>(), invalid_entry));
				   // now allocate the respective memory
  unsigned int *new_colnums = new unsigned int[nonzero_elements];


				   // reserve temporary storage to
				   // store the entries of one row
  vector<unsigned int> tmp_entries (max_row_length);
  
				   // Traverse all rows
  for (unsigned int line=0; line<rows; ++line)
    {
				       // copy used entries, break if
				       // first unused entry is reached
      row_length = 0;
      for (unsigned int j=rowstart[line]; j<rowstart[line+1]; ++j,++row_length)
	if (colnums[j] != invalid_entry)
	  tmp_entries[row_length] = colnums[j];
	else
	  break;
				       // now #rowstart# is
				       // the number of entries in
				       // this line

				       // for square matrices, the
				       // first entry in each row
				       // is the diagonal one. In
				       // this case only sort the
				       // remaining entries, otherwise
				       // sort all
      sort ((rows==cols) ? &tmp_entries[1] : &tmp_entries[0],
	    &tmp_entries[row_length]);

				       // insert column numbers
				       // into the new field
      for (unsigned int j=0; j<row_length; ++j)
	new_colnums[next_free_entry++] = tmp_entries[j];

				       // note new start of this and
				       // the next row
      rowstart[line] = next_row_start;
      next_row_start = next_free_entry;

				       // some internal checks
      Assert ((rows!=cols) ||
	      (new_colnums[rowstart[line]] == line),
	      ExcInternalError());
				       // assert that the first entry
				       // does not show up in
				       // the remaining ones and that
				       // the remaining ones are unique
				       // among themselves (this handles
				       // both cases, quadratic and
				       // rectangular matrices)
				       //
				       // the only exception here is
				       // if the row contains no
				       // entries at all
      Assert ((rowstart[line] == next_row_start)
	      ||
	      (find (&new_colnums[rowstart[line]+1],
		     &new_colnums[next_row_start],
		     new_colnums[rowstart[line]]) ==
	       &new_colnums[next_row_start]),
	      ExcInternalError());
      Assert ((rowstart[line] == next_row_start)
	      ||
	      (adjacent_find(&new_colnums[rowstart[line]+1],
			     &new_colnums[next_row_start]) ==
	       &new_colnums[next_row_start]),
	      ExcInternalError());
    };

				   // assert that we have used all
				   // allocated space, no more and no
				   // less
  Assert (next_free_entry == nonzero_elements,
	  ExcInternalError());

				   // set iterator-past-the-end
  rowstart[rows] = next_row_start;

				   // set colnums to the newly
				   // allocated array and delete the
				   // old one
  delete[] colnums;
  colnums = new_colnums;

				   // store the size
  max_vec_len = nonzero_elements;
  
  compressed = true;
};



bool
SparsityPattern::empty () const
{
				   // let's try to be on the safe side of
				   // life by using multiple possibilities in
				   // the check for emptiness... (sorry for
				   // this kludge -- emptying matrices and
				   // freeing memory was not present in the
				   // original implementation and I don't
				   // know at how many places I missed
				   // something in adding it, so I try to
				   // be cautious. wb)
  if ((rowstart==0) || (rows==0) || (cols==0))
    {
      Assert (rowstart==0, ExcInternalError());
      Assert (rows==0, ExcInternalError());
      Assert (cols==0, ExcInternalError());
      Assert (colnums==0, ExcInternalError());
      Assert (max_vec_len==0, ExcInternalError());

      return true;
    };
  return false;
};



unsigned int
SparsityPattern::max_entries_per_row () const 
{
				   // if compress() has not yet been
				   // called, we can get the maximum
				   // number of elements per row using
				   // the stored value
  if (!compressed)
    return max_row_length;

				   // if compress() was called, we
				   // use a better algorithm which
				   // gives us a sharp bound
  unsigned int m = 0;
  for (unsigned int i=1; i<rows; ++i)
    m = max (m, rowstart[i]-rowstart[i-1]);

  return m;
};



unsigned int
SparsityPattern::operator () (const unsigned int i,
			      const unsigned int j) const
{
  Assert ((rowstart!=0) && (colnums!=0), ExcEmptyObject());  
  Assert (i<rows, ExcInvalidIndex(i,rows));
  Assert (j<cols, ExcInvalidIndex(j,cols));
  Assert (compressed, ExcNotCompressed());

				   // check first entry separately, since
				   // for square matrices this is
				   // the diagonal entry (check only
				   // if a first entry exists)
  if (rowstart[i] != rowstart[i+1]) 
    {
      if (j == colnums[rowstart[i]])
	return rowstart[i];
    }
  else
				     // no first entry exists for this
				     // line
    return invalid_entry;

				   // all other entries are sorted, so
				   // we can use a binary seach algorithm
				   //
				   // note that the entries are only sorted
				   // upon compression, so this would fail
				   // for non-compressed sparsity patterns;
				   // however, that is why the Assertion is
				   // at the top of this function, so it
				   // may not be called for noncompressed
				   // structures.
  const unsigned int * const p = lower_bound (&colnums[rowstart[i]+1],
					      &colnums[rowstart[i+1]],
					      j);
  if ((*p == j) &&
      (p != &colnums[rowstart[i+1]]))
    return (p - &colnums[0]);
  else
    return invalid_entry;
}



void
SparsityPattern::add (const unsigned int i,
		      const unsigned int j)
{
  Assert ((rowstart!=0) && (colnums!=0), ExcEmptyObject());  
  Assert (i<rows, ExcInvalidIndex(i,rows));
  Assert (j<cols, ExcInvalidIndex(j,cols));
  Assert (compressed==false, ExcMatrixIsCompressed());

  for (unsigned int k=rowstart[i]; k<rowstart[i+1]; k++)
    {
				       // entry already exists
      if (colnums[k] == j) return;
				       // empty entry found, put new
				       // entry here
      if (colnums[k] == invalid_entry)
	{
	  colnums[k] = j;
	  return;
	};
    };

				   // if we came thus far, something went
				   // wrong: there was not enough space
				   // in this line
  Assert (false, ExcNotEnoughSpace(i, rowstart[i+1]-rowstart[i]));
}



void
SparsityPattern::print_gnuplot (ostream &out) const
{
  Assert ((rowstart!=0) && (colnums!=0), ExcEmptyObject());  
  for (unsigned int i=0; i<rows; ++i)
    for (unsigned int j=rowstart[i]; j<rowstart[i+1]; ++j)
      if (colnums[j] != invalid_entry)
					 // while matrix entries are
					 // usually written (i,j),
					 // with i vertical and j
					 // horizontal, gnuplot output
					 // is x-y, that is we have to
					 // exchange the order of
					 // output
	out << colnums[j] << " " << -static_cast<signed int>(i) << endl;

  AssertThrow (out, ExcIO());
}



unsigned int
SparsityPattern::bandwidth () const
{
  Assert ((rowstart!=0) && (colnums!=0), ExcEmptyObject());  
  unsigned int b=0;
  for (unsigned int i=0; i<rows; ++i)
    for (unsigned int j=rowstart[i]; j<rowstart[i+1]; ++j)
      if (colnums[j] != invalid_entry)
	{
	  if (static_cast<unsigned int>(abs(static_cast<int>(i-colnums[j]))) > b)
	    b = abs(static_cast<signed int>(i-colnums[j]));
	}
      else
					 // leave if at the end of
					 // the entries of this line
	break;
  return b;
};

