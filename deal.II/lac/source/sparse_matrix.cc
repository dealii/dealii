// $Id$

// This file was once part of the DEAL Library
// DEAL is Copyright(1995) by
// Roland Becker, Guido Kanschat, Franz-Theo Suttmeier
// Revised, modified and extended by Wolfgang Bangerth, 1998, 1999


#include <lac/sparsematrix.h>
#include <lac/sparsematrix.templates.h>
#include <lac/ivector.h>

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <numeric>


SparseMatrixStruct::SparseMatrixStruct () :
		max_dim(0),
		max_vec_len(0),
		rowstart(0),
		colnums(0)
{
  reinit (0,0,0);
};



SparseMatrixStruct::SparseMatrixStruct (const SparseMatrixStruct &s) :
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



SparseMatrixStruct::SparseMatrixStruct (const unsigned int m,
					const unsigned int n,
					const unsigned int max_per_row) 
		: max_dim(0),
		  max_vec_len(0),
		  rowstart(0),
		  colnums(0)
{
  reinit (m,n,max_per_row);
};



SparseMatrixStruct::SparseMatrixStruct (const unsigned int n,
					const unsigned int max_per_row)
		: max_dim(0),
		  max_vec_len(0),
		  rowstart(0),
		  colnums(0)
{
  reinit (n,n,max_per_row);
};


SparseMatrixStruct::SparseMatrixStruct (const SparseMatrixStruct &original,
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
      const int * const
	original_row_start = &original.colnums[original.rowstart[row]] + 1;
				       // thw following requires that
				       // #original# be compressed since
				       // otherwise there might be -1's
      const int * const
	original_row_end   = &original.colnums[original.rowstart[row+1]];

      const int * const
	original_last_before_side_diagonals = lower_bound (original_row_start,
							   original_row_end,
							   static_cast<int>
							   (row
							    -extra_off_diagonals));
      const int * const
	original_first_after_side_diagonals = upper_bound (original_row_start,
							   original_row_end,
							   static_cast<int>
							   (row
							    +extra_off_diagonals));

      int * next_free_slot = &colnums[rowstart[row]] + 1;

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



SparseMatrixStruct::~SparseMatrixStruct ()
{
  if (rowstart != 0)  delete[] rowstart;
  if (colnums != 0)   delete[] colnums;
}




void
SparseMatrixStruct::reinit (const unsigned int m,
			    const unsigned int n,
			    const unsigned int max_per_row)
{
  Assert ((max_per_row>0) || ((m==0) && (n==0)), ExcInvalidNumber(max_per_row));
  rows = m;
  cols = n;
  vec_len = m * max_per_row;
  max_row_len = max_per_row;

				   // delete empty matrices
  if ((m==0) || (n==0))
    {
      if (rowstart)  delete[] rowstart;
      if (colnums)   delete[] colnums;
      rowstart = 0;
      colnums = 0;
      max_vec_len = vec_len = max_dim = rows = cols = 0;
      compressed = false;
      return;
    };
  
  if (rows > max_dim)
    {
      if (rowstart) delete[] rowstart;
      max_dim = rows;
      rowstart = new unsigned int[max_dim+1];
    };
  
  if (vec_len > max_vec_len)
    {
      if (colnums) delete[] colnums;
      max_vec_len = vec_len;
      colnums = new int[max_vec_len];
    };
  
  for (unsigned int i=0; i<=rows; i++)
    rowstart[i] = i * max_per_row;
  fill_n (&colnums[0], vec_len, -1);

  if (rows == cols)
    for (unsigned int i=0;i<rows;i++)
      colnums[rowstart[i]] = i;

  compressed = false;
}


void
SparseMatrixStruct::compress ()
{
  Assert ((rowstart!=0) && (colnums!=0), ExcEmptyObject());
  
  if (compressed) return;
  unsigned int next_free_entry = 0,
		next_row_start = 0,
		    row_length = 0;

				   // reserve temporary storage to
				   // store the entries of one wor
  int *tmp_entries = new int[max_row_len];

				   // Traverse all rows
  for (unsigned int line=0; line<rows; ++line)
    {
				       // copy used entries, break if
				       // first unused entry is reached
      row_length = 0;
      for (unsigned int j=rowstart[line]; j<rowstart[line+1]; ++j,++row_length)
	if (colnums[j] != -1)
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

				       // Re-insert column numbers
				       // into the field
      for (unsigned int j=0; j<row_length; ++j)
	colnums[next_free_entry++] = tmp_entries[j];

				       // note new start of this and
				       // the next row
      rowstart[line] = next_row_start;
      next_row_start = next_free_entry;

				       // some internal checks
      Assert ((rows!=cols) ||
	      (colnums[rowstart[line]] == static_cast<signed int>(line)),
	      ExcInternalError());
				       // assert that the first entry
				       // does not show up in
				       // the remaining ones and that
				       // the remaining ones are unique
				       // among themselves (this handles
				       // both cases, quadratic and
				       // rectangular matrices)
      Assert (find (&colnums[rowstart[line]+1],
		    &colnums[next_row_start],
		    colnums[rowstart[line]]) ==
	      &colnums[next_row_start],
	      ExcInternalError());
      Assert (adjacent_find(&colnums[rowstart[line]+1],
			    &colnums[next_row_start]) ==
	      &colnums[next_row_start],
	      ExcInternalError());
    };
  
  vec_len = rowstart[rows] = next_row_start;
  compressed = true;

  delete[] tmp_entries;
};



bool
SparseMatrixStruct::empty () const {
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
      Assert (vec_len==0, ExcInternalError());
      Assert (max_vec_len==0, ExcInternalError());
      Assert (vec_len==0, ExcInternalError());

      return true;
    };
  return false;
};



int
SparseMatrixStruct::operator () (const unsigned int i, const unsigned int j) const
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
      if (static_cast<signed int>(j) == colnums[rowstart[i]])
	return rowstart[i];
    }
  else
				     // no first entry exists for this
				     // line
    return -1;

				   // all other entries are sorted, so
				   // we can use a binary seach algorithm
  const int * const p = lower_bound (&colnums[rowstart[i]+1],
				     &colnums[rowstart[i+1]],
				     static_cast<signed int>(j));
  if ((*p == static_cast<signed int>(j)) &&
      (p != &colnums[rowstart[i+1]]))
    return (p - &colnums[0]);
  else
    return -1;
}


void
SparseMatrixStruct::add (const unsigned int i, const unsigned int j)
{
  Assert ((rowstart!=0) && (colnums!=0), ExcEmptyObject());  
  Assert (i<rows, ExcInvalidIndex(i,rows));
  Assert (j<cols, ExcInvalidIndex(j,cols));
  Assert (compressed==false, ExcMatrixIsCompressed());

  for (unsigned int k=rowstart[i]; k<rowstart[i+1]; k++)
    {
				       // entry already exists
      if (colnums[k] == (signed int)j) return;
				       // empty entry found, put new
				       // entry here
      if (colnums[k] == -1)
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
SparseMatrixStruct::add_matrix (const unsigned int n, const int* rowcols)
{
  Assert ((rowstart!=0) && (colnums!=0), ExcEmptyObject());  
  for (unsigned int i=0; i<n; ++i)
    for (unsigned int j=0; j<n; ++j)
      add(rowcols[i], rowcols[j]);
}



void
SparseMatrixStruct::add_matrix (const unsigned int m, const unsigned int n,
				const int* rows, const int* cols)
{
  Assert ((rowstart!=0) && (colnums!=0), ExcEmptyObject());  
  for (unsigned i=0; i<m; ++i)
    for (unsigned j=0; j<n; ++j)
      add(rows[i], cols[j]);
}



void
SparseMatrixStruct::add_matrix (const iVector& rowcols)
{
  Assert ((rowstart!=0) && (colnums!=0), ExcEmptyObject());  
  unsigned int i,j;
  for (i=0;i<rowcols.n();i++)
    for (j=0;j<rowcols.n();j++)
      add(rowcols(i), rowcols(j));
}



void
SparseMatrixStruct::add_matrix (const iVector& rows, const iVector& cols)
{
  Assert ((rowstart!=0) && (colnums!=0), ExcEmptyObject());  
  unsigned int i,j;
  for (i=0;i<rows.n();i++)
    for (j=0;j<cols.n();j++)
      add(rows(i), cols(j));
}



void
SparseMatrixStruct::print_gnuplot (ostream &out) const
{
  Assert ((rowstart!=0) && (colnums!=0), ExcEmptyObject());  
  for (unsigned int i=0; i<rows; ++i)
    for (unsigned int j=rowstart[i]; j<rowstart[i+1]; ++j)
      if (colnums[j]>=0)
	out << i << " " << -colnums[j] << endl;

  AssertThrow (out, ExcIO());
}



unsigned int
SparseMatrixStruct::bandwidth () const
{
  Assert ((rowstart!=0) && (colnums!=0), ExcEmptyObject());  
  unsigned int b=0;
  for (unsigned int i=0; i<rows; ++i)
    for (unsigned int j=rowstart[i]; j<rowstart[i+1]; ++j)
      if (colnums[j]>=0) 
	{
	  if (static_cast<unsigned int>(abs(static_cast<int>(i-colnums[j]))) > b)
	    b = abs(static_cast<int>(i-colnums[j]));
	}
      else
					 // leave if at the end of
					 // the entries of this line
	break;
  return b;
};



unsigned int
SparseMatrixStruct::n_nonzero_elements () const {
  Assert ((rowstart!=0) && (colnums!=0), ExcEmptyObject());  
  Assert (compressed, ExcNotCompressed());
  return colnums[rows]-colnums[0];
};








// explicit instantiations for "float"
template class SparseMatrix<float>;

template SparseMatrix<float> & SparseMatrix<float>::copy_from (const SparseMatrix<float> &);
template SparseMatrix<float> & SparseMatrix<float>::copy_from (const SparseMatrix<double> &);

template void SparseMatrix<float>::add_scaled (const float, const SparseMatrix<float> &);
template void SparseMatrix<float>::add_scaled (const float, const SparseMatrix<double> &);

template void SparseMatrix<float>::vmult (Vector<float> &, const Vector<float> &) const;
template void SparseMatrix<float>::vmult (Vector<double> &, const Vector<double> &) const;

template void SparseMatrix<float>::Tvmult (Vector<float> &, const Vector<float> &) const;
template void SparseMatrix<float>::Tvmult (Vector<double> &, const Vector<double> &) const;

template double SparseMatrix<float>::matrix_norm (const Vector<float> &) const;
template double SparseMatrix<float>::matrix_norm (const Vector<double> &) const;

template double SparseMatrix<float>::residual (Vector<float> &,
					       const Vector<float> &,
					       const Vector<float> &) const;
template double SparseMatrix<float>::residual (Vector<double> &,
					       const Vector<double> &,
					       const Vector<double> &) const;

template void SparseMatrix<float>::precondition (Vector<float> &,
						 const Vector<float> &) const;
template void SparseMatrix<float>::precondition (Vector<double> &,
						 const Vector<double> &) const;

template void SparseMatrix<float>::precondition_SSOR (Vector<float> &,
						      const Vector<float> &,
						      float) const;
template void SparseMatrix<float>::precondition_SSOR (Vector<double> &,
						      const Vector<double> &,
						      float) const;

template void SparseMatrix<float>::precondition_SOR (Vector<float> &,
						     const Vector<float> &,
						     float) const;
template void SparseMatrix<float>::precondition_SOR (Vector<double> &,
						     const Vector<double> &,
						     float) const;

template void SparseMatrix<float>::precondition_Jacobi (Vector<float> &,
							const Vector<float> &,
							float) const;
template void SparseMatrix<float>::precondition_Jacobi (Vector<double> &,
							const Vector<double> &,
							float) const;

template void SparseMatrix<float>::SOR (Vector<float> &,
					float) const;
template void SparseMatrix<float>::SOR (Vector<double> &,
					float) const;

template void SparseMatrix<float>::SSOR (Vector<float> &,
					 float) const;
template void SparseMatrix<float>::SSOR (Vector<double> &,
					 float) const;



// explicit instantiations for "double"
template class SparseMatrix<double>;

template SparseMatrix<double> & SparseMatrix<double>::copy_from (const SparseMatrix<float> &);
template SparseMatrix<double> & SparseMatrix<double>::copy_from (const SparseMatrix<double> &);

template void SparseMatrix<double>::add_scaled (const double, const SparseMatrix<float> &);
template void SparseMatrix<double>::add_scaled (const double, const SparseMatrix<double> &);

template void SparseMatrix<double>::vmult (Vector<float> &, const Vector<float> &) const;
template void SparseMatrix<double>::vmult (Vector<double> &, const Vector<double> &) const;

template void SparseMatrix<double>::Tvmult (Vector<float> &, const Vector<float> &) const;
template void SparseMatrix<double>::Tvmult (Vector<double> &, const Vector<double> &) const;

template double SparseMatrix<double>::matrix_norm (const Vector<float> &) const;
template double SparseMatrix<double>::matrix_norm (const Vector<double> &) const;

template double SparseMatrix<double>::residual (Vector<float> &,
						const Vector<float> &,
						const Vector<float> &) const;
template double SparseMatrix<double>::residual (Vector<double> &,
						const Vector<double> &,
						const Vector<double> &) const;

template void SparseMatrix<double>::precondition (Vector<float> &,
						  const Vector<float> &) const;
template void SparseMatrix<double>::precondition (Vector<double> &,
						  const Vector<double> &) const;

template void SparseMatrix<double>::precondition_SSOR (Vector<float> &,
						       const Vector<float> &,
						       double) const;
template void SparseMatrix<double>::precondition_SSOR (Vector<double> &,
						       const Vector<double> &,
						       double) const;

template void SparseMatrix<double>::precondition_SOR (Vector<float> &,
						      const Vector<float> &,
						      double) const;
template void SparseMatrix<double>::precondition_SOR (Vector<double> &,
						      const Vector<double> &,
						      double) const;

template void SparseMatrix<double>::precondition_Jacobi (Vector<float> &,
							 const Vector<float> &,
							 double) const;
template void SparseMatrix<double>::precondition_Jacobi (Vector<double> &,
							 const Vector<double> &,
							 double) const;

template void SparseMatrix<double>::SOR (Vector<float> &,
					 double) const;
template void SparseMatrix<double>::SOR (Vector<double> &,
					 double) const;

template void SparseMatrix<double>::SSOR (Vector<float> &,
					  double) const;
template void SparseMatrix<double>::SSOR (Vector<double> &,
					  double) const;
