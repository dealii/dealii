//----------------------------  sparse_ilu.templates.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparse_ilu.templates.h  ---------------------------
#ifndef __deal2__sparse_ilu_templates_h
#define __deal2__sparse_ilu_templates_h



#include <lac/vector.h>
#include <lac/sparse_ilu.h>

#include <algorithm>
#include <cmath>


template <typename number>
SparseILU<number>::SparseILU () 
{};



template <typename number>
SparseILU<number>::SparseILU (const SparsityPattern &sparsity) :
		SparseMatrix<number> (sparsity)
{};



template <typename number>
void SparseILU<number>::reinit ()
{
  SparseMatrix<number>::reinit ();
};



template <typename number>
void SparseILU<number>::reinit (const SparsityPattern &sparsity)
{
  SparseMatrix<number>::reinit (sparsity);
};



template <typename number>
template <typename somenumber>
void SparseILU<number>::decompose (const SparseMatrix<somenumber> &matrix,
				   const double                    strengthen_diagonal)
{
  Assert (matrix.m()==matrix.n(), ExcMatrixNotSquare ());
  Assert (m()==n(),               ExcMatrixNotSquare ());
  Assert (matrix.m()==m(),        ExcSizeMismatch(matrix.m(), m()));
  
  Assert (strengthen_diagonal>=0, ExcInvalidStrengthening (strengthen_diagonal));


				   // first thing: copy over all elements
				   // of @p{matrix} to the present object
				   //
				   // note that some elements in this
				   // matrix may not be in @p{matrix},
				   // so we need to preset our matrix
				   // by zeroes.
  if (true)
    {
				       // preset the elements
      std::fill_n (&global_entry(0),
		   n_nonzero_elements(),
		   0);

				       // note: pointers to the sparsity
				       // pattern of the old matrix!
      const unsigned int * const rowstart_indices
	= matrix.get_sparsity_pattern().get_rowstart_indices();
      const unsigned int * const column_numbers
	= matrix.get_sparsity_pattern().get_column_numbers();
      
      for (unsigned int row=0; row<m(); ++row)
	for (const unsigned int * col = &column_numbers[rowstart_indices[row]];
	     col != &column_numbers[rowstart_indices[row+1]]; ++col)
	  set (row, *col, matrix.global_entry(col-column_numbers));
    };

  if (strengthen_diagonal > 0)
    for (unsigned int row=0; row<m(); ++row)
      {
					 // get the length of the row
					 // (without the diagonal element)
	const unsigned int rowlength = get_sparsity_pattern().get_rowstart_indices()[row+1]
				       -get_sparsity_pattern().get_rowstart_indices()[row]
				       -1;
	
					 // get the global index of the first
					 // non-diagonal element in this row
	const unsigned int rowstart
	  = get_sparsity_pattern().get_rowstart_indices()[row] + 1;
	number * const diagonal_element = &global_entry(rowstart-1);

	number rowsum = 0;
	for (unsigned int global_index=rowstart;
	     global_index<rowstart+rowlength; ++global_index)
	  rowsum += fabs(global_entry(global_index));

	*diagonal_element += strengthen_diagonal * rowsum;
      };


				   // now work only on this
				   // matrix
  const SparsityPattern             &sparsity = get_sparsity_pattern();
  const unsigned int * const rowstart_indices = sparsity.get_rowstart_indices();
  const unsigned int * const column_numbers   = sparsity.get_column_numbers();
  
/*
  PSEUDO-ALGORITHM
  (indices=0..N-1)
  
  for i=1..N-1
    a[i-1,i-1] = a[i-1,i-1]^{-1}

    for k=0..i-1
      a[i,k] = a[i,k] * a[k,k]

      for j=k+1..N-1
        if (a[i,j] exists & a[k,j] exists)
          a[i,j] -= a[i,k] * a[k,j]
*/


				   // i := row
  for (unsigned int row=1; row<m(); ++row)
    {
				       // invert diagonal element of the
				       // previous row. this is a hack,
				       // which is possible since this
				       // element is not needed any more
				       // in the process of decomposition
				       // and since it makes the backward
				       // step when applying the decomposition
				       // significantly faster
      AssertThrow((global_entry(rowstart_indices[row-1]) !=0),
		  ExcDivideByZero());
      
      global_entry (rowstart_indices[row-1])
	= 1./global_entry (rowstart_indices[row-1]);

				       // let k run over all lower-left
				       // elements of row i; skip
				       // diagonal element at start
      const unsigned int * first_of_row
	= &column_numbers[rowstart_indices[row]+1];
      const unsigned int * first_after_diagonal
	= std::lower_bound (&column_numbers[rowstart_indices[row]+1],
			    &column_numbers[rowstart_indices[row+1]],
			    row);

				       // k := *col_ptr
      for (const unsigned int * col_ptr = first_of_row; col_ptr!=first_after_diagonal; ++col_ptr)
	{
	  const unsigned int global_index_ik = col_ptr-column_numbers;
	  global_entry(global_index_ik) *= diag_element(*col_ptr);

					   // now do the inner loop over
					   // j. note that we need to do
					   // it in the right order, i.e.
					   // taking into account that the
					   // columns are sorted within each
					   // row correctly, but excluding
					   // the main diagonal entry
	  bool left_of_diagonal = true;
	  for (const unsigned int * j = col_ptr+1;
	       j<&column_numbers[rowstart_indices[row+1]];
	       ++j)
	    {
					       // note: this inner loop could
					       // be made considerable faster
					       // if we consulted the row
					       // with number *col_ptr,
					       // instead of always asking
					       // sparsity(*col_ptr,*j),
					       // since we traverse this
					       // row linearly. I just didn't
					       // have the time to figure out
					       // the details.
	      
					       // check whether we have just
					       // traversed the diagonal
					       //
					       // note that j should never point
					       // to the diagonal itself!
	      if (left_of_diagonal && (*j > row))
		{
		  Assert (*j != row, ExcInternalError());
		  
		  left_of_diagonal = false;
						   // a[i,i] -= a[i,k]*a[k,i]
		  const int global_index_ki = sparsity(*col_ptr,row);

		  if (global_index_ki != -1)
		    diag_element(row) -= global_entry(global_index_ik) *
					 global_entry(global_index_ki);
		  
		};
	      
       	      const int global_index_ij = j - &column_numbers[0],
			global_index_kj = sparsity(*col_ptr,*j);
	      if ((global_index_ij != -1) &&
		  (global_index_kj != -1))
		global_entry(global_index_ij) -= global_entry(global_index_ik) *
						 global_entry(global_index_kj);
	    };
	};
    };

/*
  OLD CODE, rather crude first implementation with an algorithm taken
  from 'W. Hackbusch, G. Wittum: Incomplete Decompositions (ILU)-
  Algorithms, Theory, and Applications', page 6.
  
  for (unsigned int k=0; k<m()-1; ++k)
    for (unsigned int i=k+1; i<m(); ++i)
      {
					 // get the global index
					 // of the element (i,k)
	const int global_index_ik = get_sparsity_pattern()(i,k);

					 // if this element is zero,
					 // then we continue with the
					 // next i, since e would be
					 // zero and nothing would happen
					 // in this loop
	if (global_index_ik == -1)
	  continue;
	
	const number e = global_entry(global_index_ik) / diag_element(k);
	global_entry(global_index_ik) = e;

	for (unsigned int j=k+1; j<m(); ++j)
	  {
					     // find out about a_kj
					     // if this does not exist,
					     // then the updates within
					     // this innermost loop would
					     // be zero, invariable of the
					     // fact of whether a_ij is a
					     // nonzero or a zero element
	    const int global_index_kj = get_sparsity_pattern()(k,j);
	    if (global_index_kj == -1)
	      continue;

	    const int global_index_ij = get_sparsity_pattern()(i,j);
	    if (global_index_ij != -1)
	      global_entry(global_index_ij) -= e*global_entry(global_index_kj);
	    else
	      diag_element(i) -= e*global_entry(global_index_kj);
	  };
      };
*/      
};



template <typename number>
template <typename somenumber>
void SparseILU<number>::apply_decomposition (Vector<somenumber>       &dst,
					     const Vector<somenumber> &src) const 
{
  Assert (dst.size() == src.size(), ExcSizeMismatch(dst.size(), src.size()));
  Assert (dst.size() == m(), ExcSizeMismatch(dst.size(), m()));
  
  const unsigned int N=dst.size();
  const unsigned int * const rowstart_indices = get_sparsity_pattern().get_rowstart_indices();
  const unsigned int * const column_numbers   = get_sparsity_pattern().get_column_numbers();
				   // solve LUx=b in two steps:
				   // first Ly = b, then
				   //       Ux = y
				   //
				   // first a forward solve. since
				   // the diagonal values of L are
				   // one, there holds
				   // y_i = b_i
				   //       - sum_{j=0}^{i-1} L_{ij}y_j
				   // we split the y_i = b_i off and
				   // perform it at the outset of the
				   // loop
  dst = src;
  for (unsigned int row=0; row<N; ++row)
    {
				       // get start of this row. skip the
				       // diagonal element
      const unsigned int * const rowstart = &column_numbers[rowstart_indices[row]+1];
				       // find the position where the part
				       // right of the diagonal starts
      const unsigned int * const first_after_diagonal
	= std::lower_bound (rowstart,
			    &column_numbers[rowstart_indices[row+1]],
			    row);
      
      for (const unsigned int * col=rowstart; col!=first_after_diagonal; ++col)
	dst(row) -= global_entry (col-column_numbers) * dst(*col);
    };

				   // now the backward solve. same
				   // procedure, but we need not set
				   // dst before, since this is already
				   // done.
				   //
				   // note that we need to scale now,
				   // since the diagonal is not zero
				   // now
  for (int row=N-1; row>=0; --row)
    {
				       // get end of this row
      const unsigned int * const rowend = &column_numbers[rowstart_indices[row+1]];
				       // find the position where the part
				       // right of the diagonal starts
      const unsigned int * const first_after_diagonal
	= std::lower_bound (&column_numbers[rowstart_indices[row]+1],
			    &column_numbers[rowstart_indices[row+1]],
			    static_cast<unsigned int>(row));
      
      for (const unsigned int * col=first_after_diagonal; col!=rowend; ++col)
	dst(row) -= global_entry (col-column_numbers) * dst(*col);

				       // scale by the diagonal element.
				       // note that the diagonal element
				       // was stored inverted
      dst(row) *= diag_element(row);
    };
};



template <typename number>
unsigned int
SparseILU<number>::memory_consumption () const
{
  return SparseMatrix<number>::memory_consumption ();
};



/*----------------------------   sparse_ilu.templates.h     ---------------------------*/

#endif
/*----------------------------   sparse_ilu.templates.h     ---------------------------*/
