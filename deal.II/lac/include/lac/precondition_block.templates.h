//----------------------------  precondition_block.templates.h  ---------------------------
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
//----------------------------  precondition_block.templates.h  ---------------------------
#ifndef __deal2__precondition_block_templates_h
#define __deal2__precondition_block_templates_h


#include <lac/precondition_block.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>


template <typename number, typename inverse_type>
PreconditionBlock<number,inverse_type>::PreconditionBlock ():
		blocksize(0),
		A(0)  {};


template <typename number, typename inverse_type>
PreconditionBlock<number,inverse_type>::~PreconditionBlock ()
{
  if (inverse.size()!=0)
    inverse.erase(inverse.begin(), inverse.end());
}


template <typename number, typename inverse_type>
void PreconditionBlock<number,inverse_type>::clear ()
{
  if (inverse.size()!=0)
    inverse.erase(inverse.begin(), inverse.end());
  blocksize=0;
}


template <typename number, typename inverse_type>
void PreconditionBlock<number,inverse_type>::use_matrix(
  const SparseMatrix<number> &M)
{
  A = &M;
}


template <typename number, typename inverse_type>
void PreconditionBlock<number,inverse_type>::set_block_size(unsigned int bsize) {
  blocksize=bsize;
}


template <typename number, typename inverse_type>
unsigned int PreconditionBlock<number,inverse_type>::block_size() const {
  return blocksize;
}


template <typename number, typename inverse_type>
void PreconditionBlock<number,inverse_type>::invert_diagblocks()
{
  Assert (A!=0, ExcNoMatrixGivenToUse());
  const SparseMatrix<number> &M=*A;
  Assert (M.m() == M.n(), ExcMatrixNotSquare());
  Assert (inverse.size()==0, ExcInverseMatricesAlreadyExist());

  Assert (blocksize!=0, ExcBlockSizeNotSet());
  Assert (M.m()%blocksize==0, ExcWrongBlockSize(blocksize, M.m()));

  const unsigned int n_cells = M.m()/blocksize;

				   // cell_row, cell_column are the
				   // numbering of the blocks (cells).
				   // row_cell, column_cell are the local
				   // numbering of the unknowns in the
				   // blocks.
				   // row, column are the global numbering
				   // of the unkowns.

				   // set the #inverse# array to the right
				   // size. we could do it like this:
				   // inverse = vector<>(n_cells,FullMatrix<>())
				   // but this would involve copying many
				   // FullMatrix objects.
				   //
				   // the following is a neat trick which
				   // avoids copying
  if (true)
    {
      vector<FullMatrix<inverse_type> > tmp(n_cells,
					    FullMatrix<inverse_type>(blocksize));
      inverse.swap (tmp);
    };
  
  FullMatrix<inverse_type> M_cell(blocksize);
  
  for (unsigned int cell=0, row=0; cell<n_cells; ++cell)
    {
      for (unsigned int row_cell=0; row_cell<blocksize; ++row_cell, ++row)
	for (unsigned int column_cell=0, column=cell*blocksize;
	     column_cell<blocksize; ++column_cell, ++column)
	  M_cell(row_cell,column_cell)=M(row,column);
//      try
//	{
      if (blocksize <=4)
	{  
	  inverse[cell].invert(M_cell);
	}
      else
	{
	  M_cell.gauss_jordan();
	  inverse[cell]=M_cell;
	}
//      }
//      catch (ExcNotImplemented &)
    }
}


/*--------------------- PreconditionBlockSOR -----------------------*/

template<typename number, typename inverse_type>
PreconditionBlockSOR<number,inverse_type>::PreconditionBlockSOR(const number omega):
		omega(omega)  {}


template<typename number, typename inverse_type>
PreconditionBlockSOR<number,inverse_type>::~PreconditionBlockSOR(){}


template <typename number, typename inverse_type>
template <typename number2>
void PreconditionBlockSOR<number,inverse_type>::operator() (Vector<number2>       &dst,
							    const Vector<number2> &src) const
{
  Assert(A!=0, ExcNoMatrixGivenToUse());
  const SparseMatrix<number> &M=*A;
  Assert (M.m() == M.n(), ExcMatrixNotSquare());
  Assert (blocksize!=0, ExcBlockSizeNotSet());
  Assert (M.m()%blocksize==0, ExcWrongBlockSize(blocksize, M.m()));
  const unsigned int n_cells=M.m()/blocksize;
  Assert (inverse.size()==0 || inverse.size()==n_cells,
	  ExcWrongNumberOfInverses(inverse.size(), n_cells));

  const SparsityPattern &spars    = M.get_sparsity_pattern();
  const unsigned int    *rowstart = spars.get_rowstart_indices();
  const unsigned int    *columns  = spars.get_column_numbers();

  Vector<number2> b_cell(blocksize), x_cell(blocksize);

				       // cell_row, cell_column are the
				       // numbering of the blocks (cells).
				       // row_cell, column_cell are the local
				       // numbering of the unknowns in the
				       // blocks.
				       // row, column are the global numbering
				       // of the unkowns.
  unsigned int row, column, row_cell, begin_diag_block=0;
  number2 b_cell_row;

  if (inverse.size()==0)
    {
      FullMatrix<number> M_cell(blocksize);
      for (unsigned int cell=0; cell<n_cells; ++cell)
	{
	  for (row=cell*blocksize, row_cell=0; row_cell<blocksize; ++row_cell, ++row)
	    {
	      b_cell_row=src(row);
	      for (unsigned int j=rowstart[row]; j<rowstart[row+1]; ++j)
		if ((column=columns[j]) < begin_diag_block)
		    b_cell_row -= M.global_entry(j) * dst(column);
	      b_cell(row_cell)=b_cell_row;
	      for (unsigned int column_cell=0, column=cell*blocksize;
		   column_cell<blocksize; ++column_cell, ++column)
		  M_cell(row_cell,column_cell)=M(row,column);
	    }
	  M_cell.householder(b_cell);
	  M_cell.backward(x_cell,b_cell);
					   // distribute x_cell to dst
	  for (row=cell*blocksize, row_cell=0; row_cell<blocksize; ++row_cell, ++row)
	    dst(row)=omega*x_cell(row_cell);
	  
	  begin_diag_block+=blocksize;
	}
    }
  else
    for (unsigned int cell=0; cell<n_cells; ++cell)
      {
	for (row=cell*blocksize, row_cell=0; row_cell<blocksize; ++row_cell, ++row)
	  {
	    b_cell_row=src(row);
	    for (unsigned int j=rowstart[row]; j<rowstart[row+1]; ++j)
	      if ((column=columns[j]) < begin_diag_block)
		{
		  b_cell_row -= M.global_entry(j) * dst(column);
		}
	    b_cell(row_cell)=b_cell_row;
	  }
	inverse[cell].vmult(x_cell, b_cell);
					 // distribute x_cell to dst
	for (row=cell*blocksize, row_cell=0; row_cell<blocksize; ++row_cell, ++row)
	  dst(row)=omega*x_cell(row_cell);
	
	begin_diag_block+=blocksize;
      }
}


#endif
