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


#include <base/exceptions.h>
#include <base/logstream.h>
#include <lac/precondition_block.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>


template <typename number, typename inverse_type>
PreconditionBlock<number,inverse_type>::PreconditionBlock ():
		blocksize(0),
		A(0),
		same_diagonal(false)
{};


template <typename number, typename inverse_type>
PreconditionBlock<number,inverse_type>::~PreconditionBlock ()
{
  if (_inverse.size()!=0)
    _inverse.erase(_inverse.begin(), _inverse.end());
}


template <typename number, typename inverse_type>
void PreconditionBlock<number,inverse_type>::clear ()
{
  if (_inverse.size()!=0)
    _inverse.erase(_inverse.begin(), _inverse.end());
  blocksize     = 0;
  same_diagonal = false;
  A = 0;
}


template <typename number, typename inverse_type>
void PreconditionBlock<number,inverse_type>::initialize (const SparseMatrix<number> &M,
							 unsigned int bsize)
{
  clear();
  Assert (M.m() == M.n(), ExcMatrixNotSquare());
  A = &M;
  Assert (bsize>0, ExcIndexRange(bsize, 1, M.m()));
  Assert (A->m()%bsize==0, ExcWrongBlockSize(bsize, A->m()));
  blocksize=bsize;
}


template <typename number, typename inverse_type>
const FullMatrix<inverse_type>&
PreconditionBlock<number,inverse_type>::inverse(unsigned int i) const
{
  if (same_diagonal)
    return _inverse[0];
  
  Assert (i < _inverse.size(), ExcIndexRange(i,0,_inverse.size()));
  return _inverse[i];
}


template <typename number, typename inverse_type>
unsigned int PreconditionBlock<number,inverse_type>::block_size() const
{
  return blocksize;
}


template <typename number, typename inverse_type>
void
PreconditionBlock<number,inverse_type>::set_same_diagonal()
{
  Assert(_inverse.size()==0, ExcInverseMatricesAlreadyExist());
  same_diagonal = true;
}


template <typename number, typename inverse_type>
bool
PreconditionBlock<number,inverse_type>::inverses_ready() const
{
  return (_inverse.size() != 0);
}


template <typename number, typename inverse_type>
void PreconditionBlock<number,inverse_type>::invert_diagblocks()
{
  Assert (A!=0, ExcNotInitialized());
  Assert (blocksize!=0, ExcNotInitialized());

  const SparseMatrix<number> &M=*A;
  Assert (_inverse.size()==0, ExcInverseMatricesAlreadyExist());

  const unsigned int n_cells = M.m()/blocksize;

  FullMatrix<inverse_type> M_cell(blocksize);

  if (same_diagonal)
    {
      deallog << "PreconditionBlock uses only one diagonal block" << endl;
				       // Invert only the first block
				       // This is a copy of the code in the
				       // 'else' part, stripped of the outer loop
      _inverse.resize(1);

      for (unsigned int row_cell=0; row_cell<blocksize; ++row_cell)
	for (unsigned int column_cell=0; column_cell<blocksize; ++column_cell)
	  M_cell(row_cell,column_cell)=M.el(row_cell,column_cell);

      _inverse[0].invert(M_cell);
    }
  else
    {
				       // cell_row, cell_column are the
				       // numbering of the blocks (cells).
				       // row_cell, column_cell are the local
				       // numbering of the unknowns in the
				       // blocks.
				       // row, column are the global numbering
				       // of the unkowns.

				       // set the @p{_inverse} array to the right
				       // size. we could do it like this:
				       // _inverse = vector<>(n_cells,FullMatrix<>())
				       // but this would involve copying many
				       // FullMatrix objects.
				       //
				       // the following is a neat trick which
				       // avoids copying
      vector<FullMatrix<inverse_type> > tmp(n_cells,
					    FullMatrix<inverse_type>(blocksize));
      _inverse.swap (tmp);

      M_cell.clear ();
      
      for (unsigned int cell=0, row=0; cell<n_cells; ++cell)
	{
	  for (unsigned int row_cell=0; row_cell<blocksize; ++row_cell, ++row)
	    for (unsigned int column_cell=0, column=cell*blocksize;
		 column_cell<blocksize; ++column_cell, ++column)
	      M_cell(row_cell,column_cell)=M.el(row,column);

	  _inverse[cell].invert(M_cell);
	}
    }
}


/*--------------------- PreconditionBlockJacobi -----------------------*/


template <typename number, typename inverse_type>
template <typename number2>
void PreconditionBlockJacobi<number,inverse_type>
::vmult (Vector<number2>       &dst,
	 const Vector<number2> &src) const
{
				   // introduce the following typedef
				   // since in the use of exceptions,
				   // strict C++ requires us to
				   // specify them fully as they are
				   // from a template dependent base
				   // class. thus, we'd have to write
				   // PreconditionBlock<number,inverse_type>::ExcNoMatrixGivenToUse,
				   // which is lengthy, but also poses
				   // some problems to the
				   // preprocessor due to the comma in
				   // the template arg list. we could
				   // then wrap the whole thing into
				   // parentheses, but that creates a
				   // parse error for gcc for the
				   // exceptions that do not take
				   // args...
  typedef PreconditionBlock<number,inverse_type> BaseClass;
  Assert(A!=0, ExcNotInitialized());
  
  const SparseMatrix<number> &M=*A;
  const unsigned int n_cells=M.m()/blocksize;

  Vector<number2> b_cell(blocksize), x_cell(blocksize);

				       // cell_row, cell_column are the
				       // numbering of the blocks (cells).
				       // row_cell, column_cell are the local
				       // numbering of the unknowns in the
				       // blocks.
				       // row, column are the global numbering
				       // of the unkowns.
  unsigned int row, row_cell, begin_diag_block=0;

  if (!inverses_ready())
    {
      FullMatrix<number> M_cell(blocksize);
      for (unsigned int cell=0; cell<n_cells; ++cell)
	{
	  for (row=cell*blocksize, row_cell=0; row_cell<blocksize; ++row_cell, ++row)
	    {
	      b_cell(row_cell)=src(row);
	      for (unsigned int column_cell=0, column=cell*blocksize;
		   column_cell<blocksize; ++column_cell, ++column)
		M_cell(row_cell,column_cell)=M(row,column);
	    }
	  M_cell.householder(b_cell);
	  M_cell.backward(x_cell,b_cell);
					   // distribute x_cell to dst
	  for (row=cell*blocksize, row_cell=0; row_cell<blocksize; ++row_cell, ++row)
	    dst(row)=x_cell(row_cell);
	  
	  begin_diag_block+=blocksize;
	}
    }
  else
    for (unsigned int cell=0; cell<n_cells; ++cell)
      {
	for (row=cell*blocksize, row_cell=0; row_cell<blocksize; ++row_cell, ++row)
	  {
	    b_cell(row_cell)=src(row);
	  }
	inverse(cell).vmult(x_cell, b_cell);
					 // distribute x_cell to dst
	for (row=cell*blocksize, row_cell=0; row_cell<blocksize; ++row_cell, ++row)
	  dst(row)=x_cell(row_cell);
	
	begin_diag_block+=blocksize;
      }
}


template <typename number, typename inverse_type>
template <typename number2>
void PreconditionBlockJacobi<number,inverse_type>
::Tvmult (Vector<number2>       &dst,
	 const Vector<number2> &src) const
{
  vmult(dst, src);
}


/*--------------------- PreconditionBlockSOR -----------------------*/

template<typename number, typename inverse_type>
PreconditionBlockSOR<number,inverse_type>::PreconditionBlockSOR(const number omega):
		omega(omega)
{}


template<typename number, typename inverse_type>
PreconditionBlockSOR<number,inverse_type>::~PreconditionBlockSOR()
{}



template<typename number, typename inverse_type>
void
PreconditionBlockSOR<number,inverse_type>::set_omega(number om)
{
  omega = om;
}


//TODO: implement Tvmult

template <typename number, typename inverse_type>
template <typename number2>
void PreconditionBlockSOR<number,inverse_type>::vmult (Vector<number2>       &dst,
						       const Vector<number2> &src) const
{
				   // introduce the following typedef
				   // since in the use of exceptions,
				   // strict C++ requires us to
				   // specify them fully as they are
				   // from a template dependent base
				   // class. thus, we'd have to write
				   // PreconditionBlock<number,inverse_type>::ExcNoMatrixGivenToUse,
				   // which is lengthy, but also poses
				   // some problems to the
				   // preprocessor due to the comma in
				   // the template arg list. we could
				   // then wrap the whole thing into
				   // parentheses, but that creates a
				   // parse error for gcc for the
				   // exceptions that do not take
				   // args...
  typedef PreconditionBlock<number,inverse_type> BaseClass;

  Assert (A!=0, ExcNotInitialized());
  
  const SparseMatrix<number> &M=*A;
  const unsigned int n_cells=M.m()/blocksize;

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

  if (!inverses_ready())
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
	inverse(cell).vmult(x_cell, b_cell);
					 // distribute x_cell to dst
	for (row=cell*blocksize, row_cell=0; row_cell<blocksize; ++row_cell, ++row)
	  dst(row)=omega*x_cell(row_cell);
	
	begin_diag_block+=blocksize;
      }
}


#endif
