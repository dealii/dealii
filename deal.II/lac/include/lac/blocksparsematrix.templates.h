/*----------------------------   blocksparsematrix.templates.h     --------------------*/
/*      $Id$                 */
/*           Ralf Hartmann, University of Heidelberg                                   */
/*----------------------------   blocksparsematrix.templates.h     ---------------------*/



#include <lac/blocksparsematrix.h>
#include <lac/vector.h>


template <typename number, typename inverse_type>
BlockSparseMatrix<number, inverse_type>::BlockSparseMatrix ():
		blocksize(0) {};


template <typename number, typename inverse_type>
BlockSparseMatrix<number, inverse_type>::~BlockSparseMatrix ()
{
  if (inverse.size()!=0)
    inverse.erase(inverse.begin(), inverse.end());
}


template <typename number, typename inverse_type>
void BlockSparseMatrix<number, inverse_type>::reinit ()
{
  if (inverse.size()!=0)
    inverse.erase(inverse.begin(), inverse.end());
  blocksize=0;
  SparseMatrix<number>::reinit ();
}


template <typename number, typename inverse_type>
void BlockSparseMatrix<number, inverse_type>::reinit (const SparseMatrixStruct &sparsity)
{
  if (inverse.size()!=0)
    inverse.erase(inverse.begin(), inverse.end());
  blocksize=0;
  SparseMatrix<number>::reinit (sparsity);
}


template <typename number, typename inverse_type>
void BlockSparseMatrix<number, inverse_type>::clear ()
{
  SparseMatrix<number>::clear();
  if (inverse.size()!=0)
    inverse.erase(inverse.begin(), inverse.end());
  blocksize=0;
}

template <typename number, typename inverse_type>
void BlockSparseMatrix<number, inverse_type>::set_block_size(unsigned int bsize) {
  blocksize=bsize;
}


template <typename number, typename inverse_type>
unsigned int BlockSparseMatrix<number, inverse_type>::block_size() const {
  return blocksize;
}


template <typename number, typename inverse_type>
template <typename number2>
void BlockSparseMatrix<number, inverse_type>::precondition_BlockSOR (Vector<number2> &dst, const Vector<number2> &src, const number) const
{
  Assert (m() == n(), ExcMatrixNotSquare());
  Assert (blocksize!=0, ExcBlockSizeNotSet());
  Assert (m()%blocksize==0, ExcWrongBlockSize(blocksize, m()));
  unsigned int n_cells=m()/blocksize;
  Assert (inverse.size()==0 || inverse.size()==n_cells,
	  ExcWrongNumberOfInverses(inverse.size(), n_cells));

  const SparseMatrixStruct &spars=get_sparsity_pattern();
  const unsigned int *rowstart = spars.get_rowstart_indices();
  const int *columns = spars.get_column_numbers();

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
		if ((column=static_cast<unsigned int>(columns[j]))
		    < begin_diag_block)
		    b_cell_row -= global_entry(j) * dst(column);
	      b_cell(row_cell)=b_cell_row;
	      for (unsigned int column_cell=0, column=cell*blocksize;
		   column_cell<blocksize; ++column_cell, ++column)
		  M_cell(row_cell,column_cell)=(*this)(row,column);
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
	    b_cell_row=src(row);
	    for (unsigned int j=rowstart[row]; j<rowstart[row+1]; ++j)
	      if ((column=static_cast<unsigned int>(columns[j])) < begin_diag_block)
		{
		  b_cell_row -= global_entry(j) * dst(column);
		}
	    b_cell(row_cell)=b_cell_row;
	  }
	inverse[cell].vmult(x_cell, b_cell);
					 // distribute x_cell to dst
	for (row=cell*blocksize, row_cell=0; row_cell<blocksize; ++row_cell, ++row)
	  dst(row)=x_cell(row_cell);
	
	begin_diag_block+=blocksize;
      }
}



template <typename number, typename inverse_type>
void BlockSparseMatrix<number, inverse_type>::invert_diagblocks()
{
  Assert (m() == n(), ExcMatrixNotSquare());
  Assert (inverse.size()==0, ExcInverseMatricesAlreadyExist());

  Assert (blocksize!=0, ExcBlockSizeNotSet());
  Assert (m()%blocksize==0, ExcWrongBlockSize(blocksize, m()));

  unsigned int n_cells=m()/blocksize;

				   // cell_row, cell_column are the
				   // numbering of the blocks (cells).
				   // row_cell, column_cell are the local
				   // numbering of the unknowns in the
				   // blocks.
				   // row, column are the global numbering
				   // of the unkowns.

  inverse = vector<FullMatrix<inverse_type> > (n_cells, FullMatrix<inverse_type>(blocksize));
  FullMatrix<inverse_type> M_cell(blocksize);
  
  for (unsigned int cell=0, row=0; cell<n_cells; ++cell)
    {
      for (unsigned int row_cell=0; row_cell<blocksize; ++row_cell, ++row)
	for (unsigned int column_cell=0, column=cell*blocksize;
	     column_cell<blocksize; ++column_cell, ++column)
	  M_cell(row_cell,column_cell)=(*this)(row,column);
      if (blocksize<=4)
	inverse[cell].invert(M_cell);
      else
	{
	  M_cell.gauss_jordan();
	  inverse[cell]=M_cell;
	}
    }
}



/*----------------------------   blocksparsematrix.templates.h     ---------------------*/
