/*----------------------------   dblocksmatrix.cc     ---------------------------*/
/*      $Id$                 */
/*----------------------------   dblocksmatrix.cc     ---------------------------*/

#include <lac/dblocksmatrix.h>
#include <lac/dvector.h>


dBlockSMatrix::dBlockSMatrix ():
		blocksize(0) {};

dBlockSMatrix::~dBlockSMatrix ()
{
  if (inverse.size()!=0)
    inverse.erase(inverse.begin(), inverse.end());
}


void dBlockSMatrix::reinit ()
{
  if (inverse.size()!=0)
    inverse.erase(inverse.begin(), inverse.end());
  blocksize=0;
  dSMatrix::reinit ();
}


void dBlockSMatrix::reinit (const dSMatrixStruct &sparsity)
{
  if (inverse.size()!=0)
    inverse.erase(inverse.begin(), inverse.end());
  blocksize=0;
  dSMatrix::reinit (sparsity);
}


void dBlockSMatrix::clear ()
{
  dSMatrix::clear();
  if (inverse.size()!=0)
    inverse.erase(inverse.begin(), inverse.end());
  blocksize=0;
}


void dBlockSMatrix::set_block_size(unsigned int bsize) {
  blocksize=bsize;
}



unsigned int dBlockSMatrix::block_size() const {
  return blocksize;
}



void dBlockSMatrix::precondition_BlockSOR (dVector &dst, const dVector &src) const
{
  Assert (m() == n(), ExcMatrixNotSquare());
  Assert (blocksize!=0, ExcBlockSizeNotSet());
  Assert (m()%blocksize==0, ExcWrongBlockSize(blocksize, m()));
  unsigned int n_cells=m()/blocksize;
  Assert (inverse.size()==0 || inverse.size()==n_cells,
	  ExcWrongInverses(inverse.size(), n_cells));

  const dSMatrixStruct &spars=get_sparsity_pattern();
  const unsigned int *rowstart = spars.get_rowstart_indices();
  const int *columns = spars.get_column_numbers();

  dVector b_cell(blocksize), x_cell(blocksize);

				       // cell_row, cell_column are the
				       // numbering of the blocks (cells).
				       // row_cell, column_cell are the local
				       // numbering of the unknowns in the
				       // blocks.
				       // row, column are the global numbering
				       // of the unkowns.
  unsigned int row, column, row_cell, begin_diag_block=0;
  double b_cell_row;

  if (inverse.size()==0)
    {
      dFMatrix M_cell(blocksize);
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


void dBlockSMatrix::invert_diagblocks()
{
  Assert (m() == n(), ExcMatrixNotSquare());
  Assert (inverse.size()==0, ExcInverseMatricesAlreadyExist());

  Assert (blocksize!=0, ExcBlockSizeNotSet());
  Assert (m()%blocksize==0, ExcWrongBlockSize(blocksize, m()));
  unsigned int n_cells=m()/blocksize;

  inverse.insert(inverse.begin(), n_cells, dFMatrix(blocksize));
  
				   // cell_row, cell_column are the
				   // numbering of the blocks (cells).
				   // row_cell, column_cell are the local
				   // numbering of the unknowns in the
				   // blocks.
				   // row, column are the global numbering
				   // of the unkowns.
  dFMatrix M_cell(blocksize);

  for (unsigned int cell=0, row=0; cell<n_cells; ++cell)
    {
      for (unsigned int row_cell=0; row_cell<blocksize; ++row_cell, ++row)
	for (unsigned int column_cell=0, column=cell*blocksize;
	     column_cell<blocksize; ++column_cell, ++column)
	    M_cell(row_cell,column_cell)=(*this)(row,column);

				       // perhaps #dFMatrix::invert# should
				       // be change such that it calls
				       // #gauss_jordan()# automatically
				       // if blocksize > 4
      if (blocksize<=4)
	inverse[cell].invert(M_cell);
      else
	{
	  M_cell.gauss_jordan();
	  inverse[cell]=M_cell;
	}
    }
}



/*----------------------------   dblocksmatrix.cc     ---------------------------*/
