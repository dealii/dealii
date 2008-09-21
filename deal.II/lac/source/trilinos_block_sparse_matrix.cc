//---------------------------------------------------------------------------
//    $Id: trilinos_block_sparse_matrix.cc 15631 2008-01-17 23:47:31Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <lac/trilinos_block_sparse_matrix.h>


#ifdef DEAL_II_USE_TRILINOS

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  BlockSparseMatrix::BlockSparseMatrix ()
  {}



  BlockSparseMatrix::~BlockSparseMatrix ()
  {}



  BlockSparseMatrix &
  BlockSparseMatrix::operator = (const BlockSparseMatrix &m)
  {
    BaseClass::operator = (m);

    return *this;
  }



  BlockSparseMatrix &
  BlockSparseMatrix::operator = (const double d)
  {
    for (unsigned int r=0; r<this->n_block_rows(); ++r)
      for (unsigned int c=0; c<this->n_block_cols(); ++c)
        this->block(r,c) = d;

    return *this;
  }



  void
  BlockSparseMatrix::
  reinit (const unsigned int n_block_rows,
          const unsigned int n_block_columns)
  {
                                     // first delete previous content of
                                     // the subobjects array
    for (unsigned int r=0; r<this->n_block_rows(); ++r)
      for (unsigned int c=0; c<this->n_block_cols(); ++c)
        {
          BlockType *p = this->sub_objects[r][c];
          this->sub_objects[r][c] = 0;
          delete p;
        }
  
    this->sub_objects.reinit (0,0);

                                     // then resize. set sizes of blocks to
                                     // zero. user will later have to call
                                     // collect_sizes for this
    this->sub_objects.reinit (n_block_rows,
                              n_block_columns);
    this->row_block_indices.reinit (n_block_rows, 0);
    this->column_block_indices.reinit (n_block_columns, 0);

                                     // and reinitialize the blocks
    for (unsigned int r=0; r<this->n_block_rows(); ++r)
      for (unsigned int c=0; c<this->n_block_cols(); ++c)
        {
          BlockType *p = new BlockType();
          this->sub_objects[r][c] = p;
        }
  }



  void
  BlockSparseMatrix::
  reinit (const std::vector<Epetra_Map> &input_maps,
	  const BlockSparsityPattern    &block_sparsity_pattern)
  {
    Assert (input_maps.size() == block_sparsity_pattern.n_block_rows(),
	    ExcDimensionMismatch (input_maps.size(),
				  block_sparsity_pattern.n_block_rows()));
    Assert (input_maps.size() == block_sparsity_pattern.n_block_cols(),
	    ExcDimensionMismatch (input_maps.size(),
				  block_sparsity_pattern.n_block_cols()));
    
    const unsigned int n_block_rows = input_maps.size();

    if (input_maps[0].Comm().MyPID()==0)
      {
	Assert (n_block_rows == block_sparsity_pattern.n_block_rows(),
		ExcDimensionMismatch (n_block_rows,
				      block_sparsity_pattern.n_block_rows()));
	Assert (n_block_rows == block_sparsity_pattern.n_block_cols(),
		ExcDimensionMismatch (n_block_rows,
				      block_sparsity_pattern.n_block_cols()));
      }

				     // Call the other basic reinit function ...
    reinit (n_block_rows, n_block_rows);
	
				     // ... and then assign the correct
				     // data to the blocks.
    for (unsigned int r=0; r<this->n_block_rows(); ++r)
      for (unsigned int c=0; c<this->n_block_cols(); ++c)
        {
          this->block(r,c).reinit(input_maps[r],input_maps[c],
				  block_sparsity_pattern.block(r,c));
        }
	
    collect_sizes();
  }



  void
  BlockSparseMatrix::
  reinit (const BlockSparsityPattern    &block_sparsity_pattern)
  {
    Assert (block_sparsity_pattern.n_block_rows() ==
	    block_sparsity_pattern.n_block_cols(),
	    ExcDimensionMismatch (block_sparsity_pattern.n_block_rows(),
				  block_sparsity_pattern.n_block_cols()));
    Assert (block_sparsity_pattern.n_rows() ==
	    block_sparsity_pattern.n_cols(),
	    ExcDimensionMismatch (block_sparsity_pattern.n_rows(),
				  block_sparsity_pattern.n_cols()));
    
				     // produce a dummy local map and pass it
				     // off to the other function
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
    Epetra_MpiComm    trilinos_communicator (MPI_COMM_WORLD);
#else
    Epetra_SerialComm trilinos_communicator;
#endif

    std::vector<Epetra_Map> input_maps;
    for (unsigned int i=0; i<block_sparsity_pattern.n_block_rows(); ++i)
      input_maps.push_back (Epetra_Map(block_sparsity_pattern.block(i,0).n_rows(),
				       0,
				       trilinos_communicator));

    reinit (input_maps, block_sparsity_pattern);
  }
  


  void
  BlockSparseMatrix::
  reinit (const std::vector<Epetra_Map>             &input_maps,
	  const ::dealii::BlockSparseMatrix<double> &dealii_block_sparse_matrix,
	  const double                               drop_tolerance)
  {
    const unsigned int n_block_rows = input_maps.size();
    
    Assert (n_block_rows == dealii_block_sparse_matrix.n_block_rows(),
	    ExcDimensionMismatch (n_block_rows,
				  dealii_block_sparse_matrix.n_block_rows()));
    Assert (n_block_rows == dealii_block_sparse_matrix.n_block_cols(),
	    ExcDimensionMismatch (n_block_rows,
				  dealii_block_sparse_matrix.n_block_cols()));

				     // Call the other basic reinit function ...
    reinit (n_block_rows, n_block_rows);
	
				     // ... and then assign the correct
				     // data to the blocks.
    for (unsigned int r=0; r<this->n_block_rows(); ++r)
      for (unsigned int c=0; c<this->n_block_cols(); ++c)
        {
          this->block(r,c).reinit(input_maps[r],input_maps[c],
				  dealii_block_sparse_matrix.block(r,c),
				  drop_tolerance);
        }

    collect_sizes();
  }



  void
  BlockSparseMatrix::compress()
  {
    for (unsigned int r=0; r<this->n_block_rows(); ++r)
      for (unsigned int c=0; c<this->n_block_cols(); ++c)
	this->block(r,c).compress();
  }



  void
  BlockSparseMatrix::collect_sizes ()
  {
    compress();
    BaseClass::collect_sizes ();
  }



  unsigned int
  BlockSparseMatrix::n_nonzero_elements () const
  {
    unsigned int n_nonzero = 0;
    for (unsigned int rows = 0; rows<this->n_block_rows(); ++rows)
      for (unsigned int cols = 0; cols<this->n_block_cols(); ++cols)
	n_nonzero += this->block(rows,cols).n_nonzero_elements();

    return n_nonzero;
  }

}


DEAL_II_NAMESPACE_CLOSE

#endif
