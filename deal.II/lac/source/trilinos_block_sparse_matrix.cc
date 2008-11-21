//---------------------------------------------------------------------------
//    $Id$
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

#include <lac/block_sparsity_pattern.h>

#ifdef DEAL_II_USE_TRILINOS

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  BlockSparseMatrix::BlockSparseMatrix ()
  {}



  BlockSparseMatrix::~BlockSparseMatrix ()
  {
  				   // delete previous content of
				   // the subobjects array
    clear ();
  }



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
    clear ();

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

	  Assert (this->sub_objects[r][c] == 0,
		  ExcInternalError());
          this->sub_objects[r][c] = p;
        }
  }



  template <typename BlockSparsityType>  
  void
  BlockSparseMatrix::
  reinit (const std::vector<Epetra_Map> &input_maps,
	  const BlockSparsityType       &block_sparsity_pattern)
  {
    Assert (input_maps.size() == block_sparsity_pattern.n_block_rows(),
	    ExcDimensionMismatch (input_maps.size(),
				  block_sparsity_pattern.n_block_rows()));
    Assert (input_maps.size() == block_sparsity_pattern.n_block_cols(),
	    ExcDimensionMismatch (input_maps.size(),
				  block_sparsity_pattern.n_block_cols()));
    
    const unsigned int n_block_rows = input_maps.size();

    Assert (n_block_rows == block_sparsity_pattern.n_block_rows(),
	    ExcDimensionMismatch (n_block_rows,
				  block_sparsity_pattern.n_block_rows()));
    Assert (n_block_rows == block_sparsity_pattern.n_block_cols(),
	    ExcDimensionMismatch (n_block_rows,
				  block_sparsity_pattern.n_block_cols()));

    
				     // Call the other basic reinit function, ...
    reinit (block_sparsity_pattern.n_block_rows(),
	    block_sparsity_pattern.n_block_cols());

				     // ... set the correct sizes, ...
    this->row_block_indices    = block_sparsity_pattern.get_row_indices();
    this->column_block_indices = block_sparsity_pattern.get_column_indices();
	
				     // ... and then assign the correct
				     // data to the blocks.
    for (unsigned int r=0; r<this->n_block_rows(); ++r)
      for (unsigned int c=0; c<this->n_block_cols(); ++c)
        {
	  this->sub_objects[r][c]->reinit (input_maps[r], input_maps[c],
					   block_sparsity_pattern.block(r,c));
        }
  }



  template <typename BlockSparsityType>
  void
  BlockSparseMatrix::
  reinit (const BlockSparsityType &block_sparsity_pattern)
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
          this->sub_objects[r][c]->reinit(input_maps[r],input_maps[c],
					  dealii_block_sparse_matrix.block(r,c),
					  drop_tolerance);
        }

    collect_sizes();
  }



  void
  BlockSparseMatrix::
  reinit (const ::dealii::BlockSparseMatrix<double> &dealii_block_sparse_matrix,
	  const double                               drop_tolerance)
  {
    Assert (dealii_block_sparse_matrix.n_block_rows() ==
	    dealii_block_sparse_matrix.n_block_cols(),
	    ExcDimensionMismatch (dealii_block_sparse_matrix.n_block_rows(),
				  dealii_block_sparse_matrix.n_block_cols()));
    Assert (dealii_block_sparse_matrix.m() ==
	    dealii_block_sparse_matrix.n(),
	    ExcDimensionMismatch (dealii_block_sparse_matrix.m(),
				  dealii_block_sparse_matrix.n()));
    
				     // produce a dummy local map and pass it
				     // off to the other function
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
    Epetra_MpiComm    trilinos_communicator (MPI_COMM_WORLD);
#else
    Epetra_SerialComm trilinos_communicator;
#endif

    std::vector<Epetra_Map> input_maps;
    for (unsigned int i=0; i<dealii_block_sparse_matrix.n_block_rows(); ++i)
      input_maps.push_back (Epetra_Map(dealii_block_sparse_matrix.block(i,0).m(),
				       0,
				       trilinos_communicator));

    reinit (input_maps, dealii_block_sparse_matrix, drop_tolerance);
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



  void
  BlockSparseMatrix::set (const unsigned int   i,
			  const unsigned int   j,
			  const TrilinosScalar value)
  {
    BaseClass::set (i, j, value);
  }



  void
  BlockSparseMatrix::set (const std::vector<unsigned int>  &row_indices,
			  const std::vector<unsigned int>  &col_indices,
			  const FullMatrix<TrilinosScalar> &values)
  {
    Assert (row_indices.size() == values.m(),
	    ExcDimensionMismatch(row_indices.size(), values.m()));
    Assert (col_indices.size() == values.n(),
	    ExcDimensionMismatch(col_indices.size(), values.n()));

    set (row_indices.size(), &row_indices[0], 
	 col_indices.size(), &col_indices[0], &values(0,0));
  }



  void
  BlockSparseMatrix::set (const unsigned int                 row,
			  const std::vector<unsigned int>   &col_indices,
			  const std::vector<TrilinosScalar> &values)
  {
    Assert (col_indices.size() == values.size(),
	    ExcDimensionMismatch(col_indices.size(), values.size()));

    set (1, &row, col_indices.size(), &col_indices[0], &values[0]);
  }



  void
  BlockSparseMatrix::set (const unsigned int    n_rows,
			  const unsigned int   *row_indices,
			  const unsigned int    n_cols,
			  const unsigned int   *col_indices,
			  const TrilinosScalar *values)
  {
				   // Resize scratch arrays
    block_col_indices.resize (this->n_block_cols());
    local_row_length.resize (this->n_block_cols());	
    local_col_indices.resize (n_cols);

				   // Clear the content in local_row_length
    for (unsigned int i=0; i<this->n_block_cols(); ++i)
      local_row_length[i] = 0;

				   // Go through the column indices to find
				   // out which portions of the values
				   // should be written into which block
				   // matrix. This can be done before
				   // starting the loop over all the rows,
				   // since we assume a rectangular set of
				   // matrix data.
    {
      unsigned int current_block = 0, row_length = 0;
      block_col_indices[0] = 0;
      for (unsigned int j=0; j<n_cols; ++j, ++row_length)
	{
	  const std::pair<unsigned int, unsigned int>
	    col_index = this->column_block_indices.global_to_local(col_indices[j]);
	  local_col_indices[j] = col_index.second;
	  if (col_index.first > current_block)
	    {
	      local_row_length[current_block] = row_length;
	      row_length = 0;
	      while (col_index.first > current_block)
		current_block++;
	      block_col_indices[current_block] = j;
	    }

	  Assert (col_index.first == current_block,
		  ExcInternalError());
	}
      local_row_length[current_block] = row_length;
      Assert (current_block < this->n_block_cols(),
	      ExcInternalError());

#ifdef DEBUG
				   // If in debug mode, do a check whether
				   // the right length has been obtained.
      unsigned int length = 0;
      for (unsigned int i=0; i<this->n_block_cols(); ++i)
	length += local_row_length[i];
      Assert (length == n_cols,
	      ExcDimensionMismatch(length, n_cols));
#endif
    }

				   // Now we found out about where the
				   // individual columns should start and
				   // where we should start reading out
				   // data. Now let's write the data into
				   // the individual blocks!
    for (unsigned int i=0; i<n_rows; ++i)
      {
	const std::pair<unsigned int,unsigned int> 
	  row_index = this->row_block_indices.global_to_local (row_indices[i]);
	for (unsigned int block_col=0; block_col<n_block_cols(); ++block_col)
	  {
	    if (local_row_length[block_col] == 0)
	      continue;

	    block(row_index.first, block_col).set 
	      (1, &row_index.second, 
	       local_row_length[block_col], 
	       &local_col_indices[block_col_indices[block_col]],
	       &values[n_cols*i + block_col_indices[block_col]]);
	  }
      }
  }



  void
  BlockSparseMatrix::add (const unsigned int   i,
			  const unsigned int   j,
			  const TrilinosScalar value)
  {
				   // For adding one single element, it is
				   // faster to rely on the BaseClass add
				   // function, than doing all the strange
				   // operations below.
    BaseClass::add (i, j, value);
  }



  void
  BlockSparseMatrix::add (const std::vector<unsigned int>  &row_indices,
			  const std::vector<unsigned int>  &col_indices,
			  const FullMatrix<TrilinosScalar> &values)
  {
    Assert (row_indices.size() == values.m(),
	    ExcDimensionMismatch(row_indices.size(), values.m()));
    Assert (col_indices.size() == values.n(),
	    ExcDimensionMismatch(col_indices.size(), values.n()));

    add (row_indices.size(), &row_indices[0], 
	 col_indices.size(), &col_indices[0], &values(0,0));
  }



  void
  BlockSparseMatrix::add (const unsigned int                 row,
			  const std::vector<unsigned int>   &col_indices,
			  const std::vector<TrilinosScalar> &values)
  {
    Assert (col_indices.size() == values.size(),
	    ExcDimensionMismatch(col_indices.size(), values.size()));

    add (1, &row, col_indices.size(), &col_indices[0], &values[0]);
  }



  void
  BlockSparseMatrix::add (const unsigned int    n_rows,
			  const unsigned int   *row_indices,
			  const unsigned int    n_cols,
			  const unsigned int   *col_indices,
			  const TrilinosScalar *values)
  {
				   // TODO: Look over this to find out
				   // whether we can do that more
				   // efficiently.

				   // Resize scratch arrays
    block_col_indices.resize (this->n_block_cols());
    local_row_length.resize (this->n_block_cols());	
    local_col_indices.resize (n_cols);

				   // Clear the content in local_row_length
    for (unsigned int i=0; i<this->n_block_cols(); ++i)
      local_row_length[i] = 0;

				   // Go through the column indices to find
				   // out which portions of the values
				   // should be written into which block
				   // matrix. This can be done before
				   // starting the loop over all the rows,
				   // since we assume a rectangular set of
				   // matrix data.
    {
      unsigned int current_block = 0, row_length = 0;
      block_col_indices[0] = 0;
      for (unsigned int j=0; j<n_cols; ++j, ++row_length)
	{
	  const std::pair<unsigned int, unsigned int>
	    col_index = this->column_block_indices.global_to_local(col_indices[j]);
	  local_col_indices[j] = col_index.second;
	  if (col_index.first > current_block)
	    {
	      local_row_length[current_block] = row_length;
	      row_length = 0;
	      while (col_index.first > current_block)
		current_block++;
	      block_col_indices[current_block] = j;
	    }

	  Assert (col_index.first == current_block,
		  ExcInternalError());
	}
      local_row_length[current_block] = row_length;
      Assert (current_block < this->n_block_cols(),
	      ExcInternalError());

#ifdef DEBUG
				   // If in debug mode, do a check whether
				   // the right length has been obtained.
      unsigned int length = 0;
      for (unsigned int i=0; i<this->n_block_cols(); ++i)
	length += local_row_length[i];
      Assert (length == n_cols,
	      ExcDimensionMismatch(length, n_cols));
#endif
    }

				   // Now we found out about where the
				   // individual columns should start and
				   // where we should start reading out
				   // data. Now let's write the data into
				   // the individual blocks!
    for (unsigned int i=0; i<n_rows; ++i)
      {
	const std::pair<unsigned int,unsigned int> 
	  row_index = this->row_block_indices.global_to_local (row_indices[i]);
	for (unsigned int block_col=0; block_col<n_block_cols(); ++block_col)
	  {
	    if (local_row_length[block_col] == 0)
	      continue;

	    block(row_index.first, block_col).add 
	      (1, &row_index.second, 
	       local_row_length[block_col], 
	       &local_col_indices[block_col_indices[block_col]],
	       &values[n_cols*i + block_col_indices[block_col]]);
	  }
      }
  }



  TrilinosScalar
  BlockSparseMatrix::residual (MPI::BlockVector       &dst,
			       const MPI::BlockVector &x,
			       const MPI::BlockVector &b) const
  {
    vmult (dst, x);
    dst -= b;
    dst *= -1.;

    return dst.l2_norm();
  }



					// TODO: In the following we
					// use the same code as just
					// above six more times. Use
					// templates.
  TrilinosScalar
  BlockSparseMatrix::residual (BlockVector       &dst,
			       const BlockVector &x,
			       const BlockVector &b) const
  {
    vmult (dst, x);
    dst -= b;
    dst *= -1.;

    return dst.l2_norm();
  }



  TrilinosScalar
  BlockSparseMatrix::residual (MPI::BlockVector       &dst,
			       const MPI::Vector      &x,
			       const MPI::BlockVector &b) const
  {
    vmult (dst, x);
    dst -= b;
    dst *= -1.;

    return dst.l2_norm();
  }



  TrilinosScalar
  BlockSparseMatrix::residual (BlockVector       &dst,
			       const Vector      &x,
			       const BlockVector &b) const
  {
    vmult (dst, x);
    dst -= b;
    dst *= -1.;

    return dst.l2_norm();
  }



  TrilinosScalar
  BlockSparseMatrix::residual (MPI::Vector            &dst,
			       const MPI::BlockVector &x,
			       const MPI::Vector      &b) const
  {
    vmult (dst, x);
    dst -= b;
    dst *= -1.;

    return dst.l2_norm();
  }



  TrilinosScalar
  BlockSparseMatrix::residual (Vector            &dst,
			       const BlockVector &x,
			       const Vector      &b) const
  {
    vmult (dst, x);
    dst -= b;
    dst *= -1.;

    return dst.l2_norm();
  }



  TrilinosScalar
  BlockSparseMatrix::residual (VectorBase       &dst,
			       const VectorBase &x,
			       const VectorBase &b) const
  {
    vmult (dst, x);
    dst -= b;
    dst *= -1.;

    return dst.l2_norm();
  }





  // -------------------- explicit instantiations -----------------------
  //
  template void
  BlockSparseMatrix::reinit (const BlockSparsityPattern &);
  template void
  BlockSparseMatrix::reinit (const BlockCompressedSparsityPattern &);
  template void
  BlockSparseMatrix::reinit (const BlockCompressedSetSparsityPattern &);
  template void
  BlockSparseMatrix::reinit (const BlockCompressedSimpleSparsityPattern &);


  template void
  BlockSparseMatrix::reinit (const std::vector<Epetra_Map> &,
			     const BlockSparsityPattern    &);
  template void
  BlockSparseMatrix::reinit (const std::vector<Epetra_Map> &,
			     const BlockCompressedSparsityPattern &);
  template void
  BlockSparseMatrix::reinit (const std::vector<Epetra_Map> &,
			     const BlockCompressedSetSparsityPattern &);
  template void
  BlockSparseMatrix::reinit (const std::vector<Epetra_Map> &,
			     const BlockCompressedSimpleSparsityPattern &);

}


DEAL_II_NAMESPACE_CLOSE

#endif
