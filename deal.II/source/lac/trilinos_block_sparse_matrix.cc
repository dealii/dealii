//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <deal.II/lac/trilinos_block_sparse_matrix.h>

#ifdef DEAL_II_USE_TRILINOS

#  include <deal.II/lac/block_sparse_matrix.h>
#  include <deal.II/lac/block_sparsity_pattern.h>

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
  reinit (const std::vector<Epetra_Map> &parallel_partitioning,
	  const BlockSparsityType       &block_sparsity_pattern)
  {
    Assert (parallel_partitioning.size() == block_sparsity_pattern.n_block_rows(),
	    ExcDimensionMismatch (parallel_partitioning.size(),
				  block_sparsity_pattern.n_block_rows()));
    Assert (parallel_partitioning.size() == block_sparsity_pattern.n_block_cols(),
	    ExcDimensionMismatch (parallel_partitioning.size(),
				  block_sparsity_pattern.n_block_cols()));

    const unsigned int n_block_rows = parallel_partitioning.size();

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
	  this->sub_objects[r][c]->reinit (parallel_partitioning[r],
					   parallel_partitioning[c],
					   block_sparsity_pattern.block(r,c));
        }
  }



  template <typename BlockSparsityType>
  void
  BlockSparseMatrix::
  reinit (const std::vector<IndexSet> &parallel_partitioning,
	  const BlockSparsityType     &block_sparsity_pattern,
	  const MPI_Comm              &communicator)
  {
    std::vector<Epetra_Map> epetra_maps;
    for (unsigned int i=0; i<block_sparsity_pattern.n_block_rows(); ++i)
      epetra_maps.push_back
	(parallel_partitioning[i].make_trilinos_map(communicator, false));

    reinit (epetra_maps, block_sparsity_pattern);

  }



  template <typename BlockSparsityType>
  void
  BlockSparseMatrix::
  reinit (const BlockSparsityType &block_sparsity_pattern)
  {
    std::vector<Epetra_Map> parallel_partitioning;
    for (unsigned int i=0; i<block_sparsity_pattern.n_block_rows(); ++i)
      parallel_partitioning.push_back
	(Epetra_Map(block_sparsity_pattern.block(i,0).n_rows(),
		    0,
		    Utilities::Trilinos::comm_self()));

    reinit (parallel_partitioning, block_sparsity_pattern);
  }



  template <>
  void
  BlockSparseMatrix::
  reinit (const BlockSparsityPattern    &block_sparsity_pattern)
  {

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
	  this->sub_objects[r][c]->reinit (block_sparsity_pattern.block(r,c));
        }
  }



  void
  BlockSparseMatrix::
  reinit (const std::vector<Epetra_Map>             &parallel_partitioning,
	  const ::dealii::BlockSparseMatrix<double> &dealii_block_sparse_matrix,
	  const double                               drop_tolerance)
  {
    const unsigned int n_block_rows = parallel_partitioning.size();

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
          this->sub_objects[r][c]->reinit(parallel_partitioning[r],
					  parallel_partitioning[c],
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
    Epetra_MpiComm    trilinos_communicator (MPI_COMM_SELF);
#else
    Epetra_SerialComm trilinos_communicator;
#endif

    std::vector<Epetra_Map> parallel_partitioning;
    for (unsigned int i=0; i<dealii_block_sparse_matrix.n_block_rows(); ++i)
      parallel_partitioning.push_back (Epetra_Map(dealii_block_sparse_matrix.block(i,0).m(),
				       0,
				       trilinos_communicator));

    reinit (parallel_partitioning, dealii_block_sparse_matrix, drop_tolerance);
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
  BlockSparseMatrix::reinit (const dealii::BlockSparsityPattern &);
  template void
  BlockSparseMatrix::reinit (const dealii::BlockCompressedSparsityPattern &);
  template void
  BlockSparseMatrix::reinit (const dealii::BlockCompressedSetSparsityPattern &);
  template void
  BlockSparseMatrix::reinit (const dealii::BlockCompressedSimpleSparsityPattern &);


  template void
  BlockSparseMatrix::reinit (const std::vector<Epetra_Map> &,
			     const dealii::BlockSparsityPattern    &);
  template void
  BlockSparseMatrix::reinit (const std::vector<Epetra_Map> &,
			     const dealii::BlockCompressedSparsityPattern &);
  template void
  BlockSparseMatrix::reinit (const std::vector<Epetra_Map> &,
			     const dealii::BlockCompressedSetSparsityPattern &);
  template void
  BlockSparseMatrix::reinit (const std::vector<Epetra_Map> &,
			     const dealii::BlockCompressedSimpleSparsityPattern &);

}


DEAL_II_NAMESPACE_CLOSE

#endif
