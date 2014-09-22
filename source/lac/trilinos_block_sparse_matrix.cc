// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/lac/trilinos_block_sparse_matrix.h>

#ifdef DEAL_II_WITH_TRILINOS

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
  reinit (const size_type n_block_rows,
          const size_type n_block_columns)
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
    for (size_type r=0; r<this->n_block_rows(); ++r)
      for (size_type c=0; c<this->n_block_cols(); ++c)
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
          const BlockSparsityType       &block_sparsity_pattern,
          const bool                   exchange_data)
  {
    Assert (parallel_partitioning.size() == block_sparsity_pattern.n_block_rows(),
            ExcDimensionMismatch (parallel_partitioning.size(),
                                  block_sparsity_pattern.n_block_rows()));
    Assert (parallel_partitioning.size() == block_sparsity_pattern.n_block_cols(),
            ExcDimensionMismatch (parallel_partitioning.size(),
                                  block_sparsity_pattern.n_block_cols()));

    const size_type n_block_rows = parallel_partitioning.size();

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
    for (size_type r=0; r<this->n_block_rows(); ++r)
      for (size_type c=0; c<this->n_block_cols(); ++c)
        {
          this->sub_objects[r][c]->reinit (parallel_partitioning[r],
                                           parallel_partitioning[c],
                                           block_sparsity_pattern.block(r,c),
                                           exchange_data);
        }
  }



  template <typename BlockSparsityType>
  void
  BlockSparseMatrix::
  reinit (const std::vector<IndexSet> &parallel_partitioning,
          const BlockSparsityType     &block_sparsity_pattern,
          const MPI_Comm              &communicator,
          const bool                   exchange_data)
  {
    std::vector<Epetra_Map> epetra_maps;
    for (size_type i=0; i<block_sparsity_pattern.n_block_rows(); ++i)
      epetra_maps.push_back
      (parallel_partitioning[i].make_trilinos_map(communicator, false));

    reinit (epetra_maps, block_sparsity_pattern, exchange_data);

  }



  template <typename BlockSparsityType>
  void
  BlockSparseMatrix::
  reinit (const BlockSparsityType &block_sparsity_pattern)
  {
    std::vector<Epetra_Map> parallel_partitioning;
    for (size_type i=0; i<block_sparsity_pattern.n_block_rows(); ++i)
      parallel_partitioning.push_back
      (Epetra_Map(static_cast<TrilinosWrappers::types::int_type>(block_sparsity_pattern.block(i,0).n_rows()),
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
    for (size_type r=0; r<this->n_block_rows(); ++r)
      for (size_type c=0; c<this->n_block_cols(); ++c)
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
    const size_type n_block_rows = parallel_partitioning.size();

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
    for (size_type r=0; r<this->n_block_rows(); ++r)
      for (size_type c=0; c<this->n_block_cols(); ++c)
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
#ifdef DEAL_II_WITH_MPI
    Epetra_MpiComm    trilinos_communicator (MPI_COMM_SELF);
#else
    Epetra_SerialComm trilinos_communicator;
#endif

    std::vector<Epetra_Map> parallel_partitioning;
    for (size_type i=0; i<dealii_block_sparse_matrix.n_block_rows(); ++i)
      parallel_partitioning.push_back (Epetra_Map(static_cast<TrilinosWrappers::types::int_type>(dealii_block_sparse_matrix.block(i,0).m()),
                                                  0,
                                                  trilinos_communicator));

    reinit (parallel_partitioning, dealii_block_sparse_matrix, drop_tolerance);
  }





  void
  BlockSparseMatrix::collect_sizes ()
  {
    compress();
    BaseClass::collect_sizes ();
  }



  BlockSparseMatrix::size_type
  BlockSparseMatrix::n_nonzero_elements () const
  {
    size_type n_nonzero = 0;
    for (size_type rows = 0; rows<this->n_block_rows(); ++rows)
      for (size_type cols = 0; cols<this->n_block_cols(); ++cols)
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
                             const dealii::BlockSparsityPattern &,
                             const bool);
  template void
  BlockSparseMatrix::reinit (const std::vector<Epetra_Map> &,
                             const dealii::BlockCompressedSparsityPattern &,
                             const bool);
  template void
  BlockSparseMatrix::reinit (const std::vector<Epetra_Map> &,
                             const dealii::BlockCompressedSetSparsityPattern &,
                             const bool);
  template void
  BlockSparseMatrix::reinit (const std::vector<Epetra_Map> &,
                             const dealii::BlockCompressedSimpleSparsityPattern &,
                             const bool);

  template void
  BlockSparseMatrix::reinit (const std::vector<IndexSet> &,
                             const dealii::BlockCompressedSimpleSparsityPattern &,
                             const MPI_Comm &,
                             const bool);

}


DEAL_II_NAMESPACE_CLOSE

#endif
