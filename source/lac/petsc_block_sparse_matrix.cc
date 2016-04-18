// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2014 by the deal.II authors
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

#include <deal.II/lac/petsc_block_sparse_matrix.h>

#ifdef DEAL_II_WITH_PETSC

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
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
          this->sub_objects[r][c] = p;
        }
  }

  std::vector<BlockSparseMatrix::size_type >
  BlockSparseMatrix::
  locally_owned_domain_sizes() const
  {
    std::vector< size_type > index_sets;

    for ( unsigned int i=0; i<this->n_block_cols(); ++i)
      index_sets.push_back(this->block(0,i).n());

    return index_sets;
  }

  std::vector<BlockSparseMatrix::size_type >
  BlockSparseMatrix::
  locally_owned_range_sizes() const
  {
    std::vector< size_type > index_sets;

    for ( unsigned int i=0; i<this->n_block_rows(); ++i)
      index_sets.push_back(this->block(i,0).m());

    return index_sets;
  }

  void
  BlockSparseMatrix::collect_sizes ()
  {
    BaseClass::collect_sizes ();
  }

}


DEAL_II_NAMESPACE_CLOSE

#endif
