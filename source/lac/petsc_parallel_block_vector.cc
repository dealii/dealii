// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
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

#include <deal.II/lac/petsc_parallel_block_vector.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/lac/petsc_block_vector.h>

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  namespace MPI
  {
    typedef types::global_dof_index size_type;

    BlockVector &
    BlockVector::operator = (const PETScWrappers::BlockVector &v)
    {
      Assert (v.get_block_indices() == this->get_block_indices(),
              ExcNonMatchingBlockVectors());

      for (unsigned int i=0; i<this->n_blocks(); ++i)
        this->block(i) = v.block(i);

      return *this;
    }


    void
    BlockVector::reinit (const unsigned int num_blocks)
    {
      std::vector<size_type> block_sizes (num_blocks, 0);
      this->block_indices.reinit (block_sizes);
      if (this->components.size() != this->n_blocks())
        this->components.resize(this->n_blocks());

      for (unsigned int i=0; i<this->n_blocks(); ++i)
        components[i].reinit (MPI_COMM_SELF, 0, 0);

      collect_sizes();
    }
  }

}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
