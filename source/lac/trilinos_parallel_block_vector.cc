// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II authors
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


#include <deal.II/lac/trilinos_parallel_block_vector.h>

#ifdef DEAL_II_WITH_TRILINOS

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  namespace MPI
  {
    BlockVector::BlockVector (const std::vector<Epetra_Map> &parallel_partitioning)
    {
      reinit (parallel_partitioning, false);
    }



    bool
    BlockVector::is_compressed () const
    {
      bool compressed = true;
      for (unsigned int row=0; row<n_blocks(); ++row)
        if (block(row).is_compressed() == false)
          {
            compressed = false;
            break;
          }

      return compressed;
    }
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS
