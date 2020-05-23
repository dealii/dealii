// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/base/mpi_consensus_algorithms.templates.h>

DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  namespace MPI
  {
    namespace ConsensusAlgorithms
    {
      template class Process<unsigned int, unsigned int>;

      template class NBX<unsigned int, unsigned int>;

      template class PEX<unsigned int, unsigned int>;

      template class Selector<unsigned int, unsigned int>;


      template class Process<
        std::pair<types::global_dof_index, types::global_dof_index>,
        unsigned int>;

      template class Selector<
        std::pair<types::global_dof_index, types::global_dof_index>,
        unsigned int>;

      template class NBX<
        std::pair<types::global_dof_index, types::global_dof_index>,
        unsigned int>;

      template class PEX<
        std::pair<types::global_dof_index, types::global_dof_index>,
        unsigned int>;

    } // namespace ConsensusAlgorithms
  }   // end of namespace MPI
} // end of namespace Utilities

DEAL_II_NAMESPACE_CLOSE
