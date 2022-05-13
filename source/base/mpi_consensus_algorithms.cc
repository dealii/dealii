// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2021 by the deal.II authors
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

#include <boost/serialization/utility.hpp>


DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  namespace MPI
  {
    namespace ConsensusAlgorithms
    {
      template class Process<std::vector<unsigned int>,
                             std::vector<unsigned int>>;

      template class Interface<std::vector<unsigned int>,
                               std::vector<unsigned int>>;

      template class NBX<std::vector<unsigned int>, std::vector<unsigned int>>;

      template class PEX<std::vector<unsigned int>, std::vector<unsigned int>>;

      template class Serial<std::vector<unsigned int>,
                            std::vector<unsigned int>>;

      template class Selector<std::vector<unsigned int>,
                              std::vector<unsigned int>>;


      template class Process<std::vector<std::pair<types::global_dof_index,
                                                   types::global_dof_index>>,
                             std::vector<unsigned int>>;

      template class Interface<std::vector<std::pair<types::global_dof_index,
                                                     types::global_dof_index>>,
                               std::vector<unsigned int>>;

      template class Selector<std::vector<std::pair<types::global_dof_index,
                                                    types::global_dof_index>>,
                              std::vector<unsigned int>>;

      template class NBX<std::vector<std::pair<types::global_dof_index,
                                               types::global_dof_index>>,
                         std::vector<unsigned int>>;

      template class Serial<std::vector<std::pair<types::global_dof_index,
                                                  types::global_dof_index>>,
                            std::vector<unsigned int>>;

      template class PEX<std::vector<std::pair<types::global_dof_index,
                                               types::global_dof_index>>,
                         std::vector<unsigned int>>;

#ifdef DEAL_II_WITH_64BIT_INDICES
      template class Process<std::vector<types::global_dof_index>,
                             std::vector<unsigned int>>;

      template class Interface<std::vector<types::global_dof_index>,
                               std::vector<unsigned int>>;

      template class NBX<std::vector<types::global_dof_index>,
                         std::vector<unsigned int>>;

      template class Serial<std::vector<types::global_dof_index>,
                            std::vector<unsigned int>>;

      template class PEX<std::vector<types::global_dof_index>,
                         std::vector<unsigned int>>;

      template class Selector<std::vector<types::global_dof_index>,
                              std::vector<unsigned int>>;
#endif

      template class Process<std::vector<char>, std::vector<char>>;

      template class Interface<std::vector<char>, std::vector<char>>;

      template class NBX<std::vector<char>, std::vector<char>>;

      template class PEX<std::vector<char>, std::vector<char>>;

      template class Serial<std::vector<char>, std::vector<char>>;

      template class Selector<std::vector<char>, std::vector<char>>;

    } // namespace ConsensusAlgorithms
  }   // end of namespace MPI
} // end of namespace Utilities

DEAL_II_NAMESPACE_CLOSE
