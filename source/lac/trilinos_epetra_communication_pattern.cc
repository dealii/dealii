// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/lac/trilinos_epetra_communication_pattern.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/base/index_set.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  include <Epetra_Map.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#  include <memory>

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace EpetraWrappers
  {
    CommunicationPattern::CommunicationPattern(
      const IndexSet &locally_owned_indices,
      const IndexSet &ghost_indices,
      const MPI_Comm  communicator)
    {
      // virtual functions called in constructors and destructors never use the
      // override in a derived class
      // for clarity be explicit on which function is called
      CommunicationPattern::reinit(locally_owned_indices,
                                   ghost_indices,
                                   communicator);
    }



    void
    CommunicationPattern::reinit(const IndexSet &locally_owned_indices,
                                 const IndexSet &ghost_indices,
                                 const MPI_Comm  communicator)
    {
      comm = std::make_shared<const MPI_Comm>(communicator);

      Epetra_Map vector_space_vector_map =
        locally_owned_indices.make_trilinos_map(*comm, false);
      Epetra_Map read_write_vector_map =
        ghost_indices.make_trilinos_map(*comm, true);

      // Target map is read_write_vector_map
      // Source map is vector_space_vector_map. This map must have uniquely
      // owned GID.
      importer = std::make_unique<Epetra_Import>(read_write_vector_map,
                                                 vector_space_vector_map);
    }



    MPI_Comm
    CommunicationPattern::get_mpi_communicator() const
    {
      return *comm;
    }



    const Epetra_Import &
    CommunicationPattern::get_epetra_import() const
    {
      return *importer;
    }
  } // namespace EpetraWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif
