// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include "deal.II/base/memory_space.h"

#include "deal.II/lac/trilinos_tpetra_types.h"
#include <deal.II/lac/trilinos_tpetra_communication_pattern.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA

#  include <deal.II/base/index_set.h>

#  include <Tpetra_Map.hpp>

#  include <memory>

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace TpetraWrappers
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
      comm = Teuchos::rcpFromUndefRef(communicator);

      auto vector_space_vector_map =
        locally_owned_indices.make_tpetra_map_rcp(*comm, false);
      auto read_write_vector_map =
        ghost_indices.make_tpetra_map_rcp(*comm, true);

      // Target map is read_write_vector_map
      // Source map is vector_space_vector_map. This map must have uniquely
      // owned GID.
      tpetra_import =
        Teuchos::rcp(new Tpetra::Import<int, types::signed_global_dof_index>(
          read_write_vector_map, vector_space_vector_map));
      tpetra_export =
        Teuchos::rcp(new Tpetra::Export<int, types::signed_global_dof_index>(
          read_write_vector_map, vector_space_vector_map));
    }



    MPI_Comm
    CommunicationPattern::get_mpi_communicator() const
    {
      return *comm;
    }



    const Tpetra::Import<int, types::signed_global_dof_index> &
    CommunicationPattern::get_tpetra_import() const
    {
      return *tpetra_import;
    }



    Teuchos::RCP<TpetraTypes::ImportType<dealii::MemorySpace::Host>>
    CommunicationPattern::get_tpetra_import_rcp() const
    {
      return tpetra_import;
    }



    const Tpetra::Export<int, types::signed_global_dof_index> &
    CommunicationPattern::get_tpetra_export() const
    {
      return *tpetra_export;
    }



    Teuchos::RCP<TpetraTypes::ExportType<dealii::MemorySpace::Host>>
    CommunicationPattern::get_tpetra_export_rcp() const
    {
      return tpetra_export;
    }
  } // namespace TpetraWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif
