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
    template <typename MemorySpace>
    CommunicationPattern<MemorySpace>::CommunicationPattern(
      const IndexSet &locally_owned_indices,
      const IndexSet &ghost_indices,
      const MPI_Comm  communicator)
    {
      // virtual functions called in constructors and destructors never use the
      // override in a derived class
      // for clarity be explicit on which function is called
      CommunicationPattern<MemorySpace>::reinit(locally_owned_indices,
                                                ghost_indices,
                                                communicator);
    }



    template <typename MemorySpace>
    void
    CommunicationPattern<MemorySpace>::reinit(
      const IndexSet &locally_owned_indices,
      const IndexSet &ghost_indices,
      const MPI_Comm  communicator)
    {
      comm = Teuchos::rcpFromUndefRef(communicator);

      auto vector_space_vector_map =
        locally_owned_indices
          .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(
            *comm, false);
      auto read_write_vector_map =
        ghost_indices
          .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(
            *comm, true);

      // Target map is read_write_vector_map
      // Source map is vector_space_vector_map. This map must have uniquely
      // owned GID.
      tpetra_import =
        Teuchos::rcp(new Tpetra::Import<int,
                                        types::signed_global_dof_index,
                                        TpetraTypes::NodeType<MemorySpace>>(
          read_write_vector_map, vector_space_vector_map));
      tpetra_export =
        Teuchos::rcp(new Tpetra::Export<int,
                                        types::signed_global_dof_index,
                                        TpetraTypes::NodeType<MemorySpace>>(
          read_write_vector_map, vector_space_vector_map));
    }



    template <typename MemorySpace>
    MPI_Comm
    CommunicationPattern<MemorySpace>::get_mpi_communicator() const
    {
      return *comm;
    }



    template <typename MemorySpace>
    const Tpetra::Import<int,
                         types::signed_global_dof_index,
                         TpetraTypes::NodeType<MemorySpace>> &
    CommunicationPattern<MemorySpace>::get_tpetra_import() const
    {
      return *tpetra_import;
    }


    template <typename MemorySpace>
    Teuchos::RCP<TpetraTypes::ImportType<MemorySpace>>
    CommunicationPattern<MemorySpace>::get_tpetra_import_rcp() const
    {
      return tpetra_import;
    }



    template <typename MemorySpace>
    const Tpetra::Export<int,
                         types::signed_global_dof_index,
                         TpetraTypes::NodeType<MemorySpace>> &
    CommunicationPattern<MemorySpace>::get_tpetra_export() const
    {
      return *tpetra_export;
    }



    template <typename MemorySpace>
    Teuchos::RCP<TpetraTypes::ExportType<MemorySpace>>
    CommunicationPattern<MemorySpace>::get_tpetra_export_rcp() const
    {
      return tpetra_export;
    }
  } // namespace TpetraWrappers
} // namespace LinearAlgebra


template class LinearAlgebra::TpetraWrappers::CommunicationPattern<
  MemorySpace::Default>;
template class LinearAlgebra::TpetraWrappers::CommunicationPattern<
  MemorySpace::Host>;

DEAL_II_NAMESPACE_CLOSE

#endif
