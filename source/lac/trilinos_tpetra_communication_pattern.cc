// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
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

#include <deal.II/lac/trilinos_tpetra_communication_pattern.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA

#  ifdef DEAL_II_WITH_MPI

#    include <deal.II/base/index_set.h>

#    include <Tpetra_Map.hpp>

#    include <memory>

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
    CommunicationPattern::CommunicationPattern(
      const IndexSet &vector_space_vector_index_set,
      const IndexSet &read_write_vector_index_set,
      const MPI_Comm &communicator)
    {
      // virtual functions called in constructors and destructors never use the
      // override in a derived class
      // for clarity be explicit on which function is called
      CommunicationPattern::reinit(vector_space_vector_index_set,
                                   read_write_vector_index_set,
                                   communicator);
    }



    void
    CommunicationPattern::reinit(const IndexSet &vector_space_vector_index_set,
                                 const IndexSet &read_write_vector_index_set,
                                 const MPI_Comm &communicator)
    {
      comm = std::make_shared<const MPI_Comm>(communicator);

      auto vector_space_vector_map =
        Teuchos::rcp(new Tpetra::Map<int, types::global_dof_index>(
          vector_space_vector_index_set.make_tpetra_map(*comm, false)));
      auto read_write_vector_map =
        Teuchos::rcp(new Tpetra::Map<int, types::global_dof_index>(
          read_write_vector_index_set.make_tpetra_map(*comm, true)));

      // Target map is read_write_vector_map
      // Source map is vector_space_vector_map. This map must have uniquely
      // owned GID.
      tpetra_import =
        std::make_unique<Tpetra::Import<int, types::global_dof_index>>(
          read_write_vector_map, vector_space_vector_map);
      tpetra_export =
        std::make_unique<Tpetra::Export<int, types::global_dof_index>>(
          read_write_vector_map, vector_space_vector_map);
    }



    const MPI_Comm &
    CommunicationPattern::get_mpi_communicator() const
    {
      return *comm;
    }



    const Tpetra::Import<int, types::global_dof_index> &
    CommunicationPattern::get_tpetra_import() const
    {
      return *tpetra_import;
    }



    const Tpetra::Export<int, types::global_dof_index> &
    CommunicationPattern::get_tpetra_export() const
    {
      return *tpetra_export;
    }
  } // namespace TpetraWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#  endif

#endif
