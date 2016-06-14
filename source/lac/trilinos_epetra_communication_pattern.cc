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

#include <deal.II/lac/trilinos_epetra_communication_pattern.h>

#ifdef DEAL_II_WITH_TRILINOS

#ifdef DEAL_II_WITH_MPI

#include <deal.II/base/index_set.h>
#include <deal.II/base/std_cxx11/shared_ptr.h>
#include "Epetra_Map.h"

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace EpetraWrappers
  {
    CommunicationPattern::CommunicationPattern(const IndexSet &vector_space_vector_index_set,
                                               const IndexSet &read_write_vector_index_set,
                                               const MPI_Comm &communicator)
    {
      reinit(vector_space_vector_index_set, read_write_vector_index_set, communicator);
    }



    void CommunicationPattern::reinit(const IndexSet &vector_space_vector_index_set,
                                      const IndexSet &read_write_vector_index_set,
                                      const MPI_Comm &communicator)
    {
      comm = std_cxx11::make_shared<const MPI_Comm>(communicator);

      Epetra_Map vector_space_vector_map = vector_space_vector_index_set.make_trilinos_map(*comm,
                                           false);
      Epetra_Map read_write_vector_map = read_write_vector_index_set.make_trilinos_map(*comm,
                                         true);

      // Target map is read_write_vector_map
      // Source map is vector_space_vector_map. This map must have uniquely
      // owned GID.
      import.reset(new Epetra_Import(read_write_vector_map, vector_space_vector_map));
    }



    const MPI_Comm &CommunicationPattern::get_mpi_communicator() const
    {
      return *comm;
    }



    const Epetra_Import &CommunicationPattern::get_epetra_import() const
    {
      return *import;
    }
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif

#endif
