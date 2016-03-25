// ---------------------------------------------------------------------
//
// Copyright (C) 2015-2016 by the deal.II authors
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

#ifndef dealii__trilinos_epetra_communication_pattern_h
#define dealii__trilinos_epetra_communication_pattern_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_TRILINOS

#ifdef DEAL_II_WITH_MPI

#include <deal.II/base/std_cxx11/shared_ptr.h>
#include <deal.II/lac/communication_pattern_base.h>
#include "Epetra_Import.h"

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace EpetraWrappers
  {
    /**
     * This class implements a wrapper to Trilinos Import.
     */
    class CommunicationPattern : public CommunicationPatternBase
    {
    public:
      /**
       * Reinitialize the communication pattern. The first argument @p
       * vector_space_vector_index_set is the index set associated to a
       * VectorSpaceVector object. The second argument @p
       * read_write_vector_index_set is the index set associated to a
       * ReadWriteVector object.
       */
      CommunicationPattern(const IndexSet &vector_space_vector_index_set,
                           const IndexSet &read_write_vector_index_set,
                           const MPI_Comm &communicator);

      /**
       * Reinitialize the object.
       */
      void reinit(const IndexSet &vector_space_vector_index_set,
                  const IndexSet &read_write_vector_index_set,
                  const MPI_Comm &communicator);

      /**
       * Return the underlying MPI communicator.
       */
      const MPI_Comm &get_mpi_communicator() const;

      /**
       * Return the underlying Epetra_Import object.
       */
      const Epetra_Import &get_epetra_import() const;

    private:
      /**
       * Shared pointer to the MPI communicator used.
       */
      std_cxx11::shared_ptr<const MPI_Comm> comm;

      /**
       * Shared pointer to the Epetra_Import object used.
       */
      std_cxx11::shared_ptr<Epetra_Import> import;
    };
  } // end of namespace EpetraWrappers
} // end of namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif

#endif

#endif
