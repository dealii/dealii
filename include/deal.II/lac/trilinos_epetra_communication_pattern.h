// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II authors
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

#ifndef dealii_trilinos_epetra_communication_pattern_h
#define dealii_trilinos_epetra_communication_pattern_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_TRILINOS

#  ifdef DEAL_II_WITH_MPI

#    include <deal.II/base/communication_pattern_base.h>

#    include <Epetra_Import.h>

#    include <memory>

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace EpetraWrappers
  {
    /**
     * This class implements a wrapper to a Trilinos Epetra_Import object,
     * for use in places where a Utilities::MPI::CommunicationPatternBase object
     * is required.
     */
    class CommunicationPattern : public Utilities::MPI::CommunicationPatternBase
    {
    public:
      /**
       * Initialize the communication pattern. The first argument @p
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
      virtual void
      reinit(const IndexSet &vector_space_vector_index_set,
             const IndexSet &read_write_vector_index_set,
             const MPI_Comm &communicator) override;

      /**
       * Return the underlying MPI communicator.
       */
      virtual const MPI_Comm &
      get_mpi_communicator() const override;

      /**
       * Return the underlying Epetra_Import object.
       */
      const Epetra_Import &
      get_epetra_import() const;

    private:
      /**
       * Shared pointer to the MPI communicator used.
       */
      std::shared_ptr<const MPI_Comm> comm;

      /**
       * Shared pointer to the Epetra_Import object used.
       */
      std::unique_ptr<Epetra_Import> importer;
    };
  } // end of namespace EpetraWrappers
} // end of namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#  endif

#endif

#endif
