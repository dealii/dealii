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

#ifndef dealii_trilinos_tpetra_communication_pattern_h
#define dealii_trilinos_tpetra_communication_pattern_h


#include <deal.II/base/config.h>

#if defined(DEAL_II_TRILINOS_WITH_TPETRA) && defined(DEAL_II_WITH_MPI)

#  include <deal.II/base/communication_pattern_base.h>

#  include <Tpetra_Export.hpp>
#  include <Tpetra_Import.hpp>

#  include <memory>

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
    /**
     * This class implements a wrapper to Tpetra::Import and Tpetra::Export.
     */
    class CommunicationPattern : public Utilities::MPI::CommunicationPatternBase
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
       * Return the underlying Tpetra::Import object.
       */
      const Tpetra::Import<int, types::global_dof_index> &
      get_tpetra_import() const;

      /**
       * Return the underlying Tpetra::Export object.
       */
      const Tpetra::Export<int, types::global_dof_index> &
      get_tpetra_export() const;

    private:
      /**
       * Shared pointer to the MPI communicator used.
       */
      std::shared_ptr<const MPI_Comm> comm;

      /**
       * Shared pointer to the Tpetra::Import object used.
       */
      std::unique_ptr<Tpetra::Import<int, types::global_dof_index>>
        tpetra_import;

      /**
       * Shared pointer to the Tpetra::Export object used.
       */
      std::unique_ptr<Tpetra::Export<int, types::global_dof_index>>
        tpetra_export;
    };
  } // end of namespace TpetraWrappers
} // end of namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif

#endif
