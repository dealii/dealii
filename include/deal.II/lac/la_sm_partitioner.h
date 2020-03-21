// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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

#ifndef dealii_la_sm_partitioner_h
#define dealii_la_sm_partitioner_h


#include <deal.II/base/config.h>

#include <deal.II/lac/communication_pattern_base.h>

DEAL_II_NAMESPACE_OPEN


namespace LinearAlgebra
{
  namespace SharedMPI
  {
    template <typename Number = double>
    class Partitioner : public LinearAlgebra::CommunicationPatternBase
    {
    public:
      Partitioner(const MPI_Comm &comm, const MPI_Comm &comm_sm)
        : comm(comm)
        , comm_sm(comm_sm)
      {}

      const MPI_Comm &
      get_mpi_communicator() const override
      {
        return comm;
      }

      void
      reinit(const IndexSet &indexset_has,
             const IndexSet &indexset_want,
             const MPI_Comm &communicator) override
      {
        (void)indexset_has;
        (void)indexset_want;
        (void)communicator;
      }

    private:
      const MPI_Comm &comm;
      const MPI_Comm &comm_sm;
    };

  } // end of namespace SharedMPI
} // end of namespace LinearAlgebra


DEAL_II_NAMESPACE_CLOSE

#endif