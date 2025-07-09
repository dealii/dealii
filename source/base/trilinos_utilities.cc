// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/mpi.h>
#include <deal.II/base/trilinos_utilities.h>

#ifdef DEAL_II_WITH_TRILINOS
DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  ifdef DEAL_II_WITH_MPI
#    include <Epetra_MpiComm.h>
#    include <Teuchos_DefaultComm.hpp>
#  endif
#  include <Epetra_SerialComm.h>
#  include <Teuchos_RCP.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS
#endif

DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
#ifdef DEAL_II_WITH_TRILINOS
  namespace Trilinos
  {
    const Epetra_Comm &
    comm_world()
    {
#  ifdef DEAL_II_WITH_MPI
      static Teuchos::RCP<Epetra_MpiComm> communicator =
        Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD), true);
#  else
      static Teuchos::RCP<Epetra_SerialComm> communicator =
        Teuchos::rcp(new Epetra_SerialComm(), true);
#  endif

      return *communicator;
    }



    const Teuchos::RCP<const Teuchos::Comm<int>> &
    tpetra_comm_self()
    {
#  ifdef DEAL_II_WITH_MPI
      static auto communicator = Teuchos::RCP<const Teuchos::Comm<int>>(
        new Teuchos::MpiComm<int>(MPI_COMM_SELF));
#  else
      static auto communicator =
        Teuchos::RCP<const Teuchos::Comm<int>>(new Teuchos::Comm<int>());
#  endif

      return communicator;
    }



    const Epetra_Comm &
    comm_self()
    {
#  ifdef DEAL_II_WITH_MPI
      static Teuchos::RCP<Epetra_MpiComm> communicator =
        Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_SELF), true);
#  else
      static Teuchos::RCP<Epetra_SerialComm> communicator =
        Teuchos::rcp(new Epetra_SerialComm(), true);
#  endif

      return *communicator;
    }



    Epetra_Comm *
    duplicate_communicator(const Epetra_Comm &communicator)
    {
#  ifdef DEAL_II_WITH_MPI

      // see if the communicator is in fact a
      // parallel MPI communicator; if so,
      // return a duplicate of it
      const Epetra_MpiComm *mpi_comm =
        dynamic_cast<const Epetra_MpiComm *>(&communicator);
      if (mpi_comm != nullptr)
        return new Epetra_MpiComm(
          Utilities::MPI::duplicate_communicator(mpi_comm->GetMpiComm()));
#  endif

      // if we don't support MPI, or if the
      // communicator in question was in fact
      // not an MPI communicator, return a
      // copy of the same object again
      Assert(dynamic_cast<const Epetra_SerialComm *>(&communicator) != nullptr,
             ExcInternalError());
      return new Epetra_SerialComm(
        dynamic_cast<const Epetra_SerialComm &>(communicator));
    }



    void
    destroy_communicator(Epetra_Comm &communicator)
    {
      // save the communicator, reset the map, and delete the communicator if
      // this whole thing was created as an MPI communicator
#  ifdef DEAL_II_WITH_MPI
      Epetra_MpiComm *mpi_comm = dynamic_cast<Epetra_MpiComm *>(&communicator);
      if (mpi_comm != nullptr)
        {
          MPI_Comm comm = mpi_comm->GetMpiComm();
          *mpi_comm     = Epetra_MpiComm(MPI_COMM_SELF);

          Utilities::MPI::free_communicator(comm);
        }
#  endif
    }



    unsigned int
    get_n_mpi_processes(const Epetra_Comm &mpi_communicator)
    {
      return mpi_communicator.NumProc();
    }


    unsigned int
    get_this_mpi_process(const Epetra_Comm &mpi_communicator)
    {
      return static_cast<unsigned int>(mpi_communicator.MyPID());
    }



    Epetra_Map
    duplicate_map(const Epetra_BlockMap &map, const Epetra_Comm &comm)
    {
      if (map.LinearMap() == true)
        {
          // each processor stores a
          // contiguous range of
          // elements in the
          // following constructor
          // call
          return Epetra_Map(map.NumGlobalElements(),
                            map.NumMyElements(),
                            map.IndexBase(),
                            comm);
        }
      else
        {
          // the range is not
          // contiguous
          return Epetra_Map(map.NumGlobalElements(),
                            map.NumMyElements(),
                            map.MyGlobalElements(),
                            0,
                            comm);
        }
    }
  } // namespace Trilinos
#endif


#ifdef DEAL_II_TRILINOS_WITH_TPETRA
  namespace Trilinos
  {
    MPI_Comm
    teuchos_comm_to_mpi_comm(
      const Teuchos::RCP<const Teuchos::Comm<int>> &teuchos_comm)
    {
      MPI_Comm out;
#  ifdef DEAL_II_WITH_MPI
      // Cast from Teuchos::Comm<int> to Teuchos::MpiComm<int>.
      const Teuchos::MpiComm<int> *mpi_comm =
        dynamic_cast<const Teuchos::MpiComm<int> *>(teuchos_comm.get());
      Assert(mpi_comm != nullptr, ExcInternalError());
      // From the Teuchos::MpiComm<int> object we can extract
      // the MPI_Comm object via the getRawMpiComm() function.
      out = *(mpi_comm->getRawMpiComm())();
#  else
      out = MPI_COMM_SELF;
#  endif
      return out;
    }
  }    // namespace Trilinos
#endif // DEAL_II_TRILINOS_WITH_TPETRA

} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE
