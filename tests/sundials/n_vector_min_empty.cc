// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Check N_VMin when one MPI rank owns no entries. Test both choices of
// owning rank and positive and negative minima.

#include <deal.II/base/mpi.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/sundials/n_vector.h>
#include <deal.II/sundials/n_vector.templates.h>

#include "../tests.h"

#if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
static SUNContext context;
#endif


template <typename VectorType>
void
check()
{
  const auto rank    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const auto n_ranks = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  for (unsigned int owner = 0; owner < n_ranks; ++owner)
    {
      IndexSet owned(3);
      if (rank == owner)
        owned.add_range(0, 3);
      VectorType vector(owned, MPI_COMM_WORLD);
      for (const double minimum : {2.0, -3.0})
        {
          if (rank == owner)
            {
              vector[0] = minimum + 4.0;
              vector[1] = minimum;
              vector[2] = minimum + 1.0;
            }
          vector.compress(VectorOperation::insert);
          auto view = SUNDIALS::internal::make_nvector_view(vector
#if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                                                            ,
                                                            context
#endif
          );
          AssertThrow(N_VMin(view) == minimum, ExcInternalError());
        }
    }
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    initlog();
#if DEAL_II_SUNDIALS_VERSION_GTE(7, 0, 0)
  AssertThrow(SUNContext_Create(MPI_COMM_WORLD, &context) == 0,
              ExcInternalError());
#elif DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
  MPI_Comm communicator = MPI_COMM_WORLD;
  AssertThrow(SUNContext_Create(&communicator, &context) == 0,
              ExcInternalError());
#endif

#ifdef DEAL_II_WITH_PETSC
  check<PETScWrappers::MPI::Vector>();
#endif
  check<LinearAlgebra::distributed::Vector<double>>();
#ifdef DEAL_II_WITH_TRILINOS
  check<TrilinosWrappers::MPI::Vector>();
#endif

#if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
  SUNContext_Free(&context);
#endif
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
}
