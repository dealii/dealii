// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2022 by the deal.II authors
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

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/tools.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>

#include <deal.II/numerics/vector_tools.h>

#include "ecl.h"


// Like ecl_01.cc but running with shared-memory MPI enabled.

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  mpi_initlog();

  // test without shared memory
  {
    test<2, 2, 3, double, VectorizedArray<double>>(0, 6, false, MPI_COMM_SELF);
  }

  // test with shared memory
  {
    const auto rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

    MPI_Comm subcommunicator;
    MPI_Comm_split_type(MPI_COMM_WORLD,
                        MPI_COMM_TYPE_SHARED,
                        rank,
                        MPI_INFO_NULL,
                        &subcommunicator);

    test<2, 2, 3, double, VectorizedArray<double>>(0,
                                                   6,
                                                   false,
                                                   subcommunicator);

    MPI_Comm_free(&subcommunicator);
  }
}
