// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
