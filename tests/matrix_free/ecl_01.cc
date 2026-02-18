// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2021 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

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


// Compare the evaluation of a SIP Laplace operator with face- and element-
// centric loops on a hypercube, hypershell, and hyperball.

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  mpi_initlog();
  test<2, 2, 3, double, VectorizedArray<double>>(0);
  test<2, 2, 3, double, VectorizedArray<double>>(1);
  test<2, 2, 3, double, VectorizedArray<double>>(2);
}
