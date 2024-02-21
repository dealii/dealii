// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/lac/la_parallel_vector.h>

#include "fe_tools_extrapolate_common.h"

// check FETools::extrapolate on distributed triangulations
// for LinearAlgebra::distributed::Vector<double>

template <int dim>
void
check(const FiniteElement<dim> &fe1,
      const FiniteElement<dim> &fe2,
      const std::string        &name)
{
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "Checking " << name << " in " << dim << "d:" << std::endl;

  // call main function in .cc files
  check_this_dealii<dim, LinearAlgebra::distributed::Vector<double>>(fe1, fe2);
}

#define CHECK(EL1, deg1, EL2, deg2, dim)                \
  {                                                     \
    FE_##EL1<dim> fe1(deg1);                            \
    FE_##EL2<dim> fe2(deg2);                            \
    check(fe1, fe2, #EL1 #deg1 " against " #EL2 #deg2); \
  }

#define CHECK_ALL(EL1, deg1, EL2, deg2) \
  {                                     \
    CHECK(EL1, deg1, EL2, deg2, 2);     \
    CHECK(EL1, deg1, EL2, deg2, 3);     \
  }


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll log;

  CHECK_ALL(Q, 1, Q, 1);
  CHECK_ALL(Q, 1, Q, 2);
  CHECK_ALL(Q, 2, Q, 1);
  CHECK_ALL(Q, 2, Q, 2);

  CHECK_ALL(DGQ, 0, DGQ, 0);
  CHECK_ALL(DGQ, 0, DGQ, 1);
  CHECK_ALL(DGQ, 1, DGQ, 0);
  CHECK_ALL(DGQ, 1, DGQ, 1);
}
