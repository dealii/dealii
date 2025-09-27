// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check VectorTools::project_parallel() for matrix-free quadrature data for
// FE_Q on a mesh with hanging nodes and specified fe_component

#include "../tests.h"

#include "project_parallel_qpmf_common.h"

template <int dim>
void
test()
{
  FE_Q<dim> fe1(1);
  FE_Q<dim> fe2(2);
  FE_Q<dim> fe4(4);
  test_with_hanging_nodes<1, 2, dim>({{&fe1, &fe2}}, 1, 0); // take first
  test_with_hanging_nodes<2, 3, dim>({{&fe1, &fe2}}, 2, 1); // take second
  test_with_hanging_nodes<2, 6, dim>(
    {{&fe2, &fe4}}, 2, 0); // take first with non-standard n_q_points_1d
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);

  initlog();

  test<2>();
  test<3>();
}
