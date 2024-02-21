// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2020 by the deal.II authors
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
// FE_Q on a mesh without hanging nodes.


#include "../tests.h"

#include "project_parallel_qpmf_common.h"

template <int dim>
void
test()
{
  test_no_hanging_nodes<1, 2, dim>(FE_Q<dim>(1), 1);
  test_no_hanging_nodes<2, 3, dim>(FE_Q<dim>(2), 2);
  test_no_hanging_nodes<3, 4, dim>(FE_Q<dim>(3), 3);
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
