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


// check VectorTools::project_parallel() of quadrature data for
// FE_Q on a mesh without hanging nodes.

#include "../tests.h"

#include "project_parallel_qp_common.h"

namespace LA
{
  using namespace ::LinearAlgebraTrilinos;
}

template <int dim>
void
test()
{
  for (unsigned int p = 1; p < 7 - dim; ++p)
    test_no_hanging_nodes<LA::MPI::Vector>(FE_Q<dim>(p), p);
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();

  test<2>();
  test<3>();
}
