// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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


// check VectorTools::project_parallel() of quadrature data for
// FE_Q on a mesh with hanging nodes.

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
    test_with_hanging_nodes<LA::MPI::Vector>(FE_Q<dim>(p), p);
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();

  test<2>();
  test<3>();
}
