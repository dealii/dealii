// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// check VectorTools::project_parallel() for matrix-free quadrature data for
// FE_Q on a mesh with hanging nodes.

#include "project_parallel_qpmf_common.h"

template <int dim>
void test ()
{
  test_with_hanging_nodes<1, 2,dim> (FE_Q<dim>(1), 1);
  test_with_hanging_nodes<2, 3,dim> (FE_Q<dim>(2), 2);
  test_with_hanging_nodes<3, 4,dim> (FE_Q<dim>(3), 3);
}

int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv,
                                                       numbers::invalid_unsigned_int);

  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  test<2>();
  test<3>();
}

