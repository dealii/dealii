// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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


// Test that we can assign LA::d::Vector from one MemorySpace to another one


#include <deal.II/base/exceptions.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <cmath>

#include "../tests.h"


using namespace dealii;

template <typename Number>
void
test()
{
  const unsigned int                                            size = 100;
  LinearAlgebra::distributed::Vector<Number, MemorySpace::Host> vec_ref(size);
  for (unsigned int i = 0; i < size; ++i)
    vec_ref[i] = i;

  // Assignment from Host to CUDA
  LinearAlgebra::distributed::Vector<Number, MemorySpace::CUDA> vec_dev;
  vec_dev = vec_ref;
  vec_dev *= 2;

  // Assignment from CUDA to HOST
  LinearAlgebra::distributed::Vector<Number, MemorySpace::Host> vec_host;
  vec_host = vec_dev;

  vec_ref *= 2;

  for (unsigned int i = 0; i < size; ++i)
    AssertThrow(std::fabs(vec_ref[i] - vec_host[i]) < 1e-12,
                ExcInternalError());
}


int
main(int argc, char *argv[])
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  test<float>();
  test<double>();

  deallog << "OK" << std::endl;

  return 0;
}
