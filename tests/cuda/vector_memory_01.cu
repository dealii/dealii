// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2017 by the deal.II authors
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


// Test that we can successfully fill a GrowingVectorMemory pool
// with LinearAlgebra::distributed::Vector<Number, MemorySpace::CUDA> objects.
// Partially copied from lac/vector_memory.cc


#include <deal.II/base/exceptions.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector_memory.h>

#include "../tests.h"


using namespace dealii;

template <typename VectorType>
void
test_stat()
{
  GrowingVectorMemory<VectorType> mem(1, true);
  VectorType *                    v1 = mem.alloc();
  VectorType *                    v2 = mem.alloc();
  VectorType *                    v3 = mem.alloc();
  VectorType *                    v4 = mem.alloc();
  VectorType *                    v5 = mem.alloc();
  VectorType *                    v6 = mem.alloc();
  v1->reinit(5);
  v2->reinit(5);
  v3->reinit(5);
  v4->reinit(5);
  v5->reinit(5);
  v6->reinit(5);
  mem.free(v1);
  mem.free(v2);
  mem.free(v3);
  mem.free(v4);
  mem.free(v5);
  mem.free(v6);
  v1 = mem.alloc();
  mem.free(v1);
  v1 = mem.alloc();
  mem.free(v1);
  v1 = mem.alloc();
  mem.free(v1);
  v1 = mem.alloc();
  mem.free(v1);
}


int
main(int argc, char *argv[])
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  test_stat<LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA>>();
  test_stat<LinearAlgebra::distributed::Vector<float, MemorySpace::CUDA>>();

  return 0;
}
