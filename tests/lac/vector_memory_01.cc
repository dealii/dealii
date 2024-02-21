// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test that we can successfully fill a GrowingVectorMemory pool
// with LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>
// objects. Partially copied from lac/vector_memory.cc


#include <deal.II/base/exceptions.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector_memory.h>

#include "../tests.h"



template <typename VectorType>
void
test_stat()
{
  GrowingVectorMemory<VectorType> mem(1, true);
  VectorType                     *v1 = mem.alloc();
  VectorType                     *v2 = mem.alloc();
  VectorType                     *v3 = mem.alloc();
  VectorType                     *v4 = mem.alloc();
  VectorType                     *v5 = mem.alloc();
  VectorType                     *v6 = mem.alloc();
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

  test_stat<LinearAlgebra::distributed::Vector<double, MemorySpace::Default>>();
  test_stat<LinearAlgebra::distributed::Vector<float, MemorySpace::Default>>();

  return 0;
}
