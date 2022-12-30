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


// Test that memory leaks are detected correctly for a GrowingVectorMemory pool
// with LinearAlgebra::distributed<Number, MemorySpace::CUDA> objects.
// Partially copied from lac/vector_memory.cc


#include <deal.II/base/exceptions.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector_memory.h>

#include "../tests.h"


using namespace dealii;

template <typename VectorType>
void
test_leak()
{
  GrowingVectorMemory<VectorType> mem;
  VectorType *                    v = mem.alloc();
  v->reinit(5);
}

int
main(int argc, char *argv[])
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);
  deal_II_exceptions::disable_abort_on_exception();

  try
    {
      test_leak<
        LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA>>();
      test_leak<LinearAlgebra::distributed::Vector<float, MemorySpace::CUDA>>();
    }
  catch (const StandardExceptions::ExcMemoryLeak &e)
    {
      deallog << "Exception: " << e.get_exc_name() << std::endl;
    }

  return 0;
}
