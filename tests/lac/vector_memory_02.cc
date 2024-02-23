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


// Test that memory leaks are detected correctly for a GrowingVectorMemory pool
// with LinearAlgebra::distributed<Number, MemorySpace::Default> objects.
// Partially copied from lac/vector_memory.cc


#include <deal.II/base/exceptions.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector_memory.h>

#include "../tests.h"



template <typename VectorType>
void
test_leak()
{
  GrowingVectorMemory<VectorType> mem;
  VectorType                     *v = mem.alloc();
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
        LinearAlgebra::distributed::Vector<double, MemorySpace::Default>>();
      test_leak<
        LinearAlgebra::distributed::Vector<float, MemorySpace::Default>>();
    }
  catch (const StandardExceptions::ExcMemoryLeak &e)
    {
      deallog << "Exception: " << e.get_exc_name() << std::endl;
    }

  return 0;
}
