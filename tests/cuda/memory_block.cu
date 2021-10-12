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


// Check MemoryBlock

#include <deal.II/base/memory_block.h>

#include <iostream>

#include "../tests.h"

struct MultiplyFunctor
{
  MultiplyFunctor(ArrayView<double, MemorySpace::CUDA> &memory_block_view)
    : mb_view(memory_block_view)
  {}

  DEAL_II_CUDA_HOST_DEV
  void
  operator()(int i) const
  {
    mb_view[i] *= 2.;
  }

  ArrayView<double, MemorySpace::CUDA> mb_view;
};

struct CheckFunctor
{
  CheckFunctor(ArrayView<double, MemorySpace::Host> &memory_block_view,
               std::vector<double>                   ref)
    : mb_view(memory_block_view)
    , reference(ref)
  {}

  void
  operator()(unsigned int i) const
  {
    AssertThrow(mb_view[i] == reference[i], ExcInternalError());
  }

  ArrayView<double, MemorySpace::Host> mb_view;
  std::vector<double>                  reference;
};

void
test()
{
  // Check that operator=0 works
  MemoryBlock<double, MemorySpace::Host> memory_block_host(10);
  memory_block_host = 0.;
  std::vector<double>                  zeros(10, 0.);
  ArrayView<double, MemorySpace::Host> mb_host_view(memory_block_host.data(),
                                                    memory_block_host.size());
  CheckFunctor                         check_functor_zero(mb_host_view, zeros);
  Utilities::MemorySpace::for_each_index(MemorySpace::Host{},
                                         memory_block_host.size(),
                                         check_functor_zero);

  // Check reinit
  std::vector<double> vector(5);
  for (unsigned int i = 0; i < vector.size(); ++i)
    vector[i] = 2. * i;
  ArrayView<double, MemorySpace::Host>   vector_view(vector);
  MemoryBlock<double, MemorySpace::CUDA> memory_block_dev(vector_view);
  AssertThrow(memory_block_dev.size() == vector.size(), ExcInternalError());
  memory_block_host.reinit(memory_block_dev);
  AssertThrow(memory_block_dev.size() == memory_block_host.size(),
              ExcInternalError());
  mb_host_view.reinit(memory_block_host.data(), memory_block_host.size());
  CheckFunctor check_functor_set(mb_host_view, vector);
  Utilities::MemorySpace::for_each_index(MemorySpace::Host{},
                                         memory_block_host.size(),
                                         check_functor_set);


  // Check access on GPU
  ArrayView<double, MemorySpace::CUDA> mb_dev_view(memory_block_dev.data(),
                                                   memory_block_dev.size());
  MultiplyFunctor                      multiply_functor(mb_dev_view);
  Utilities::MemorySpace::for_each_index(MemorySpace::CUDA{},
                                         memory_block_dev.size(),
                                         multiply_functor);
  memory_block_host.reinit(memory_block_dev);
  for (unsigned int i = 0; i < vector.size(); ++i)
    vector[i] *= 2.;
  mb_host_view.reinit(memory_block_host.data(), memory_block_host.size());
  CheckFunctor check_functor_multiply(mb_host_view, vector);
  Utilities::MemorySpace::for_each_index(MemorySpace::Host{},
                                         memory_block_host.size(),
                                         check_functor_multiply);

  deallog << "OK" << std::endl;
}

int
main(int argc, char *argv[])
{
  initlog();

  init_cuda();

  test();
  return 0;
}
