// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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


// check that we detect creating ArrayView objects wit the wrong memory space.

#include <deal.II/base/array_view.h>

#include "../tests.h"

int
main(int argc, char **argv)
{
  deal_II_exceptions::disable_abort_on_exception();

  initlog();

  init_cuda();

  std::vector<unsigned int>                                 dummy_host(2);
  std::unique_ptr<unsigned int[], void (*)(unsigned int *)> dummy_cuda(
    Utilities::CUDA::allocate_device_data<unsigned int>(2),
    Utilities::CUDA::delete_device_data<unsigned int>);

  deallog << "Testing host ArrayView with host memory" << std::endl;
  ArrayView<unsigned int, MemorySpace::Host> view_1(dummy_host);

  deallog << "Testing device ArrayView with host memory" << std::endl;
  try
    {
      ArrayView<unsigned int, MemorySpace::CUDA> view_2(dummy_host);
    }
  catch (const ExceptionBase &exc)
    {
      deallog << exc.what() << std::endl;
    }

  deallog << "Testing host ArrayView with device memory" << std::endl;
  try
    {
      ArrayView<unsigned int, MemorySpace::Host> view_3(dummy_cuda.get(), 2);
    }
  catch (const ExceptionBase &exc)
    {
      deallog << exc.what() << std::endl;
    }

  deallog << "Testing device ArrayView with device memory" << std::endl;
  ArrayView<unsigned int, MemorySpace::CUDA> view_4(dummy_cuda.get(), 2);

  return 0;
}
