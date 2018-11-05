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


// check that we detect that accessing CUDA memory in an ArrayView object
// is not allowed.

#include <deal.II/base/array_view.h>

#include "../tests.h"

int
main(int argc, char **argv)
{
  deal_II_exceptions::disable_abort_on_exception();

  initlog();

  init_cuda();

  std::unique_ptr<unsigned int[], void (*)(unsigned int *)> dummy_cuda(
    Utilities::CUDA::allocate_device_data<unsigned int>(2),
    Utilities::CUDA::delete_device_data<unsigned int>);

  try
    {
      ArrayView<unsigned int, MemorySpace::CUDA> view(dummy_cuda.get(), 2);
      const auto                                 dummy = view[0];
    }
  catch (const ExceptionBase &exc)
    {
      deallog << exc.what() << std::endl;
    }

  return 0;
}
