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


// check that we detect that accessing memory in MemorySpace::Default using an
// ArrayView object is not allowed.

#include <deal.II/base/array_view.h>

#include "../tests.h"

int
main(int argc, char **argv)
{
  deal_II_exceptions::disable_abort_on_exception();

  initlog();

  Kokkos::ScopeGuard                                               guard;
  Kokkos::View<unsigned int *, MemorySpace::Default::kokkos_space> dummy(
    "dummy", 2);

  try
    {
      ArrayView<unsigned int, MemorySpace::Default> view(dummy.data(), 2);
      const auto                                    dummy = view[0];
    }
  catch (const ExceptionBase &exc)
    {
      deallog << exc.what() << std::endl;
    }

  return 0;
}
