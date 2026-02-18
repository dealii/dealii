// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2019 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// check that exception macros work on device

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <Kokkos_Core.hpp>

#include "../tests.h"

DEAL_II_HOST_DEVICE
void
test(int dim)
{
  Assert(dim == 2, ExcMessage("a message"));

  if (dim == 2)
    {
      return;
    }
  else
    DEAL_II_NOT_IMPLEMENTED();

  DEAL_II_ASSERT_UNREACHABLE();
}

class Functor
{
public:
  DEAL_II_HOST_DEVICE void
  operator()(const long /*n*/) const
  {
    test(2);
  }
};

int
main()
{
  initlog();

  Kokkos::initialize();

  Functor f;
  Kokkos::parallel_for("single", 1, f);

  Kokkos::finalize();

  deallog << "OK" << std::endl;
  return 0;
}
