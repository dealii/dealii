// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check that Tensor works on device

#include <deal.II/base/config.h>

#include "deal.II/base/symmetric_tensor.h"
#include <deal.II/base/exceptions.h>

#include <Kokkos_Core.hpp>

#include "../tests.h"

DEAL_II_HOST_DEVICE
void
test_tensor()
{
  constexpr const int dim = 2;

  Tensor<2, dim> t1;
  Tensor<2, dim> t2;
  t1       = 0;
  t2[0][1] = 1.0;
  t1 += t2;
  t1 *= 2.0;
  Assert(t1.norm() == 2.0, ExcInternalError());
}

class Functor
{
public:
  DEAL_II_HOST_DEVICE void
  operator()(const long n) const
  {
    test_tensor();
  }
};



int
main()
{
  initlog();

  Kokkos::initialize();

  Functor f;
  Kokkos::parallel_for("test", 1, f);

  Kokkos::finalize();

  deallog << "OK" << std::endl;
  return 0;
}
