// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test ScalarFunctionFromFunctionObject with functions with and without
// time dependency.

#include <deal.II/base/function.h>

#include "../tests.h"


template <int dim>
void
check1()
{
  ScalarFunctionFromFunctionObject<dim> object(
    [](const auto &p) { return p.norm(); });

  for (unsigned int i = 0; i < 10; ++i)
    {
      Point<dim> p;
      for (unsigned int d = 0; d < dim; ++d)
        p[d] = i + d;

      AssertThrow(object.value(p) == p.norm(), ExcInternalError());
    }

  deallog << "OK" << std::endl;
}


template <int dim>
void
check2()
{
  ScalarFunctionFromFunctionObject<dim> object(
    [](const auto t, const auto &p) { return p.norm() * t; });

  for (unsigned int i = 0; i < 10; ++i)
    {
      object.set_time(i);

      Point<dim> p;
      for (unsigned int d = 0; d < dim; ++d)
        p[d] = i + d;

      AssertThrow(object.value(p) == (p.norm() * static_cast<double>(i)),
                  ExcInternalError());
    }

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  check1<1>();
  check1<2>();
  check1<3>();

  check2<1>();
  check2<2>();
  check2<3>();
}
