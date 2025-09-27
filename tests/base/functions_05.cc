// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test ScalarFunctionFromFunctionObject

#include <deal.II/base/function.h>

#include "../tests.h"


template <int dim>
void
check1()
{
  ScalarFunctionFromFunctionObject<dim> object(&Point<dim>::norm);

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
  Point<dim> q;
  for (unsigned int d = 0; d < dim; ++d)
    q[d] = d;

  ScalarFunctionFromFunctionObject<dim> object(
    std::bind(&Point<dim>::distance, q, std::placeholders::_1));

  for (unsigned int i = 0; i < 10; ++i)
    {
      Point<dim> p;
      for (unsigned int d = 0; d < dim; ++d)
        p[d] = i + d;

      AssertThrow(object.value(p) == q.distance(p), ExcInternalError());
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
