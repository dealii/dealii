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

      DEAL_II_AssertThrow(object.value(p) == p.norm(), ExcInternalError());
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

      DEAL_II_AssertThrow(object.value(p) == q.distance(p), ExcInternalError());
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
