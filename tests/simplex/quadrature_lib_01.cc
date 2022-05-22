// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2022 by the deal.II authors
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


// Test QGaussSimplex: output its quadrature points and weights.


#include <deal.II/base/quadrature_lib.h>

#include "../tests.h"

using namespace dealii;

template <int dim>
void
test(const unsigned int n_points)
{
  QGaussSimplex<dim> quad(n_points);

  for (unsigned int q = 0; q < quad.size(); ++q)
    {
      deallog << quad.point(q) << ' ';
      deallog << quad.weight(q) << ' ';
      deallog << std::endl;
    }
}

int
main()
{
  initlog();

  {
    deallog.push("1d-1");
    test<1>(1 /*n_points*/);
    deallog.pop();
  }
  {
    deallog.push("1d-2");
    test<1>(2);
    deallog.pop();
  }
  {
    deallog.push("1d-3");
    test<1>(3);
    deallog.pop();
  }
  {
    deallog.push("2d-1");
    test<2>(1);
    deallog.pop();
  }
  {
    deallog.push("2d-3");
    test<2>(2);
    deallog.pop();
  }
  {
    deallog.push("2d-7");
    test<2>(3);
    deallog.pop();
  }

  {
    deallog.push("3d-1");
    test<3>(1);
    deallog.pop();
  }
  {
    deallog.push("3d-4");
    test<3>(2);
    deallog.pop();
  }
  {
    deallog.push("3d-10");
    test<3>(3);
    deallog.pop();
  }
}
