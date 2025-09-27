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


// synthesize the laplacian_list function out of repeated calls to
// laplacian

#include <deal.II/base/function.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"


template <int dim>
class F : public Function<dim>
{
public:
  double
  laplacian(const Point<dim> &p, const unsigned int c) const
  {
    Assert(c == 0, ExcInternalError());
    return p.norm();
  }
};


template <int dim>
void
check()
{
  std::vector<Point<dim>> points;
  for (unsigned int i = 0; i < 10; ++i)
    {
      Point<dim> p;
      for (unsigned int d = 0; d < dim; ++d)
        p[d] = d;
      points.push_back(p);
    }

  F<dim>              f;
  std::vector<double> laplacians(10);
  f.laplacian_list(points, laplacians);

  for (unsigned int i = 0; i < 10; ++i)
    AssertThrow(points[i].norm() == laplacians[i], ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  check<1>();
  check<2>();
  check<3>();
}
