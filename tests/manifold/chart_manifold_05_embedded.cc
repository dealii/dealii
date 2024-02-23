// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test direction vector of flat manifold without periodicity, where the
// flat manifold is implemented as a ChartManifold with identity
// pull-back and push-forward
//
// make the chart higher dimensional

#include <deal.II/grid/manifold.h>

#include "../tests.h"


template <int dim, int spacedim>
class MyFlatManifold : public ChartManifold<dim, spacedim, spacedim>
{
public:
  virtual std::unique_ptr<Manifold<dim, spacedim>>
  clone() const override
  {
    return std::unique_ptr<Manifold<dim, spacedim>>(new MyFlatManifold());
  }

  virtual Point<spacedim>
  pull_back(const Point<spacedim> &space_point) const override
  {
    Point<spacedim> p;
    for (unsigned int d = 0; d < spacedim; ++d)
      p[d] = space_point[d];
    return p;
  }


  virtual Point<spacedim>
  push_forward(const Point<spacedim> &chart_point) const override
  {
    Point<spacedim> p;
    for (unsigned int d = 0; d < spacedim; ++d)
      p[d] = chart_point[d];
    return p;
  }

  virtual DerivativeForm<1, spacedim, spacedim>
  push_forward_gradient(const Point<spacedim> &chart_point) const override
  {
    DerivativeForm<1, spacedim, spacedim> x;
    for (unsigned int d = 0; d < spacedim; ++d)
      x[d][d] = 1;
    return x;
  }
};



// Helper function
template <int dim, int spacedim>
void
test()
{
  deallog << "Testing dim=" << dim << ", spacedim=" << spacedim << std::endl;

  MyFlatManifold<dim, spacedim> manifold;

  Point<spacedim> x1, x2;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      x1[d] = 0.1;
      x2[d] = 0.9;
    }

  // get the connecting vector between these two points. because we
  // have no periodicity, this should simply be the vector with
  // components all equal to 0.8
  deallog << manifold.get_tangent_vector(x1, x2) << std::endl;
}

int
main()
{
  initlog();

  test<1, 1>();
  test<1, 2>();
  test<2, 2>();

  return 0;
}
