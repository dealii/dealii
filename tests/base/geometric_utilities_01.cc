// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check GeometricUtilities::Coordinates::to_spherical and
// GeometricUtilities::Coordinates::from_spherical.

#include <deal.II/base/geometric_utilities.h>

#include "../tests.h"


using namespace dealii::GeometricUtilities::Coordinates;

template <unsigned int dim, typename T>
void
print(T point1, T point2)
{
  deallog << std::endl << "Point 1: ";
  for (unsigned int i = 0; i < dim; ++i)
    deallog << point1[i] << ' ';

  deallog << std::endl << "Point 2: ";
  for (unsigned int i = 0; i < dim; ++i)
    deallog << point2[i] << ' ';

  deallog << std::endl;
}

template <int dim>
void
test()
{
  Assert(dim > 1, ExcNotImplemented());

  const Point<dim> origin;

  std::array<double, dim> sorigin;
  for (unsigned int d = 0; d < dim; ++d)
    sorigin[d] = 0.;

  Point<dim> one;
  for (unsigned int d = 0; d < dim; ++d)
    one[d] = 1.;

  std::array<double, dim> sone;
  sone[0] = std::sqrt(1. * dim);
  sone[1] = numbers::PI / 4;
  if (dim == 3)
    sone[2] = std::acos(1 / std::sqrt(3.));

  print<dim>(to_spherical(origin), sorigin);
  print<dim>(origin, from_spherical(sorigin));

  print<dim>(to_spherical(one), sone);
  print<dim>(one, from_spherical(sone));
}

void
test3d()
{
  const dealii::Point<3> x(1, 0, 0);
  const dealii::Point<3> y(0, 1, 0);
  const dealii::Point<3> z(0, 0, 1);

  std::array<double, 3> sx;
  sx[0] = 1.;
  sx[1] = 0;
  sx[2] = numbers::PI / 2;
  std::array<double, 3> sy;
  sy[0] = 1;
  sy[1] = numbers::PI / 2;
  sy[2] = numbers::PI / 2;
  std::array<double, 3> sz;
  sz[0] = 1.;
  sz[1] = 0.;
  sz[2] = 0.;

  print<3>(x, from_spherical(sx));
  print<3>(y, from_spherical(sy));
  print<3>(z, from_spherical(sz));

  print<3>(to_spherical(x), sx);
  print<3>(to_spherical(y), sy);
  print<3>(to_spherical(z), sz);

  const Point<3>        dateline(0, -1, 0);
  std::array<double, 3> sdateline;
  sdateline[0] = 1.;
  sdateline[1] = 3 * numbers::PI / 2;
  sdateline[2] = numbers::PI / 2;

  print<3>(dateline, from_spherical(sdateline));
  print<3>(to_spherical(dateline), sdateline);
}


int
main()
{
  initlog();

  test<2>();
  test<3>();
  test3d();

  return 0;
}
