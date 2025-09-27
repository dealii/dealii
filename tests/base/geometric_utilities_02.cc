// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2021 by the deal.II authors
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

#include <deal.II/base/exceptions.h>
#include <deal.II/base/geometric_utilities.h>

#include "../tests.h"


DeclException3(DifferentComponent,
               int,
               double,
               double,
               << arg1 << "-th component is different: " << arg2
               << "!=" << arg3);


template <int dim>
void
test()
{
  for (double r = 0.1; r < 10; r += 0.35)
    for (double theta = 0; theta < 2 * numbers::PI; theta += numbers::PI / 3.)
      // Note: for phi=0 or Pi, \theta can be arbitrary and thereby we can't
      // have one-to-one correspondence
      for (double phi = 0.01; phi <= numbers::PI; phi += numbers::PI / 4.)
        {
          std::array<double, dim> sp;
          sp[0]        = r;
          sp[1]        = theta;
          sp[2]        = phi;
          Point<dim> p = GeometricUtilities::Coordinates::from_spherical(sp);
          const std::array<double, dim> sp2 =
            GeometricUtilities::Coordinates::to_spherical(p);
          for (unsigned int i = 0; i < dim; ++i)
            AssertThrow(std::fabs(sp[i] - sp2[i]) <= std::fabs(sp[i]) * 1e-10,
                        DifferentComponent(i, sp[i], sp2[i]));
        }

  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();

  test<3>();

  return 0;
}
