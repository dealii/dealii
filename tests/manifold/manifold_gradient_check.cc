// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Check push_forward_gradient() of the chart manifolds in manifold_lib.cc
// against central finite differences of push_forward() at generic chart
// points. CylindricalManifold used to miss the factor r in the derivative
// with respect to the angle.

#include <deal.II/grid/manifold_lib.h>

#include <cmath>
#include <string>

#include "../tests.h"


template <int chartdim, int spacedim, typename ManifoldType>
void
check(const std::string     &name,
      const ManifoldType    &manifold,
      const Point<chartdim> &chart_point)
{
  const auto   grad    = manifold.push_forward_gradient(chart_point);
  const double h       = 1e-6;
  double       max_err = 0.;
  double       max_val = 0.;
  bool         finite  = true;

  for (unsigned int d = 0; d < chartdim; ++d)
    {
      Point<chartdim> cp = chart_point, cm = chart_point;
      cp[d] += h;
      cm[d] -= h;
      const Tensor<1, spacedim> fd =
        (manifold.push_forward(cp) - manifold.push_forward(cm)) / (2. * h);
      for (unsigned int i = 0; i < spacedim; ++i)
        {
          finite  = finite && std::isfinite(grad[i][d]);
          max_err = std::max(max_err, std::abs(grad[i][d] - fd[i]));
          max_val = std::max(max_val, std::abs(fd[i]));
        }
    }
  const bool ok = finite && std::isfinite(max_err) && max_val > 0. &&
                  max_err <= 1e-6 * std::max(1., max_val);

  // Only print the verdict: the finite-difference error itself is of the
  // order of 1e-10 and depends on the platform.
  deallog << name << ": " << (ok ? "ok" : "WRONG") << std::endl;
}



int
main()
{
  initlog();

  check<2, 2>("PolarManifold<2>",
              PolarManifold<2>(Point<2>(0.3, -0.2)),
              Point<2>(1.7, 0.6));
  check<3, 3>("PolarManifold<3>",
              PolarManifold<3>(Point<3>(0.1, 0.2, -0.3)),
              Point<3>(1.7, 0.6, 1.1));
  check<3, 3>("CylindricalManifold<3>(axis 2)",
              CylindricalManifold<3>(2),
              Point<3>(2.0, 0.3, 0.5));
  check<3, 3>("CylindricalManifold<3>(direction, point)",
              CylindricalManifold<3>(Tensor<1, 3>({0.3, 0.4, 0.866}),
                                     Point<3>(0.2, -0.1, 0.4)),
              Point<3>(1.5, 2.1, -0.7));
  check<2, 2>("EllipticalManifold<2>",
              EllipticalManifold<2>(Point<2>(0.1, 0.2),
                                    Tensor<1, 2>({0.8, 0.6}),
                                    0.6),
              Point<2>(1.4, 0.9));
  check<3, 3>("TorusManifold<3>",
              TorusManifold<3>(2.0, 0.5),
              Point<3>(0.4, 0.7, 1.3));
}
