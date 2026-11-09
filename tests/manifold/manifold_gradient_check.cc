// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2016 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Checks push_forward_gradient() of the chart manifolds in manifold_lib.cc
// against central finite differences of push_forward() at generic chart points.

#include <deal.II/grid/manifold_lib.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

using namespace dealii;

template <int chartdim, int spacedim, typename M>
bool
check(const std::string &name, const M &m, const Point<chartdim> &chart)
{
  const auto   grad = m.push_forward_gradient(chart);
  const double h    = 1e-6;
  double       max_err = 0, max_val = 0;
  for (unsigned int d = 0; d < chartdim; ++d)
    {
      Point<chartdim> cp = chart, cm = chart;
      cp[d] += h;
      cm[d] -= h;
      const Tensor<1, spacedim> fd = (m.push_forward(cp) - m.push_forward(cm)) / (2 * h);
      for (unsigned int i = 0; i < spacedim; ++i)
        {
          max_err = std::max(max_err, std::abs(grad[i][d] - fd[i]));
          max_val = std::max(max_val, std::abs(fd[i]));
        }
    }
  bool finite = true;
  for (unsigned int i = 0; i < spacedim; ++i)
    for (unsigned int d = 0; d < chartdim; ++d)
      finite = finite && std::isfinite(grad[i][d]);
  const bool ok = finite && std::isfinite(max_err) && max_val > 0 && max_err <= 1e-6 * std::max(1.0, max_val);
  std::cout << std::left << std::setw(46) << name << " max |grad - FD| = " << std::setw(12) << max_err
            << (ok ? "  ok" : "  WRONG") << "\n";
  return ok;
}

int
main()
{
  bool ok = true;
  ok &= check<2, 2>("PolarManifold<2>", PolarManifold<2>(Point<2>(0.3, -0.2)), Point<2>(1.7, 0.6));
  ok &= check<3, 3>("PolarManifold<3> (spherical chart)", PolarManifold<3>(Point<3>(0.1, 0.2, -0.3)),
                    Point<3>(1.7, 0.6, 1.1));
  ok &= check<3, 3>("CylindricalManifold<3>(axis 2)", CylindricalManifold<3>(2), Point<3>(2.0, 0.3, 0.5));
  ok &= check<3, 3>("CylindricalManifold<3>(direction, point)",
                    CylindricalManifold<3>(Tensor<1, 3>({0.3, 0.4, 0.866}), Point<3>(0.2, -0.1, 0.4)),
                    Point<3>(1.5, 2.1, -0.7));
  ok &= check<2, 2>("EllipticalManifold<2>", EllipticalManifold<2>(Point<2>(0.1, 0.2), Tensor<1, 2>({0.8, 0.6}), 0.6),
                    Point<2>(1.4, 0.9)); // chart (c, theta), c > 1
  ok &= check<3, 3>("TorusManifold<3>(R = 2, r = 0.5)", TorusManifold<3>(2.0, 0.5), Point<3>(0.4, 0.7, 1.3));
  std::cout << "\n" << (ok ? "all ok" : "ERRORS FOUND") << "\n";
  return 0;
}
