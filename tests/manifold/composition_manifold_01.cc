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


// Test the combination of simple ChartManifolds: parabolic + translation

#include "../tests.h"


// all include files you need here
#include <deal.II/grid/composition_manifold.h>
#include <deal.II/grid/manifold_lib.h>


int
main()
{
  initlog();

  const int dim = 2, spacedim = 2;

  FunctionManifold<dim, spacedim, 1>        F("x;x^2", "x");
  FunctionManifold<dim, spacedim, spacedim> G("x;y+1", "x;y-1");

  CompositionManifold<dim, spacedim, 1, dim, dim, spacedim> manifold(F, G);

  // Chart points.
  Point<1> cp[2];
  cp[1][0] = 1.0;

  // Spacedim points
  std::vector<Point<spacedim>> sp(2);

  // Weights
  std::vector<double> w(2);

  sp[0] = manifold.push_forward(cp[0]);
  sp[1] = manifold.push_forward(cp[1]);

  for (unsigned int d = 0; d < 2; ++d)
    if (cp[d].distance(manifold.pull_back(sp[d])) > 1e-10)
      deallog << "Error!" << std::endl;

  unsigned int n_intermediates = 16;

  deallog << "P0: " << sp[0] << ", P1: " << sp[1] << std::endl;

  for (unsigned int i = 0; i < n_intermediates + 1; ++i)
    {
      w[0] = 1.0 - (double)i / ((double)n_intermediates);
      w[1] = 1.0 - w[0];

      Point<spacedim> ip =
        manifold.get_new_point(make_array_view(sp), make_array_view(w));
      Tensor<1, spacedim> t1 = manifold.get_tangent_vector(ip, sp[0]);
      Tensor<1, spacedim> t2 = manifold.get_tangent_vector(ip, sp[1]);

      deallog << "P: " << ip << ", T(P, P0): " << t1 << ", T(P, P1): " << t2
              << std::endl;
    }


  return 0;
}
