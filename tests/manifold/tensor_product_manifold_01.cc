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


// Test TensorProductManifold

#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tensor_product_manifold.h>

#include "../tests.h"


void
test1()
{
  const int dim = 2, spacedim = 2 + 1;

  FunctionManifold<1, 2, 1> F("x;x^2", "x");
  FunctionManifold<1, 1, 1> G("1.0+2*x", "0.5*(x-1.0)");

  TensorProductManifold<1, 1, 2, 1, 1, 1, 1> manifold(F, G);

  // Chart points.
  Point<2> cp[2];
  cp[1][0] = 1.0;
  cp[1][1] = 0.5;

  // Spacedim points
  std::vector<Point<spacedim>> sp(2);

  // Weights
  std::vector<double> w(2);

  sp[0] = manifold.push_forward(cp[0]);
  sp[1] = manifold.push_forward(cp[1]);

  for (unsigned int d = 0; d < 2; ++d)
    if (cp[d].distance(manifold.pull_back(sp[d])) > 1e-10)
      deallog << "Error! " << cp[d] << "->" << sp[d] << "->"
              << manifold.pull_back(sp[d]) << std::endl;

  unsigned int n_intermediates = 8;

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
}



int
main()
{
  initlog();

  test1();


  return 0;
}
