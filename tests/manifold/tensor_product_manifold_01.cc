// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// Test TensorProductManifold

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>

#include <deal.II/grid/tensor_product_manifold.h>
#include <deal.II/grid/manifold_lib.h>


void test1()
{
  const int dim=2, spacedim=2+1;

  FunctionManifold<1,2,1>    F("x;x^2", "x");
  FunctionManifold<1,1,1>    G("1.0+2*x", "0.5*(x-1.0)");

  TensorProductManifold<1, 1,2,1, 1,1,1> manifold(F, G);

  // Chart points.
  Point<2> cp[2];
  cp[1][0] = 1.0;
  cp[1][1] = 0.5;

  // Spacedim points
  std::vector<Point<spacedim> > sp(2);

  // Weights
  std::vector<double> w(2);

  sp[0] = manifold.push_forward(cp[0]);
  sp[1] = manifold.push_forward(cp[1]);

  for (unsigned int d=0; d<2; ++d)
    if (cp[d].distance(manifold.pull_back(sp[d])) > 1e-10)
      deallog << "Error! "
              << cp[d] << "->" << sp[d] << "->" << manifold.pull_back(sp[d])
              << std::endl;

  unsigned int n_intermediates = 8;

  deallog << "P0: " << sp[0]
          << ", P1: " << sp[1] << std::endl;

  for (unsigned int i=0; i<n_intermediates+1; ++i)
    {
      w[0] = 1.0-(double)i/((double)n_intermediates);
      w[1] = 1.0 - w[0];

      Point<spacedim> ip = manifold.get_new_point(Quadrature<spacedim>(sp, w));
      Tensor<1,spacedim> t1 = manifold.get_tangent_vector(ip, sp[0]);
      Tensor<1,spacedim> t2 = manifold.get_tangent_vector(ip, sp[1]);

      deallog << "P: " << ip
              << ", T(P, P0): " << t1
              << ", T(P, P1): " << t2 << std::endl;

    }
}



int main ()
{
  initlog();

  test1();


  return 0;
}

