// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2018 by the deal.II authors
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



// Test a simple parabolic manifold, including gradients and tangent vector

#include "../tests.h"


// all include files you need here
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>

// Helper function
template <int dim, int spacedim>
void
test(unsigned int ref=1)
{
  deallog << "Testing dim " << dim
          << ", spacedim " << spacedim << std::endl;

  // Here the only allowed axis is z. In cylinder the default is x.
  std::string push_forward_expression;
  std::string pull_back_expression;

  switch (spacedim)
    {
    case 2:
      push_forward_expression = "x; x^2";
      pull_back_expression = "x";
      break;
    case 3:
      push_forward_expression = "x; x^2; 0";
      pull_back_expression = "x";
      break;
    default:
      Assert(false, ExcInternalError());
    }

  FunctionManifold<dim,spacedim,1> manifold(push_forward_expression,
                                            pull_back_expression);

  // Two points and two weights
  std::vector<Point<spacedim> > p(2);
  p[1][0] = 1.0;
  p[1][1] = 1.0;
  std::vector<double> w(2);

  unsigned int n_intermediates = 16;

  deallog << "P0: " << p[0]
          << ", P1: " << p[1] << std::endl;

  for (unsigned int i=0; i<n_intermediates+1; ++i)
    {
      w[0] = 1.0-(double)i/((double)n_intermediates);
      w[1] = 1.0 - w[0];

      Point<spacedim> ip = manifold.get_new_point(make_array_view(p),
                                                  make_array_view(w));
      Tensor<1,spacedim> t1 = manifold.get_tangent_vector(ip, p[0]);
      Tensor<1,spacedim> t2 = manifold.get_tangent_vector(ip, p[1]);

      deallog << "P: " << ip
              << ", T(P, P0): " << t1
              << ", T(P, P1): " << t2 << std::endl;

    }
}

int
main ()
{
  initlog();


  test<2,2>();
  test<2,3>();
  test<3,3>();

  return 0;
}
