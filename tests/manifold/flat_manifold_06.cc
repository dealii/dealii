// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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

// Test direction vector of flat manifold with periodicity

#include "../tests.h"
#include <deal.II/grid/manifold.h>

// Helper function
template <int dim, int spacedim>
void
test()
{
  deallog << "Testing dim=" << dim << ", spacedim=" << spacedim << std::endl;

  // make the domain periodic in the first direction with periodicity 1.1
  Tensor<1, spacedim> periodicity;
  periodicity[0] = 1.1;
  FlatManifold<dim, spacedim> manifold(periodicity);

  Point<spacedim> x1, x2;
  for(unsigned int d = 0; d < spacedim; ++d)
    {
      x1[d] = 0.1;
      x2[d] = 0.9;
    }

  // get the connecting vector between these two points. because we
  // have no periodicity, this should simply be the vector with
  // components all equal to 0.8 except for the first, which ought to be -0.3
  deallog << manifold.get_tangent_vector(x1, x2) << std::endl;

  // then also test the opposite direction
  deallog << manifold.get_tangent_vector(x2, x1) << std::endl;
}

int
main()
{
  initlog();

  test<1, 1>();
  test<1, 2>();
  test<1, 3>();
  test<2, 2>();
  test<2, 3>();
  test<3, 3>();

  return 0;
}
