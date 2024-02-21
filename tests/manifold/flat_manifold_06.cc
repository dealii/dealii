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


// Test direction vector of flat manifold with periodicity

#include <deal.II/grid/manifold.h>

#include "../tests.h"


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
  for (unsigned int d = 0; d < spacedim; ++d)
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
