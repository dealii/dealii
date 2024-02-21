// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check get_tangent_vector for spherical manifold, on simple points.

#include <deal.II/base/utilities.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <numeric>

#include "../tests.h"

template <int dim, int spacedim>
void
test()
{
  deallog << "dim=" << dim << ", spacedim=" << spacedim << std::endl;

  Point<spacedim>                           center;
  static const PolarManifold<dim, spacedim> manifold(center);

  // Go from 0,1 to 1,0
  Point<spacedim> p0, p1;
  p0[1] = 1.0;
  p1[0] = 1.0;

  Tensor<1, spacedim> T = manifold.get_tangent_vector(p0, p1);

  deallog << "P0      : " << p0 << std::endl;
  deallog << "P1      : " << p1 << std::endl;
  deallog << "T(P0-P1): " << T << std::endl;
  deallog << "Error   : " << T.norm() - numbers::PI / 2 << std::endl;
}


int
main()
{
  initlog();

  test<2, 2>();
  test<2, 3>();
  test<3, 3>();
}
