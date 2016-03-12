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

// Check get_tangent_vector for spherical manifold, on simple points.

#include "../tests.h"

#include <deal.II/base/utilities.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/mapping_q_generic.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/manifold_lib.h>

#include <numeric>

template<int dim, int spacedim>
void test()
{
  deallog << "dim=" << dim << ", spacedim=" << spacedim << std::endl;

  Point<spacedim> center;
  static  const SphericalManifold<dim,spacedim> manifold(center);

  // Go from 0,1 to 1,0
  Point<spacedim> p0,p1;
  p0[1] = 1.0;
  p1[0] = 1.0;

  Tensor<1,spacedim> T = manifold.get_tangent_vector(p0, p1);

  deallog << "P0      : " << p0 << std::endl;
  deallog << "P1      : " << p1 << std::endl;
  deallog << "T(P0-P1): " << T << std::endl;
  deallog << "Error   : " << T.norm() - numbers::PI/2 << std::endl;
}


int
main()
{
  std::ofstream logfile ("output");
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  test<2,2>();
  test<2,3>();

  test<3,3>();

  return 0;
}



