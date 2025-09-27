// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// like _01, but use a quadratic mapping. since we now map line segments to
// curves, the normal vectors at different quadrature points should no longer
// be parallel

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
test(unsigned int degree)
{
  Triangulation<dim - 1, dim> mesh;
  GridGenerator::hyper_cube(mesh);

  QGauss<dim - 1>        quadrature(dim == 2 ? 3 : 2);
  MappingQ<dim - 1, dim> mapping(degree);
  Point<dim>             p;

  // Try to project a point on the
  // surface
  for (unsigned int i = 0; i < dim; ++i)
    p[i] = .2;

  Point<dim - 1> q =
    mapping.transform_real_to_unit_cell(mesh.begin_active(), p);

  deallog << "Mapping Q(" << degree << "): P: " << p << ", on unit: " << q
          << std::endl;
}



int
main()
{
  initlog();

  test<2>(1);
  test<2>(2);

  test<3>(1);
  test<3>(2);

  return 0;
}
