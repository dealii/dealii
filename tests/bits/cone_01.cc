// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check CylindricalManifold and GridGenerator::truncated_cone

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_c1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
check()
{
  AssertThrow(false, ExcNotImplemented());
}

template <>
void
check<2>()
{
  constexpr int      dim = 2;
  Triangulation<dim> triangulation;
  GridGenerator::truncated_cone(triangulation);

  triangulation.refine_global(2);

  GridOut().write_gnuplot(triangulation, deallog.get_file_stream());
}

template <>
void
check<3>()
{
  constexpr int      dim = 3;
  Triangulation<dim> triangulation;
  GridGenerator::truncated_cone(triangulation);
  static const CylindricalManifold<dim> boundary;
  triangulation.set_manifold(0, boundary);

  triangulation.refine_global(1);

  GridOut().write_gnuplot(triangulation, deallog.get_file_stream());
}


int
main()
{
  initlog();

  check<2>();
  check<3>();
}
