// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2018 by the deal.II authors
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

  triangulation.refine_global(2);

  GridOut().write_gnuplot(triangulation, deallog.get_file_stream());
}


int
main()
{
  initlog();

  check<2>();
  check<3>();
}
