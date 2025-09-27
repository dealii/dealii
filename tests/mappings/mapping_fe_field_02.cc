// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q1_eulerian.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test()
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria);

  FESystem<dim, spacedim>   fe(FE_Q<dim, spacedim>(1), spacedim);
  DoFHandler<dim, spacedim> dh(tria);

  dh.distribute_dofs(fe);

  deallog << "dim, spacedim: " << dim << ", " << spacedim << std::endl
          << "cells: " << tria.n_active_cells() << ", dofs: " << dh.n_dofs()
          << std::endl;

  // Create a Mapping
  Vector<double> map_vector(dh.n_dofs());
  VectorTools::get_position_vector(dh, map_vector);
  MappingFEField<dim, spacedim> mapping(dh, map_vector);

  tria.refine_global(1);

  dh.distribute_dofs(fe);

  deallog << "After refine:" << std::endl
          << "cells: " << tria.n_active_cells() << ", dofs: " << dh.n_dofs()
          << std::endl;
}

int
main()
{
  initlog();
  test<1, 1>();
  test<1, 2>();
  test<2, 2>();
  test<2, 3>();
  test<3, 3>();
}
