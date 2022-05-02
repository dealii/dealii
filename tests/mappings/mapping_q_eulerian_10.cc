// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// Check that it is possible to perform two consecutive distribute dofs
// when using MappingQEulerian

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q1_eulerian.h>
#include <deal.II/fe/mapping_q_eulerian.h>

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
  Vector<double>                                  map_vector(dh.n_dofs());
  MappingQEulerian<dim, Vector<double>, spacedim> mapping(1, dh, map_vector);

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
