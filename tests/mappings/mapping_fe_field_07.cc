// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// check that MappingFEField is equivalent to MappingQ on a curved
// shell mesh

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test()
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_shell(tria, Point<dim>(), 0.5, 1.);

  tria.refine_global(2);

  FESystem<dim, spacedim>   fe(FE_Q<dim, spacedim>(2), spacedim);
  DoFHandler<dim, spacedim> dh(tria);

  dh.distribute_dofs(fe);

  deallog << "dim, spacedim: " << dim << ", " << spacedim << std::endl
          << "cells: " << tria.n_active_cells() << ", dofs: " << dh.n_dofs()
          << std::endl;

  // Create a Mapping
  Vector<double> map_vector(dh.n_dofs());
  VectorTools::get_position_vector(dh, map_vector);
  MappingFEField<dim, spacedim, Vector<double>> mapping(dh, map_vector);
  MappingQ<dim>                                 mapping_ref(fe.degree);

  QGauss<dim>   quad(1);
  FEValues<dim> fe_values_ref(mapping_ref, fe, quad, update_quadrature_points);
  FEValues<dim> fe_values(mapping, fe, quad, update_quadrature_points);

  for (const auto &cell : tria.active_cell_iterators())
    {
      fe_values_ref.reinit(cell);
      fe_values.reinit(cell);

      if (fe_values_ref.quadrature_point(0).distance(
            fe_values.quadrature_point(0)) > 1e-12)
        deallog << "Mapped point should be "
                << fe_values_ref.quadrature_point(0) << " and is "
                << fe_values.quadrature_point(0) << std::endl;
    }
  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();
  test<2, 2>();
  test<3, 3>();
}
