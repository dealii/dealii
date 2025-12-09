// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check
//   DoFTools::map_dofs_to_support_points (const Mapping<dim> &,
//                   const DoFHandler<dim> &,
//                   std::vector<Point<dim> > &, ComponentMask& mask)
// with a system of finite elements containing a finite element which does not
// have support points.

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



template <int dim>
void
test()
{
  Triangulation<dim> triangulation;

  GridGenerator::hyper_cube(triangulation, 0, 1);
  triangulation.refine_global(4 - dim);

  FESystem<dim> fe_velocity(FE_Q<dim>(2), dim);
  FESystem<dim> fe(fe_velocity, 1, FE_DGP<dim>(1), 1);

  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  const MappingQ1<dim> mapping;

  FEValuesExtractors::Vector velocities(0);
  ComponentMask              velocity_mask = fe.component_mask(velocities);

  const FiniteElement<dim> &velocity_fe = fe.base_element(0);

  deallog << "Total components in FESystem: " << fe.n_components() << std::endl;
  deallog << "Base element 0 components: " << velocity_fe.n_components()
          << std::endl;
  deallog << "Velocity mask: " << velocity_mask << std::endl;

  ComponentMask pressure_mask =
    fe.component_mask(FEValuesExtractors::Scalar(dim));
  deallog << "Pressure mask: " << pressure_mask << std::endl;

  const std::map<types::global_dof_index, Point<dim>> &support_points =
    DoFTools::map_dofs_to_support_points(mapping, dof_handler, velocity_mask);

  deallog << "Number of support points: " << support_points.size() << std::endl;
  for (const auto &[dof_idx, point] : support_points)
    deallog << dof_idx << " -> " << point << std::endl;
}


int
main(int argc, char **argv)
{
  initlog();

  deallog.push("1d");
  test<1>();
  deallog.pop();
  deallog.push("2d");
  test<2>();
  deallog.pop();
  deallog.push("3d");
  test<3>();
  deallog.pop();
}
