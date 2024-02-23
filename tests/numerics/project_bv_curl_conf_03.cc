// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------
// The test checks that project_boundary_values_curl_conforming
// works for FESystems with FE_Nedelec of different order.

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


int
main(int argc, char **argv)
{
  initlog();
  const unsigned int dim = 3;

  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, 0, 1, true);
  MappingQ<dim>             mapping(1);
  AffineConstraints<double> constraints;

  // Testing cases with only Nedelec Elements in FE_System
  {
    deallog << "Equal order Nedelec degree 0 \t";
    // Equal-order elements
    FESystem<dim>   fe_system(FE_Nedelec<dim>(0), 1, FE_Nedelec<dim>(0), 1);
    DoFHandler<dim> dof_handler(triangulation);
    dof_handler.distribute_dofs(fe_system);

    constraints.clear();
    VectorTools::project_boundary_values_curl_conforming_l2(
      dof_handler,
      dim,
      Functions::ZeroFunction<dim>(dim + dim),
      0,
      constraints,
      mapping);
    constraints.close();
    deallog << "OK" << std::endl;
  }

  {
    deallog << "Mixed order Nedelec (degree 1,degree 0) \t";
    // Elements of mixed order
    FESystem<dim>   fe_system(FE_Nedelec<dim>(1), 1, FE_Nedelec<dim>(0), 1);
    DoFHandler<dim> dof_handler(triangulation);
    dof_handler.distribute_dofs(fe_system);

    constraints.clear();
    VectorTools::project_boundary_values_curl_conforming_l2(
      dof_handler,
      dim,
      Functions::ZeroFunction<dim>(dim + dim),
      0,
      constraints,
      mapping);
    constraints.close();
    deallog << "OK" << std::endl;
  }

  return 0;
}
