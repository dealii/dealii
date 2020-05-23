// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2019 by the deal.II authors
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



// check VectorTools::project_boundary_values_curl_conforming_l2 for a
// combination element that has edge DoFs other than those just from
// the Nedelec base element

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



int
main()
{
  initlog();

  const unsigned int dim = 3;


  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, 0, 1, true);
  MappingQ<dim>             mapping(1);
  AffineConstraints<double> constraints;

  FESystem<dim>   fe_system(FE_Q<dim>(2), 1, FE_Nedelec<dim>(0), 1);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe_system);

  constraints.clear();
  VectorTools::project_boundary_values_curl_conforming_l2(
    dof_handler,
    1,
    Functions::ZeroFunction<dim>(dim + 1),
    0,
    constraints,
    mapping);
  constraints.close();

  constraints.print(deallog.get_file_stream());
}
