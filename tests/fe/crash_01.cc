// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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


// Test project_boundary_values_curl_conforming_l2

#include <deal.II/base/function.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/vector_tools.h>

#include <string>
#include <vector>

#include "../tests.h"


template <int dim>
void
test(FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl;
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, -1.0, 1.0);
  tria.refine_global(1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  AffineConstraints<double>    constraints;
  Functions::ZeroFunction<dim> boundary_values(fe.n_components());
  VectorTools::project_boundary_values_curl_conforming_l2(
    dof_handler,
    0,
    boundary_values,
    0,
    constraints,
    StaticMappingQ1<dim>::mapping);
}



int
main()
{
  initlog();
  deallog << std::setprecision(7);

  FE_Nedelec<3> fe1(0); // works
  test<3>(fe1);
  FESystem<3> fe2(FE_Nedelec<3>(0), 2); // crash
  test<3>(fe2);

  return 0;
}
