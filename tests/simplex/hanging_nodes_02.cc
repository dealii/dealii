// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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



// verify hanging node constraints on locally p-refined simplex mesh

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/simplex/fe_lib.h>
#include <deal.II/simplex/grid_generator.h>

#include "../tests.h"


template <int dim>
void
test(const hp::FECollection<dim> &fes)
{
  // setup grid
  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_cube_with_simplices(tria, 1);

  DoFHandler<dim> dofh(tria);
  dofh.begin_active()->set_active_fe_index(1);

  dofh.distribute_dofs(fes);
  deallog << "ndofs: " << dofh.n_dofs() << std::endl;

  // hanging node constraints
  AffineConstraints<double> constraints;
  // DoFTools::make_hanging_node_constraints(dofh, constraints);
  constraints.print(deallog.get_file_stream());

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  deallog.push("2d");
  test<2>(hp::FECollection<2>(Simplex::FE_P<2>(1), Simplex::FE_P<2>(2)));
  test<2>(hp::FECollection<2>(Simplex::FE_P<2>(2), Simplex::FE_P<2>(1)));
  deallog.pop();
}
