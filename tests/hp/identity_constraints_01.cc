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



// Create a 2x1x1 grid with FE_Q elements with degrees 2 & 4 assigned.
// Verify that we do not unify dofs on lines in 3D.
// On each of the four lines on the interface, the central dofs are identical
// and will be treated with constraints.
//
// +----+----+
// | Q2 | Q4 |
// +----+----+

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/lac/affine_constraints.h>

#include "../tests.h"

#include "../test_grids.h"


template <int dim>
void
test()
{
  Triangulation<dim> tria;
  TestGrids::hyper_line(tria, 2);

  hp::FECollection<dim> fe;
  fe.push_back(FE_Q<dim>(4));
  fe.push_back(FE_Q<dim>(2));

  DoFHandler<dim> dh(tria);
  dh.begin_active()->set_active_fe_index(1);
  dh.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dh, constraints);
  constraints.close();

  deallog << "Total constraints:          " << constraints.n_constraints()
          << std::endl
          << "  Inhomogenous constraints: " << constraints.n_inhomogeneities()
          << std::endl
          << "  Identity constraints:     " << constraints.n_identities()
          << std::endl;
}


int
main()
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
