// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test the function VectorTools::create_point_source_vector for hp.



#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim>
void
check()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  hp::FECollection<dim> fe_collection;

  for (unsigned int i = 1; i <= tria.n_active_cells(); ++i)
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(i), dim));

  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe_collection);
  Point<dim> orientation;
  Point<dim> p(tria.begin_active()->center());

  for (unsigned int i = 0; i < dim; ++i)
    orientation[i] = i;

  Vector<double> vector(dof.n_dofs());

  VectorTools::create_point_source_vector(dof, p, orientation, vector);

  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    deallog << vector(i) << std::endl;
}



int
main()
{
  initlog();
  deallog << std::setprecision(2);
  deallog << std::fixed;
  deallog.push("1d");
  check<1>();
  deallog.pop();
  deallog.push("2d");
  check<2>();
  deallog.pop();
  deallog.push("3d");
  check<3>();
  deallog.pop();
}
