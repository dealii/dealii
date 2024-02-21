// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that VectorTools::interpolate works for FE_Nothing
// elements. this is another confirmation that we had previously fixed
// a bug again reported by Krishna Garikipati but that had previously
// been addressed by ensuring that
// FE_Nothing::has_unit_support_points() returns true

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <vector>

#include "../tests.h"


template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(3);

  hp::FECollection<dim> hp_fe;
  hp_fe.push_back(FESystem<dim>(FE_Q<dim>(2), 1, FE_Nothing<dim>(), 1));
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(hp_fe);

  Vector<double> interpolant(dof_handler.n_dofs());
  Vector<float>  error(triangulation.n_active_cells());

  // interpolate the function
  VectorTools::interpolate(dof_handler,
                           Functions::ZeroFunction<dim>(
                             hp_fe[0].n_components()),
                           interpolant);
  deallog << interpolant.l2_norm() << std::endl;
}



int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test<1>();
  test<2>();
  test<3>();
}
