// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// check that VectorTools::interpolate works for FE_Nothing
// elements. this is another confirmation that we had previously fixed
// a bug again reported by Krishna Garikipati but that had previously
// been addressed by ensuring that
// FE_Nothing::has_unit_support_points() returns true

#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>

#include <deal.II/dofs/dof_accessor.h>
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
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/numerics/vector_tools.h>

#include <vector>

template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(3);

  hp::FECollection<dim> hp_fe;
  hp_fe.push_back(FESystem<dim>(FE_Q<dim>(2), 1, FE_Nothing<dim>(), 1));
  hp::DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(hp_fe);

  Vector<double> interpolant(dof_handler.n_dofs());
  Vector<float>  error(triangulation.n_active_cells());

  // interpolate the function
  VectorTools::interpolate(
    dof_handler,
    Functions::ZeroFunction<dim>(hp_fe[0].n_components()),
    interpolant);
  deallog << interpolant.l2_norm() << std::endl;
}

int
main()
{
  std::ofstream logfile("output");
  logfile.precision(3);
  deallog << std::setprecision(3);

  deallog.attach(logfile);

  test<1>();
  test<2>();
  test<3>();
}
