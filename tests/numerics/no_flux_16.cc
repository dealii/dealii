// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2022 by the deal.II authors
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



// add tests for compute_no_normal_flux_constraints_on_level()

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/vector_tools_constraints.h>

#include "../tests.h"

template <int dim>
void
run(const Triangulation<dim> &          triangulation,
    const std::set<types::boundary_id> &no_flux_boundary)
{
  FESystem<dim>   fe(FE_Q<dim>(1), dim);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  dof_handler.distribute_mg_dofs();

  MGConstrainedDoFs mg_constrained_dofs;
  mg_constrained_dofs.initialize(dof_handler);

  const unsigned int n_levels = triangulation.n_global_levels();

  MappingQ<dim> mapping(4);

  for (unsigned int level = 0; level < n_levels; ++level)
    {
      AffineConstraints<double> user_level_constraints;

      VectorTools::compute_no_normal_flux_constraints_on_level(
        dof_handler,
        mg_constrained_dofs,
        level,
        0,
        no_flux_boundary,
        user_level_constraints,
        mapping);

      user_level_constraints.print(deallog.get_file_stream());

      deallog.get_file_stream() << std::flush;
      user_level_constraints.close();

      deallog << "Level " << level << " OK" << std::endl;
    }
}


int
main()
{
  initlog();
  deallog.get_file_stream().precision(7);
  deallog.get_file_stream().setf(std::ios::fixed);

  {
    const unsigned int dim = 2;
    Triangulation<dim> triangulation(
      Triangulation<dim>::limit_level_difference_at_vertices);
    GridGenerator::quarter_hyper_ball(triangulation);
    triangulation.refine_global(1);
    std::set<types::boundary_id> no_flux_boundary{0, 1};
    run<dim>(triangulation, no_flux_boundary);
  }
  {
    const unsigned int dim = 3;
    Triangulation<dim> triangulation(
      Triangulation<dim>::limit_level_difference_at_vertices);
    GridGenerator::quarter_hyper_ball(triangulation);
    triangulation.refine_global(1);
    std::set<types::boundary_id> no_flux_boundary{0, 1};
    run<dim>(triangulation, no_flux_boundary);
  }
}
