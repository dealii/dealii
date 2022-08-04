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

bool
compare_constraints(const AffineConstraints<double> &constraints1,
                    const AffineConstraints<double> &constraints2)
{
  if (constraints1.n_constraints() != constraints2.n_constraints())
    return false;

  for (auto line : constraints1.get_lines())
    {
      const unsigned int line_index = line.index;
      const std::vector<std::pair<types::global_dof_index, double>>
        *constraint_entries_1 = constraints1.get_constraint_entries(line_index);
      const std::vector<std::pair<types::global_dof_index, double>>
        *constraint_entries_2 = constraints2.get_constraint_entries(line_index);
      if (constraint_entries_1->size() != constraint_entries_2->size() ||
          (constraints1.get_inhomogeneity(line_index) !=
           constraints2.get_inhomogeneity(line_index)))
        return false;
      for (unsigned int i = 0; i < constraint_entries_1->size(); ++i)
        {
          if ((*constraint_entries_1)[i].first !=
              (*constraint_entries_2)[i].first)
            return false;
          if ((*constraint_entries_1)[i].second !=
              (*constraint_entries_2)[i].second)
            return false;
        }
    }

  return true;
}

template <int dim>
void
constraints_on_active_cells(
  const Triangulation<dim> &          tria,
  const std::set<types::boundary_id> &no_normal_flux_boundaries,
  AffineConstraints<double> &         constraints)
{
  SphericalManifold<dim> spherical;

  MappingQ<dim> mapping(4);

  FESystem<dim>   fe(FE_Q<dim>(1), dim);
  DoFHandler<dim> dofh(tria);

  dofh.distribute_dofs(fe);

  VectorTools::compute_no_normal_flux_constraints(
    dofh, 0, no_normal_flux_boundaries, constraints, mapping);
}

template <int dim>
void
constraints_on_levels(
  const Triangulation<dim> &                triangulation,
  const std::set<types::boundary_id> &      no_flux_boundary,
  MGLevelObject<AffineConstraints<double>> &level_constraints)
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

      const IndexSet &refinement_edge_indices =
        mg_constrained_dofs.get_refinement_edge_indices(level);

      VectorTools::compute_no_normal_flux_constraints(dof_handler,
                                                      0,
                                                      no_flux_boundary,
                                                      user_level_constraints,
                                                      mapping,
                                                      refinement_edge_indices,
                                                      level);

      // user_level_constraints.print(deallog.get_file_stream());

      // deallog.get_file_stream() << std::flush;
      user_level_constraints.close();
      level_constraints[level].copy_from(user_level_constraints);
    }
}

template <int dim>
void
run(Triangulation<dim> &         triangulation,
    std::set<types::boundary_id> no_flux_boundary)
{
  const unsigned int                     n_levels = 2;
  std::vector<AffineConstraints<double>> ref_constraints(n_levels);

  constraints_on_active_cells<dim>(triangulation,
                                   no_flux_boundary,
                                   ref_constraints[0]);
  for (unsigned int level = 1; level < n_levels; ++level)
    {
      triangulation.refine_global(1);
      constraints_on_active_cells<dim>(triangulation,
                                       no_flux_boundary,
                                       ref_constraints[level]);
    }

  MGLevelObject<AffineConstraints<double>> mg_level_constraints(0,
                                                                n_levels - 1);
  constraints_on_levels<dim>(triangulation,
                             no_flux_boundary,
                             mg_level_constraints);

  deallog << " dim " << dim << std::endl;
  for (unsigned int i = 0; i < n_levels; ++i)
    {
      if (compare_constraints(ref_constraints[i], mg_level_constraints[i]))
        deallog << "Level " << i << " OK" << std::endl;
      else
        deallog << "Level " << i << " failed" << std::endl;
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
    std::set<types::boundary_id> no_flux_boundary{0, 1};
    run(triangulation, no_flux_boundary);
  }
  {
    const unsigned int dim = 3;
    Triangulation<dim> triangulation(
      Triangulation<dim>::limit_level_difference_at_vertices);
    GridGenerator::quarter_hyper_ball(triangulation);
    std::set<types::boundary_id> no_flux_boundary{0, 1};
    run(triangulation, no_flux_boundary);
  }
}
