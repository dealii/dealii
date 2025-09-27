// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test compute_no_normal_flux_constraints_on_lefel().

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
  const double tol = 1.e-15;
  if (constraints1.n_constraints() != constraints2.n_constraints())
    {
      deallog << "n_constraints_1 = " << constraints1.n_constraints()
              << "n_constraints_2 = " << constraints2.n_constraints()
              << std::endl;
      return false;
    }

  auto line1    = (constraints1.get_lines()).begin();
  auto line2    = (constraints2.get_lines()).begin();
  auto line_end = (constraints1.get_lines()).end();

  for (; line1 != line_end; ++line1, ++line2)
    {
      const unsigned int line_index  = line1->index;
      const unsigned int line_index2 = line2->index;
      if (line_index != line_index2)
        {
          deallog << " line_index_1 = " << line_index
                  << " line_index_2 = " << line_index2 << std::endl;
          return false;
        }

      const std::vector<std::pair<types::global_dof_index, double>>
        *constraint_entries_1 = constraints1.get_constraint_entries(line_index);
      const std::vector<std::pair<types::global_dof_index, double>>
        *constraint_entries_2 = constraints2.get_constraint_entries(line_index);
      if (constraint_entries_1 == nullptr && constraint_entries_2 == nullptr)
        {
          return true;
        }
      else if (constraint_entries_1 != nullptr &&
               constraint_entries_2 != nullptr)
        {
          if (constraint_entries_1->size() != constraint_entries_2->size() ||
              std::abs(constraints1.get_inhomogeneity(line_index) -
                       constraints2.get_inhomogeneity(line_index)) > tol)
            {
              deallog << " constraints1.get_inhomogeneity(line_index) "
                      << constraints1.get_inhomogeneity(line_index)
                      << " constraints2.get_inhomogeneity(line_index) "
                      << constraints2.get_inhomogeneity(line_index) << std::endl
                      << " constraint_entries_1->size() "
                      << constraint_entries_1->size()
                      << " constraint_entries_2->size() "
                      << constraint_entries_2->size() << std::endl;
              return false;
            }
          for (unsigned int i = 0; i < constraint_entries_1->size(); ++i)
            {
              if ((*constraint_entries_1)[i].first !=
                    (*constraint_entries_2)[i].first ||
                  std::abs((*constraint_entries_1)[i].second -
                           (*constraint_entries_2)[i].second) > tol)
                {
                  deallog << " (*constraint_entries_1)[i].first "
                          << (*constraint_entries_1)[i].first
                          << " (*constraint_entries_1)[i].second "
                          << (*constraint_entries_1)[i].second << " diff: "
                          << (*constraint_entries_1)[i].first -
                               (*constraint_entries_2)[i].first
                          << std::endl
                          << " (*constraint_entries_2)[i].first "
                          << (*constraint_entries_2)[i].first
                          << " (*constraint_entries_2)[i].second "
                          << (*constraint_entries_2)[i].second << " diff: "
                          << (*constraint_entries_1)[i].first -
                               (*constraint_entries_2)[i].first
                          << std::endl;
                  return false;
                }
            }
        }
      else
        {
          return false;
        }
    }

  return true;
}

template <int dim>
void
get_constraints_on_active_cells(
  const Triangulation<dim>           &tria,
  const std::set<types::boundary_id> &no_normal_flux_boundaries,
  AffineConstraints<double>          &constraints)
{
  MappingQ<dim> mapping(4);

  FESystem<dim>   fe(FE_Q<dim>(1), dim);
  DoFHandler<dim> dofh(tria);

  dofh.distribute_dofs(fe);

  VectorTools::compute_no_normal_flux_constraints(
    dofh, 0, no_normal_flux_boundaries, constraints, mapping);
  constraints.close();
  // deallog<<" ***active cell constraints:
  // "<<constraints.n_constraints()<<std::endl;
  // constraints.print(deallog.get_file_stream());
}

template <int dim>
void
get_constraints_on_levels(
  const Triangulation<dim>                 &triangulation,
  const std::set<types::boundary_id>       &no_flux_boundary,
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

      VectorTools::compute_no_normal_flux_constraints_on_level(
        dof_handler,
        0,
        no_flux_boundary,
        user_level_constraints,
        mapping,
        refinement_edge_indices,
        level);
      user_level_constraints.close();
      level_constraints[level].copy_from(user_level_constraints);
      // deallog<<" ***level cell constraints:
      // "<<user_level_constraints.n_constraints()<<std::endl;
      // user_level_constraints.print(deallog.get_file_stream());
    }
}

template <int dim>
void
run(Triangulation<dim>          &triangulation,
    std::set<types::boundary_id> no_flux_boundary)
{
  const unsigned int                     n_levels = 2;
  std::vector<AffineConstraints<double>> ref_constraints(n_levels);

  get_constraints_on_active_cells<dim>(triangulation,
                                       no_flux_boundary,
                                       ref_constraints[0]);
  for (unsigned int level = 1; level < n_levels; ++level)
    {
      triangulation.refine_global(1);
      get_constraints_on_active_cells<dim>(triangulation,
                                           no_flux_boundary,
                                           ref_constraints[level]);
    }

  MGLevelObject<AffineConstraints<double>> mg_level_constraints(0, n_levels);
  get_constraints_on_levels<dim>(triangulation,
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
    deallog << " quarter_hyper_ball: " << std::endl;
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
  {
    deallog << " hyper_shell: " << std::endl;
    const unsigned int dim = 2;
    Triangulation<dim> triangulation(
      Triangulation<dim>::limit_level_difference_at_vertices);
    GridGenerator::hyper_shell(triangulation, Point<dim>(), 0.5, 1.);
    std::set<types::boundary_id> no_flux_boundary{0};
    run(triangulation, no_flux_boundary);
  }
  {
    const unsigned int dim = 3;
    Triangulation<dim> triangulation(
      Triangulation<dim>::limit_level_difference_at_vertices);
    GridGenerator::hyper_shell(triangulation, Point<dim>(), 0.5, 1.);
    std::set<types::boundary_id> no_flux_boundary{0};
    run(triangulation, no_flux_boundary);
  }
  {
    deallog << " half_hyper_shell: " << std::endl;
    const unsigned int dim = 2;
    Triangulation<dim> triangulation(
      Triangulation<dim>::limit_level_difference_at_vertices);
    GridGenerator::half_hyper_shell(triangulation, Point<dim>(), 0.5, 1.);
    std::set<types::boundary_id> no_flux_boundary{0};
    run(triangulation, no_flux_boundary);
  }
  {
    const unsigned int dim = 3;
    Triangulation<dim> triangulation(
      Triangulation<dim>::limit_level_difference_at_vertices);
    GridGenerator::half_hyper_shell(triangulation, Point<dim>(), 0.5, 1.);
    std::set<types::boundary_id> no_flux_boundary{0};
    run(triangulation, no_flux_boundary);
  }
  {
    deallog << " quarter_hyper_shell: " << std::endl;
    const unsigned int dim = 2;
    Triangulation<dim> triangulation(
      Triangulation<dim>::limit_level_difference_at_vertices);
    GridGenerator::quarter_hyper_shell(triangulation, Point<dim>(), 0.5, 1.);
    std::set<types::boundary_id> no_flux_boundary{0};
    run(triangulation, no_flux_boundary);
  }
  {
    const unsigned int dim = 3;
    Triangulation<dim> triangulation(
      Triangulation<dim>::limit_level_difference_at_vertices);
    GridGenerator::quarter_hyper_shell(triangulation, Point<dim>(), 0.5, 1.);
    std::set<types::boundary_id> no_flux_boundary{0};
    run(triangulation, no_flux_boundary);
  }
  {
    deallog << " hyper_cube:" << std::endl;
    const unsigned int           dim = 2;
    std::set<types::boundary_id> no_flux_boundary;
    for (unsigned int i = 0; i < GeometryInfo<dim>::faces_per_cell; ++i)
      {
        no_flux_boundary.insert(i);
        Triangulation<dim> triangulation(
          Triangulation<dim>::limit_level_difference_at_vertices);
        GridGenerator::hyper_cube(triangulation, 0., 1., true);
        run(triangulation, no_flux_boundary);
      }
  }
  {
    const unsigned int           dim = 3;
    std::set<types::boundary_id> no_flux_boundary;
    for (unsigned int i = 0; i < GeometryInfo<dim>::faces_per_cell; ++i)
      {
        no_flux_boundary.insert(i);
        Triangulation<dim> triangulation(
          Triangulation<dim>::limit_level_difference_at_vertices);
        GridGenerator::hyper_cube(triangulation, 0., 1., true);
        run(triangulation, no_flux_boundary);
      }
  }
  {
    // test adaptive refinement
    // the output file is generated by the old version of
    // compute_no_normal_flux_constraints() and copied here.
    deallog << " adaptive hyper_shell:" << std::endl;

    const unsigned int dim = 2;
    FESystem<dim>      fe(FE_Q<dim>(2), dim);

    const MappingQ<dim> mapping(4);
    Triangulation<dim>  triangulation(
      Triangulation<dim>::limit_level_difference_at_vertices);
    GridGenerator::hyper_shell(triangulation, Point<dim>(), 0.5, 1.);
    std::set<types::boundary_id> no_flux_boundary{0};
    triangulation.begin_active()->set_refine_flag();
    triangulation.execute_coarsening_and_refinement();

    DoFHandler<dim> dofh(triangulation);
    dofh.distribute_dofs(fe);

    AffineConstraints<double> constraints;

    VectorTools::compute_no_normal_flux_constraints(
      dofh, 0, no_flux_boundary, constraints, mapping);

    constraints.print(deallog.get_file_stream());
  }
}
