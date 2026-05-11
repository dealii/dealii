// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// Test based on a performance benchmark of DynamicSparsityPattern. Prints some
// essential statistics, but can be turned back into a performance benchmark
// (with details on memory usage and wall clock time).

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/numerics/vector_tools_boundary.h>
#include <deal.II/numerics/vector_tools_constraints.h>

#include <chrono>

#include "../tests.h"

using namespace dealii;

static bool run_as_benchmark = false;

template <typename SparsityType>
void
print_quartiles(const SparsityType &dsp)
{
  if (run_as_benchmark)
    deallog << "mem  = " << dsp.memory_consumption() / 1024 / 1024 << " MB"
            << std::endl;

  std::vector<std::size_t> row_lengths(dsp.n_rows());
  for (std::size_t row_no = 0; row_no < dsp.n_rows(); ++row_no)
    row_lengths[row_no] = dsp.row_length(row_no);

  std::sort(row_lengths.begin(), row_lengths.end());
  deallog << "q0   = " << row_lengths.front() << std::endl;
  deallog << "q1   = " << row_lengths[row_lengths.size() / 4] << std::endl;
  deallog << "q2   = " << row_lengths[row_lengths.size() / 2] << std::endl;
  deallog << "q3   = " << row_lengths[3 * row_lengths.size() / 4] << std::endl;
  deallog << "high = " << row_lengths.back() << std::endl;
}

int
main()
{
  initlog();

  // benchmark based on step-40:
  {
    Triangulation<3> tria;
    GridGenerator::hyper_cube(tria);
    tria.refine_global(3);

    for (unsigned int i = 0; i < 3; ++i)
      {
        const auto diameter = tria.begin_active()->diameter();
        for (auto &cell : tria.active_cell_iterators())
          {
            const auto center = cell->barycenter();

            if (std::abs(0.5 + 0.25 * std::sin(4.0 * numbers::PI * center[0]) -
                         center[1]) < diameter)
              cell->set_refine_flag();
          }

        tria.execute_coarsening_and_refinement();
      }

    DoFHandler<3> dof_handler(tria);
    FE_Q<3>       fe(2);
    dof_handler.distribute_dofs(fe);

    const auto locally_owned_dofs = dof_handler.locally_owned_dofs();
    const auto locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);
    AffineConstraints<double> constraints(locally_owned_dofs,
                                          locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<3>(),
                                             constraints);
    constraints.close();

    std::chrono::milliseconds best_time(1000000);
    const std::size_t         upper_bound = run_as_benchmark ? 32 : 1;
    for (unsigned int i = 0; i < upper_bound; ++i)
      {
        DynamicSparsityPattern dsp(locally_relevant_dofs);
        const auto             t0 = std::chrono::high_resolution_clock::now();
        DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
        const auto t1 = std::chrono::high_resolution_clock::now();

        const auto current_time =
          std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);

        if (current_time < best_time)
          best_time = current_time;

        if (i == upper_bound - 1)
          print_quartiles(dsp);
      }
    deallog << "number of DoFs    = " << dof_handler.n_dofs() << std::endl;
    deallog << "dofs per cell     = " << fe.n_dofs_per_cell() << std::endl;
    if (run_as_benchmark)
      deallog << "best time (of " << upper_bound << ") = " << best_time.count()
              << "ms" << std::endl;
  }

  // benchmark based on step-69:
  {
    Triangulation<3> tria;
    GridGenerator::hyper_cube(tria);
    tria.refine_global(6);

    DoFHandler<3> dof_handler(tria);
    FE_Q<3>       fe(1);
    dof_handler.distribute_dofs(fe);
    const auto locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);

    std::chrono::milliseconds best_time(1000000);
    const std::size_t         upper_bound = run_as_benchmark ? 32 : 1;
    for (unsigned int i = 0; i < upper_bound; ++i)
      {
        DynamicSparsityPattern dsp(locally_relevant_dofs);

        const auto t0 = std::chrono::high_resolution_clock::now();
        std::vector<types::global_dof_index> dof_indices(fe.n_dofs_per_cell());
        for (const auto &cell : dof_handler.active_cell_iterators())
          {
            cell->get_dof_indices(dof_indices);

            for (const auto dof : dof_indices)
              dsp.add_entries(dof, dof_indices.begin(), dof_indices.end());
          }
        const auto t1 = std::chrono::high_resolution_clock::now();

        const auto current_time =
          std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);

        if (current_time < best_time)
          best_time = current_time;

        if (i == upper_bound - 1)
          print_quartiles(dsp);
      }
    deallog << "number of DoFs    = " << dof_handler.n_dofs() << std::endl;
    deallog << "dofs per cell     = " << fe.n_dofs_per_cell() << std::endl;
    if (run_as_benchmark)
      deallog << "best time (of " << upper_bound << ") = " << best_time.count()
              << "ms" << std::endl;
  }

  // benchmark based on step-31 and step-32 (velocity only):
  {
    Triangulation<3> tria;
    GridGenerator::hyper_shell(tria, Point<3>(), 1.0, 2.0, 96, true);
    tria.refine_global(1);

    for (unsigned int i = 0; i < 2; ++i)
      {
        const auto diameter = tria.begin_active()->diameter();
        for (auto &cell : tria.active_cell_iterators())
          {
            const auto center = cell->barycenter();

            for (unsigned int d = 0; d < 3; ++d)
              if (std::abs(center[d]) < diameter)
                cell->set_refine_flag();

            for (const auto &face : cell->face_iterators())
              if (face->boundary_id() == 1)
                cell->set_refine_flag();
          }

        tria.execute_coarsening_and_refinement();
      }

    DoFHandler<3> dof_handler(tria);
    FESystem<3>   fe(FE_Q<3>(2), 3);
    dof_handler.distribute_dofs(fe);

    DoFRenumbering::component_wise(dof_handler);

    const auto locally_owned_dofs = dof_handler.locally_owned_dofs();
    const auto locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);
    AffineConstraints<double> constraints(locally_owned_dofs,
                                          locally_relevant_dofs);

    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    const FEValuesExtractors::Vector velocity_components(0);
    VectorTools::interpolate_boundary_values(
      dof_handler,
      1,
      Functions::ZeroFunction<3>(fe.n_components()),
      constraints,
      fe.component_mask(velocity_components));
    const std::set<types::boundary_id> no_normal_flux_boundaries = {0};
    VectorTools::compute_no_normal_flux_constraints(dof_handler,
                                                    0,
                                                    no_normal_flux_boundaries,
                                                    constraints);
    constraints.close();

    std::chrono::milliseconds best_time(1000000);
    std::size_t               upper_bound = run_as_benchmark ? 8 : 1;
    for (unsigned int i = 0; i < upper_bound; ++i)
      {
        DynamicSparsityPattern dsp(locally_relevant_dofs);

        const auto t0 = std::chrono::high_resolution_clock::now();
        DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
        const auto t1 = std::chrono::high_resolution_clock::now();

        const auto current_time =
          std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);

        if (current_time < best_time)
          best_time = current_time;

        if (i == upper_bound - 1)
          print_quartiles(dsp);
      }
    deallog << "number of DoFs    = " << dof_handler.n_dofs() << std::endl;
    deallog << "dofs per cell     = " << fe.n_dofs_per_cell() << std::endl;
    if (run_as_benchmark)
      deallog << "best time (of " << upper_bound << ")  = " << best_time.count()
              << "ms" << std::endl;
  }

  // benchmark based on step-31 and step-32 (Stokes):
  {
    Triangulation<3> tria;
    GridGenerator::hyper_shell(tria, Point<3>(), 1.0, 2.0, 96, true);
    tria.refine_global(1);

    for (unsigned int i = 0; i < 2; ++i)
      {
        const auto diameter = tria.begin_active()->diameter();
        for (auto &cell : tria.active_cell_iterators())
          {
            const auto center = cell->barycenter();

            for (unsigned int d = 0; d < 3; ++d)
              if (std::abs(center[d]) < diameter)
                cell->set_refine_flag();

            for (const auto &face : cell->face_iterators())
              if (face->boundary_id() == 1)
                cell->set_refine_flag();
          }

        tria.execute_coarsening_and_refinement();
      }

    DoFHandler<3> dof_handler(tria);
    FESystem<3>   fe(FE_Q<3>(2) ^ 3, FE_Q<3>(1));
    dof_handler.distribute_dofs(fe);

    const std::vector<unsigned int> sub_blocks{0, 0, 0, 1};
    DoFRenumbering::component_wise(dof_handler, sub_blocks);

    const auto locally_owned_dofs = dof_handler.locally_owned_dofs();
    const auto locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);
    AffineConstraints<double> constraints(locally_owned_dofs,
                                          locally_relevant_dofs);

    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    const FEValuesExtractors::Vector velocity_components(0);
    VectorTools::interpolate_boundary_values(
      dof_handler,
      1,
      Functions::ZeroFunction<3>(fe.n_components()),
      constraints,
      fe.component_mask(velocity_components));
    const std::set<types::boundary_id> no_normal_flux_boundaries = {0};
    VectorTools::compute_no_normal_flux_constraints(dof_handler,
                                                    0,
                                                    no_normal_flux_boundaries,
                                                    constraints);
    constraints.close();

    const std::vector<types::global_dof_index> dofs_per_block =
      DoFTools::count_dofs_per_fe_block(dof_handler, sub_blocks);

    const types::global_dof_index n_u = dofs_per_block[0],
                                  n_p = dofs_per_block[1];
    std::vector<IndexSet> partitioning(2);
    partitioning[0] = complete_index_set(n_u);
    partitioning[1] = complete_index_set(n_p);

    Table<2, DoFTools::Coupling> coupling(fe.n_components(), fe.n_components());
    for (unsigned int c = 0; c < fe.n_components(); ++c)
      for (unsigned int d = 0; d < fe.n_components(); ++d)
        if (c == 3 && d == 3)
          coupling[c][d] = DoFTools::none;
        else
          coupling[c][d] = DoFTools::always;

    std::chrono::milliseconds best_time(1000000);
    const std::size_t         upper_bound = run_as_benchmark ? 8 : 1;
    for (unsigned int i = 0; i < 8; ++i)
      {
        BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);

        const auto t0 = std::chrono::high_resolution_clock::now();
        DoFTools::make_sparsity_pattern(
          dof_handler, coupling, dsp, constraints, false);
        const auto t1 = std::chrono::high_resolution_clock::now();

        const auto current_time =
          std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);

        if (current_time < best_time)
          best_time = current_time;

        if (i == upper_bound - 1)
          print_quartiles(dsp);
      }
    deallog << "number of DoFs    = " << dof_handler.n_dofs() << std::endl;
    deallog << "dofs per cell     = " << fe.n_dofs_per_cell() << std::endl;
    if (run_as_benchmark)
      deallog << "best time (of " << upper_bound << ")  = " << best_time.count()
              << "ms" << std::endl;
  }
}
