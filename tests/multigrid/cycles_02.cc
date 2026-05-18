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


/*
 * Test multigrid cycles for 1D Poisson problem with unit source term and zero
 * Dirchlet boundary.
 */

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"

using namespace dealii;

int
main()
{
  initlog();

  const unsigned int dim = 1;

  for (unsigned int n_global_refinements = 0; n_global_refinements < 10;
       ++n_global_refinements)
    {
      Triangulation<dim> triangulation(
        Triangulation<dim>::MeshSmoothing::limit_level_difference_at_vertices);
      GridGenerator::hyper_cube(triangulation, 0.0, 1.0);
      triangulation.refine_global(n_global_refinements);

      FE_Q<dim>       fe(1);
      DoFHandler<dim> dof_handler(triangulation);
      dof_handler.distribute_dofs(fe);
      dof_handler.distribute_mg_dofs();

      AffineConstraints<double> constraints;
      DoFTools::make_hanging_node_constraints(dof_handler, constraints);
      VectorTools::interpolate_boundary_values(dof_handler,
                                               0,
                                               Functions::ZeroFunction<dim>(),
                                               constraints);
      VectorTools::interpolate_boundary_values(dof_handler,
                                               1,
                                               Functions::ZeroFunction<dim>(),
                                               constraints);
      constraints.close();

      DynamicSparsityPattern dsp(dof_handler.n_dofs());
      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);

      SparsityPattern      sparsity;
      SparseMatrix<double> system_matrix;
      sparsity.copy_from(dsp);
      system_matrix.reinit(sparsity);

      Vector<double> rhs(dof_handler.n_dofs());

      QGauss<dim>   quadrature(fe.degree + 1);
      FEValues<dim> fe_values(
        fe, quadrature, update_values | update_gradients | update_JxW_values);

      const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
      const unsigned int n_q           = quadrature.size();

      FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
      Vector<double>     cell_rhs(dofs_per_cell);
      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          cell_matrix = 0;
          cell_rhs    = 0;

          fe_values.reinit(cell);

          for (unsigned int q = 0; q < n_q; ++q)
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  cell_matrix(i, j) += fe_values.shape_grad(i, q) *
                                       fe_values.shape_grad(j, q) *
                                       fe_values.JxW(q);

                cell_rhs(i) +=
                  fe_values.shape_value(i, q) * 1.0 * fe_values.JxW(q);
              }

          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(
            cell_matrix, cell_rhs, local_dof_indices, system_matrix, rhs);
        }

      const unsigned int n_levels = triangulation.n_global_levels();

      MGConstrainedDoFs mg_constrained_dofs;
      mg_constrained_dofs.initialize(dof_handler);
      mg_constrained_dofs.make_zero_boundary_constraints(dof_handler, {0, 1});

      MGLevelObject<AffineConstraints<double>> mg_constraints(0, n_levels - 1);
      MGLevelObject<SparsityPattern>           mg_sparsities(0, n_levels - 1);
      MGLevelObject<SparseMatrix<double>>      mg_matrices(0, n_levels - 1);

      for (unsigned int level = 0; level < n_levels; ++level)
        {
          mg_constrained_dofs.merge_constraints(
            mg_constraints[level], level, true, true, true, true);

          DynamicSparsityPattern dsp_level(dof_handler.n_dofs(level),
                                           dof_handler.n_dofs(level));

          MGTools::make_sparsity_pattern(dof_handler, dsp_level, level);

          mg_sparsities[level].copy_from(dsp_level);
          mg_matrices[level].reinit(mg_sparsities[level]);
        }

      for (const auto &cell : dof_handler.cell_iterators())
        if (cell->level_subdomain_id() != numbers::invalid_subdomain_id)
          {
            cell_matrix = 0;
            fe_values.reinit(cell);

            for (unsigned int q = 0; q < n_q; ++q)
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  cell_matrix(i, j) += fe_values.shape_grad(i, q) *
                                       fe_values.shape_grad(j, q) *
                                       fe_values.JxW(q);

            cell->get_mg_dof_indices(local_dof_indices);
            mg_constraints[cell->level()].distribute_local_to_global(
              cell_matrix, local_dof_indices, mg_matrices[cell->level()]);
          }

      mg::Matrix<Vector<double>> mg_matrix(mg_matrices);

      using SmootherPreconditioner = PreconditionJacobi<SparseMatrix<double>>;

      typename SmootherPreconditioner::AdditionalData smoother_data;
      smoother_data.relaxation = 0.7;

      mg::SmootherRelaxation<SmootherPreconditioner, Vector<double>> smoother;
      smoother.initialize(mg_matrices, smoother_data);
      smoother.set_steps(3);
      smoother.set_symmetric(true);


      mg::SmootherRelaxation<SmootherPreconditioner, Vector<double>>
        coarse_smoother;
      coarse_smoother.initialize(mg_matrices, 1.0);
      coarse_smoother.set_steps(1);

      MGCoarseGridApplySmoother<Vector<double>> coarse_grid_solver;
      coarse_grid_solver.initialize(coarse_smoother);

      MGTransferPrebuilt<Vector<double>> transfer(mg_constrained_dofs);
      transfer.build(dof_handler);

      const auto run = [&](const auto cycle_type) {
        Multigrid<Vector<double>> mg(mg_matrix,
                                     coarse_grid_solver,
                                     transfer,
                                     smoother,
                                     smoother,
                                     0,
                                     numbers::invalid_unsigned_int,
                                     cycle_type);

        PreconditionMG<dim, Vector<double>, MGTransferPrebuilt<Vector<double>>>
          preconditioner(dof_handler, mg, transfer);

        ReductionControl                 solver_control(1000, 1e-20, 1e-6);
        SolverRichardson<Vector<double>> solver(solver_control);

        Vector<double> solution(dof_handler.n_dofs());
        solver.solve(system_matrix, solution, rhs, preconditioner);
        constraints.distribute(solution);

        return solver_control.last_step();
      };

      deallog << n_global_refinements << " "
              << run(Multigrid<Vector<double>>::Cycle::v_cycle) << " "
              << run(Multigrid<Vector<double>>::Cycle::w_cycle) << " "
              << run(Multigrid<Vector<double>>::Cycle::f_cycle) << std::endl;
    }
}
