// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


/**
 * Check that everything works also for locally-refined grids.
 *
 */

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include "multigrid_util.h"


template <int dim>
void
create_quadrant(Triangulation<dim> &tria, const unsigned int n_refinements)
{
  // Taken from @cite munch2022gc, according to the description in A FLEXIBLE,
  // PARALLEL, ADAPTIVE GEOMETRIC MULTIGRID METHOD FOR FEM (Clevenger, Heister,
  // Kanschat, Kronbichler): https://arxiv.org/pdf/1904.03317.pdf

  GridGenerator::hyper_cube(tria, -1.0, +1.0);

  if (n_refinements == 0)
    return;

  tria.refine_global(1);

  for (unsigned int i = 1; i < n_refinements; i++)
    {
      for (auto cell : tria.active_cell_iterators())
        if (cell->is_locally_owned())
          {
            bool flag = true;
            for (int d = 0; d < dim; d++)
              if (cell->center()[d] > 0.0)
                flag = false;
            if (flag)
              cell->set_refine_flag();
          }
      tria.execute_coarsening_and_refinement();
    }

  AssertDimension(tria.n_global_levels() - 1, n_refinements);
}


template <int dim, typename Number = double>
void
test(const unsigned int n_refinements,
     const unsigned int fe_degree_fine,
     const bool         do_simplex_mesh)
{
  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  const unsigned int min_level = 0;
  const unsigned int max_level = n_refinements;

  MGLevelObject<std::shared_ptr<parallel::distributed::Triangulation<dim>>>
                                           triangulations(min_level, max_level);
  MGLevelObject<DoFHandler<dim>>           dof_handlers(min_level, max_level);
  MGLevelObject<AffineConstraints<Number>> constraints(min_level, max_level);
  MGLevelObject<std::unique_ptr<Mapping<dim>>> mappings(min_level, max_level);
  MGLevelObject<std::shared_ptr<MGTwoLevelTransferNonNested<dim, VectorType>>>
                                          transfers(min_level, max_level);
  MGLevelObject<Operator<dim, 1, Number>> operators(min_level, max_level);


  // set up levels
  for (auto l = min_level; l <= max_level; ++l)
    {
      triangulations[l] =
        std::make_shared<parallel::distributed::Triangulation<dim>>(
          MPI_COMM_WORLD);
      auto &dof_handler = dof_handlers[l];
      auto &constraint  = constraints[l];
      auto &op          = operators[l];

      std::unique_ptr<FiniteElement<dim>> fe;
      std::unique_ptr<Quadrature<dim>>    quad;
      std::unique_ptr<Mapping<dim>>       mapping;

      if (do_simplex_mesh)
        {
          fe      = std::make_unique<FE_SimplexP<dim>>(fe_degree_fine);
          quad    = std::make_unique<QGaussSimplex<dim>>(fe_degree_fine + 1);
          mapping = std::make_unique<MappingFE<dim>>(FE_SimplexP<dim>(1));
        }
      else
        {
          fe      = std::make_unique<FE_Q<dim>>(fe_degree_fine);
          quad    = std::make_unique<QGauss<dim>>(fe_degree_fine + 1);
          mapping = std::make_unique<MappingQ<dim>>(1);
        }

      mappings[l] = mapping->clone();

      // set up triangulation

      create_quadrant(*triangulations[l], 2 * l + 1);

      // set up dofhandler
      dof_handler.reinit(*triangulations[l]);
      dof_handler.distribute_dofs(*fe);

      // set up constraints
      const IndexSet locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(dof_handler);
      constraint.reinit(dof_handler.locally_owned_dofs(),
                        locally_relevant_dofs);
      VectorTools::interpolate_boundary_values(
        *mapping, dof_handler, 0, Functions::ZeroFunction<dim>(), constraint);
      DoFTools::make_hanging_node_constraints(dof_handler, constraint);
      constraint.close();

      // set up operator
      op.reinit(*mapping, dof_handler, *quad, constraint);
    }

  // set up transfer operator
  for (unsigned int l = min_level; l < max_level; ++l)
    {
      transfers[l + 1] =
        std::make_shared<MGTwoLevelTransferNonNested<dim, VectorType>>();
      transfers[l + 1]->reinit(dof_handlers[l + 1],
                               dof_handlers[l],
                               *mappings[l + 1],
                               *mappings[l],
                               constraints[l + 1],
                               constraints[l]);
    }

  MGTransferGlobalCoarsening<dim, VectorType> transfer(
    transfers,
    [&](const auto l, auto &vec) { operators[l].initialize_dof_vector(vec); });


  GMGParameters mg_data; // TODO

  VectorType dst, src;
  operators[max_level].initialize_dof_vector(dst);
  operators[max_level].initialize_dof_vector(src);

  operators[max_level].rhs(src);

  ReductionControl solver_control(
    mg_data.maxiter, mg_data.abstol, mg_data.reltol, false, false);

  mg_solve(solver_control,
           dst,
           src,
           mg_data,
           dof_handlers[max_level],
           operators[max_level],
           operators,
           transfer);

  deallog << dim << ' ' << fe_degree_fine << ' ' << n_refinements << ' '
          << (do_simplex_mesh ? "tri " : "quad") << ' '
          << solver_control.last_step() << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.precision(8);

  const unsigned int n_refinements = 3;
  const unsigned int degree        = 4;
  test<2>(n_refinements, degree, false /*quadrilateral*/);
}
