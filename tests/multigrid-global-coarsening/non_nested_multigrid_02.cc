// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
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
 * Test global-coarsening multigrid for a uniformly refined mesh both for
 * simplex and hypercube mesh.
 */

#include <deal.II/grid/grid_tools.h>

#include "multigrid_util.h"

template <int dim, typename Number = double>
void
test(const unsigned int n_refinements,
     const unsigned int fe_degree_fine,
     const double       factor)
{
  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  const unsigned int min_level = 0;
  const unsigned int max_level = n_refinements;

  MGLevelObject<Triangulation<dim>>        triangulations(min_level, max_level);
  MGLevelObject<DoFHandler<dim>>           dof_handlers(min_level, max_level);
  MGLevelObject<AffineConstraints<Number>> constraints(min_level, max_level);
  MGLevelObject<MappingQ1<dim>>            mappings(min_level, max_level);
  MGLevelObject<std::shared_ptr<MGTwoLevelTransferNonNested<dim, VectorType>>>
                                          transfers(min_level, max_level);
  MGLevelObject<Operator<dim, 1, Number>> operators(min_level, max_level);


  // set up levels
  for (auto l = min_level; l <= max_level; ++l)
    {
      auto &tria        = triangulations[l];
      auto &dof_handler = dof_handlers[l];
      auto &constraint  = constraints[l];
      auto &mapping     = mappings[l];
      auto &op          = operators[l];

      std::unique_ptr<FiniteElement<dim>> fe =
        std::make_unique<FE_Q<dim>>(fe_degree_fine);
      std::unique_ptr<Quadrature<dim>> quad =
        std::make_unique<QGauss<dim>>(fe_degree_fine + 1);
      std::unique_ptr<Mapping<dim>> _mapping =
        std::make_unique<MappingFE<dim>>(FE_Q<dim>(1));


      // set up triangulation
      GridGenerator::subdivided_hyper_cube(tria, 2);
      tria.refine_global(l);
      GridTools::distort_random(factor,
                                tria,
                                true,
                                boost::random::mt19937::default_seed);

      // set up dofhandler
      dof_handler.reinit(tria);
      dof_handler.distribute_dofs(*fe);

      // set up constraints
      const IndexSet locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(dof_handler);
      constraint.reinit(dof_handler.locally_owned_dofs(),
                        locally_relevant_dofs);
      VectorTools::interpolate_boundary_values(
        *_mapping, dof_handler, 0, Functions::ZeroFunction<dim>(), constraint);
      constraint.close();

      // set up operator
      op.reinit(*_mapping, dof_handler, *quad, constraint);
    }

  // set up transfer operator
  for (unsigned int l = min_level; l < max_level; ++l)
    {
      transfers[l + 1] =
        std::make_shared<MGTwoLevelTransferNonNested<dim, VectorType>>();
      transfers[l + 1]->reinit(dof_handlers[l + 1],
                               dof_handlers[l],
                               mappings[l + 1],
                               mappings[l],
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
          << "quad" << ' ' << solver_control.last_step() << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.precision(8);

  const double factor = .36;
  for (unsigned int n_refinements = 2; n_refinements <= 4; ++n_refinements)
    for (unsigned int degree = 1; degree <= 4; ++degree)
      test<2>(n_refinements, degree, factor);
}
