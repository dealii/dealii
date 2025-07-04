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
 * Test p-multigrid so that it also works on any refinement level without
 * hanging nodes.
 */

#include "multigrid_util.h"

template <int dim, typename Number = double>
void
test(const unsigned int n_refinements,
     const unsigned int level,
     const unsigned int fe_degree_fine)
{
  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  parallel::distributed::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    Triangulation<dim>::MeshSmoothing::none,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  GridGenerator::subdivided_hyper_cube(tria, 2);
  tria.refine_global(n_refinements);

  const auto level_degrees =
    MGTransferGlobalCoarseningTools::create_polynomial_coarsening_sequence(
      fe_degree_fine,
      MGTransferGlobalCoarseningTools::PolynomialCoarseningSequenceType::
        bisect);

  const unsigned int min_level = 0;
  const unsigned int max_level = level_degrees.size() - 1;

  MGLevelObject<DoFHandler<dim>> dof_handlers(min_level, max_level, tria);
  MGLevelObject<AffineConstraints<Number>> constraints(min_level, max_level);
  MGLevelObject<MGTwoLevelTransfer<dim, VectorType>> transfers(min_level,
                                                               max_level);
  MGLevelObject<Operator<dim, 1, Number>> operators(min_level, max_level);

  std::unique_ptr<Mapping<dim>> mapping_;

  // set up levels
  for (auto l = min_level; l <= max_level; ++l)
    {
      auto &dof_handler = dof_handlers[l];
      auto &constraint  = constraints[l];
      auto &op          = operators[l];

      std::unique_ptr<FiniteElement<dim>> fe;
      std::unique_ptr<Quadrature<dim>>    quad;
      std::unique_ptr<Mapping<dim>>       mapping;

      fe      = std::make_unique<FE_Q<dim>>(level_degrees[l]);
      quad    = std::make_unique<QGauss<dim>>(level_degrees[l] + 1);
      mapping = std::make_unique<MappingFE<dim>>(FE_Q<dim>(1));

      if (l == max_level)
        mapping_ = mapping->clone();

      // set up dofhandler
      dof_handler.distribute_dofs(*fe);
      dof_handler.distribute_mg_dofs();

      // set up constraints
      const std::set<types::boundary_id> dirichlet_boundary = {0};
      MGConstrainedDoFs                  mg_constrained_dofs;
      mg_constrained_dofs.initialize(dof_handler);
      mg_constrained_dofs.make_zero_boundary_constraints(dof_handler,
                                                         dirichlet_boundary);

      const IndexSet relevant_dofs =
        DoFTools::extract_locally_relevant_level_dofs(dof_handler, level);
      constraint.reinit(dof_handler.locally_owned_mg_dofs(level),
                        relevant_dofs);
      constraint.add_lines(mg_constrained_dofs.get_boundary_indices(level));
      constraint.close();

      constraint.print(std::cout);

      // set up operator
      op.reinit(*mapping, dof_handler, *quad, constraint, level);
    }

  // set up transfer operator
  for (unsigned int l = min_level; l < max_level; ++l)
    transfers[l + 1].reinit(dof_handlers[l + 1],
                            dof_handlers[l],
                            constraints[l + 1],
                            constraints[l],
                            level,
                            level);

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

  constraints[max_level].distribute(dst);

  deallog << dim << ' ' << fe_degree_fine << ' ' << n_refinements << ' '
          << level << ' ' << dst.size() << ' ' << solver_control.last_step()
          << std::endl;

  return;

  static unsigned int counter = 0;

  MGLevelObject<VectorType> results(min_level, max_level);

  transfer.interpolate_to_mg(dof_handlers[max_level], results, dst);

  for (unsigned int l = min_level; l <= max_level; ++l)
    {
      DataOut<dim> data_out;

      data_out.attach_dof_handler(dof_handlers[l]);
      data_out.add_data_vector(
        results[l],
        "solution",
        DataOut_DoFData<dim, dim>::DataVectorType::type_dof_data);
      data_out.build_patches(*mapping_, 2);

      std::ofstream output("test." + std::to_string(dim) + "." +
                           std::to_string(counter) + "." + std::to_string(l) +
                           ".vtk");
      data_out.write_vtk(output);
    }

  counter++;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.precision(8);

  for (unsigned int n_refinements = 2; n_refinements <= 4; ++n_refinements)
    for (unsigned int degree = 2; degree <= 4; ++degree)
      for (unsigned int level = 0; level <= n_refinements; ++level)
        test<2>(n_refinements, level, degree);
}
