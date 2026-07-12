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


/**
 *
 * Like matrix_free_multigrid_h_02.cc but for a system of finite elements.
 */

#include "multigrid_util.h"

using namespace dealii;

template <int dim, unsigned int fe_degree, typename Number = double>
void
test(const unsigned int n_refinements)
{
  using VectorType =
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>;
  using HostVectorType =
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Host>;

  using MatrixType   = Portable::LaplaceOperator<dim, fe_degree, dim, Number>;
  using SmootherType = PreconditionChebyshev<MatrixType, VectorType>;
  using TransferType = Portable::MGTwoLevelTransfer<dim, VectorType>;

  parallel::distributed::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    Triangulation<dim>::MeshSmoothing::none,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  GridGenerator::subdivided_hyper_cube(tria, 2);
  tria.refine_global(n_refinements);

  const unsigned int min_level = 0;
  const unsigned int max_level = n_refinements;

  MGLevelObject<AffineConstraints<Number>> constraints(min_level, max_level);
  MGLevelObject<TransferType>              transfers(min_level, max_level);
  MGLevelObject<MatrixType>                operators(min_level, max_level);
  MGLevelObject<SmootherType>              smoothers(min_level, max_level);

  std::unique_ptr<FiniteElement<dim>> fe;
  fe = std::make_unique<FE_Q<dim>>(fe_degree);

  // set up dofhandler
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(FESystem<dim>(*fe, dim));
  dof_handler.distribute_mg_dofs();

  const std::set<types::boundary_id> dirichlet_boundary = {0};
  MGConstrainedDoFs                  mg_constrained_dofs;
  mg_constrained_dofs.initialize(dof_handler);
  mg_constrained_dofs.make_zero_boundary_constraints(dof_handler,
                                                     dirichlet_boundary);

  AffineConstraints<Number> constraints_fine;

  // set up levels
  for (auto l = min_level; l <= max_level; ++l)
    {
      auto &constraint = constraints[l];
      auto &op         = operators[l];

      // set up constraints
      const IndexSet relevant_dofs =
        DoFTools::extract_locally_relevant_level_dofs(dof_handler, l);
      constraint.reinit(dof_handler.locally_owned_mg_dofs(l), relevant_dofs);
      constraint.add_lines(mg_constrained_dofs.get_boundary_indices(l));
      constraint.close();

      if (l == max_level)
        constraints_fine.copy_from(constraint);

      // set up operator
      op.reinit(dof_handler, constraint, l);

      typename SmootherType::AdditionalData smoother_data;
      if (l > 0)
        {
          smoother_data.smoothing_range     = 15.;
          smoother_data.degree              = 5;
          smoother_data.eig_cg_n_iterations = 10;
        }
      else
        {
          smoother_data.smoothing_range     = 1e-3;
          smoother_data.degree              = numbers::invalid_unsigned_int;
          smoother_data.eig_cg_n_iterations = op.m();
        }

      smoother_data.constraints.copy_from(constraint);
      op.compute_diagonal();
      smoother_data.preconditioner = op.get_matrix_diagonal_inverse();

      smoothers[l].initialize(op, smoother_data);
    }

  // set up transfer operator
  for (unsigned int l = min_level; l < max_level; ++l)
    {
      transfers[l + 1].reinit_geometric_transfer(
        dof_handler, dof_handler, constraints[l + 1], constraints[l], l + 1, l);
    }

  VectorType src, dst;

  operators[max_level].initialize_dof_vector(src);
  operators[max_level].initialize_dof_vector(dst);

  operators[max_level].rhs(src, constraints_fine);

  const unsigned int maxiter = 10000;
  const double       abstol  = 1e-20;
  const double       reltol  = 1e-4;

  ReductionControl solver_control(maxiter, abstol, reltol, false, false);

  Portable::
    VCycleMultigrid<dim, MatrixType, VectorType, TransferType, SmootherType>
      mg_solve(operators, transfers, smoothers);

  mg_solve.solve_cg(solver_control, dst, src);

  deallog << dim << ' ' << fe_degree << ' ' << n_refinements << ' '
          << solver_control.last_step() << std::endl;
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.precision(8);

  constexpr unsigned int max_fe_degree = 4;


  for (unsigned int n_refinements = 2; n_refinements <= 4; ++n_refinements)
    {
      test<2, 1>(n_refinements);
      test<2, 2>(n_refinements);
      test<2, 3>(n_refinements);
      test<2, 4>(n_refinements);
    }
}
