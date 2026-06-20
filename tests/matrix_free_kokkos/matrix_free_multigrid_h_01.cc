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
 * Test global-coarsening multigrid for a uniformly refined mesh both for
 *  hypercube mesh.
 */

#include "multigrid_util.h"

using namespace dealii;

template <int dim, unsigned int fe_degree, typename Number = double>
void
test(const unsigned int n_refinements)
{
  using VectorType =
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>;

  using MatrixType   = Portable::LaplaceOperator<dim, fe_degree, Number>;
  using SmootherType = PreconditionChebyshev<MatrixType, VectorType>;
  using TransferType = Portable::MGTwoLevelTransfer<dim, VectorType>;

  const unsigned int min_level = 0;
  const unsigned int max_level = n_refinements;

  MGLevelObject<Triangulation<dim>>        triangulations(min_level, max_level);
  MGLevelObject<DoFHandler<dim>>           dof_handlers(min_level, max_level);
  MGLevelObject<AffineConstraints<Number>> constraints(min_level, max_level);
  MGLevelObject<TransferType>              transfers(min_level, max_level);
  MGLevelObject<MatrixType>                operators(min_level, max_level);
  MGLevelObject<SmootherType>              smoothers(min_level, max_level);

  AffineConstraints<Number> constraints_fine;

  // set up levels
  for (auto l = min_level; l <= max_level; ++l)
    {
      auto &tria        = triangulations[l];
      auto &dof_handler = dof_handlers[l];
      auto &constraint  = constraints[l];
      auto &op          = operators[l];

      std::unique_ptr<FiniteElement<dim>> fe;

      fe = std::make_unique<FE_Q<dim>>(fe_degree);

      // set up triangulation
      GridGenerator::subdivided_hyper_cube(tria, 2);
      tria.refine_global(l);

      // set up dofhandler
      dof_handler.reinit(tria);
      dof_handler.distribute_dofs(*fe);

      // set up constraints
      const IndexSet locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(dof_handler);
      constraint.reinit(dof_handler.locally_owned_dofs(),
                        locally_relevant_dofs);
      VectorTools::interpolate_boundary_values(dof_handler,
                                               0,
                                               Functions::ZeroFunction<dim>(),
                                               constraint);
      constraint.close();

      if (l == max_level)
        constraints_fine.copy_from(constraint);

      // set up operator
      op.reinit(dof_handler, constraint);

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
    transfers[l + 1].reinit_geometric_transfer(dof_handlers[l + 1],
                                               dof_handlers[l],
                                               constraints[l + 1],
                                               constraints[l]);

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
