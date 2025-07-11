// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// test the mumps sparse direct solver in parallel on a
// TrilinosWrappers::SparseMatrix

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <cstdlib>
#include <fstream>
#include <iostream>

#include "../tests/tests.h"

template <typename MatrixType, typename VectorType>
void
solve_and_check(const MatrixType &M,
                const VectorType &rhs,
                const VectorType &solution)
{
  SparseDirectMUMPS::AdditionalData data;
  data.output_details = false;
  data.symmetric      = true;
  data.posdef         = true;
  SparseDirectMUMPS solver(data, M.get_mpi_communicator());
  solver.initialize(M);
  VectorType dst(rhs);
  solver.vmult(dst, rhs);
  dst -= solution;
  Assert(dst.l2_norm() < 1e-8, ExcInternalError());
}

template <int dim, typename MatrixType, typename VectorType>
void
test()
{
  deallog << dim << 'd' << std::endl;

  MPI_Comm mpi_communicator = MPI_COMM_WORLD;

  parallel::distributed::Triangulation<dim> tria(mpi_communicator);
  GridGenerator::hyper_cube(tria, -1, 1);
  tria.refine_global(1);

  // destroy the uniformity of the matrix by
  // refining one cell
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.refine_global(8 - 2 * dim);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  deallog << "Number of dofs = " << dof_handler.n_dofs() << std::endl;

  auto locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);
  auto                   locally_owned_dofs = dof_handler.locally_owned_dofs();
  DynamicSparsityPattern dsp(locally_relevant_dofs);
  AffineConstraints<double> constraints;

  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
  SparsityTools::distribute_sparsity_pattern(dsp,
                                             locally_owned_dofs,
                                             mpi_communicator,
                                             locally_relevant_dofs);

  MatrixType B;
  B.reinit(locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);

  QGauss<dim> qr(2);
  MatrixTools::create_mass_matrix(dof_handler, qr, B);
  B.compress(VectorOperation::add);

  SparseDirectMUMPS::AdditionalData data;
  data.output_details   = true;
  data.error_statistics = true;
  SparseDirectMUMPS Binv(data, B.get_mpi_communicator());
  Binv.initialize(B);

  // for a number of different solution
  // vectors, make up a matching rhs vector
  // and check what the solver finds
  for (unsigned int i = 0; i < 3; ++i)
    {
      VectorType solution(dof_handler.locally_owned_dofs(), mpi_communicator);
      VectorType x(dof_handler.locally_owned_dofs(), mpi_communicator);
      VectorType b(dof_handler.locally_owned_dofs(), mpi_communicator);

      for (const types::global_dof_index idx : dof_handler.locally_owned_dofs())
        solution(idx) = idx + idx * (i + 1) * (i + 1);

      solution.compress(VectorOperation::insert);

      B.vmult(b, solution);

      Binv.vmult(x, b);

      x -= solution;
      deallog << "relative norm distance = " << x.l2_norm() / solution.l2_norm()
              << std::endl;
      deallog << "absolute norms = " << x.l2_norm() << ' ' << solution.l2_norm()
              << std::endl;
      Assert(x.l2_norm() / solution.l2_norm() < 1e-8, ExcInternalError());

      // check with also the posdef option
      solve_and_check<MatrixType, VectorType>(B, b, solution);
    }
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);

  MPILogInitAll log(true);

  deallog << "Trilinos matrices" << std::endl;
  test<2, TrilinosWrappers::SparseMatrix, TrilinosWrappers::MPI::Vector>();
  test<3, TrilinosWrappers::SparseMatrix, TrilinosWrappers::MPI::Vector>();


  deallog << "PETSc matrices" << std::endl;
  test<2, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>();
  test<3, PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector>();
}
