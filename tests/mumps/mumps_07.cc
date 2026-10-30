// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2001 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// test the mumps sparse direct solver in parallel on
// TrilinosWrappers::BlockSparseMatrix and PETScWrappers::MPI::BlockSparseMatrix
// with their respective block vectors

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/petsc_block_sparse_matrix.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <cstdlib>
#include <fstream>
#include <iostream>

#include "../tests/tests.h"


// Test MUMPS with block matrices and block vectors.

template <int dim, typename MatrixType, typename VectorType>
void
test()
{
  deallog << dim << 'd' << std::endl;

  MPI_Comm mpi_communicator = MPI_COMM_WORLD;

  parallel::distributed::Triangulation<dim> tria(mpi_communicator);
  GridGenerator::hyper_cube(tria, -1, 1);
  tria.refine_global(1);

  // Destroy the uniformity of the mesh by refining one cell
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.refine_global(8 - 2 * dim);

  // Use a FESystem with 2 components to create a block structure
  const FESystem<dim> fe(FE_Q<dim>(1), 2);
  DoFHandler<dim>     dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  // Renumber in block structure
  DoFRenumbering::component_wise(dof_handler);

  const std::vector<types::global_dof_index> dofs_per_component =
    DoFTools::count_dofs_per_fe_component(dof_handler);
  const unsigned int n_u = dofs_per_component[0];
  const unsigned int n_v = dofs_per_component[1];

  deallog << "Number of dofs = " << dof_handler.n_dofs() << " (" << n_u << " + "
          << n_v << ")" << std::endl;

  const auto locally_owned_dofs = dof_handler.locally_owned_dofs();
  const auto locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  // Create block index sets
  std::vector<IndexSet> locally_owned_partitioning(2);
  std::vector<IndexSet> locally_relevant_partitioning(2);

  locally_owned_partitioning[0] = locally_owned_dofs.get_view(0, n_u);
  locally_owned_partitioning[1] = locally_owned_dofs.get_view(n_u, n_u + n_v);

  locally_relevant_partitioning[0] = locally_relevant_dofs.get_view(0, n_u);
  locally_relevant_partitioning[1] =
    locally_relevant_dofs.get_view(n_u, n_u + n_v);

  // Build block sparsity pattern
  BlockDynamicSparsityPattern dsp(2, 2);
  dsp.block(0, 0).reinit(n_u, n_u);
  dsp.block(0, 1).reinit(n_u, n_v);
  dsp.block(1, 0).reinit(n_v, n_u);
  dsp.block(1, 1).reinit(n_v, n_v);
  dsp.collect_sizes();

  AffineConstraints<double> constraints;
  constraints.close();

  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);

  SparsityTools::distribute_sparsity_pattern(dsp,
                                             locally_owned_dofs,
                                             mpi_communicator,
                                             locally_relevant_dofs);

  // Create block matrix
  MatrixType B;
  B.reinit(locally_owned_partitioning, dsp, mpi_communicator);

  QGauss<dim>   qr(2);
  FEValues<dim> fe_values(
    fe, qr, update_values | update_JxW_values | update_quadrature_points);

  const FEValuesExtractors::Scalar comp_0(0);
  const FEValuesExtractors::Scalar comp_1(1);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  const unsigned int n_q_points    = qr.size();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);
        cell_matrix = 0;

        for (unsigned int q = 0; q < n_q_points; ++q)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              cell_matrix(i, j) +=
                (fe_values[comp_0].value(i, q) * fe_values[comp_0].value(j, q) +
                 fe_values[comp_1].value(i, q) *
                   fe_values[comp_1].value(j, q)) *
                fe_values.JxW(q);

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix,
                                               local_dof_indices,
                                               B);
      }
  B.compress(VectorOperation::add);

  // Initialize MUMPS with the block matrix
  SparseDirectMUMPS::AdditionalData data;
  data.output_details = false;
  data.symmetric      = true;
  data.posdef         = true;
  SparseDirectMUMPS Binv(data, B.get_mpi_communicator());
  Binv.initialize(B);

  // same as mumps_07
  for (unsigned int i = 0; i < 3; ++i)
    {
      VectorType solution(locally_owned_partitioning, mpi_communicator);
      VectorType x(locally_owned_partitioning, mpi_communicator);
      VectorType b(locally_owned_partitioning, mpi_communicator);

      for (const types::global_dof_index idx : dof_handler.locally_owned_dofs())
        {
          if (idx < n_u)
            solution.block(0)(idx) = idx + idx * (i + 1) * (i + 1);
          else
            solution.block(1)(idx - n_u) =
              (idx - n_u) + (idx - n_u) * (i + 1) * (i + 1);
        }
      solution.compress(VectorOperation::insert);

      B.vmult(b, solution);

      Binv.vmult(x, b);

      x -= solution;
      deallog << "relative norm distance = " << x.l2_norm() / solution.l2_norm()
              << std::endl;
      deallog << "absolute norms = " << x.l2_norm() << ' ' << solution.l2_norm()
              << std::endl;
      Assert(x.l2_norm() / solution.l2_norm() < 1e-8, ExcInternalError());
    }
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);

  MPILogInitAll log(true);

  deallog << "Trilinos block matrices and vectors" << std::endl;
  test<2,
       TrilinosWrappers::BlockSparseMatrix,
       TrilinosWrappers::MPI::BlockVector>();
  test<3,
       TrilinosWrappers::BlockSparseMatrix,
       TrilinosWrappers::MPI::BlockVector>();

  deallog << "PETSc block matrices and vectors" << std::endl;
  test<2,
       PETScWrappers::MPI::BlockSparseMatrix,
       PETScWrappers::MPI::BlockVector>();
  test<3,
       PETScWrappers::MPI::BlockSparseMatrix,
       PETScWrappers::MPI::BlockVector>();
}
