// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// Test TrilinosPayload vmult and Tvmult operations for MPI vectors,
// specifically under conditions where the transpose flag is set
// This is the parallel version of linear_operator_13.cc

#include <deal.II/base/function.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/packaged_operation.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_linear_operator.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
void
build_matrix_vector(TrilinosWrappers::BlockSparseMatrix &matrix,
                    TrilinosWrappers::MPI::BlockVector & vector,
                    const FE_Q<dim> &                    fe_test,
                    const FE_Q<dim> &                    fe_trial)
{
  deallog.push("build_matrix_vector");

  MPI_Comm mpi_communicator(MPI_COMM_WORLD);

  // Configure block system
  // Block data
  const unsigned int        n_blocks = 2;
  std::vector<unsigned int> block_component(n_blocks);
  block_component[0] = 0;
  block_component[1] = 1;
  std::vector<types::global_dof_index> dofs_per_block(n_blocks);

  // DoF index data
  std::vector<IndexSet> all_locally_owned_dofs;
  IndexSet              locally_owned_dofs;
  IndexSet              locally_relevant_dofs;
  std::vector<IndexSet> locally_owned_partitioning;
  std::vector<IndexSet> locally_relevant_partitioning;

  // Initialise
  const FESystem<dim>                       fe(fe_test, 1, fe_trial, 1);
  parallel::distributed::Triangulation<dim> triangulation(
    mpi_communicator,
    typename Triangulation<dim>::MeshSmoothing(
      Triangulation<dim>::smoothing_on_refinement |
      Triangulation<dim>::smoothing_on_coarsening));
  QGauss<dim>               quadrature_formula(fe_trial.degree + 1);
  DoFHandler<dim>           dof_handler(triangulation);
  AffineConstraints<double> constraints;

  const unsigned int dofs_per_cell = fe.dofs_per_cell;

  // Make grid
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(2);

  // Setup system
  dof_handler.distribute_dofs(fe);
  DoFRenumbering::component_wise(dof_handler, block_component);
  dofs_per_block =
    DoFTools::count_dofs_per_fe_block(dof_handler, block_component);

  locally_owned_dofs = dof_handler.locally_owned_dofs();
  locally_owned_partitioning.push_back(
    locally_owned_dofs.get_view(0, dofs_per_block[0]));
  locally_owned_partitioning.push_back(
    locally_owned_dofs.get_view(dofs_per_block[0],
                                dofs_per_block[0] + dofs_per_block[1]));
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
  locally_relevant_partitioning.push_back(
    locally_relevant_dofs.get_view(0, dofs_per_block[0]));
  locally_relevant_partitioning.push_back(
    locally_relevant_dofs.get_view(dofs_per_block[0],
                                   dofs_per_block[0] + dofs_per_block[1]));

  constraints.clear();
  constraints.reinit(locally_relevant_dofs);
  constraints.close();

  // See
  // https://www.dealii.org/developer/doxygen/deal.II/step_32.html#TheBoussinesqFlowProblemsetupfunctions
  Table<2, DoFTools::Coupling> coupling(n_blocks, n_blocks);
  coupling.fill(DoFTools::always);
  TrilinosWrappers::BlockSparsityPattern dsp(locally_owned_partitioning,
                                             locally_owned_partitioning,
                                             locally_relevant_partitioning,
                                             mpi_communicator);
  DoFTools::make_sparsity_pattern(dof_handler,
                                  coupling,
                                  dsp,
                                  constraints,
                                  false,
                                  Utilities::MPI::this_mpi_process(
                                    mpi_communicator));
  dsp.compress();

  matrix.clear();
  matrix.reinit(dsp);
  vector.reinit(locally_owned_partitioning, mpi_communicator);

  // Assemble system: Mass matrix and constraint RHS vector
  FEValues<dim>      fe_values(fe,
                          quadrature_formula,
                          update_values | update_JxW_values);
  const unsigned int n_q_points = quadrature_formula.size();
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned() == false)
        continue;

      fe_values.reinit(cell);
      cell_matrix = 0;
      cell_rhs    = 0;

      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              // Globally symmetric contributions, but the off-diagonal
              // blocks are non-square
              // This is useful for checking implementation of transpose
              // operator
              cell_matrix(i, j) +=
                (fe_values.shape_value(i, q_point) *
                 fe_values.shape_value(j, q_point) * fe_values.JxW(q_point));

            // Non-trivial vector contribution
            cell_rhs(i) += (fe_values.shape_value(i, q_point) * 1.0 *
                            fe_values.JxW(q_point));
          }

      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(
        cell_matrix, cell_rhs, local_dof_indices, matrix, vector);
    }
  matrix.compress(VectorOperation::add);
  vector.compress(VectorOperation::add);

  deallog.pop();
}

void
evaluate_ops(const TrilinosWrappers::BlockSparseMatrix &matrix,
             const TrilinosWrappers::MPI::BlockVector & vector)
{
  const double                                   tol = 1e-12;
  typedef dealii::TrilinosWrappers::SparseMatrix MatrixType;
  typedef dealii::TrilinosWrappers::MPI::Vector  VectorType;
  typedef dealii::TrilinosWrappers::internal::LinearOperatorImplementation::
    TrilinosPayload                        PayloadType;
  typedef typename PayloadType::VectorType PayloadVectorType;
  typedef dealii::types::global_dof_index  size_type;

  deallog.push("System info");
  {
    deallog << "Matrix frobenius norm" << std::endl
            << "block(0,0): " << matrix.block(0, 0).frobenius_norm()
            << std::endl
            << "block(0,1): " << matrix.block(0, 1).frobenius_norm()
            << std::endl
            << "block(1,0): " << matrix.block(1, 0).frobenius_norm()
            << std::endl
            << "block(1,1): " << matrix.block(1, 1).frobenius_norm()
            << std::endl
            << "Vector L1 norm" << std::endl
            << "block(0): " << vector.block(0).l1_norm() << std::endl
            << "block(1): " << vector.block(1).l1_norm() << std::endl
            << "Vector L2 norm" << std::endl
            << "block(0): " << vector.block(0).l2_norm() << std::endl
            << "block(1): " << vector.block(1).l2_norm() << std::endl
            << "Vector Linfty norm" << std::endl
            << "block(0): " << vector.block(0).linfty_norm() << std::endl
            << "block(1): " << vector.block(1).linfty_norm() << std::endl
            << std::endl;
  }
  deallog.pop();

  deallog.push("Matrix TW::Vector");
  {
    {
      deallog.push("vmult");

      const MatrixType &A = matrix.block(1, 0);
      const VectorType &b = vector.block(0);
      const VectorType &r = vector.block(1);
      Assert(A.frobenius_norm() > 0.0, ExcInternalError());
      Assert(b.l2_norm() > 0.0, ExcInternalError());
      deallog << "System size: "
              << "  A.m(): " << A.m() << "  A.n(): " << A.n()
              << "  b.size(): " << b.size() << "  r.size(): " << r.size()
              << std::endl;

      VectorType out_ref(r);
      VectorType out_lo(r);
      VectorType out_lo_pyld(r);

      const auto lo_A = linear_operator<VectorType, VectorType, PayloadType>(A);

      // First test the standard operation to get a reference result
      A.vmult(out_ref, b);
      deallog << "Reference result norm squared: " << out_ref.norm_sqr()
              << std::endl;

      // Next we check the native operation of the LinearOperator
      out_lo = lo_A * b; // lo_A.vmult(out_lo,b);
      deallog << "LinearOperator result norm squared: " << out_lo.norm_sqr()
              << std::endl;
      const VectorType diff1 = out_lo - out_ref;
      Assert(
        std::sqrt(diff1.norm_sqr()) < tol,
        ExcMessage(
          "LinearOperator vmult operation does not match reference result"));

      // Lastly we test functionality added by the Payload
      // These Trilinos vectors are composed of a set of pointers to elements
      // in the deal.II itself, and so don't own their own elements.
      // So when we perform an operation on them, we perform the
      // operation on the deal.II vector itself.
      // For how the Epetra_MultiVector's are initialised, see
      // TrilinosWrappers::SparseMatrix::vmult
      out_lo_pyld                  = 0.0;
      const size_type o_local_size = out_lo_pyld.end() - out_lo_pyld.begin();
      AssertDimension(o_local_size,
                      static_cast<size_type>(
                        lo_A.OperatorRangeMap().NumMyPoints()));
      PayloadVectorType tril_out_lo_pyld(View,
                                         lo_A.OperatorRangeMap(),
                                         const_cast<TrilinosScalar *>(
                                           out_lo_pyld.begin()),
                                         o_local_size,
                                         1);
      const size_type   b_local_size = b.end() - b.begin();
      AssertDimension(b_local_size,
                      static_cast<size_type>(
                        lo_A.OperatorDomainMap().NumMyPoints()));
      PayloadVectorType tril_b_pyld(View,
                                    lo_A.OperatorDomainMap(),
                                    const_cast<TrilinosScalar *>(b.begin()),
                                    b_local_size,
                                    1);

      deallog << "LinearOperator status:"
              << "  UseTranspose: " << lo_A.UseTranspose() << std::endl;
      lo_A.Apply(tril_b_pyld, tril_out_lo_pyld);
      deallog << "LinearOperator payload result norm squared: "
              << out_lo_pyld.norm_sqr() << std::endl;
      const VectorType diff2 = out_lo_pyld - out_ref;
      Assert(
        std::sqrt(diff2.norm_sqr()) < tol,
        ExcMessage(
          "LinearOperator payload vmult operation does not match reference result"));

      deallog.pop();
    }

    {
      deallog.push("Tvmult");

      const MatrixType &A = matrix.block(1, 0);
      const VectorType &b = vector.block(1);
      const VectorType &r = vector.block(0);
      Assert(A.frobenius_norm() > 0.0, ExcInternalError());
      Assert(b.l2_norm() > 0.0, ExcInternalError());
      deallog << "System size: "
              << "  A.m(): " << A.m() << "  A.n(): " << A.n()
              << "  b.size(): " << b.size() << "  r.size(): " << r.size()
              << std::endl;

      VectorType out_ref(r);
      VectorType out_lo(r);
      VectorType out_lo_pyld(r);

      // First test the standard operation to get a reference result
      A.Tvmult(out_ref, b);
      deallog << "Reference result norm squared: " << out_ref.norm_sqr()
              << std::endl;

      // Next we check the native operation of the LinearOperator
      const auto lo_A_T =
        transpose_operator<VectorType, VectorType, PayloadType>(A);
      out_lo = lo_A_T * b; // lo_A.Tvmult(out_lo,b);
      deallog << "LinearOperator result norm squared: " << out_lo.norm_sqr()
              << std::endl;
      const VectorType diff1 = out_lo - out_ref;
      Assert(
        std::sqrt(diff1.norm_sqr()) < tol,
        ExcMessage(
          "LinearOperator Tvmult operation does not match reference result"));

      // Lastly we test functionality added by the Payload
      out_lo_pyld                  = 0.0;
      const size_type o_local_size = out_lo_pyld.end() - out_lo_pyld.begin();
      AssertDimension(o_local_size,
                      static_cast<size_type>(
                        lo_A_T.OperatorRangeMap().NumMyPoints()));
      PayloadVectorType tril_out_lo_pyld(View,
                                         lo_A_T.OperatorRangeMap(),
                                         const_cast<TrilinosScalar *>(
                                           out_lo_pyld.begin()),
                                         o_local_size,
                                         1);
      const size_type   b_local_size = b.end() - b.begin();
      AssertDimension(b_local_size,
                      static_cast<size_type>(
                        lo_A_T.OperatorDomainMap().NumMyPoints()));
      PayloadVectorType tril_b_pyld(View,
                                    lo_A_T.OperatorDomainMap(),
                                    const_cast<TrilinosScalar *>(b.begin()),
                                    b_local_size,
                                    1);

      deallog << "LinearOperator status:"
              << "  UseTranspose: " << lo_A_T.UseTranspose() << std::endl;
      lo_A_T.Apply(tril_b_pyld, tril_out_lo_pyld);
      deallog << "LinearOperator payload result norm squared: "
              << out_lo_pyld.norm_sqr() << std::endl;
      const VectorType diff2 = out_lo_pyld - out_ref;
      Assert(
        std::sqrt(diff2.norm_sqr()) < tol,
        ExcMessage(
          "LinearOperator payload Tvmult operation does not match reference result"));

      deallog.pop();
    }

    {
      deallog.push("Composite vmult");

      const MatrixType &A = matrix.block(1, 0);
      const VectorType &b = vector.block(0);
      const VectorType &r = vector.block(0);
      const VectorType &i = vector.block(1);
      Assert(A.frobenius_norm() > 0.0, ExcInternalError());
      Assert(b.l2_norm() > 0.0, ExcInternalError());
      deallog << "System size: "
              << "  A.m(): " << A.m() << "  A.n(): " << A.n()
              << "  b.size(): " << b.size() << "  r.size(): " << r.size()
              << std::endl;

      VectorType out_ref(r);
      VectorType out_lo(r);
      VectorType out_lo_pyld(r);

      // First test the standard operation to get a reference result
      {
        VectorType int_lo_pyld(i); // intermediate solution
        A.vmult(int_lo_pyld, b);
        A.Tvmult(out_ref, int_lo_pyld);
      }
      deallog << "Reference result norm squared: " << out_ref.norm_sqr()
              << std::endl;

      // Next we check the native operation of the LinearOperator
      const auto lo_A = linear_operator<VectorType, VectorType, PayloadType>(A);
      const auto lo_A_T =
        transpose_operator<VectorType, VectorType, PayloadType>(A);
      out_lo = (lo_A_T * lo_A) * b;
      deallog << "LinearOperator result norm squared: " << out_lo.norm_sqr()
              << std::endl;
      const VectorType diff1 = out_lo - out_ref;
      Assert(
        std::sqrt(diff1.norm_sqr()) < tol,
        ExcMessage(
          "LinearOperator composite vmult operation does not match reference result"));

      // Lastly we test functionality added by the Payload
      const auto lo_A_T_x_lo_A = lo_A_T * lo_A; // Construct composite operator
      out_lo_pyld              = 0.0;
      const size_type o_local_size = out_lo_pyld.end() - out_lo_pyld.begin();
      AssertDimension(o_local_size,
                      static_cast<size_type>(
                        lo_A_T_x_lo_A.OperatorRangeMap().NumMyPoints()));
      PayloadVectorType tril_out_lo_pyld(View,
                                         lo_A_T_x_lo_A.OperatorRangeMap(),
                                         const_cast<TrilinosScalar *>(
                                           out_lo_pyld.begin()),
                                         o_local_size,
                                         1);
      const size_type   b_local_size = b.end() - b.begin();
      AssertDimension(b_local_size,
                      static_cast<size_type>(
                        lo_A_T_x_lo_A.OperatorDomainMap().NumMyPoints()));
      PayloadVectorType tril_b_pyld(View,
                                    lo_A_T_x_lo_A.OperatorDomainMap(),
                                    const_cast<TrilinosScalar *>(b.begin()),
                                    b_local_size,
                                    1);

      deallog << "LinearOperator status:"
              << "  UseTranspose: " << lo_A_T_x_lo_A.UseTranspose()
              << std::endl;
      lo_A_T_x_lo_A.Apply(tril_b_pyld, tril_out_lo_pyld);
      deallog << "LinearOperator payload result norm squared: "
              << out_lo_pyld.norm_sqr() << std::endl;
      const VectorType diff2 = out_lo_pyld - out_ref;
      Assert(
        std::sqrt(diff2.norm_sqr()) < tol,
        ExcMessage(
          "LinearOperator payload composite vmult operation does not match reference result"));

      deallog.pop();
    }

    {
      deallog.push("Composite mult Tvmult");

      const MatrixType &A = matrix.block(1, 0);
      const VectorType &b = vector.block(1);
      const VectorType &r = vector.block(1);
      const VectorType &i = vector.block(0);
      Assert(A.frobenius_norm() > 0.0, ExcInternalError());
      Assert(b.l2_norm() > 0.0, ExcInternalError());
      deallog << "System size: "
              << "  A.m(): " << A.m() << "  A.n(): " << A.n()
              << "  b.size(): " << b.size() << "  r.size(): " << r.size()
              << std::endl;

      VectorType out_ref(r);
      VectorType out_lo(r);
      VectorType out_lo_pyld(r);

      // First test the standard operation to get a reference result
      {
        VectorType int_lo_pyld(i); // intermediate
        A.Tvmult(int_lo_pyld, b);
        A.vmult(out_ref, int_lo_pyld);
      }
      deallog << "Reference result norm squared: " << out_ref.norm_sqr()
              << std::endl;

      // Next we check the native operation of the LinearOperator
      const auto lo_A = linear_operator<VectorType, VectorType, PayloadType>(A);
      const auto lo_A_T =
        transpose_operator<VectorType, VectorType, PayloadType>(A);
      out_lo = (lo_A * lo_A_T) * b;
      deallog << "LinearOperator result norm squared: " << out_lo.norm_sqr()
              << std::endl;
      const VectorType diff = out_lo - out_ref;
      Assert(
        std::sqrt(diff.norm_sqr()) < tol,
        ExcMessage(
          "LinearOperator composite Tvmult operation does not match reference result"));

      // Lastly we test functionality added by the Payload
      const auto lo_A_x_lo_A_T =
        transpose_operator(lo_A * lo_A_T); // Construct composite operator
      out_lo_pyld                  = 0.0;
      const size_type o_local_size = out_lo_pyld.end() - out_lo_pyld.begin();
      AssertDimension(o_local_size,
                      static_cast<size_type>(
                        lo_A_x_lo_A_T.OperatorRangeMap().NumMyPoints()));
      PayloadVectorType tril_out_lo_pyld(View,
                                         lo_A_x_lo_A_T.OperatorRangeMap(),
                                         const_cast<TrilinosScalar *>(
                                           out_lo_pyld.begin()),
                                         o_local_size,
                                         1);
      const size_type   b_local_size = b.end() - b.begin();
      AssertDimension(b_local_size,
                      static_cast<size_type>(
                        lo_A_x_lo_A_T.OperatorDomainMap().NumMyPoints()));
      PayloadVectorType tril_b_pyld(View,
                                    lo_A_x_lo_A_T.OperatorDomainMap(),
                                    const_cast<TrilinosScalar *>(b.begin()),
                                    b_local_size,
                                    1);

      deallog << "LinearOperator status:"
              << "  UseTranspose: " << lo_A_x_lo_A_T.UseTranspose()
              << std::endl;
      lo_A_x_lo_A_T.Apply(tril_b_pyld, tril_out_lo_pyld);
      deallog << "LinearOperator payload result norm squared: "
              << out_lo_pyld.norm_sqr() << std::endl;
      const VectorType diff2 = out_lo_pyld - out_ref;
      Assert(
        std::sqrt(diff2.norm_sqr()) < tol,
        ExcMessage(
          "LinearOperator payload composite vmult operation does not match reference result"));

      deallog.pop();
    }
  }
  deallog.pop();
  deallog << "Matrix TW::Vector OK" << std::endl;
}

int
main(int argc, char *argv[])
{
  const int dim = 2;

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.depth_console(0);
  deallog << std::setprecision(10);

  FE_Q<dim> fe_test_1(1);
  FE_Q<dim> fe_test_2(2);
  FE_Q<dim> fe_trial(1);

  TrilinosWrappers::BlockSparseMatrix A;
  TrilinosWrappers::MPI::BlockVector  b;

  deallog.push("Square");
  build_matrix_vector<dim>(A, b, fe_test_1, fe_trial);
  evaluate_ops(A, b);
  deallog.pop();

  deallog << std::endl << std::endl;

  deallog.push("Non-square");
  build_matrix_vector<dim>(A, b, fe_test_2, fe_trial);
  evaluate_ops(A, b);
  deallog.pop();
}
