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

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/numbers.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparse_matrix_ez.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/matrix_scaling.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_matrix_ez.h>
#include <deal.II/lac/vector.h>

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  class SparseMatrix;
  namespace MPI
  {
    class SparseMatrix;
    class Vector;
  } // namespace MPI
} // namespace TrilinosWrappers
namespace PETScWrappers
{
  namespace MPI
  {
    class SparseMatrix;
    class Vector;
  } // namespace MPI
} // namespace PETScWrappers

#include "lac/matrix_scaling.inst"

// Check for unsupported matrix types
template <typename Matrix>
constexpr bool
is_supported_matrix()
{
  return std::is_same_v<Matrix, dealii::SparseMatrix<double>> ||
         std::is_same_v<Matrix, dealii::SparseMatrix<float>> ||
         std::is_same_v<Matrix, dealii::FullMatrix<double>> ||
         std::is_same_v<Matrix, dealii::FullMatrix<float>> ||
         std::is_same_v<Matrix, dealii::BlockSparseMatrix<double>> ||
         std::is_same_v<Matrix, dealii::BlockSparseMatrix<float>> ||
         std::is_same_v<Matrix, dealii::SparseMatrixEZ<double>> ||
         std::is_same_v<Matrix, dealii::SparseMatrixEZ<float>>;
}



MatrixScaling::MatrixScaling(const AdditionalData &control)
  : control(control)
  , row_scaling()
  , column_scaling()
  , converged(false)
{}



template <class Matrix>
void
MatrixScaling::scale_matrix(Matrix &matrix)
{
  // Check for unsupported matrix types
  if constexpr (is_supported_matrix<Matrix>())
    {
      const auto n_rows = matrix.m();
      const auto n_cols = matrix.n();

      converged = false;

      row_scaling.reinit(n_rows);
      column_scaling.reinit(n_cols);

      row_scaling    = 1.0;
      column_scaling = 1.0;

      switch (control.algorithm)
        {
          case MatrixScaling::AdditionalData::ScalingAlgorithm::sinkhorn_knopp:
            {
              sk_scaling(matrix,
                         control.sinkhorn_knopp_parameters.max_iterations);
            }
            break;
          case MatrixScaling::AdditionalData::ScalingAlgorithm::
            l1_linf_symmetry_preserving:
            {
              if (control.l1linf_parameters.start_inf_norm_steps > 0 &&
                  !converged)
                linfty_scaling(matrix,
                               control.l1linf_parameters.start_inf_norm_steps);
              if (control.l1linf_parameters.l1_norm_steps > 0 && !converged)
                l1_scaling(matrix, control.l1linf_parameters.l1_norm_steps);
              if (control.l1linf_parameters.end_inf_norm_steps > 0 &&
                  !converged)
                linfty_scaling(matrix,
                               control.l1linf_parameters.end_inf_norm_steps);
            }
            break;
        }
    }
  else if constexpr (std::is_same_v<Matrix, TrilinosWrappers::SparseMatrix> ||
                     std::is_same_v<Matrix, PETScWrappers::MPI::SparseMatrix>)
    {
      locally_owned_rows = matrix.locally_owned_range_indices();
      if (matrix.n() == matrix.m())
        locally_owned_cols = locally_owned_rows;
      else
        locally_owned_cols =
          Utilities::MPI::create_evenly_distributed_partitioning(
            matrix.get_mpi_communicator(), matrix.n());

      converged = false;
      ghost_columns.clear();
      ghost_columns.set_size(matrix.n());

      row_scaling.reinit(locally_owned_rows.n_elements());
      column_scaling.reinit(locally_owned_cols.n_elements());

      row_scaling    = 1.0;
      column_scaling = 1.0;

      // Identify ghost columns
      for (auto row = locally_owned_rows.begin();
           row != locally_owned_rows.end();
           ++row)
        {
          for (auto it = matrix.begin(*row); it != matrix.end(*row); ++it)
            {
              if (!locally_owned_cols.is_element(it->column()))
                ghost_columns.add_index(it->column());
            }
        }
      partitioner.reinit(locally_owned_cols,
                         ghost_columns,
                         matrix.get_mpi_communicator());

      ghost_column_owners =
        Utilities::MPI::compute_index_owner(locally_owned_cols,
                                            ghost_columns,
                                            matrix.get_mpi_communicator());

      switch (control.algorithm)
        {
          case MatrixScaling::AdditionalData::ScalingAlgorithm::sinkhorn_knopp:
            {
              sk_scaling(matrix,
                         control.sinkhorn_knopp_parameters.max_iterations);
            }
            break;
          case MatrixScaling::AdditionalData::ScalingAlgorithm::
            l1_linf_symmetry_preserving:
            {
              if (control.l1linf_parameters.start_inf_norm_steps > 0 &&
                  !converged)
                linfty_scaling(matrix,
                               control.l1linf_parameters.start_inf_norm_steps);
              if (control.l1linf_parameters.l1_norm_steps > 0 && !converged)
                l1_scaling(matrix, control.l1linf_parameters.l1_norm_steps);
              if (control.l1linf_parameters.end_inf_norm_steps > 0 &&
                  !converged)
                linfty_scaling(matrix,
                               control.l1linf_parameters.end_inf_norm_steps);
            }
            break;
        }
    }
  else
    {
      (void)matrix;
      Assert(false, ExcNotImplemented());
    }
}



template <class Matrix, class VectorType>
void
MatrixScaling::scale_linear_system(Matrix &matrix, VectorType &rhs)
{
  AssertDimension(matrix.m(), rhs.size());
  AssertDimension(matrix.m(), matrix.n());

  if constexpr (is_supported_matrix<Matrix>() ||
                std::is_same_v<Matrix, TrilinosWrappers::SparseMatrix> ||
                std::is_same_v<Matrix, PETScWrappers::MPI::SparseMatrix>)
    scale_matrix(matrix);
  else
    Assert(false, ExcNotImplemented());

  if constexpr (std::is_same_v<VectorType, dealii::Vector<double>> ||
                std::is_same_v<VectorType, dealii::Vector<float>>)
    rhs.scale(row_scaling);
  else if constexpr (std::is_same_v<VectorType, dealii::BlockVector<double>> ||
                     std::is_same_v<VectorType, dealii::BlockVector<float>>)
    {
      for (unsigned int i = 0; i < rhs.size(); ++i)
        rhs[i] *= row_scaling[i];
    }
  else if constexpr (std::is_same_v<VectorType,
                                    TrilinosWrappers::MPI::Vector> ||
                     std::is_same_v<VectorType, PETScWrappers::MPI::Vector>)
    {
      Assert(matrix.locally_owned_range_indices() ==
               rhs.locally_owned_elements(),
             ExcMessage("Matrix and vector must have the same partitioning"));
      std::vector<double> local_updates(row_scaling.size());
      for (auto i : locally_owned_rows)
        {
          auto local_idx           = partitioner.global_to_local(i);
          local_updates[local_idx] = rhs[i] * row_scaling[local_idx];
        }
      rhs.set(locally_owned_rows.get_index_vector(), local_updates);
      rhs.compress(VectorOperation::insert);
    }
  else
    Assert(false, ExcNotImplemented());
}



template <class VectorType>
void
MatrixScaling::scale_system_solution(VectorType &sol) const
{
  if constexpr (std::is_same_v<VectorType, dealii::Vector<double>> ||
                std::is_same_v<VectorType, dealii::Vector<float>>)
    {
      AssertDimension(sol.size(), column_scaling.size());
      sol.scale(column_scaling);
    }
  else if constexpr (std::is_same_v<VectorType, dealii::BlockVector<double>> ||
                     std::is_same_v<VectorType, dealii::BlockVector<float>>)
    {
      AssertDimension(sol.size(), column_scaling.size());
      for (unsigned int i = 0; i < sol.size(); ++i)
        sol[i] *= column_scaling[i];
    }
  else if constexpr (std::is_same_v<VectorType,
                                    TrilinosWrappers::MPI::Vector> ||
                     std::is_same_v<VectorType, PETScWrappers::MPI::Vector>)
    {
      Assert(locally_owned_cols == sol.locally_owned_elements(),
             ExcMessage("Matrix and vector must have the same partitioning"));
      std::vector<double> local_updates(column_scaling.size());
      for (auto i : locally_owned_cols)
        {
          auto local_idx           = partitioner.global_to_local(i);
          local_updates[local_idx] = sol[i] * column_scaling[local_idx];
        }
      sol.set(locally_owned_cols.get_index_vector(), local_updates);
      sol.compress(VectorOperation::insert);
    }
  else
    Assert(false, ExcNotImplemented());
}



const Vector<double> &
MatrixScaling::get_row_scaling() const
{
  return row_scaling;
}



const Vector<double> &
MatrixScaling::get_column_scaling() const
{
  return column_scaling;
}



bool
MatrixScaling::is_converged() const
{
  return converged;
}



template <typename Number>
bool
MatrixScaling::check_convergence(Vector<Number> row_col_norm,
                                 std::string    norm_type) const
{
  row_col_norm.add(-1.0);
  if (norm_type == "l1")
    return row_col_norm.l1_norm() < control.scaling_tolerance;
  else if (norm_type == "linfty")
    return row_col_norm.linfty_norm() < control.scaling_tolerance;
  else
    Assert(false, ExcNotImplemented());
}



template <typename Number>
bool
MatrixScaling::check_convergence(Vector<Number> row_norm,
                                 Vector<Number> col_norm,
                                 std::string    norm_type) const
{
  row_norm.add(-1.0);
  col_norm.add(-1.0);
  if (norm_type == "l1")
    return (row_norm.l1_norm() < control.scaling_tolerance &&
            col_norm.l1_norm() < control.scaling_tolerance);
  else if (norm_type == "linfty")
    return (row_norm.linfty_norm() < control.scaling_tolerance &&
            col_norm.linfty_norm() < control.scaling_tolerance);
  else
    Assert(false, ExcNotImplemented());
}



#ifdef DEAL_II_WITH_MPI
bool
MatrixScaling::check_convergence(Vector<double>    local_row_col_norm,
                                 const std::string norm_type,
                                 const MPI_Comm    mpi_communicator) const
{
  local_row_col_norm.add(-1.0);
  bool local_not_converged;

  if (norm_type == "l1")
    local_not_converged =
      !(local_row_col_norm.l1_norm() < control.scaling_tolerance);
  else if (norm_type == "linfty")
    local_not_converged =
      !(local_row_col_norm.linfty_norm() < control.scaling_tolerance);

  bool any_not_converged =
    Utilities::MPI::logical_or(local_not_converged, mpi_communicator);

  return !any_not_converged;
}



bool
MatrixScaling::check_convergence(Vector<double>    local_row_norm,
                                 Vector<double>    local_col_norm,
                                 const std::string norm_type,
                                 const MPI_Comm    mpi_communicator) const
{
  local_row_norm.add(-1.0);
  local_col_norm.add(-1.0);
  bool local_not_converged;

  if (norm_type == "l1")
    local_not_converged =
      !(local_row_norm.l1_norm() < control.scaling_tolerance &&
        local_col_norm.l1_norm() < control.scaling_tolerance);
  else if (norm_type == "linfty")
    local_not_converged =
      !(local_row_norm.linfty_norm() < control.scaling_tolerance &&
        local_col_norm.linfty_norm() < control.scaling_tolerance);

  bool any_not_converged =
    Utilities::MPI::logical_or(local_not_converged, mpi_communicator);

  return !any_not_converged;
}



void
MatrixScaling::send_prepare_col_norms(
  std::map<unsigned int,
           std::vector<std::pair<types::global_dof_index, double>>> &send_data,
  const std::map<types::global_dof_index, double> &partial_column_norms,
  Vector<double>                                  &local_col_norms)
{
  for (const auto &[col_idx, norm_value] : partial_column_norms)
    {
      if (locally_owned_cols.is_element(col_idx))
        {
          unsigned int local_idx     = partitioner.global_to_local(col_idx);
          local_col_norms[local_idx] = norm_value;
        }
      else
        {
          auto         ghost_pos   = ghost_columns.index_within_set(col_idx);
          unsigned int target_rank = ghost_column_owners[ghost_pos];

          send_data[target_rank].emplace_back(col_idx, norm_value);
        }
    }
}



void
MatrixScaling::send_prepare_updated_col_norms(
  std::map<unsigned int,
           std::vector<std::pair<types::global_dof_index, double>>>
    &send_column_norms,
  const std::map<unsigned int,
                 std::vector<std::pair<types::global_dof_index, double>>>
                       &received_data,
  const Vector<double> &local_col_norms)
{
  // For each rank that requested column norm data from me,
  // send them the corresponding column norms updated
  for (const auto &[sender_rank, pairs] : received_data)
    {
      for (const auto &[global_col, norm_contribution] : pairs)
        {
          // I received norm data for global_col, so sender_rank needs
          // my scaling value for global_col
          if (locally_owned_cols.is_element(global_col))
            {
              unsigned int local_idx = partitioner.global_to_local(global_col);
              send_column_norms[sender_rank].emplace_back(
                global_col, local_col_norms[local_idx]);
            }
        }
    }
}
#endif



template <class Matrix>
void
MatrixScaling::l1_scaling(Matrix &matrix, const unsigned int nsteps)
{
  if constexpr (is_supported_matrix<Matrix>())
    {
      using Number = typename Matrix::value_type;

      Vector<Number> row_norms(matrix.m());
      Vector<Number> col_norms(matrix.n());

      for (unsigned int i = 0; i < nsteps; i++)
        {
          row_norms = 0;
          col_norms = 0;

          for (unsigned int row = 0; row < matrix.m(); ++row)
            {
              for (auto it = matrix.begin(row); it != matrix.end(row); ++it)
                {
                  row_norms[row] += std::abs(it->value());
                  col_norms[it->column()] += std::abs(it->value());
                }
            }

          if (check_convergence(row_norms, col_norms, "l1"))
            {
              converged = true;
              break;
            }

          for (unsigned int row = 0; row < matrix.m(); ++row)
            {
              for (auto it = matrix.begin(row); it != matrix.end(row); ++it)
                {
                  matrix.set(row,
                             it->column(),
                             it->value() / std::sqrt(row_norms[row] *
                                                     col_norms[it->column()]));
                }
              row_scaling[row] /= std::sqrt(row_norms[row]);
            }
          for (unsigned int col = 0; col < matrix.n(); ++col)
            {
              column_scaling[col] /= std::sqrt(col_norms[col]);
            }
        }
    }
  else if constexpr (std::is_same_v<Matrix, TrilinosWrappers::SparseMatrix> ||
                     std::is_same_v<Matrix, PETScWrappers::MPI::SparseMatrix>)
    {
      Vector<double> local_row_norms(locally_owned_rows.n_elements());
      Vector<double> local_col_norms(locally_owned_cols.n_elements());
      std::map<types::global_dof_index, double>
        partial_column_norms; // column -> local norm

      for (unsigned int i = 0; i < nsteps; i++)
        {
          local_row_norms = 0;
          local_col_norms = 0;
          partial_column_norms.clear();

          for (auto row = locally_owned_rows.begin();
               row != locally_owned_rows.end();
               ++row)
            {
              auto local_row_idx = partitioner.global_to_local(*row);
              for (auto it = matrix.begin(*row); it != matrix.end(*row); ++it)
                {
                  local_row_norms[local_row_idx] += std::abs(it->value());
                  partial_column_norms[it->column()] += std::abs(it->value());
                }
            }
          // Communicate partial column norms
          std::map<unsigned int,
                   std::vector<std::pair<types::global_dof_index, double>>>
            send_data;

          send_prepare_col_norms(send_data,
                                 partial_column_norms,
                                 local_col_norms);

          auto received_data =
            Utilities::MPI::some_to_some(matrix.get_mpi_communicator(),
                                         send_data);

          // Process received data and fill local_col_norms
          for (const auto &[sender_rank, pairs] : received_data)
            {
              for (const auto &[global_col, contribution] : pairs)
                {
                  unsigned int local_idx =
                    partitioner.global_to_local(global_col);
                  local_col_norms[local_idx] += contribution;
                }
            }

          if (check_convergence(local_row_norms,
                                local_col_norms,
                                "l1",
                                matrix.get_mpi_communicator()))
            {
              converged = true;
              break;
            }

          for (unsigned int i = 0; i < local_row_norms.size(); ++i)
            row_scaling[i] /= std::sqrt(local_row_norms[i]);
          for (unsigned int i = 0; i < local_col_norms.size(); ++i)
            column_scaling[i] /= std::sqrt(local_col_norms[i]);

          // Communicate column norm values to all ranks that need them
          std::map<unsigned int,
                   std::vector<std::pair<types::global_dof_index, double>>>
            send_column_norms;

          send_prepare_updated_col_norms(send_column_norms,
                                         received_data,
                                         local_col_norms);

          auto received_column_norms =
            Utilities::MPI::some_to_some(matrix.get_mpi_communicator(),
                                         send_column_norms);

          std::map<types::global_dof_index, double> ghost_column_norms_lookup;
          for (const auto &[sender_rank, pairs] : received_column_norms)
            {
              for (const auto &[col_id, norm_val] : pairs)
                {
                  ghost_column_norms_lookup[col_id] = norm_val;
                }
            }

          // scale the matrix
          for (auto row = locally_owned_rows.begin();
               row != locally_owned_rows.end();
               ++row)
            {
              auto local_row_idx = partitioner.global_to_local(*row);
              for (auto it = matrix.begin(*row); it != matrix.end(*row); ++it)
                {
                  if (locally_owned_cols.is_element(it->column()))
                    {
                      unsigned int local_col_idx =
                        partitioner.global_to_local(it->column());
                      matrix.set(*row,
                                 it->column(),
                                 it->value() /
                                   std::sqrt(local_col_norms[local_col_idx] *
                                             local_row_norms[local_row_idx]));
                    }
                  else
                    {
                      double col_norms =
                        ghost_column_norms_lookup[it->column()];

                      matrix.set(*row,
                                 it->column(),
                                 it->value() /
                                   std::sqrt(col_norms *
                                             local_row_norms[local_row_idx]));
                    }
                  // This is atrocious, PETSc needs it, this NEEDS to be changed
                  if constexpr (std::is_same_v<
                                  Matrix,
                                  PETScWrappers::MPI::SparseMatrix>)
                    matrix.compress(VectorOperation::insert);
                }
            }
          if constexpr (std::is_same_v<Matrix, TrilinosWrappers::SparseMatrix>)
            matrix.compress(VectorOperation::insert);
        }
    }
}



template <class Matrix>
void
MatrixScaling::linfty_scaling(Matrix &matrix, const unsigned int nsteps)
{
  if constexpr (is_supported_matrix<Matrix>())
    {
      using Number = typename Matrix::value_type;

      Vector<Number> row_norms(matrix.m());
      Vector<Number> col_norms(matrix.n());

      for (unsigned int i = 0; i < nsteps; i++)
        {
          row_norms = 0;
          col_norms = 0;

          for (unsigned int row = 0; row < matrix.m(); ++row)
            {
              for (auto it = matrix.begin(row); it != matrix.end(row); ++it)
                {
                  row_norms[row] =
                    std::max(row_norms[row], std::abs(it->value()));
                  col_norms[it->column()] =
                    std::max(col_norms[it->column()], std::abs(it->value()));
                }
            }

          if (check_convergence(row_norms, col_norms, "linfty"))
            {
              converged = true;
              break;
            }

          for (unsigned int row = 0; row < matrix.m(); ++row)
            {
              for (auto it = matrix.begin(row); it != matrix.end(row); ++it)
                {
                  matrix.set(row,
                             it->column(),
                             it->value() / std::sqrt(row_norms[row] *
                                                     col_norms[it->column()]));
                }
              row_scaling[row] /= std::sqrt(row_norms[row]);
            }
          for (unsigned int col = 0; col < matrix.n(); ++col)
            {
              column_scaling[col] /= std::sqrt(col_norms[col]);
            }
        }
    }
  else if constexpr (std::is_same_v<Matrix, TrilinosWrappers::SparseMatrix> ||
                     std::is_same_v<Matrix, PETScWrappers::MPI::SparseMatrix>)
    {
      Vector<double> local_row_norms(locally_owned_rows.n_elements());
      Vector<double> local_col_norms(locally_owned_cols.n_elements());
      std::map<types::global_dof_index, double>
        partial_column_norms; // column -> local norm

      for (unsigned int i = 0; i < nsteps; i++)
        {
          local_row_norms = 0;
          local_col_norms = 0;
          partial_column_norms.clear();

          for (auto row = locally_owned_rows.begin();
               row != locally_owned_rows.end();
               ++row)
            {
              auto local_row_idx = partitioner.global_to_local(*row);
              for (auto it = matrix.begin(*row); it != matrix.end(*row); ++it)
                {
                  local_row_norms[local_row_idx] =
                    std::max(local_row_norms[local_row_idx],
                             std::abs(it->value()));
                  partial_column_norms[it->column()] =
                    std::max(partial_column_norms[it->column()],
                             std::abs(it->value()));
                }
            }
          // Communicate partial column norms
          std::map<unsigned int,
                   std::vector<std::pair<types::global_dof_index, double>>>
            send_data;

          send_prepare_col_norms(send_data,
                                 partial_column_norms,
                                 local_col_norms);

          auto received_data =
            Utilities::MPI::some_to_some(matrix.get_mpi_communicator(),
                                         send_data);

          // Process received data and fill local_col_norms
          for (const auto &[sender_rank, pairs] : received_data)
            {
              for (const auto &[global_col, contribution] : pairs)
                {
                  unsigned int local_idx =
                    partitioner.global_to_local(global_col);
                  local_col_norms[local_idx] =
                    std::max(local_col_norms[local_idx], contribution);
                }
            }

          if (check_convergence(local_row_norms,
                                local_col_norms,
                                "linfty",
                                matrix.get_mpi_communicator()))
            {
              converged = true;
              break;
            }

          for (unsigned int i = 0; i < local_row_norms.size(); ++i)
            row_scaling[i] /= std::sqrt(local_row_norms[i]);
          for (unsigned int i = 0; i < local_col_norms.size(); ++i)
            column_scaling[i] /= std::sqrt(local_col_norms[i]);

          // Communicate column norm values to all ranks that need them
          std::map<unsigned int,
                   std::vector<std::pair<types::global_dof_index, double>>>
            send_column_norms;

          send_prepare_updated_col_norms(send_column_norms,
                                         received_data,
                                         local_col_norms);

          auto received_column_norms =
            Utilities::MPI::some_to_some(matrix.get_mpi_communicator(),
                                         send_column_norms);

          std::map<types::global_dof_index, double> ghost_column_norms_lookup;
          for (const auto &[sender_rank, pairs] : received_column_norms)
            {
              for (const auto &[col_id, norm_val] : pairs)
                {
                  ghost_column_norms_lookup[col_id] = norm_val;
                }
            }

          // scale the matrix
          for (auto row = locally_owned_rows.begin();
               row != locally_owned_rows.end();
               ++row)
            {
              auto local_row_idx = partitioner.global_to_local(*row);
              for (auto it = matrix.begin(*row); it != matrix.end(*row); ++it)
                {
                  if (locally_owned_cols.is_element(it->column()))
                    {
                      unsigned int local_col_idx =
                        partitioner.global_to_local(it->column());
                      matrix.set(*row,
                                 it->column(),
                                 it->value() /
                                   std::sqrt(local_col_norms[local_col_idx] *
                                             local_row_norms[local_row_idx]));
                    }
                  else
                    {
                      double col_norms =
                        ghost_column_norms_lookup[it->column()];

                      matrix.set(*row,
                                 it->column(),
                                 it->value() /
                                   std::sqrt(col_norms *
                                             local_row_norms[local_row_idx]));
                    }
                  // This is atrocious, PETSc needs it, this NEEDS to be changed
                  if constexpr (std::is_same_v<
                                  Matrix,
                                  PETScWrappers::MPI::SparseMatrix>)
                    matrix.compress(VectorOperation::insert);
                }
            }
          if constexpr (std::is_same_v<Matrix, TrilinosWrappers::SparseMatrix>)
            matrix.compress(VectorOperation::insert);
        }
    }
}



template <class Matrix>
void
MatrixScaling::sk_scaling(Matrix &matrix, const unsigned int nsteps)
{
  if constexpr (is_supported_matrix<Matrix>())
    {
      using Number = typename Matrix::value_type;

      Vector<Number> row_norms(matrix.m());
      Vector<Number> col_norms(matrix.n());

      switch (control.sinkhorn_knopp_parameters.norm_type)
        {
          case MatrixScaling::AdditionalData::SKParameters::NormType::l1:
            {
              // Row_norms_0 to start the procedure
              row_norms = 0;
              for (unsigned int row = 0; row < matrix.m(); ++row)
                for (auto it = matrix.begin(row); it != matrix.end(row); ++it)
                  row_norms[row] += std::abs(it->value());

              for (unsigned int i = 0; i < nsteps; i++)
                {
                  // Row step
                  col_norms = 0;
                  for (unsigned int row = 0; row < matrix.m(); ++row)
                    {
                      for (auto it = matrix.begin(row); it != matrix.end(row);
                           ++it)
                        {
                          matrix.set(row,
                                     it->column(),
                                     it->value() / row_norms[row]);
                          col_norms[it->column()] += std::abs(it->value());
                        }
                      row_scaling[row] /= row_norms[row];
                    }

                  if (check_convergence(col_norms, "l1"))
                    {
                      converged = true;
                      break;
                    }

                  // Column step
                  row_norms = 0;
                  for (unsigned int row = 0; row < matrix.m(); ++row)
                    {
                      for (auto it = matrix.begin(row); it != matrix.end(row);
                           ++it)
                        {
                          matrix.set(row,
                                     it->column(),
                                     it->value() / col_norms[it->column()]);
                          row_norms[row] += std::abs(it->value());
                        }
                    }
                  for (unsigned int col = 0; col < matrix.n(); ++col)
                    {
                      column_scaling[col] /= col_norms[col];
                    }

                  if (check_convergence(row_norms, "l1"))
                    {
                      converged = true;
                      break;
                    }
                }
            }
            break;
          case MatrixScaling::AdditionalData::SKParameters::NormType::
            l_infinity:
            {
              // Row_norms_0 to start the procedure
              row_norms = 0;
              for (unsigned int row = 0; row < matrix.m(); ++row)
                for (auto it = matrix.begin(row); it != matrix.end(row); ++it)
                  row_norms[row] +=
                    std::max(row_norms[row], std::abs(it->value()));

              for (unsigned int i = 0; i < nsteps; i++)
                {
                  // Row step
                  col_norms = 0;
                  for (unsigned int row = 0; row < matrix.m(); ++row)
                    {
                      for (auto it = matrix.begin(row); it != matrix.end(row);
                           ++it)
                        {
                          matrix.set(row,
                                     it->column(),
                                     it->value() / row_norms[row]);
                          col_norms[it->column()] =
                            std::max(col_norms[it->column()],
                                     std::abs(it->value()));
                        }
                      row_scaling[row] /= row_norms[row];
                    }

                  if (check_convergence(col_norms, "linfty"))
                    {
                      converged = true;
                      break;
                    }

                  // Column step
                  row_norms = 0;
                  for (unsigned int row = 0; row < matrix.m(); ++row)
                    {
                      for (auto it = matrix.begin(row); it != matrix.end(row);
                           ++it)
                        {
                          matrix.set(row,
                                     it->column(),
                                     it->value() / col_norms[it->column()]);
                          row_norms[row] +=
                            std::max(row_norms[row], std::abs(it->value()));
                        }
                    }
                  for (unsigned int col = 0; col < matrix.n(); ++col)
                    {
                      column_scaling[col] /= col_norms[col];
                    }

                  if (check_convergence(row_norms, "linfty"))
                    {
                      converged = true;
                      break;
                    }
                }
            }
            break;
        }
    }
  else if constexpr (std::is_same_v<Matrix, TrilinosWrappers::SparseMatrix> ||
                     std::is_same_v<Matrix, PETScWrappers::MPI::SparseMatrix>)
    {
      Vector<double> local_row_norms(locally_owned_rows.n_elements());
      Vector<double> local_col_norms(locally_owned_cols.n_elements());
      std::map<types::global_dof_index, double>
        partial_column_norms; // column -> local norm

      switch (control.sinkhorn_knopp_parameters.norm_type)
        {
          case MatrixScaling::AdditionalData::SKParameters::NormType::l1:
            {
              // Row_norms_0 to start the procedure
              local_row_norms = 0;
              for (auto row = locally_owned_rows.begin();
                   row != locally_owned_rows.end();
                   ++row)
                {
                  auto local_row_idx = partitioner.global_to_local(*row);
                  for (auto it = matrix.begin(*row); it != matrix.end(*row);
                       ++it)
                    local_row_norms[local_row_idx] += std::abs(it->value());
                }

              for (unsigned int i = 0; i < nsteps; i++)
                {
                  // Row step
                  local_col_norms = 0;
                  partial_column_norms.clear();

                  for (auto row = locally_owned_rows.begin();
                       row != locally_owned_rows.end();
                       ++row)
                    {
                      auto local_row_idx = partitioner.global_to_local(*row);
                      for (auto it = matrix.begin(*row); it != matrix.end(*row);
                           ++it)
                        {
                          partial_column_norms[it->column()] += std::abs(
                            it->value() / local_row_norms[local_row_idx]);
                          matrix.set(*row,
                                     it->column(),
                                     it->value() /
                                       local_row_norms[local_row_idx]);
                          // This is atrocious, PETSc needs it, this NEEDS to be
                          // changed
                          if constexpr (std::is_same_v<
                                          Matrix,
                                          PETScWrappers::MPI::SparseMatrix>)
                            matrix.compress(VectorOperation::insert);
                        }
                      row_scaling[local_row_idx] /=
                        local_row_norms[local_row_idx];
                    }
                  if constexpr (std::is_same_v<Matrix,
                                               TrilinosWrappers::SparseMatrix>)
                    matrix.compress(VectorOperation::insert);

                  // Communicate partial column norms
                  std::map<
                    unsigned int,
                    std::vector<std::pair<types::global_dof_index, double>>>
                    send_data;

                  send_prepare_col_norms(send_data,
                                         partial_column_norms,
                                         local_col_norms);

                  auto received_data =
                    Utilities::MPI::some_to_some(matrix.get_mpi_communicator(),
                                                 send_data);

                  // Process received data and fill local_col_norms
                  for (const auto &[sender_rank, pairs] : received_data)
                    {
                      for (const auto &[global_col, contribution] : pairs)
                        {
                          unsigned int local_idx =
                            partitioner.global_to_local(global_col);
                          local_col_norms[local_idx] += contribution;
                        }
                    }

                  // Convergence check only on columns
                  if (check_convergence(local_col_norms,
                                        "l1",
                                        matrix.get_mpi_communicator()))
                    {
                      converged = true;
                      break;
                    }

                  // Column step
                  local_row_norms = 0;

                  // Communicate column norm values to all ranks that need them
                  std::map<
                    unsigned int,
                    std::vector<std::pair<types::global_dof_index, double>>>
                    send_column_norms;

                  send_prepare_updated_col_norms(send_column_norms,
                                                 received_data,
                                                 local_col_norms);

                  auto received_column_norms =
                    Utilities::MPI::some_to_some(matrix.get_mpi_communicator(),
                                                 send_column_norms);

                  std::map<types::global_dof_index, double>
                    ghost_column_norms_lookup;
                  for (const auto &[sender_rank, pairs] : received_column_norms)
                    {
                      for (const auto &[col_id, norm_val] : pairs)
                        {
                          ghost_column_norms_lookup[col_id] = norm_val;
                        }
                    }

                  for (auto row = locally_owned_rows.begin();
                       row != locally_owned_rows.end();
                       ++row)
                    {
                      auto local_row_idx = partitioner.global_to_local(*row);
                      for (auto it = matrix.begin(*row); it != matrix.end(*row);
                           ++it)
                        {
                          if (locally_owned_cols.is_element(it->column()))
                            {
                              unsigned int local_col_idx =
                                partitioner.global_to_local(it->column());

                              local_row_norms[local_row_idx] += std::abs(
                                it->value() / local_col_norms[local_col_idx]);

                              matrix.set(*row,
                                         it->column(),
                                         it->value() /
                                           local_col_norms[local_col_idx]);
                            }
                          else
                            {
                              double col_norms =
                                ghost_column_norms_lookup[it->column()];

                              local_row_norms[local_row_idx] +=
                                std::abs(it->value() / col_norms);

                              matrix.set(*row,
                                         it->column(),
                                         it->value() / col_norms);
                            }
                          // This is atrocious, PETSc needs it, this NEEDS to be
                          // changed
                          if constexpr (std::is_same_v<
                                          Matrix,
                                          PETScWrappers::MPI::SparseMatrix>)
                            matrix.compress(VectorOperation::insert);
                        }
                    }
                  if constexpr (std::is_same_v<Matrix,
                                               TrilinosWrappers::SparseMatrix>)
                    matrix.compress(VectorOperation::insert);
                  for (unsigned int i = 0; i < local_col_norms.size(); ++i)
                    column_scaling[i] /= local_col_norms[i];

                  if (check_convergence(local_row_norms,
                                        "l1",
                                        matrix.get_mpi_communicator()))
                    {
                      converged = true;
                      break;
                    }
                }
            }
            break;
          case MatrixScaling::AdditionalData::SKParameters::NormType::
            l_infinity:
            {
              // Row_norms_0 to start the procedure
              local_row_norms = 0;
              for (auto row = locally_owned_rows.begin();
                   row != locally_owned_rows.end();
                   ++row)
                {
                  auto local_row_idx = partitioner.global_to_local(*row);
                  for (auto it = matrix.begin(*row); it != matrix.end(*row);
                       ++it)
                    local_row_norms[local_row_idx] =
                      std::max(local_row_norms[local_row_idx],
                               std::abs(it->value()));
                }

              for (unsigned int i = 0; i < nsteps; i++)
                {
                  // Row step
                  local_col_norms = 0;
                  partial_column_norms.clear();

                  for (auto row = locally_owned_rows.begin();
                       row != locally_owned_rows.end();
                       ++row)
                    {
                      auto local_row_idx = partitioner.global_to_local(*row);
                      for (auto it = matrix.begin(*row); it != matrix.end(*row);
                           ++it)
                        {
                          partial_column_norms[it->column()] =
                            std::max(partial_column_norms[it->column()],
                                     std::abs(it->value() /
                                              local_row_norms[local_row_idx]));
                          matrix.set(*row,
                                     it->column(),
                                     it->value() /
                                       local_row_norms[local_row_idx]);

                          // This is atrocious, PETSc needs it, this NEEDS to be
                          // changed
                          if constexpr (std::is_same_v<
                                          Matrix,
                                          PETScWrappers::MPI::SparseMatrix>)
                            matrix.compress(VectorOperation::insert);
                        }
                      row_scaling[local_row_idx] /=
                        local_row_norms[local_row_idx];
                    }
                  if constexpr (std::is_same_v<Matrix,
                                               TrilinosWrappers::SparseMatrix>)
                    matrix.compress(VectorOperation::insert);

                  // Communicate partial column norms
                  std::map<
                    unsigned int,
                    std::vector<std::pair<types::global_dof_index, double>>>
                    send_data;

                  send_prepare_col_norms(send_data,
                                         partial_column_norms,
                                         local_col_norms);

                  auto received_data =
                    Utilities::MPI::some_to_some(matrix.get_mpi_communicator(),
                                                 send_data);

                  // Process received data and fill local_col_norms
                  for (const auto &[sender_rank, pairs] : received_data)
                    {
                      for (const auto &[global_col, contribution] : pairs)
                        {
                          unsigned int local_idx =
                            partitioner.global_to_local(global_col);
                          local_col_norms[local_idx] =
                            std::max(local_col_norms[local_idx], contribution);
                        }
                    }

                  if (check_convergence(local_col_norms,
                                        "linfty",
                                        matrix.get_mpi_communicator()))
                    {
                      converged = true;
                      break;
                    }

                  // Column step
                  local_row_norms = 0;

                  // Communicate column norm values to all ranks that need them
                  std::map<
                    unsigned int,
                    std::vector<std::pair<types::global_dof_index, double>>>
                    send_column_norms;

                  send_prepare_updated_col_norms(send_column_norms,
                                                 received_data,
                                                 local_col_norms);

                  auto received_column_norms =
                    Utilities::MPI::some_to_some(matrix.get_mpi_communicator(),
                                                 send_column_norms);

                  std::map<types::global_dof_index, double>
                    ghost_column_norms_lookup;
                  for (const auto &[sender_rank, pairs] : received_column_norms)
                    {
                      for (const auto &[col_id, norm_val] : pairs)
                        {
                          ghost_column_norms_lookup[col_id] = norm_val;
                        }
                    }

                  for (auto row = locally_owned_rows.begin();
                       row != locally_owned_rows.end();
                       ++row)
                    {
                      auto local_row_idx = partitioner.global_to_local(*row);
                      for (auto it = matrix.begin(*row); it != matrix.end(*row);
                           ++it)
                        {
                          if (locally_owned_cols.is_element(it->column()))
                            {
                              unsigned int local_col_idx =
                                partitioner.global_to_local(it->column());

                              local_row_norms[local_row_idx] =
                                std::max(local_row_norms[local_row_idx],
                                         std::abs(
                                           it->value() /
                                           local_col_norms[local_col_idx]));
                              matrix.set(*row,
                                         it->column(),
                                         it->value() /
                                           local_col_norms[local_col_idx]);
                            }
                          else
                            {
                              double col_norms =
                                ghost_column_norms_lookup[it->column()];

                              local_row_norms[local_row_idx] =
                                std::max(local_row_norms[local_row_idx],
                                         std::abs(it->value() / col_norms));

                              matrix.set(*row,
                                         it->column(),
                                         it->value() / col_norms);
                            }
                          // This is atrocious, PETSc needs it, this NEEDS to be
                          // changed
                          if constexpr (std::is_same_v<
                                          Matrix,
                                          PETScWrappers::MPI::SparseMatrix>)
                            matrix.compress(VectorOperation::insert);
                        }
                    }
                  if constexpr (std::is_same_v<Matrix,
                                               TrilinosWrappers::SparseMatrix>)
                    matrix.compress(VectorOperation::insert);
                  for (unsigned int i = 0; i < local_col_norms.size(); ++i)
                    column_scaling[i] /= local_col_norms[i];

                  if (check_convergence(local_row_norms,
                                        "linfty",
                                        matrix.get_mpi_communicator()))
                    {
                      converged = true;
                      break;
                    }
                }
            }
            break;
        }
    }
}



#define InstantiateMatrixScaling(MATRIX, VECTOR)                        \
  template void MatrixScaling::scale_matrix(MATRIX &);                  \
  template void MatrixScaling::scale_linear_system(MATRIX &, VECTOR &); \
  template void MatrixScaling::scale_system_solution(VECTOR &) const;

#ifdef DEAL_II_WITH_TRILINOS
InstantiateMatrixScaling(TrilinosWrappers::SparseMatrix,
                         TrilinosWrappers::MPI::Vector)
#endif

#ifdef DEAL_II_WITH_PETSC
  InstantiateMatrixScaling(PETScWrappers::MPI::SparseMatrix,
                           PETScWrappers::MPI::Vector)
#endif


    DEAL_II_NAMESPACE_CLOSE
