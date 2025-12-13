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

#include <boost/serialization/utility.hpp>

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

// Check if the matrix type is sequential
template <typename Matrix>
constexpr bool
is_sequential_matrix()
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



MatrixScaling::AdditionalData::SKParameters::SKParameters(
  const NormType     norm_type,
  const unsigned int max_iterations)
  : max_iterations(max_iterations)
  , norm_type(norm_type)
{}



MatrixScaling::AdditionalData::l1linfParameters::l1linfParameters(
  const unsigned int start_inf_norm_steps,
  const unsigned int l1_norm_steps,
  const unsigned int end_inf_norm_steps)
  : start_inf_norm_steps(start_inf_norm_steps)
  , l1_norm_steps(l1_norm_steps)
  , end_inf_norm_steps(end_inf_norm_steps)
{}



MatrixScaling::AdditionalData::AdditionalData(
  const double           scaling_tolerance,
  const ScalingAlgorithm alg,
  const SKParameters     sk_params,
  const l1linfParameters l1linf_params)
  : scaling_tolerance(scaling_tolerance)
  , algorithm(alg)
  , sinkhorn_knopp_parameters(sk_params)
  , l1linf_parameters(l1linf_params)
{}



MatrixScaling::MatrixScaling(const AdditionalData &control)
  : control(control)
  , row_scaling()
  , column_scaling()
{}



template <class Matrix>
bool
MatrixScaling::find_scaling_and_scale_matrix(Matrix &matrix)
{
  bool converged = false;

  if constexpr (is_sequential_matrix<Matrix>())
    {
      const auto n_rows = matrix.m();
      const auto n_cols = matrix.n();

      row_scaling.reinit(n_rows);
      column_scaling.reinit(n_cols);

      row_scaling    = 1.0;
      column_scaling = 1.0;

      switch (control.algorithm)
        {
          case MatrixScaling::AdditionalData::ScalingAlgorithm::sinkhorn_knopp:
            {
              converged =
                do_sk_scaling(matrix,
                              control.sinkhorn_knopp_parameters.max_iterations);

              break;
            }

          case MatrixScaling::AdditionalData::ScalingAlgorithm::
            l1_linf_symmetry_preserving:
            {
              if (control.l1linf_parameters.start_inf_norm_steps > 0 &&
                  !converged)
                converged = do_linfty_scaling(
                  matrix, control.l1linf_parameters.start_inf_norm_steps);
              if (control.l1linf_parameters.l1_norm_steps > 0 && !converged)
                converged =
                  do_l1_scaling(matrix,
                                control.l1linf_parameters.l1_norm_steps);
              if (control.l1linf_parameters.end_inf_norm_steps > 0 &&
                  !converged)
                converged = do_linfty_scaling(
                  matrix, control.l1linf_parameters.end_inf_norm_steps);

              break;
            }

          default:
            DEAL_II_ASSERT_UNREACHABLE();
        }
    }
  else if constexpr (std::is_same_v<Matrix, TrilinosWrappers::SparseMatrix> ||
                     std::is_same_v<Matrix, PETScWrappers::MPI::SparseMatrix>)
    {
      locally_owned_rows = matrix.locally_owned_range_indices();

      // If the matrix is square then use the same partitioning for
      // columns and rows to easily scale a linear system, otherwise create a
      // new balanced partitioning for the columns. We do this because
      // distributed matrices do not really distribute columns, but only rows:
      // we want to avoid saving the whole column scaling on each MPI rank
      if (matrix.n() == matrix.m())
        locally_owned_cols = locally_owned_rows;
      else
        locally_owned_cols =
          Utilities::MPI::create_evenly_distributed_partitioning(
            matrix.get_mpi_communicator(), matrix.n());

      ghost_columns.clear();
      ghost_columns.set_size(matrix.n());

      row_scaling.reinit(locally_owned_rows.n_elements());
      column_scaling.reinit(locally_owned_cols.n_elements());

      row_scaling    = 1.0;
      column_scaling = 1.0;

      // Identify ghost columns
      for (const auto row : locally_owned_rows)
        {
          for (auto it = matrix.begin(row); it != matrix.end(row); ++it)
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
              converged =
                do_sk_scaling(matrix,
                              control.sinkhorn_knopp_parameters.max_iterations);

              break;
            }

          case MatrixScaling::AdditionalData::ScalingAlgorithm::
            l1_linf_symmetry_preserving:
            {
              if (control.l1linf_parameters.start_inf_norm_steps > 0 &&
                  !converged)
                converged = do_linfty_scaling(
                  matrix, control.l1linf_parameters.start_inf_norm_steps);
              if (control.l1linf_parameters.l1_norm_steps > 0 && !converged)
                converged =
                  do_l1_scaling(matrix,
                                control.l1linf_parameters.l1_norm_steps);
              if (control.l1linf_parameters.end_inf_norm_steps > 0 &&
                  !converged)
                converged = do_linfty_scaling(
                  matrix, control.l1linf_parameters.end_inf_norm_steps);

              break;
            }

          default:
            DEAL_II_ASSERT_UNREACHABLE();
        }
    }
  else
    {
      (void)matrix;
      Assert(false, ExcNotImplemented());
    }

  return converged;
}



template <class Matrix, class VectorType>
bool
MatrixScaling::find_scaling_and_scale_linear_system(Matrix     &matrix,
                                                    VectorType &rhs)
{
  AssertDimension(matrix.m(), rhs.size());
  AssertDimension(matrix.m(), matrix.n());

  bool converged = false;

  converged = find_scaling_and_scale_matrix(matrix);

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

  return converged;
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



template <typename Number>
bool
MatrixScaling::check_convergence(
  const Vector<Number>                     &row_col_norm,
  const MatrixScaling::ConvergenceNormType &norm_type) const
{
  Number convergence_norm = 0;
  switch (norm_type)
    {
      case ConvergenceNormType::l1:
        {
          for (const auto &val : row_col_norm)
            convergence_norm += std::abs(val - Number(1.0));

          return convergence_norm < control.scaling_tolerance;
        }

      case ConvergenceNormType::l_infty:
        {
          for (const auto &val : row_col_norm)
            convergence_norm =
              std::max(convergence_norm, std::abs(val - Number(1.0)));

          return convergence_norm < control.scaling_tolerance;
        }

      default:
        {
          DEAL_II_ASSERT_UNREACHABLE();
          return false;
        }
    }
}



template <typename Number>
bool
MatrixScaling::check_convergence(
  const Vector<Number>                     &row_norm,
  const Vector<Number>                     &col_norm,
  const MatrixScaling::ConvergenceNormType &norm_type) const
{
  Number convergence_row_norm = 0;
  Number convergence_col_norm = 0;
  switch (norm_type)
    {
      case ConvergenceNormType::l1:
        {
          for (const auto &val : row_norm)
            convergence_row_norm += std::abs(val - Number(1.0));
          for (const auto &val : col_norm)
            convergence_col_norm += std::abs(val - Number(1.0));

          return (convergence_row_norm < control.scaling_tolerance &&
                  convergence_col_norm < control.scaling_tolerance);
        }

      case ConvergenceNormType::l_infty:
        {
          for (const auto &val : row_norm)
            convergence_row_norm =
              std::max(convergence_row_norm, std::abs(val - Number(1.0)));
          for (const auto &val : col_norm)
            convergence_col_norm =
              std::max(convergence_col_norm, std::abs(val - Number(1.0)));

          return (convergence_row_norm < control.scaling_tolerance &&
                  convergence_col_norm < control.scaling_tolerance);
        }

      default:
        {
          DEAL_II_ASSERT_UNREACHABLE();
          return false;
        }
    }
}



bool
MatrixScaling::check_convergence(
  const Vector<double>                     &local_row_col_norm,
  const MatrixScaling::ConvergenceNormType &norm_type,
  const MPI_Comm                            mpi_communicator) const
{
  double convergence_norm = 0;
  bool   local_not_converged;

  switch (norm_type)
    {
      case ConvergenceNormType::l1:
        {
          for (const auto &val : local_row_col_norm)
            convergence_norm += std::abs(val - 1.0);

          local_not_converged = !(convergence_norm < control.scaling_tolerance);
          break;
        }

      case ConvergenceNormType::l_infty:
        {
          for (const auto &val : local_row_col_norm)
            convergence_norm = std::max(convergence_norm, std::abs(val - 1.0));

          local_not_converged = !(convergence_norm < control.scaling_tolerance);
          break;
        }

      default:
        {
          DEAL_II_ASSERT_UNREACHABLE();
          return false;
        }
    }

  const bool any_not_converged =
    Utilities::MPI::logical_or(local_not_converged, mpi_communicator);

  return !any_not_converged;
}



bool
MatrixScaling::check_convergence(
  const Vector<double>                     &local_row_norm,
  const Vector<double>                     &local_col_norm,
  const MatrixScaling::ConvergenceNormType &norm_type,
  const MPI_Comm                            mpi_communicator) const
{
  double convergence_row_norm = 0;
  double convergence_col_norm = 0;
  bool   local_not_converged;

  switch (norm_type)
    {
      case ConvergenceNormType::l1:
        {
          for (const auto &val : local_row_norm)
            convergence_row_norm += std::abs(val - 1.0);
          for (const auto &val : local_col_norm)
            convergence_col_norm += std::abs(val - 1.0);

          local_not_converged =
            !(convergence_row_norm < control.scaling_tolerance &&
              convergence_col_norm < control.scaling_tolerance);
          break;
        }

      case ConvergenceNormType::l_infty:
        {
          for (const auto &val : local_row_norm)
            convergence_row_norm =
              std::max(convergence_row_norm, std::abs(val - 1.0));
          for (const auto &val : local_col_norm)
            convergence_col_norm =
              std::max(convergence_col_norm, std::abs(val - 1.0));

          local_not_converged =
            !(convergence_row_norm < control.scaling_tolerance &&
              convergence_col_norm < control.scaling_tolerance);
          break;
        }

      default:
        {
          DEAL_II_ASSERT_UNREACHABLE();
          return false;
        }
    }

  const bool any_not_converged =
    Utilities::MPI::logical_or(local_not_converged, mpi_communicator);

  return !any_not_converged;
}



void
MatrixScaling::send_prepare_col_norms(
  const std::map<types::global_dof_index, double> &partial_column_norms,
  std::map<unsigned int,
           std::vector<std::pair<types::global_dof_index, double>>> &send_data,
  Vector<double> &local_col_norms)
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
  const std::map<unsigned int,
                 std::vector<std::pair<types::global_dof_index, double>>>
                       &received_data,
  const Vector<double> &local_col_norms,
  std::map<unsigned int,
           std::vector<std::pair<types::global_dof_index, double>>>
    &send_column_norms)
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



template <class Matrix>
bool
MatrixScaling::do_l1_scaling(Matrix &matrix, const unsigned int nsteps)
{
  if constexpr (is_sequential_matrix<Matrix>())
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

          if (check_convergence(row_norms,
                                col_norms,
                                MatrixScaling::ConvergenceNormType::l1))
            return true;

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

          for (const auto row : locally_owned_rows)
            {
              auto local_row_idx = partitioner.global_to_local(row);
              for (auto it = matrix.begin(row); it != matrix.end(row); ++it)
                {
                  local_row_norms[local_row_idx] += std::abs(it->value());
                  partial_column_norms[it->column()] += std::abs(it->value());
                }
            }
          // Communicate partial column norms
          std::map<unsigned int,
                   std::vector<std::pair<types::global_dof_index, double>>>
            send_data;

          send_prepare_col_norms(partial_column_norms,
                                 send_data,
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
                                MatrixScaling::ConvergenceNormType::l1,
                                matrix.get_mpi_communicator()))
            return true;


          for (unsigned int i = 0; i < local_row_norms.size(); ++i)
            row_scaling[i] /= std::sqrt(local_row_norms[i]);
          for (unsigned int i = 0; i < local_col_norms.size(); ++i)
            column_scaling[i] /= std::sqrt(local_col_norms[i]);

          // Communicate column norm values to all ranks that need them
          std::map<unsigned int,
                   std::vector<std::pair<types::global_dof_index, double>>>
            send_column_norms;

          send_prepare_updated_col_norms(received_data,
                                         local_col_norms,
                                         send_column_norms);

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
          for (const auto row : locally_owned_rows)
            {
              auto local_row_idx = partitioner.global_to_local(row);
              auto row_size      = matrix.row_length(row);
              std::vector<double>                  values(row_size);
              std::vector<types::global_dof_index> columns(row_size);
              unsigned int                         idx = 0;
              for (auto it = matrix.begin(row); it != matrix.end(row); ++it)
                {
                  if (locally_owned_cols.is_element(it->column()))
                    {
                      unsigned int local_col_idx =
                        partitioner.global_to_local(it->column());
                      columns[idx] = it->column();
                      values[idx] =
                        it->value() / std::sqrt(local_col_norms[local_col_idx] *
                                                local_row_norms[local_row_idx]);
                      ++idx;
                    }
                  else
                    {
                      double col_norms =
                        ghost_column_norms_lookup[it->column()];
                      columns[idx] = it->column();
                      values[idx] =
                        it->value() /
                        std::sqrt(col_norms * local_row_norms[local_row_idx]);
                      ++idx;
                    }
                }
              matrix.set(row, columns, values);
              if constexpr (std::is_same_v<Matrix,
                                           PETScWrappers::MPI::SparseMatrix>)
                matrix.compress(VectorOperation::insert);
            }
          if constexpr (std::is_same_v<Matrix, TrilinosWrappers::SparseMatrix>)
            matrix.compress(VectorOperation::insert);
        }
    }
  return false;
}



template <class Matrix>
bool
MatrixScaling::do_linfty_scaling(Matrix &matrix, const unsigned int nsteps)
{
  if constexpr (is_sequential_matrix<Matrix>())
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

          if (check_convergence(row_norms,
                                col_norms,
                                MatrixScaling::ConvergenceNormType::l_infty))
            return true;


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

          for (const auto row : locally_owned_rows)
            {
              auto local_row_idx = partitioner.global_to_local(row);
              for (auto it = matrix.begin(row); it != matrix.end(row); ++it)
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

          send_prepare_col_norms(partial_column_norms,
                                 send_data,
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
                                MatrixScaling::ConvergenceNormType::l_infty,
                                matrix.get_mpi_communicator()))
            return true;

          for (unsigned int i = 0; i < local_row_norms.size(); ++i)
            row_scaling[i] /= std::sqrt(local_row_norms[i]);
          for (unsigned int i = 0; i < local_col_norms.size(); ++i)
            column_scaling[i] /= std::sqrt(local_col_norms[i]);

          // Communicate column norm values to all ranks that need them
          std::map<unsigned int,
                   std::vector<std::pair<types::global_dof_index, double>>>
            send_column_norms;

          send_prepare_updated_col_norms(received_data,
                                         local_col_norms,
                                         send_column_norms);

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
          for (const auto row : locally_owned_rows)
            {
              auto local_row_idx = partitioner.global_to_local(row);
              auto row_size      = matrix.row_length(row);
              std::vector<double>                  values(row_size);
              std::vector<types::global_dof_index> columns(row_size);
              unsigned int                         idx = 0;
              for (auto it = matrix.begin(row); it != matrix.end(row); ++it)
                {
                  if (locally_owned_cols.is_element(it->column()))
                    {
                      unsigned int local_col_idx =
                        partitioner.global_to_local(it->column());
                      columns[idx] = it->column();
                      values[idx] =
                        it->value() / std::sqrt(local_col_norms[local_col_idx] *
                                                local_row_norms[local_row_idx]);
                      ++idx;
                    }
                  else
                    {
                      double col_norms =
                        ghost_column_norms_lookup[it->column()];
                      columns[idx] = it->column();
                      values[idx] =
                        it->value() /
                        std::sqrt(col_norms * local_row_norms[local_row_idx]);
                      ++idx;
                    }
                }
              matrix.set(row, columns, values);
              if constexpr (std::is_same_v<Matrix,
                                           PETScWrappers::MPI::SparseMatrix>)
                matrix.compress(VectorOperation::insert);
            }
          if constexpr (std::is_same_v<Matrix, TrilinosWrappers::SparseMatrix>)
            matrix.compress(VectorOperation::insert);
        }
    }
  return false;
}



template <class Matrix>
bool
MatrixScaling::do_sk_scaling(Matrix &matrix, const unsigned int nsteps)
{
  if constexpr (is_sequential_matrix<Matrix>())
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

                  if (check_convergence(col_norms,
                                        MatrixScaling::ConvergenceNormType::l1))
                    return true;

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

                  if (check_convergence(row_norms,
                                        MatrixScaling::ConvergenceNormType::l1))
                    return true;
                }
              break;
            }

          case MatrixScaling::AdditionalData::SKParameters::NormType::l_infty:
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

                  if (check_convergence(
                        col_norms, MatrixScaling::ConvergenceNormType::l_infty))
                    return true;

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

                  if (check_convergence(
                        row_norms, MatrixScaling::ConvergenceNormType::l_infty))
                    return true;
                }
              break;
            }

          default:
            DEAL_II_ASSERT_UNREACHABLE();
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
              for (const auto row : locally_owned_rows)
                {
                  auto local_row_idx = partitioner.global_to_local(row);
                  for (auto it = matrix.begin(row); it != matrix.end(row); ++it)
                    local_row_norms[local_row_idx] += std::abs(it->value());
                }

              for (unsigned int i = 0; i < nsteps; i++)
                {
                  // Row step
                  local_col_norms = 0;
                  partial_column_norms.clear();

                  for (const auto row : locally_owned_rows)
                    {
                      auto local_row_idx = partitioner.global_to_local(row);
                      auto row_size      = matrix.row_length(row);
                      std::vector<double>                  values(row_size);
                      std::vector<types::global_dof_index> columns(row_size);
                      unsigned int                         idx = 0;
                      for (auto it = matrix.begin(row); it != matrix.end(row);
                           ++it)
                        {
                          partial_column_norms[it->column()] += std::abs(
                            it->value() / local_row_norms[local_row_idx]);
                          columns[idx] = it->column();
                          values[idx] =
                            it->value() / local_row_norms[local_row_idx];
                          ++idx;
                        }
                      matrix.set(row, columns, values);
                      if constexpr (std::is_same_v<
                                      Matrix,
                                      PETScWrappers::MPI::SparseMatrix>)
                        matrix.compress(VectorOperation::insert);

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

                  send_prepare_col_norms(partial_column_norms,
                                         send_data,
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
                                        MatrixScaling::ConvergenceNormType::l1,
                                        matrix.get_mpi_communicator()))
                    return true;

                  // Column step
                  local_row_norms = 0;

                  // Communicate column norm values to all ranks that need them
                  std::map<
                    unsigned int,
                    std::vector<std::pair<types::global_dof_index, double>>>
                    send_column_norms;

                  send_prepare_updated_col_norms(received_data,
                                                 local_col_norms,
                                                 send_column_norms);

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

                  for (const auto row : locally_owned_rows)
                    {
                      auto local_row_idx = partitioner.global_to_local(row);
                      auto row_size      = matrix.row_length(row);
                      std::vector<double>                  values(row_size);
                      std::vector<types::global_dof_index> columns(row_size);
                      unsigned int                         idx = 0;
                      for (auto it = matrix.begin(row); it != matrix.end(row);
                           ++it)
                        {
                          if (locally_owned_cols.is_element(it->column()))
                            {
                              unsigned int local_col_idx =
                                partitioner.global_to_local(it->column());

                              local_row_norms[local_row_idx] += std::abs(
                                it->value() / local_col_norms[local_col_idx]);

                              columns[idx] = it->column();
                              values[idx] =
                                it->value() / local_col_norms[local_col_idx];
                              ++idx;
                            }
                          else
                            {
                              double col_norms =
                                ghost_column_norms_lookup[it->column()];

                              local_row_norms[local_row_idx] +=
                                std::abs(it->value() / col_norms);

                              columns[idx] = it->column();
                              values[idx]  = it->value() / col_norms;
                              ++idx;
                            }
                        }
                      matrix.set(row, columns, values);
                      if constexpr (std::is_same_v<
                                      Matrix,
                                      PETScWrappers::MPI::SparseMatrix>)
                        matrix.compress(VectorOperation::insert);
                    }
                  if constexpr (std::is_same_v<Matrix,
                                               TrilinosWrappers::SparseMatrix>)
                    matrix.compress(VectorOperation::insert);
                  for (unsigned int i = 0; i < local_col_norms.size(); ++i)
                    column_scaling[i] /= local_col_norms[i];

                  if (check_convergence(local_row_norms,
                                        MatrixScaling::ConvergenceNormType::l1,
                                        matrix.get_mpi_communicator()))
                    return true;
                }
              break;
            }

          case MatrixScaling::AdditionalData::SKParameters::NormType::l_infty:
            {
              // Row_norms_0 to start the procedure
              local_row_norms = 0;
              for (const auto row : locally_owned_rows)
                {
                  auto local_row_idx = partitioner.global_to_local(row);
                  for (auto it = matrix.begin(row); it != matrix.end(row); ++it)
                    local_row_norms[local_row_idx] =
                      std::max(local_row_norms[local_row_idx],
                               std::abs(it->value()));
                }

              for (unsigned int i = 0; i < nsteps; i++)
                {
                  // Row step
                  local_col_norms = 0;
                  partial_column_norms.clear();

                  for (const auto row : locally_owned_rows)
                    {
                      auto local_row_idx = partitioner.global_to_local(row);
                      auto row_size      = matrix.row_length(row);
                      std::vector<double>                  values(row_size);
                      std::vector<types::global_dof_index> columns(row_size);
                      unsigned int                         idx = 0;
                      for (auto it = matrix.begin(row); it != matrix.end(row);
                           ++it)
                        {
                          partial_column_norms[it->column()] =
                            std::max(partial_column_norms[it->column()],
                                     std::abs(it->value() /
                                              local_row_norms[local_row_idx]));
                          columns[idx] = it->column();
                          values[idx] =
                            it->value() / local_row_norms[local_row_idx];
                          ++idx;
                        }
                      matrix.set(row, columns, values);
                      if constexpr (std::is_same_v<
                                      Matrix,
                                      PETScWrappers::MPI::SparseMatrix>)
                        matrix.compress(VectorOperation::insert);

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

                  send_prepare_col_norms(partial_column_norms,
                                         send_data,
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

                  if (check_convergence(
                        local_col_norms,
                        MatrixScaling::ConvergenceNormType::l_infty,
                        matrix.get_mpi_communicator()))
                    return true;

                  // Column step
                  local_row_norms = 0;

                  // Communicate column norm values to all ranks that need them
                  std::map<
                    unsigned int,
                    std::vector<std::pair<types::global_dof_index, double>>>
                    send_column_norms;

                  send_prepare_updated_col_norms(received_data,
                                                 local_col_norms,
                                                 send_column_norms);

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

                  for (const auto row : locally_owned_rows)
                    {
                      auto local_row_idx = partitioner.global_to_local(row);
                      auto row_size      = matrix.row_length(row);
                      std::vector<double>                  values(row_size);
                      std::vector<types::global_dof_index> columns(row_size);
                      unsigned int                         idx = 0;
                      for (auto it = matrix.begin(row); it != matrix.end(row);
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

                              columns[idx] = it->column();
                              values[idx] =
                                it->value() / local_col_norms[local_col_idx];
                              ++idx;
                            }
                          else
                            {
                              double col_norms =
                                ghost_column_norms_lookup[it->column()];

                              local_row_norms[local_row_idx] =
                                std::max(local_row_norms[local_row_idx],
                                         std::abs(it->value() / col_norms));

                              columns[idx] = it->column();
                              values[idx]  = it->value() / col_norms;
                              ++idx;
                            }
                        }
                      matrix.set(row, columns, values);
                      if constexpr (std::is_same_v<
                                      Matrix,
                                      PETScWrappers::MPI::SparseMatrix>)
                        matrix.compress(VectorOperation::insert);
                    }
                  if constexpr (std::is_same_v<Matrix,
                                               TrilinosWrappers::SparseMatrix>)
                    matrix.compress(VectorOperation::insert);
                  for (unsigned int i = 0; i < local_col_norms.size(); ++i)
                    column_scaling[i] /= local_col_norms[i];

                  if (check_convergence(
                        local_row_norms,
                        MatrixScaling::ConvergenceNormType::l_infty,
                        matrix.get_mpi_communicator()))
                    return true;
                }
              break;
            }

          default:
            DEAL_II_ASSERT_UNREACHABLE();
        }
    }
  return false;
}



#define InstantiateMatrixScaling(MATRIX, VECTOR)                               \
  template bool MatrixScaling::find_scaling_and_scale_matrix(MATRIX &);        \
  template bool MatrixScaling::find_scaling_and_scale_linear_system(MATRIX &,  \
                                                                    VECTOR &); \
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
