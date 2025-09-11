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


MatrixScaling::MatrixScaling(const ScalingControl &control)
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
          case MatrixScaling::ScalingControl::ScalingAlgorithm::sinkhorn_knopp:
            {
              sk_scaling(matrix,
                         control.sinkhorn_knopp_parameters.max_iterations);
            }
            break;
          case MatrixScaling::ScalingControl::ScalingAlgorithm::
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

  if constexpr (is_supported_matrix<Matrix>())
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
  else
    Assert(false, ExcNotImplemented());
}



template <class VectorType>
void
MatrixScaling::scale_system_solution(VectorType &sol) const
{
  AssertDimension(sol.size(), column_scaling.size());

  if constexpr (std::is_same_v<VectorType, dealii::Vector<double>> ||
                std::is_same_v<VectorType, dealii::Vector<float>>)
    sol.scale(column_scaling);
  else if constexpr (std::is_same_v<VectorType, dealii::BlockVector<double>> ||
                     std::is_same_v<VectorType, dealii::BlockVector<float>>)
    {
      for (unsigned int i = 0; i < sol.size(); ++i)
        sol[i] *= column_scaling[i];
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



template <class Matrix>
void
MatrixScaling::l1_scaling(Matrix &matrix, const unsigned int nsteps)
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

      // Check convergence. Subtract and add again. Alternative is copying every
      // time the whole vectors and subtracting -1
      row_norms.add(-1.0);
      col_norms.add(-1.0);
      if (row_norms.l1_norm() < control.scaling_tolerance &&
          col_norms.l1_norm() < control.scaling_tolerance)
        {
          converged = true;
          break;
        }
      row_norms.add(1.0);
      col_norms.add(1.0);

      for (unsigned int row = 0; row < matrix.m(); ++row)
        {
          for (auto it = matrix.begin(row); it != matrix.end(row); ++it)
            {
              // it->value() /=
              //   std::sqrt(row_norms[row] * col_norms[it->column()]);
              matrix.set(row,
                         it->column(),
                         it->value() /
                           std::sqrt(row_norms[row] * col_norms[it->column()]));
            }
          row_scaling[row] /= std::sqrt(row_norms[row]);
        }
      for (unsigned int col = 0; col < matrix.n(); ++col)
        {
          column_scaling[col] /= std::sqrt(col_norms[col]);
        }
    }
}



template <class Matrix>
void
MatrixScaling::linfty_scaling(Matrix &matrix, const unsigned int nsteps)
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
              row_norms[row] = std::max(row_norms[row], std::abs(it->value()));
              col_norms[it->column()] =
                std::max(col_norms[it->column()], std::abs(it->value()));
            }
        }

      // Check convergence. Sum and add again. Alternative is copying every time
      // the whole vectors and subtractin -1
      row_norms.add(-1.0);
      col_norms.add(-1.0);
      if (row_norms.linfty_norm() < control.scaling_tolerance &&
          col_norms.linfty_norm() < control.scaling_tolerance)
        {
          converged = true;
          break;
        }
      row_norms.add(1.0);
      col_norms.add(1.0);

      for (unsigned int row = 0; row < matrix.m(); ++row)
        {
          for (auto it = matrix.begin(row); it != matrix.end(row); ++it)
            {
              // it->value() /=
              //   std::sqrt(row_norms[row] * col_norms[it->column()]);
              matrix.set(row,
                         it->column(),
                         it->value() /
                           std::sqrt(row_norms[row] * col_norms[it->column()]));
            }
          row_scaling[row] /= std::sqrt(row_norms[row]);
        }
      for (unsigned int col = 0; col < matrix.n(); ++col)
        {
          column_scaling[col] /= std::sqrt(col_norms[col]);
        }
    }
}



template <class Matrix>
void
MatrixScaling::sk_scaling(Matrix &matrix, const unsigned int nsteps)
{
  using Number = typename Matrix::value_type;

  Vector<Number> row_norms(matrix.m());
  Vector<Number> col_norms(matrix.n());

  switch (control.sinkhorn_knopp_parameters.norm_type)
    {
      case MatrixScaling::ScalingControl::SKParameters::NormType::l1:
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
                  for (auto it = matrix.begin(row); it != matrix.end(row); ++it)
                    {
                      // it->value() /= row_norms[row];
                      matrix.set(row,
                                 it->column(),
                                 it->value() / row_norms[row]);
                      col_norms[it->column()] += std::abs(it->value());
                    }
                  row_scaling[row] /= row_norms[row];
                }
              // Here exit condition only on columns, rows are already
              // normalized
              col_norms.add(-1.0);
              if (col_norms.l1_norm() < control.scaling_tolerance)
                {
                  converged = true;
                  break;
                }
              col_norms.add(1.0);

              // Column step
              row_norms = 0;
              for (unsigned int row = 0; row < matrix.m(); ++row)
                {
                  for (auto it = matrix.begin(row); it != matrix.end(row); ++it)
                    {
                      // it->value() /= col_norms[it->column()];
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
              // Here exit condition only on rows, columns are already
              // normalized
              row_norms.add(-1.0);
              if (row_norms.l1_norm() < control.scaling_tolerance)
                {
                  converged = true;
                  break;
                }
              row_norms.add(1.0);
            }
        }
        break;
      case MatrixScaling::ScalingControl::SKParameters::NormType::l_infinity:
        {
          // Row_norms_0 to start the procedure
          row_norms = 0;
          for (unsigned int row = 0; row < matrix.m(); ++row)
            for (auto it = matrix.begin(row); it != matrix.end(row); ++it)
              row_norms[row] += std::max(row_norms[row], std::abs(it->value()));

          for (unsigned int i = 0; i < nsteps; i++)
            {
              // Row step
              col_norms = 0;
              for (unsigned int row = 0; row < matrix.m(); ++row)
                {
                  for (auto it = matrix.begin(row); it != matrix.end(row); ++it)
                    {
                      // it->value() /= row_norms[row];
                      matrix.set(row,
                                 it->column(),
                                 it->value() / row_norms[row]);
                      col_norms[it->column()] =
                        std::max(col_norms[it->column()],
                                 std::abs(it->value()));
                    }
                  row_scaling[row] /= row_norms[row];
                }
              // Here exit condition only on columns, rows are already
              // normalized
              col_norms.add(-1.0);
              if (col_norms.linfty_norm() < control.scaling_tolerance)
                {
                  converged = true;
                  break;
                }
              col_norms.add(1.0);
              // Column step
              row_norms = 0;
              for (unsigned int row = 0; row < matrix.m(); ++row)
                {
                  for (auto it = matrix.begin(row); it != matrix.end(row); ++it)
                    {
                      // it->value() /= col_norms[it->column()];
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
              // Here exit condition only on rows, columns are already
              // normalized
              row_norms.add(-1.0);
              if (row_norms.linfty_norm() < control.scaling_tolerance)
                {
                  converged = true;
                  break;
                }
              row_norms.add(1.0);
            }
        }
        break;
    }
}



DEAL_II_NAMESPACE_CLOSE
