// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test PreconditionChebyshev for different optimization levels. The
// test is similar to precondition_relaxation_01.


#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/matrix_tools.h>

#include "../tests.h"

template <typename VectorType>
class MyDiagonalMatrix
{
public:
  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    diagonal_matrix.vmult(dst, src);
  }

  void
  Tvmult(VectorType &dst, const VectorType &src) const
  {
    diagonal_matrix.vmult(dst, src);
  }

  VectorType &
  get_vector()
  {
    return diagonal_matrix.get_vector();
  }

private:
  DiagonalMatrix<VectorType> diagonal_matrix;
};

template <typename VectorType>
class MyDiagonalMatrixWithPreAndPost
{
public:
  void
  vmult(VectorType       &dst,
        const VectorType &src,
        const std::function<void(const unsigned int, const unsigned int)>
          &operation_before_matrix_vector_product = {},
        const std::function<void(const unsigned int, const unsigned int)>
          &operation_after_matrix_vector_product = {}) const
  {
    if (operation_before_matrix_vector_product)
      operation_before_matrix_vector_product(0, src.size());

    if (operation_before_matrix_vector_product)
      diagonal_matrix.vmult_add(dst, src);
    else
      diagonal_matrix.vmult(dst, src);

    if (operation_after_matrix_vector_product)
      operation_after_matrix_vector_product(0, src.size());
  }

  VectorType &
  get_vector()
  {
    return diagonal_matrix.get_vector();
  }

private:
  DiagonalMatrix<VectorType> diagonal_matrix;
};

template <typename SparseMatrixType>
class MySparseMatrix : public EnableObserverPointer
{
public:
  MySparseMatrix(const SparseMatrixType &sparse_matrix)
    : sparse_matrix(sparse_matrix)
  {}

  template <typename VectorType>
  void
  vmult(VectorType       &dst,
        const VectorType &src,
        const std::function<void(const unsigned int, const unsigned int)>
          &operation_before_matrix_vector_product = {},
        const std::function<void(const unsigned int, const unsigned int)>
          &operation_after_matrix_vector_product = {}) const
  {
    if (operation_before_matrix_vector_product)
      operation_before_matrix_vector_product(0, src.size());

    sparse_matrix.vmult(dst, src);

    if (operation_after_matrix_vector_product)
      operation_after_matrix_vector_product(0, src.size());
  }

  types::global_dof_index
  m() const
  {
    return sparse_matrix.m();
  }

  double
  el(unsigned int i, unsigned int j) const
  {
    return sparse_matrix.el(i, j);
  }

private:
  const SparseMatrixType &sparse_matrix;
};


template <typename PreconditionerType, typename VectorType>
std::tuple<double, double>
test(const PreconditionerType &preconditioner, const VectorType &src)
{
  VectorType dst;
  dst.reinit(src);

  dst = 1.0;
  preconditioner.vmult(dst, src);
  const double norm_0 = dst.l2_norm();

  dst = 1.0;
  preconditioner.step(dst, src);
  const double norm_1 = dst.l2_norm();

  return std::tuple<double, double>{norm_0, norm_1};
}



int
main()
{
  initlog();

  using Number     = double;
  using VectorType = Vector<Number>;
  using MatrixType = SparseMatrix<Number>;

  const unsigned int dim    = 2;
  const unsigned int degree = 1;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(3);

  FE_Q<dim>      fe(degree);
  QGauss<dim>    quad(degree + 1);
  MappingQ1<dim> mapping;

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  MatrixType system_matrix;

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);

  SparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from(dsp);
  system_matrix.reinit(sparsity_pattern);

  MatrixCreator::create_laplace_matrix(mapping,
                                       dof_handler,
                                       quad,
                                       system_matrix);

  VectorType diagonal(dof_handler.n_dofs());
  VectorType src(dof_handler.n_dofs());

  for (const auto &entry : system_matrix)
    if (entry.row() == entry.column())
      diagonal[entry.row()] = 1.0 / entry.value();

  src = 1.0;

  std::vector<unsigned int> n_iterationss{1, 2, 3};

  for (const auto n_iterations : n_iterationss)
    {
      std::vector<std::tuple<double, double>> results;

      {
        // 1) Test PreconditionChebyshev + DiagonalMatrix and matrix with
        // pre/post
        using PreconditionerType = DiagonalMatrix<VectorType>;

        using MyMatrixType = MySparseMatrix<MatrixType>;

        MyMatrixType my_system_matrix(system_matrix);

        PreconditionChebyshev<MyMatrixType, VectorType, PreconditionerType>
          preconditioner;

        PreconditionChebyshev<MyMatrixType, VectorType, PreconditionerType>::
          AdditionalData ad;
        ad.degree         = n_iterations;
        ad.preconditioner = std::make_shared<PreconditionerType>();
        ad.preconditioner->get_vector() = diagonal;

        preconditioner.initialize(my_system_matrix, ad);

        results.emplace_back(test(preconditioner, src));
      }

      {
        // 2) Test PreconditionChebyshev + DiagonalMatrix and matrix without
        // pre/post
        using PreconditionerType = DiagonalMatrix<VectorType>;

        PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>
          preconditioner;

        PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::
          AdditionalData ad;
        ad.degree         = n_iterations;
        ad.preconditioner = std::make_shared<PreconditionerType>();
        ad.preconditioner->get_vector() = diagonal;

        preconditioner.initialize(system_matrix, ad);

        results.emplace_back(test(preconditioner, src));
      }

      {
        // 3) Test PreconditionChebyshev + arbitrary preconditioner and matrix
        // without pre/post
        using PreconditionerType = MyDiagonalMatrix<VectorType>;

        PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>
          preconditioner;

        PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::
          AdditionalData ad;
        ad.degree         = n_iterations;
        ad.preconditioner = std::make_shared<PreconditionerType>();
        ad.preconditioner->get_vector() = diagonal;

        preconditioner.initialize(system_matrix, ad);

        results.emplace_back(test(preconditioner, src));
      }

      {
        // 4) Test PreconditionChebyshev + preconditioner with pre/post and
        // matrix without pre/post
        using PreconditionerType = MyDiagonalMatrixWithPreAndPost<VectorType>;

        PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>
          preconditioner;

        PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::
          AdditionalData ad;
        ad.degree         = n_iterations;
        ad.preconditioner = std::make_shared<PreconditionerType>();
        ad.preconditioner->get_vector() = diagonal;

        preconditioner.initialize(system_matrix, ad);

        results.emplace_back(test(preconditioner, src));
      }

      {
        // 5) Test PreconditionChebyshev + preconditioner with pre/post and
        // matrix with pre/post
        using PreconditionerType = MyDiagonalMatrixWithPreAndPost<VectorType>;

        using MyMatrixType = MySparseMatrix<MatrixType>;

        MyMatrixType my_system_matrix(system_matrix);

        PreconditionChebyshev<MyMatrixType, VectorType, PreconditionerType>
          preconditioner;

        PreconditionChebyshev<MyMatrixType, VectorType, PreconditionerType>::
          AdditionalData ad;
        ad.degree         = n_iterations;
        ad.preconditioner = std::make_shared<PreconditionerType>();
        ad.preconditioner->get_vector() = diagonal;

        preconditioner.initialize(my_system_matrix, ad);

        results.emplace_back(test(preconditioner, src));
      }

      if (std::all_of(results.begin(), results.end(), [&](const auto &a) {
            if (std::abs(std::get<0>(a) - std::get<0>(results[0])) > 1e-6)
              return false;

            if (std::abs(std::get<1>(a) - std::get<1>(results[0])) > 1e-6)
              return false;

            return true;
          }))
        deallog << "OK!" << std::endl;
      else
        deallog << "ERROR!" << std::endl;
    }
}
