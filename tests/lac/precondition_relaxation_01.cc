// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2022 by the deal.II authors
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


// Test PreconditionRelaxation for different Jacobi preconditioners.


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



template <typename PreconditionerType, typename VectorType>
std::tuple<double, double, double, double>
test(const PreconditionerType &preconditioner, const VectorType &src)
{
  VectorType dst;
  dst.reinit(src);

  preconditioner.vmult(dst, src);
  const double norm_0 = dst.l2_norm();

  preconditioner.step(dst, src);
  const double norm_1 = dst.l2_norm();

  preconditioner.Tvmult(dst, src);
  const double norm_2 = dst.l2_norm();

  preconditioner.Tstep(dst, src);
  const double norm_3 = dst.l2_norm();

  return std::tuple<double, double, double, double>{norm_0,
                                                    norm_1,
                                                    norm_2,
                                                    norm_3};
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

  std::vector<double>       relaxations{1.0, 0.9, 1.1};
  std::vector<unsigned int> n_iterationss{1, 2, 3};

  for (const auto relaxation : relaxations)
    for (const auto n_iterations : n_iterationss)
      {
        std::vector<std::tuple<double, double, double, double>> results;

        {
          // Test PreconditionJacobi
          PreconditionJacobi<MatrixType> preconditioner;

          PreconditionJacobi<MatrixType>::AdditionalData ad;
          ad.relaxation   = relaxation;
          ad.n_iterations = n_iterations;

          preconditioner.initialize(system_matrix, ad);

          results.emplace_back(test(preconditioner, src));
        }

        {
          // Test PreconditionRelaxation + DiagonalMatrix: optimized path
          using PreconditionerType = DiagonalMatrix<VectorType>;

          PreconditionRelaxation<MatrixType, PreconditionerType> preconditioner;

          PreconditionRelaxation<MatrixType, PreconditionerType>::AdditionalData
            ad;
          ad.relaxation     = relaxation;
          ad.n_iterations   = n_iterations;
          ad.preconditioner = std::make_shared<PreconditionerType>();
          ad.preconditioner->get_vector() = diagonal;

          preconditioner.initialize(system_matrix, ad);

          results.emplace_back(test(preconditioner, src));
        }

        {
          // Test PreconditionRelaxation + wrapper around DiagonalMatrix:
          // optimized path cannot be used
          using PreconditionerType = MyDiagonalMatrix<VectorType>;

          PreconditionRelaxation<MatrixType, PreconditionerType> preconditioner;

          PreconditionRelaxation<MatrixType, PreconditionerType>::AdditionalData
            ad;
          ad.relaxation     = relaxation;
          ad.n_iterations   = n_iterations;
          ad.preconditioner = std::make_shared<PreconditionerType>();
          ad.preconditioner->get_vector() = diagonal;

          preconditioner.initialize(system_matrix, ad);

          results.emplace_back(test(preconditioner, src));
        }

        if (std::equal(results.begin(),
                       results.end(),
                       results.begin(),
                       [](const auto &a, const auto &b) {
                         if (std::abs(std::get<0>(a) - std::get<0>(b)) > 1e-6)
                           return false;

                         if (std::abs(std::get<1>(a) - std::get<1>(b)) > 1e-6)
                           return false;

                         if (std::abs(std::get<2>(a) - std::get<2>(b)) > 1e-6)
                           return false;

                         if (std::abs(std::get<3>(a) - std::get<3>(b)) > 1e-6)
                           return false;

                         return true;
                       }))
          deallog << "OK!" << std::endl;
        else
          deallog << "ERROR!" << std::endl;
      }
}
