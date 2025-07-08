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



// this tests the correctness of using the SolverCG class with interleaved
// vector operations, making use of the data locality options available by
// MatrixFree loops.

#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/solver_cg.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <iostream>

#include "../tests.h"

#include "matrix_vector_mf.h"



template <int dim, typename number = double>
class HelmholtzOperator : public EnableObserverPointer
{
public:
  using value_type = number;
  using VectorType = LinearAlgebra::distributed::Vector<number>;

  HelmholtzOperator(const MatrixFree<dim, number> &matrix_free)
    : n_calls_vmult(0)
    , data(matrix_free)
  {}

  void
  print_n_calls_special()
  {
    if (n_calls_vmult > 0)
      deallog << "Number of calls to special vmult: " << n_calls_vmult
              << std::endl;
  }

  void
  initialize(const Mapping<dim> &mapping, DoFHandler<dim> &dof_handler)
  {
    n_calls_vmult = 0;
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    const std::function<
      void(const MatrixFree<dim, typename VectorType::value_type> &,
           VectorType &,
           const VectorType &,
           const std::pair<unsigned int, unsigned int> &)>
      function = helmholtz_operator<dim, -1, VectorType, 0>;
    data.cell_loop(function, dst, src, true);
    for (const unsigned int i : data.get_constrained_dofs())
      dst.local_element(i) = src.local_element(i);
  }

  void
  vmult(VectorType       &dst,
        const VectorType &src,
        const std::function<void(const unsigned int, const unsigned int)>
          &operation_before_loop,
        const std::function<void(const unsigned int, const unsigned int)>
          &operation_after_loop) const
  {
    ++n_calls_vmult;
    const std::function<
      void(const MatrixFree<dim, typename VectorType::value_type> &,
           VectorType &,
           const VectorType &,
           const std::pair<unsigned int, unsigned int> &)>
      function = helmholtz_operator<dim, -1, VectorType, 0>;
    data.cell_loop(
      function, dst, src, operation_before_loop, operation_after_loop);
    for (const unsigned int i : data.get_constrained_dofs())
      dst.local_element(i) = src.local_element(i);
  }

private:
  mutable unsigned int    n_calls_vmult;
  MatrixFree<dim, number> data;
};


// Preconditioner class without any apply or apply_to_subrange function, as
// opposed to deal.II's DiagonalMatrix
template <typename Number>
struct MyDiagonalMatrix
{
  MyDiagonalMatrix(const LinearAlgebra::distributed::Vector<Number> &vec)
    : vec(vec)
  {}

  void
  vmult(LinearAlgebra::distributed::Vector<Number>       &dst,
        const LinearAlgebra::distributed::Vector<Number> &src) const
  {
    dst = src;
    dst.scale(vec);
  }

  const LinearAlgebra::distributed::Vector<Number> &vec;
};



// Preconditioner class only proving an apply_to_subrange function,
template <typename Number>
struct DiagonalMatrixSubrange
{
  DiagonalMatrixSubrange(const LinearAlgebra::distributed::Vector<Number> &vec)
    : vec(vec)
  {}

  void
  vmult(LinearAlgebra::distributed::Vector<Number>       &dst,
        const LinearAlgebra::distributed::Vector<Number> &src) const
  {
    dst = src;
    dst.scale(vec);
  }

  void
  apply_to_subrange(const unsigned int begin_range,
                    const unsigned int end_range,
                    const Number      *src_pointer_to_current_range,
                    Number            *dst_pointer_to_current_range) const
  {
    AssertIndexRange(begin_range,
                     vec.locally_owned_elements().n_elements() + 1);
    AssertIndexRange(end_range, vec.locally_owned_elements().n_elements() + 1);

    const Number      *diagonal_entry = vec.begin() + begin_range;
    const unsigned int length         = end_range - begin_range;

    DEAL_II_OPENMP_SIMD_PRAGMA
    for (unsigned int i = 0; i < length; ++i)
      dst_pointer_to_current_range[i] =
        diagonal_entry[i] * src_pointer_to_current_range[i];
  }

  const LinearAlgebra::distributed::Vector<Number> &vec;
};



template <int dim, typename number>
void
test(const unsigned int fe_degree)
{
  using VectorType = LinearAlgebra::distributed::Vector<number>;

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(6 - dim);

  MappingQ1<dim>  mapping;
  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);

  typename MatrixFree<dim, number>::AdditionalData addit_data;
  addit_data.tasks_parallel_scheme =
    MatrixFree<dim, number>::AdditionalData::none;

  IndexSet                  relevant_set;
  AffineConstraints<double> constraints;
  {
    constraints.clear();
    relevant_set = DoFTools::extract_locally_relevant_dofs(dof);
    constraints.reinit(dof.locally_owned_dofs(), relevant_set);
    constraints.close();
  }

  DoFRenumbering::matrix_free_data_locality(dof, constraints, addit_data);

  {
    constraints.clear();
    relevant_set = DoFTools::extract_locally_relevant_dofs(dof);
    constraints.reinit(dof.locally_owned_dofs(), relevant_set);
    constraints.close();
  }

  const QGauss<1>         quadrature(dof.get_fe().degree + 1);
  MatrixFree<dim, number> mf_data;
  mf_data.reinit(mapping, dof, constraints, quadrature, addit_data);

  MatrixFreeTest<dim, -1, number, VectorType> mf(mf_data);
  VectorType                                  rhs, sol;
  mf_data.initialize_dof_vector(rhs);
  mf_data.initialize_dof_vector(sol);
  rhs = 1. / std::sqrt(rhs.size());
  DiagonalMatrix<VectorType> preconditioner;
  mf_data.initialize_dof_vector(preconditioner.get_vector());
  mf.vmult(preconditioner.get_vector(), rhs);
  for (number &a : preconditioner.get_vector())
    if (a != 0.)
      a = 1. / a;
    else
      a = 0.;

  // Step 1: solve with CG solver for a matrix that does not support the
  // interleaving operation
  {
    deallog << "CG solver without interleaving support" << std::endl;
    SolverControl        control(200, 1e-2 * rhs.l2_norm());
    SolverCG<VectorType> solver(control);
    solver.solve(mf, sol, rhs, preconditioner);
    deallog << "Norm of the solution: " << sol.l2_norm() << std::endl;
  }

  // Step 2: solve with CG solver for a matrix that does support the
  // interleaving operation
  {
    deallog << "CG solver with interleaving support" << std::endl;
    sol = 0;
    HelmholtzOperator<dim, number> matrix(mf_data);
    SolverControl                  control(200, 1e-2 * rhs.l2_norm());
    SolverCG<VectorType>           solver(control);
    solver.solve(matrix, sol, rhs, preconditioner);
    deallog << "Norm of the solution: " << sol.l2_norm() << std::endl;
    matrix.print_n_calls_special();
  }

  // Step 3: solve with CG solver for a matrix that does support but a
  // preconditioner that does not support the interleaving operation
  {
    deallog << "CG solver with matrix with interleaving support but "
            << "no preconditioner" << std::endl;
    sol = 0;
    HelmholtzOperator<dim, number> matrix(mf_data);
    MyDiagonalMatrix<number>       simple_diagonal(preconditioner.get_vector());
    SolverControl                  control(200, 1e-2 * rhs.l2_norm());
    SolverCG<VectorType>           solver(control);
    solver.solve(matrix, sol, rhs, simple_diagonal);
    deallog << "Norm of the solution: " << sol.l2_norm() << std::endl;
    matrix.print_n_calls_special();
  }

  // Step 4: solve with CG solver for a matrix that does support the
  // interleaving operation and a preconditioner that only supports the
  // apply_to_subrange function
  {
    deallog << "CG solver with interleaving support and "
            << "preconditioner working on subrange" << std::endl;
    sol = 0;
    HelmholtzOperator<dim, number> matrix(mf_data);
    DiagonalMatrixSubrange<number> diagonal_subrange(
      preconditioner.get_vector());
    SolverControl        control(200, 1e-2 * rhs.l2_norm());
    SolverCG<VectorType> solver(control);
    solver.solve(matrix, sol, rhs, diagonal_subrange);
    deallog << "Norm of the solution: " << sol.l2_norm() << std::endl;
    matrix.print_n_calls_special();
  }
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  mpi_initlog();
  deallog << std::setprecision(8);

  test<2, double>(3);
  test<3, double>(4);
  test<3, double>(3);
}
