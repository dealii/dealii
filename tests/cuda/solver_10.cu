// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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

// Check that dealii::SolverRelaxation works with CUDAWrappers::SparseMatrix

#include <deal.II/base/cuda.h>
#include <deal.II/base/exceptions.h>

#include <deal.II/lac/cuda_sparse_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_relaxation.h>
#include <deal.II/lac/vector.h>

#include "../testmatrix.h"
#include "../tests.h"


template <typename MatrixType>
class RelaxationOperator
{
public:
  RelaxationOperator(const MatrixType &system_matrix_,
                     const MatrixType &inverse_diagonal_matrix_)
    : system_matrix(system_matrix_)
    , inverse_diagonal_matrix(inverse_diagonal_matrix_)
  {}

  template <typename VectorType>
  void
  step(VectorType &u, const VectorType &v) const
  {
    // u = u - omega*inverse_diagonal_matrix*(system_matrix*u-v)
    const double omega = 1.;
    VectorType   tmp_1(v.size());
    system_matrix.vmult(tmp_1, u);
    tmp_1 -= v;
    VectorType tmp_2(u.size());
    inverse_diagonal_matrix.vmult(tmp_2, tmp_1);
    tmp_2 *= omega;
    u -= tmp_2;
  }

  template <typename VectorType>
  void
  Tstep(VectorType &u, const VectorType &v) const
  {
    AssertThrow(false, ExcNotImplemented());
  }

private:
  const MatrixType &system_matrix;
  const MatrixType &inverse_diagonal_matrix;
};


void
test(Utilities::CUDA::Handle &cuda_handle)
{
  // Build the sparse matrix on the host
  const unsigned int   problem_size = 10;
  unsigned int         size         = (problem_size - 1) * (problem_size - 1);
  FDMatrix             testproblem(problem_size, problem_size);
  SparsityPattern      structure(size, size, 5);
  SparseMatrix<double> A;
  testproblem.five_point_structure(structure);
  structure.compress();
  A.reinit(structure);
  testproblem.five_point(A);
  SparseMatrix<double> A_diagonal_inverse;
  A_diagonal_inverse.reinit(structure);
  for (unsigned int i = 0; i < size; ++i)
    A_diagonal_inverse(i, i) = 1. / A(i, i);

  // Solve on the host
  RelaxationOperator<SparseMatrix<double>> relaxation_operator(
    A, A_diagonal_inverse);
  SolverControl      control(1000, 1.e-3);
  SolverRelaxation<> relaxation_host(control);
  Vector<double>     sol_host(size);
  Vector<double>     rhs_host(size);
  for (unsigned int i = 0; i < size; ++i)
    rhs_host[i] = static_cast<double>(i);
  relaxation_host.solve(A, sol_host, rhs_host, relaxation_operator);

  // Solve on the device
  CUDAWrappers::SparseMatrix<double> A_dev(cuda_handle, A);
  CUDAWrappers::SparseMatrix<double> A_diagonal_inverse_dev(cuda_handle,
                                                            A_diagonal_inverse);
  RelaxationOperator<CUDAWrappers::SparseMatrix<double>>
                                              relaxation_operator_dev(A_dev, A_diagonal_inverse_dev);
  LinearAlgebra::CUDAWrappers::Vector<double> sol_dev(size);
  LinearAlgebra::CUDAWrappers::Vector<double> rhs_dev(size);
  LinearAlgebra::ReadWriteVector<double>      rw_vector(size);
  for (unsigned int i = 0; i < size; ++i)
    rw_vector[i] = static_cast<double>(i);
  rhs_dev.import(rw_vector, VectorOperation::insert);
  SolverRelaxation<LinearAlgebra::CUDAWrappers::Vector<double>> relaxation_dev(
    control);
  relaxation_dev.solve(A_dev, sol_dev, rhs_dev, relaxation_operator_dev);

  // Check the result
  rw_vector.import(sol_dev, VectorOperation::insert);
  for (unsigned int i = 0; i < size; ++i)
    AssertThrow(std::fabs(rw_vector[i] - sol_host[i]) < 1e-8,
                ExcInternalError());
}

int
main()
{
  initlog();
  deallog.depth_console(10);

  Utilities::CUDA::Handle cuda_handle;
  test(cuda_handle);

  deallog << "OK" << std::endl;

  return 0;
}
