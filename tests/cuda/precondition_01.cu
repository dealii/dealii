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

// Check that dealii::SolverCG works with CUDAWrappers::SparseMatrix
// and PreconditionIC

#include <deal.II/base/cuda.h>
#include <deal.II/base/exceptions.h>

#include <deal.II/lac/cuda_precondition.h>
#include <deal.II/lac/cuda_sparse_matrix.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>

#include "../tests.h"

#include "../testmatrix.h"

template <typename Number>
void
test(Utilities::CUDA::Handle &cuda_handle)
{
  // Build the sparse matrix on the host
  const unsigned int   problem_size = 10;
  unsigned int         size         = (problem_size - 1) * (problem_size - 1);
  FDMatrix             testproblem(problem_size, problem_size);
  SparsityPattern      structure(size, size, 5);
  SparseMatrix<Number> A;
  testproblem.five_point_structure(structure);
  structure.compress();
  A.reinit(structure);
  testproblem.five_point(A);

  // Solve on the device
  CUDAWrappers::SparseMatrix<Number>          A_dev(cuda_handle, A);
  LinearAlgebra::CUDAWrappers::Vector<Number> sol_dev(size);
  LinearAlgebra::CUDAWrappers::Vector<Number> rhs_dev(size);
  LinearAlgebra::ReadWriteVector<Number>      rw_vector(size);
  for (unsigned int i = 0; i < size; ++i)
    rw_vector[i] = static_cast<Number>(i);
  rhs_dev.import(rw_vector, VectorOperation::insert);
  const Number  tolerance = 1000. * std::numeric_limits<Number>::epsilon();
  SolverControl control(100, tolerance);
  SolverCG<LinearAlgebra::CUDAWrappers::Vector<Number>> cg_dev(control);

  CUDAWrappers::PreconditionIC<Number> prec_ic(cuda_handle);
  prec_ic.initialize(A_dev);

  cg_dev.solve(A_dev, sol_dev, rhs_dev, prec_ic);

  // Check the result
  LinearAlgebra::CUDAWrappers::Vector<Number> residual(size);
  A_dev.residual(residual, sol_dev, rhs_dev);
  Assert(residual.l2_norm() < 20 * tolerance, ExcInternalError());
  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();
  deallog << std::setprecision(10);
  deallog.depth_console(0);

  init_cuda();

  Utilities::CUDA::Handle cuda_handle;
  deallog << "Testing float" << std::endl;
  test<float>(cuda_handle);
  deallog << "Testing double" << std::endl;
  test<double>(cuda_handle);

  deallog << "OK" << std::endl;

  return 0;
}
