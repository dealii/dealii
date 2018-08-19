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
// and PreconditionILU

#include <deal.II/base/cuda.h>
#include <deal.II/base/exceptions.h>

#include <deal.II/lac/cuda_precondition.h>
#include <deal.II/lac/cuda_sparse_matrix.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>

#include "../testmatrix.h"
#include "../tests.h"

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
  A.print(std::cout);

  // Solve on the device
  CUDAWrappers::SparseMatrix<Number>          A_dev(cuda_handle, A);
  LinearAlgebra::CUDAWrappers::Vector<Number> sol_dev(size);
  LinearAlgebra::CUDAWrappers::Vector<Number> rhs_dev(size);
  LinearAlgebra::ReadWriteVector<Number>      rw_vector(size);
  for (unsigned int i = 0; i < size; ++i)
    rw_vector[i] = static_cast<Number>(i);
  rhs_dev.import(rw_vector, VectorOperation::insert);
  SolverControl                                         control(100, 1.e-10);
  SolverCG<LinearAlgebra::CUDAWrappers::Vector<Number>> cg_dev(control);

  CUDAWrappers::PreconditionILU<Number> prec_ilu(cuda_handle);
  prec_ilu.initialize(A_dev);

  cg_dev.solve(A_dev, sol_dev, rhs_dev, prec_ilu);

  // Check the result
  rw_vector.import(sol_dev, VectorOperation::insert);
  for (unsigned int i = 0; i < size; ++i)
    deallog << rw_vector[i] << std::endl;
}

int
main()
{
  initlog();
  deallog.depth_console(0);

  Utilities::CUDA::Handle cuda_handle;
  test<float>(cuda_handle);
  test<double>(cuda_handle);

  deallog << "OK" << std::endl;

  return 0;
}
