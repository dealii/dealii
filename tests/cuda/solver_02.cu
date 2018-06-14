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

// Check that CUDA direct solvers work

#include <deal.II/base/cuda.h>

#include <deal.II/lac/cuda_solver_direct.h>
#include <deal.II/lac/cuda_sparse_matrix.h>
#include <deal.II/lac/cuda_vector.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/vector.h>

#include "../testmatrix.h"
#include "../tests.h"


void
test(Utilities::CUDA::Handle &cuda_handle)
{
  // Create the matrix on the host.
  dealii::SparsityPattern                sparsity_pattern;
  dealii::SparseMatrix<double>           matrix;
  unsigned int const                     size = 30;
  std::vector<std::vector<unsigned int>> column_indices(size);
  for (unsigned int i = 0; i < size; ++i)
    {
      unsigned int j_max = std::min(size, i + 2);
      unsigned int j_min = (i == 0) ? 0 : i - 1;
      for (unsigned int j = j_min; j < j_max; ++j)
        column_indices[i].emplace_back(j);
    }
  sparsity_pattern.copy_from(size,
                             size,
                             column_indices.begin(),
                             column_indices.end());
  matrix.reinit(sparsity_pattern);
  for (unsigned int i = 0; i < size; ++i)
    {
      unsigned int j_max = std::min(size - 1, i + 1);
      unsigned int j_min = (i == 0) ? 0 : i - 1;
      matrix.set(i, j_min, -1.);
      matrix.set(i, j_max, -1.);
      matrix.set(i, i, 4.);
    }

  // Generate a random solution and then compute the rhs
  dealii::Vector<double> sol_ref(size);
  for (auto &val : sol_ref)
    val = random_value(5., 15.);

  dealii::Vector<double> rhs(size);
  matrix.vmult(rhs, sol_ref);

  // Move the matrix and the rhs to the host
  CUDAWrappers::SparseMatrix<double> matrix_dev(cuda_handle, matrix);

  LinearAlgebra::CUDAWrappers::Vector<double> rhs_dev(size);
  LinearAlgebra::ReadWriteVector<double>      rhs_host(size);
  std::copy(rhs.begin(), rhs.end(), rhs_host.begin());
  rhs_dev.import(rhs_host, VectorOperation::insert);

  LinearAlgebra::CUDAWrappers::Vector<double> solution_dev(size);
  const std::array<std::string, 3>            solver_names{"Cholesky",
                                                "LU_dense",
                                                "LU_host"};

  for (auto solver_type : solver_names)
    {
      // Solve on the device
      CUDAWrappers::SolverDirect<double>::AdditionalData data(solver_type);
      SolverControl                                      solver_control;

      CUDAWrappers::SolverDirect<double> solver(cuda_handle,
                                                solver_control,
                                                data);
      solver.solve(matrix_dev, solution_dev, rhs_dev);

      // Move the result back to the host
      LinearAlgebra::ReadWriteVector<double> solution_host(size);
      solution_host.import(solution_dev, VectorOperation::insert);

      // Check the result
      for (unsigned int i = 0; i < size; ++i)
        AssertThrow(std::abs(solution_host[i] - sol_ref[i]) < 1e-12,
                    ExcInternalError());
      deallog << solver_type << std::endl;
    }
}

int
main()
{
  initlog();
  deallog.depth_console(0);

  Utilities::CUDA::Handle cuda_handle;
  test(cuda_handle);

  deallog << "OK" << std::endl;

  return 0;
}
