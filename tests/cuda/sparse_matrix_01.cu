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

// Check multiplications and norms

#include <deal.II/base/cuda.h>
#include <deal.II/base/exceptions.h>

#include <deal.II/lac/cuda_sparse_matrix.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

#include "../testmatrix.h"


void
check_matrix(SparseMatrix<double> const &        A,
             CUDAWrappers::SparseMatrix<double> &A_dev)
{
  cudaError_t cuda_error_code;
  double *    val_dev          = nullptr;
  int *       column_index_dev = nullptr;
  int *       row_ptr_dev      = nullptr;
  std::tie(val_dev, column_index_dev, row_ptr_dev, std::ignore, std::ignore) =
    A_dev.get_cusparse_matrix();

  int                 nnz = A_dev.n_nonzero_elements();
  std::vector<double> val_host(nnz);
  cuda_error_code = cudaMemcpy(&val_host[0],
                               val_dev,
                               nnz * sizeof(double),
                               cudaMemcpyDeviceToHost);
  AssertCuda(cuda_error_code);

  std::vector<int> column_index_host(nnz);
  cuda_error_code = cudaMemcpy(&column_index_host[0],
                               column_index_dev,
                               nnz * sizeof(int),
                               cudaMemcpyDeviceToHost);
  AssertCuda(cuda_error_code);

  int const        n_rows = A_dev.m() + 1;
  std::vector<int> row_ptr_host(n_rows + 1);
  cuda_error_code = cudaMemcpy(&row_ptr_host[0],
                               row_ptr_dev,
                               (A_dev.m() + 1) * sizeof(int),
                               cudaMemcpyDeviceToHost);
  AssertCuda(cuda_error_code);

  for (int i = 0; i < n_rows; ++i)
    for (int j = row_ptr_host[i]; j < row_ptr_host[i + 1]; ++j)
      AssertThrow(std::abs(val_host[j] - A(i, column_index_host[j])) < 1e-15,
                  ExcInternalError());
}

void
check_vector(Vector<double> const &                        a,
             LinearAlgebra::ReadWriteVector<double> const &b)
{
  unsigned int size = a.size();
  for (unsigned int i = 0; i < size; ++i)
    AssertThrow(std::abs(a[i] - b[i]) < 1e-15, ExcInternalError());
}

void
test(Utilities::CUDA::Handle &cuda_handle)
{
  // Build the sparse matrix on the host
  const unsigned int   size = 10;
  unsigned int         dim  = (size - 1) * (size - 1);
  FDMatrix             testproblem(size, size);
  SparsityPattern      structure(dim, dim, 5);
  SparseMatrix<double> A;
  testproblem.five_point_structure(structure);
  structure.compress();
  A.reinit(structure);
  testproblem.five_point(A, true);

  // Create the sparse matrix on the device
  CUDAWrappers::SparseMatrix<double> A_dev(cuda_handle, A);
  check_matrix(A, A_dev);

  AssertDimension(A.m(), A_dev.m());
  AssertDimension(A.n(), A_dev.n());

  // Multiply by a constant
  A *= 2.;
  A_dev *= 2.;
  check_matrix(A, A_dev);

  // Divide by a constant
  A /= 2.;
  A_dev /= 2.;
  check_matrix(A, A_dev);

  // Matrix-vector multiplication
  const unsigned int vector_size = A.n();
  Vector<double>     dst(vector_size);
  Vector<double>     src(vector_size);
  for (unsigned int i = 0; i < vector_size; ++i)
    src[i] = i;
  A.vmult(dst, src);
  LinearAlgebra::CUDAWrappers::Vector<double> dst_dev(vector_size);
  LinearAlgebra::CUDAWrappers::Vector<double> src_dev(vector_size);
  LinearAlgebra::ReadWriteVector<double>      read_write(vector_size);
  for (unsigned int i = 0; i < vector_size; ++i)
    read_write[i] = i;
  src_dev.import(read_write, VectorOperation::insert);
  A_dev.vmult(dst_dev, src_dev);
  read_write.import(dst_dev, VectorOperation::insert);
  check_vector(dst, read_write);

  // Transpose matrix-vector multiplication
  A.Tvmult(dst, src);
  A_dev.Tvmult(dst_dev, src_dev);
  read_write.import(dst_dev, VectorOperation::insert);
  check_vector(dst, read_write);

  // Matrix-vector multiplication and add
  A.vmult_add(dst, src);
  A_dev.vmult_add(dst_dev, src_dev);
  read_write.import(dst_dev, VectorOperation::insert);
  check_vector(dst, read_write);

  // Transpose matrix-vector multiplication and add
  A.Tvmult_add(dst, src);
  A_dev.Tvmult_add(dst_dev, src_dev);
  read_write.import(dst_dev, VectorOperation::insert);
  check_vector(dst, read_write);

  // Matrix norm square
  double value      = A.matrix_norm_square(src);
  double value_host = A_dev.matrix_norm_square(src_dev);
  AssertThrow(std::abs(value - value_host) < 1e-15, ExcInternalError());

  // Matrix scalar product (reuse dst and src but they are both input)
  value      = A.matrix_scalar_product(dst, src);
  value_host = A_dev.matrix_scalar_product(dst_dev, src_dev);
  AssertThrow(std::abs(value - value_host) < 1e-15, ExcInternalError());

  // Compute the residual
  Vector<double> b(src);
  for (unsigned int i = 0; i < vector_size; ++i)
    {
      b[i]          = i;
      src[i]        = i;
      read_write[i] = i;
    }
  LinearAlgebra::CUDAWrappers::Vector<double> b_dev(vector_size);
  b_dev.import(read_write, VectorOperation::insert);
  src_dev.import(read_write, VectorOperation::insert);
  value      = A.residual(dst, src, b);
  value_host = A_dev.residual(dst_dev, src_dev, b_dev);
  AssertThrow(std::abs(value - value_host) < 1e-15, ExcInternalError());
  read_write.import(dst_dev, VectorOperation::insert);
  check_vector(dst, read_write);

  // Compute L1 norm
  value      = A.l1_norm();
  value_host = A_dev.l1_norm();
  AssertThrow(std::abs(value - value_host) < 1e-15, ExcInternalError());

  // Compute Linfty norm
  value      = A.linfty_norm();
  value_host = A_dev.linfty_norm();
  AssertThrow(std::abs(value - value_host) < 1e-15, ExcInternalError());

  // Compute Frobenius norm
  value      = A.frobenius_norm();
  value_host = A_dev.frobenius_norm();
  AssertThrow(std::abs(value - value_host) < 1e-15, ExcInternalError());

  // Compute L1 norm second test
  SparsityPattern sparsity_pattern(vector_size, vector_size, 3);
  for (unsigned int i = 0; i < vector_size; ++i)
    {
      sparsity_pattern.add(i, 0);
      sparsity_pattern.add(i, i);
      if (i < vector_size - 1)
        sparsity_pattern.add(i, i + 1);
    }
  sparsity_pattern.compress();
  SparseMatrix<double> B(sparsity_pattern);
  for (unsigned int i = 0; i < vector_size; ++i)
    {
      B.set(i, 0, 1);
      B.set(i, i, 1);
      if (i < vector_size - 1)
        B.set(i, i + 1, 1);
    }
  CUDAWrappers::SparseMatrix<double> B_dev(cuda_handle, B);
  value      = B.l1_norm();
  value_host = B_dev.l1_norm();
  AssertThrow(std::abs(value - value_host) < 1e-15, ExcInternalError());
}

int
main()
{
  initlog();
  deallog.depth_console(0);

  init_cuda();

  Utilities::CUDA::Handle cuda_handle;

  test(cuda_handle);

  deallog << "OK" << std::endl;

  return 0;
}
