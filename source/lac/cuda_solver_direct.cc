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

#include <deal.II/lac/cuda_solver_direct.h>

DEAL_II_NAMESPACE_OPEN

namespace CUDAWrappers
{
  namespace
  {
    void
    cusparsecsr2dense(cusparseHandle_t           cusparse_handle,
                      const SparseMatrix<float> &matrix,
                      float *                    dense_matrix_dev)
    {
      auto cusparse_matrix = matrix.get_cusparse_matrix();

      const cusparseStatus_t cusparse_error_code =
        cusparseScsr2dense(cusparse_handle,
                           matrix.m(),
                           matrix.n(),
                           std::get<3>(cusparse_matrix),
                           std::get<0>(cusparse_matrix),
                           std::get<2>(cusparse_matrix),
                           std::get<1>(cusparse_matrix),
                           dense_matrix_dev,
                           matrix.m());
      AssertCusparse(cusparse_error_code);
    }



    void
    cusparsecsr2dense(cusparseHandle_t            cusparse_handle,
                      const SparseMatrix<double> &matrix,
                      double *                    dense_matrix_dev)
    {
      auto cusparse_matrix = matrix.get_cusparse_matrix();

      const cusparseStatus_t cusparse_error_code =
        cusparseDcsr2dense(cusparse_handle,
                           matrix.m(),
                           matrix.n(),
                           std::get<3>(cusparse_matrix),
                           std::get<0>(cusparse_matrix),
                           std::get<2>(cusparse_matrix),
                           std::get<1>(cusparse_matrix),
                           dense_matrix_dev,
                           matrix.m());
      AssertCusparse(cusparse_error_code);
    }



    void
    cusolverDngetrf_buffer_size(cusolverDnHandle_t cusolver_dn_handle,
                                int                m,
                                int                n,
                                float *            dense_matrix_dev,
                                int &              workspace_size)
    {
      const cusolverStatus_t cusolver_error_code = cusolverDnSgetrf_bufferSize(
        cusolver_dn_handle, m, n, dense_matrix_dev, m, &workspace_size);
      AssertCusolver(cusolver_error_code);
    }



    void
    cusolverDngetrf_buffer_size(cusolverDnHandle_t cusolver_dn_handle,
                                int                m,
                                int                n,
                                double *           dense_matrix_dev,
                                int &              workspace_size)
    {
      const cusolverStatus_t cusolver_error_code = cusolverDnDgetrf_bufferSize(
        cusolver_dn_handle, m, n, dense_matrix_dev, m, &workspace_size);
      AssertCusolver(cusolver_error_code);
    }



    void
    cusolverDngetrf(cusolverDnHandle_t cusolver_dn_handle,
                    int                m,
                    int                n,
                    float *            dense_matrix_dev,
                    float *            workspace_dev,
                    int *              pivot_dev,
                    int *              info_dev)
    {
      const cusolverStatus_t cusolver_error_code =
        cusolverDnSgetrf(cusolver_dn_handle,
                         m,
                         n,
                         dense_matrix_dev,
                         m,
                         workspace_dev,
                         pivot_dev,
                         info_dev);
      AssertCusolver(cusolver_error_code);
    }



    void
    cusolverDngetrf(cusolverDnHandle_t cusolver_dn_handle,
                    int                m,
                    int                n,
                    double *           dense_matrix_dev,
                    double *           workspace_dev,
                    int *              pivot_dev,
                    int *              info_dev)
    {
      const cusolverStatus_t cusolver_error_code =
        cusolverDnDgetrf(cusolver_dn_handle,
                         m,
                         n,
                         dense_matrix_dev,
                         m,
                         workspace_dev,
                         pivot_dev,
                         info_dev);
      AssertCusolver(cusolver_error_code);
    }



    void
    cusolverDngetrs(cusolverDnHandle_t cusolver_dn_handle,
                    int                m,
                    float *            dense_matrix_dev,
                    int *              pivot_dev,
                    float *            b,
                    int *              info_dev)
    {
      const int              n_rhs = 1;
      const cusolverStatus_t cusolver_error_code =
        cusolverDnSgetrs(cusolver_dn_handle,
                         CUBLAS_OP_N,
                         m,
                         n_rhs,
                         dense_matrix_dev,
                         m,
                         pivot_dev,
                         b,
                         m,
                         info_dev);
      AssertCusolver(cusolver_error_code);
    }



    void
    cusolverDngetrs(cusolverDnHandle_t cusolver_dn_handle,
                    int                m,
                    double *           dense_matrix_dev,
                    int *              pivot_dev,
                    double *           b,
                    int *              info_dev)
    {
      const int              n_rhs = 1;
      const cusolverStatus_t cusolver_error_code =
        cusolverDnDgetrs(cusolver_dn_handle,
                         CUBLAS_OP_N,
                         m,
                         n_rhs,
                         dense_matrix_dev,
                         m,
                         pivot_dev,
                         b,
                         m,
                         info_dev);
      AssertCusolver(cusolver_error_code);
    }



    void
    cusolverSpcsrlsvluHost(cusolverSpHandle_t cusolver_sp_handle,
                           const unsigned int n_rows,
                           const unsigned int nnz,
                           cusparseMatDescr_t descr,
                           const float *      val_host,
                           const int *        row_ptr_host,
                           const int *        column_index_host,
                           const float *      b_host,
                           float *            x_host)
    {
      int                    singularity = 0;
      const cusolverStatus_t cusolver_error_code =
        cusolverSpScsrlsvluHost(cusolver_sp_handle,
                                n_rows,
                                nnz,
                                descr,
                                val_host,
                                row_ptr_host,
                                column_index_host,
                                b_host,
                                0.,
                                1,
                                x_host,
                                &singularity);
      AssertCusolver(cusolver_error_code);
      Assert(singularity == -1, ExcMessage("Coarse matrix is singular"));
    }



    void
    cusolverSpcsrlsvluHost(cusolverSpHandle_t cusolver_sp_handle,
                           const unsigned int n_rows,
                           unsigned int       nnz,
                           cusparseMatDescr_t descr,
                           const double *     val_host,
                           const int *        row_ptr_host,
                           const int *        column_index_host,
                           const double *     b_host,
                           double *           x_host)
    {
      int                    singularity = 0;
      const cusolverStatus_t cusolver_error_code =
        cusolverSpDcsrlsvluHost(cusolver_sp_handle,
                                n_rows,
                                nnz,
                                descr,
                                val_host,
                                row_ptr_host,
                                column_index_host,
                                b_host,
                                0.,
                                1,
                                x_host,
                                &singularity);
      AssertCusolver(cusolver_error_code);
      Assert(singularity == -1, ExcMessage("Coarse matrix is singular"));
    }



    void
    cholesky_factorization(cusolverSpHandle_t         cusolver_sp_handle,
                           const SparseMatrix<float> &matrix,
                           const float *              b,
                           float *                    x)
    {
      auto cusparse_matrix = matrix.get_cusparse_matrix();
      int  singularity     = 0;

      const cusolverStatus_t cusolver_error_code =
        cusolverSpScsrlsvchol(cusolver_sp_handle,
                              matrix.m(),
                              matrix.n_nonzero_elements(),
                              std::get<3>(cusparse_matrix),
                              std::get<0>(cusparse_matrix),
                              std::get<2>(cusparse_matrix),
                              std::get<1>(cusparse_matrix),
                              b,
                              0.,
                              0,
                              x,
                              &singularity);
      AssertCusolver(cusolver_error_code);
      Assert(singularity == -1, ExcMessage("Coarse matrix is not SPD"));
    }



    void
    cholesky_factorization(cusolverSpHandle_t          cusolver_sp_handle,
                           const SparseMatrix<double> &matrix,
                           const double *              b,
                           double *                    x)
    {
      auto cusparse_matrix = matrix.get_cusparse_matrix();
      int  singularity     = 0;

      const cusolverStatus_t cusolver_error_code =
        cusolverSpDcsrlsvchol(cusolver_sp_handle,
                              matrix.m(),
                              matrix.n_nonzero_elements(),
                              std::get<3>(cusparse_matrix),
                              std::get<0>(cusparse_matrix),
                              std::get<2>(cusparse_matrix),
                              std::get<1>(cusparse_matrix),
                              b,
                              0.,
                              0,
                              x,
                              &singularity);
      AssertCusolver(cusolver_error_code);
      Assert(singularity == -1, ExcMessage("Coarse matrix is not SPD"));
    }



    template <typename Number>
    void
    lu_factorization(cusparseHandle_t            cusparse_handle,
                     cusolverDnHandle_t          cusolver_dn_handle,
                     const SparseMatrix<Number> &matrix,
                     const Number *              b_dev,
                     Number *                    x_dev)
    {
      // Change the format of the matrix from sparse to dense
      unsigned int const m = matrix.m();
      unsigned int const n = matrix.n();
      Assert(m == n, ExcMessage("The matrix is not square"));
      Number *dense_matrix_dev;
      Utilities::CUDA::malloc(dense_matrix_dev, m * n);

      // Change the format of matrix to dense
      cusparsecsr2dense(cusparse_handle, matrix, dense_matrix_dev);

      // Create the working space
      int workspace_size = 0;
      cusolverDngetrf_buffer_size(
        cusolver_dn_handle, m, n, dense_matrix_dev, workspace_size);
      Assert(workspace_size > 0, ExcMessage("No workspace was allocated"));
      Number *workspace_dev;
      Utilities::CUDA::malloc(workspace_dev, workspace_size);

      // LU factorization
      int *pivot_dev;
      Utilities::CUDA::malloc(pivot_dev, m);
      int *info_dev;
      Utilities::CUDA::malloc(info_dev, 1);

      cusolverDngetrf(cusolver_dn_handle,
                      m,
                      n,
                      dense_matrix_dev,
                      workspace_dev,
                      pivot_dev,
                      info_dev);

#ifdef DEBUG
      int         info = 0;
      cudaError_t cuda_error_code_debug =
        cudaMemcpy(&info, info_dev, sizeof(int), cudaMemcpyDeviceToHost);
      AssertCuda(cuda_error_code_debug);
      Assert(info == 0,
             ExcMessage("There was a problem during the LU factorization"));
#endif

      // Solve Ax = b
      cudaError_t cuda_error_code =
        cudaMemcpy(x_dev, b_dev, m * sizeof(Number), cudaMemcpyDeviceToDevice);
      AssertCuda(cuda_error_code);
      cusolverDngetrs(
        cusolver_dn_handle, m, dense_matrix_dev, pivot_dev, x_dev, info_dev);
#ifdef DEBUG
      cuda_error_code =
        cudaMemcpy(&info, info_dev, sizeof(int), cudaMemcpyDeviceToHost);
      AssertCuda(cuda_error_code);
      Assert(info == 0, ExcMessage("There was a problem during the LU solve"));
#endif

      // Free the memory allocated
      Utilities::CUDA::free(dense_matrix_dev);
      Utilities::CUDA::free(workspace_dev);
      Utilities::CUDA::free(pivot_dev);
      Utilities::CUDA::free(info_dev);
    }



    template <typename Number>
    void
    lu_factorization(cusolverSpHandle_t          cusolver_sp_handle,
                     const SparseMatrix<Number> &matrix,
                     const Number *              b_dev,
                     Number *                    x_dev)
    {
      // cuSOLVER does not support LU factorization of sparse matrix on the
      // device, so we need to move everything to the host first and then back
      // to the host.
      const unsigned int  nnz    = matrix.n_nonzero_elements();
      const unsigned int  n_rows = matrix.m();
      std::vector<Number> val_host(nnz);
      std::vector<int>    column_index_host(nnz);
      std::vector<int>    row_ptr_host(n_rows + 1);
      auto                cusparse_matrix = matrix.get_cusparse_matrix();
      Utilities::CUDA::copy_to_host(std::get<0>(cusparse_matrix), val_host);
      Utilities::CUDA::copy_to_host(std::get<1>(cusparse_matrix),
                                    column_index_host);
      Utilities::CUDA::copy_to_host(std::get<2>(cusparse_matrix), row_ptr_host);
      std::vector<Number> b_host(n_rows);
      Utilities::CUDA::copy_to_host(b_dev, b_host);
      std::vector<Number> x_host(n_rows);
      Utilities::CUDA::copy_to_host(x_dev, x_host);

      cusolverSpcsrlsvluHost(cusolver_sp_handle,
                             n_rows,
                             nnz,
                             std::get<3>(cusparse_matrix),
                             val_host.data(),
                             row_ptr_host.data(),
                             column_index_host.data(),
                             b_host.data(),
                             x_host.data());

      // Move the solution back to the device
      Utilities::CUDA::copy_to_dev(x_host, x_dev);
    }
  } // namespace



  template <typename Number>
  SolverDirect<Number>::AdditionalData::AdditionalData(
    const std::string &solver_type)
    : solver_type(solver_type)
  {}



  template <typename Number>
  SolverDirect<Number>::SolverDirect(const Utilities::CUDA::Handle &handle,
                                     SolverControl &                cn,
                                     const AdditionalData &         data)
    : cuda_handle(handle)
    , solver_control(cn)
    , additional_data(data.solver_type)
  {}



  template <typename Number>
  SolverControl &
  SolverDirect<Number>::control() const
  {
    return solver_control;
  }



  template <typename Number>
  void
  SolverDirect<Number>::solve(
    const SparseMatrix<Number> &                       A,
    LinearAlgebra::CUDAWrappers::Vector<Number> &      x,
    const LinearAlgebra::CUDAWrappers::Vector<Number> &b)
  {
    if (additional_data.solver_type == "Cholesky")
      cholesky_factorization(cuda_handle.cusolver_sp_handle,
                             A,
                             b.get_values(),
                             x.get_values());
    else if (additional_data.solver_type == "LU_dense")
      lu_factorization(cuda_handle.cusparse_handle,
                       cuda_handle.cusolver_dn_handle,
                       A,
                       b.get_values(),
                       x.get_values());
    else if (additional_data.solver_type == "LU_host")
      lu_factorization(cuda_handle.cusolver_sp_handle,
                       A,
                       b.get_values(),
                       x.get_values());
    else
      AssertThrow(false,
                  ExcMessage("The provided solver name " +
                             additional_data.solver_type + " is invalid."));

    // Force the SolverControl object to report convergence
    solver_control.check(0, 0);
  }


  // Explicit Instanationation
  template class SolverDirect<float>;
  template class SolverDirect<double>;
} // namespace CUDAWrappers

DEAL_II_NAMESPACE_CLOSE
