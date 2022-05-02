// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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

#include <deal.II/base/cuda_size.h>
#include <deal.II/base/exceptions.h>

#include <deal.II/lac/cuda_atomic.h>
#include <deal.II/lac/cuda_sparse_matrix.h>

#ifdef DEAL_II_WITH_CUDA

#  include <cusparse.h>

DEAL_II_NAMESPACE_OPEN

namespace CUDAWrappers
{
  namespace internal
  {
    template <typename Number>
    __global__ void
    scale(Number *                                       val,
          const Number                                   a,
          const typename SparseMatrix<Number>::size_type N)
    {
      const typename SparseMatrix<Number>::size_type idx =
        threadIdx.x + blockIdx.x * blockDim.x;
      if (idx < N)
        val[idx] *= a;
    }



    void
    create_sp_mat_descr(int                   m,
                        int                   n,
                        int                   nnz,
                        const float *         A_val_dev,
                        const int *           A_row_ptr_dev,
                        const int *           A_column_index_dev,
                        cusparseSpMatDescr_t &sp_descr)
    {
      cusparseStatus_t error_code = cusparseCreateCsr(
        &sp_descr,
        m,
        n,
        nnz,
        reinterpret_cast<void *>(const_cast<int *>(A_row_ptr_dev)),
        reinterpret_cast<void *>(const_cast<int *>(A_column_index_dev)),
        reinterpret_cast<void *>(const_cast<float *>(A_val_dev)),
        CUSPARSE_INDEX_32I,
        CUSPARSE_INDEX_32I,
        CUSPARSE_INDEX_BASE_ZERO,
        CUDA_R_32F);
      AssertCusparse(error_code);
    }



    void
    create_sp_mat_descr(int                   m,
                        int                   n,
                        int                   nnz,
                        const double *        A_val_dev,
                        const int *           A_row_ptr_dev,
                        const int *           A_column_index_dev,
                        cusparseSpMatDescr_t &sp_descr)
    {
      cusparseStatus_t error_code = cusparseCreateCsr(
        &sp_descr,
        m,
        n,
        nnz,
        reinterpret_cast<void *>(const_cast<int *>(A_row_ptr_dev)),
        reinterpret_cast<void *>(const_cast<int *>(A_column_index_dev)),
        reinterpret_cast<void *>(const_cast<double *>(A_val_dev)),
        CUSPARSE_INDEX_32I,
        CUSPARSE_INDEX_32I,
        CUSPARSE_INDEX_BASE_ZERO,
        CUDA_R_64F);
      AssertCusparse(error_code);
    }



    void
    csrmv(cusparseHandle_t           handle,
          bool                       transpose,
          int                        m,
          int                        n,
          const cusparseSpMatDescr_t sp_descr,
          const float *              x,
          bool                       add,
          float *                    y)
    {
      float               alpha = 1.;
      float               beta  = add ? 1. : 0.;
      cusparseOperation_t cusparse_operation =
        transpose ? CUSPARSE_OPERATION_TRANSPOSE :
                    CUSPARSE_OPERATION_NON_TRANSPOSE;

      // Move the data to cuSPARSE data type
      cusparseDnVecDescr_t x_cuvec;
      cusparseStatus_t     error_code =
        cusparseCreateDnVec(&x_cuvec,
                            n,
                            reinterpret_cast<void *>(const_cast<float *>(x)),
                            CUDA_R_32F);
      AssertCusparse(error_code);

      cusparseDnVecDescr_t y_cuvec;
      error_code =
        cusparseCreateDnVec(&y_cuvec,
                            m,
                            reinterpret_cast<void *>(const_cast<float *>(y)),
                            CUDA_R_32F);
      AssertCusparse(error_code);

      // This function performs y = alpha*op(A)*x + beta*y
      size_t buffer_size = 0;
      error_code         = cusparseSpMV_bufferSize(handle,
                                           cusparse_operation,
                                           &alpha,
                                           sp_descr,
                                           x_cuvec,
                                           &beta,
                                           y_cuvec,
                                           CUDA_R_32F,
                                           CUSPARSE_MV_ALG_DEFAULT,
                                           &buffer_size);

      void *      buffer          = nullptr;
      cudaError_t cuda_error_code = cudaMalloc(&buffer, buffer_size);
      AssertCuda(cuda_error_code);

      // execute SpMV
      error_code = cusparseSpMV(handle,
                                cusparse_operation,
                                &alpha,
                                sp_descr,
                                x_cuvec,
                                &beta,
                                y_cuvec,
                                CUDA_R_32F,
                                CUSPARSE_MV_ALG_DEFAULT,
                                buffer);
      AssertCusparse(error_code);

      cuda_error_code = cudaFree(buffer);
      AssertCuda(cuda_error_code);
      error_code = cusparseDestroyDnVec(x_cuvec);
      AssertCusparse(error_code);
      error_code = cusparseDestroyDnVec(y_cuvec);
      AssertCusparse(error_code);
    }



    void
    csrmv(cusparseHandle_t           handle,
          bool                       transpose,
          int                        m,
          int                        n,
          const cusparseSpMatDescr_t sp_descr,
          const double *             x,
          bool                       add,
          double *                   y)
    {
      double              alpha = 1.;
      double              beta  = add ? 1. : 0.;
      cusparseOperation_t cusparse_operation =
        transpose ? CUSPARSE_OPERATION_TRANSPOSE :
                    CUSPARSE_OPERATION_NON_TRANSPOSE;

      // Move the data to cuSPARSE data type
      cusparseDnVecDescr_t x_cuvec;
      cusparseStatus_t     error_code =
        cusparseCreateDnVec(&x_cuvec,
                            n,
                            reinterpret_cast<void *>(const_cast<double *>(x)),
                            CUDA_R_64F);
      AssertCusparse(error_code);

      cusparseDnVecDescr_t y_cuvec;
      error_code =
        cusparseCreateDnVec(&y_cuvec,
                            m,
                            reinterpret_cast<void *>(const_cast<double *>(y)),
                            CUDA_R_64F);
      AssertCusparse(error_code);

      // This function performs y = alpha*op(A)*x + beta*y
      size_t buffer_size = 0;
      error_code         = cusparseSpMV_bufferSize(handle,
                                           cusparse_operation,
                                           &alpha,
                                           sp_descr,
                                           x_cuvec,
                                           &beta,
                                           y_cuvec,
                                           CUDA_R_64F,
                                           CUSPARSE_MV_ALG_DEFAULT,
                                           &buffer_size);

      void *      buffer          = nullptr;
      cudaError_t cuda_error_code = cudaMalloc(&buffer, buffer_size);
      AssertCuda(cuda_error_code);

      // execute SpMV
      error_code = cusparseSpMV(handle,
                                cusparse_operation,
                                &alpha,
                                sp_descr,
                                x_cuvec,
                                &beta,
                                y_cuvec,
                                CUDA_R_64F,
                                CUSPARSE_MV_ALG_DEFAULT,
                                buffer);
      AssertCusparse(error_code);

      cuda_error_code = cudaFree(buffer);
      AssertCuda(cuda_error_code);
      error_code = cusparseDestroyDnVec(x_cuvec);
      AssertCusparse(error_code);
      error_code = cusparseDestroyDnVec(y_cuvec);
      AssertCusparse(error_code);
    }



    template <typename Number>
    __global__ void
    l1_norm(const typename SparseMatrix<Number>::size_type n_rows,
            const Number *                                 val_dev,
            const int *                                    column_index_dev,
            const int *                                    row_ptr_dev,
            Number *                                       sums)
    {
      const typename SparseMatrix<Number>::size_type row =
        threadIdx.x + blockIdx.x * blockDim.x;

      if (row < n_rows)
        {
          for (int j = row_ptr_dev[row]; j < row_ptr_dev[row + 1]; ++j)
            atomicAdd(&sums[column_index_dev[j]], abs(val_dev[j]));
        }
    }



    template <typename Number>
    __global__ void
    linfty_norm(const typename SparseMatrix<Number>::size_type n_rows,
                const Number *                                 val_dev,
                const int *                                    column_index_dev,
                const int *                                    row_ptr_dev,
                Number *                                       sums)
    {
      const typename SparseMatrix<Number>::size_type row =
        threadIdx.x + blockIdx.x * blockDim.x;

      if (row < n_rows)
        {
          sums[row] = (Number)0.;
          for (int j = row_ptr_dev[row]; j < row_ptr_dev[row + 1]; ++j)
            sums[row] += abs(val_dev[j]);
        }
    }
  } // namespace internal



  template <typename Number>
  SparseMatrix<Number>::SparseMatrix()
    : nnz(0)
    , n_rows(0)
    , val_dev(nullptr, Utilities::CUDA::delete_device_data<Number>)
    , column_index_dev(nullptr, Utilities::CUDA::delete_device_data<int>)
    , row_ptr_dev(nullptr, Utilities::CUDA::delete_device_data<int>)
    , descr(nullptr)
    , sp_descr(nullptr)
  {}



  template <typename Number>
  SparseMatrix<Number>::SparseMatrix(
    Utilities::CUDA::Handle &             handle,
    const ::dealii::SparseMatrix<Number> &sparse_matrix_host)
    : val_dev(nullptr, Utilities::CUDA::delete_device_data<Number>)
    , column_index_dev(nullptr, Utilities::CUDA::delete_device_data<int>)
    , row_ptr_dev(nullptr, Utilities::CUDA::delete_device_data<int>)
    , descr(nullptr)
    , sp_descr(nullptr)
  {
    reinit(handle, sparse_matrix_host);
  }



  template <typename Number>
  SparseMatrix<Number>::SparseMatrix(CUDAWrappers::SparseMatrix<Number> &&other)
    : cusparse_handle(other.cusparse_handle)
    , nnz(other.nnz)
    , n_rows(other.n_rows)
    , n_cols(other.n_cols)
    , val_dev(std::move(other.val_dev))
    , column_index_dev(std::move(other.column_index_dev))
    , row_ptr_dev(std::move(other.row_ptr_dev))
    , descr(other.descr)
    , sp_descr(other.sp_descr)
  {
    other.nnz      = 0;
    other.n_rows   = 0;
    other.n_cols   = 0;
    other.descr    = nullptr;
    other.sp_descr = nullptr;
  }



  template <typename Number>
  SparseMatrix<Number>::~SparseMatrix<Number>()
  {
    if (descr != nullptr)
      {
        const cusparseStatus_t cusparse_error_code =
          cusparseDestroyMatDescr(descr);
        AssertNothrowCusparse(cusparse_error_code);
        descr = nullptr;
      }

    if (sp_descr != nullptr)
      {
        const cusparseStatus_t cusparse_error_code =
          cusparseDestroySpMat(sp_descr);
        AssertNothrowCusparse(cusparse_error_code);
        sp_descr = nullptr;
      }

    nnz    = 0;
    n_rows = 0;
  }



  template <typename Number>
  SparseMatrix<Number> &
  SparseMatrix<Number>::operator=(SparseMatrix<Number> &&other)
  {
    cusparse_handle  = other.cusparse_handle;
    nnz              = other.nnz;
    n_rows           = other.n_rows;
    n_cols           = other.n_cols;
    val_dev          = std::move(other.val_dev);
    column_index_dev = std::move(other.column_index_dev);
    row_ptr_dev      = std::move(other.row_ptr_dev);
    descr            = other.descr;
    sp_descr         = other.sp_descr;

    other.nnz      = 0;
    other.n_rows   = 0;
    other.n_cols   = 0;
    other.descr    = nullptr;
    other.sp_descr = nullptr;

    return *this;
  }



  template <typename Number>
  void
  SparseMatrix<Number>::reinit(
    Utilities::CUDA::Handle &             handle,
    const ::dealii::SparseMatrix<Number> &sparse_matrix_host)
  {
    cusparse_handle                  = handle.cusparse_handle;
    nnz                              = sparse_matrix_host.n_nonzero_elements();
    n_rows                           = sparse_matrix_host.m();
    n_cols                           = sparse_matrix_host.n();
    unsigned int const  row_ptr_size = n_rows + 1;
    std::vector<Number> val;
    val.reserve(nnz);
    std::vector<int> column_index;
    column_index.reserve(nnz);
    std::vector<int> row_ptr(row_ptr_size, 0);

    // dealii::SparseMatrix stores the diagonal first in each row so we need to
    // do some reordering
    for (int row = 0; row < n_rows; ++row)
      {
        auto         p_end   = sparse_matrix_host.end(row);
        unsigned int counter = 0;
        for (auto p = sparse_matrix_host.begin(row); p != p_end; ++p)
          {
            val.emplace_back(p->value());
            column_index.emplace_back(p->column());
            ++counter;
          }
        row_ptr[row + 1] = row_ptr[row] + counter;

        // Sort the elements in the row
        unsigned int const offset     = row_ptr[row];
        int const          diag_index = column_index[offset];
        Number             diag_elem  = sparse_matrix_host.diag_element(row);
        unsigned int       pos        = 1;
        while ((column_index[offset + pos] < row) && (pos < counter))
          {
            val[offset + pos - 1]          = val[offset + pos];
            column_index[offset + pos - 1] = column_index[offset + pos];
            ++pos;
          }
        val[offset + pos - 1]          = diag_elem;
        column_index[offset + pos - 1] = diag_index;
      }

    // Copy the elements to the gpu
    val_dev.reset(Utilities::CUDA::allocate_device_data<Number>(nnz));
    cudaError_t error_code = cudaMemcpy(val_dev.get(),
                                        val.data(),
                                        nnz * sizeof(Number),
                                        cudaMemcpyHostToDevice);
    AssertCuda(error_code);

    // Copy the column indices to the gpu
    column_index_dev.reset(Utilities::CUDA::allocate_device_data<int>(nnz));
    AssertCuda(error_code);
    error_code = cudaMemcpy(column_index_dev.get(),
                            column_index.data(),
                            nnz * sizeof(int),
                            cudaMemcpyHostToDevice);
    AssertCuda(error_code);

    // Copy the row pointer to the gpu
    row_ptr_dev.reset(Utilities::CUDA::allocate_device_data<int>(row_ptr_size));
    AssertCuda(error_code);
    error_code = cudaMemcpy(row_ptr_dev.get(),
                            row_ptr.data(),
                            row_ptr_size * sizeof(int),
                            cudaMemcpyHostToDevice);
    AssertCuda(error_code);

    // Create the matrix descriptor
    cusparseStatus_t cusparse_error_code = cusparseCreateMatDescr(&descr);
    AssertCusparse(cusparse_error_code);
    cusparse_error_code =
      cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
    AssertCusparse(cusparse_error_code);
    cusparse_error_code =
      cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);
    AssertCusparse(cusparse_error_code);

    // Create the sparse matrix descriptor
    internal::create_sp_mat_descr(n_rows,
                                  n_cols,
                                  nnz,
                                  val_dev.get(),
                                  row_ptr_dev.get(),
                                  column_index_dev.get(),
                                  sp_descr);
  }



  template <typename Number>
  SparseMatrix<Number> &
  SparseMatrix<Number>::operator*=(const Number factor)
  {
    AssertIsFinite(factor);
    const int n_blocks = 1 + (nnz - 1) / block_size;
    internal::scale<Number>
      <<<n_blocks, block_size>>>(val_dev.get(), factor, nnz);
    AssertCudaKernel();

    return *this;
  }



  template <typename Number>
  SparseMatrix<Number> &
  SparseMatrix<Number>::operator/=(const Number factor)
  {
    AssertIsFinite(factor);
    Assert(factor != Number(0.), ExcZero());
    const int n_blocks = 1 + (nnz - 1) / block_size;
    internal::scale<Number>
      <<<n_blocks, block_size>>>(val_dev.get(), 1. / factor, nnz);
    AssertCudaKernel();

    return *this;
  }



  template <typename Number>
  void
  SparseMatrix<Number>::vmult(
    LinearAlgebra::CUDAWrappers::Vector<Number> &      dst,
    const LinearAlgebra::CUDAWrappers::Vector<Number> &src) const
  {
    internal::csrmv(cusparse_handle,
                    false,
                    n_rows,
                    n_cols,
                    sp_descr,
                    src.get_values(),
                    false,
                    dst.get_values());
  }



  template <typename Number>
  void
  SparseMatrix<Number>::Tvmult(
    LinearAlgebra::CUDAWrappers::Vector<Number> &      dst,
    const LinearAlgebra::CUDAWrappers::Vector<Number> &src) const
  {
    internal::csrmv(cusparse_handle,
                    true,
                    n_rows,
                    n_cols,
                    sp_descr,
                    src.get_values(),
                    false,
                    dst.get_values());
  }



  template <typename Number>
  void
  SparseMatrix<Number>::vmult_add(
    LinearAlgebra::CUDAWrappers::Vector<Number> &      dst,
    const LinearAlgebra::CUDAWrappers::Vector<Number> &src) const
  {
    internal::csrmv(cusparse_handle,
                    false,
                    n_rows,
                    n_cols,
                    sp_descr,
                    src.get_values(),
                    true,
                    dst.get_values());
  }



  template <typename Number>
  void
  SparseMatrix<Number>::Tvmult_add(
    LinearAlgebra::CUDAWrappers::Vector<Number> &      dst,
    const LinearAlgebra::CUDAWrappers::Vector<Number> &src) const
  {
    internal::csrmv(cusparse_handle,
                    true,
                    n_rows,
                    n_cols,
                    sp_descr,
                    src.get_values(),
                    true,
                    dst.get_values());
  }



  template <typename Number>
  Number
  SparseMatrix<Number>::matrix_norm_square(
    const LinearAlgebra::CUDAWrappers::Vector<Number> &v) const
  {
    LinearAlgebra::CUDAWrappers::Vector<Number> tmp = v;
    vmult(tmp, v);

    return v * tmp;
  }



  template <typename Number>
  Number
  SparseMatrix<Number>::matrix_scalar_product(
    const LinearAlgebra::CUDAWrappers::Vector<Number> &u,
    const LinearAlgebra::CUDAWrappers::Vector<Number> &v) const
  {
    LinearAlgebra::CUDAWrappers::Vector<Number> tmp = v;
    vmult(tmp, v);

    return u * tmp;
  }



  template <typename Number>
  Number
  SparseMatrix<Number>::residual(
    LinearAlgebra::CUDAWrappers::Vector<Number> &      dst,
    const LinearAlgebra::CUDAWrappers::Vector<Number> &x,
    const LinearAlgebra::CUDAWrappers::Vector<Number> &b) const
  {
    vmult(dst, x);
    dst.sadd(-1., 1., b);

    return dst.l2_norm();
  }



  template <typename Number>
  Number
  SparseMatrix<Number>::l1_norm() const
  {
    LinearAlgebra::CUDAWrappers::Vector<real_type> column_sums(n_cols);
    const int n_blocks = 1 + (nnz - 1) / block_size;
    internal::l1_norm<Number>
      <<<n_blocks, block_size>>>(n_rows,
                                 val_dev.get(),
                                 column_index_dev.get(),
                                 row_ptr_dev.get(),
                                 column_sums.get_values());
    AssertCudaKernel();

    return column_sums.linfty_norm();
  }



  template <typename Number>
  Number
  SparseMatrix<Number>::linfty_norm() const
  {
    LinearAlgebra::CUDAWrappers::Vector<real_type> row_sums(n_rows);
    const int n_blocks = 1 + (nnz - 1) / block_size;
    internal::linfty_norm<Number>
      <<<n_blocks, block_size>>>(n_rows,
                                 val_dev.get(),
                                 column_index_dev.get(),
                                 row_ptr_dev.get(),
                                 row_sums.get_values());
    AssertCudaKernel();

    return row_sums.linfty_norm();
  }



  template <typename Number>
  Number
  SparseMatrix<Number>::frobenius_norm() const
  {
    LinearAlgebra::CUDAWrappers::Vector<real_type> matrix_values(nnz);
    cudaError_t cuda_error = cudaMemcpy(matrix_values.get_values(),
                                        val_dev.get(),
                                        nnz * sizeof(Number),
                                        cudaMemcpyDeviceToDevice);

    return matrix_values.l2_norm();
  }



  template <typename Number>
  std::tuple<Number *, int *, int *, cusparseMatDescr_t, cusparseSpMatDescr_t>
  SparseMatrix<Number>::get_cusparse_matrix() const
  {
    return std::make_tuple(val_dev.get(),
                           column_index_dev.get(),
                           row_ptr_dev.get(),
                           descr,
                           sp_descr);
  }



  template class SparseMatrix<float>;
  template class SparseMatrix<double>;
} // namespace CUDAWrappers
DEAL_II_NAMESPACE_CLOSE

#endif
