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

#include <deal.II/base/cuda.h>
#include <deal.II/base/exceptions.h>

#include <deal.II/lac/cuda_sparse_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/vector.h>

#include <memory>

#include "../testmatrix.h"
#include "../tests.h"

DEAL_II_NAMESPACE_OPEN

namespace CUDAWrappers
{
  /** \addtogroup CUDAWrappers
   *  @{
   */

  /**
   * Template wrapper for cusparse<t>csrilu02.
   * (https://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-csrilu02).
   *  function performs the solve phase of the incomplete-LU factorization with
   * 0 fill-in and no pivoting.
   */
  template <typename Number>
  cusparseStatus_t
  cusparseXcsrilu02(cusparseHandle_t         handle,
                    int                      m,
                    int                      nnz,
                    const cusparseMatDescr_t descrA,
                    Number *                 csrValA_valM,
                    const int *              csrRowPtrA,
                    const int *              csrColIndA,
                    csrilu02Info_t           info,
                    cusparseSolvePolicy_t    policy,
                    void *                   pBuffer)
  {
    AssertThrow(false, ExcNotImplemented());
    return CUSPARSE_STATUS_INVALID_VALUE;
  }

  template <>
  cusparseStatus_t
  cusparseXcsrilu02<float>(cusparseHandle_t         handle,
                           int                      m,
                           int                      nnz,
                           const cusparseMatDescr_t descrA,
                           float *                  csrValA_valM,
                           const int *              csrRowPtrA,
                           const int *              csrColIndA,
                           csrilu02Info_t           info,
                           cusparseSolvePolicy_t    policy,
                           void *                   pBuffer)
  {
    return cusparseScsrilu02(handle,
                             m,
                             nnz,
                             descrA,
                             csrValA_valM,
                             csrRowPtrA,
                             csrColIndA,
                             info,
                             policy,
                             pBuffer);
  }

  template <>
  cusparseStatus_t
  cusparseXcsrilu02<double>(cusparseHandle_t         handle,
                            int                      m,
                            int                      nnz,
                            const cusparseMatDescr_t descrA,
                            double *                 csrValA_valM,
                            const int *              csrRowPtrA,
                            const int *              csrColIndA,
                            csrilu02Info_t           info,
                            cusparseSolvePolicy_t    policy,
                            void *                   pBuffer)
  {
    return cusparseDcsrilu02(handle,
                             m,
                             nnz,
                             descrA,
                             csrValA_valM,
                             csrRowPtrA,
                             csrColIndA,
                             info,
                             policy,
                             pBuffer);
  }

  template <>
  cusparseStatus_t
  cusparseXcsrilu02<cuComplex>(cusparseHandle_t         handle,
                               int                      m,
                               int                      nnz,
                               const cusparseMatDescr_t descrA,
                               cuComplex *              csrValA_valM,
                               const int *              csrRowPtrA,
                               const int *              csrColIndA,
                               csrilu02Info_t           info,
                               cusparseSolvePolicy_t    policy,
                               void *                   pBuffer)
  {
    return cusparseCcsrilu02(handle,
                             m,
                             nnz,
                             descrA,
                             csrValA_valM,
                             csrRowPtrA,
                             csrColIndA,
                             info,
                             policy,
                             pBuffer);
  }

  template <>
  cusparseStatus_t
  cusparseXcsrilu02<cuDoubleComplex>(cusparseHandle_t         handle,
                                     int                      m,
                                     int                      nnz,
                                     const cusparseMatDescr_t descrA,
                                     cuDoubleComplex *        csrValA_valM,
                                     const int *              csrRowPtrA,
                                     const int *              csrColIndA,
                                     csrilu02Info_t           info,
                                     cusparseSolvePolicy_t    policy,
                                     void *                   pBuffer)
  {
    return cusparseZcsrilu02(handle,
                             m,
                             nnz,
                             descrA,
                             csrValA_valM,
                             csrRowPtrA,
                             csrColIndA,
                             info,
                             policy,
                             pBuffer);
  }



  /**
   * Template wrapper for cusparse<t>csrilu02_analysis.
   * (https://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-csrilu02_analysis).
   * This function performs the analysis phase of the incomplete-LU
   * factorization with 0 fill-in and no pivoting.
   */
  template <typename Number>
  cusparseStatus_t
  cusparseXcsrilu02_analysis(cusparseHandle_t         handle,
                             int                      m,
                             int                      nnz,
                             const cusparseMatDescr_t descrA,
                             const Number *           csrValA,
                             const int *              csrRowPtrA,
                             const int *              csrColIndA,
                             csrilu02Info_t           info,
                             cusparseSolvePolicy_t    policy,
                             void *                   pBuffer)
  {
    AssertThrow(false, ExcNotImplemented());
    return CUSPARSE_STATUS_INVALID_VALUE;
  }

  template <>
  cusparseStatus_t
  cusparseXcsrilu02_analysis<float>(cusparseHandle_t         handle,
                                    int                      m,
                                    int                      nnz,
                                    const cusparseMatDescr_t descrA,
                                    const float *            csrValA,
                                    const int *              csrRowPtrA,
                                    const int *              csrColIndA,
                                    csrilu02Info_t           info,
                                    cusparseSolvePolicy_t    policy,
                                    void *                   pBuffer)
  {
    return cusparseScsrilu02_analysis(handle,
                                      m,
                                      nnz,
                                      descrA,
                                      csrValA,
                                      csrRowPtrA,
                                      csrColIndA,
                                      info,
                                      policy,
                                      pBuffer);
  }

  template <>
  cusparseStatus_t
  cusparseXcsrilu02_analysis<double>(cusparseHandle_t         handle,
                                     int                      m,
                                     int                      nnz,
                                     const cusparseMatDescr_t descrA,
                                     const double *           csrValA,
                                     const int *              csrRowPtrA,
                                     const int *              csrColIndA,
                                     csrilu02Info_t           info,
                                     cusparseSolvePolicy_t    policy,
                                     void *                   pBuffer)
  {
    return cusparseDcsrilu02_analysis(handle,
                                      m,
                                      nnz,
                                      descrA,
                                      csrValA,
                                      csrRowPtrA,
                                      csrColIndA,
                                      info,
                                      policy,
                                      pBuffer);
  }

  template <>
  cusparseStatus_t
  cusparseXcsrilu02_analysis<cuComplex>(cusparseHandle_t         handle,
                                        int                      m,
                                        int                      nnz,
                                        const cusparseMatDescr_t descrA,
                                        const cuComplex *        csrValA,
                                        const int *              csrRowPtrA,
                                        const int *              csrColIndA,
                                        csrilu02Info_t           info,
                                        cusparseSolvePolicy_t    policy,
                                        void *                   pBuffer)
  {
    return cusparseCcsrilu02_analysis(handle,
                                      m,
                                      nnz,
                                      descrA,
                                      csrValA,
                                      csrRowPtrA,
                                      csrColIndA,
                                      info,
                                      policy,
                                      pBuffer);
  }

  template <>
  cusparseStatus_t
  cusparseXcsrilu02_analysis<cuDoubleComplex>(cusparseHandle_t         handle,
                                              int                      m,
                                              int                      nnz,
                                              const cusparseMatDescr_t descrA,
                                              const cuDoubleComplex *  csrValA,
                                              const int *           csrRowPtrA,
                                              const int *           csrColIndA,
                                              csrilu02Info_t        info,
                                              cusparseSolvePolicy_t policy,
                                              void *                pBuffer)
  {
    return cusparseZcsrilu02_analysis(handle,
                                      m,
                                      nnz,
                                      descrA,
                                      csrValA,
                                      csrRowPtrA,
                                      csrColIndA,
                                      info,
                                      policy,
                                      pBuffer);
  }



  /**
   * Template wrapper for cusparse<t>csrilu02_bufferSize.
   * (https://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-csrilu02_bufferSize).
   * This function returns size of the buffer used in computing the
   * incomplete-LU factorization with 0 fill-in and no pivoting.
   */
  template <typename Number>
  cusparseStatus_t
  cusparseXcsrilu02_bufferSize(cusparseHandle_t         handle,
                               int                      m,
                               int                      nnz,
                               const cusparseMatDescr_t descrA,
                               Number *                 csrValA,
                               const int *              csrRowPtrA,
                               const int *              csrColIndA,
                               csrilu02Info_t           info,
                               int *                    pBufferSizeInBytes)
  {
    AssertThrow(false, ExcNotImplemented());
    return CUSPARSE_STATUS_INVALID_VALUE;
  }

  template <>
  cusparseStatus_t
  cusparseXcsrilu02_bufferSize<float>(cusparseHandle_t         handle,
                                      int                      m,
                                      int                      nnz,
                                      const cusparseMatDescr_t descrA,
                                      float *                  csrValA,
                                      const int *              csrRowPtrA,
                                      const int *              csrColIndA,
                                      csrilu02Info_t           info,
                                      int *pBufferSizeInBytes)
  {
    return cusparseScsrilu02_bufferSize(handle,
                                        m,
                                        nnz,
                                        descrA,
                                        csrValA,
                                        csrRowPtrA,
                                        csrColIndA,
                                        info,
                                        pBufferSizeInBytes);
  }

  template <>
  cusparseStatus_t
  cusparseXcsrilu02_bufferSize<double>(cusparseHandle_t         handle,
                                       int                      m,
                                       int                      nnz,
                                       const cusparseMatDescr_t descrA,
                                       double *                 csrValA,
                                       const int *              csrRowPtrA,
                                       const int *              csrColIndA,
                                       csrilu02Info_t           info,
                                       int *pBufferSizeInBytes)
  {
    return cusparseDcsrilu02_bufferSize(handle,
                                        m,
                                        nnz,
                                        descrA,
                                        csrValA,
                                        csrRowPtrA,
                                        csrColIndA,
                                        info,
                                        pBufferSizeInBytes);
  }


  template <>
  cusparseStatus_t
  cusparseXcsrilu02_bufferSize<cuComplex>(cusparseHandle_t         handle,
                                          int                      m,
                                          int                      nnz,
                                          const cusparseMatDescr_t descrA,
                                          cuComplex *              csrValA,
                                          const int *              csrRowPtrA,
                                          const int *              csrColIndA,
                                          csrilu02Info_t           info,
                                          int *pBufferSizeInBytes)
  {
    return cusparseCcsrilu02_bufferSize(handle,
                                        m,
                                        nnz,
                                        descrA,
                                        csrValA,
                                        csrRowPtrA,
                                        csrColIndA,
                                        info,
                                        pBufferSizeInBytes);
  }

  template <>
  cusparseStatus_t
  cusparseXcsrilu02_bufferSize<cuDoubleComplex>(cusparseHandle_t         handle,
                                                int                      m,
                                                int                      nnz,
                                                const cusparseMatDescr_t descrA,
                                                cuDoubleComplex *csrValA,
                                                const int *      csrRowPtrA,
                                                const int *      csrColIndA,
                                                csrilu02Info_t   info,
                                                int *pBufferSizeInBytes)
  {
    return cusparseZcsrilu02_bufferSize(handle,
                                        m,
                                        nnz,
                                        descrA,
                                        csrValA,
                                        csrRowPtrA,
                                        csrColIndA,
                                        info,
                                        pBufferSizeInBytes);
  }


  /**
   * Template wrapper for cusparse<t>csric02
   * (https://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-csric02).
   * This function performs the solve phase of the computing the
   * incomplete-Cholesky factorization with 0 fill-in and no pivoting.
   */
  template <typename Number>
  cusparseStatus_t
  cusparseXcsric02(cusparseHandle_t         handle,
                   int                      m,
                   int                      nnz,
                   const cusparseMatDescr_t descrA,
                   Number *                 csrValA_valM,
                   const int *              csrRowPtrA,
                   const int *              csrColIndA,
                   csric02Info_t            info,
                   cusparseSolvePolicy_t    policy,
                   void *                   pBuffer)
  {
    AssertThrow(false, ExcNotImplemented());
    return CUSPARSE_STATUS_INVALID_VALUE;
  }

  template <>
  cusparseStatus_t
  cusparseXcsric02<float>(cusparseHandle_t         handle,
                          int                      m,
                          int                      nnz,
                          const cusparseMatDescr_t descrA,
                          float *                  csrValA_valM,
                          const int *              csrRowPtrA,
                          const int *              csrColIndA,
                          csric02Info_t            info,
                          cusparseSolvePolicy_t    policy,
                          void *                   pBuffer)
  {
    return cusparseScsric02(handle,
                            m,
                            nnz,
                            descrA,
                            csrValA_valM,
                            csrRowPtrA,
                            csrColIndA,
                            info,
                            policy,
                            pBuffer);
  }

  template <>
  cusparseStatus_t
  cusparseXcsric02<double>(cusparseHandle_t         handle,
                           int                      m,
                           int                      nnz,
                           const cusparseMatDescr_t descrA,
                           double *                 csrValA_valM,
                           const int *              csrRowPtrA,
                           const int *              csrColIndA,
                           csric02Info_t            info,
                           cusparseSolvePolicy_t    policy,
                           void *                   pBuffer)
  {
    return cusparseDcsric02(handle,
                            m,
                            nnz,
                            descrA,
                            csrValA_valM,
                            csrRowPtrA,
                            csrColIndA,
                            info,
                            policy,
                            pBuffer);
  }

  template <>
  cusparseStatus_t
  cusparseXcsric02<cuComplex>(cusparseHandle_t         handle,
                              int                      m,
                              int                      nnz,
                              const cusparseMatDescr_t descrA,
                              cuComplex *              csrValA_valM,
                              const int *              csrRowPtrA,
                              const int *              csrColIndA,
                              csric02Info_t            info,
                              cusparseSolvePolicy_t    policy,
                              void *                   pBuffer)
  {
    return cusparseCcsric02(handle,
                            m,
                            nnz,
                            descrA,
                            csrValA_valM,
                            csrRowPtrA,
                            csrColIndA,
                            info,
                            policy,
                            pBuffer);
  }

  template <>
  cusparseStatus_t
  cusparseXcsric02<cuDoubleComplex>(cusparseHandle_t         handle,
                                    int                      m,
                                    int                      nnz,
                                    const cusparseMatDescr_t descrA,
                                    cuDoubleComplex *        csrValA_valM,
                                    const int *              csrRowPtrA,
                                    const int *              csrColIndA,
                                    csric02Info_t            info,
                                    cusparseSolvePolicy_t    policy,
                                    void *                   pBuffer)
  {
    return cusparseZcsric02(handle,
                            m,
                            nnz,
                            descrA,
                            csrValA_valM,
                            csrRowPtrA,
                            csrColIndA,
                            info,
                            policy,
                            pBuffer);
  }


  /**
   * Template wrapper for cusparse<t>csrsv2_solve
   *(https://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-csrsv2_solve).
   * This function performs the solve phase of csrsv2, a new sparse triangular
   *linear system op(A)*y = alpha*x.
   */
  template <typename Number>
  cusparseStatus_t
  cusparseXcsrsv2_solve(cusparseHandle_t         handle,
                        cusparseOperation_t      transA,
                        int                      m,
                        int                      nnz,
                        const Number *           alpha,
                        const cusparseMatDescr_t descra,
                        const Number *           csrValA,
                        const int *              csrRowPtrA,
                        const int *              csrColIndA,
                        csrsv2Info_t             info,
                        const Number *           x,
                        Number *                 y,
                        cusparseSolvePolicy_t    policy,
                        void *                   pBuffer)
  {
    AssertThrow(false, ExcNotImplemented());
    return CUSPARSE_STATUS_INVALID_VALUE;
  }

  template <>
  cusparseStatus_t
  cusparseXcsrsv2_solve<float>(cusparseHandle_t         handle,
                               cusparseOperation_t      transA,
                               int                      m,
                               int                      nnz,
                               const float *            alpha,
                               const cusparseMatDescr_t descra,
                               const float *            csrValA,
                               const int *              csrRowPtrA,
                               const int *              csrColIndA,
                               csrsv2Info_t             info,
                               const float *            x,
                               float *                  y,
                               cusparseSolvePolicy_t    policy,
                               void *                   pBuffer)
  {
    return cusparseScsrsv2_solve(handle,
                                 transA,
                                 m,
                                 nnz,
                                 alpha,
                                 descra,
                                 csrValA,
                                 csrRowPtrA,
                                 csrColIndA,
                                 info,
                                 x,
                                 y,
                                 policy,
                                 pBuffer);
  }

  template <>
  cusparseStatus_t
  cusparseXcsrsv2_solve<double>(cusparseHandle_t         handle,
                                cusparseOperation_t      transA,
                                int                      m,
                                int                      nnz,
                                const double *           alpha,
                                const cusparseMatDescr_t descra,
                                const double *           csrValA,
                                const int *              csrRowPtrA,
                                const int *              csrColIndA,
                                csrsv2Info_t             info,
                                const double *           x,
                                double *                 y,
                                cusparseSolvePolicy_t    policy,
                                void *                   pBuffer)
  {
    return cusparseDcsrsv2_solve(handle,
                                 transA,
                                 m,
                                 nnz,
                                 alpha,
                                 descra,
                                 csrValA,
                                 csrRowPtrA,
                                 csrColIndA,
                                 info,
                                 x,
                                 y,
                                 policy,
                                 pBuffer);
  }

  template <>
  cusparseStatus_t
  cusparseXcsrsv2_solve<cuComplex>(cusparseHandle_t         handle,
                                   cusparseOperation_t      transA,
                                   int                      m,
                                   int                      nnz,
                                   const cuComplex *        alpha,
                                   const cusparseMatDescr_t descra,
                                   const cuComplex *        csrValA,
                                   const int *              csrRowPtrA,
                                   const int *              csrColIndA,
                                   csrsv2Info_t             info,
                                   const cuComplex *        x,
                                   cuComplex *              y,
                                   cusparseSolvePolicy_t    policy,
                                   void *                   pBuffer)
  {
    return cusparseCcsrsv2_solve(handle,
                                 transA,
                                 m,
                                 nnz,
                                 alpha,
                                 descra,
                                 csrValA,
                                 csrRowPtrA,
                                 csrColIndA,
                                 info,
                                 x,
                                 y,
                                 policy,
                                 pBuffer);
  }

  template <>
  cusparseStatus_t
  cusparseXcsrsv2_solve<cuDoubleComplex>(cusparseHandle_t         handle,
                                         cusparseOperation_t      transA,
                                         int                      m,
                                         int                      nnz,
                                         const cuDoubleComplex *  alpha,
                                         const cusparseMatDescr_t descra,
                                         const cuDoubleComplex *  csrValA,
                                         const int *              csrRowPtrA,
                                         const int *              csrColIndA,
                                         csrsv2Info_t             info,
                                         const cuDoubleComplex *  x,
                                         cuDoubleComplex *        y,
                                         cusparseSolvePolicy_t    policy,
                                         void *                   pBuffer)
  {
    return cusparseZcsrsv2_solve(handle,
                                 transA,
                                 m,
                                 nnz,
                                 alpha,
                                 descra,
                                 csrValA,
                                 csrRowPtrA,
                                 csrColIndA,
                                 info,
                                 x,
                                 y,
                                 policy,
                                 pBuffer);
  }


  /**
   * Template wrapper for cusparse<t>csrsv2_analysis
   * (https://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-csrsv2_analysis).
   * This function performs the analysis phase of csrsv2, a new sparse
   * triangular linear system op(A)*y = alpha*x.
   */
  template <typename Number>
  cusparseStatus_t
  cusparseXcsrsv2_analysis(cusparseHandle_t         handle,
                           cusparseOperation_t      transA,
                           int                      m,
                           int                      nnz,
                           const cusparseMatDescr_t descrA,
                           const Number *           csrValA,
                           const int *              csrRowPtrA,
                           const int *              csrColIndA,
                           csrsv2Info_t             info,
                           cusparseSolvePolicy_t    policy,
                           void *                   pBuffer)
  {
    AssertThrow(false, ExcNotImplemented());
    return CUSPARSE_STATUS_INVALID_VALUE;
  }

  template <>
  cusparseStatus_t
  cusparseXcsrsv2_analysis<float>(cusparseHandle_t         handle,
                                  cusparseOperation_t      transA,
                                  int                      m,
                                  int                      nnz,
                                  const cusparseMatDescr_t descrA,
                                  const float *            csrValA,
                                  const int *              csrRowPtrA,
                                  const int *              csrColIndA,
                                  csrsv2Info_t             info,
                                  cusparseSolvePolicy_t    policy,
                                  void *                   pBuffer)
  {
    return cusparseScsrsv2_analysis(handle,
                                    transA,
                                    m,
                                    nnz,
                                    descrA,
                                    csrValA,
                                    csrRowPtrA,
                                    csrColIndA,
                                    info,
                                    policy,
                                    pBuffer);
  }

  template <>
  cusparseStatus_t
  cusparseXcsrsv2_analysis<double>(cusparseHandle_t         handle,
                                   cusparseOperation_t      transA,
                                   int                      m,
                                   int                      nnz,
                                   const cusparseMatDescr_t descrA,
                                   const double *           csrValA,
                                   const int *              csrRowPtrA,
                                   const int *              csrColIndA,
                                   csrsv2Info_t             info,
                                   cusparseSolvePolicy_t    policy,
                                   void *                   pBuffer)
  {
    return cusparseDcsrsv2_analysis(handle,
                                    transA,
                                    m,
                                    nnz,
                                    descrA,
                                    csrValA,
                                    csrRowPtrA,
                                    csrColIndA,
                                    info,
                                    policy,
                                    pBuffer);
  }

  template <>
  cusparseStatus_t
  cusparseXcsrsv2_analysis<cuComplex>(cusparseHandle_t         handle,
                                      cusparseOperation_t      transA,
                                      int                      m,
                                      int                      nnz,
                                      const cusparseMatDescr_t descrA,
                                      const cuComplex *        csrValA,
                                      const int *              csrRowPtrA,
                                      const int *              csrColIndA,
                                      csrsv2Info_t             info,
                                      cusparseSolvePolicy_t    policy,
                                      void *                   pBuffer)
  {
    return cusparseCcsrsv2_analysis(handle,
                                    transA,
                                    m,
                                    nnz,
                                    descrA,
                                    csrValA,
                                    csrRowPtrA,
                                    csrColIndA,
                                    info,
                                    policy,
                                    pBuffer);
  }

  template <>
  cusparseStatus_t
  cusparseXcsrsv2_analysis<cuDoubleComplex>(cusparseHandle_t         handle,
                                            cusparseOperation_t      transA,
                                            int                      m,
                                            int                      nnz,
                                            const cusparseMatDescr_t descrA,
                                            const cuDoubleComplex *  csrValA,
                                            const int *              csrRowPtrA,
                                            const int *              csrColIndA,
                                            csrsv2Info_t             info,
                                            cusparseSolvePolicy_t    policy,
                                            void *                   pBuffer)
  {
    return cusparseZcsrsv2_analysis(handle,
                                    transA,
                                    m,
                                    nnz,
                                    descrA,
                                    csrValA,
                                    csrRowPtrA,
                                    csrColIndA,
                                    info,
                                    policy,
                                    pBuffer);
  }



  /**
   * Template wrapper for cusparse<t>csric02_analysis
   * (https://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-csric02_analysis).
   * This function performs the analysis phase of the incomplete-Cholesky
   * factorization with 0 fill-in and no pivoting.
   */
  template <typename Number>
  cusparseStatus_t
  cusparseXcsric02_analysis(cusparseHandle_t         handle,
                            int                      m,
                            int                      nnz,
                            const cusparseMatDescr_t descrA,
                            const Number *           csrValA,
                            const int *              csrRowPtrA,
                            const int *              csrColIndA,
                            csric02Info_t            info,
                            cusparseSolvePolicy_t    policy,
                            void *                   pBuffer)
  {
    AssertThrow(false, ExcNotImplemented());
    return CUSPARSE_STATUS_INVALID_VALUE;
  }

  template <>
  cusparseStatus_t
  cusparseXcsric02_analysis<float>(cusparseHandle_t         handle,
                                   int                      m,
                                   int                      nnz,
                                   const cusparseMatDescr_t descrA,
                                   const float *            csrValA,
                                   const int *              csrRowPtrA,
                                   const int *              csrColIndA,
                                   csric02Info_t            info,
                                   cusparseSolvePolicy_t    policy,
                                   void *                   pBuffer)
  {
    return cusparseScsric02_analysis(handle,
                                     m,
                                     nnz,
                                     descrA,
                                     csrValA,
                                     csrRowPtrA,
                                     csrColIndA,
                                     info,
                                     policy,
                                     pBuffer);
  }

  template <>
  cusparseStatus_t
  cusparseXcsric02_analysis<double>(cusparseHandle_t         handle,
                                    int                      m,
                                    int                      nnz,
                                    const cusparseMatDescr_t descrA,
                                    const double *           csrValA,
                                    const int *              csrRowPtrA,
                                    const int *              csrColIndA,
                                    csric02Info_t            info,
                                    cusparseSolvePolicy_t    policy,
                                    void *                   pBuffer)
  {
    return cusparseDcsric02_analysis(handle,
                                     m,
                                     nnz,
                                     descrA,
                                     csrValA,
                                     csrRowPtrA,
                                     csrColIndA,
                                     info,
                                     policy,
                                     pBuffer);
  }

  template <>
  cusparseStatus_t
  cusparseXcsric02_analysis<cuComplex>(cusparseHandle_t         handle,
                                       int                      m,
                                       int                      nnz,
                                       const cusparseMatDescr_t descrA,
                                       const cuComplex *        csrValA,
                                       const int *              csrRowPtrA,
                                       const int *              csrColIndA,
                                       csric02Info_t            info,
                                       cusparseSolvePolicy_t    policy,
                                       void *                   pBuffer)
  {
    return cusparseCcsric02_analysis(handle,
                                     m,
                                     nnz,
                                     descrA,
                                     csrValA,
                                     csrRowPtrA,
                                     csrColIndA,
                                     info,
                                     policy,
                                     pBuffer);
  }

  template <>
  cusparseStatus_t
  cusparseXcsric02_analysis<cuDoubleComplex>(cusparseHandle_t         handle,
                                             int                      m,
                                             int                      nnz,
                                             const cusparseMatDescr_t descrA,
                                             const cuDoubleComplex *  csrValA,
                                             const int *           csrRowPtrA,
                                             const int *           csrColIndA,
                                             csric02Info_t         info,
                                             cusparseSolvePolicy_t policy,
                                             void *                pBuffer)
  {
    return cusparseZcsric02_analysis(handle,
                                     m,
                                     nnz,
                                     descrA,
                                     csrValA,
                                     csrRowPtrA,
                                     csrColIndA,
                                     info,
                                     policy,
                                     pBuffer);
  }


  /**
   * Template wrapper for cusparse<t>csrsv2_bufferSize
   * (https://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-csrsv2_bufferSize).
   * This function returns the size of the buffer used in csrsv2, a new sparse
   * triangular linear system op(A)*y = alpha*x.
   */
  template <typename Number>
  cusparseStatus_t
  cusparseXcsrsv2_bufferSize(cusparseHandle_t         handle,
                             cusparseOperation_t      transA,
                             int                      m,
                             int                      nnz,
                             const cusparseMatDescr_t descrA,
                             Number *                 csrValA,
                             const int *              csrRowPtrA,
                             const int *              csrColIndA,
                             csrsv2Info_t             info,
                             int *                    pBufferSizeInBytes)
  {
    AssertThrow(false, ExcNotImplemented());
    return CUSPARSE_STATUS_INVALID_VALUE;
  }

  template <>
  cusparseStatus_t
  cusparseXcsrsv2_bufferSize<float>(cusparseHandle_t         handle,
                                    cusparseOperation_t      transA,
                                    int                      m,
                                    int                      nnz,
                                    const cusparseMatDescr_t descrA,
                                    float *                  csrValA,
                                    const int *              csrRowPtrA,
                                    const int *              csrColIndA,
                                    csrsv2Info_t             info,
                                    int *                    pBufferSizeInBytes)
  {
    return cusparseScsrsv2_bufferSize(handle,
                                      transA,
                                      m,
                                      nnz,
                                      descrA,
                                      csrValA,
                                      csrRowPtrA,
                                      csrColIndA,
                                      info,
                                      pBufferSizeInBytes);
  }

  template <>
  cusparseStatus_t
  cusparseXcsrsv2_bufferSize<double>(cusparseHandle_t         handle,
                                     cusparseOperation_t      transA,
                                     int                      m,
                                     int                      nnz,
                                     const cusparseMatDescr_t descrA,
                                     double *                 csrValA,
                                     const int *              csrRowPtrA,
                                     const int *              csrColIndA,
                                     csrsv2Info_t             info,
                                     int *pBufferSizeInBytes)
  {
    return cusparseDcsrsv2_bufferSize(handle,
                                      transA,
                                      m,
                                      nnz,
                                      descrA,
                                      csrValA,
                                      csrRowPtrA,
                                      csrColIndA,
                                      info,
                                      pBufferSizeInBytes);
  }

  template <>
  cusparseStatus_t
  cusparseXcsrsv2_bufferSize<cuComplex>(cusparseHandle_t         handle,
                                        cusparseOperation_t      transA,
                                        int                      m,
                                        int                      nnz,
                                        const cusparseMatDescr_t descrA,
                                        cuComplex *              csrValA,
                                        const int *              csrRowPtrA,
                                        const int *              csrColIndA,
                                        csrsv2Info_t             info,
                                        int *pBufferSizeInBytes)
  {
    return cusparseCcsrsv2_bufferSize(handle,
                                      transA,
                                      m,
                                      nnz,
                                      descrA,
                                      csrValA,
                                      csrRowPtrA,
                                      csrColIndA,
                                      info,
                                      pBufferSizeInBytes);
  }

  template <>
  cusparseStatus_t
  cusparseXcsrsv2_bufferSize<cuDoubleComplex>(cusparseHandle_t         handle,
                                              cusparseOperation_t      transA,
                                              int                      m,
                                              int                      nnz,
                                              const cusparseMatDescr_t descrA,
                                              cuDoubleComplex *        csrValA,
                                              const int *  csrRowPtrA,
                                              const int *  csrColIndA,
                                              csrsv2Info_t info,
                                              int *        pBufferSizeInBytes)
  {
    return cusparseZcsrsv2_bufferSize(handle,
                                      transA,
                                      m,
                                      nnz,
                                      descrA,
                                      csrValA,
                                      csrRowPtrA,
                                      csrColIndA,
                                      info,
                                      pBufferSizeInBytes);
  }



  /**
   * Template wrapper for cusparse<t>csric02_bufferSize
   * (https://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-csric02_bufferSize).
   *This function returns size of buffer used in computing the
   *incomplete-Cholesky factorization with 0 fill-in and no pivoting.
   */
  template <typename Number>
  cusparseStatus_t
  cusparseXcsric02_bufferSize(cusparseHandle_t         handle,
                              int                      m,
                              int                      nnz,
                              const cusparseMatDescr_t descrA,
                              Number *                 csrValA,
                              const int *              csrRowPtrA,
                              const int *              csrColIndA,
                              csric02Info_t            info,
                              int *                    pBufferSizeInBytes)
  {
    AssertThrow(false, ExcNotImplemented());
    return CUSPARSE_STATUS_INVALID_VALUE;
  }

  template <>
  cusparseStatus_t
  cusparseXcsric02_bufferSize<float>(cusparseHandle_t         handle,
                                     int                      m,
                                     int                      nnz,
                                     const cusparseMatDescr_t descrA,
                                     float *                  csrValA,
                                     const int *              csrRowPtrA,
                                     const int *              csrColIndA,
                                     csric02Info_t            info,
                                     int *pBufferSizeInBytes)
  {
    return cusparseScsric02_bufferSize(handle,
                                       m,
                                       nnz,
                                       descrA,
                                       csrValA,
                                       csrRowPtrA,
                                       csrColIndA,
                                       info,
                                       pBufferSizeInBytes);
  }

  template <>
  cusparseStatus_t
  cusparseXcsric02_bufferSize<double>(cusparseHandle_t         handle,
                                      int                      m,
                                      int                      nnz,
                                      const cusparseMatDescr_t descrA,
                                      double *                 csrValA,
                                      const int *              csrRowPtrA,
                                      const int *              csrColIndA,
                                      csric02Info_t            info,
                                      int *pBufferSizeInBytes)
  {
    return cusparseDcsric02_bufferSize(handle,
                                       m,
                                       nnz,
                                       descrA,
                                       csrValA,
                                       csrRowPtrA,
                                       csrColIndA,
                                       info,
                                       pBufferSizeInBytes);
  }

  template <>
  cusparseStatus_t
  cusparseXcsric02_bufferSize<cuComplex>(cusparseHandle_t         handle,
                                         int                      m,
                                         int                      nnz,
                                         const cusparseMatDescr_t descrA,
                                         cuComplex *              csrValA,
                                         const int *              csrRowPtrA,
                                         const int *              csrColIndA,
                                         csric02Info_t            info,
                                         int *pBufferSizeInBytes)
  {
    return cusparseCcsric02_bufferSize(handle,
                                       m,
                                       nnz,
                                       descrA,
                                       csrValA,
                                       csrRowPtrA,
                                       csrColIndA,
                                       info,
                                       pBufferSizeInBytes);
  }

  template <>
  cusparseStatus_t
  cusparseXcsric02_bufferSize<cuDoubleComplex>(cusparseHandle_t         handle,
                                               int                      m,
                                               int                      nnz,
                                               const cusparseMatDescr_t descrA,
                                               cuDoubleComplex *        csrValA,
                                               const int *   csrRowPtrA,
                                               const int *   csrColIndA,
                                               csric02Info_t info,
                                               int *         pBufferSizeInBytes)
  {
    return cusparseZcsric02_bufferSize(handle,
                                       m,
                                       nnz,
                                       descrA,
                                       csrValA,
                                       csrRowPtrA,
                                       csrColIndA,
                                       info,
                                       pBufferSizeInBytes);
  }
  /**
   * @}
   */
} // namespace CUDAWrappers

DEAL_II_NAMESPACE_CLOSE

namespace
{
  template <typename Number>
  void
  delete_device_vector(Number *device_ptr) noexcept
  {
    const cudaError_t error_code = cudaFree(device_ptr);
    (void)error_code;
    AssertNothrow(error_code == cudaSuccess,
                  dealii::ExcCudaError(cudaGetErrorString(error_code)));
  }
  template <typename Number>
  Number *
  allocate_device_vector(const std::size_t size)
  {
    Number *device_ptr;
    Utilities::CUDA::malloc(device_ptr, size);
    return device_ptr;
  }
} // namespace

namespace dealii
{
  namespace CUDAWrappers
  {
    /**
     * This class implements an incomplete Cholesky factorization (IC)
     * preconditioner for @em symmetric CUDAWrappers::SparseMatrix matrices.
     *
     * The implementation closely follows the one documented in the cuSPARSE
     * documentation
     * (https://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-csrilu02).
     *
     * @note Instantiations for this template are provided for <tt>@<float@> and
     * @<double@></tt>.
     *
     * @ingroup Preconditioners CUDAWrappers
     * @author Daniel Arndt
     * @date 2018
     */
    template <typename Number>
    class PreconditionILU
    {
    public:
      /**
       * Declare the type for container size.
       */
      using size_type = int;

      /**
       * Standardized data struct to pipe additional flags to the
       * preconditioner.
       */
      struct AdditionalData
      {
        /**
         * Constructor. cuSPARSE allows to compute and use level information.
         * According to the documentation it is this might improve performance.
         * It is suggested to try both options.
         */
        AdditionalData(bool use_level_analysis = true);

        /**
         * Flag that determines if level informations are used when creating and
         * applying the preconditioner. See the documentation for
         * cusparseSolvePolicy_t at
         * https://docs.nvidia.com/cuda/cusparse/index.html#cusparsesolvepolicy_t
         * for more information.
         */
        bool use_level_analysis;
      };

      /**
       * Constructor.
       */
      PreconditionILU(const Utilities::CUDA::Handle &handle);

      /**
       * The copy constructor is deleted.
       */
      PreconditionILU(const PreconditionILU<Number> &) = delete;

      /**
       * The copy assignment operator is deleted.
       */
      PreconditionILU &
      operator=(const PreconditionILU<Number> &) = delete;

      /**
       * Destructor. Free all resources that were initialized in this class.
       */
      ~PreconditionILU();

      /**
       * Initialize this object. In particular, the given matrix is copied to be
       * modified in-place. For the underlying sparsity pattern pointers are
       * stored. Specifically, this means
       * that the current object can only be used reliably as long as @p matrix is valid
       * and has not been changed since calling this function.
       *
       * The @p additional_data determines if level information are used.
       */
      void
      initialize(const SparseMatrix<Number> &matrix,
                 const AdditionalData &additional_data = AdditionalData());

      /**
       * Apply the preconditioner.
       */
      void
      vmult(LinearAlgebra::CUDAWrappers::Vector<Number> &      dst,
            const LinearAlgebra::CUDAWrappers::Vector<Number> &src) const;

      /**
       * Apply the preconditioner. Since the preconditioner is symmetric, this
       * is the same as vmult().
       */
      void
      Tvmult(LinearAlgebra::CUDAWrappers::Vector<Number> &      dst,
             const LinearAlgebra::CUDAWrappers::Vector<Number> &src) const;

      /**
       *  Return the dimension of the codomain (or range) space. Note that the
       * matrix is square and has dimension $m \times m$.
       *
       * @note This function should only be called if the preconditioner has been
       * initialized.
       */
      size_type
      m() const;

      /**
       *  Return the dimension of the codomain (or range) space. Note that the
       * matrix is square and has dimension $m \times m$.
       *
       * @note This function should only be called if the preconditioner has been
       * initialized.
       */
      size_type
      n() const;

    private:
      /**
       * cuSPARSE handle used to call cuSPARSE functions.
       */
      cusparseHandle_t cusparse_handle;

      /**
       * cuSPARSE description of the sparse matrix $M=LU$.
       */
      cusparseMatDescr_t descr_M;

      /**
       * cuSPARSE description of the lower triangular matrix $L$.
       */
      cusparseMatDescr_t descr_L;

      /**
       * cuSPARSE description of the upper triangular matrix $U$.
       */
      cusparseMatDescr_t descr_U;

      /**
       * Solve and analysis structure for $M=LL^T$.
       */
      csrilu02Info_t info_M;

      /**
       * Solve and analysis structure for the lower triangular matrix $L$.
       */
      csrsv2Info_t info_L;

      /**
       * Solve and analysis structure for the upper triangular matrix $U$.
       */
      csrsv2Info_t info_U;

      /**
       * Pointer to the values (on the device) of the computed preconditioning
       * matrix.
       */
      std::unique_ptr<Number[], void (*)(Number *)> P_val_dev;

      /**
       * Pointer to the row pointer (on the device) of the sparse matrix this
       * object was initialized with.
       */
      const int *P_row_ptr_dev;

      /**
       * Pointer to the column indices (on the device) of the sparse matrix this
       * object was initialized with.
       */
      const int *P_column_index_dev;

      /**
       * Pointer to the value (on the device) for a temporary (helper) vector
       * used in vmult().
       */
      std::unique_ptr<Number[], void (*)(Number *)> tmp_dev;

      /**
       *
       */
      std::unique_ptr<void, void (*)(void *)> buffer_dev;

      /**
       * Determine if level information should be generated for the lower
       * triangular matrix $L$. This value can be modified through an
       * AdditionalData object.
       */
      cusparseSolvePolicy_t policy_L;

      /**
       * Determine if level information should be generated for the upper
       * triangular matrix $L^T$. This value can be modified through an
       * AdditionalData object.
       */
      cusparseSolvePolicy_t policy_U;

      /**
       * Determine if level information should be generated for $M=LL^T$. This
       * value can be modified through an AdditionalData object.
       */
      cusparseSolvePolicy_t policy_M;

      /**
       * The number of rows is the same as for the matrix this object has been
       * initialized with.
       */
      int n_rows;

      /**
       * The number of non-zero elements is the same as for the matrix this
       * object has been initialized with.
       */
      int n_nonzero_elements;
    };

    template <typename Number>
    PreconditionILU<Number>::AdditionalData::AdditionalData(
      bool use_level_analysis_)
      : use_level_analysis(use_level_analysis_)
    {}



    template <typename Number>
    PreconditionILU<Number>::PreconditionILU(
      const Utilities::CUDA::Handle &handle)
      : cusparse_handle(handle.cusparse_handle)
      , P_val_dev(nullptr, delete_device_vector<Number>)
      , P_row_ptr_dev(nullptr)
      , P_column_index_dev(nullptr)
      , tmp_dev(nullptr, delete_device_vector<Number>)
      , buffer_dev(nullptr, delete_device_vector<void>)
      , policy_L(CUSPARSE_SOLVE_POLICY_USE_LEVEL)
      , policy_U(CUSPARSE_SOLVE_POLICY_USE_LEVEL)
      , policy_M(CUSPARSE_SOLVE_POLICY_USE_LEVEL)
      , n_rows(0)
      , n_nonzero_elements(0)
    {
      cusparseStatus_t status;
      // step 1: create a descriptor which contains
      // - matrix M is base-0
      // - matrix L is base-0
      // - matrix L is lower triangular
      // - matrix L has unit diagonal
      // - matrix U is base-0
      // - matrix U is upper triangular
      // - matrix U has non-unit diagonal
      status = cusparseCreateMatDescr(&descr_M);
      AssertCusparse(status);
      status = cusparseSetMatIndexBase(descr_M, CUSPARSE_INDEX_BASE_ZERO);
      AssertCusparse(status);
      status = cusparseSetMatType(descr_M, CUSPARSE_MATRIX_TYPE_GENERAL);
      AssertCusparse(status);

      status = cusparseCreateMatDescr(&descr_L);
      AssertCusparse(status);
      status = cusparseSetMatIndexBase(descr_L, CUSPARSE_INDEX_BASE_ZERO);
      AssertCusparse(status);
      status = cusparseSetMatType(descr_L, CUSPARSE_MATRIX_TYPE_GENERAL);
      AssertCusparse(status);
      status = cusparseSetMatFillMode(descr_L, CUSPARSE_FILL_MODE_LOWER);
      AssertCusparse(status);
      status = cusparseSetMatDiagType(descr_L, CUSPARSE_DIAG_TYPE_UNIT);
      AssertCusparse(status);

      status = cusparseCreateMatDescr(&descr_U);
      AssertCusparse(status);
      status = cusparseSetMatIndexBase(descr_U, CUSPARSE_INDEX_BASE_ZERO);
      AssertCusparse(status);
      status = cusparseSetMatType(descr_U, CUSPARSE_MATRIX_TYPE_GENERAL);
      AssertCusparse(status);
      status = cusparseSetMatFillMode(descr_U, CUSPARSE_FILL_MODE_UPPER);
      AssertCusparse(status);
      status = cusparseSetMatDiagType(descr_U, CUSPARSE_DIAG_TYPE_NON_UNIT);
      AssertCusparse(status);

      // step 2: create a empty info structure
      // we need one info for csrilu02 and two info's for csrsv2
      status = cusparseCreateCsrilu02Info(&info_M);
      AssertCusparse(status);
      status = cusparseCreateCsrsv2Info(&info_L);
      AssertCusparse(status);
      status = cusparseCreateCsrsv2Info(&info_U);
      AssertCusparse(status);
    }

    template <typename Number>
    PreconditionILU<Number>::~PreconditionILU()
    {
      // step 8: free resources
      cusparseStatus_t status = cusparseDestroyMatDescr(descr_M);
      AssertNothrowCusparse(status);

      status = cusparseDestroyMatDescr(descr_L);
      AssertNothrowCusparse(status);

      status = cusparseDestroyMatDescr(descr_U);
      AssertNothrowCusparse(status);

      status = cusparseDestroyCsrilu02Info(info_M);
      AssertNothrowCusparse(status);

      status = cusparseDestroyCsrsv2Info(info_L);
      AssertNothrowCusparse(status);

      status = cusparseDestroyCsrsv2Info(info_U);
      AssertNothrowCusparse(status);
    }



    template <typename Number>
    void
    PreconditionILU<Number>::initialize(const SparseMatrix<Number> &A,
                                        const AdditionalData &additional_data)
    {
      if (additional_data.use_level_analysis)
        {
          policy_L = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
          policy_U = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
          policy_M = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
        }
      else
        {
          policy_L = CUSPARSE_SOLVE_POLICY_NO_LEVEL;
          policy_U = CUSPARSE_SOLVE_POLICY_NO_LEVEL;
          policy_M = CUSPARSE_SOLVE_POLICY_NO_LEVEL;
        }


      n_rows             = A.m();
      n_nonzero_elements = A.n_nonzero_elements();
      AssertDimension(A.m(), A.n());

      const auto          cusparse_matrix = A.get_cusparse_matrix();
      const Number *const A_val_dev       = std::get<0>(cusparse_matrix);

      // create a copy of the matrix entries
      P_val_dev.reset(allocate_device_vector<Number>(n_nonzero_elements));
      cudaError_t cuda_status            = cudaMemcpy(P_val_dev.get(),
                                           A_val_dev,
                                           n_nonzero_elements * sizeof(Number),
                                           cudaMemcpyDeviceToDevice);
      P_column_index_dev                 = std::get<1>(cusparse_matrix);
      P_row_ptr_dev                      = std::get<2>(cusparse_matrix);
      const cusparseMatDescr_t mat_descr = std::get<3>(cusparse_matrix);

      // initializa an internal buffer we need later on
      tmp_dev.reset(allocate_device_vector<Number>(n_rows));

      // step 3: query how much memory used in csrilu02 and csrsv2, and allocate
      // the buffer
      int              BufferSize_M;
      cusparseStatus_t status = cusparseXcsrilu02_bufferSize(cusparse_handle,
                                                             n_rows,
                                                             n_nonzero_elements,
                                                             descr_M,
                                                             P_val_dev.get(),
                                                             P_row_ptr_dev,
                                                             P_column_index_dev,
                                                             info_M,
                                                             &BufferSize_M);
      AssertCusparse(status);

      int BufferSize_L;
      status = cusparseXcsrsv2_bufferSize(cusparse_handle,
                                          CUSPARSE_OPERATION_NON_TRANSPOSE,
                                          n_rows,
                                          n_nonzero_elements,
                                          descr_L,
                                          P_val_dev.get(),
                                          P_row_ptr_dev,
                                          P_column_index_dev,
                                          info_L,
                                          &BufferSize_L);
      AssertCusparse(status);

      int BufferSize_U;
      status = cusparseXcsrsv2_bufferSize(cusparse_handle,
                                          CUSPARSE_OPERATION_NON_TRANSPOSE,
                                          n_rows,
                                          n_nonzero_elements,
                                          descr_U,
                                          P_val_dev.get(),
                                          P_row_ptr_dev,
                                          P_column_index_dev,
                                          info_U,
                                          &BufferSize_U);
      AssertCusparse(status);

      const int BufferSize =
        std::max(BufferSize_M, std::max(BufferSize_L, BufferSize_U));
      // workaround: since allocate_device_vector needs a type, we pass char
      // which is required to have size 1.
      buffer_dev.reset(static_cast<void *>(
        allocate_device_vector<char>(BufferSize / sizeof(char))));

      // step 4: perform analysis of incomplete Cholesky on M
      //         perform analysis of triangular solve on L
      //         perform analysis of triangular solve on U
      // The lower(upper) triangular part of M has the same sparsity pattern as
      // L(U), we can do analysis of csrilu0 and csrsv2 simultaneously.

      status = cusparseXcsrilu02_analysis(cusparse_handle,
                                          n_rows,
                                          n_nonzero_elements,
                                          descr_M,
                                          P_val_dev.get(),
                                          P_row_ptr_dev,
                                          P_column_index_dev,
                                          info_M,
                                          policy_M,
                                          buffer_dev.get());
      AssertCusparse(status);

      int structural_zero;
      status =
        cusparseXcsrilu02_zeroPivot(cusparse_handle, info_M, &structural_zero);
      AssertCusparse(status);

      status = cusparseXcsrsv2_analysis(cusparse_handle,
                                        CUSPARSE_OPERATION_NON_TRANSPOSE,
                                        n_rows,
                                        n_nonzero_elements,
                                        descr_L,
                                        P_val_dev.get(),
                                        P_row_ptr_dev,
                                        P_column_index_dev,
                                        info_L,
                                        policy_L,
                                        buffer_dev.get());
      AssertCusparse(status);

      status = cusparseXcsrsv2_analysis(cusparse_handle,
                                        CUSPARSE_OPERATION_NON_TRANSPOSE,
                                        n_rows,
                                        n_nonzero_elements,
                                        descr_U,
                                        P_val_dev.get(),
                                        P_row_ptr_dev,
                                        P_column_index_dev,
                                        info_U,
                                        policy_U,
                                        buffer_dev.get());

      // step 5: M = L * U
      status = cusparseXcsrilu02(cusparse_handle,
                                 n_rows,
                                 n_nonzero_elements,
                                 descr_M,
                                 P_val_dev.get(),
                                 P_row_ptr_dev,
                                 P_column_index_dev,
                                 info_M,
                                 policy_M,
                                 buffer_dev.get());
      AssertCusparse(status);

      int numerical_zero;
      status =
        cusparseXcsrilu02_zeroPivot(cusparse_handle, info_M, &numerical_zero);
      AssertCusparse(status);
    }



    template <typename Number>
    void
    PreconditionILU<Number>::vmult(
      LinearAlgebra::CUDAWrappers::Vector<Number> &      dst,
      const LinearAlgebra::CUDAWrappers::Vector<Number> &src) const
    {
      Assert(P_val_dev != nullptr, ExcNotInitialized());
      Assert(P_row_ptr_dev != nullptr, ExcNotInitialized());
      Assert(P_column_index_dev != nullptr, ExcNotInitialized());
      AssertDimension(dst.size(), static_cast<unsigned int>(n_rows));
      AssertDimension(src.size(), static_cast<unsigned int>(n_rows));
      Assert(tmp_dev != nullptr, ExcInternalError());

      const Number *const src_dev = src.get_values();
      Number *const       dst_dev = dst.get_values();

      // step 6: solve L*z = alpha*x
      const Number     alpha = 1.;
      cusparseStatus_t status =
        cusparseXcsrsv2_solve(cusparse_handle,
                              CUSPARSE_OPERATION_NON_TRANSPOSE,
                              n_rows,
                              n_nonzero_elements,
                              &alpha,
                              descr_L,
                              P_val_dev.get(),
                              P_row_ptr_dev,
                              P_column_index_dev,
                              info_L,
                              src_dev,
                              tmp_dev.get(),
                              policy_L,
                              buffer_dev.get());
      AssertCusparse(status);

      // step 7: solve U*y = alpha*z
      status = cusparseXcsrsv2_solve(cusparse_handle,
                                     CUSPARSE_OPERATION_NON_TRANSPOSE,
                                     n_rows,
                                     n_nonzero_elements,
                                     &alpha,
                                     descr_U,
                                     P_val_dev.get(),
                                     P_row_ptr_dev,
                                     P_column_index_dev,
                                     info_U,
                                     tmp_dev.get(),
                                     dst_dev,
                                     policy_U,
                                     buffer_dev.get());
      AssertCusparse(status);
    }



    template <typename Number>
    void
    PreconditionILU<Number>::Tvmult(
      LinearAlgebra::CUDAWrappers::Vector<Number> &      dst,
      const LinearAlgebra::CUDAWrappers::Vector<Number> &src) const
    {
      // the constructed preconditioner is symmetric
      vmult(dst, src);
    }

    template <typename Number>
    PreconditionILU<Number>::size_type
    PreconditionILU<Number>::m() const
    {
      return n_rows;
    }


    template <typename Number>
    PreconditionILU<Number>::size_type
    PreconditionILU<Number>::n() const
    {
      return n_rows;
    }



    template <typename Number>
    void
    apply_preconditioner(const SparseMatrix<Number> &A,
                         const cusparseHandle_t      cusparse_handle,
                         LinearAlgebra::CUDAWrappers::Vector<Number> &      dst,
                         const LinearAlgebra::CUDAWrappers::Vector<Number> &src)
    {
      const Number *const    src_dev = src.get_values();
      Number *               dst_dev = dst.get_values();
      const cusparseHandle_t handle  = cusparse_handle;

      const auto       cusparse_matrix    = A.get_cusparse_matrix();
      Number *         A_val_dev          = std::get<0>(cusparse_matrix);
      const int *const A_row_ptr_dev      = std::get<2>(cusparse_matrix);
      const int *const A_column_index_dev = std::get<1>(cusparse_matrix);
      const cusparseMatDescr_t mat_descr  = std::get<3>(cusparse_matrix);

      const unsigned int n_rows             = A.m();
      const unsigned int n_nonzero_elements = A.n_nonzero_elements();

      AssertDimension(dst.size(), src.size());
      AssertDimension(A.m(), src.size());
      AssertDimension(A.n(), src.size());

      std::unique_ptr<Number[], void (*)(Number *)> tmp_dev(
        allocate_device_vector<Number>(dst.size()),
        delete_device_vector<Number>);

      // Suppose that A is a m x m sparse matrix represented by CSR format,
      // Assumption:
      // - handle is already created by cusparseCreate(),
      // - (A_row_ptr_dev, A_column_index_dev, A_val_dev) is CSR of A on device
      // memory,
      // - src_dev is right hand side vector on device memory,
      // - dst_dev is solution vector on device memory.
      // - tmp_dev is intermediate result on device memory.

      cusparseMatDescr_t          descr_M = mat_descr;
      cusparseMatDescr_t          descr_L = mat_descr;
      cusparseMatDescr_t          descr_U = mat_descr;
      csrilu02Info_t              info_M  = 0;
      csrsv2Info_t                info_L  = 0;
      csrsv2Info_t                info_U  = 0;
      int                         BufferSize_M;
      int                         BufferSize_L;
      int                         BufferSize_U;
      int                         BufferSize;
      void *                      buffer_dev = 0;
      int                         structural_zero;
      int                         numerical_zero;
      const double                alpha    = 1.;
      const cusparseSolvePolicy_t policy_M = CUSPARSE_SOLVE_POLICY_NO_LEVEL;
      const cusparseSolvePolicy_t policy_L = CUSPARSE_SOLVE_POLICY_NO_LEVEL;
      const cusparseSolvePolicy_t policy_U = CUSPARSE_SOLVE_POLICY_USE_LEVEL;

      // step 1: create a descriptor which contains
      // - matrix M is base-0
      // - matrix L is base-0
      // - matrix L is lower triangular
      // - matrix L has unit diagonal
      // - matrix U is base-0
      // - matrix U is upper triangular
      // - matrix U has non-unit diagonal
      cusparseStatus_t status = cusparseCreateMatDescr(&descr_M);
      AssertCusparse(status);
      status = cusparseSetMatIndexBase(descr_M, CUSPARSE_INDEX_BASE_ZERO);
      AssertCusparse(status);
      status = cusparseSetMatType(descr_M, CUSPARSE_MATRIX_TYPE_GENERAL);
      AssertCusparse(status);

      status = cusparseCreateMatDescr(&descr_L);
      AssertCusparse(status);
      status = cusparseSetMatIndexBase(descr_L, CUSPARSE_INDEX_BASE_ZERO);
      AssertCusparse(status);
      status = cusparseSetMatType(descr_L, CUSPARSE_MATRIX_TYPE_GENERAL);
      AssertCusparse(status);
      status = cusparseSetMatFillMode(descr_L, CUSPARSE_FILL_MODE_LOWER);
      AssertCusparse(status);
      status = cusparseSetMatDiagType(descr_L, CUSPARSE_DIAG_TYPE_UNIT);
      AssertCusparse(status);

      status = cusparseCreateMatDescr(&descr_U);
      AssertCusparse(status);
      status = cusparseSetMatIndexBase(descr_U, CUSPARSE_INDEX_BASE_ZERO);
      AssertCusparse(status);
      status = cusparseSetMatType(descr_U, CUSPARSE_MATRIX_TYPE_GENERAL);
      AssertCusparse(status);
      status = cusparseSetMatFillMode(descr_U, CUSPARSE_FILL_MODE_UPPER);
      AssertCusparse(status);
      status = cusparseSetMatDiagType(descr_U, CUSPARSE_DIAG_TYPE_NON_UNIT);
      AssertCusparse(status);

      // step 2: create a empty info structure
      // we need one info for csrilu02 and two info's for csrsv2
      status = cusparseCreateCsrilu02Info(&info_M);
      AssertCusparse(status);
      status = cusparseCreateCsrsv2Info(&info_L);
      AssertCusparse(status);
      status = cusparseCreateCsrsv2Info(&info_U);
      AssertCusparse(status);

      // step 3: query how much memory used in csrilu02 and csrsv2, and allocate
      // the buffer
      status = cusparseXcsrilu02_bufferSize(handle,
                                            n_rows,
                                            n_nonzero_elements,
                                            descr_M,
                                            A_val_dev,
                                            A_row_ptr_dev,
                                            A_column_index_dev,
                                            info_M,
                                            &BufferSize_M);
      AssertCusparse(status);

      status = cusparseXcsrsv2_bufferSize(handle,
                                          CUSPARSE_OPERATION_NON_TRANSPOSE,
                                          n_rows,
                                          n_nonzero_elements,
                                          descr_L,
                                          A_val_dev,
                                          A_row_ptr_dev,
                                          A_column_index_dev,
                                          info_L,
                                          &BufferSize_L);
      AssertCusparse(status);

      status = cusparseXcsrsv2_bufferSize(handle,
                                          CUSPARSE_OPERATION_NON_TRANSPOSE,
                                          n_rows,
                                          n_nonzero_elements,
                                          descr_U,
                                          A_val_dev,
                                          A_row_ptr_dev,
                                          A_column_index_dev,
                                          info_U,
                                          &BufferSize_U);
      AssertCusparse(status);

      BufferSize = max(BufferSize_M, max(BufferSize_L, BufferSize_U));

      // Buffer returned by cudaMalloc is automatically aligned to 128 bytes.
      cudaMalloc((void **)&buffer_dev, BufferSize);

      // step 4: perform analysis of incomplete Cholesky on M
      //         perform analysis of triangular solve on L
      //         perform analysis of triangular solve on U
      // The lower(upper) triangular part of M has the same sparsity pattern as
      // L(U), we can do analysis of csrilu0 and csrsv2 simultaneously.

      status = cusparseXcsrilu02_analysis(handle,
                                          n_rows,
                                          n_nonzero_elements,
                                          descr_M,
                                          A_val_dev,
                                          A_row_ptr_dev,
                                          A_column_index_dev,
                                          info_M,
                                          policy_M,
                                          buffer_dev);
      status = cusparseXcsrilu02_zeroPivot(handle, info_M, &structural_zero);
      AssertCusparse(status);
      if (CUSPARSE_STATUS_ZERO_PIVOT == status)
        {
          printf("A(%d,%d) is missing\n", structural_zero, structural_zero);
        }

      status = cusparseXcsrsv2_analysis(handle,
                                        CUSPARSE_OPERATION_NON_TRANSPOSE,
                                        n_rows,
                                        n_nonzero_elements,
                                        descr_L,
                                        A_val_dev,
                                        A_row_ptr_dev,
                                        A_column_index_dev,
                                        info_L,
                                        policy_L,
                                        buffer_dev);
      AssertCusparse(status);

      status = cusparseXcsrsv2_analysis(handle,
                                        CUSPARSE_OPERATION_NON_TRANSPOSE,
                                        n_rows,
                                        n_nonzero_elements,
                                        descr_U,
                                        A_val_dev,
                                        A_row_ptr_dev,
                                        A_column_index_dev,
                                        info_U,
                                        policy_U,
                                        buffer_dev);
      AssertCusparse(status);

      // step 5: M = L * U
      status = cusparseXcsrilu02(handle,
                                 n_rows,
                                 n_nonzero_elements,
                                 descr_M,
                                 A_val_dev,
                                 A_row_ptr_dev,
                                 A_column_index_dev,
                                 info_M,
                                 policy_M,
                                 buffer_dev);
      status = cusparseXcsrilu02_zeroPivot(handle, info_M, &numerical_zero);
      AssertCusparse(status);
      if (CUSPARSE_STATUS_ZERO_PIVOT == status)
        {
          printf("U(%d,%d) is zero\n", numerical_zero, numerical_zero);
        }

      // step 6: solve L*z = x
      status = cusparseXcsrsv2_solve(handle,
                                     CUSPARSE_OPERATION_NON_TRANSPOSE,
                                     n_rows,
                                     n_nonzero_elements,
                                     &alpha,
                                     descr_L,
                                     A_val_dev,
                                     A_row_ptr_dev,
                                     A_column_index_dev,
                                     info_L,
                                     src_dev,
                                     tmp_dev.get(),
                                     policy_L,
                                     buffer_dev);
      AssertCusparse(status);

      // step 7: solve U*y = z
      status = cusparseXcsrsv2_solve(handle,
                                     CUSPARSE_OPERATION_NON_TRANSPOSE,
                                     n_rows,
                                     n_nonzero_elements,
                                     &alpha,
                                     descr_U,
                                     A_val_dev,
                                     A_row_ptr_dev,
                                     A_column_index_dev,
                                     info_U,
                                     tmp_dev.get(),
                                     dst_dev,
                                     policy_U,
                                     buffer_dev);
      AssertCusparse(status);

      // step 8: free resources
      cudaFree(buffer_dev);
      status = cusparseDestroyMatDescr(descr_M);
      AssertCusparse(status);
      status = cusparseDestroyMatDescr(descr_L);
      AssertCusparse(status);
      status = cusparseDestroyMatDescr(descr_U);
      AssertCusparse(status);
      status = cusparseDestroyCsrilu02Info(info_M);
      AssertCusparse(status);
      status = cusparseDestroyCsrsv2Info(info_L);
      AssertCusparse(status);
      status = cusparseDestroyCsrsv2Info(info_U);
      AssertCusparse(status);
    }
  } // namespace CUDAWrappers
} // namespace dealii

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
  A.print(std::cout);

  // Solve on the host
  PreconditionIdentity prec_no;
  SolverControl        control(100, 1.e-10);
  SolverCG<>           cg_host(control);
  Vector<double>       sol_host(size);
  Vector<double>       rhs_host(size);
  for (unsigned int i = 0; i < size; ++i)
    rhs_host[i] = static_cast<double>(i);
  cg_host.solve(A, sol_host, rhs_host, prec_no);

  // Solve on the device
  CUDAWrappers::SparseMatrix<double>          A_dev(cuda_handle, A);
  LinearAlgebra::CUDAWrappers::Vector<double> sol_dev(size);
  LinearAlgebra::CUDAWrappers::Vector<double> rhs_dev(size);
  LinearAlgebra::ReadWriteVector<double>      rw_vector(size);
  for (unsigned int i = 0; i < size; ++i)
    rw_vector[i] = static_cast<double>(i);
  rhs_dev.import(rw_vector, VectorOperation::insert);
  SolverCG<LinearAlgebra::CUDAWrappers::Vector<double>> cg_dev(control);

  A_dev.print(std::cout);
  A_dev.print_formatted(std::cout);
  CUDAWrappers::PreconditionILU<double>    prec_double(cuda_handle);
  CUDAWrappers::PreconditionILU<float>     prec_float(cuda_handle);
  CUDAWrappers::PreconditionILU<cuComplex> prec_complex_float(cuda_handle);
  CUDAWrappers::PreconditionILU<cuDoubleComplex> prec_complex_double(
    cuda_handle);

  // apply_preconditioner(A_dev, cuda_handle.cusparse_handle, sol_dev, rhs_dev);
  // A_dev.print_formatted(std::cout);
  prec_double.initialize(A_dev);
  // A_dev.print_formatted(std::cout);
  // prec_double.vmult(sol_dev, rhs_dev);
  // A_dev.print_formatted(std::cout);
  cg_dev.solve(A_dev, sol_dev, rhs_dev, prec_double);

  // Check the result
  rw_vector.import(sol_dev, VectorOperation::insert);
  for (unsigned int i = 0; i < size; ++i)
    deallog << rw_vector[i] << " " << sol_host[i] << std::endl;
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
