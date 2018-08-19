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

#include <deal.II/lac/cuda/precondition.h>

DEAL_II_NAMESPACE_OPEN

namespace
{
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

  namespace CUDAWrappers
  {
    template <typename Number>
    PreconditionIC<Number>::AdditionalData::AdditionalData(
      bool use_level_analysis_)
      : use_level_analysis(use_level_analysis_)
    {}



    template <typename Number>
    PreconditionIC<Number>::PreconditionIC(
      const Utilities::CUDA::Handle &handle)
      : cusparse_handle(handle.cusparse_handle)
      , P_val_dev(nullptr, delete_device_vector<Number>)
      , P_row_ptr_dev(nullptr)
      , P_column_index_dev(nullptr)
      , tmp_dev(nullptr, delete_device_vector<Number>)
      , buffer_dev(nullptr, delete_device_vector<void>)
      , policy_L(CUSPARSE_SOLVE_POLICY_USE_LEVEL)
      , policy_Lt(CUSPARSE_SOLVE_POLICY_USE_LEVEL)
      , policy_M(CUSPARSE_SOLVE_POLICY_USE_LEVEL)
      , n_rows(0)
      , n_nonzero_elements(0)
    {
      cusparseStatus_t status;
      // step 1: create a descriptor which contains
      // - matrix M is base-0
      // - matrix L is base-0
      // - matrix L is lower triangular
      // - matrix L has non-unit diagonal
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
      status = cusparseSetMatDiagType(descr_L, CUSPARSE_DIAG_TYPE_NON_UNIT);
      AssertCusparse(status);

      // step 2: create a empty info structure
      // we need one info for csric02 and two info's for csrsv2
      status = cusparseCreateCsric02Info(&info_M);
      AssertCusparse(status);
      status = cusparseCreateCsrsv2Info(&info_L);
      AssertCusparse(status);
      status = cusparseCreateCsrsv2Info(&info_Lt);
      AssertCusparse(status);
    }



    template <typename Number>
    PreconditionIC<Number>::~PreconditionIC()
    {
      // step 8: free resources
      cusparseStatus_t status = cusparseDestroyMatDescr(descr_M);
      AssertNothrowCusparse(status);

      status = cusparseDestroyMatDescr(descr_L);
      AssertNothrowCusparse(status);

      status = cusparseDestroyCsric02Info(info_M);
      AssertNothrowCusparse(status);

      status = cusparseDestroyCsrsv2Info(info_L);
      AssertNothrowCusparse(status);

      status = cusparseDestroyCsrsv2Info(info_Lt);
      AssertNothrowCusparse(status);
    }



    template <typename Number>
    void
    PreconditionIC<Number>::initialize(const SparseMatrix<Number> &A,
                                       const AdditionalData &additional_data)
    {
      if (additional_data.use_level_analysis)
        {
          policy_L  = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
          policy_Lt = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
          policy_M  = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
        }
      else
        {
          policy_L  = CUSPARSE_SOLVE_POLICY_NO_LEVEL;
          policy_Lt = CUSPARSE_SOLVE_POLICY_NO_LEVEL;
          policy_M  = CUSPARSE_SOLVE_POLICY_NO_LEVEL;
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

      // step 3: query how much memory used in csric02 and csrsv2, and allocate
      // the buffer
      int              BufferSize_M;
      cusparseStatus_t status = cusparseXcsric02_bufferSize(cusparse_handle,
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

      int BufferSize_Lt;
      status = cusparseXcsrsv2_bufferSize(cusparse_handle,
                                          CUSPARSE_OPERATION_TRANSPOSE,
                                          n_rows,
                                          n_nonzero_elements,
                                          descr_L,
                                          P_val_dev.get(),
                                          P_row_ptr_dev,
                                          P_column_index_dev,
                                          info_Lt,
                                          &BufferSize_Lt);
      AssertCusparse(status);

      const int BufferSize =
        std::max(BufferSize_M, std::max(BufferSize_L, BufferSize_Lt));
      // workaround: since allocate_device_vector needs a type, we pass char
      // which is required to have size 1.
      buffer_dev.reset(static_cast<void *>(
        allocate_device_vector<char>(BufferSize / sizeof(char))));

      // step 4: perform analysis of incomplete Cholesky on M
      //         perform analysis of triangular solve on L
      //         perform analysis of triangular solve on L'
      // The lower triangular part of M has the same sparsity pattern as L, so
      // we can do analysis of csric02 and csrsv2 simultaneously.

      status = cusparseXcsric02_analysis(cusparse_handle,
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
        cusparseXcsric02_zeroPivot(cusparse_handle, info_M, &structural_zero);
      AssertCusparse(status);

      status = cusparseXcsrsv2_analysis(cusparse_handle,
                                        CUSPARSE_OPERATION_TRANSPOSE,
                                        n_rows,
                                        n_nonzero_elements,
                                        descr_L,
                                        P_val_dev.get(),
                                        P_row_ptr_dev,
                                        P_column_index_dev,
                                        info_Lt,
                                        policy_Lt,
                                        buffer_dev.get());
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

      // step 5: M = L * L'
      status = cusparseXcsric02(cusparse_handle,
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
        cusparseXcsric02_zeroPivot(cusparse_handle, info_M, &numerical_zero);
      AssertCusparse(status);
    }



    template <typename Number>
    void
    PreconditionIC<Number>::vmult(
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
      const double     alpha = 1.;
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

      // step 7: solve L'*y = alpha*z
      status = cusparseXcsrsv2_solve(cusparse_handle,
                                     CUSPARSE_OPERATION_TRANSPOSE,
                                     n_rows,
                                     n_nonzero_elements,
                                     &alpha,
                                     descr_L,
                                     P_val_dev.get(),
                                     P_row_ptr_dev,
                                     P_column_index_dev,
                                     info_Lt,
                                     tmp_dev.get(),
                                     dst_dev,
                                     policy_Lt,
                                     buffer_dev.get());
      AssertCusparse(status);
    }



    template <typename Number>
    void
    PreconditionIC<Number>::Tvmult(
      LinearAlgebra::CUDAWrappers::Vector<Number> &      dst,
      const LinearAlgebra::CUDAWrappers::Vector<Number> &src) const
    {
      // the constructed preconditioner is symmetric
      vmult(dst, src);
    }



    template <typename Number>
    PreconditionIC<Number>::size_type
    PreconditionIC<Number>::m() const
    {
      return n_rows;
    }



    template <typename Number>
    PreconditionIC<Number>::size_type
    PreconditionIC<Number>::n() const
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
      csric02Info_t               info_M  = 0;
      csrsv2Info_t                info_L  = 0;
      csrsv2Info_t                info_Lt = 0;
      int                         BufferSize_M;
      int                         BufferSize_L;
      int                         BufferSize_Lt;
      int                         BufferSize;
      void *                      buffer_dev = 0;
      int                         structural_zero;
      int                         numerical_zero;
      const double                alpha     = 1.;
      const cusparseSolvePolicy_t policy_M  = CUSPARSE_SOLVE_POLICY_NO_LEVEL;
      const cusparseSolvePolicy_t policy_L  = CUSPARSE_SOLVE_POLICY_NO_LEVEL;
      const cusparseSolvePolicy_t policy_Lt = CUSPARSE_SOLVE_POLICY_USE_LEVEL;

      cusparseStatus_t status;
      // step 1: create a descriptor which contains
      // - matrix M is base-0
      // - matrix L is base-0
      // - matrix L is lower triangular
      // - matrix L has non-unit diagonal
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
      status = cusparseSetMatDiagType(descr_L, CUSPARSE_DIAG_TYPE_NON_UNIT);
      AssertCusparse(status);

      // step 2: create a empty info structure
      // we need one info for csric02 and two info's for csrsv2
      status = cusparseCreateCsric02Info(&info_M);
      AssertCusparse(status);
      status = cusparseCreateCsrsv2Info(&info_L);
      AssertCusparse(status);
      status = cusparseCreateCsrsv2Info(&info_Lt);
      AssertCusparse(status);

      // step 3: query how much memory used in csric02 and csrsv2, and allocate
      // the buffer
      status = cusparseXcsric02_bufferSize(handle,
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
                                          CUSPARSE_OPERATION_TRANSPOSE,
                                          n_rows,
                                          n_nonzero_elements,
                                          descr_L,
                                          A_val_dev,
                                          A_row_ptr_dev,
                                          A_column_index_dev,
                                          info_Lt,
                                          &BufferSize_Lt);
      AssertCusparse(status);

      BufferSize = max(BufferSize_M, max(BufferSize_L, BufferSize_Lt));

      // buffer_dev returned by cudaMalloc is automatically aligned to 128
      // bytes.
      cudaError_t status_cuda = cudaMalloc((void **)&buffer_dev, BufferSize);
      Assert(cudaSuccess == status_cuda, ExcInternalError());

      // step 4: perform analysis of incomplete Cholesky on M
      //         perform analysis of triangular solve on L
      //         perform analysis of triangular solve on L'
      // The lower triangular part of M has the same sparsity pattern as L, so
      // we can do analysis of csric02 and csrsv2 simultaneously.

      status = cusparseXcsric02_analysis(handle,
                                         n_rows,
                                         n_nonzero_elements,
                                         descr_M,
                                         A_val_dev,
                                         A_row_ptr_dev,
                                         A_column_index_dev,
                                         info_M,
                                         policy_M,
                                         buffer_dev);
      AssertCusparse(status);
      status = cusparseXcsric02_zeroPivot(handle, info_M, &structural_zero);
      if (CUSPARSE_STATUS_ZERO_PIVOT == status)
        {
          printf("A(%d,%d) is missing\n", structural_zero, structural_zero);
        }

      status = cusparseXcsrsv2_analysis(handle,
                                        CUSPARSE_OPERATION_TRANSPOSE,
                                        n_rows,
                                        n_nonzero_elements,
                                        descr_L,
                                        A_val_dev,
                                        A_row_ptr_dev,
                                        A_column_index_dev,
                                        info_Lt,
                                        policy_Lt,
                                        buffer_dev);
      AssertCusparse(status);

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

      // step 5: M = L * L'
      status = cusparseXcsric02(handle,
                                n_rows,
                                n_nonzero_elements,
                                descr_M,
                                A_val_dev,
                                A_row_ptr_dev,
                                A_column_index_dev,
                                info_M,
                                policy_M,
                                buffer_dev);
      AssertCusparse(status);
      status = cusparseXcsric02_zeroPivot(handle, info_M, &numerical_zero);
      if (CUSPARSE_STATUS_ZERO_PIVOT == status)
        {
          printf("L(%d,%d) is zero\n", numerical_zero, numerical_zero);
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

      // step 7: solve L'*y = z
      status = cusparseXcsrsv2_solve(handle,
                                     CUSPARSE_OPERATION_TRANSPOSE,
                                     n_rows,
                                     n_nonzero_elements,
                                     &alpha,
                                     descr_L,
                                     A_val_dev,
                                     A_row_ptr_dev,
                                     A_column_index_dev,
                                     info_Lt,
                                     tmp_dev.get(),
                                     dst_dev,
                                     policy_Lt,
                                     buffer_dev);
      AssertCusparse(status);

      // step 8: free resources
      status_cuda = cudaFree(buffer_dev);
      AssertCuda(status_cuda);
      status = cusparseDestroyMatDescr(descr_M);
      AssertCusparse(status);
      status = cusparseDestroyMatDescr(descr_L);
      AssertCusparse(status);
      status = cusparseDestroyCsric02Info(info_M);
      AssertCusparse(status);
      status = cusparseDestroyCsrsv2Info(info_L);
      AssertCusparse(status);
      status = cusparseDestroyCsrsv2Info(info_Lt);
      AssertCusparse(status);
    }



    // explicit instantiations
    template class PreconditionIC<float>;
    template class PreconditionIC<double>;
  } // namespace CUDAWrappers

DEAL_II_NAMESPACE_CLOSE
