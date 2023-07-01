// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2023 by the deal.II authors
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

#include <deal.II/lac/cuda_precondition.h>
#include <deal.II/lac/cuda_sparse_matrix.h>

DEAL_II_NAMESPACE_OPEN

namespace
{
  /**
   * Template wrapper for cusparse<t>csrilu02.
   * (https://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-csrilu02).
   * function performs the solve phase of the incomplete-LU factorization with
   * 0 fill-in and no pivoting.
   */
  template <typename Number>
  cusparseStatus_t
  cusparseXcsrilu02(cusparseHandle_t /*handle*/,
                    int /*m*/,
                    int /*nnz*/,
                    const cusparseMatDescr_t /*descrA*/,
                    Number * /*csrValA_valM*/,
                    const int * /*csrRowPtrA*/,
                    const int * /*csrColIndA*/,
                    csrilu02Info_t /*info*/,
                    cusparseSolvePolicy_t /*policy*/,
                    void * /*pBuffer*/)
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

  /*
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
  */



  /**
   * Template wrapper for cusparse<t>csrilu02_analysis.
   * (https://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-csrilu02_analysis).
   * This function performs the analysis phase of the incomplete-LU
   * factorization with 0 fill-in and no pivoting.
   */
  template <typename Number>
  cusparseStatus_t
  cusparseXcsrilu02_analysis(cusparseHandle_t /*handle*/,
                             int /*m*/,
                             int /*nnz*/,
                             const cusparseMatDescr_t /*descrA*/,
                             const Number * /*csrValA*/,
                             const int * /*csrRowPtrA*/,
                             const int * /*csrColIndA*/,
                             csrilu02Info_t /*info*/,
                             cusparseSolvePolicy_t /*policy*/,
                             void * /*pBuffer*/)
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

  /*
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
*/



  /**
   * Template wrapper for cusparse<t>csrilu02_bufferSize.
   * (https://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-csrilu02_bufferSize).
   * This function returns size of the buffer used in computing the
   * incomplete-LU factorization with 0 fill-in and no pivoting.
   */
  template <typename Number>
  cusparseStatus_t
  cusparseXcsrilu02_bufferSize(cusparseHandle_t /*handle*/,
                               int /*m*/,
                               int /*nnz*/,
                               const cusparseMatDescr_t /*descrA*/,
                               Number * /*csrValA*/,
                               const int * /*csrRowPtrA*/,
                               const int * /*csrColIndA*/,
                               csrilu02Info_t /*info*/,
                               int * /*pBufferSizeInBytes*/)
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

  /*
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
    cusparseXcsrilu02_bufferSize<cuDoubleComplex>(cusparseHandle_t handle, int
    m, int                      nnz, const cusparseMatDescr_t descrA,
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
*/



  /**
   * Template wrapper for cusparse<t>csric02
   * (https://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-csric02).
   * This function performs the solve phase of the computing the
   * incomplete-Cholesky factorization with 0 fill-in and no pivoting.
   */
  template <typename Number>
  cusparseStatus_t
  cusparseXcsric02(cusparseHandle_t /*handle*/,
                   int /*m*/,
                   int /*nnz*/,
                   const cusparseMatDescr_t /*descrA*/,
                   Number * /*csrValA_valM*/,
                   const int * /*csrRowPtrA*/,
                   const int * /*csrColIndA*/,
                   csric02Info_t /*info*/,
                   cusparseSolvePolicy_t /*policy*/,
                   void * /*pBuffer*/)
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

  /*
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
    */



  /**
   * Template wrapper for cusparse<t>csrsv2_solve
   * (https://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-csrsv2_solve).
   * This function performs the solve phase of csrsv2, a new sparse triangular
   * linear system op(A)*y = alpha*x.
   */
  template <typename Number>
  cusparseStatus_t
  cusparseXcsrsv2_solve(cusparseHandle_t /*handle*/,
                        cusparseOperation_t /*transA*/,
                        int /*m*/,
                        int /*nnz*/,
                        const Number * /*alpha*/,
                        const cusparseMatDescr_t /*descra*/,
                        const Number * /*csrValA*/,
                        const int * /*csrRowPtrA*/,
                        const int * /*csrColIndA*/,
                        csrsv2Info_t /*info*/,
                        const Number * /*x*/,
                        Number * /*y*/,
                        cusparseSolvePolicy_t /*policy*/,
                        void * /*pBuffer*/)
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

  /*
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
*/



  /**
   * Template wrapper for cusparse<t>csrsv2_analysis
   * (https://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-csrsv2_analysis).
   * This function performs the analysis phase of csrsv2, a new sparse
   * triangular linear system op(A)*y = alpha*x.
   */
  template <typename Number>
  cusparseStatus_t
  cusparseXcsrsv2_analysis(cusparseHandle_t /*handle*/,
                           cusparseOperation_t /*transA*/,
                           int /*m*/,
                           int /*nnz*/,
                           const cusparseMatDescr_t /*descrA*/,
                           const Number * /*csrValA*/,
                           const int * /*csrRowPtrA*/,
                           const int * /*csrColIndA*/,
                           csrsv2Info_t /*info*/,
                           cusparseSolvePolicy_t /*policy*/,
                           void * /*pBuffer*/)
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

  /*
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
*/


  /**
   * Template wrapper for cusparse<t>csric02_analysis
   * (https://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-csric02_analysis).
   * This function performs the analysis phase of the incomplete-Cholesky
   * factorization with 0 fill-in and no pivoting.
   */
  template <typename Number>
  cusparseStatus_t
  cusparseXcsric02_analysis(cusparseHandle_t /*handle*/,
                            int /*m*/,
                            int /*nnz*/,
                            const cusparseMatDescr_t /*descrA*/,
                            const Number * /*csrValA*/,
                            const int * /*csrRowPtrA*/,
                            const int * /*csrColIndA*/,
                            csric02Info_t /*info*/,
                            cusparseSolvePolicy_t /*policy*/,
                            void * /*pBuffer*/)
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

  /*
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
*/



  /**
   * Template wrapper for cusparse<t>csrsv2_bufferSize
   * (https://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-csrsv2_bufferSize).
   * This function returns the size of the buffer used in csrsv2, a new sparse
   * triangular linear system op(A)*y = alpha*x.
   */
  template <typename Number>
  cusparseStatus_t
  cusparseXcsrsv2_bufferSize(cusparseHandle_t /*handle*/,
                             cusparseOperation_t /*transA*/,
                             int /*m*/,
                             int /*nnz*/,
                             const cusparseMatDescr_t /*descrA*/,
                             Number * /*csrValA*/,
                             const int * /*csrRowPtrA*/,
                             const int * /*csrColIndA*/,
                             csrsv2Info_t /*info*/,
                             int * /*pBufferSizeInBytes*/)
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

  /*
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
*/



  /**
   * Template wrapper for cusparse<t>csric02_bufferSize
   * (https://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-csric02_bufferSize).
   * This function returns size of buffer used in computing the
   * incomplete-Cholesky factorization with 0 fill-in and no pivoting.
   */
  template <typename Number>
  cusparseStatus_t
  cusparseXcsric02_bufferSize(cusparseHandle_t /*handle*/,
                              int /*m*/,
                              int /*nnz*/,
                              const cusparseMatDescr_t /*descrA*/,
                              Number * /*csrValA*/,
                              const int * /*csrRowPtrA*/,
                              const int * /*csrColIndA*/,
                              csric02Info_t /*info*/,
                              int * /*pBufferSizeInBytes*/)
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

  /*
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
  */
} // namespace

namespace CUDAWrappers
{
  template <typename Number>
  PreconditionIC<Number>::AdditionalData::AdditionalData(
    bool use_level_analysis_)
    : use_level_analysis(use_level_analysis_)
  {}



  template <typename Number>
  PreconditionIC<Number>::PreconditionIC(const Utilities::CUDA::Handle &handle)
    : cusparse_handle(handle.cusparse_handle)
    , P_val_dev(nullptr, Utilities::CUDA::delete_device_data<Number>)
    , P_row_ptr_dev(nullptr)
    , P_column_index_dev(nullptr)
    , tmp_dev(nullptr, Utilities::CUDA::delete_device_data<Number>)
    , buffer_dev(nullptr, Utilities::CUDA::delete_device_data<void>)
    , policy_L(CUSPARSE_SOLVE_POLICY_USE_LEVEL)
    , policy_Lt(CUSPARSE_SOLVE_POLICY_USE_LEVEL)
    , policy_M(CUSPARSE_SOLVE_POLICY_USE_LEVEL)
    , n_rows(0)
    , n_nonzero_elements(0)
  {
    // step 1: create a descriptor which contains
    // - matrix M is base-0
    // - matrix L is base-0
    // - matrix L is lower triangular
    // - matrix L has non-unit diagonal
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

    matrix_pointer = &A;
    const Number *A_val_dev;
    std::tie(A_val_dev,
             P_column_index_dev,
             P_row_ptr_dev,
             std::ignore,
             std::ignore) = A.get_cusparse_matrix();

    // create a copy of the matrix entries since the algorithm works in-place.
    P_val_dev.reset(
      Utilities::CUDA::allocate_device_data<Number>(n_nonzero_elements));
    cudaError_t cuda_status = cudaMemcpy(P_val_dev.get(),
                                         A_val_dev,
                                         n_nonzero_elements * sizeof(Number),
                                         cudaMemcpyDeviceToDevice);
    AssertCuda(cuda_status);

    // initialize an internal buffer we need later on
    tmp_dev.reset(Utilities::CUDA::allocate_device_data<Number>(n_rows));

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
    // workaround: since allocate_device_data needs a type, we pass char
    // which is required to have size 1.
    buffer_dev.reset(static_cast<void *>(
      Utilities::CUDA::allocate_device_data<char>(BufferSize / sizeof(char))));

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
    const Number     alpha = internal::NumberType<Number>::value(1.);
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
  PreconditionILU<Number>::AdditionalData::AdditionalData(
    bool use_level_analysis_)
    : use_level_analysis(use_level_analysis_)
  {}



  template <typename Number>
  PreconditionILU<Number>::PreconditionILU(
    const Utilities::CUDA::Handle &handle)
    : cusparse_handle(handle.cusparse_handle)
    , P_val_dev(nullptr, Utilities::CUDA::delete_device_data<Number>)
    , P_row_ptr_dev(nullptr)
    , P_column_index_dev(nullptr)
    , tmp_dev(nullptr, Utilities::CUDA::delete_device_data<Number>)
    , buffer_dev(nullptr, Utilities::CUDA::delete_device_data<void>)
    , policy_L(CUSPARSE_SOLVE_POLICY_USE_LEVEL)
    , policy_U(CUSPARSE_SOLVE_POLICY_USE_LEVEL)
    , policy_M(CUSPARSE_SOLVE_POLICY_USE_LEVEL)
    , n_rows(0)
    , n_nonzero_elements(0)
  {
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

    matrix_pointer     = &A;
    n_rows             = A.m();
    n_nonzero_elements = A.n_nonzero_elements();
    AssertDimension(A.m(), A.n());

    const Number *A_val_dev;
    std::tie(A_val_dev,
             P_column_index_dev,
             P_row_ptr_dev,
             std::ignore,
             std::ignore) = A.get_cusparse_matrix();

    // create a copy of the matrix entries since the algorithm works in-place.
    P_val_dev.reset(
      Utilities::CUDA::allocate_device_data<Number>(n_nonzero_elements));
    cudaError_t cuda_status = cudaMemcpy(P_val_dev.get(),
                                         A_val_dev,
                                         n_nonzero_elements * sizeof(Number),
                                         cudaMemcpyDeviceToDevice);
    AssertCuda(cuda_status);

    // initialize an internal buffer we need later on
    tmp_dev.reset(Utilities::CUDA::allocate_device_data<Number>(n_rows));

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
    // workaround: since allocate_device_data needs a type, we pass char
    // which is required to have size 1.
    buffer_dev.reset(static_cast<void *>(
      Utilities::CUDA::allocate_device_data<char>(BufferSize / sizeof(char))));

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
    const Number     alpha = internal::NumberType<Number>::value(1.);
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
    LinearAlgebra::CUDAWrappers::Vector<Number> & /*dst*/,
    const LinearAlgebra::CUDAWrappers::Vector<Number> & /*src*/) const
  {
    Assert(false, ExcNotImplemented());
  }



  // explicit instantiations
  template class PreconditionIC<float>;
  template class PreconditionIC<double>;
  // template class PreconditionIC<cuComplex>;
  // template class PreconditionIC<cuDoubleComplex>;
  template class PreconditionILU<float>;
  template class PreconditionILU<double>;
  // template class PreconditionILU<cuComplex>;
  // template class PreconditionILU<cuDoubleComplex>;
} // namespace CUDAWrappers

DEAL_II_NAMESPACE_CLOSE
