// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
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

#ifndef dealii_cuda_precondition_h
#define dealii_cuda_precondition_h

#include <deal.II/base/config.h>

#include <deal.II/base/cuda.h>
#include <deal.II/base/smartpointer.h>

#include <memory>

#ifdef DEAL_II_COMPILER_CUDA_AWARE

DEAL_II_NAMESPACE_OPEN

// forward-definition
#  ifndef DOXYGEN
namespace LinearAlgebra
{
  namespace CUDAWrappers
  {
    template <typename Number>
    class Vector;
  }
} // namespace LinearAlgebra
#  endif

namespace CUDAWrappers
{
  // forward definition
  template <typename Number>
  class SparseMatrix;

  /**
   * This class implements an incomplete Cholesky factorization (IC)
   * preconditioner for @em symmetric CUDAWrappers::SparseMatrix matrices.
   *
   * The implementation closely follows the one documented in the cuSPARSE
   * documentation
   * (https://docs.nvidia.com/cuda/cusparse/index.html#cusparse-lt-t-gt-csric02).
   *
   * @note Instantiations for this template are provided for <tt>@<float@> and
   * @<double@></tt>.
   *
   * @ingroup Preconditioners CUDAWrappers
   * @author Daniel Arndt
   * @date 2018
   */
  template <typename Number>
  class PreconditionIC
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
       * According to the documentation this might improve performance.
       * It is suggested to try both options.
       */
      AdditionalData(bool use_level_analysis = true);

      /**
       * Flag that determines if level information is used when creating and
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
    PreconditionIC(const Utilities::CUDA::Handle &handle);

    /**
     * The copy constructor is deleted.
     */
    PreconditionIC(const PreconditionIC<Number> &) = delete;

    /**
     * The copy assignment operator is deleted.
     */
    PreconditionIC &
    operator=(const PreconditionIC<Number> &) = delete;

    /**
     * Destructor. Free all resources that were initialized in this class.
     */
    ~PreconditionIC();

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
               const AdditionalData &      additional_data = AdditionalData());

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
     * Return the dimension of the codomain (or range) space. Note that the
     * matrix is square and has dimension $m \times m$.
     *
     * @note This function should only be called if the preconditioner has been
     * initialized.
     */
    size_type
    m() const;

    /**
     * Return the dimension of the codomain (or range) space. Note that the
     * matrix is square and has dimension $n \times n$.
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
     * cuSPARSE description of the sparse matrix $M=LL^T$.
     */
    cusparseMatDescr_t descr_M;

    /**
     * cuSPARSE description of the lower triangular matrix $L$.
     */
    cusparseMatDescr_t descr_L;

    /**
     * Solve and analysis structure for $M=LL^T$.
     */
    csric02Info_t info_M;

    /**
     * Solve and analysis structure for the lower triangular matrix $L$.
     */
    csrsv2Info_t info_L;

    /**
     * Solve and analysis structure for the upper triangular matrix $L^T$.
     */
    csrsv2Info_t info_Lt;

    /**
     * Pointer to the matrix this object was initialized with.
     */
    SmartPointer<const SparseMatrix<Number>> matrix_pointer;

    /**
     * Pointer to the values (on the device) of the computed preconditioning
     * matrix.
     */
    std::unique_ptr<Number[], void (*)(Number *)> P_val_dev;

    /**
     * Pointer to the row pointer (on the device) of the sparse matrix this
     * object was initialized with. Guarded by matrix_pointer.
     */
    const int *P_row_ptr_dev;

    /**
     * Pointer to the column indices (on the device) of the sparse matrix this
     * object was initialized with. Guarded by matrix_pointer.
     */
    const int *P_column_index_dev;

    /**
     * Pointer to the value (on the device) for a temporary (helper) vector
     * used in vmult().
     */
    std::unique_ptr<Number[], void (*)(Number *)> tmp_dev;

    /**
     * Pointer to an internal buffer (on the device) that is used for
     * computing the decomposition.
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
    cusparseSolvePolicy_t policy_Lt;

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

  /**
   * This class implements an incomplete LU factorization preconditioner for
   * CUDAWrappers::SparseMatrix matrices.
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
       *  to the documentation this might improve performance.
       * It is suggested to try both options.
       */
      AdditionalData(bool use_level_analysis = true);

      /**
       * Flag that determines if level information is used when creating and
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
               const AdditionalData &      additional_data = AdditionalData());

    /**
     * Apply the preconditioner.
     */
    void
    vmult(LinearAlgebra::CUDAWrappers::Vector<Number> &      dst,
          const LinearAlgebra::CUDAWrappers::Vector<Number> &src) const;

    /**
     * Apply the transposed preconditioner. Not yet implemented.
     */
    void
    Tvmult(LinearAlgebra::CUDAWrappers::Vector<Number> &      dst,
           const LinearAlgebra::CUDAWrappers::Vector<Number> &src) const;

    /**
     * Return the dimension of the codomain (or range) space. Note that the
     * matrix is square and has dimension $m \times m$.
     *
     * @note This function should only be called if the preconditioner has been
     * initialized.
     */
    size_type
    m() const;

    /**
     * Return the dimension of the codomain (or range) space. Note that the
     * matrix is square and has dimension $n \times n$.
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
     * Solve and analysis structure for $M=LU$.
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
     * Pointer to the matrix this object was initialized with.
     */
    SmartPointer<const SparseMatrix<Number>> matrix_pointer;

    /**
     * Pointer to the values (on the device) of the computed preconditioning
     * matrix.
     */
    std::unique_ptr<Number[], void (*)(Number *)> P_val_dev;

    /**
     * Pointer to the row pointer (on the device) of the sparse matrix this
     * object was initialized with. Guarded by matrix_pointer.
     */
    const int *P_row_ptr_dev;

    /**
     * Pointer to the column indices (on the device) of the sparse matrix this
     * object was initialized with. Guarded by matrix_pointer.
     */
    const int *P_column_index_dev;

    /**
     * Pointer to the value (on the device) for a temporary (helper) vector
     * used in vmult().
     */
    std::unique_ptr<Number[], void (*)(Number *)> tmp_dev;

    /**
     * Pointer to an internal buffer (on the device) that is used for
     * computing the decomposition.
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
     * triangular matrix $U$. This value can be modified through an
     * AdditionalData object.
     */
    cusparseSolvePolicy_t policy_U;

    /**
     * Determine if level information should be generated for $M=LU$. This
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

  /*--------------------------- inline functions ----------------------------*/

#  ifndef DOXYGEN
  template <typename Number>
  inline typename PreconditionIC<Number>::size_type
  PreconditionIC<Number>::m() const
  {
    return n_rows;
  }



  template <typename Number>
  inline typename PreconditionIC<Number>::size_type
  PreconditionIC<Number>::n() const
  {
    return n_rows;
  }



  template <typename Number>
  inline typename PreconditionILU<Number>::size_type
  PreconditionILU<Number>::m() const
  {
    return n_rows;
  }



  template <typename Number>
  inline typename PreconditionILU<Number>::size_type
  PreconditionILU<Number>::n() const
  {
    return n_rows;
  }
#  endif // DOXYGEN

} // namespace CUDAWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_CUDA

#endif // dealii_cuda_precondition_h
