// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_cuda_solver_direct_h
#define dealii_cuda_solver_direct_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_CUDA
#  include <deal.II/base/cuda.h>

#  include <deal.II/lac/cuda_sparse_matrix.h>
#  include <deal.II/lac/cuda_vector.h>
#  include <deal.II/lac/solver_control.h>

DEAL_II_NAMESPACE_OPEN

namespace CUDAWrappers
{
  /**
   * Direct solvers. These solvers call cuSOLVER underneath.
   *
   * @note Instantiations for this template are provided for <tt>@<float@></tt>
   * and <tt>@<double@></tt>.
   *
   * @ingroup CUDAWrappers
   */
  template <typename Number>
  class SolverDirect
  {
  public:
    /**
     * Struct for additional settings for SolverDirect.
     */
    struct AdditionalData
    {
      /**
       * Set the additional data field to the desired solver.
       */
      explicit AdditionalData(const std::string &solver_type = "LU_dense");

      /**
       * Set the solver type. Possibilities are:
       * <ul>
       * <li> "Cholesky" which performs a Cholesky decomposition on the
       * @ref GlossDevice "device"
       * </li>
       * <li> "LU_dense" which converts the sparse matrix to a dense
       * matrix and uses LU factorization </li>
       * <li> "LU_host" which uses LU factorization on the host </li>
       * </ul>
       */
      std::string solver_type;
    };

    /**
     * Constructor. Takes the solver control object and creates the solver.
     */
    SolverDirect(const Utilities::CUDA::Handle &handle,
                 SolverControl                 &cn,
                 const AdditionalData          &data = AdditionalData());

    /**
     * Destructor.
     */
    virtual ~SolverDirect() = default;

    /**
     * Solve the linear system <tt>Ax=b</tt>.
     */
    void
    solve(const SparseMatrix<Number>                        &A,
          LinearAlgebra::CUDAWrappers::Vector<Number>       &x,
          const LinearAlgebra::CUDAWrappers::Vector<Number> &b);

    /**
     * Access to object that controls convergence.
     */
    SolverControl &
    control() const;

  private:
    /**
     * Handle
     */
    const Utilities::CUDA::Handle &cuda_handle;

    /**
     * Reference to the object that controls convergence of the iterative
     * solver. In fact, for these CUDA wrappers, cuSOLVER and cuSPARSE do so
     * themselves, but we copy the data from this object before starting the
     * solution process, and copy the data back into it afterwards.
     */
    SolverControl &solver_control;

    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;
  };
} // namespace CUDAWrappers

DEAL_II_NAMESPACE_CLOSE

#endif

#endif
