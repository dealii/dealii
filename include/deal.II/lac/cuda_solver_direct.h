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
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

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
   * @author Bruno Turcksin
   * @date 2018
   */
  template <typename Number>
  class SolverDirect
  {
  public:
    struct AdditionalData
    {
      /**
       * Set the additional data field to the desired solver.
       */
      explicit AdditionalData(const std::string& solver_type = "LU_dense");

      /**
       * Set the solver type. Possibilities are:
       * <ul>
       * <li> "Cholesky" which performs a Cholesky decomposition on the device </li>
       * <li> "LU_dense" which converts the sparse matrix to a dense matrix and
       * uses LU factorization </li>
       * <li> "LU_host" which uses LU factorization on the host </li>
       * </ul>
       */
      std::string solver_type;
    };

    /**
     * Constructor. Takes the solver control object and creates the solver.
     */
    SolverDirect(const Utilities::CUDA::Handle& handle,
                 SolverControl&                 cn,
                 const AdditionalData&          data = AdditionalData());

    /**
     * Destructor.
     */
    virtual ~SolverDirect() = default;

    /**
     * Solve the linear system <tt>Ax=b</tt>.
     */
    void
    solve(const SparseMatrix<Number>&                        A,
          LinearAlgebra::CUDAWrappers::Vector<Number>&       x,
          const LinearAlgebra::CUDAWrappers::Vector<Number>& b);

    /**
     * Access to object that controls convergence.
     */
    SolverControl&
    control() const;

  private:
    /**
     * Handle
     */
    const Utilities::CUDA::Handle& cuda_handle;

    /**
     * Reference to the object that controls convergence of the iterative
     * solver. In fact, for these Trilinos wrappers, Trilinos does so itself,
     * but we copy the data from this object before starting the solution
     * process, and copy the data back into it afterwards.
     */
    SolverControl& solver_control;

    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;
  };
} // namespace CUDAWrappers

DEAL_II_NAMESPACE_CLOSE

#endif

#endif
