//-----------------------------------------------------------
//
//    Copyright (C) 2021 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//-----------------------------------------------------------

#ifndef dealii_sundials_sunlinsol_wrapper_h
#define dealii_sundials_sunlinsol_wrapper_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SUNDIALS
#  if DEAL_II_SUNDIALS_VERSION_GTE(4, 0, 0)

#    include <sundials/sundials_linearsolver.h>

#    include <functional>
#    include <memory>

DEAL_II_NAMESPACE_OPEN

namespace SUNDIALS
{
#    ifndef DOXYGEN
  // forward declarations
  namespace internal
  {
    template <typename VectorType>
    struct LinearSolverContent;
  }
#    endif

  /**
   * A linear operator that wraps SUNDIALS functionality.
   */
  template <typename VectorType>
  struct SundialsOperator
  {
    /**
     * Apply this LinearOperator to @p src and store the result in @p dst.
     */
    void
    vmult(VectorType &dst, const VectorType &src) const;

    /**
     * Constructor.
     *
     * @param A_data Data required by @p a_times_fn
     * @param a_times_fn A function pointer to the function that computes A*v
     */
    SundialsOperator(void *A_data, ATimesFn a_times_fn);

  private:
    /**
     * Data necessary to evaluate a_times_fn.
     */
    void *A_data;

    /**
     * %Function pointer declared by SUNDIALS to evaluate the matrix vector
     * product.
     */
    ATimesFn a_times_fn;
  };



  /**
   * A linear operator that wraps preconditioner functionality as specified by
   * SUNDIALS. The vmult() function solves the preconditioner equation $Px=b$,
   * i.e., it computes $x=P^{-1}b$.
   */
  template <typename VectorType>
  struct SundialsPreconditioner
  {
    /**
     * Apply the wrapped preconditioner, i.e., solve $Px=b$ where $x$ is the
     * @p dst vector and $b$ the @p src vector.
     *
     * @param dst Result vector of the preconditioner application
     * @param src Target vector of the preconditioner application
     */
    void
    vmult(VectorType &dst, const VectorType &src) const;

    /**
     * Constructor.
     *
     * @param P_data Data required by @p p_solve_fn
     * @param p_solve_fn A function pointer to the function that computes A*v
     * @param tol Tolerance, that an iterative solver should use to judge
     *   convergence
     */
    SundialsPreconditioner(void *P_data, PSolveFn p_solve_fn, double tol);

  private:
    /**
     * Data necessary to calls p_solve_fn
     */
    void *P_data;

    /**
     * %Function pointer to a function that computes the preconditioner
     * application.
     */
    PSolveFn p_solve_fn;

    /**
     * Potential tolerance to use in the internal solve of the preconditioner
     * equation.
     */
    double tol;
  };

  /**
   * Type of function objects to interface with SUNDIALS linear solvers
   *
   * This function type encapsulates the action of solving $P^{-1}Ax=P^{-1}b$.
   * The LinearOperator @p op encapsulates the matrix vector product $Ax$ and
   * the LinearOperator @p prec encapsulates the application of the
   * preconditioner $P^{-1}z$.
   * The user can specify function objects of this type to attach custom linear
   * solver routines to SUNDIALS. The two LinearOperators @p op and @p prec are
   * built internally by SUNDIALS based on user settings. The parameters are
   * interpreted as follows:
   *
   * @param[in] op A LinearOperator that applies the matrix vector product
   * @param[in] prec A LinearOperator that applies the preconditioner
   * @param[out] x The output solution vector
   * @param[in] b The right-hand side
   * @param[in] tol Tolerance for the iterative solver
   *
   * This function should return:
   * - 0: Success
   * - >0: Recoverable error, ARKode will reattempt the solution and call this
   *       function again.
   * - <0: Unrecoverable error, the computation will be aborted and an
   *       assertion will be thrown.
   */
  template <typename VectorType>
  using LinearSolveFunction =
    std::function<int(SundialsOperator<VectorType> &      op,
                      SundialsPreconditioner<VectorType> &prec,
                      VectorType &                        x,
                      const VectorType &                  b,
                      double                              tol)>;

  namespace internal
  {
    /*!
     * Attach wrapper functions to SUNDIALS' linear solver interface. We pretend
     * that the user-supplied linear solver is matrix-free, even though it can
     * be matrix-based. This way SUNDIALS does not need to understand our matrix
     * types.
     */
    template <typename VectorType>
    class LinearSolverWrapper
    {
    public:
      explicit LinearSolverWrapper(LinearSolveFunction<VectorType> lsolve);

      ~LinearSolverWrapper();

      /**
       * Implicit conversion to SUNLinearSolver.
       */
      operator SUNLinearSolver();

    private:
      SUNLinearSolver                                  sun_linear_solver;
      std::unique_ptr<LinearSolverContent<VectorType>> content;
    };
  } // namespace internal
} // namespace SUNDIALS

DEAL_II_NAMESPACE_CLOSE

#  endif
#endif
#endif
