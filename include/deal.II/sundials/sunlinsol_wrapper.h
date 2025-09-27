// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_sundials_sunlinsol_wrapper_h
#define dealii_sundials_sunlinsol_wrapper_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SUNDIALS

#  include <sundials/sundials_linearsolver.h>

#  include <exception>
#  include <functional>
#  include <memory>

DEAL_II_NAMESPACE_OPEN

namespace SUNDIALS
{
#  ifndef DOXYGEN
  // forward declarations
  namespace internal
  {
    template <typename VectorType>
    struct LinearSolverContent;
  }
#  endif

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

#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    /**
     * Constructor.
     *
     * @param A_data Data required by @p a_times_fn
     * @param a_times_fn A function pointer to the function that computes A*v
     * @param linsol_ctx The context object used to set up the linear solver and all vectors
     */
    SundialsOperator(void       *A_data,
                     SUNATimesFn a_times_fn,
                     SUNContext  linsol_ctx);
#  else
    /**
     * Constructor.
     *
     * @param A_data Data required by @p a_times_fn
     * @param a_times_fn A function pointer to the function that computes A*v
     */
    SundialsOperator(void *A_data, ATimesFn a_times_fn);
#  endif

  private:
    /**
     * Data necessary to evaluate a_times_fn.
     */
    void *A_data;

#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    /**
     * %Function pointer declared by SUNDIALS to evaluate the matrix vector
     * product.
     */
    SUNATimesFn a_times_fn;

    /**
     * Context object used for SUNDIALS logging.
     */
    SUNContext linsol_ctx;
#  else
    /**
     * %Function pointer declared by SUNDIALS to evaluate the matrix vector
     * product.
     */
    ATimesFn a_times_fn;
#  endif
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

#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    /**
     * Constructor.
     *
     * @param P_data Data required by @p p_solve_fn
     * @param p_solve_fn A function pointer to the function that computes A*v
     * @param linsol_ctx The context object used to set up the linear solver and all vectors
     * @param tol Tolerance, that an iterative solver should use to judge
     *   convergence
     *
     * @note This function is only available with SUNDIALS 6.0.0 and later.
     * 6.0.0. If you are using an earlier version of SUNDIALS then you need to
     * use the other constructor.
     */
    SundialsPreconditioner(void       *P_data,
                           SUNPSolveFn p_solve_fn,
                           SUNContext  linsol_ctx,
                           double      tol);
#  else
    /**
     * Constructor.
     *
     * @param P_data Data required by @p p_solve_fn
     * @param p_solve_fn A function pointer to the function that computes A*v
     * @param tol Tolerance, that an iterative solver should use to judge
     *   convergence
     *
     * @note This function is only available with versions of SUNDIALS prior to
     * 6.0.0. If you are using SUNDIALS 6 or later then you need to use the
     * other constructor.
     */
    SundialsPreconditioner(void *P_data, PSolveFn p_solve_fn, double tol);
#  endif

  private:
    /**
     * Data necessary to calls p_solve_fn
     */
    void *P_data;

#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    /**
     * %Function pointer to a function that computes the preconditioner
     * application.
     */
    SUNPSolveFn p_solve_fn;

    /**
     * Context object used for SUNDIALS logging.
     */
    SUNContext linsol_ctx;
#  else
    /**
     * %Function pointer to a function that computes the preconditioner
     * application.
     */
    PSolveFn p_solve_fn;
#  endif

    /**
     * Potential tolerance to use in the internal solve of the preconditioner
     * equation.
     */
    double tol;
  };

  /**
   * Type of function objects to interface with SUNDIALS' linear solvers
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
   *
   * @note This variable represents a
   * @ref GlossUserProvidedCallBack "user provided callback".
   * See there for a description of how to deal with errors and other
   * requirements and conventions. In particular, ARKode can deal
   * with "recoverable" errors in some circumstances, so callbacks
   * can throw exceptions of type RecoverableUserCallbackError.
   */
  template <typename VectorType>
  using LinearSolveFunction =
    std::function<void(SundialsOperator<VectorType>       &op,
                       SundialsPreconditioner<VectorType> &prec,
                       VectorType                         &x,
                       const VectorType                   &b,
                       double                              tol)>;

  namespace internal
  {
    /**
     * Attach wrapper functions to SUNDIALS' linear solver interface. We pretend
     * that the user-supplied linear solver is matrix-free, even though it can
     * be matrix-based. This way SUNDIALS does not need to understand our matrix
     * types.
     */
    template <typename VectorType>
    class LinearSolverWrapper
    {
    public:
      explicit LinearSolverWrapper(
        const LinearSolveFunction<VectorType> &lsolve,
        std::exception_ptr                    &pending_exception
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
        ,
        SUNContext linsol_ctx
#  endif
      );

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

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif
#endif
