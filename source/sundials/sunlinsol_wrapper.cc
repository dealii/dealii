// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/config.h>

#include <deal.II/sundials/sunlinsol_wrapper.h>

#ifdef DEAL_II_WITH_SUNDIALS

#  include <deal.II/base/exceptions.h>

#  include <deal.II/lac/block_vector.h>
#  include <deal.II/lac/la_parallel_block_vector.h>
#  include <deal.II/lac/la_parallel_vector.h>
#  include <deal.II/lac/vector.h>
#  ifdef DEAL_II_WITH_TRILINOS
#    include <deal.II/lac/trilinos_parallel_block_vector.h>
#    include <deal.II/lac/trilinos_vector.h>
#  endif
#  ifdef DEAL_II_WITH_PETSC
#    include <deal.II/lac/petsc_block_vector.h>
#    include <deal.II/lac/petsc_vector.h>
#  endif

#  include <deal.II/sundials/n_vector.h>
#  include <deal.II/sundials/sundials_types.h>

#  include <sundials/sundials_iterative.h>
#  include <sundials/sundials_linearsolver.h>

DEAL_II_NAMESPACE_OPEN



namespace SUNDIALS
{
  DeclException1(ExcSundialsSolverError,
                 int,
                 << "One of the SUNDIALS linear solver internal"
                 << " functions returned a negative error code: " << arg1
                 << ". Please consult SUNDIALS manual.");

#  define AssertSundialsSolver(code) \
    Assert(code >= 0, ExcSundialsSolverError(code))

  namespace
  {
    /**
     * A function that calls the function object given by its first argument
     * with the set of arguments following at the end. If the call returns
     * regularly, the current function returns zero to indicate success. If the
     * call fails with an exception of type RecoverableUserCallbackError, then
     * the current function returns 1 to indicate that the called function
     * object thought the error it encountered is recoverable. If the call fails
     * with any other exception, then the current function returns with an error
     * code of -1. In each of the last two cases, the exception thrown by `f`
     * is captured and `eptr` is set to the exception. In case of success,
     * `eptr` is set to `nullptr`.
     */
    template <typename F, typename... Args>
    int
    call_and_possibly_capture_exception(const F            &f,
                                        std::exception_ptr &eptr,
                                        Args &&...args)
    {
      // See whether there is already something in the exception pointer
      // variable. This can only happen if we had previously had
      // a recoverable exception, and the underlying library actually
      // did recover successfully. In that case, we can abandon the
      // exception previously thrown. If eptr contains anything other,
      // then we really don't know how that could have happened, and
      // should probably bail out:
      if (eptr)
        {
          try
            {
              std::rethrow_exception(eptr);
            }
          catch (const RecoverableUserCallbackError &)
            {
              // ok, ignore, but reset the pointer
              eptr = nullptr;
            }
          catch (...)
            {
              // uh oh:
              AssertThrow(false, ExcInternalError());
            }
        }

      // Call the function and if that succeeds, return zero:
      try
        {
          f(std::forward<Args>(args)...);
          eptr = nullptr;
          return 0;
        }
      // If the call failed with a recoverable error, then
      // ignore the exception for now (but store a pointer to it)
      // and return a positive return value (+1). If the underlying
      // implementation manages to recover
      catch (const RecoverableUserCallbackError &)
        {
          eptr = std::current_exception();
          return 1;
        }
      // For any other exception, capture the exception and
      // return -1:
      catch (const std::exception &)
        {
          eptr = std::current_exception();
          return -1;
        }
    }
  } // namespace


  namespace internal
  {
    /**
     * storage for internal content of the linear solver wrapper
     */
    template <typename VectorType>
    struct LinearSolverContent
    {
      LinearSolverContent(std::exception_ptr &pending_exception)
        : a_times_fn(nullptr)
        , preconditioner_setup(nullptr)
        , preconditioner_solve(nullptr)
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
        , linsol_ctx(nullptr)
#  endif
        , P_data(nullptr)
        , A_data(nullptr)
        , pending_exception(pending_exception)
      {}

#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
      SUNATimesFn a_times_fn;
      SUNPSetupFn preconditioner_setup;
      SUNPSolveFn preconditioner_solve;

      SUNContext linsol_ctx;
#  else
      ATimesFn a_times_fn;
      PSetupFn preconditioner_setup;
      PSolveFn preconditioner_solve;
#  endif

      LinearSolveFunction<VectorType> lsolve;

      void *P_data;
      void *A_data;

      /**
       * A reference to a location where we can store exceptions, should they
       * be thrown by a linear solver object.
       */
      std::exception_ptr &pending_exception;
    };
  } // namespace internal

  namespace
  {
    using namespace SUNDIALS::internal;

    /**
     * Access our LinearSolverContent from the generic content of the
     * SUNLinearSolver @p ls.
     */
    template <typename VectorType>
    LinearSolverContent<VectorType> *
    access_content(SUNLinearSolver ls)
    {
      Assert(ls->content != nullptr, ExcInternalError());
      return static_cast<LinearSolverContent<VectorType> *>(ls->content);
    }



    SUNLinearSolver_Type
    arkode_linsol_get_type(SUNLinearSolver)
    {
      return SUNLINEARSOLVER_ITERATIVE;
    }



    template <typename VectorType>
    int
    arkode_linsol_solve(SUNLinearSolver LS,
                        SUNMatrix /*ignored*/,
                        N_Vector           x,
                        N_Vector           b,
                        SUNDIALS::realtype tol)
    {
      LinearSolverContent<VectorType> *content = access_content<VectorType>(LS);

      auto *src_b = unwrap_nvector_const<VectorType>(b);
      auto *dst_x = unwrap_nvector<VectorType>(x);

      SundialsOperator<VectorType> op(content->A_data,
                                      content->a_times_fn
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                                      ,
                                      content->linsol_ctx
#  endif
      );

      SundialsPreconditioner<VectorType> preconditioner(
        content->P_data,
        content->preconditioner_solve,
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
        content->linsol_ctx,
#  endif
        tol);

      return call_and_possibly_capture_exception(content->lsolve,
                                                 content->pending_exception,
                                                 op,
                                                 preconditioner,
                                                 *dst_x,
                                                 *src_b,
                                                 tol);
    }



    template <typename VectorType>
    int
    arkode_linsol_setup(SUNLinearSolver LS, SUNMatrix /*ignored*/)
    {
      LinearSolverContent<VectorType> *content = access_content<VectorType>(LS);

      if (content->preconditioner_setup)
        return content->preconditioner_setup(content->P_data);
      return 0;
    }



    template <typename VectorType>
    int
    arkode_linsol_initialize(SUNLinearSolver)
    {
      // this method is currently only provided because SUNDIALS 4.0.0 requires
      // it - no user-set action is implemented so far
      return 0;
    }


    template <typename VectorType>
    int
    arkode_linsol_set_a_times(SUNLinearSolver LS,
                              void           *A_data,
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                              SUNATimesFn ATimes
#  else
                              ATimesFn        ATimes
#  endif
    )
    {
      LinearSolverContent<VectorType> *content = access_content<VectorType>(LS);

      content->A_data     = A_data;
      content->a_times_fn = ATimes;
      return 0;
    }



    template <typename VectorType>
    int
    arkode_linsol_set_preconditioner(SUNLinearSolver LS,
                                     void           *P_data,
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                                     SUNPSetupFn p_setup,
                                     SUNPSolveFn p_solve
#  else
                                     PSetupFn p_setup,
                                     PSolveFn p_solve
#  endif
    )
    {
      LinearSolverContent<VectorType> *content = access_content<VectorType>(LS);

      content->P_data               = P_data;
      content->preconditioner_setup = p_setup;
      content->preconditioner_solve = p_solve;
      return 0;
    }
  } // namespace



  template <typename VectorType>
  internal::LinearSolverWrapper<VectorType>::LinearSolverWrapper(
    const LinearSolveFunction<VectorType> &lsolve,
    std::exception_ptr                    &pending_exception
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    ,
    SUNContext linsol_ctx
#  endif
    )
    : content(
        std::make_unique<LinearSolverContent<VectorType>>(pending_exception))
  {
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    sun_linear_solver = SUNLinSolNewEmpty(linsol_ctx);
#  else
    sun_linear_solver = SUNLinSolNewEmpty();
#  endif

    sun_linear_solver->ops->gettype    = arkode_linsol_get_type;
    sun_linear_solver->ops->solve      = arkode_linsol_solve<VectorType>;
    sun_linear_solver->ops->setup      = arkode_linsol_setup<VectorType>;
    sun_linear_solver->ops->initialize = arkode_linsol_initialize<VectorType>;
    sun_linear_solver->ops->setatimes  = arkode_linsol_set_a_times<VectorType>;
    sun_linear_solver->ops->setpreconditioner =
      arkode_linsol_set_preconditioner<VectorType>;

    content->lsolve = lsolve;
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    content->linsol_ctx = linsol_ctx;
#  endif
    sun_linear_solver->content = content.get();
  }



  namespace internal
  {
    template <typename VectorType>
    LinearSolverWrapper<VectorType>::~LinearSolverWrapper()
    {
      SUNLinSolFreeEmpty(sun_linear_solver);
    }



    template <typename VectorType>
    LinearSolverWrapper<VectorType>::operator SUNLinearSolver()
    {
      return sun_linear_solver;
    }
  } // namespace internal



#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
  template <typename VectorType>
  SundialsOperator<VectorType>::SundialsOperator(void       *A_data,
                                                 SUNATimesFn a_times_fn,
                                                 SUNContext  linsol_ctx)
    : A_data(A_data)
    , a_times_fn(a_times_fn)
    , linsol_ctx(linsol_ctx)

  {
    Assert(a_times_fn != nullptr, ExcInternalError());
    Assert(linsol_ctx != nullptr, ExcInternalError());
  }
#  else
  template <typename VectorType>
  SundialsOperator<VectorType>::SundialsOperator(void    *A_data,
                                                 ATimesFn a_times_fn)
    : A_data(A_data)
    , a_times_fn(a_times_fn)

  {
    Assert(a_times_fn != nullptr, ExcInternalError());
  }
#  endif



  template <typename VectorType>
  void
  SundialsOperator<VectorType>::vmult(VectorType       &dst,
                                      const VectorType &src) const
  {
    auto sun_dst = internal::make_nvector_view(dst
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                                               ,
                                               linsol_ctx
#  endif
    );
    auto sun_src = internal::make_nvector_view(src
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                                               ,
                                               linsol_ctx
#  endif
    );
    int status = a_times_fn(A_data, sun_src, sun_dst);
    (void)status;
    AssertSundialsSolver(status);
  }



#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
  template <typename VectorType>
  SundialsPreconditioner<VectorType>::SundialsPreconditioner(
    void       *P_data,
    SUNPSolveFn p_solve_fn,
    SUNContext  linsol_ctx,
    double      tol)
    : P_data(P_data)
    , p_solve_fn(p_solve_fn)
    , linsol_ctx(linsol_ctx)
    , tol(tol)
  {}
#  else
  template <typename VectorType>
  SundialsPreconditioner<VectorType>::SundialsPreconditioner(
    void    *P_data,
    PSolveFn p_solve_fn,
    double   tol)
    : P_data(P_data)
    , p_solve_fn(p_solve_fn)
    , tol(tol)
  {}
#  endif



  template <typename VectorType>
  void
  SundialsPreconditioner<VectorType>::vmult(VectorType       &dst,
                                            const VectorType &src) const
  {
    // apply identity preconditioner if nothing else specified
    if (!p_solve_fn)
      {
        dst = src;
        return;
      }

    auto sun_dst = internal::make_nvector_view(dst
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                                               ,
                                               linsol_ctx
#  endif
    );
    auto sun_src = internal::make_nvector_view(src
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                                               ,
                                               linsol_ctx
#  endif
    );
    // for custom preconditioners no distinction between left and right
    // preconditioning is made
    int status =
      p_solve_fn(P_data, sun_src, sun_dst, tol, 0 /*precondition_type*/);
    (void)status;
    AssertSundialsSolver(status);
  }

  template struct SundialsOperator<Vector<double>>;
  template struct SundialsOperator<BlockVector<double>>;
  template struct SundialsOperator<LinearAlgebra::distributed::Vector<double>>;
  template struct SundialsOperator<
    LinearAlgebra::distributed::BlockVector<double>>;

  template struct SundialsPreconditioner<Vector<double>>;
  template struct SundialsPreconditioner<BlockVector<double>>;
  template struct SundialsPreconditioner<
    LinearAlgebra::distributed::Vector<double>>;
  template struct SundialsPreconditioner<
    LinearAlgebra::distributed::BlockVector<double>>;

  template class internal::LinearSolverWrapper<Vector<double>>;
  template class internal::LinearSolverWrapper<BlockVector<double>>;
  template class internal::LinearSolverWrapper<
    LinearAlgebra::distributed::Vector<double>>;
  template class internal::LinearSolverWrapper<
    LinearAlgebra::distributed::BlockVector<double>>;

#  ifdef DEAL_II_WITH_MPI
#    ifdef DEAL_II_WITH_TRILINOS
  template struct SundialsOperator<TrilinosWrappers::MPI::Vector>;
  template struct SundialsOperator<TrilinosWrappers::MPI::BlockVector>;

  template struct SundialsPreconditioner<TrilinosWrappers::MPI::Vector>;
  template struct SundialsPreconditioner<TrilinosWrappers::MPI::BlockVector>;

  template class internal::LinearSolverWrapper<TrilinosWrappers::MPI::Vector>;
  template class internal::LinearSolverWrapper<
    TrilinosWrappers::MPI::BlockVector>;
#    endif // DEAL_II_WITH_TRILINOS

#    ifdef DEAL_II_WITH_PETSC
#      ifndef PETSC_USE_COMPLEX

  template struct SundialsOperator<PETScWrappers::MPI::Vector>;
  template struct SundialsOperator<PETScWrappers::MPI::BlockVector>;

  template struct SundialsPreconditioner<PETScWrappers::MPI::Vector>;
  template struct SundialsPreconditioner<PETScWrappers::MPI::BlockVector>;

  template class internal::LinearSolverWrapper<PETScWrappers::MPI::Vector>;
  template class internal::LinearSolverWrapper<PETScWrappers::MPI::BlockVector>;
#      endif // PETSC_USE_COMPLEX
#    endif   // DEAL_II_WITH_PETSC

#  endif // DEAL_II_WITH_MPI
} // namespace SUNDIALS
DEAL_II_NAMESPACE_CLOSE

#endif
