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

#include <deal.II/base/config.h>

#include <deal.II/sundials/sunlinsol_wrapper.h>

#ifdef DEAL_II_WITH_SUNDIALS
#  if DEAL_II_SUNDIALS_VERSION_GTE(4, 0, 0)

#    include <deal.II/base/exceptions.h>

#    include <deal.II/lac/block_vector.h>
#    include <deal.II/lac/la_parallel_block_vector.h>
#    include <deal.II/lac/la_parallel_vector.h>
#    include <deal.II/lac/vector.h>
#    ifdef DEAL_II_WITH_TRILINOS
#      include <deal.II/lac/trilinos_parallel_block_vector.h>
#      include <deal.II/lac/trilinos_vector.h>
#    endif
#    ifdef DEAL_II_WITH_PETSC
#      include <deal.II/lac/petsc_block_vector.h>
#      include <deal.II/lac/petsc_vector.h>
#    endif

#    include <deal.II/sundials/n_vector.h>

#    if DEAL_II_SUNDIALS_VERSION_LT(5, 0, 0)
#      include <deal.II/sundials/sunlinsol_newempty.h>
#    endif

DEAL_II_NAMESPACE_OPEN



namespace SUNDIALS
{
  DeclException1(ExcSundialsSolverError,
                 int,
                 << "One of the SUNDIALS linear solver internal"
                 << " functions returned a negative error code: " << arg1
                 << ". Please consult SUNDIALS manual.");

#    define AssertSundialsSolver(code) \
      Assert(code >= 0, ExcSundialsSolverError(code))

  namespace internal
  {
    /**
     * storage for internal content of the linear solver wrapper
     */
    template <typename VectorType>
    struct LinearSolverContent
    {
      ATimesFn a_times_fn;
      PSetupFn preconditioner_setup;
      PSolveFn preconditioner_solve;

      LinearSolveFunction<VectorType> lsolve;

      void *P_data;
      void *A_data;
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



    SUNLinearSolver_Type arkode_linsol_get_type(SUNLinearSolver)
    {
      return SUNLINEARSOLVER_ITERATIVE;
    }



    template <typename VectorType>
    int
    arkode_linsol_solve(SUNLinearSolver LS,
                        SUNMatrix /*ignored*/,
                        N_Vector x,
                        N_Vector b,
                        realtype tol)
    {
      auto content = access_content<VectorType>(LS);

      auto *src_b = unwrap_nvector_const<VectorType>(b);
      auto *dst_x = unwrap_nvector<VectorType>(x);

      SundialsOperator<VectorType> op(content->A_data, content->a_times_fn);

      SundialsPreconditioner<VectorType> preconditioner(
        content->P_data, content->preconditioner_solve, tol);

      return content->lsolve(op, preconditioner, *dst_x, *src_b, tol);
    }



    template <typename VectorType>
    int
    arkode_linsol_setup(SUNLinearSolver LS, SUNMatrix /*ignored*/)
    {
      auto content = access_content<VectorType>(LS);
      if (content->preconditioner_setup)
        return content->preconditioner_setup(content->P_data);
      return 0;
    }



    template <typename VectorType>
    int arkode_linsol_initialize(SUNLinearSolver)
    {
      // this method is currently only provided because SUNDIALS 4.0.0 requires
      // it - no user-set action is implemented so far
      return 0;
    }



    template <typename VectorType>
    int
    arkode_linsol_set_a_times(SUNLinearSolver LS, void *A_data, ATimesFn ATimes)
    {
      auto content        = access_content<VectorType>(LS);
      content->A_data     = A_data;
      content->a_times_fn = ATimes;
      return 0;
    }



    template <typename VectorType>
    int
    arkode_linsol_set_preconditioner(SUNLinearSolver LS,
                                     void *          P_data,
                                     PSetupFn        p_setup,
                                     PSolveFn        p_solve)
    {
      auto content                  = access_content<VectorType>(LS);
      content->P_data               = P_data;
      content->preconditioner_setup = p_setup;
      content->preconditioner_solve = p_solve;
      return 0;
    }
  } // namespace



  template <typename VectorType>
  internal::LinearSolverWrapper<VectorType>::LinearSolverWrapper(
    LinearSolveFunction<VectorType> lsolve)
    : content(std::make_unique<LinearSolverContent<VectorType>>())
  {
    sun_linear_solver                  = SUNLinSolNewEmpty();
    sun_linear_solver->ops->gettype    = arkode_linsol_get_type;
    sun_linear_solver->ops->solve      = arkode_linsol_solve<VectorType>;
    sun_linear_solver->ops->setup      = arkode_linsol_setup<VectorType>;
    sun_linear_solver->ops->initialize = arkode_linsol_initialize<VectorType>;
    sun_linear_solver->ops->setatimes  = arkode_linsol_set_a_times<VectorType>;
    sun_linear_solver->ops->setpreconditioner =
      arkode_linsol_set_preconditioner<VectorType>;

    content->lsolve            = lsolve;
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



  template <typename VectorType>
  SundialsOperator<VectorType>::SundialsOperator(void *   A_data,
                                                 ATimesFn a_times_fn)
    : A_data(A_data)
    , a_times_fn(a_times_fn)

  {
    Assert(a_times_fn != nullptr, ExcInternalError());
  }



  template <typename VectorType>
  void
  SundialsOperator<VectorType>::vmult(VectorType &      dst,
                                      const VectorType &src) const
  {
    auto sun_dst = internal::make_nvector_view(dst);
    auto sun_src = internal::make_nvector_view(src);
    int  status  = a_times_fn(A_data, sun_src, sun_dst);
    (void)status;
    AssertSundialsSolver(status);
  }



  template <typename VectorType>
  SundialsPreconditioner<VectorType>::SundialsPreconditioner(
    void *   P_data,
    PSolveFn p_solve_fn,
    double   tol)
    : P_data(P_data)
    , p_solve_fn(p_solve_fn)
    , tol(tol)
  {}



  template <typename VectorType>
  void
  SundialsPreconditioner<VectorType>::vmult(VectorType &      dst,
                                            const VectorType &src) const
  {
    // apply identity preconditioner if nothing else specified
    if (!p_solve_fn)
      {
        dst = src;
        return;
      }

    auto sun_dst = internal::make_nvector_view(dst);
    auto sun_src = internal::make_nvector_view(src);
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

#    ifdef DEAL_II_WITH_MPI
#      ifdef DEAL_II_WITH_TRILINOS
  template struct SundialsOperator<TrilinosWrappers::MPI::Vector>;
  template struct SundialsOperator<TrilinosWrappers::MPI::BlockVector>;

  template struct SundialsPreconditioner<TrilinosWrappers::MPI::Vector>;
  template struct SundialsPreconditioner<TrilinosWrappers::MPI::BlockVector>;

  template class internal::LinearSolverWrapper<TrilinosWrappers::MPI::Vector>;
  template class internal::LinearSolverWrapper<
    TrilinosWrappers::MPI::BlockVector>;
#      endif // DEAL_II_WITH_TRILINOS

#      ifdef DEAL_II_WITH_PETSC
#        ifndef PETSC_USE_COMPLEX

  template struct SundialsOperator<PETScWrappers::MPI::Vector>;
  template struct SundialsOperator<PETScWrappers::MPI::BlockVector>;

  template struct SundialsPreconditioner<PETScWrappers::MPI::Vector>;
  template struct SundialsPreconditioner<PETScWrappers::MPI::BlockVector>;

  template class internal::LinearSolverWrapper<PETScWrappers::MPI::Vector>;
  template class internal::LinearSolverWrapper<PETScWrappers::MPI::BlockVector>;
#        endif // PETSC_USE_COMPLEX
#      endif   // DEAL_II_WITH_PETSC

#    endif // DEAL_II_WITH_MPI
} // namespace SUNDIALS
DEAL_II_NAMESPACE_CLOSE

#  endif
#endif
