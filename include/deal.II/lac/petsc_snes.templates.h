// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_petsc_snes_templates_h
#define dealii_petsc_snes_templates_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PETSC
#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/mpi_stub.h>

#  include <deal.II/lac/petsc_precondition.h>
#  include <deal.II/lac/petsc_snes.h>

#  include <petscdm.h>
#  include <petscerror.h>
#  include <petscsnes.h>

DEAL_II_NAMESPACE_OPEN

// Shorthand notation for PETSc error codes.
// This is used in deal.II code to raise exceptions
#  define AssertPETSc(code)                          \
    do                                               \
      {                                              \
        PetscErrorCode ierr = (code);                \
        AssertThrow(ierr == 0, ExcPETScError(ierr)); \
      }                                              \
    while (false)

// Macro to wrap PETSc inside callbacks.
// This is used to raise "PETSc" exceptions, i.e.
// start a cascade of errors inside PETSc
#  ifndef PetscCall
#    define PetscCall(code)             \
      do                                \
        {                               \
          PetscErrorCode ierr = (code); \
          CHKERRQ(ierr);                \
        }                               \
      while (false)
#    define undefPetscCall
#  endif

namespace PETScWrappers
{
  /**
   * A function that calls the function object given by its first argument
   * with the set of arguments following at the end. If the call returns
   * regularly, the current function returns zero to indicate success. If
   * the call fails with an exception, then the current function returns with
   * an error code of -1. In that case, the exception thrown by @p f is
   * captured and @p eptr is set to the exception. In case of success,
   * @p eptr is set to `nullptr`.
   *
   * If the user callback fails with a recoverable exception, then (i) if
   * @p recoverable_action is set, execute it, eat the exception, and return
   * zero; (ii) if @p recoverable_action is an empty function object, store the
   * exception and return -1.
   */
  template <typename F, typename... Args>
  int
  call_and_possibly_capture_snes_exception(
    const F                     &f,
    std::exception_ptr          &eptr,
    const std::function<void()> &recoverable_action,
    Args &&...args)
  {
    // See whether there is already something in the exception pointer
    // variable. There is no reason why this should be so, and
    // we should probably bail out:
    AssertThrow(eptr == nullptr, ExcInternalError());

    // Call the function and if that succeeds, return zero:
    try
      {
        f(std::forward<Args>(args)...);
        eptr = nullptr;
        return 0;
      }
    // In case of a recoverable exception call the action if present:
    catch (const RecoverableUserCallbackError &)
      {
        if (recoverable_action)
          {
            // recover and eat exception
            recoverable_action();
            eptr = nullptr;
            return 0;
          }
        else
          {
            // No action provided.
            // capture exception and return -1
            eptr = std::current_exception();
            return -1;
          }
      }
    // In case of an exception, capture the exception and
    // return -1:
    catch (...)
      {
        eptr = std::current_exception();
        return -1;
      }
  }



  template <typename VectorType, typename PMatrixType, typename AMatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<PMatrixType> ||
               std::constructible_from<
                 PMatrixType,
                 Mat>)&&(concepts::is_dealii_petsc_matrix_type<AMatrixType> ||
                         std::constructible_from<AMatrixType, Mat>))
  NonlinearSolver<VectorType, PMatrixType, AMatrixType>::NonlinearSolver(
    const NonlinearSolverData &data,
    const MPI_Comm             mpi_comm)
    : pending_exception(nullptr)
  {
    AssertPETSc(SNESCreate(mpi_comm, &snes));
    AssertPETSc(SNESSetApplicationContext(snes, this));
    reinit(data);
  }



  template <typename VectorType, typename PMatrixType, typename AMatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<PMatrixType> ||
               std::constructible_from<
                 PMatrixType,
                 Mat>)&&(concepts::is_dealii_petsc_matrix_type<AMatrixType> ||
                         std::constructible_from<AMatrixType, Mat>))
  NonlinearSolver<VectorType, PMatrixType, AMatrixType>::operator SNES() const
  {
    return snes;
  }



  template <typename VectorType, typename PMatrixType, typename AMatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<PMatrixType> ||
               std::constructible_from<
                 PMatrixType,
                 Mat>)&&(concepts::is_dealii_petsc_matrix_type<AMatrixType> ||
                         std::constructible_from<AMatrixType, Mat>))
  SNES NonlinearSolver<VectorType, PMatrixType, AMatrixType>::petsc_snes()
  {
    return snes;
  }



  template <typename VectorType, typename PMatrixType, typename AMatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<PMatrixType> ||
               std::constructible_from<
                 PMatrixType,
                 Mat>)&&(concepts::is_dealii_petsc_matrix_type<AMatrixType> ||
                         std::constructible_from<AMatrixType, Mat>))
  inline MPI_Comm NonlinearSolver<VectorType, PMatrixType, AMatrixType>::
    get_mpi_communicator() const
  {
    return PetscObjectComm(reinterpret_cast<PetscObject>(snes));
  }



  template <typename VectorType, typename PMatrixType, typename AMatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<PMatrixType> ||
               std::constructible_from<
                 PMatrixType,
                 Mat>)&&(concepts::is_dealii_petsc_matrix_type<AMatrixType> ||
                         std::constructible_from<AMatrixType, Mat>))
  NonlinearSolver<VectorType, PMatrixType, AMatrixType>::~NonlinearSolver()
  {
    AssertPETSc(SNESDestroy(&snes));

    Assert(pending_exception == nullptr, ExcInternalError());
  }



  template <typename VectorType, typename PMatrixType, typename AMatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<PMatrixType> ||
               std::constructible_from<
                 PMatrixType,
                 Mat>)&&(concepts::is_dealii_petsc_matrix_type<AMatrixType> ||
                         std::constructible_from<AMatrixType, Mat>))
  void NonlinearSolver<VectorType, PMatrixType, AMatrixType>::reinit()
  {
    AssertPETSc(SNESReset(snes));
    // By default PETSc does not check for Jacobian errors in optimized
    // mode. Here we do it unconditionally.
#  if DEAL_II_PETSC_VERSION_GTE(3, 11, 0)
    AssertPETSc(SNESSetCheckJacobianDomainError(snes, PETSC_TRUE));
#  endif
  }



  template <typename VectorType, typename PMatrixType, typename AMatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<PMatrixType> ||
               std::constructible_from<
                 PMatrixType,
                 Mat>)&&(concepts::is_dealii_petsc_matrix_type<AMatrixType> ||
                         std::constructible_from<AMatrixType, Mat>))
  void NonlinearSolver<VectorType, PMatrixType, AMatrixType>::reinit(
    const NonlinearSolverData &data)
  {
    reinit();

    // Solver type
    if (data.snes_type.size())
      AssertPETSc(SNESSetType(snes, data.snes_type.c_str()));

    // Linesearch type
    if (data.snes_linesearch_type.size())
      {
        SNESLineSearch linesearch;
        AssertPETSc(SNESGetLineSearch(snes, &linesearch));
        AssertPETSc(
          SNESLineSearchSetType(linesearch, data.snes_linesearch_type.c_str()));
      }

    // Options prefix
    if (data.options_prefix.size())
      AssertPETSc(SNESSetOptionsPrefix(snes, data.options_prefix.c_str()));

    // Solver tolerances
    PetscReal atol = data.absolute_tolerance > 0.0 ?
                       data.absolute_tolerance :
                       static_cast<PetscReal>(PETSC_DEFAULT);
    PetscReal rtol = data.relative_tolerance > 0.0 ?
                       data.relative_tolerance :
                       static_cast<PetscReal>(PETSC_DEFAULT);
    PetscReal stol = data.step_tolerance > 0.0 ?
                       data.step_tolerance :
                       static_cast<PetscReal>(PETSC_DEFAULT);

    // Maximum number of iterations and function evaluations.
    PetscInt maxit = data.maximum_non_linear_iterations >= 0 ?
                       data.maximum_non_linear_iterations :
                       static_cast<PetscInt>(PETSC_DEFAULT);
    PetscInt maxfe = data.max_n_function_evaluations >= 0 ?
                       data.max_n_function_evaluations :
                       static_cast<PetscInt>(PETSC_DEFAULT);
    AssertPETSc(SNESSetTolerances(snes, atol, rtol, stol, maxit, maxfe));
  }



  template <typename VectorType, typename PMatrixType, typename AMatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<PMatrixType> ||
               std::constructible_from<
                 PMatrixType,
                 Mat>)&&(concepts::is_dealii_petsc_matrix_type<AMatrixType> ||
                         std::constructible_from<AMatrixType, Mat>))
  void NonlinearSolver<VectorType, PMatrixType, AMatrixType>::set_matrix(
    PMatrixType &P)
  {
    this->A = nullptr;
    this->P = &P;
  }



  template <typename VectorType, typename PMatrixType, typename AMatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<PMatrixType> ||
               std::constructible_from<
                 PMatrixType,
                 Mat>)&&(concepts::is_dealii_petsc_matrix_type<AMatrixType> ||
                         std::constructible_from<AMatrixType, Mat>))
  void NonlinearSolver<VectorType, PMatrixType, AMatrixType>::set_matrices(
    AMatrixType &A,
    PMatrixType &P)
  {
    this->A = &A;
    this->P = &P;
  }



  template <typename VectorType, typename PMatrixType, typename AMatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<PMatrixType> ||
               std::constructible_from<
                 PMatrixType,
                 Mat>)&&(concepts::is_dealii_petsc_matrix_type<AMatrixType> ||
                         std::constructible_from<AMatrixType, Mat>))
  unsigned int NonlinearSolver<VectorType, PMatrixType, AMatrixType>::solve(
    VectorType &x)
  {
    const auto snes_function =
      [](SNES snes, Vec x, Vec f, void *ctx) -> PetscErrorCode {
      PetscFunctionBeginUser;
      auto user = static_cast<NonlinearSolver *>(ctx);

      VectorType xdealii(x);
      VectorType fdealii(f);
      const int  lineno = __LINE__;
      const int  err    = call_and_possibly_capture_snes_exception(
        user->residual,
        user->pending_exception,
        [snes]() -> void { AssertPETSc(SNESSetFunctionDomainError(snes)); },
        xdealii,
        fdealii);
      if (err)
        return PetscError(
          PetscObjectComm((PetscObject)snes),
          lineno + 1,
          "residual",
          __FILE__,
          PETSC_ERR_LIB,
          PETSC_ERROR_INITIAL,
          "Failure in snes_function from dealii::PETScWrappers::NonlinearSolver");
      petsc_increment_state_counter(f);
      PetscFunctionReturn(PETSC_SUCCESS);
    };

    const auto snes_jacobian =
      [](SNES snes, Vec x, Mat A, Mat P, void *ctx) -> PetscErrorCode {
      PetscFunctionBeginUser;
      auto user = static_cast<NonlinearSolver *>(ctx);

      VectorType  xdealii(x);
      AMatrixType Adealii(A);
      PMatrixType Pdealii(P);

      const int lineno = __LINE__;
      const int err    = call_and_possibly_capture_snes_exception(
        user->jacobian,
        user->pending_exception,
        [snes]() -> void { snes_set_jacobian_domain_error(snes); },
        xdealii,
        Adealii,
        Pdealii);
      if (err)
        return PetscError(
          PetscObjectComm((PetscObject)snes),
          lineno + 1,
          "jacobian",
          __FILE__,
          PETSC_ERR_LIB,
          PETSC_ERROR_INITIAL,
          "Failure in snes_jacobian from dealii::PETScWrappers::NonlinearSolver");
      petsc_increment_state_counter(P);

      // Handle the Jacobian-free case
      // This call allows to resample the linearization point
      // of the MFFD tangent operator
      PetscBool flg;
      PetscCall(PetscObjectTypeCompare((PetscObject)A, MATMFFD, &flg));
      if (flg)
        {
          PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
          PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
        }
      else
        petsc_increment_state_counter(A);

      PetscFunctionReturn(PETSC_SUCCESS);
    };

    const auto snes_jacobian_with_setup =
      [](SNES snes, Vec x, Mat A, Mat P, void *ctx) -> PetscErrorCode {
      PetscFunctionBeginUser;
      auto user = static_cast<NonlinearSolver *>(ctx);

      VectorType xdealii(x);

      const int lineno = __LINE__;
      const int err    = call_and_possibly_capture_snes_exception(
        user->setup_jacobian,
        user->pending_exception,
        [snes]() -> void { snes_set_jacobian_domain_error(snes); },
        xdealii);
      if (err)
        return PetscError(
          PetscObjectComm((PetscObject)snes),
          lineno + 1,
          "setup_jacobian",
          __FILE__,
          PETSC_ERR_LIB,
          PETSC_ERROR_INITIAL,
          "Failure in snes_jacobian_with_setup from dealii::PETScWrappers::NonlinearSolver");
      // The MatCopy calls below are 99% of the times dummy calls.
      // They are only used in case we get different Mats then the one we passed
      // to SNESSetJacobian.
      if (user->P)
        PetscCall(MatCopy(user->P->petsc_matrix(), P, SAME_NONZERO_PATTERN));
      if (user->A)
        PetscCall(MatCopy(user->A->petsc_matrix(), A, SAME_NONZERO_PATTERN));
      petsc_increment_state_counter(P);

      // Handle older versions of PETSc for which we cannot pass a MATSHELL
      // matrix to DMSetMatType. This has been fixed from 3.13 on.
      if (user->need_dummy_assemble)
        {
          PetscCall(MatZeroEntries(P));
          PetscCall(MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY));
          PetscCall(MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY));
        }
      user->need_dummy_assemble = false;

      // Handle the Jacobian-free case
      // This call allows to resample the linearization point
      // of the MFFD tangent operator
      PetscBool flg;
      PetscCall(PetscObjectTypeCompare((PetscObject)A, MATMFFD, &flg));
      if (flg)
        {
          PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
          PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
        }
      else
        petsc_increment_state_counter(A);
      PetscFunctionReturn(PETSC_SUCCESS);
    };

    const auto snes_monitor =
      [](SNES snes, PetscInt it, PetscReal f, void *ctx) -> PetscErrorCode {
      PetscFunctionBeginUser;
      auto user = static_cast<NonlinearSolver *>(ctx);

      Vec x;
      PetscCall(SNESGetSolution(snes, &x));
      VectorType xdealii(x);
      const int  lineno = __LINE__;
      const int  err    = call_and_possibly_capture_snes_exception(
        user->monitor, user->pending_exception, {}, xdealii, it, f);
      if (err)
        return PetscError(
          PetscObjectComm((PetscObject)snes),
          lineno + 1,
          "monitor",
          __FILE__,
          PETSC_ERR_LIB,
          PETSC_ERROR_INITIAL,
          "Failure in snes_monitor from dealii::PETScWrappers::NonlinearSolver");
      PetscFunctionReturn(PETSC_SUCCESS);
    };

    const auto snes_objective =
      [](SNES snes, Vec x, PetscReal *f, void *ctx) -> PetscErrorCode {
      PetscFunctionBeginUser;
      auto user = static_cast<NonlinearSolver *>(ctx);

      real_type  v;
      VectorType xdealii(x);
      const int  lineno = __LINE__;
      const int  err    = call_and_possibly_capture_snes_exception(
        user->energy,
        user->pending_exception,
        [snes, &v]() -> void {
          AssertPETSc(SNESSetFunctionDomainError(snes));
          v = std::numeric_limits<real_type>::quiet_NaN();
        },
        xdealii,
        v);
      if (err)
        return PetscError(
          PetscObjectComm((PetscObject)snes),
          lineno + 1,
          "energy",
          __FILE__,
          PETSC_ERR_LIB,
          PETSC_ERROR_INITIAL,
          "Failure in snes_objective from dealii::PETScWrappers::NonlinearSolver");
      *f = v;
      PetscFunctionReturn(PETSC_SUCCESS);
    };

    AssertThrow(residual,
                StandardExceptions::ExcFunctionNotProvided("residual"));

    AssertPETSc(SNESSetSolution(snes, x.petsc_vector()));

    AssertPETSc(SNESSetFunction(snes, nullptr, snes_function, this));

    if (energy)
      AssertPETSc(SNESSetObjective(snes, snes_objective, this));

    this->need_dummy_assemble = false;
    if (setup_jacobian)
      {
        AssertPETSc(SNESSetJacobian(snes,
                                    A ? A->petsc_matrix() : nullptr,
                                    P ? P->petsc_matrix() : nullptr,
                                    snes_jacobian_with_setup,
                                    this));

        // Tell PETSc to set up a MFFD operator for the linear system matrix
        if (!A)
          set_use_matrix_free(snes, true, false);

        // Do not waste memory by creating a dummy AIJ matrix inside PETSc.
        if (!P)
          {
#  if DEAL_II_PETSC_VERSION_GTE(3, 13, 0)
            DM dm;
            AssertPETSc(SNESGetDM(snes, &dm));
            AssertPETSc(DMSetMatType(dm, MATSHELL));
#  else
            this->need_dummy_assemble = true;
#  endif
          }
      }
    else
      {
        if (jacobian)
          {
            AssertPETSc(SNESSetJacobian(snes,
                                        A ? A->petsc_matrix() :
                                            (P ? P->petsc_matrix() : nullptr),
                                        P ? P->petsc_matrix() : nullptr,
                                        snes_jacobian,
                                        this));
          }
        else
          // The user did not set any Jacobian callback. PETSc default in this
          // case is to use FD and thus assemble a dense operator by finite
          // differencing the residual callbacks. Here instead we decide to
          // use a full matrix-free approach by default. This choice can always
          // be overridden from command line.
          {
            set_use_matrix_free(snes, false, true);
          }
      }

    // In case solve_with_jacobian is provided, create a shell
    // preconditioner wrapping the user call. The default internal Krylov
    // solver only applies the preconditioner. This choice
    // can be overridden by command line and users can use any other
    // Krylov method if their solve is not accurate enough.
    // Using solve_with_jacobian as a preconditioner allows users
    // to provide approximate solvers and possibly iterate on a matrix-free
    // approximation of the tangent operator.
    PreconditionShell precond(
      PetscObjectComm(reinterpret_cast<PetscObject>(snes)));
    if (solve_with_jacobian)
      {
        precond.vmult = [&](VectorBase       &indst,
                            const VectorBase &insrc) -> void {
          VectorType       dst(static_cast<const Vec &>(indst));
          const VectorType src(static_cast<const Vec &>(insrc));
          solve_with_jacobian(src, dst);
        };

        // Default Krylov solver (preconditioner only)
        KSP ksp;
        AssertPETSc(SNESGetKSP(snes, &ksp));
        AssertPETSc(KSPSetType(ksp, KSPPREONLY));
        AssertPETSc(KSPSetPC(ksp, precond.get_pc()));
      }

    // Attach user monitoring routine
    if (monitor)
      AssertPETSc(SNESMonitorSet(snes, snes_monitor, this, nullptr));

    // Allow command line customization.
    AssertPETSc(SNESSetFromOptions(snes));

    // Having set everything up, now do the actual work
    // and let PETSc solve the system.
    // Older versions of PETSc requires the solution vector specified even
    // if we specified SNESSetSolution upfront.
    //
    // If there is a pending exception, then one of the user callbacks
    // threw an exception we didn't know how to deal with
    // at the time. It is possible that PETSc managed to
    // recover anyway and in that case would have returned
    // a zero error code -- if so, just eat the exception and
    // continue on; otherwise, just rethrow the exception
    // and get outta here.
    const PetscErrorCode status = SNESSolve(snes, nullptr, x.petsc_vector());
    if (pending_exception)
      {
        try
          {
            std::rethrow_exception(pending_exception);
          }
        catch (...)
          {
            pending_exception = nullptr;
            if (status == 0)
              /* just eat the exception */;
            else
              throw;
          }
      }
    AssertPETSc(status);

    // Get the number of steps taken.
    PetscInt n_iterations;
    AssertPETSc(SNESGetIterationNumber(snes, &n_iterations));

    // Raise an exception if the solver has not converged
    SNESConvergedReason reason;
    AssertPETSc(SNESGetConvergedReason(snes, &reason));
    AssertThrow(reason > 0,
                ExcMessage("SNES solver did not converge after " +
                           std::to_string(n_iterations) +
                           " iterations with reason " +
                           SNESConvergedReasons[reason]));

    // Finally return
    return n_iterations;
  }



  template <typename VectorType, typename PMatrixType, typename AMatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<PMatrixType> ||
               std::constructible_from<
                 PMatrixType,
                 Mat>)&&(concepts::is_dealii_petsc_matrix_type<AMatrixType> ||
                         std::constructible_from<AMatrixType, Mat>))
  unsigned int NonlinearSolver<VectorType, PMatrixType, AMatrixType>::solve(
    VectorType  &x,
    PMatrixType &P)
  {
    set_matrix(P);
    return solve(x);
  }



  template <typename VectorType, typename PMatrixType, typename AMatrixType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_petsc_vector_type<VectorType> ||
     std::constructible_from<
       VectorType,
       Vec>)&&(concepts::is_dealii_petsc_matrix_type<PMatrixType> ||
               std::constructible_from<
                 PMatrixType,
                 Mat>)&&(concepts::is_dealii_petsc_matrix_type<AMatrixType> ||
                         std::constructible_from<AMatrixType, Mat>))
  unsigned int NonlinearSolver<VectorType, PMatrixType, AMatrixType>::solve(
    VectorType  &x,
    AMatrixType &A,
    PMatrixType &P)
  {
    set_matrices(A, P);
    return solve(x);
  }

} // namespace PETScWrappers

#  undef AssertPETSc
#  if defined(undefPetscCall)
#    undef PetscCall
#    undef undefPetscCall
#  endif

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC

#endif
