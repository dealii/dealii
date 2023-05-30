//-----------------------------------------------------------
//
//    Copyright (C) 2023 by the deal.II authors
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
//---------------------------------------------------------------

#ifndef dealii_petsc_ts_templates_h
#define dealii_petsc_ts_templates_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PETSC
#  include <deal.II/base/exceptions.h>

#  include <deal.II/lac/petsc_precondition.h>
#  include <deal.II/lac/petsc_ts.h>

#  include <petscdm.h>
#  include <petscts.h>

DEAL_II_NAMESPACE_OPEN

// Shorthand notation for PETSc error codes.
#  define AssertPETSc(code)                          \
    do                                               \
      {                                              \
        PetscErrorCode ierr = (code);                \
        AssertThrow(ierr == 0, ExcPETScError(ierr)); \
      }                                              \
    while (0)

namespace PETScWrappers
{
  namespace
  {
    /**
     * A function that calls the function object given by its first argument
     * with the set of arguments following at the end. If the call returns
     * regularly, the current function returns zero to indicate success. If
     * the call fails with an exception, then the current function returns with
     * an error code of -1. In that case, the exception thrown by `f` is
     * captured and `eptr` is set to the exception. In case of success,
     * `eptr` is set to `nullptr`.
     */
    template <typename F, typename... Args>
    int
    call_and_possibly_capture_exception(const F &           f,
                                        std::exception_ptr &eptr,
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
      // In case of an exception, capture the exception and
      // return -1:
      catch (...)
        {
          eptr = std::current_exception();
          return -1;
        }
    }
  } // namespace



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
  TimeStepper<VectorType, PMatrixType, AMatrixType>::TimeStepper(
    const TimeStepperData &data,
    const MPI_Comm         mpi_comm)
    : pending_exception(nullptr)
  {
    AssertPETSc(TSCreate(mpi_comm, &ts));
    AssertPETSc(TSSetApplicationContext(ts, this));
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
  TimeStepper<VectorType, PMatrixType, AMatrixType>::~TimeStepper()
  {
    AssertPETSc(TSDestroy(&ts));

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
  TimeStepper<VectorType, PMatrixType, AMatrixType>::operator TS() const
  {
    return ts;
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
  TS TimeStepper<VectorType, PMatrixType, AMatrixType>::petsc_ts()
  {
    return ts;
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
  inline MPI_Comm
    TimeStepper<VectorType, PMatrixType, AMatrixType>::get_mpi_communicator()
      const
  {
    return PetscObjectComm(reinterpret_cast<PetscObject>(ts));
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
  typename TimeStepper<VectorType, PMatrixType, AMatrixType>::real_type
    TimeStepper<VectorType, PMatrixType, AMatrixType>::get_time()
  {
    PetscReal      t;
    PetscErrorCode ierr = TSGetTime(ts, &t);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
    return t;
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
  typename TimeStepper<VectorType, PMatrixType, AMatrixType>::real_type
    TimeStepper<VectorType, PMatrixType, AMatrixType>::get_time_step()
  {
    PetscReal      dt;
    PetscErrorCode ierr = TSGetTimeStep(ts, &dt);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
    return dt;
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
  void TimeStepper<VectorType, PMatrixType, AMatrixType>::reinit()
  {
    AssertPETSc(TSReset(ts));
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
  void TimeStepper<VectorType, PMatrixType, AMatrixType>::reinit(
    const TimeStepperData &data)
  {
    reinit();

    // Solver type
    if (data.ts_type.size())
      AssertPETSc(TSSetType(ts, data.ts_type.c_str()));

    // Options prefix
    if (data.options_prefix.size())
      AssertPETSc(TSSetOptionsPrefix(ts, data.options_prefix.c_str()));

    // Time and steps limits
    AssertPETSc(TSSetTime(ts, data.initial_time));
    if (data.final_time > data.initial_time)
      ts_set_max_time(ts, data.final_time);
    if (data.initial_step_size > 0.0)
      AssertPETSc(TSSetTimeStep(ts, data.initial_step_size));
    if (data.max_steps >= 0)
      ts_set_max_steps(ts, data.max_steps);

    // Decide how to end the integration. Either stepover the final time or
    // match it.
    AssertPETSc(TSSetExactFinalTime(ts,
                                    data.match_step ?
                                      TS_EXACTFINALTIME_MATCHSTEP :
                                      TS_EXACTFINALTIME_STEPOVER));

    // Adaptive tolerances
    const PetscReal atol = data.absolute_tolerance > 0.0 ?
                             data.absolute_tolerance :
                             static_cast<PetscReal>(PETSC_DEFAULT);
    const PetscReal rtol = data.relative_tolerance > 0.0 ?
                             data.relative_tolerance :
                             static_cast<PetscReal>(PETSC_DEFAULT);
    AssertPETSc(TSSetTolerances(ts, atol, nullptr, rtol, nullptr));

    // At this point we do not know the problem size so we cannot
    // set variable tolerances for differential and algebratic equations
    // Store this value and use it during solve.
    this->need_dae_tolerances = data.ignore_algebraic_lte;

    // Adaptive time stepping
    TSAdapt tsadapt;
    AssertPETSc(TSGetAdapt(ts, &tsadapt));
    AssertPETSc(TSAdaptSetType(tsadapt, data.ts_adapt_type.c_str()));

    // As of 3.19, PETSc does not propagate options prefixes to the
    // adaptors.
    if (data.options_prefix.size())
      AssertPETSc(
        TSAdaptSetOptionsPrefix(tsadapt, data.options_prefix.c_str()));

    // Time step limits
    const PetscReal hmin = data.minimum_step_size > 0.0 ?
                             data.minimum_step_size :
                             static_cast<PetscReal>(PETSC_DEFAULT);
    const PetscReal hmax = data.maximum_step_size > 0.0 ?
                             data.maximum_step_size :
                             static_cast<PetscReal>(PETSC_DEFAULT);
    AssertPETSc(TSAdaptSetStepLimits(tsadapt, hmin, hmax));
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
  void TimeStepper<VectorType, PMatrixType, AMatrixType>::set_matrix(
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
  void TimeStepper<VectorType, PMatrixType, AMatrixType>::set_matrices(
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
  unsigned int TimeStepper<VectorType, PMatrixType, AMatrixType>::solve(
    VectorType &y)
  {
    const auto ts_ifunction =
      [](TS, PetscReal t, Vec x, Vec xdot, Vec f, void *ctx) -> PetscErrorCode {
      PetscFunctionBeginUser;
      auto user = static_cast<TimeStepper *>(ctx);

      VectorType xdealii(x);
      VectorType xdotdealii(xdot);
      VectorType fdealii(f);
      const int  err =
        call_and_possibly_capture_exception(user->implicit_function,
                                            user->pending_exception,
                                            t,
                                            xdealii,
                                            xdotdealii,
                                            fdealii);
      petsc_increment_state_counter(f);
      PetscFunctionReturn(err);
    };

    const auto ts_ijacobian =
      [](TS, PetscReal t, Vec x, Vec xdot, PetscReal s, Mat A, Mat P, void *ctx)
      -> PetscErrorCode {
      PetscFunctionBeginUser;
      auto user = static_cast<TimeStepper *>(ctx);

      VectorType  xdealii(x);
      VectorType  xdotdealii(xdot);
      AMatrixType Adealii(A);
      PMatrixType Pdealii(P);

      const int err =
        call_and_possibly_capture_exception(user->implicit_jacobian,
                                            user->pending_exception,
                                            t,
                                            xdealii,
                                            xdotdealii,
                                            s,
                                            Adealii,
                                            Pdealii);

      petsc_increment_state_counter(P);

      // Handle the Jacobian-free case
      // This call allow to resample the linearization point
      // of the MFFD tangent operator
      PetscBool flg;
      AssertPETSc(PetscObjectTypeCompare((PetscObject)A, MATMFFD, &flg));
      if (flg)
        {
          AssertPETSc(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
          AssertPETSc(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
        }
      else
        petsc_increment_state_counter(A);

      PetscFunctionReturn(err);
    };

    const auto ts_ijacobian_with_setup =
      [](TS, PetscReal t, Vec x, Vec xdot, PetscReal s, Mat A, Mat P, void *ctx)
      -> PetscErrorCode {
      PetscFunctionBeginUser;
      auto user = static_cast<TimeStepper *>(ctx);

      VectorType  xdealii(x);
      VectorType  xdotdealii(xdot);
      AMatrixType Adealii(A);
      PMatrixType Pdealii(P);

      user->A = &Adealii;
      user->P = &Pdealii;
      const int err =
        call_and_possibly_capture_exception(user->setup_jacobian,
                                            user->pending_exception,
                                            t,
                                            xdealii,
                                            xdotdealii,
                                            s);

      petsc_increment_state_counter(P);

      // Handle older versions of PETSc for which we cannot pass a MATSHELL
      // matrix to DMSetMatType. This has been fixed from 3.13 on.
      // What we need to do for older version of PETSc is instead to have
      // a zero matrix with all diagonal entries present.
      if (user->need_dummy_assemble)
        {
          AssertPETSc(MatZeroEntries(P));
          AssertPETSc(MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY));
          AssertPETSc(MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY));
          AssertPETSc(MatShift(P, 0.0));
        }

      // Handle the Jacobian-free case
      // This call allow to resample the linearization point
      // of the MFFD tangent operator
      PetscBool flg;
      AssertPETSc(PetscObjectTypeCompare((PetscObject)A, MATMFFD, &flg));
      if (flg)
        {
          AssertPETSc(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
          AssertPETSc(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
        }
      else
        petsc_increment_state_counter(A);

      PetscFunctionReturn(err);
    };

    const auto ts_rhsfunction =
      [](TS, PetscReal t, Vec x, Vec f, void *ctx) -> PetscErrorCode {
      PetscFunctionBeginUser;
      auto user = static_cast<TimeStepper *>(ctx);

      VectorType xdealii(x);
      VectorType fdealii(f);

      const int err = call_and_possibly_capture_exception(
        user->explicit_function, user->pending_exception, t, xdealii, fdealii);
      petsc_increment_state_counter(f);
      PetscFunctionReturn(err);
    };

    const auto ts_rhsjacobian =
      [](TS, PetscReal t, Vec x, Mat A, Mat P, void *ctx) -> PetscErrorCode {
      PetscFunctionBeginUser;
      auto user = static_cast<TimeStepper *>(ctx);

      VectorType  xdealii(x);
      AMatrixType Adealii(A);
      PMatrixType Pdealii(P);

      const int err =
        call_and_possibly_capture_exception(user->explicit_jacobian,
                                            user->pending_exception,
                                            t,
                                            xdealii,
                                            Adealii,
                                            Pdealii);

      petsc_increment_state_counter(P);

      // Handle older versions of PETSc for which we cannot pass a MATSHELL
      // matrix to DMSetMatType
      if (user->need_dummy_assemble)
        {
          AssertPETSc(MatZeroEntries(P));
          AssertPETSc(MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY));
          AssertPETSc(MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY));
          AssertPETSc(MatShift(P, 0.0));
        }

      // Handle the Jacobian-free case
      // This call allow to resample the linearization point
      // of the MFFD tangent operator
      PetscBool flg;
      AssertPETSc(PetscObjectTypeCompare((PetscObject)A, MATMFFD, &flg));
      if (flg)
        {
          AssertPETSc(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
          AssertPETSc(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
        }
      else
        petsc_increment_state_counter(A);

      PetscFunctionReturn(err);
    };

    const auto ts_monitor =
      [](TS, PetscInt it, PetscReal t, Vec x, void *ctx) -> PetscErrorCode {
      PetscFunctionBeginUser;
      auto user = static_cast<TimeStepper *>(ctx);

      VectorType xdealii(x);
      const int  err = call_and_possibly_capture_exception(
        user->monitor, user->pending_exception, t, xdealii, it);
      PetscFunctionReturn(err);
    };

    AssertThrow(explicit_function || implicit_function,
                StandardExceptions::ExcFunctionNotProvided(
                  "explicit_function || implicit_function"));

    AssertPETSc(TSSetSolution(ts, y.petsc_vector()));

    if (explicit_function)
      AssertPETSc(TSSetRHSFunction(ts, nullptr, ts_rhsfunction, this));

    if (implicit_function)
      AssertPETSc(TSSetIFunction(ts, nullptr, ts_ifunction, this));

    if (setup_jacobian)
      {
        AssertPETSc(TSSetIJacobian(ts,
                                   A ? A->petsc_matrix() : nullptr,
                                   P ? P->petsc_matrix() : nullptr,
                                   ts_ijacobian_with_setup,
                                   this));

        // Tell PETSc to setup a MFFD operator for the linear system matrix
        if (!A)
          set_use_matrix_free(ts, true, false);

        // Do not waste memory by creating a dummy AIJ matrix inside PETSc.
        this->need_dummy_assemble = false;
        if (!P)
          {
#  if DEAL_II_PETSC_VERSION_GTE(3, 13, 0)
            DM   dm;
            SNES snes;
            AssertPETSc(TSGetSNES(ts, &snes));
            AssertPETSc(SNESGetDM(snes, &dm));
            AssertPETSc(DMSetMatType(dm, MATSHELL));
#  else
            this->need_dummy_assemble = true;
#  endif
          }
      }
    else
      {
        if (explicit_jacobian)
          {
            AssertPETSc(TSSetRHSJacobian(ts,
                                         A ? A->petsc_matrix() :
                                             (P ? P->petsc_matrix() : nullptr),
                                         P ? P->petsc_matrix() : nullptr,
                                         ts_rhsjacobian,
                                         this));
          }

        if (implicit_jacobian)
          {
            AssertPETSc(TSSetIJacobian(ts,
                                       A ? A->petsc_matrix() :
                                           (P ? P->petsc_matrix() : nullptr),
                                       P ? P->petsc_matrix() : nullptr,
                                       ts_ijacobian,
                                       this));
          }

        // The user did not set any Jacobian callback. PETSc default in this
        // case is to use FD and thus assemble a dense operator by finite
        // differencing the residual callbacks. Here instead we decide to
        // use a full matrix-free approach by default. This choice can always
        // be overriden from command line.
        if (!explicit_jacobian && !implicit_jacobian)
          {
            set_use_matrix_free(ts, false, true);
          }
      }

    // In case solve_with_jacobian is provided, create a shell
    // preconditioner wrapping the user call. The default internal Krylov
    // solver only applies the preconditioner. This choice
    // can be overriden by command line and users can use any other
    // Krylov method if their solve is not accurate enough.
    // Using solve_with_jacobian as a preconditioner allow users
    // to provide approximate solvers and possibly iterate on a matrix-free
    // approximation of the tangent operator.
    PreconditionShell precond(
      PetscObjectComm(reinterpret_cast<PetscObject>(ts)));
    if (solve_with_jacobian)
      {
        precond.vmult = [&](VectorBase &indst, const VectorBase &insrc) -> int {
          VectorType       dst(static_cast<const Vec &>(indst));
          const VectorType src(static_cast<const Vec &>(insrc));
          return call_and_possibly_capture_exception(solve_with_jacobian,
                                                     pending_exception,
                                                     src,
                                                     dst);
        };

        // Default Krylov solver (preconditioner only)
        SNES snes;
        KSP  ksp;
        AssertPETSc(TSGetSNES(ts, &snes));
        AssertPETSc(SNESGetKSP(snes, &ksp));
        AssertPETSc(KSPSetType(ksp, KSPPREONLY));
        AssertPETSc(KSPSetPC(ksp, precond.get_pc()));
      }

    // Attach user monitoring routine.
    if (monitor)
      AssertPETSc(TSMonitorSet(ts, ts_monitor, this, nullptr));

    // Allow command line customization.
    AssertPETSc(TSSetFromOptions(ts));

    // Handle algebraic components.
    if (algebraic_components && need_dae_tolerances)
      {
        PetscReal atol, rtol;
        AssertPETSc(TSGetTolerances(ts, &atol, nullptr, &rtol, nullptr));

        Vec av, rv;
        AssertPETSc(VecDuplicate(y.petsc_vector(), &av));
        AssertPETSc(VecDuplicate(y.petsc_vector(), &rv));

        VectorBase avdealii(av);
        VectorBase rvdealii(rv);
        avdealii = atol;
        rvdealii = rtol;
        for (auto i : algebraic_components())
          {
            avdealii[i] = -1.0;
            rvdealii[i] = -1.0;
          }
        avdealii.compress(VectorOperation::insert);
        rvdealii.compress(VectorOperation::insert);
        AssertPETSc(TSSetTolerances(ts, atol, av, rtol, rv));
        AssertPETSc(VecDestroy(&av));
        AssertPETSc(VecDestroy(&rv));
      }

    // Having set everything up, now do the actual work
    // and let PETSc do the time stepping. If there is
    // a pending exception, then one of the user callbacks
    // threw an exception we didn't know how to deal with
    // at the time. It is possible that PETSc managed to
    // recover anyway and in that case would have returned
    // a zero error code -- if so, just eat the exception and
    // continue on; otherwise, just rethrow the exception
    // and get outta here.
    const int status = TSSolve(ts, nullptr);
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
    auto nt = ts_get_step_number(ts);

    // Raise an exception if the solver has not converged
    TSConvergedReason reason;
    AssertPETSc(TSGetConvergedReason(ts, &reason));
    AssertThrow(reason > 0,
                ExcMessage("TS solver did not converge after " +
                           std::to_string(nt) + " iterations with reason " +
                           TSConvergedReasons[reason]));

    // Finally return
    return nt;
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
  unsigned int TimeStepper<VectorType, PMatrixType, AMatrixType>::solve(
    VectorType & y,
    PMatrixType &P)
  {
    set_matrix(P);
    return solve(y);
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
  unsigned int TimeStepper<VectorType, PMatrixType, AMatrixType>::solve(
    VectorType & y,
    AMatrixType &A,
    PMatrixType &P)
  {
    set_matrices(A, P);
    return solve(y);
  }

} // namespace PETScWrappers

#  undef AssertPETSc

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC

#endif
