//-----------------------------------------------------------
//
//    Copyright (C) 2022 by the deal.II authors
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
//
// Author: Stefano Zampini, King Abdullah University of Science and Technology.

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
#  define AssertTS(code)                             \
    do                                               \
      {                                              \
        PetscErrorCode ierr = (code);                \
        AssertThrow(ierr == 0, ExcPETScError(ierr)); \
      }                                              \
    while (0)

// Shorthand notation for User error codes.
#  define AssertUser(code, name)                                               \
    do                                                                         \
      {                                                                        \
        int ierr = (code);                                                     \
        AssertThrow(ierr == 0,                                                 \
                    StandardExceptions::ExcFunctionNonzeroReturn(name, ierr)); \
      }                                                                        \
    while (0)

namespace PETScWrappers
{
  template <typename VectorType, typename PMatrixType, typename AMatrixType>
  TimeStepper<VectorType, PMatrixType, AMatrixType>::TimeStepper(
    const TimeStepperData &data,
    const MPI_Comm &       mpi_comm)
  {
    AssertTS(TSCreate(mpi_comm, &ts));
    AssertTS(TSSetApplicationContext(ts, this));
    reinit(data);
  }



  template <typename VectorType, typename PMatrixType, typename AMatrixType>
  TimeStepper<VectorType, PMatrixType, AMatrixType>::~TimeStepper()
  {
    AssertTS(TSDestroy(&ts));
  }



  template <typename VectorType, typename PMatrixType, typename AMatrixType>
  void
  TimeStepper<VectorType, PMatrixType, AMatrixType>::reinit()
  {
    AssertTS(TSReset(ts));
  }



  template <typename VectorType, typename PMatrixType, typename AMatrixType>
  void
  TimeStepper<VectorType, PMatrixType, AMatrixType>::reinit(
    const TimeStepperData &data)
  {
    reinit();

    // Solver type
    if (data.tstype.size())
      AssertTS(TSSetType(ts, data.tstype.c_str()));

    // Options prefix
    if (data.opts_prefix.size())
      AssertTS(TSSetOptionsPrefix(ts, data.opts_prefix.c_str()));

    // Time and steps limits
    AssertTS(TSSetTime(ts, data.initial_time));
    if (data.final_time > data.initial_time)
      AssertTS(TSSetMaxTime(ts, data.final_time));
    if (data.initial_step_size > 0.0)
      AssertTS(TSSetTimeStep(ts, data.initial_step_size));
    if (data.max_steps >= 0)
      AssertTS(TSSetMaxSteps(ts, data.max_steps));

    // Decide how to end the integration. Either stepover the final time or
    // match it.
    AssertTS(TSSetExactFinalTime(ts,
                                 data.match_step ? TS_EXACTFINALTIME_MATCHSTEP :
                                                   TS_EXACTFINALTIME_STEPOVER));

    // Adaptive tolerances
    const PetscReal atol = data.absolute_tolerance > 0.0 ?
                             data.absolute_tolerance :
                             static_cast<PetscReal>(PETSC_DEFAULT);
    const PetscReal rtol = data.relative_tolerance > 0.0 ?
                             data.relative_tolerance :
                             static_cast<PetscReal>(PETSC_DEFAULT);
    AssertTS(TSSetTolerances(ts, atol, nullptr, rtol, nullptr));

    // At this point we do not know the problem size so we cannot
    // set variable tolerances for differential and algebratic equations
    // Store this value and use it during solve.
    this->need_dae_tolerances = data.ignore_algebraic_lte;

    // Adaptive time stepping
    TSAdapt tsadapt;
    AssertTS(TSGetAdapt(ts, &tsadapt));
    AssertTS(TSAdaptSetType(tsadapt, data.tsadapttype.c_str()));

    // As of 3.19, PETSc does not propagate options prefixes to the
    // adaptors.
    if (data.opts_prefix.size())
      AssertTS(TSAdaptSetOptionsPrefix(tsadapt, data.opts_prefix.c_str()));

    // Time step limits
    const PetscReal hmin = data.minimum_step_size > 0.0 ?
                             data.minimum_step_size :
                             static_cast<PetscReal>(PETSC_DEFAULT);
    const PetscReal hmax = data.maximum_step_size > 0.0 ?
                             data.maximum_step_size :
                             static_cast<PetscReal>(PETSC_DEFAULT);
    AssertTS(TSAdaptSetStepLimits(tsadapt, hmin, hmax));
  }



  template <typename VectorType, typename PMatrixType, typename AMatrixType>
  void
  TimeStepper<VectorType, PMatrixType, AMatrixType>::reinit_matrix(
    PMatrixType &P)
  {
    this->A = nullptr;
    this->P = &P;
  }



  template <typename VectorType, typename PMatrixType, typename AMatrixType>
  void
  TimeStepper<VectorType, PMatrixType, AMatrixType>::reinit_matrices(
    AMatrixType &A,
    PMatrixType &P)
  {
    this->A = &A;
    this->P = &P;
  }



  template <typename VectorType, typename PMatrixType, typename AMatrixType>
  unsigned int
  TimeStepper<VectorType, PMatrixType, AMatrixType>::solve(VectorType &y)
  {
    auto ts_ifunction_ =
      [](TS ts, PetscReal t, Vec x, Vec xdot, Vec f, void *ctx)
      -> PetscErrorCode {
      PetscFunctionBeginUser;
      (void)ts;
      auto user = static_cast<TimeStepper *>(ctx);

      VectorType xdealii(x);
      VectorType xdotdealii(xdot);
      VectorType fdealii(f);
      AssertUser(user->implicit_function(t, xdealii, xdotdealii, fdealii),
                 "implicit_function");
      PetscFunctionReturn(0);
    };

    auto ts_ijacobian_ = [](TS        ts,
                            PetscReal t,
                            Vec       x,
                            Vec       xdot,
                            PetscReal s,
                            Mat       A,
                            Mat       P,
                            void *    ctx) -> PetscErrorCode {
      PetscFunctionBeginUser;
      (void)ts;
      auto user = static_cast<TimeStepper *>(ctx);

      VectorType  xdealii(x);
      VectorType  xdotdealii(xdot);
      AMatrixType Adealii(A);
      PMatrixType Pdealii(P);

      AssertUser(
        user->implicit_jacobian(t, xdealii, xdotdealii, s, Adealii, Pdealii),
        "implicit_jacobian");

      // Handle the Jacobian-free case
      // This call allow to resample the linearization point
      // of the MFFD tangent operator
      PetscBool flg;
      AssertTS(PetscObjectTypeCompare((PetscObject)A, MATMFFD, &flg));
      if (flg)
        {
          AssertTS(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
          AssertTS(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
        }
      PetscFunctionReturn(0);
    };

    auto ts_ijacobian_with_setup_ = [](TS        ts,
                                       PetscReal t,
                                       Vec       x,
                                       Vec       xdot,
                                       PetscReal s,
                                       Mat       A,
                                       Mat       P,
                                       void *    ctx) -> PetscErrorCode {
      PetscFunctionBeginUser;
      (void)ts;
      auto user = static_cast<TimeStepper *>(ctx);

      VectorType  xdealii(x);
      VectorType  xdotdealii(xdot);
      AMatrixType Adealii(A);
      PMatrixType Pdealii(P);

      user->A = &Adealii;
      user->P = &Pdealii;
      AssertUser(user->setup_jacobian(t, xdealii, xdotdealii, s),
                 "setup_jacobian");

      // Handle the Jacobian-free case
      // This call allow to resample the linearization point
      // of the MFFD tangent operator
      PetscBool flg;
      AssertTS(PetscObjectTypeCompare((PetscObject)A, MATMFFD, &flg));
      if (flg)
        {
          AssertTS(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
          AssertTS(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
        }
      PetscFunctionReturn(0);
    };

    auto ts_rhsfunction_ =
      [](TS ts, PetscReal t, Vec x, Vec f, void *ctx) -> PetscErrorCode {
      PetscFunctionBeginUser;
      (void)ts;
      auto user = static_cast<TimeStepper *>(ctx);

      VectorType xdealii(x);
      VectorType fdealii(f);

      AssertUser(user->explicit_function(t, xdealii, fdealii),
                 "explicit_function");
      PetscFunctionReturn(0);
    };

    auto ts_rhsjacobian_ =
      [](TS ts, PetscReal t, Vec x, Mat A, Mat P, void *ctx) -> PetscErrorCode {
      PetscFunctionBeginUser;
      (void)ts;
      auto user = static_cast<TimeStepper *>(ctx);

      VectorType  xdealii(x);
      AMatrixType Adealii(A);
      PMatrixType Pdealii(P);

      AssertUser(user->explicit_jacobian(t, xdealii, Adealii, Pdealii),
                 "explicit_jacobian");

      // Handle the Jacobian-free case
      // This call allow to resample the linearization point
      // of the MFFD tangent operator
      PetscBool flg;
      AssertTS(PetscObjectTypeCompare((PetscObject)A, MATMFFD, &flg));
      if (flg)
        {
          AssertTS(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
          AssertTS(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
        }
      PetscFunctionReturn(0);
    };

    auto ts_monitor_ =
      [](TS ts, PetscInt it, PetscReal t, Vec x, void *ctx) -> PetscErrorCode {
      PetscFunctionBeginUser;
      (void)ts;
      auto user = static_cast<TimeStepper *>(ctx);

      VectorType xdealii(x);
      AssertUser(user->monitor(t, xdealii, it), "monitor");
      PetscFunctionReturn(0);
    };

    AssertThrow(explicit_function || implicit_function,
                StandardExceptions::ExcFunctionNotProvided(
                  "explicit_function || implicit_function"));

    AssertTS(TSSetSolution(ts, y.petsc_vector()));

    if (explicit_function)
      AssertTS(TSSetRHSFunction(ts, nullptr, ts_rhsfunction_, this));

    if (implicit_function)
      AssertTS(TSSetIFunction(ts, nullptr, ts_ifunction_, this));

    if (setup_jacobian)
      {
        AssertTS(TSSetIJacobian(ts,
                                A ? A->petsc_matrix() : nullptr,
                                P ? P->petsc_matrix() : nullptr,
                                ts_ijacobian_with_setup_,
                                this));

        // Tell PETSc to setup a MFFD operator for the linear system matrix
        SNES snes;
        AssertTS(TSGetSNES(ts, &snes));
        if (!A)
          AssertTS(SNESSetUseMatrixFree(snes, PETSC_TRUE, PETSC_FALSE));

        // Do not waste memory by creating a dummy AIJ matrix inside PETSc.
        if (!P)
          {
            DM dm;
            AssertTS(SNESGetDM(snes, &dm));
            AssertTS(DMSetMatType(dm, MATSHELL));
          }
      }
    else
      {
        if (explicit_jacobian)
          {
            AssertTS(TSSetRHSJacobian(ts,
                                      A ? A->petsc_matrix() :
                                          (P ? P->petsc_matrix() : nullptr),
                                      P ? P->petsc_matrix() : nullptr,
                                      ts_rhsjacobian_,
                                      this));
          }

        if (implicit_jacobian)
          {
            AssertTS(TSSetIJacobian(ts,
                                    A ? A->petsc_matrix() :
                                        (P ? P->petsc_matrix() : nullptr),
                                    P ? P->petsc_matrix() : nullptr,
                                    ts_ijacobian_,
                                    this));
          }

        // The user did not set any Jacobian callback. PETSc default in this
        // case is to use FD and thus assemble a dense operator by finite
        // differencing the residual callbacks. Here instead we decide to
        // use a full matrix-free approach by default. This choice can always
        // be overriden from command line.
        if (!explicit_jacobian && !implicit_jacobian)
          {
            SNES snes;
            AssertTS(TSGetSNES(ts, &snes));
            AssertTS(SNESSetUseMatrixFree(snes, PETSC_FALSE, PETSC_TRUE));
          }
      }

    // In case solve_for_jacobian_system is provided, create a shell
    // preconditioner wrapping the user call. The default internal Krylov
    // solver only applies the preconditioner. This choice
    // can be overriden by command line and users can use any other
    // Krylov method if their solve is not accurate enough.
    // Using solve_for_jacobian_system as a preconditioner allow users
    // to provide approximate solvers and possibly iterate on a matrix-free
    // approximation of the tangent operator.
    PreconditionShell precond(
      PetscObjectComm(reinterpret_cast<PetscObject>(ts)));
    if (solve_for_jacobian_system)
      {
        precond.vmult = [&](VectorBase &indst, const VectorBase &insrc) -> int {
          VectorType       dst(static_cast<const Vec &>(indst));
          const VectorType src(static_cast<const Vec &>(insrc));
          return solve_for_jacobian_system(src, dst);
        };

        // Default Krylov solver (preconditioner only)
        SNES snes;
        KSP  ksp;
        AssertTS(TSGetSNES(ts, &snes));
        AssertTS(SNESGetKSP(snes, &ksp));
        AssertTS(KSPSetType(ksp, KSPPREONLY));
        AssertTS(KSPSetPC(ksp, precond.get_pc()));
      }

    // Attach user monitoring routine.
    if (monitor)
      AssertTS(TSMonitorSet(ts, ts_monitor_, this, nullptr));

    // Allow command line customization.
    AssertTS(TSSetFromOptions(ts));

    // Handle algebraic components.
    if (algebraic_components && need_dae_tolerances)
      {
        PetscReal atol, rtol;
        AssertTS(TSGetTolerances(ts, &atol, nullptr, &rtol, nullptr));

        Vec av, rv;
        AssertTS(VecDuplicate(y.petsc_vector(), &av));
        AssertTS(VecDuplicate(y.petsc_vector(), &rv));

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
        AssertTS(TSSetTolerances(ts, atol, av, rtol, rv));
        AssertTS(VecDestroy(&av));
        AssertTS(VecDestroy(&rv));
      }

    // Having set everything up, now do the actual work
    // and let PETSc do the time stepping.
    AssertTS(TSSolve(ts, nullptr));

    // Return the number of steps taken.
    PetscInt nt;
    AssertTS(TSGetStepNumber(ts, &nt));
    return nt;
  }



  template <typename VectorType, typename PMatrixType, typename AMatrixType>
  unsigned int
  TimeStepper<VectorType, PMatrixType, AMatrixType>::solve(VectorType & y,
                                                           PMatrixType &P)
  {
    reinit_matrix(P);
    return solve(y);
  }



  template <typename VectorType, typename PMatrixType, typename AMatrixType>
  unsigned int
  TimeStepper<VectorType, PMatrixType, AMatrixType>::solve(VectorType & y,
                                                           AMatrixType &A,
                                                           PMatrixType &P)
  {
    reinit_matrices(A, P);
    return solve(y);
  }

} // namespace PETScWrappers

#  undef AssertTS
#  undef AssertUser

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC

#endif
