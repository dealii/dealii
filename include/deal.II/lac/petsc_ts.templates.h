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

#ifndef dealii_petsc_ts_templates_h
#define dealii_petsc_ts_templates_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PETSC
#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/mpi_stub.h>

#  include <deal.II/lac/petsc_precondition.h>
#  include <deal.II/lac/petsc_ts.h>

#  include <petscdm.h>
#  include <petscerror.h>
#  include <petscts.h>

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
  call_and_possibly_capture_ts_exception(
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
  TimeStepper<VectorType, PMatrixType, AMatrixType>::TimeStepper(
    const TimeStepperData &data,
    const MPI_Comm         mpi_comm)
    : solve_with_jacobian_pc(mpi_comm)
    , pending_exception(nullptr)
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
    AssertPETSc(PetscObjectComposeFunction((PetscObject)ts,
                                           "__dealii_ts_resize_setup__",
                                           nullptr));
    AssertPETSc(PetscObjectComposeFunction((PetscObject)ts,
                                           "__dealii_ts_resize_transfer__",
                                           nullptr));
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
    PetscReal t;
    AssertPETSc(TSGetTime(ts, &t));
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
  unsigned int TimeStepper<VectorType, PMatrixType, AMatrixType>::
    get_step_number()
  {
    return ts_get_step_number(ts);
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
    PetscReal dt;

    AssertPETSc(TSGetTimeStep(ts, &dt));
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

    // By default we always allow recovery from failed nonlinear solves
    AssertPETSc(TSSetMaxSNESFailures(ts, -1));

    // By default we do not want PETSc to return a nonzero error code from
    // TSSolve
    AssertPETSc(TSSetErrorIfStepFails(ts, PETSC_FALSE));
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

    // Store the boolean to be passed to TSSetResize.
    this->restart_if_remesh = data.restart_if_remesh;

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
  void TimeStepper<VectorType, PMatrixType, AMatrixType>::setup_callbacks()
  {
    const auto ts_resize_setup =
      [](TS ts, PetscInt it, PetscReal t, Vec x, PetscBool *flg, void *ctx)
      -> PetscErrorCode {
      PetscFunctionBeginUser;
      auto user = static_cast<TimeStepper *>(ctx);

      VectorType xdealii(x);
      bool       resize = false;

      const int lineno = __LINE__;
      const int err    = call_and_possibly_capture_ts_exception(
        (user->decide_and_prepare_for_remeshing ?
              [&](const real_type    t,
               const unsigned int step,
               const VectorType  &y,
               bool              &resize) {
             resize = user->decide_and_prepare_for_remeshing(t, step, y);
           } :
              user->decide_for_coarsening_and_refinement),
        user->pending_exception,
        {},
        t,
        it,
        xdealii,
        resize);
      *flg = resize ? PETSC_TRUE : PETSC_FALSE;
      if (err)
        return PetscError(
          PetscObjectComm((PetscObject)ts),
          lineno + 1,
          "decide_and_prepare_for_remeshing",
          __FILE__,
          PETSC_ERR_LIB,
          PETSC_ERROR_INITIAL,
          "Failure in ts_resize_setup from dealii::PETScWrappers::TimeStepper");

      PetscFunctionReturn(PETSC_SUCCESS);
    };

    const auto ts_resize_transfer = [](TS       ts,
                                       PetscInt nv,
                                       Vec      vin[],
                                       Vec      vout[],
                                       void    *ctx) -> PetscErrorCode {
      PetscFunctionBeginUser;
      auto user = static_cast<TimeStepper *>(ctx);

      std::vector<VectorType> all_in, all_out;
      for (PetscInt i = 0; i < nv; i++)
        {
          all_in.push_back(VectorType(vin[i]));
        }

      const int lineno = __LINE__;
      const int err    = call_and_possibly_capture_ts_exception(
        (user->transfer_solution_vectors_to_new_mesh ?
              // If we can, call the new callback
           user->transfer_solution_vectors_to_new_mesh :
              // otherwise, use the old one where we just ignore the time:
           [user](const real_type /*t*/,
                  const std::vector<VectorType> &all_in,
                  std::vector<VectorType>       &all_out) {
             user->interpolate(all_in, all_out);
           }),
        user->pending_exception,
        {},
        user->get_time(),
        all_in,
        all_out);
      if (err)
        return PetscError(
          PetscObjectComm((PetscObject)ts),
          lineno + 1,
          "interpolate",
          __FILE__,
          PETSC_ERR_LIB,
          PETSC_ERROR_INITIAL,
          "Failure in ts_resize_transfer from dealii::PETScWrappers::TimeStepper");

      if (all_out.size() != all_in.size())
        SETERRQ(PetscObjectComm((PetscObject)ts),
                PETSC_ERR_LIB,
                "Broken all_out");

      for (PetscInt i = 0; i < nv; i++)
        {
          vout[i] = all_out[i].petsc_vector();
          PetscCall(PetscObjectReference((PetscObject)vout[i]));
        }

      try
        {
          user->setup_callbacks();
        }
      catch (...)
        {
          SETERRQ(PetscObjectComm((PetscObject)ts),
                  PETSC_ERR_LIB,
                  "Exception in setup_callbacks");
        }

      try
        {
          user->setup_algebraic_constraints(all_out[0]);
        }
      catch (...)
        {
          SETERRQ(PetscObjectComm((PetscObject)ts),
                  PETSC_ERR_LIB,
                  "Exception in setup_algebraic_constraints");
        }
      PetscFunctionReturn(PETSC_SUCCESS);
    };

    const auto ts_poststep_amr = [](TS ts) -> PetscErrorCode {
      PetscFunctionBeginUser;
      void *ctx;
      PetscCall(TSGetApplicationContext(ts, &ctx));
      auto user = static_cast<TimeStepper *>(ctx);

      using transfer_setup_fn =
        PetscErrorCode (*)(TS, PetscInt, PetscReal, Vec, PetscBool *, void *);
      using transfer_fn =
        PetscErrorCode (*)(TS, PetscInt, Vec[], Vec[], void *);
      transfer_setup_fn transfer_setup;
      transfer_fn       transfer;

      AssertPETSc(PetscObjectQueryFunction((PetscObject)ts,
                                           "__dealii_ts_resize_setup__",
                                           &transfer_setup));
      AssertPETSc(PetscObjectQueryFunction((PetscObject)ts,
                                           "__dealii_ts_resize_transfer__",
                                           &transfer));

      PetscBool resize = PETSC_FALSE;
      Vec       x;
      PetscCall(TSGetSolution(ts, &x));
      PetscCall(PetscObjectReference((PetscObject)x));
      PetscCall(transfer_setup(
        ts, user->get_step_number(), user->get_time(), x, &resize, user));

      if (resize)
        {
          Vec new_x;

          user->reinit();

          // Reinitialize DM. Currently it is not possible to do with public API
          ts_reset_dm(ts);

          PetscCall(transfer(ts, 1, &x, &new_x, user));
          PetscCall(TSSetSolution(ts, new_x));
          PetscCall(VecDestroy(&new_x));
        }
#  if DEAL_II_PETSC_VERSION_LT(3, 17, 0)
      // Older versions of PETSc assume that the user does not
      // change the solution vector during TSPostStep
      // We "fix" it by taking a reference to the object and
      // increment its state so that the time stepper is restarted
      // properly.
      AssertPETSc(PetscObjectCompose((PetscObject)ts,
                                     "__dealii_ts_resize_bug__",
                                     (PetscObject)x));
      petsc_increment_state_counter(x);
#  endif
      PetscCall(VecDestroy(&x));
      PetscFunctionReturn(PETSC_SUCCESS);
    };

    const auto ts_prestage = [](TS ts, PetscReal) -> PetscErrorCode {
      PetscFunctionBeginUser;
      void *ctx;
      PetscCall(TSGetApplicationContext(ts, &ctx));
      auto user               = static_cast<TimeStepper *>(ctx);
      user->error_in_function = false;

      // only in case we are using an implicit solver
      if (ts_has_snes(ts))
        {
          SNES snes;
          PetscCall(TSGetSNES(ts, &snes));
          snes_reset_domain_flags(snes);
        }
      PetscFunctionReturn(PETSC_SUCCESS);
    };

    const auto ts_poststage =
      [](TS ts, PetscReal t, PetscInt i, Vec *xx) -> PetscErrorCode {
      PetscFunctionBeginUser;
      void *ctx;
      PetscCall(TSGetApplicationContext(ts, &ctx));
      auto user = static_cast<TimeStepper *>(ctx);

      VectorType xdealii(xx[i]);

      const int lineno = __LINE__;
      const int err    = call_and_possibly_capture_ts_exception(
        // Call the user-provided callback. Use the new name if used,
        // or the legacy name otherwise.
        (user->update_constrained_components ?
              user->update_constrained_components :
              user->distribute),
        user->pending_exception,
        [user, ts]() -> void {
          user->error_in_function = true;

          SNES snes;
          AssertPETSc(TSGetSNES(ts, &snes));
          AssertPETSc(SNESSetFunctionDomainError(snes));
        },
        t,
        xdealii);
      if (err)
        return PetscError(
          PetscObjectComm((PetscObject)ts),
          lineno + 1,
          "update_constrained_components",
          __FILE__,
          PETSC_ERR_LIB,
          PETSC_ERROR_INITIAL,
          "Failure in ts_poststage from dealii::PETScWrappers::TimeStepper");
      petsc_increment_state_counter(xx[i]);
      PetscFunctionReturn(PETSC_SUCCESS);
    };

    const auto ts_functiondomainerror =
      [](TS ts, PetscReal, Vec, PetscBool *accept) -> PetscErrorCode {
      PetscFunctionBeginUser;
      void *ctx;
      AssertPETSc(TSGetApplicationContext(ts, &ctx));
      auto user = static_cast<TimeStepper *>(ctx);
      *accept   = user->error_in_function ? PETSC_FALSE : PETSC_TRUE;
      PetscFunctionReturn(PETSC_SUCCESS);
    };

    const auto ts_ifunction =
      [](TS ts, PetscReal t, Vec x, Vec xdot, Vec f, void *ctx)
      -> PetscErrorCode {
      PetscFunctionBeginUser;
      auto user = static_cast<TimeStepper *>(ctx);

      VectorType xdealii(x);
      VectorType xdotdealii(xdot);
      VectorType fdealii(f);

      const int lineno = __LINE__;
      const int err    = call_and_possibly_capture_ts_exception(
        user->implicit_function,
        user->pending_exception,
        [user, ts]() -> void {
          user->error_in_function = true;

          SNES snes;
          AssertPETSc(TSGetSNES(ts, &snes));
          AssertPETSc(SNESSetFunctionDomainError(snes));
        },
        t,
        xdealii,
        xdotdealii,
        fdealii);
      if (err)
        return PetscError(
          PetscObjectComm((PetscObject)ts),
          lineno + 1,
          "implicit_function",
          __FILE__,
          PETSC_ERR_LIB,
          PETSC_ERROR_INITIAL,
          "Failure in ts_ifunction from dealii::PETScWrappers::TimeStepper");
      petsc_increment_state_counter(f);
      PetscFunctionReturn(PETSC_SUCCESS);
    };

    const auto ts_ijacobian = [](TS        ts,
                                 PetscReal t,
                                 Vec       x,
                                 Vec       xdot,
                                 PetscReal s,
                                 Mat       A,
                                 Mat       P,
                                 void     *ctx) -> PetscErrorCode {
      PetscFunctionBeginUser;
      auto user = static_cast<TimeStepper *>(ctx);

      VectorType  xdealii(x);
      VectorType  xdotdealii(xdot);
      AMatrixType Adealii(A);
      PMatrixType Pdealii(P);

      const int lineno = __LINE__;
      const int err    = call_and_possibly_capture_ts_exception(
        user->implicit_jacobian,
        user->pending_exception,
        [ts]() -> void {
          SNES snes;
          AssertPETSc(TSGetSNES(ts, &snes));
          snes_set_jacobian_domain_error(snes);
        },
        t,
        xdealii,
        xdotdealii,
        s,
        Adealii,
        Pdealii);
      if (err)
        return PetscError(
          PetscObjectComm((PetscObject)ts),
          lineno + 1,
          "implicit_jacobian",
          __FILE__,
          PETSC_ERR_LIB,
          PETSC_ERROR_INITIAL,
          "Failure in ts_ijacobian from dealii::PETScWrappers::TimeStepper");
      petsc_increment_state_counter(P);

      // Handle the Jacobian-free case
      // This call allow to resample the linearization point
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

    const auto ts_ijacobian_with_setup = [](TS        ts,
                                            PetscReal t,
                                            Vec       x,
                                            Vec       xdot,
                                            PetscReal s,
                                            Mat       A,
                                            Mat       P,
                                            void     *ctx) -> PetscErrorCode {
      PetscFunctionBeginUser;
      auto user = static_cast<TimeStepper *>(ctx);

      VectorType xdealii(x);
      VectorType xdotdealii(xdot);

      const int lineno = __LINE__;
      const int err    = call_and_possibly_capture_ts_exception(
        user->setup_jacobian,
        user->pending_exception,
        [ts]() -> void {
          SNES snes;
          AssertPETSc(TSGetSNES(ts, &snes));
          snes_set_jacobian_domain_error(snes);
        },
        t,
        xdealii,
        xdotdealii,
        s);
      if (err)
        return PetscError(
          PetscObjectComm((PetscObject)ts),
          lineno + 1,
          "setup_jacobian",
          __FILE__,
          PETSC_ERR_LIB,
          PETSC_ERROR_INITIAL,
          "Failure in ts_ijacobian_with_setup from dealii::PETScWrappers::TimeStepper");
      // The MatCopy calls below are 99% of the times dummy calls.
      // They are only used in case we get different Mats then the one we passed
      // to TSSetIJacobian.
      if (user->P)
        PetscCall(MatCopy(user->P->petsc_matrix(), P, SAME_NONZERO_PATTERN));
      if (user->A)
        PetscCall(MatCopy(user->A->petsc_matrix(), A, SAME_NONZERO_PATTERN));
      petsc_increment_state_counter(P);

      // Handle older versions of PETSc for which we cannot pass a MATSHELL
      // matrix to DMSetMatType. This has been fixed from 3.13 on.
      // What we need to do for older version of PETSc is instead to have
      // a zero matrix with all diagonal entries present.
      // We use PETSC_MACHINE_EPSILON because some versions of PETSc skip
      // MatShift if we pass 0.0
      if (user->need_dummy_assemble)
        {
          PetscCall(MatZeroEntries(P));
          PetscCall(MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY));
          PetscCall(MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY));
          PetscCall(MatShift(P, PETSC_MACHINE_EPSILON));
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

    const auto ts_rhsfunction =
      [](TS ts, PetscReal t, Vec x, Vec f, void *ctx) -> PetscErrorCode {
      PetscFunctionBeginUser;
      auto user = static_cast<TimeStepper *>(ctx);

      VectorType xdealii(x);
      VectorType fdealii(f);

      const int lineno = __LINE__;
      const int err    = call_and_possibly_capture_ts_exception(
        user->explicit_function,
        user->pending_exception,
        [user, ts]() -> void {
          user->error_in_function = true;

          if (ts_has_snes(ts))
            {
              SNES snes;
              AssertPETSc(TSGetSNES(ts, &snes));
              AssertPETSc(SNESSetFunctionDomainError(snes));
            }
        },
        t,
        xdealii,
        fdealii);
      if (err)
        return PetscError(
          PetscObjectComm((PetscObject)ts),
          lineno + 1,
          "explicit_function",
          __FILE__,
          PETSC_ERR_LIB,
          PETSC_ERROR_INITIAL,
          "Failure in ts_rhsfunction from dealii::PETScWrappers::TimeStepper");
      petsc_increment_state_counter(f);
      PetscFunctionReturn(PETSC_SUCCESS);
    };

    const auto ts_rhsjacobian =
      [](TS ts, PetscReal t, Vec x, Mat A, Mat P, void *ctx) -> PetscErrorCode {
      PetscFunctionBeginUser;
      auto user = static_cast<TimeStepper *>(ctx);

      VectorType  xdealii(x);
      AMatrixType Adealii(A);
      PMatrixType Pdealii(P);

      const int lineno = __LINE__;
      const int err    = call_and_possibly_capture_ts_exception(
        user->explicit_jacobian,
        user->pending_exception,
        [ts]() -> void {
          SNES snes;
          AssertPETSc(TSGetSNES(ts, &snes));
          snes_set_jacobian_domain_error(snes);
        },
        t,
        xdealii,
        Adealii,
        Pdealii);
      if (err)
        return PetscError(
          PetscObjectComm((PetscObject)ts),
          lineno + 1,
          "explicit_jacobian",
          __FILE__,
          PETSC_ERR_LIB,
          PETSC_ERROR_INITIAL,
          "Failure in ts_rhsjacobian from dealii::PETScWrappers::TimeStepper");
      petsc_increment_state_counter(P);

      // Handle the Jacobian-free case
      // This call allow to resample the linearization point
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

    const auto ts_monitor =
      [](TS ts, PetscInt it, PetscReal t, Vec x, void *ctx) -> PetscErrorCode {
      PetscFunctionBeginUser;
      auto user = static_cast<TimeStepper *>(ctx);

      VectorType xdealii(x);

      const int lineno = __LINE__;
      const int err    = call_and_possibly_capture_ts_exception(
        user->monitor, user->pending_exception, {}, t, xdealii, it);
      if (err)
        return PetscError(
          PetscObjectComm((PetscObject)ts),
          lineno + 1,
          "monitor",
          __FILE__,
          PETSC_ERR_LIB,
          PETSC_ERROR_INITIAL,
          "Failure in ts_monitor from dealii::PETScWrappers::TimeStepper");
      PetscFunctionReturn(PETSC_SUCCESS);
    };

    AssertThrow(explicit_function || implicit_function,
                StandardExceptions::ExcFunctionNotProvided(
                  "explicit_function || implicit_function"));

    // Handle recoverable errors.
    this->error_in_function = false;
    AssertPETSc(TSSetPreStage(ts, ts_prestage));
    AssertPETSc(TSSetFunctionDomainError(ts, ts_functiondomainerror));

    // Problem description.
    if (explicit_function)
      AssertPETSc(TSSetRHSFunction(ts, nullptr, ts_rhsfunction, this));
    if (implicit_function)
      AssertPETSc(TSSetIFunction(ts, nullptr, ts_ifunction, this));

    // Handle Jacobians.
    this->need_dummy_assemble = false;
    if (setup_jacobian)
      {
        AssertPETSc(TSSetIJacobian(ts,
                                   A ? A->petsc_matrix() : nullptr,
                                   P ? P->petsc_matrix() : nullptr,
                                   ts_ijacobian_with_setup,
                                   this));

        // Tell PETSc to set up a MFFD operator for the linear system matrix
        if (!A)
          set_use_matrix_free(ts, true, false);

        // Do not waste memory by creating a dummy AIJ matrix inside PETSc.
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
        // be overridden from command line.
        if (!explicit_jacobian && !implicit_jacobian)
          {
            set_use_matrix_free(ts, false, true);
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
    if (solve_with_jacobian)
      {
        solve_with_jacobian_pc.vmult = [&](VectorBase       &indst,
                                           const VectorBase &insrc) -> void {
          VectorType       dst(static_cast<const Vec &>(indst));
          const VectorType src(static_cast<const Vec &>(insrc));
          solve_with_jacobian(src, dst);
        };

        // Default Krylov solver (preconditioner only)
        SNES snes;
        KSP  ksp;
        AssertPETSc(TSGetSNES(ts, &snes));
        AssertPETSc(SNESGetKSP(snes, &ksp));
        AssertPETSc(KSPSetType(ksp, KSPPREONLY));
        AssertPETSc(KSPSetPC(ksp, solve_with_jacobian_pc.get_pc()));
      }

    if (distribute || update_constrained_components)
      {
        if (update_constrained_components)
          Assert(!distribute,
                 ExcMessage(
                   "The 'distribute' callback name of the TimeStepper "
                   "class is deprecated. If you are setting the equivalent "
                   "'update_constrained_components' callback, you cannot also "
                   "set the 'distribute' callback."));
        AssertPETSc(TSSetPostStage(ts, ts_poststage));
      }

    // Attach user monitoring routine.
    if (monitor)
      AssertPETSc(TSMonitorSet(ts, ts_monitor, this, nullptr));

    // Handle AMR.
    if (decide_and_prepare_for_remeshing ||
        decide_for_coarsening_and_refinement)
      {
        if (decide_and_prepare_for_remeshing)
          Assert(
            !decide_for_coarsening_and_refinement,
            ExcMessage(
              "The 'decide_for_coarsening_and_refinement' callback name "
              "of the TimeStepper class is deprecated. If you are setting "
              "the equivalent 'decide_and_prepare_for_remeshing' callback, you "
              "cannot also set the 'decide_for_coarsening_and_refinement' "
              "callback."));

        // If we have the decide_and_prepare_for_remeshing callback
        // set, then we also need to have the callback for actually
        // transferring the solution:
        AssertThrow(interpolate || transfer_solution_vectors_to_new_mesh,
                    StandardExceptions::ExcFunctionNotProvided(
                      "transfer_solution_vectors_to_new_mesh"));

        if (transfer_solution_vectors_to_new_mesh)
          Assert(
            !interpolate,
            ExcMessage(
              "The 'interpolate' callback name "
              "of the TimeStepper class is deprecated. If you are setting "
              "the equivalent 'transfer_solution_vectors_to_new_mesh' callback, you "
              "cannot also set the 'interpolate' "
              "callback."));

#  if DEAL_II_PETSC_VERSION_GTE(3, 21, 0)
        (void)ts_poststep_amr;
        AssertPETSc(TSSetResize(ts,
                                static_cast<PetscBool>(restart_if_remesh),
                                ts_resize_setup,
                                ts_resize_transfer,
                                this));
#  else
        AssertThrow(!restart_if_remesh,
                    ExcMessage(
                      "Restart with remesh supported from PETSc 3.21."));
#    if DEAL_II_PETSC_VERSION_GTE(3, 20, 0)
        (void)ts_poststep_amr;
        AssertPETSc(TSSetResize(ts, ts_resize_setup, ts_resize_transfer, this));
#    else
        AssertPETSc(PetscObjectComposeFunction(
          (PetscObject)ts,
          "__dealii_ts_resize_setup__",
          static_cast<PetscErrorCode (*)(
            TS, PetscInt, PetscReal, Vec, PetscBool *, void *)>(
            ts_resize_setup)));
        AssertPETSc(PetscObjectComposeFunction(
          (PetscObject)ts,
          "__dealii_ts_resize_transfer__",
          static_cast<PetscErrorCode (*)(TS, PetscInt, Vec[], Vec[], void *)>(
            ts_resize_transfer)));
        AssertPETSc(TSSetPostStep(ts, ts_poststep_amr));
#    endif
#  endif
      }

    // Allow command line customization.
    AssertPETSc(TSSetFromOptions(ts));

    // By default PETSc does not check for Jacobian errors in optimized
    // mode. Here we do it unconditionally.
#  if DEAL_II_PETSC_VERSION_GTE(3, 11, 0)
    {
      SNES snes;
      AssertPETSc(TSGetSNES(ts, &snes));
      AssertPETSc(SNESSetCheckJacobianDomainError(snes, PETSC_TRUE));
    }
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
  void TimeStepper<VectorType, PMatrixType, AMatrixType>::
    setup_algebraic_constraints(const VectorType &y)
  {
    if (algebraic_components && need_dae_tolerances)
      {
        PetscReal atol, rtol;
        AssertPETSc(TSGetTolerances(ts, &atol, nullptr, &rtol, nullptr));

        // Fill up vectors with individual tolerances, using -1 for the
        // algebraic variables.
        // We need to input them with the same vector type of our solution
        // but we are not sure that the user's VectorType supports operator[]
        // We thus first create standard vectors that we know support the
        // operator, and then copy the values into the user's type operator.
        Vec av, rv;
        Vec py = const_cast<VectorType &>(y).petsc_vector();
        AssertPETSc(VecDuplicate(py, &av));
        AssertPETSc(VecDuplicate(py, &rv));

        // Standard vectors
        Vec      avT, rvT;
        PetscInt n, N;
        AssertPETSc(VecGetLocalSize(py, &n));
        AssertPETSc(VecGetSize(py, &N));
        AssertPETSc(VecCreateMPI(this->get_mpi_communicator(), n, N, &avT));
        AssertPETSc(VecDuplicate(avT, &rvT));

        // Fill-up the vectors
        VectorBase avdealii(avT);
        VectorBase rvdealii(rvT);
        avdealii = atol;
        rvdealii = rtol;
        for (auto i : algebraic_components())
          {
            avdealii[i] = -1.0;
            rvdealii[i] = -1.0;
          }
        avdealii.compress(VectorOperation::insert);
        rvdealii.compress(VectorOperation::insert);

        // Copy, set and destroy
        AssertPETSc(VecCopy(avT, av));
        AssertPETSc(VecCopy(rvT, rv));
        AssertPETSc(TSSetTolerances(ts, atol, av, rtol, rv));
        AssertPETSc(VecDestroy(&av));
        AssertPETSc(VecDestroy(&rv));
        AssertPETSc(VecDestroy(&avT));
        AssertPETSc(VecDestroy(&rvT));
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
  unsigned int TimeStepper<VectorType, PMatrixType, AMatrixType>::solve(
    VectorType &y)
  {
    // Set up PETSc data structure.
    setup_callbacks();

    // Set solution vector.
    AssertPETSc(TSSetSolution(ts, y.petsc_vector()));

    // Setup internal data structures, including the internal nonlinear
    // solver SNES if using an implicit solver.
    AssertPETSc(TSSetUp(ts));

    // Handle algebraic components. This must be done when we know the layout
    // of the solution vector.
    setup_algebraic_constraints(y);

    // Having set everything up, now do the actual work
    // and let PETSc do the time stepping. If there is
    // a pending exception, then one of the user callbacks
    // threw an exception we didn't know how to deal with
    // at the time. It is possible that PETSc managed to
    // recover anyway and in that case would have returned
    // a zero error code -- if so, just eat the exception and
    // continue on; otherwise, just rethrow the exception
    // and get outta here.
    const PetscErrorCode status = TSSolve(ts, nullptr);
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
    VectorType  &y,
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
    VectorType  &y,
    AMatrixType &A,
    PMatrixType &P)
  {
    set_matrices(A, P);
    return solve(y);
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
