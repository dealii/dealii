// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_sundials_arkode_stepper_templates_h
#define dealii_sundials_arkode_stepper_templates_h

#include <deal.II/base/config.h>

#include <deal.II/sundials/arkode_stepper.h>

#ifdef DEAL_II_WITH_SUNDIALS

#  include <deal.II/base/utilities.h>

#  include <deal.II/sundials/arkode_exception.h>
#  include <deal.II/sundials/invocation_context.h>
#  include <deal.II/sundials/n_vector.h>
#  include <deal.II/sundials/sunlinsol_wrapper.h>
#  include <deal.II/sundials/utilities.h>

#  include <arkode/arkode_arkstep.h>
#  include <sunlinsol/sunlinsol_spgmr.h>
#  include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

#  include <iostream>

DEAL_II_NAMESPACE_OPEN

namespace SUNDIALS
{
  template <typename VectorType>
  ARKStepper<VectorType>::ARKStepper(const AdditionalData &data)
    : arkode_mem(nullptr)
    , data(data)
  {}

  template <typename VectorType>
  ARKStepper<VectorType>::~ARKStepper()
  {
    Assert(arkode_mem != nullptr, ExcInternalError());

    ARKStepFree(&arkode_mem);
  }

  template <typename VectorType>
  void
  ARKStepper<VectorType>::reinit(double                      t0,
                                 const VectorType           &y0,
                                 internal::InvocationContext inv_ctx)
  {
    if (arkode_mem)
      {
        ARKStepFree(&arkode_mem);
        // Initialization is version-dependent: do that in a moment
      }

    Assert(explicit_function || implicit_function,
           ExcFunctionNotProvided("explicit_function || implicit_function"));

    auto explicit_function_callback = [](SUNDIALS::realtype tt,
                                         N_Vector           yy,
                                         N_Vector           yp,
                                         void              *user_data) -> int {
      Assert(user_data != nullptr, ExcInternalError());
      auto &callback_ctx =
        *static_cast<ARKCallbackContext<ARKStepper<VectorType>> *>(user_data);

      auto *src_yy = internal::unwrap_nvector_const<VectorType>(yy);
      auto *dst_yp = internal::unwrap_nvector<VectorType>(yp);

      return Utilities::call_and_possibly_capture_exception(
        callback_ctx.stepper->explicit_function,
        *callback_ctx.pending_exception,
        tt,
        *src_yy,
        *dst_yp);
    };

    auto implicit_function_callback = [](SUNDIALS::realtype tt,
                                         N_Vector           yy,
                                         N_Vector           yp,
                                         void              *user_data) -> int {
      Assert(user_data != nullptr, ExcInternalError());
      auto &callback_ctx =
        *static_cast<ARKCallbackContext<ARKStepper<VectorType>> *>(user_data);

      auto *src_yy = internal::unwrap_nvector_const<VectorType>(yy);
      auto *dst_yp = internal::unwrap_nvector<VectorType>(yp);

      return Utilities::call_and_possibly_capture_exception(
        callback_ctx.stepper->implicit_function,
        *callback_ctx.pending_exception,
        tt,
        *src_yy,
        *dst_yp);
    };

    // LSRKStepCreateSSP

    // just a view on the memory in solution, all write operations on yy by
    // ARKODE will automatically be mirrored to solution
    auto initial_condition_nvector =
      internal::make_nvector_view(y0
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                                  ,
                                  inv_ctx.arkode_ctx
#  endif
      );

    arkode_mem = ARKStepCreate(explicit_function ? explicit_function_callback :
                                                   ARKRhsFn(nullptr),
                               implicit_function ? implicit_function_callback :
                                                   ARKRhsFn(nullptr),
                               t0,
                               initial_condition_nvector
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                               ,
                               inv_ctx.arkode_ctx
#  endif
    );
    Assert(arkode_mem != nullptr, ExcInternalError());

    // for version 4.0.0 this call must be made before solver settings
    callback_ctx = {this, &inv_ctx.pending_exception};
    int status   = ARKStepSetUserData(arkode_mem, &callback_ctx);
    AssertARKode(status);

    setup_system_solver(y0, inv_ctx);

    setup_mass_solver(y0, inv_ctx);

    if (custom_setup)
      custom_setup(arkode_mem);
  }

  template <typename VectorType>
  void
  ARKStepper<VectorType>::setup_system_solver(
    const VectorType           &solution,
    internal::InvocationContext inv_ctx)
  {
    // no point in setting up a solver when there is no linear system
    if (!implicit_function)
      return;

    int status;
    (void)status;
    // Initialize solver
    // currently only iterative linear solvers are supported
    if (jacobian_times_vector)
      {
        SUNLinearSolver sun_linear_solver;
        if (solve_linearized_system)
          {
            linear_solver =
              std::make_unique<internal::LinearSolverWrapper<VectorType>>(
                solve_linearized_system,
                inv_ctx.pending_exception
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                ,
                inv_ctx.arkode_ctx
#  endif
              );
            sun_linear_solver = *linear_solver;
          }
        else
          {
            // use default solver from SUNDIALS
            // TODO give user options
            auto y_template = internal::make_nvector_view(solution
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                                                          ,
                                                          inv_ctx.arkode_ctx
#  endif
            );
#  if DEAL_II_SUNDIALS_VERSION_LT(6, 0, 0)
            sun_linear_solver =
              SUNLinSol_SPGMR(y_template,
                              PREC_NONE,
                              0 /*krylov subvectors, 0 uses default*/);
#  else
            sun_linear_solver =
              SUNLinSol_SPGMR(y_template,
                              SUN_PREC_NONE,
                              0 /*krylov subvectors, 0 uses default*/,
                              inv_ctx.arkode_ctx);
#  endif
          }
        status = ARKStepSetLinearSolver(arkode_mem, sun_linear_solver, nullptr);
        AssertARKode(status);

        auto jacobian_times_vector_callback = [](N_Vector           v,
                                                 N_Vector           Jv,
                                                 SUNDIALS::realtype t,
                                                 N_Vector           y,
                                                 N_Vector           fy,
                                                 void              *user_data,
                                                 N_Vector) -> int {
          Assert(user_data != nullptr, ExcInternalError());
          auto &callback_ctx =
            *static_cast<ARKCallbackContext<ARKStepper<VectorType>> *>(
              user_data);

          auto *src_v  = internal::unwrap_nvector_const<VectorType>(v);
          auto *src_y  = internal::unwrap_nvector_const<VectorType>(y);
          auto *src_fy = internal::unwrap_nvector_const<VectorType>(fy);

          auto *dst_Jv = internal::unwrap_nvector<VectorType>(Jv);

          return Utilities::call_and_possibly_capture_exception(
            callback_ctx.stepper->jacobian_times_vector,
            *callback_ctx.pending_exception,
            *src_v,
            *dst_Jv,
            t,
            *src_y,
            *src_fy);
        };

        auto jacobian_times_vector_setup_callback = [](SUNDIALS::realtype t,
                                                       N_Vector           y,
                                                       N_Vector           fy,
                                                       void *user_data) -> int {
          Assert(user_data != nullptr, ExcInternalError());
          auto &callback_ctx =
            *static_cast<ARKCallbackContext<ARKStepper<VectorType>> *>(
              user_data);

          auto *src_y  = internal::unwrap_nvector_const<VectorType>(y);
          auto *src_fy = internal::unwrap_nvector_const<VectorType>(fy);

          return Utilities::call_and_possibly_capture_exception(
            callback_ctx.stepper->jacobian_times_setup,
            *callback_ctx.pending_exception,
            t,
            *src_y,
            *src_fy);
        };
        status = ARKStepSetJacTimes(arkode_mem,
                                    jacobian_times_setup ?
                                      jacobian_times_vector_setup_callback :
                                      ARKLsJacTimesSetupFn(nullptr),
                                    jacobian_times_vector_callback);
        AssertARKode(status);
        if (jacobian_preconditioner_solve)
          {
            auto solve_with_jacobian_callback = [](SUNDIALS::realtype t,
                                                   N_Vector           y,
                                                   N_Vector           fy,
                                                   N_Vector           r,
                                                   N_Vector           z,
                                                   SUNDIALS::realtype gamma,
                                                   SUNDIALS::realtype delta,
                                                   int                lr,
                                                   void *user_data) -> int {
              Assert(user_data != nullptr, ExcInternalError());
              auto &callback_ctx =
                *static_cast<ARKCallbackContext<ARKStepper<VectorType>> *>(
                  user_data);

              auto *src_y  = internal::unwrap_nvector_const<VectorType>(y);
              auto *src_fy = internal::unwrap_nvector_const<VectorType>(fy);
              auto *src_r  = internal::unwrap_nvector_const<VectorType>(r);

              auto *dst_z = internal::unwrap_nvector<VectorType>(z);

              return Utilities::call_and_possibly_capture_exception(
                callback_ctx.stepper->jacobian_preconditioner_solve,
                *callback_ctx.pending_exception,
                t,
                *src_y,
                *src_fy,
                *src_r,
                *dst_z,
                gamma,
                delta,
                lr);
            };

            auto jacobian_solver_setup_callback =
              [](SUNDIALS::realtype  t,
                 N_Vector            y,
                 N_Vector            fy,
                 SUNDIALS::booltype  jok,
                 SUNDIALS::booltype *jcurPtr,
                 SUNDIALS::realtype  gamma,
                 void               *user_data) -> int {
              Assert(user_data != nullptr, ExcInternalError());
              auto &callback_ctx =
                *static_cast<ARKCallbackContext<ARKStepper<VectorType>> *>(
                  user_data);

              auto *src_y  = internal::unwrap_nvector_const<VectorType>(y);
              auto *src_fy = internal::unwrap_nvector_const<VectorType>(fy);

              return Utilities::call_and_possibly_capture_exception(
                callback_ctx.stepper->jacobian_preconditioner_setup,
                *callback_ctx.pending_exception,
                t,
                *src_y,
                *src_fy,
                jok,
                *jcurPtr,
                gamma);
            };

            status = ARKStepSetPreconditioner(arkode_mem,
                                              jacobian_preconditioner_setup ?
                                                jacobian_solver_setup_callback :
                                                ARKLsPrecSetupFn(nullptr),
                                              solve_with_jacobian_callback);
            AssertARKode(status);
          }
        if (data.implicit_function_is_linear)
          {
            status = ARKStepSetLinear(
              arkode_mem, data.implicit_function_is_time_independent ? 0 : 1);
            AssertARKode(status);
          }
      }
    else
      {
        auto y_template = internal::make_nvector_view(solution
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                                                      ,
                                                      inv_ctx.arkode_ctx
#  endif
        );

#  if DEAL_II_SUNDIALS_VERSION_LT(6, 0, 0)
        SUNNonlinearSolver fixed_point_solver =
          SUNNonlinSol_FixedPoint(y_template,
                                  data.anderson_acceleration_subspace);
#  else
        SUNNonlinearSolver fixed_point_solver =
          SUNNonlinSol_FixedPoint(y_template,
                                  data.anderson_acceleration_subspace,
                                  inv_ctx.arkode_ctx);
#  endif

        status = ARKStepSetNonlinearSolver(arkode_mem, fixed_point_solver);
        AssertARKode(status);
      }

    status =
      ARKStepSetMaxNonlinIters(arkode_mem, data.maximum_non_linear_iterations);
    AssertARKode(status);
  }

  template <typename VectorType>
  void
  ARKStepper<VectorType>::setup_mass_solver(const VectorType &solution,
                                            internal::InvocationContext inv_ctx)
  {
    int status;
    (void)status;

    if (mass_times_vector)
      {
        SUNLinearSolver sun_mass_linear_solver;
        if (solve_mass)
          {
            mass_solver =
              std::make_unique<internal::LinearSolverWrapper<VectorType>>(
                solve_mass,

                inv_ctx.pending_exception
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                ,
                inv_ctx.arkode_ctx
#  endif
              );
            sun_mass_linear_solver = *mass_solver;
          }
        else
          {
            auto y_template = internal::make_nvector_view(solution
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                                                          ,
                                                          inv_ctx.arkode_ctx
#  endif
            );
#  if DEAL_II_SUNDIALS_VERSION_LT(6, 0, 0)
            sun_mass_linear_solver =
              SUNLinSol_SPGMR(y_template,
                              PREC_NONE,
                              0 /*krylov subvectors, 0 uses default*/);
#  else
            sun_mass_linear_solver =
              SUNLinSol_SPGMR(y_template,
                              SUN_PREC_NONE,
                              0 /*krylov subvectors, 0 uses default*/,
                              inv_ctx.arkode_ctx);
#  endif
          }

        SUNDIALS::booltype mass_time_dependent =
          data.mass_is_time_independent ? SUNFALSE : SUNTRUE;

        status = ARKStepSetMassLinearSolver(arkode_mem,
                                            sun_mass_linear_solver,
                                            nullptr,
                                            mass_time_dependent);
        AssertARKode(status);

        auto mass_matrix_times_vector_setup_callback =
          [](SUNDIALS::realtype t, void *mtimes_data) -> int {
          Assert(mtimes_data != nullptr, ExcInternalError());
          auto &callback_ctx =
            *static_cast<ARKCallbackContext<ARKStepper<VectorType>> *>(
              mtimes_data);

          return Utilities::call_and_possibly_capture_exception(
            callback_ctx.stepper->mass_times_setup,
            *callback_ctx.pending_exception,
            t);
        };

        auto mass_matrix_times_vector_callback = [](N_Vector           v,
                                                    N_Vector           Mv,
                                                    SUNDIALS::realtype t,
                                                    void *mtimes_data) -> int {
          Assert(mtimes_data != nullptr, ExcInternalError());
          auto &callback_ctx =
            *static_cast<ARKCallbackContext<ARKStepper<VectorType>> *>(
              mtimes_data);

          auto *src_v  = internal::unwrap_nvector_const<VectorType>(v);
          auto *dst_Mv = internal::unwrap_nvector<VectorType>(Mv);

          return Utilities::call_and_possibly_capture_exception(
            callback_ctx.stepper->mass_times_vector,
            *callback_ctx.pending_exception,
            t,
            *src_v,
            *dst_Mv);
        };

        status = ARKStepSetMassTimes(arkode_mem,
                                     mass_times_setup ?
                                       mass_matrix_times_vector_setup_callback :
                                       ARKLsMassTimesSetupFn(nullptr),
                                     mass_matrix_times_vector_callback,
                                     this);
        AssertARKode(status);

        if (mass_preconditioner_solve)
          {
            auto mass_matrix_solver_setup_callback =
              [](SUNDIALS::realtype t, void *user_data) -> int {
              Assert(user_data != nullptr, ExcInternalError());
              auto &callback_ctx =
                *static_cast<ARKCallbackContext<ARKStepper<VectorType>> *>(
                  user_data);

              return Utilities::call_and_possibly_capture_exception(
                callback_ctx.stepper->mass_preconditioner_setup,
                *callback_ctx.pending_exception,
                t);
            };

            auto solve_with_mass_matrix_callback = [](SUNDIALS::realtype t,
                                                      N_Vector           r,
                                                      N_Vector           z,
                                                      SUNDIALS::realtype delta,
                                                      int                lr,
                                                      void *user_data) -> int {
              Assert(user_data != nullptr, ExcInternalError());
              auto &callback_ctx =
                *static_cast<ARKCallbackContext<ARKStepper<VectorType>> *>(
                  user_data);

              auto *src_r = internal::unwrap_nvector_const<VectorType>(r);
              auto *dst_z = internal::unwrap_nvector<VectorType>(z);

              return Utilities::call_and_possibly_capture_exception(
                callback_ctx.stepper->mass_preconditioner_solve,
                *callback_ctx.pending_exception,
                t,
                *src_r,
                *dst_z,
                delta,
                lr);
            };

            status =
              ARKStepSetMassPreconditioner(arkode_mem,
                                           mass_preconditioner_setup ?
                                             mass_matrix_solver_setup_callback :
                                             ARKLsMassPrecSetupFn(nullptr),
                                           solve_with_mass_matrix_callback);
            AssertARKode(status);
          }
      }
  }


  template <typename VectorType>
  void *
  ARKStepper<VectorType>::get_arkode_memory() const
  {
    return arkode_mem;
  }

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
