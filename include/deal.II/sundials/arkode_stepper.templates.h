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
#  include <arkode/arkode_erkstep.h>
#  if DEAL_II_SUNDIALS_VERSION_GTE(7, 2, 0)
#    include <arkode/arkode_lsrkstep.h>
#  endif
#  include <sunlinsol/sunlinsol_spgmr.h>
#  include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

#  include <iostream>

#endif // DEAL_II_WITH_SUNDIALS

DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_SUNDIALS

namespace SUNDIALS
{
  template <typename VectorType>
  ARKStepper<VectorType>::ARKStepper(const AdditionalData &data)
    : data(data)
    , arkode_mem(nullptr, [](void *mem) {
      if (mem)
        ARKStepFree(&mem);
    })
  {
#  if DEAL_II_SUNDIALS_VERSION_LT(6, 4, 0)
    AssertThrow(
      data.implicit_butcher_table.empty() &&
        data.explicit_butcher_table.empty(),
      ExcMessage(
        "Setting ARK Butcher tables by name requires SUNDIALS version 6.4.0 "
        "or later. Set the order of accuracy instead or use custom_setup "
        "callback to set the Butcher tables manually using "
        "function ARKStepSetTableNum()."));
#  endif

    AssertThrow(
      data.order == 0 || (data.implicit_butcher_table.empty() &&
                          data.explicit_butcher_table.empty()),
      ExcMessage(
        "Either the order of accuracy or the Butcher table names may be "
        "specified, but not both."));
  }



  template <typename VectorType>
  void
  ARKStepper<VectorType>::reinit(double                      t0,
                                 const VectorType           &y0,
                                 internal::InvocationContext inv_ctx)
  {
    arkode_mem.reset();

    Assert(explicit_function || implicit_function,
           ExcFunctionNotProvided("explicit_function || implicit_function"));

    const auto explicit_function_callback = [](SUNDIALS::realtype tt,
                                               N_Vector           yy,
                                               N_Vector           yp,
                                               void *user_data) -> int {
      Assert(user_data != nullptr, ExcInternalError());
      auto &callback_ctx = *static_cast<CallbackContext *>(user_data);

      auto *src_yy = internal::unwrap_nvector_const<VectorType>(yy);
      auto *dst_yp = internal::unwrap_nvector<VectorType>(yp);

      return Utilities::call_and_possibly_capture_exception(
        callback_ctx.stepper->explicit_function,
        *callback_ctx.pending_exception,
        tt,
        *src_yy,
        *dst_yp);
    };

    const auto implicit_function_callback = [](SUNDIALS::realtype tt,
                                               N_Vector           yy,
                                               N_Vector           yp,
                                               void *user_data) -> int {
      Assert(user_data != nullptr, ExcInternalError());
      auto &callback_ctx = *static_cast<CallbackContext *>(user_data);

      auto *src_yy = internal::unwrap_nvector_const<VectorType>(yy);
      auto *dst_yp = internal::unwrap_nvector<VectorType>(yp);

      return Utilities::call_and_possibly_capture_exception(
        callback_ctx.stepper->implicit_function,
        *callback_ctx.pending_exception,
        tt,
        *src_yy,
        *dst_yp);
    };

    // just a view on the memory in solution, all write operations on yy by
    // ARKODE will automatically be mirrored to solution
    const auto initial_condition_nvector =
      internal::make_nvector_view(y0
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                                  ,
                                  inv_ctx.arkode_ctx
#  endif
      );

    arkode_mem.reset(ARKStepCreate(
      explicit_function ? explicit_function_callback : ARKRhsFn(nullptr),
      implicit_function ? implicit_function_callback : ARKRhsFn(nullptr),
      t0,
      initial_condition_nvector
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
      ,
      inv_ctx.arkode_ctx
#  endif
      ));
    Assert(arkode_mem != nullptr, ExcInternalError());

    callback_ctx = {this, &inv_ctx.pending_exception};
    int status   = ARKStepSetUserData(arkode_mem.get(), &callback_ctx);
    AssertARKode(status);

    setup_system_solver(y0, inv_ctx);

    setup_mass_solver(y0, inv_ctx);

    if (data.order != 0)
      {
#  if DEAL_II_SUNDIALS_VERSION_GTE(7, 1, 0)
        status = ARKodeSetOrder(arkode_mem.get(), data.order);
#  else
        status = ARKStepSetOrder(arkode_mem.get(), data.order);
#  endif
        AssertARKode(status);
      }
    else if (!data.implicit_butcher_table.empty() ||
             !data.explicit_butcher_table.empty())
      {
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 4, 0)
        const std::string itable = data.implicit_butcher_table.empty() ?
                                     "ARKODE_DIRK_NONE" :
                                     data.implicit_butcher_table;
        const std::string etable = data.explicit_butcher_table.empty() ?
                                     "ARKODE_ERK_NONE" :
                                     data.explicit_butcher_table;
        status =
          ARKStepSetTableName(arkode_mem.get(), itable.c_str(), etable.c_str());
        AssertARKode(status);
#  endif
      }

    if (custom_setup)
      custom_setup(arkode_mem.get());
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
        status =
          ARKStepSetLinearSolver(arkode_mem.get(), sun_linear_solver, nullptr);
        AssertARKode(status);

        auto jacobian_times_vector_callback = [](N_Vector           v,
                                                 N_Vector           Jv,
                                                 SUNDIALS::realtype t,
                                                 N_Vector           y,
                                                 N_Vector           fy,
                                                 void              *user_data,
                                                 N_Vector) -> int {
          Assert(user_data != nullptr, ExcInternalError());
          auto &callback_ctx = *static_cast<CallbackContext *>(user_data);

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
          auto &callback_ctx = *static_cast<CallbackContext *>(user_data);

          auto *src_y  = internal::unwrap_nvector_const<VectorType>(y);
          auto *src_fy = internal::unwrap_nvector_const<VectorType>(fy);

          return Utilities::call_and_possibly_capture_exception(
            callback_ctx.stepper->jacobian_times_vector_setup,
            *callback_ctx.pending_exception,
            t,
            *src_y,
            *src_fy);
        };
        status = ARKStepSetJacTimes(arkode_mem.get(),
                                    jacobian_times_vector_setup ?
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
              auto &callback_ctx = *static_cast<CallbackContext *>(user_data);

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
              auto &callback_ctx = *static_cast<CallbackContext *>(user_data);

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

            status = ARKStepSetPreconditioner(arkode_mem.get(),
                                              jacobian_preconditioner_setup ?
                                                jacobian_solver_setup_callback :
                                                ARKLsPrecSetupFn(nullptr),
                                              solve_with_jacobian_callback);
            AssertARKode(status);
          }
        if (data.implicit_function_is_linear)
          {
            status =
              ARKStepSetLinear(arkode_mem.get(),
                               data.implicit_function_is_time_independent ? 0 :
                                                                            1);
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

        status =
          ARKStepSetNonlinearSolver(arkode_mem.get(), fixed_point_solver);
        AssertARKode(status);
      }

    status = ARKStepSetMaxNonlinIters(arkode_mem.get(),
                                      data.maximum_non_linear_iterations);
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

        status = ARKStepSetMassLinearSolver(arkode_mem.get(),
                                            sun_mass_linear_solver,
                                            nullptr,
                                            mass_time_dependent);
        AssertARKode(status);

        auto mass_matrix_times_vector_setup_callback =
          [](SUNDIALS::realtype t, void *mtimes_data) -> int {
          Assert(mtimes_data != nullptr, ExcInternalError());
          auto &callback_ctx = *static_cast<CallbackContext *>(mtimes_data);

          return Utilities::call_and_possibly_capture_exception(
            callback_ctx.stepper->mass_times_vector_setup,
            *callback_ctx.pending_exception,
            t);
        };

        auto mass_matrix_times_vector_callback = [](N_Vector           v,
                                                    N_Vector           Mv,
                                                    SUNDIALS::realtype t,
                                                    void *mtimes_data) -> int {
          Assert(mtimes_data != nullptr, ExcInternalError());
          auto &callback_ctx = *static_cast<CallbackContext *>(mtimes_data);

          auto *src_v  = internal::unwrap_nvector_const<VectorType>(v);
          auto *dst_Mv = internal::unwrap_nvector<VectorType>(Mv);

          return Utilities::call_and_possibly_capture_exception(
            callback_ctx.stepper->mass_times_vector,
            *callback_ctx.pending_exception,
            t,
            *src_v,
            *dst_Mv);
        };

        status = ARKStepSetMassTimes(arkode_mem.get(),
                                     mass_times_vector_setup ?
                                       mass_matrix_times_vector_setup_callback :
                                       ARKLsMassTimesSetupFn(nullptr),
                                     mass_matrix_times_vector_callback,
                                     &callback_ctx);
        AssertARKode(status);

        if (mass_preconditioner_solve)
          {
            auto mass_matrix_solver_setup_callback =
              [](SUNDIALS::realtype t, void *user_data) -> int {
              Assert(user_data != nullptr, ExcInternalError());
              auto &callback_ctx = *static_cast<CallbackContext *>(user_data);

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
              auto &callback_ctx = *static_cast<CallbackContext *>(user_data);

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
              ARKStepSetMassPreconditioner(arkode_mem.get(),
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
    return arkode_mem.get();
  }


  // ---------------------------------------------------------------------------
  // ERKStepper implementation
  // ---------------------------------------------------------------------------

  template <typename VectorType>
  ERKStepper<VectorType>::ERKStepper(const AdditionalData &data)
    : data(data)
    , arkode_mem(nullptr, [](void *mem) {
      if (mem)
        ERKStepFree(&mem);
    })
  {
#  if DEAL_II_SUNDIALS_VERSION_LT(6, 4, 0)
    AssertThrow(
      data.explicit_butcher_table.empty(),
      ExcMessage(
        "Setting a ERK Butcher table by name requires SUNDIALS version 6.4.0 "
        "or later. Set the order of accuracy instead or use custom_setup "
        "callback to set the Butcher tables manually using "
        "function ERKStepSetTableNum()."));
#  endif

    AssertThrow(
      data.order == 0 || data.explicit_butcher_table.empty(),
      ExcMessage(
        "Either the order of accuracy or the Butcher table name may be "
        "specified, but not both."));
  }



  template <typename VectorType>
  void
  ERKStepper<VectorType>::reinit(double                      t0,
                                 const VectorType           &y0,
                                 internal::InvocationContext inv_ctx)
  {
    arkode_mem.reset();

    Assert(explicit_function, ExcFunctionNotProvided("explicit_function"));

    auto explicit_function_callback = [](SUNDIALS::realtype tt,
                                         N_Vector           yy,
                                         N_Vector           yp,
                                         void              *user_data) -> int {
      Assert(user_data != nullptr, ExcInternalError());
      auto &callback_ctx = *static_cast<CallbackContext *>(user_data);

      auto *src_yy = internal::unwrap_nvector_const<VectorType>(yy);
      auto *dst_yp = internal::unwrap_nvector<VectorType>(yp);

      return Utilities::call_and_possibly_capture_exception(
        callback_ctx.stepper->explicit_function,
        *callback_ctx.pending_exception,
        tt,
        *src_yy,
        *dst_yp);
    };

    auto initial_condition_nvector =
      internal::make_nvector_view(y0
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                                  ,
                                  inv_ctx.arkode_ctx
#  endif
      );

    arkode_mem.reset(ERKStepCreate(explicit_function_callback,
                                   t0,
                                   initial_condition_nvector
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                                   ,
                                   inv_ctx.arkode_ctx
#  endif
                                   ));
    Assert(arkode_mem != nullptr, ExcInternalError());

    callback_ctx = {this, &inv_ctx.pending_exception};
    int status   = ERKStepSetUserData(arkode_mem.get(), &callback_ctx);
    AssertARKode(status);

    if (data.order != 0)
      {
#  if DEAL_II_SUNDIALS_VERSION_GTE(7, 1, 0)
        status = ARKodeSetOrder(arkode_mem.get(), data.order);
#  else
        status = ERKStepSetOrder(arkode_mem.get(), data.order);
#  endif
        AssertARKode(status);
      }
    else if (!data.explicit_butcher_table.empty())
      {
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 4, 0)
        status = ERKStepSetTableName(arkode_mem.get(),
                                     data.explicit_butcher_table.c_str());
        AssertARKode(status);
#  endif
      }

    if (custom_setup)
      custom_setup(arkode_mem.get());
  }


  template <typename VectorType>
  void *
  ERKStepper<VectorType>::get_arkode_memory() const
  {
    return arkode_mem.get();
  }


#  if DEAL_II_SUNDIALS_VERSION_GTE(7, 2, 0)

  // --------------------------------------------------------------------------
  // LSRKStepperSTS
  // --------------------------------------------------------------------------

  template <typename VectorType>
  LSRKStepperSTS<VectorType>::LSRKStepperSTS(const AdditionalData &data)
    : data(data)
    , arkode_mem(nullptr, [](void *mem) {
      if (mem)
        ARKodeFree(&mem);
    })
  {}


  template <typename VectorType>
  void
  LSRKStepperSTS<VectorType>::reinit(double                      t0,
                                     const VectorType           &y0,
                                     internal::InvocationContext inv_ctx)
  {
    arkode_mem.reset();

    Assert(explicit_function, ExcFunctionNotProvided("explicit_function"));
    Assert(dominant_eigenvalue_function,
           ExcFunctionNotProvided("dominant_eigenvalue_function"));

    auto explicit_function_callback = [](SUNDIALS::realtype tt,
                                         N_Vector           yy,
                                         N_Vector           yp,
                                         void              *user_data) -> int {
      Assert(user_data != nullptr, ExcInternalError());
      auto &ctx    = *static_cast<CallbackContext *>(user_data);
      auto *src_yy = internal::unwrap_nvector_const<VectorType>(yy);
      auto *dst_yp = internal::unwrap_nvector<VectorType>(yp);
      return Utilities::call_and_possibly_capture_exception(
        ctx.stepper->explicit_function,
        *ctx.pending_exception,
        tt,
        *src_yy,
        *dst_yp);
    };

    auto dom_eig_callback = [](SUNDIALS::realtype  tt,
                               N_Vector            yy,
                               N_Vector            fn,
                               SUNDIALS::realtype *lambdaR,
                               SUNDIALS::realtype *lambdaI,
                               void               *user_data,
                               N_Vector /*temp1*/,
                               N_Vector /*temp2*/,
                               N_Vector /*temp3*/) -> int {
      Assert(user_data != nullptr, ExcInternalError());
      auto &ctx    = *static_cast<CallbackContext *>(user_data);
      auto *src_yy = internal::unwrap_nvector_const<VectorType>(yy);
      auto *src_fn = internal::unwrap_nvector_const<VectorType>(fn);

      std::complex<double> lambda;

      const int ret = Utilities::call_and_possibly_capture_exception(
        [&]() {
          lambda =
            ctx.stepper->dominant_eigenvalue_function(tt, *src_yy, *src_fn);
        },
        *ctx.pending_exception);

      *lambdaR = static_cast<SUNDIALS::realtype>(lambda.real());
      *lambdaI = static_cast<SUNDIALS::realtype>(lambda.imag());
      return ret;
    };

    auto initial_condition_nvector =
      internal::make_nvector_view(y0, inv_ctx.arkode_ctx);

    arkode_mem.reset(LSRKStepCreateSTS(explicit_function_callback,
                                       t0,
                                       initial_condition_nvector,
                                       inv_ctx.arkode_ctx));
    Assert(arkode_mem != nullptr, ExcInternalError());

    callback_ctx = {this, &inv_ctx.pending_exception};
    int status   = ARKodeSetUserData(arkode_mem.get(), &callback_ctx);
    AssertARKode(status);

    status = LSRKStepSetDomEigFn(arkode_mem.get(), dom_eig_callback);
    AssertARKode(status);

    if (data.dom_eig_frequency >= 0)
      {
        status =
          LSRKStepSetDomEigFrequency(arkode_mem.get(), data.dom_eig_frequency);
        AssertARKode(status);
      }

    if (data.max_num_stages > 0)
      {
        status = LSRKStepSetMaxNumStages(arkode_mem.get(),
                                         static_cast<int>(data.max_num_stages));
        AssertARKode(status);
      }

    if (!data.method_name.empty())
      {
        status = LSRKStepSetSTSMethodByName(arkode_mem.get(),
                                            data.method_name.c_str());
        AssertARKode(status);
      }

    if (custom_setup)
      custom_setup(arkode_mem.get());
  }


  template <typename VectorType>
  void *
  LSRKStepperSTS<VectorType>::get_arkode_memory() const
  {
    return arkode_mem.get();
  }


  // --------------------------------------------------------------------------
  // LSRKStepperSSP
  // --------------------------------------------------------------------------

  template <typename VectorType>
  LSRKStepperSSP<VectorType>::LSRKStepperSSP(const AdditionalData &data)
    : data(data)
    , arkode_mem(nullptr, [](void *mem) {
      if (mem)
        ARKodeFree(&mem);
    })
  {}


  template <typename VectorType>
  void
  LSRKStepperSSP<VectorType>::reinit(double                      t0,
                                     const VectorType           &y0,
                                     internal::InvocationContext inv_ctx)
  {
    arkode_mem.reset();

    Assert(explicit_function, ExcFunctionNotProvided("explicit_function"));

    auto explicit_function_callback = [](SUNDIALS::realtype tt,
                                         N_Vector           yy,
                                         N_Vector           yp,
                                         void              *user_data) -> int {
      Assert(user_data != nullptr, ExcInternalError());
      auto &ctx    = *static_cast<CallbackContext *>(user_data);
      auto *src_yy = internal::unwrap_nvector_const<VectorType>(yy);
      auto *dst_yp = internal::unwrap_nvector<VectorType>(yp);
      return Utilities::call_and_possibly_capture_exception(
        ctx.stepper->explicit_function,
        *ctx.pending_exception,
        tt,
        *src_yy,
        *dst_yp);
    };

    auto initial_condition_nvector =
      internal::make_nvector_view(y0, inv_ctx.arkode_ctx);

    arkode_mem.reset(LSRKStepCreateSSP(explicit_function_callback,
                                       t0,
                                       initial_condition_nvector,
                                       inv_ctx.arkode_ctx));
    Assert(arkode_mem != nullptr, ExcInternalError());

    callback_ctx = {this, &inv_ctx.pending_exception};
    int status   = ARKodeSetUserData(arkode_mem.get(), &callback_ctx);
    AssertARKode(status);

    if (!data.method_name.empty())
      {
        status = LSRKStepSetSSPMethodByName(arkode_mem.get(),
                                            data.method_name.c_str());
        AssertARKode(status);
      }

    if (data.num_stages > 0)
      {
        status = LSRKStepSetNumSSPStages(arkode_mem.get(),
                                         static_cast<int>(data.num_stages));
        AssertARKode(status);
      }

    if (custom_setup)
      custom_setup(arkode_mem.get());
  }


  template <typename VectorType>
  void *
  LSRKStepperSSP<VectorType>::get_arkode_memory() const
  {
    return arkode_mem.get();
  }

#  endif // DEAL_II_SUNDIALS_VERSION_GTE(7, 2, 0)

} // namespace SUNDIALS

#endif // DEAL_II_WITH_SUNDIALS

DEAL_II_NAMESPACE_CLOSE

#endif
