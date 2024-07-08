// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2024 by the deal.II authors
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

#include <deal.II/sundials/arkode.h>

#ifdef DEAL_II_WITH_SUNDIALS

#  include <deal.II/base/discrete_time.h>
#  include <deal.II/base/utilities.h>

#  include <deal.II/lac/block_vector.h>
#  include <deal.II/lac/la_parallel_block_vector.h>
#  include <deal.II/lac/la_parallel_vector.h>
#  ifdef DEAL_II_WITH_TRILINOS
#    include <deal.II/lac/trilinos_parallel_block_vector.h>
#    include <deal.II/lac/trilinos_vector.h>
#  endif
#  ifdef DEAL_II_WITH_PETSC
#    include <deal.II/lac/petsc_block_vector.h>
#    include <deal.II/lac/petsc_vector.h>
#  endif

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
  ARKode<VectorType>::ARKode(const AdditionalData &data)
    : ARKode(data, MPI_COMM_SELF)
  {}


  template <typename VectorType>
  ARKode<VectorType>::ARKode(const AdditionalData &data,
                             const MPI_Comm        mpi_comm)
    : data(data)
    , arkode_mem(nullptr)
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    , arkode_ctx(nullptr)
#  endif
    , mpi_communicator(mpi_comm)
    , last_end_time(data.initial_time)
    , pending_exception(nullptr)
  {
    set_functions_to_trigger_an_assert();

    // SUNDIALS will always duplicate communicators if we provide them. This
    // can cause problems if SUNDIALS is configured with MPI and we pass along
    // MPI_COMM_SELF in a serial application as MPI won't be
    // initialized. Hence, work around that by just not providing a
    // communicator in that case.
#  if DEAL_II_SUNDIALS_VERSION_GTE(7, 0, 0)
    const int status =
      SUNContext_Create(mpi_communicator == MPI_COMM_SELF ? SUN_COMM_NULL :
                                                            mpi_communicator,
                        &arkode_ctx);
    (void)status;
    AssertARKode(status);
#  elif DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    const int status =
      SUNContext_Create(mpi_communicator == MPI_COMM_SELF ? nullptr :
                                                            &mpi_communicator,
                        &arkode_ctx);
    (void)status;
    AssertARKode(status);
#  endif
  }



  template <typename VectorType>
  ARKode<VectorType>::~ARKode()
  {
    ARKStepFree(&arkode_mem);

#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    const int status = SUNContext_Free(&arkode_ctx);
    (void)status;
    AssertARKode(status);
#  endif

    Assert(pending_exception == nullptr, ExcInternalError());
  }



  template <typename VectorType>
  unsigned int
  ARKode<VectorType>::solve_ode(VectorType &solution)
  {
    DiscreteTime time(data.initial_time,
                      data.final_time,
                      data.initial_step_size);

    return do_evolve_time(solution, time, /* force_solver_restart = */ true);
  }



  template <typename VectorType>
  unsigned int
  ARKode<VectorType>::solve_ode_incrementally(VectorType  &solution,
                                              const double intermediate_time,
                                              const bool   reset_solver)
  {
    AssertThrow(
      intermediate_time > last_end_time,
      dealii::ExcMessage(
        "The requested intermediate time is smaller than the last requested "
        "intermediate time."));

    const bool   do_reset = reset_solver || arkode_mem == nullptr;
    DiscreteTime time(last_end_time, intermediate_time, data.initial_step_size);
    return do_evolve_time(solution, time, do_reset);
  }



  template <typename VectorType>
  unsigned int
  ARKode<VectorType>::do_evolve_time(VectorType   &solution,
                                     DiscreteTime &time,
                                     const bool    do_reset)
  {
    if (do_reset)
      {
        reset(time.get_current_time(), time.get_next_step_size(), solution);
        if (output_step)
          output_step(time.get_current_time(),
                      solution,
                      time.get_step_number());
      }
    else
      {
        // If we don't do a full reset then we still need to fix the end time.
        // In SUNDIALS 6 and later, SUNDIALS will not do timesteps if the
        // current time is past the set end point (i.e., ARKStepEvolve will
        // return ARK_TSTOP_RETURN).
        const int status = ARKStepSetStopTime(arkode_mem, time.get_end_time());
        (void)status;
        AssertARKode(status);
      }

    auto solution_nvector = internal::make_nvector_view(solution
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                                                        ,
                                                        arkode_ctx
#  endif
    );

    while (!time.is_at_end())
      {
        time.set_desired_next_step_size(data.output_period);

        // Having set up all of the ancillary things, finally call the main
        // ARKode function. Once we return, check what happened:
        // - If we have a pending recoverable exception, ignore it if SUNDIAL's
        //   return code was zero -- in that case, SUNDIALS managed to indeed
        //   recover and we no longer need the exception
        // - If we have any other exception, rethrow it
        // - If no exception, test that SUNDIALS really did successfully return
        Assert(pending_exception == nullptr, ExcInternalError());
        double     actual_next_time;
        const auto status = ARKStepEvolve(arkode_mem,
                                          time.get_next_time(),
                                          solution_nvector,
                                          &actual_next_time,
                                          ARK_NORMAL);
        if (pending_exception)
          {
            try
              {
                std::rethrow_exception(pending_exception);
              }
            catch (const RecoverableUserCallbackError &exc)
              {
                pending_exception = nullptr;
                if (status == 0)
                  /* just eat the exception */;
                else
                  throw;
              }
            catch (...)
              {
                pending_exception = nullptr;
                throw;
              }
          }

        AssertARKode(status);

        // Then reflect this time advancement in our own DiscreteTime object:
        time.set_next_step_size(actual_next_time - time.get_current_time());
        time.advance_time();

        // Finally check whether resets or output calls are desired at this
        // time:
        while (solver_should_restart(time.get_current_time(), solution))
          reset(time.get_current_time(),
                time.get_previous_step_size(),
                solution);

        if (output_step)
          output_step(time.get_current_time(),
                      solution,
                      time.get_step_number());
      }
    last_end_time = time.get_current_time();

    long int   n_steps;
    const auto status = ARKStepGetNumSteps(arkode_mem, &n_steps);
    (void)status;
    AssertARKode(status);

    return n_steps;
  }



  template <typename VectorType>
  void
  ARKode<VectorType>::reset(const double      current_time,
                            const double      current_time_step,
                            const VectorType &solution)
  {
    last_end_time = current_time;
    int status;
    (void)status;

#  if DEAL_II_SUNDIALS_VERSION_GTE(7, 0, 0)
    status = SUNContext_Free(&arkode_ctx);
    AssertARKode(status);

    // Same comment applies as in class constructor:
    status =
      SUNContext_Create(mpi_communicator == MPI_COMM_SELF ? SUN_COMM_NULL :
                                                            mpi_communicator,
                        &arkode_ctx);
    AssertARKode(status);
#  elif DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    status = SUNContext_Free(&arkode_ctx);
    AssertARKode(status);

    // Same comment applies as in class constructor:
    status =
      SUNContext_Create(mpi_communicator == MPI_COMM_SELF ? nullptr :
                                                            &mpi_communicator,
                        &arkode_ctx);
    AssertARKode(status);
#  endif

    if (arkode_mem)
      {
        ARKStepFree(&arkode_mem);
        // Initialization is version-dependent: do that in a moment
      }

    // just a view on the memory in solution, all write operations on yy by
    // ARKODE will automatically be mirrored to solution
    auto initial_condition_nvector = internal::make_nvector_view(solution
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                                                                 ,
                                                                 arkode_ctx
#  endif
    );

    Assert(explicit_function || implicit_function,
           ExcFunctionNotProvided("explicit_function || implicit_function"));

    auto explicit_function_callback = [](SUNDIALS::realtype tt,
                                         N_Vector           yy,
                                         N_Vector           yp,
                                         void              *user_data) -> int {
      Assert(user_data != nullptr, ExcInternalError());
      ARKode<VectorType> &solver =
        *static_cast<ARKode<VectorType> *>(user_data);

      auto *src_yy = internal::unwrap_nvector_const<VectorType>(yy);
      auto *dst_yp = internal::unwrap_nvector<VectorType>(yp);

      return Utilities::call_and_possibly_capture_exception(
        solver.explicit_function,
        solver.pending_exception,
        tt,
        *src_yy,
        *dst_yp);
    };


    auto implicit_function_callback = [](SUNDIALS::realtype tt,
                                         N_Vector           yy,
                                         N_Vector           yp,
                                         void              *user_data) -> int {
      Assert(user_data != nullptr, ExcInternalError());
      ARKode<VectorType> &solver =
        *static_cast<ARKode<VectorType> *>(user_data);

      auto *src_yy = internal::unwrap_nvector_const<VectorType>(yy);
      auto *dst_yp = internal::unwrap_nvector<VectorType>(yp);

      return Utilities::call_and_possibly_capture_exception(
        solver.implicit_function,
        solver.pending_exception,
        tt,
        *src_yy,
        *dst_yp);
    };

    arkode_mem = ARKStepCreate(explicit_function ? explicit_function_callback :
                                                   ARKRhsFn(nullptr),
                               implicit_function ? implicit_function_callback :
                                                   ARKRhsFn(nullptr),
                               current_time,
                               initial_condition_nvector
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                               ,
                               arkode_ctx
#  endif
    );
    Assert(arkode_mem != nullptr, ExcInternalError());

    if (get_local_tolerances)
      {
        const auto abs_tols = internal::make_nvector_view(get_local_tolerances()
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                                                            ,
                                                          arkode_ctx
#  endif
        );
        status =
          ARKStepSVtolerances(arkode_mem, data.relative_tolerance, abs_tols);
        AssertARKode(status);
      }
    else
      {
        status = ARKStepSStolerances(arkode_mem,
                                     data.relative_tolerance,
                                     data.absolute_tolerance);
        AssertARKode(status);
      }

    // for version 4.0.0 this call must be made before solver settings
    status = ARKStepSetUserData(arkode_mem, this);
    AssertARKode(status);

    setup_system_solver(solution);

    setup_mass_solver(solution);

    status = ARKStepSetInitStep(arkode_mem, current_time_step);
    AssertARKode(status);

    status = ARKStepSetStopTime(arkode_mem, data.final_time);
    AssertARKode(status);

    status = ARKStepSetOrder(arkode_mem, data.maximum_order);
    AssertARKode(status);

    if (custom_setup)
      custom_setup(arkode_mem);
  }



  template <typename VectorType>
  void
  ARKode<VectorType>::setup_system_solver(const VectorType &solution)
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
                pending_exception
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                ,
                arkode_ctx
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
                                                          arkode_ctx
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
                              arkode_ctx);
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
          ARKode<VectorType> &solver =
            *static_cast<ARKode<VectorType> *>(user_data);

          auto *src_v  = internal::unwrap_nvector_const<VectorType>(v);
          auto *src_y  = internal::unwrap_nvector_const<VectorType>(y);
          auto *src_fy = internal::unwrap_nvector_const<VectorType>(fy);

          auto *dst_Jv = internal::unwrap_nvector<VectorType>(Jv);

          return Utilities::call_and_possibly_capture_exception(
            solver.jacobian_times_vector,
            solver.pending_exception,
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
          ARKode<VectorType> &solver =
            *static_cast<ARKode<VectorType> *>(user_data);

          auto *src_y  = internal::unwrap_nvector_const<VectorType>(y);
          auto *src_fy = internal::unwrap_nvector_const<VectorType>(fy);

          return Utilities::call_and_possibly_capture_exception(
            solver.jacobian_times_setup,
            solver.pending_exception,
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
              ARKode<VectorType> &solver =
                *static_cast<ARKode<VectorType> *>(user_data);

              auto *src_y  = internal::unwrap_nvector_const<VectorType>(y);
              auto *src_fy = internal::unwrap_nvector_const<VectorType>(fy);
              auto *src_r  = internal::unwrap_nvector_const<VectorType>(r);

              auto *dst_z = internal::unwrap_nvector<VectorType>(z);

              return Utilities::call_and_possibly_capture_exception(
                solver.jacobian_preconditioner_solve,
                solver.pending_exception,
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
              ARKode<VectorType> &solver =
                *static_cast<ARKode<VectorType> *>(user_data);

              auto *src_y  = internal::unwrap_nvector_const<VectorType>(y);
              auto *src_fy = internal::unwrap_nvector_const<VectorType>(fy);

              return Utilities::call_and_possibly_capture_exception(
                solver.jacobian_preconditioner_setup,
                solver.pending_exception,
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
                                                      arkode_ctx
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
                                  arkode_ctx);
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
  ARKode<VectorType>::setup_mass_solver(const VectorType &solution)
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

                pending_exception
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                ,
                arkode_ctx
#  endif
              );
            sun_mass_linear_solver = *mass_solver;
          }
        else
          {
            auto y_template = internal::make_nvector_view(solution
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
                                                          ,
                                                          arkode_ctx
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
                              arkode_ctx);
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
          ARKode<VectorType> &solver =
            *static_cast<ARKode<VectorType> *>(mtimes_data);

          return Utilities::call_and_possibly_capture_exception(
            solver.mass_times_setup, solver.pending_exception, t);
        };

        auto mass_matrix_times_vector_callback = [](N_Vector           v,
                                                    N_Vector           Mv,
                                                    SUNDIALS::realtype t,
                                                    void *mtimes_data) -> int {
          Assert(mtimes_data != nullptr, ExcInternalError());
          ARKode<VectorType> &solver =
            *static_cast<ARKode<VectorType> *>(mtimes_data);

          auto *src_v  = internal::unwrap_nvector_const<VectorType>(v);
          auto *dst_Mv = internal::unwrap_nvector<VectorType>(Mv);

          return Utilities::call_and_possibly_capture_exception(
            solver.mass_times_vector,
            solver.pending_exception,
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
              ARKode<VectorType> &solver =
                *static_cast<ARKode<VectorType> *>(user_data);

              return Utilities::call_and_possibly_capture_exception(
                solver.mass_preconditioner_setup, solver.pending_exception, t);
            };

            auto solve_with_mass_matrix_callback = [](SUNDIALS::realtype t,
                                                      N_Vector           r,
                                                      N_Vector           z,
                                                      SUNDIALS::realtype delta,
                                                      int                lr,
                                                      void *user_data) -> int {
              Assert(user_data != nullptr, ExcInternalError());
              ARKode<VectorType> &solver =
                *static_cast<ARKode<VectorType> *>(user_data);

              auto *src_r = internal::unwrap_nvector_const<VectorType>(r);
              auto *dst_z = internal::unwrap_nvector<VectorType>(z);

              return Utilities::call_and_possibly_capture_exception(
                solver.mass_preconditioner_solve,
                solver.pending_exception,
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
  void
  ARKode<VectorType>::set_functions_to_trigger_an_assert()
  {
    solver_should_restart = [](const double, VectorType &) -> bool {
      return false;
    };
  }



  template <typename VectorType>
  void *
  ARKode<VectorType>::get_arkode_memory() const
  {
    return arkode_mem;
  }

  template class ARKode<Vector<double>>;
  template class ARKode<BlockVector<double>>;

  template class ARKode<LinearAlgebra::distributed::Vector<double>>;
  template class ARKode<LinearAlgebra::distributed::BlockVector<double>>;

#  ifdef DEAL_II_WITH_MPI

#    ifdef DEAL_II_WITH_TRILINOS
  template class ARKode<TrilinosWrappers::MPI::Vector>;
  template class ARKode<TrilinosWrappers::MPI::BlockVector>;

#    endif // DEAL_II_WITH_TRILINOS

#    ifdef DEAL_II_WITH_PETSC
#      ifndef PETSC_USE_COMPLEX
  template class ARKode<PETScWrappers::MPI::Vector>;
  template class ARKode<PETScWrappers::MPI::BlockVector>;
#      endif // PETSC_USE_COMPLEX
#    endif   // DEAL_II_WITH_PETSC

#  endif // DEAL_II_WITH_MPI

} // namespace SUNDIALS

DEAL_II_NAMESPACE_CLOSE

#endif
