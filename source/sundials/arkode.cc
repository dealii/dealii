// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



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

#  include <arkode/arkode_arkstep.h>
#  include <sunlinsol/sunlinsol_spgmr.h>
#  include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

#  include <iostream>

DEAL_II_NAMESPACE_OPEN

namespace SUNDIALS
{
  namespace
  {
    template <typename VectorType>
    int
    explicit_function_callback(realtype tt,
                               N_Vector yy,
                               N_Vector yp,
                               void *   user_data)
    {
      Assert(user_data != nullptr, ExcInternalError());
      ARKode<VectorType> &solver =
        *static_cast<ARKode<VectorType> *>(user_data);

      auto *src_yy = internal::unwrap_nvector_const<VectorType>(yy);
      auto *dst_yp = internal::unwrap_nvector<VectorType>(yp);

      return solver.explicit_function(tt, *src_yy, *dst_yp);
    }



    template <typename VectorType>
    int
    implicit_function_callback(realtype tt,
                               N_Vector yy,
                               N_Vector yp,
                               void *   user_data)
    {
      Assert(user_data != nullptr, ExcInternalError());
      ARKode<VectorType> &solver =
        *static_cast<ARKode<VectorType> *>(user_data);

      auto *src_yy = internal::unwrap_nvector_const<VectorType>(yy);
      auto *dst_yp = internal::unwrap_nvector<VectorType>(yp);

      return solver.implicit_function(tt, *src_yy, *dst_yp);
    }



    template <typename VectorType>
    int
    jacobian_times_vector_callback(N_Vector v,
                                   N_Vector Jv,
                                   realtype t,
                                   N_Vector y,
                                   N_Vector fy,
                                   void *   user_data,
                                   N_Vector)
    {
      Assert(user_data != nullptr, ExcInternalError());
      ARKode<VectorType> &solver =
        *static_cast<ARKode<VectorType> *>(user_data);

      auto *src_v  = internal::unwrap_nvector_const<VectorType>(v);
      auto *src_y  = internal::unwrap_nvector_const<VectorType>(y);
      auto *src_fy = internal::unwrap_nvector_const<VectorType>(fy);

      auto *dst_Jv = internal::unwrap_nvector<VectorType>(Jv);

      return solver.jacobian_times_vector(*src_v, *dst_Jv, t, *src_y, *src_fy);
    }



    template <typename VectorType>
    int
    jacobian_times_vector_setup_callback(realtype t,
                                         N_Vector y,
                                         N_Vector fy,
                                         void *   user_data)
    {
      Assert(user_data != nullptr, ExcInternalError());
      ARKode<VectorType> &solver =
        *static_cast<ARKode<VectorType> *>(user_data);

      auto *src_y  = internal::unwrap_nvector_const<VectorType>(y);
      auto *src_fy = internal::unwrap_nvector_const<VectorType>(fy);

      return solver.jacobian_times_setup(t, *src_y, *src_fy);
    }



    template <typename VectorType>
    int
    solve_with_jacobian_callback(realtype t,
                                 N_Vector y,
                                 N_Vector fy,
                                 N_Vector r,
                                 N_Vector z,
                                 realtype gamma,
                                 realtype delta,
                                 int      lr,
                                 void *   user_data)
    {
      Assert(user_data != nullptr, ExcInternalError());
      ARKode<VectorType> &solver =
        *static_cast<ARKode<VectorType> *>(user_data);

      auto *src_y  = internal::unwrap_nvector_const<VectorType>(y);
      auto *src_fy = internal::unwrap_nvector_const<VectorType>(fy);
      auto *src_r  = internal::unwrap_nvector_const<VectorType>(r);

      auto *dst_z = internal::unwrap_nvector<VectorType>(z);

      return solver.jacobian_preconditioner_solve(
        t, *src_y, *src_fy, *src_r, *dst_z, gamma, delta, lr);
    }



    template <typename VectorType>
    int
    jacobian_solver_setup_callback(realtype     t,
                                   N_Vector     y,
                                   N_Vector     fy,
                                   booleantype  jok,
                                   booleantype *jcurPtr,
                                   realtype     gamma,
                                   void *       user_data)
    {
      Assert(user_data != nullptr, ExcInternalError());
      ARKode<VectorType> &solver =
        *static_cast<ARKode<VectorType> *>(user_data);

      auto *src_y  = internal::unwrap_nvector_const<VectorType>(y);
      auto *src_fy = internal::unwrap_nvector_const<VectorType>(fy);

      return solver.jacobian_preconditioner_setup(
        t, *src_y, *src_fy, jok, *jcurPtr, gamma);
    }



    template <typename VectorType>
    int
    mass_matrix_times_vector_callback(N_Vector v,
                                      N_Vector Mv,
                                      realtype t,
                                      void *   mtimes_data)
    {
      Assert(mtimes_data != nullptr, ExcInternalError());
      ARKode<VectorType> &solver =
        *static_cast<ARKode<VectorType> *>(mtimes_data);

      auto *src_v  = internal::unwrap_nvector_const<VectorType>(v);
      auto *dst_Mv = internal::unwrap_nvector<VectorType>(Mv);

      return solver.mass_times_vector(t, *src_v, *dst_Mv);
    }



    template <typename VectorType>
    int
    mass_matrix_times_vector_setup_callback(realtype t, void *mtimes_data)
    {
      Assert(mtimes_data != nullptr, ExcInternalError());
      ARKode<VectorType> &solver =
        *static_cast<ARKode<VectorType> *>(mtimes_data);

      return solver.mass_times_setup(t);
    }



    template <typename VectorType>
    int
    solve_with_mass_matrix_callback(realtype t,
                                    N_Vector r,
                                    N_Vector z,
                                    realtype delta,
                                    int      lr,
                                    void *   user_data)
    {
      Assert(user_data != nullptr, ExcInternalError());
      ARKode<VectorType> &solver =
        *static_cast<ARKode<VectorType> *>(user_data);

      auto *src_r = internal::unwrap_nvector_const<VectorType>(r);
      auto *dst_z = internal::unwrap_nvector<VectorType>(z);

      return solver.mass_preconditioner_solve(t, *src_r, *dst_z, delta, lr);
    }



    template <typename VectorType>
    int
    mass_matrix_solver_setup_callback(realtype t, void *user_data)
    {
      Assert(user_data != nullptr, ExcInternalError());
      ARKode<VectorType> &solver =
        *static_cast<ARKode<VectorType> *>(user_data);

      return solver.mass_preconditioner_setup(t);
    }
  } // namespace


  template <typename VectorType>
  ARKode<VectorType>::ARKode(const AdditionalData &data)
    : ARKode(data, MPI_COMM_SELF)
  {}


  template <typename VectorType>
  ARKode<VectorType>::ARKode(const AdditionalData &data,
                             const MPI_Comm &      mpi_comm)
    : data(data)
    , arkode_mem(nullptr)
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    , arkode_ctx(nullptr)
#  endif
    , mpi_communicator(mpi_comm)
    , last_end_time(data.initial_time)
  {
    set_functions_to_trigger_an_assert();

#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    // SUNDIALS will always duplicate communicators if we provide them. This
    // can cause problems if SUNDIALS is configured with MPI and we pass along
    // MPI_COMM_SELF in a serial application as MPI won't be
    // initialized. Hence, work around that by just not providing a
    // communicator in that case.
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
  ARKode<VectorType>::solve_ode_incrementally(VectorType & solution,
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
  int
  ARKode<VectorType>::do_evolve_time(VectorType &  solution,
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

        // Let ARKode advance time by one period:
        double     actual_next_time;
        const auto status = ARKStepEvolve(arkode_mem,
                                          time.get_next_time(),
                                          solution_nvector,
                                          &actual_next_time,
                                          ARK_NORMAL);
        (void)status;
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
    return time.get_step_number();
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

#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
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

    arkode_mem = ARKStepCreate(
      explicit_function ? &explicit_function_callback<VectorType> : nullptr,
      implicit_function ? &implicit_function_callback<VectorType> : nullptr,
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
                solve_linearized_system
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
                              PREC_NONE,
                              0 /*krylov subvectors, 0 uses default*/,
                              arkode_ctx);
#  endif
          }
        status = ARKStepSetLinearSolver(arkode_mem, sun_linear_solver, nullptr);
        AssertARKode(status);
        status = ARKStepSetJacTimes(
          arkode_mem,
          jacobian_times_setup ?
            jacobian_times_vector_setup_callback<VectorType> :
            nullptr,
          jacobian_times_vector_callback<VectorType>);
        AssertARKode(status);
        if (jacobian_preconditioner_solve)
          {
            status = ARKStepSetPreconditioner(
              arkode_mem,
              jacobian_preconditioner_setup ?
                jacobian_solver_setup_callback<VectorType> :
                nullptr,
              solve_with_jacobian_callback<VectorType>);
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
                solve_mass
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
                              PREC_NONE,
                              0 /*krylov subvectors, 0 uses default*/,
                              arkode_ctx);
#  endif
          }
        booleantype mass_time_dependent =
          data.mass_is_time_independent ? SUNFALSE : SUNTRUE;
        status = ARKStepSetMassLinearSolver(arkode_mem,
                                            sun_mass_linear_solver,
                                            nullptr,
                                            mass_time_dependent);
        AssertARKode(status);
        status = ARKStepSetMassTimes(
          arkode_mem,
          mass_times_setup ?
            mass_matrix_times_vector_setup_callback<VectorType> :
            nullptr,
          mass_matrix_times_vector_callback<VectorType>,
          this);
        AssertARKode(status);

        if (mass_preconditioner_solve)
          {
            status = ARKStepSetMassPreconditioner(
              arkode_mem,
              mass_preconditioner_setup ?
                mass_matrix_solver_setup_callback<VectorType> :
                nullptr,
              solve_with_mass_matrix_callback<VectorType>);
            AssertARKode(status);
          }
      }
  }



  template <typename VectorType>
  void
  ARKode<VectorType>::set_functions_to_trigger_an_assert()
  {
    reinit_vector = [](VectorType &) {
      AssertThrow(false, ExcFunctionNotProvided("reinit_vector"));
    };

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
