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


#ifndef dealii_sundials_arkode_templates_h
#define dealii_sundials_arkode_templates_h

#include <deal.II/base/config.h>

#include <deal.II/sundials/arkode.h>

#ifdef DEAL_II_WITH_SUNDIALS

#  include <deal.II/base/discrete_time.h>

#  include <deal.II/sundials/arkode_exception.h>
#  include <deal.II/sundials/arkode_stepper.h>

#  include <arkode/arkode_arkstep.h>

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
    , stepper(std::make_shared<ARKStepper<VectorType>>(
        typename ARKStepper<VectorType>::AdditionalData(
          data.maximum_non_linear_iterations,
          data.implicit_function_is_linear,
          data.implicit_function_is_time_independent,
          data.mass_is_time_independent,
          data.anderson_acceleration_subspace)))
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

    const bool do_reset =
      reset_solver || stepper->get_arkode_memory() == nullptr;
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
        const int status =
          ARKStepSetStopTime(stepper->get_arkode_memory(), time.get_end_time());
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
        const auto status = ARKStepEvolve(stepper->get_arkode_memory(),
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
    const auto status =
      ARKStepGetNumSteps(stepper->get_arkode_memory(), &n_steps);
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

    stepper->reinit(
      current_time, solution, internal::InvocationContext {
        pending_exception
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
          ,
          arkode_ctx
#  endif
      });

    auto *arkode_mem = stepper->get_arkode_memory();

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
    return stepper->get_arkode_memory();
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
