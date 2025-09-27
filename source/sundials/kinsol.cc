// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
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

#include <deal.II/sundials/kinsol.h>

#ifdef DEAL_II_WITH_SUNDIALS

#  include <deal.II/base/scope_exit.h>
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
#  include <deal.II/sundials/utilities.h>

// Make sure we #include the SUNDIALS config file...
#  include <sundials/sundials_config.h>
// ...before the rest of the SUNDIALS files:
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
#    include <kinsol/kinsol_ls.h>
#  else
#    include <kinsol/kinsol_direct.h>
#  endif
#  include <kinsol/kinsol.h>
#  include <sunlinsol/sunlinsol_dense.h>
#  include <sunmatrix/sunmatrix_dense.h>

#  include <iomanip>
#  include <iostream>

DEAL_II_NAMESPACE_OPEN

namespace SUNDIALS
{
  template <typename VectorType>
  KINSOL<VectorType>::AdditionalData::AdditionalData(
    const SolutionStrategy         &strategy,
    const unsigned int              maximum_non_linear_iterations,
    const double                    function_tolerance,
    const double                    step_tolerance,
    const bool                      no_init_setup,
    const unsigned int              maximum_setup_calls,
    const double                    maximum_newton_step,
    const double                    dq_relative_error,
    const unsigned int              maximum_beta_failures,
    const unsigned int              anderson_subspace_size,
    const OrthogonalizationStrategy anderson_qr_orthogonalization)
    : strategy(strategy)
    , maximum_non_linear_iterations(maximum_non_linear_iterations)
    , function_tolerance(function_tolerance)
    , step_tolerance(step_tolerance)
    , no_init_setup(no_init_setup)
    , maximum_setup_calls(maximum_setup_calls)
    , maximum_newton_step(maximum_newton_step)
    , dq_relative_error(dq_relative_error)
    , maximum_beta_failures(maximum_beta_failures)
    , anderson_subspace_size(anderson_subspace_size)
    , anderson_qr_orthogonalization(anderson_qr_orthogonalization)
  {}



  template <typename VectorType>
  void
  KINSOL<VectorType>::AdditionalData::add_parameters(ParameterHandler &prm)
  {
    static std::string strategy_str("newton");
    prm.add_parameter("Solution strategy",
                      strategy_str,
                      "Choose among newton|linesearch|fixed_point|picard",
                      Patterns::Selection(
                        "newton|linesearch|fixed_point|picard"));
    prm.add_action("Solution strategy", [&](const std::string &value) {
      if (value == "newton")
        strategy = newton;
      else if (value == "linesearch")
        strategy = linesearch;
      else if (value == "fixed_point")
        strategy = fixed_point;
      else if (value == "picard")
        strategy = picard;
      else
        DEAL_II_ASSERT_UNREACHABLE();
    });
    prm.add_parameter("Maximum number of nonlinear iterations",
                      maximum_non_linear_iterations);
    prm.add_parameter("Function norm stopping tolerance", function_tolerance);
    prm.add_parameter("Scaled step stopping tolerance", step_tolerance);

    prm.enter_subsection("Newton parameters");
    prm.add_parameter("No initial matrix setup", no_init_setup);
    prm.add_parameter("Maximum iterations without matrix setup",
                      maximum_setup_calls);
    prm.add_parameter("Maximum allowable scaled length of the Newton step",
                      maximum_newton_step);
    prm.add_parameter("Relative error for different quotient computation",
                      dq_relative_error);
    prm.leave_subsection();

    prm.enter_subsection("Linesearch parameters");
    prm.add_parameter("Maximum number of beta-condition failures",
                      maximum_beta_failures);
    prm.leave_subsection();


    prm.enter_subsection("Fixed point and Picard parameters");
    prm.add_parameter("Anderson acceleration subspace size",
                      anderson_subspace_size);

    static std::string orthogonalization_str("modified_gram_schmidt");
    prm.add_parameter(
      "Anderson QR orthogonalization",
      orthogonalization_str,
      "Choose among modified_gram_schmidt|inverse_compact|"
      "classical_gram_schmidt|delayed_classical_gram_schmidt",
      Patterns::Selection(
        "modified_gram_schmidt|inverse_compact|classical_gram_schmidt|"
        "delayed_classical_gram_schmidt"));
    prm.add_action("Anderson QR orthogonalization",
                   [&](const std::string &value) {
                     if (value == "modified_gram_schmidt")
                       anderson_qr_orthogonalization = modified_gram_schmidt;
                     else if (value == "inverse_compact")
                       anderson_qr_orthogonalization = inverse_compact;
                     else if (value == "classical_gram_schmidt")
                       anderson_qr_orthogonalization = classical_gram_schmidt;
                     else if (value == "delayed_classical_gram_schmidt")
                       anderson_qr_orthogonalization =
                         delayed_classical_gram_schmidt;
                     else
                       DEAL_II_ASSERT_UNREACHABLE();
                   });
    prm.leave_subsection();
  }



  template <typename VectorType>
  KINSOL<VectorType>::KINSOL(const AdditionalData &data)
    : KINSOL(data, MPI_COMM_SELF)
  {}



  template <typename VectorType>
  KINSOL<VectorType>::KINSOL(const AdditionalData &data,
                             const MPI_Comm        mpi_comm)
    : data(data)
    , mpi_communicator(mpi_comm)
    , kinsol_mem(nullptr)
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    , kinsol_ctx(nullptr)
#  endif
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
                        &kinsol_ctx);
    (void)status;
    AssertKINSOL(status);
#  elif DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    const int status =
      SUNContext_Create(mpi_communicator == MPI_COMM_SELF ? nullptr :
                                                            &mpi_communicator,
                        &kinsol_ctx);
    (void)status;
    AssertKINSOL(status);
#  endif
  }



  template <typename VectorType>
  KINSOL<VectorType>::~KINSOL()
  {
    KINFree(&kinsol_mem);
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    const int status = SUNContext_Free(&kinsol_ctx);
    (void)status;
    AssertKINSOL(status);
#  endif

    Assert(pending_exception == nullptr, ExcInternalError());
  }



  template <typename VectorType>
  unsigned int
  KINSOL<VectorType>::solve(VectorType &initial_guess_and_solution)
  {
    // Make sure we have what we need
    if (data.strategy == AdditionalData::fixed_point)
      {
        Assert(iteration_function,
               ExcFunctionNotProvided("iteration_function"));
      }
    else
      {
        Assert(residual, ExcFunctionNotProvided("residual"));
        Assert(solve_with_jacobian,
               ExcFunctionNotProvided("solve_with_jacobian"));
      }

    // Create a new solver object:
    int status = 0;
    (void)status;

    KINFree(&kinsol_mem);
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    status = SUNContext_Free(&kinsol_ctx);
    AssertKINSOL(status);
#  endif


#  if DEAL_II_SUNDIALS_VERSION_GTE(7, 0, 0)
    // Same comment applies as in class constructor:
    status =
      SUNContext_Create(mpi_communicator == MPI_COMM_SELF ? SUN_COMM_NULL :
                                                            mpi_communicator,
                        &kinsol_ctx);
    AssertKINSOL(status);

    kinsol_mem = KINCreate(kinsol_ctx);
#  elif DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    // Same comment applies as in class constructor:
    status =
      SUNContext_Create(mpi_communicator == MPI_COMM_SELF ? nullptr :
                                                            &mpi_communicator,
                        &kinsol_ctx);
    AssertKINSOL(status);

    kinsol_mem = KINCreate(kinsol_ctx);
#  else
    kinsol_mem = KINCreate();
#  endif

    status = KINSetUserData(kinsol_mem, static_cast<void *>(this));
    AssertKINSOL(status);

    // helper function to create N_Vectors compatible with different versions
    // of SUNDIALS
    const auto make_compatible_nvector_view =
#  if DEAL_II_SUNDIALS_VERSION_LT(6, 0, 0)
      [](auto &v) { return internal::make_nvector_view(v); };
#  else
      [this](auto &v) { return internal::make_nvector_view(v, kinsol_ctx); };
#  endif


    VectorType ones;
    // Prepare constant vector for scaling f or u only if we need it
    if (!get_function_scaling || !get_solution_scaling)
      {
        reinit_vector(ones);
        ones = 1.0;
      }

    auto u_scale = make_compatible_nvector_view(
      get_solution_scaling ? get_solution_scaling() : ones);
    auto f_scale = make_compatible_nvector_view(
      get_function_scaling ? get_function_scaling() : ones);

    auto solution = make_compatible_nvector_view(initial_guess_and_solution);

    // This must be called before KINSetMAA
    status = KINSetNumMaxIters(kinsol_mem, data.maximum_non_linear_iterations);
    AssertKINSOL(status);

    // From the manual: this must be called BEFORE KINInit
    status = KINSetMAA(kinsol_mem, data.anderson_subspace_size);
    AssertKINSOL(status);

#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    // From the manual: this must be called BEFORE KINInit
    status = KINSetOrthAA(kinsol_mem, data.anderson_qr_orthogonalization);
    AssertKINSOL(status);
#  else
    AssertThrow(
      data.anderson_qr_orthogonalization ==
        AdditionalData::modified_gram_schmidt,
      ExcMessage(
        "You specified an orthogonalization strategy for QR factorization "
        "different from the default (modified Gram-Schmidt) but the installed "
        "SUNDIALS version does not support this feature. Either choose the "
        "default or install a SUNDIALS version >= 6.0.0."));
#  endif

    if (data.strategy == AdditionalData::fixed_point)
      status = KINInit(
        kinsol_mem,
        /* wrap up the iteration_function() callback: */
        [](N_Vector yy, N_Vector FF, void *user_data) -> int {
          KINSOL<VectorType> &solver =
            *static_cast<KINSOL<VectorType> *>(user_data);

          auto src_yy = internal::unwrap_nvector_const<VectorType>(yy);
          auto dst_FF = internal::unwrap_nvector<VectorType>(FF);

          Assert(solver.iteration_function, ExcInternalError());

          const int err = Utilities::call_and_possibly_capture_exception(
            solver.iteration_function,
            solver.pending_exception,
            *src_yy,
            *dst_FF);

          return err;
        },
        solution);
    else
      status = KINInit(
        kinsol_mem,
        /* wrap up the residual() callback: */
        [](N_Vector yy, N_Vector FF, void *user_data) -> int {
          KINSOL<VectorType> &solver =
            *static_cast<KINSOL<VectorType> *>(user_data);

          auto src_yy = internal::unwrap_nvector_const<VectorType>(yy);
          auto dst_FF = internal::unwrap_nvector<VectorType>(FF);

          Assert(solver.residual, ExcInternalError());

          const int err = Utilities::call_and_possibly_capture_exception(
            solver.residual, solver.pending_exception, *src_yy, *dst_FF);

          return err;
        },
        solution);
    AssertKINSOL(status);

    status = KINSetFuncNormTol(kinsol_mem, data.function_tolerance);
    AssertKINSOL(status);

    status = KINSetScaledStepTol(kinsol_mem, data.step_tolerance);
    AssertKINSOL(status);

    status = KINSetMaxSetupCalls(kinsol_mem, data.maximum_setup_calls);
    AssertKINSOL(status);

    status = KINSetNoInitSetup(kinsol_mem, data.no_init_setup);
    AssertKINSOL(status);

    status = KINSetMaxNewtonStep(kinsol_mem, data.maximum_newton_step);
    AssertKINSOL(status);

    status = KINSetMaxBetaFails(kinsol_mem, data.maximum_beta_failures);
    AssertKINSOL(status);

    status = KINSetRelErrFunc(kinsol_mem, data.dq_relative_error);
    AssertKINSOL(status);

    SUNMatrix       J  = nullptr;
    SUNLinearSolver LS = nullptr;

    // user assigned a function object to the solver slot
    if (solve_with_jacobian)
      {
        // Set the operations we care for in the sun_linear_solver object
        // and attach it to the KINSOL object. The functions that will get
        // called do not actually receive the KINSOL object, just the LS
        // object, so we have to store a pointer to the current
        // object in the LS object
#  if DEAL_II_SUNDIALS_VERSION_LT(6, 0, 0)
        LS = SUNLinSolNewEmpty();
#  else
        LS = SUNLinSolNewEmpty(kinsol_ctx);
#  endif
        LS->content = this;

        LS->ops->gettype =
          [](SUNLinearSolver /*ignored*/) -> SUNLinearSolver_Type {
          return SUNLINEARSOLVER_MATRIX_ITERATIVE;
        };

        LS->ops->free = [](SUNLinearSolver LS) -> int {
          if (LS->content)
            {
              LS->content = nullptr;
            }
          if (LS->ops)
            {
              std::free(LS->ops);
              LS->ops = nullptr;
            }
          std::free(LS);
          LS = nullptr;
          return 0;
        };

        LS->ops->solve = [](SUNLinearSolver LS,
                            SUNMatrix /*ignored*/,
                            N_Vector           x,
                            N_Vector           b,
                            SUNDIALS::realtype tol) -> int {
          // Receive the object that describes the linear solver and
          // unpack the pointer to the KINSOL object from which we can then
          // get the 'reinit' and 'solve' functions.
          const KINSOL<VectorType> &solver =
            *static_cast<const KINSOL<VectorType> *>(LS->content);

          Assert(solver.solve_with_jacobian, ExcInternalError());

          auto src_b = internal::unwrap_nvector_const<VectorType>(b);
          auto dst_x = internal::unwrap_nvector<VectorType>(x);

          const int err = Utilities::call_and_possibly_capture_exception(
            solver.solve_with_jacobian,
            solver.pending_exception,
            *src_b,
            *dst_x,
            tol);

          return err;
        };

        // Even though we don't use it, KINSOL still wants us to set some
        // kind of matrix object for the nonlinear solver. This is because
        // if we don't set it, it won't call the functions that set up
        // the matrix object (i.e., the argument to the 'KINSetJacFn'
        // function below).
#  if DEAL_II_SUNDIALS_VERSION_LT(6, 0, 0)
        J = SUNMatNewEmpty();
#  else
        J  = SUNMatNewEmpty(kinsol_ctx);
#  endif
        J->content = this;

        J->ops->getid = [](SUNMatrix /*ignored*/) -> SUNMatrix_ID {
          return SUNMATRIX_CUSTOM;
        };

        J->ops->destroy = [](SUNMatrix A) {
          if (A->content)
            {
              A->content = nullptr;
            }
          if (A->ops)
            {
              std::free(A->ops);
              A->ops = nullptr;
            }
          std::free(A);
          A = nullptr;
        };

        // Now set the linear system and Jacobian objects in the solver:
        status = KINSetLinearSolver(kinsol_mem, LS, J);
        AssertKINSOL(status);

        // Finally, if we were given a set-up function, tell KINSOL about
        // it as well. The manual says that this must happen *after*
        // calling KINSetLinearSolver
        if (!setup_jacobian)
          setup_jacobian = [](const VectorType &, const VectorType &) {
            return 0;
          };
        status = KINSetJacFn(
          kinsol_mem,
          [](N_Vector u,
             N_Vector f,
             SUNMatrix /* ignored */,
             void *user_data,
             N_Vector /* tmp1 */,
             N_Vector /* tmp2 */) {
            // Receive the object that describes the linear solver and
            // unpack the pointer to the KINSOL object from which we can then
            // get the 'setup' function.
            const KINSOL<VectorType> &solver =
              *static_cast<const KINSOL<VectorType> *>(user_data);

            auto ycur = internal::unwrap_nvector_const<VectorType>(u);
            auto fcur = internal::unwrap_nvector<VectorType>(f);

            // Call the user-provided setup function with these arguments:
            return Utilities::call_and_possibly_capture_exception(
              solver.setup_jacobian, solver.pending_exception, *ycur, *fcur);
          });
        AssertKINSOL(status);
      }

    // Right before calling the main KINSol() call, allow expert users to
    // perform additional setup operations on the KINSOL object.
    if (custom_setup)
      custom_setup(kinsol_mem);


    // Having set up all of the ancillary things, finally call the main KINSol
    // function. Once we return, check what happened:
    // - If we have a pending recoverable exception, ignore it if SUNDIAL's
    //   return code was zero -- in that case, SUNDIALS managed to indeed
    //   recover and we no longer need the exception
    // - If we have any other exception, rethrow it
    // - If no exception, test that SUNDIALS really did successfully return
    //
    // This all creates difficult exit paths from this function. We have to
    // do some manual clean ups to get rid of the explicitly created
    // temporary objects of this class. To avoid having to repeat the clean-up
    // code on each exit path, we package it up and put the code into a
    // ScopeExit object that is executed automatically on each such path
    // out of this function.
    Assert(pending_exception == nullptr, ExcInternalError());
    status = KINSol(kinsol_mem, solution, data.strategy, u_scale, f_scale);

    ScopeExit upon_exit([this, &J, &LS]() mutable {
      if (J != nullptr)
        SUNMatDestroy(J);
      if (LS != nullptr)
        SUNLinSolFree(LS);
      KINFree(&kinsol_mem);
    });

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
    // It is of course also possible that KINSOL experienced
    // convergence issues even if the user-side callbacks
    // succeeded. In that case, we also want to throw an exception
    // that can be caught by the user -- whether that's actually
    // useful to determine a different course of action (i.e., whether
    // the user side can do something to recover the ability to
    // converge) is a separate matter that we need not decide
    // here. (One could imagine this happening in a time or load
    // stepping procedure where re-starting with a smaller time step
    // or load step could help.)
    AssertThrow(status >= 0, ExcKINSOLError(status));

    long nniters;
    status = KINGetNumNonlinSolvIters(kinsol_mem, &nniters);
    AssertKINSOL(status);

    return static_cast<unsigned int>(nniters);
  }



  template <typename VectorType>
  void
  KINSOL<VectorType>::set_functions_to_trigger_an_assert()
  {
    reinit_vector = [](VectorType &) {
      AssertThrow(false, ExcFunctionNotProvided("reinit_vector"));
    };
  }

  template class KINSOL<Vector<double>>;
  template class KINSOL<BlockVector<double>>;

  template class KINSOL<LinearAlgebra::distributed::Vector<double>>;
  template class KINSOL<LinearAlgebra::distributed::BlockVector<double>>;

#  ifdef DEAL_II_WITH_MPI

#    ifdef DEAL_II_WITH_TRILINOS
  template class KINSOL<TrilinosWrappers::MPI::Vector>;
  template class KINSOL<TrilinosWrappers::MPI::BlockVector>;
#    endif

#    ifdef DEAL_II_WITH_PETSC
#      ifndef PETSC_USE_COMPLEX
  template class KINSOL<PETScWrappers::MPI::Vector>;
  template class KINSOL<PETScWrappers::MPI::BlockVector>;
#      endif
#    endif

#  endif

} // namespace SUNDIALS

DEAL_II_NAMESPACE_CLOSE

#endif
