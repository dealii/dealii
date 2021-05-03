//-----------------------------------------------------------
//
//    Copyright (C) 2017 - 2020 by the deal.II authors
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
//-----------------------------------------------------------


#include <deal.II/base/config.h>

#include <deal.II/sundials/kinsol.h>

#ifdef DEAL_II_WITH_SUNDIALS

#  include <deal.II/base/utilities.h>

#  include <deal.II/lac/block_vector.h>
#  ifdef DEAL_II_WITH_TRILINOS
#    include <deal.II/lac/trilinos_parallel_block_vector.h>
#    include <deal.II/lac/trilinos_vector.h>
#  endif
#  ifdef DEAL_II_WITH_PETSC
#    include <deal.II/lac/petsc_block_vector.h>
#    include <deal.II/lac/petsc_vector.h>
#  endif

#  include <deal.II/sundials/copy.h>

// Make sure we #include the SUNDIALS config file...
#  include <sundials/sundials_config.h>
// ...before the rest of the SUNDIALS files:
#  include <kinsol/kinsol_direct.h>
#  include <sunlinsol/sunlinsol_dense.h>
#  include <sunmatrix/sunmatrix_dense.h>

#  include <iomanip>
#  include <iostream>

DEAL_II_NAMESPACE_OPEN

namespace SUNDIALS
{
  using namespace internal;


  template <typename VectorType>
  KINSOL<VectorType>::AdditionalData::AdditionalData(
    const SolutionStrategy &strategy,
    const unsigned int      maximum_non_linear_iterations,
    const double            function_tolerance,
    const double            step_tolerance,
    const bool              no_init_setup,
    const unsigned int      maximum_setup_calls,
    const double            maximum_newton_step,
    const double            dq_relative_error,
    const unsigned int      maximum_beta_failures,
    const unsigned int      anderson_subspace_size)
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
        Assert(false, ExcInternalError());
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
    prm.leave_subsection();
  }


  namespace
  {
    template <typename VectorType>
    int
    residual_or_iteration_callback(N_Vector yy, N_Vector FF, void *user_data)
    {
      KINSOL<VectorType> &solver =
        *static_cast<KINSOL<VectorType> *>(user_data);
      GrowingVectorMemory<VectorType> mem;

      typename VectorMemory<VectorType>::Pointer src_yy(mem);
      solver.reinit_vector(*src_yy);

      typename VectorMemory<VectorType>::Pointer dst_FF(mem);
      solver.reinit_vector(*dst_FF);

      copy(*src_yy, yy);

      int err = 0;
      if (solver.residual)
        err = solver.residual(*src_yy, *dst_FF);
      else if (solver.iteration_function)
        err = solver.iteration_function(*src_yy, *dst_FF);
      else
        Assert(false, ExcInternalError());

      copy(FF, *dst_FF);

      return err;
    }



#  if DEAL_II_SUNDIALS_VERSION_LT(4, 1, 0)
    template <typename VectorType>
    int
    setup_jacobian_callback(KINMem kinsol_mem)
    {
      KINSOL<VectorType> &solver =
        *static_cast<KINSOL<VectorType> *>(kinsol_mem->kin_user_data);
      GrowingVectorMemory<VectorType> mem;

      typename VectorMemory<VectorType>::Pointer src_ycur(mem);
      solver.reinit_vector(*src_ycur);

      typename VectorMemory<VectorType>::Pointer src_fcur(mem);
      solver.reinit_vector(*src_fcur);

      copy(*src_ycur, kinsol_mem->kin_uu);
      copy(*src_fcur, kinsol_mem->kin_fval);

      int err = solver.setup_jacobian(*src_ycur, *src_fcur);
      return err;
    }



    template <typename VectorType>
    int
    solve_with_jacobian_callback(KINMem    kinsol_mem,
                                 N_Vector  x,
                                 N_Vector  b,
                                 realtype *sJpnorm,
                                 realtype *sFdotJp)
    {
      KINSOL<VectorType> &solver =
        *static_cast<KINSOL<VectorType> *>(kinsol_mem->kin_user_data);
      GrowingVectorMemory<VectorType> mem;

      typename VectorMemory<VectorType>::Pointer src_ycur(mem);
      solver.reinit_vector(*src_ycur);

      typename VectorMemory<VectorType>::Pointer src_fcur(mem);
      solver.reinit_vector(*src_fcur);

      copy(*src_ycur, kinsol_mem->kin_uu);
      copy(*src_fcur, kinsol_mem->kin_fval);

      typename VectorMemory<VectorType>::Pointer src(mem);
      solver.reinit_vector(*src);

      typename VectorMemory<VectorType>::Pointer dst(mem);
      solver.reinit_vector(*dst);

      copy(*src, b);

      int err = solver.solve_jacobian_system(*src_ycur, *src_fcur, *src, *dst);
      copy(x, *dst);

      *sJpnorm = N_VWL2Norm(b, kinsol_mem->kin_fscale);
      N_VProd(b, kinsol_mem->kin_fscale, b);
      N_VProd(b, kinsol_mem->kin_fscale, b);
      *sFdotJp = N_VDotProd(kinsol_mem->kin_fval, b);

      return err;
    }

#  else // SUNDIALS 5.0 or later

    template <typename VectorType>
    int
    setup_jacobian_callback(N_Vector u,
                            N_Vector f,
                            SUNMatrix /* ignored */,
                            void *user_data,
                            N_Vector /* tmp1 */,
                            N_Vector /* tmp2 */)
    {
      // Receive the object that describes the linear solver and
      // unpack the pointer to the KINSOL object from which we can then
      // get the 'setup' function.
      const KINSOL<VectorType> &solver =
        *static_cast<const KINSOL<VectorType> *>(user_data);

      // Allocate temporary (deal.II-type) vectors into which to copy the
      // N_vectors
      GrowingVectorMemory<VectorType>            mem;
      typename VectorMemory<VectorType>::Pointer ycur(mem);
      typename VectorMemory<VectorType>::Pointer fcur(mem);
      solver.reinit_vector(*ycur);
      solver.reinit_vector(*fcur);

      copy(*ycur, u);
      copy(*fcur, f);

      // Call the user-provided setup function with these arguments:
      solver.setup_jacobian(*ycur, *fcur);

      return 0;
    }



    template <typename VectorType>
    int
    solve_with_jacobian_callback(SUNLinearSolver LS,
                                 SUNMatrix /*ignored*/,
                                 N_Vector x,
                                 N_Vector b,
                                 realtype tol)
    {
      // Receive the object that describes the linear solver and
      // unpack the pointer to the KINSOL object from which we can then
      // get the 'reinit' and 'solve' functions.
      const KINSOL<VectorType> &solver =
        *static_cast<const KINSOL<VectorType> *>(LS->content);

      // This is where we have to make a decision about which of the two
      // signals to call. Let's first check the more modern one:
      if (solver.solve_with_jacobian)
        {
          // Allocate temporary (deal.II-type) vectors into which to copy the
          // N_vectors
          GrowingVectorMemory<VectorType>            mem;
          typename VectorMemory<VectorType>::Pointer src_b(mem);
          typename VectorMemory<VectorType>::Pointer dst_x(mem);

          solver.reinit_vector(*src_b);
          solver.reinit_vector(*dst_x);

          copy(*src_b, b);

          const int err = solver.solve_with_jacobian(*src_b, *dst_x, tol);

          copy(x, *dst_x);

          return err;
        }
      else
        {
          // User has not provided the modern callback, so the fact that we are
          // here means that they must have given us something for the old
          // signal. Check this.
          Assert(solver.solve_jacobian_system, ExcInternalError());

          // Allocate temporary (deal.II-type) vectors into which to copy the
          // N_vectors
          GrowingVectorMemory<VectorType>            mem;
          typename VectorMemory<VectorType>::Pointer src_ycur(mem);
          typename VectorMemory<VectorType>::Pointer src_fcur(mem);
          typename VectorMemory<VectorType>::Pointer src_b(mem);
          typename VectorMemory<VectorType>::Pointer dst_x(mem);

          solver.reinit_vector(*src_b);
          solver.reinit_vector(*dst_x);

          copy(*src_b, b);

          // Call the user-provided setup function with these arguments. Note
          // that Sundials 4.x and later no longer provide values for
          // src_ycur and src_fcur, and so we simply pass dummy vector in.
          // These vectors will have zero lengths because we don't reinit them
          // above.
          const int err =
            solver.solve_jacobian_system(*src_ycur, *src_fcur, *src_b, *dst_x);

          copy(x, *dst_x);

          return err;
        }
    }

#  endif
  } // namespace



  template <typename VectorType>
  KINSOL<VectorType>::KINSOL(const AdditionalData &data,
                             const MPI_Comm &      mpi_comm)
    : data(data)
    , kinsol_mem(nullptr)
    , solution(nullptr)
    , u_scale(nullptr)
    , f_scale(nullptr)
    , communicator(is_serial_vector<VectorType>::value ?
                     MPI_COMM_SELF :
                     Utilities::MPI::duplicate_communicator(mpi_comm))
  {
    set_functions_to_trigger_an_assert();
  }



  template <typename VectorType>
  KINSOL<VectorType>::~KINSOL()
  {
    if (kinsol_mem)
      KINFree(&kinsol_mem);

#  ifdef DEAL_II_WITH_MPI
    if (is_serial_vector<VectorType>::value == false)
      {
        const int ierr = MPI_Comm_free(&communicator);
        (void)ierr;
        AssertNothrow(ierr == MPI_SUCCESS, ExcMPI(ierr));
      }
#  endif
  }



  template <typename VectorType>
  unsigned int
  KINSOL<VectorType>::solve(VectorType &initial_guess_and_solution)
  {
    unsigned int system_size = initial_guess_and_solution.size();

    // The solution is stored in
    // solution. Here we take only a
    // view of it.
#  ifdef DEAL_II_WITH_MPI
    if (is_serial_vector<VectorType>::value == false)
      {
        const IndexSet is = initial_guess_and_solution.locally_owned_elements();
        const unsigned int local_system_size = is.n_elements();

        solution =
          N_VNew_Parallel(communicator, local_system_size, system_size);

        u_scale = N_VNew_Parallel(communicator, local_system_size, system_size);
        N_VConst_Parallel(1.e0, u_scale);

        f_scale = N_VNew_Parallel(communicator, local_system_size, system_size);
        N_VConst_Parallel(1.e0, f_scale);
      }
    else
#  endif
      {
        Assert(is_serial_vector<VectorType>::value,
               ExcInternalError(
                 "Trying to use a serial code with a parallel vector."));
        solution = N_VNew_Serial(system_size);
        u_scale  = N_VNew_Serial(system_size);
        N_VConst_Serial(1.e0, u_scale);
        f_scale = N_VNew_Serial(system_size);
        N_VConst_Serial(1.e0, f_scale);
      }

    if (get_solution_scaling)
      copy(u_scale, get_solution_scaling());

    if (get_function_scaling)
      copy(f_scale, get_function_scaling());

    copy(solution, initial_guess_and_solution);

    if (kinsol_mem)
      KINFree(&kinsol_mem);

    kinsol_mem = KINCreate();

    int status =
      KINInit(kinsol_mem, residual_or_iteration_callback<VectorType>, solution);
    (void)status;
    AssertKINSOL(status);

    status = KINSetUserData(kinsol_mem, static_cast<void *>(this));
    AssertKINSOL(status);

    status = KINSetNumMaxIters(kinsol_mem, data.maximum_non_linear_iterations);
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

    status = KINSetMAA(kinsol_mem, data.anderson_subspace_size);
    AssertKINSOL(status);

    status = KINSetRelErrFunc(kinsol_mem, data.dq_relative_error);
    AssertKINSOL(status);

    SUNMatrix       J  = nullptr;
    SUNLinearSolver LS = nullptr;

    if (solve_jacobian_system ||
        solve_with_jacobian) // user assigned a function object to the solver
                             // slot
      {
/* interface up to and including 4.0 */
#  if DEAL_II_SUNDIALS_VERSION_LT(4, 1, 0)
        auto KIN_mem        = static_cast<KINMem>(kinsol_mem);
        KIN_mem->kin_lsolve = solve_with_jacobian_callback<VectorType>;
        if (setup_jacobian) // user assigned a function object to the Jacobian
          // set-up slot
          KIN_mem->kin_lsetup = setup_jacobian_callback<VectorType>;

/* interface up to and including 4.1 */
#  elif DEAL_II_SUNDIALS_VERSION_LT(5, 0, 0)

        // deal.II does not currently have support for KINSOL in
        // SUNDIALS 4.1. One could write this and update this section,
        // but it does not seem worthwhile spending the time to
        // interface with an old version of SUNDIAL given that the
        // code below supports modern SUNDIAL versions just fine.
        Assert(false, ExcNotImplemented());

#  else /* interface starting with SUNDIALS 5.0 */
        // Set the operations we care for in the sun_linear_solver object
        // and attach it to the KINSOL object. The functions that will get
        // called do not actually receive the KINSOL object, just the LS
        // object, so we have to store a pointer to the current
        // object in the LS object
        LS          = SUNLinSolNewEmpty();
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
              free(LS->ops);
              LS->ops = nullptr;
            }
          free(LS);
          LS = nullptr;
          return 0;
        };

        LS->ops->solve = solve_with_jacobian_callback<VectorType>;

        // Even though we don't use it, KINSOL still wants us to set some
        // kind of matrix object for the nonlinear solver. This is because
        // if we don't set it, it won't call the functions that set up
        // the matrix object (i.e., the argument to the 'KINSetJacFn'
        // function below).
        J          = SUNMatNewEmpty();
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
              free(A->ops);
              A->ops = nullptr;
            }
          free(A);
          A = nullptr;
        };

        // Now set the linear system and Jacobian objects in the solver:
        status = KINSetLinearSolver(kinsol_mem, LS, J);
        AssertKINSOL(status);

        // Finally, if we were given a set-up function, tell KINSOL about
        // it as well. The manual says that this must happen *after*
        // calling KINSetLinearSolver
        if (setup_jacobian)
          {
            status =
              KINSetJacFn(kinsol_mem, &setup_jacobian_callback<VectorType>);
            AssertKINSOL(status);
          }
#  endif
      }
    else
      {
        J      = SUNDenseMatrix(system_size, system_size);
        LS     = SUNDenseLinearSolver(u_scale, J);
        status = KINDlsSetLinearSolver(kinsol_mem, LS, J);
        AssertKINSOL(status);
      }

    if (data.strategy == AdditionalData::newton ||
        data.strategy == AdditionalData::linesearch)
      Assert(residual, ExcFunctionNotProvided("residual"));

    if (data.strategy == AdditionalData::fixed_point ||
        data.strategy == AdditionalData::picard)
      Assert(iteration_function, ExcFunctionNotProvided("iteration_function"));

    // call to KINSol
    status = KINSol(kinsol_mem, solution, data.strategy, u_scale, f_scale);
    AssertKINSOL(status);

    copy(initial_guess_and_solution, solution);

    // Free the vectors which are no longer used.
#  ifdef DEAL_II_WITH_MPI
    if (is_serial_vector<VectorType>::value == false)
      {
        N_VDestroy_Parallel(solution);
        N_VDestroy_Parallel(u_scale);
        N_VDestroy_Parallel(f_scale);
      }
    else
#  endif
      {
        N_VDestroy_Serial(solution);
        N_VDestroy_Serial(u_scale);
        N_VDestroy_Serial(f_scale);
      }

    long nniters;
    status = KINGetNumNonlinSolvIters(kinsol_mem, &nniters);
    AssertKINSOL(status);

    SUNMatDestroy(J);
    SUNLinSolFree(LS);
    KINFree(&kinsol_mem);

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
