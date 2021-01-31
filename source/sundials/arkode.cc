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

#include <deal.II/sundials/arkode.h>

#ifdef DEAL_II_WITH_SUNDIALS

#  include <deal.II/base/discrete_time.h>
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

#  include <deal.II/sundials/n_vector.h>

#  if DEAL_II_SUNDIALS_VERSION_LT(4, 0, 0)
#    include <arkode/arkode_impl.h>
#    include <sundials/sundials_config.h>
#  else
#    include <arkode/arkode_arkstep.h>
#    include <sunlinsol/sunlinsol_spgmr.h>
#    include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#    if DEAL_II_SUNDIALS_VERSION_LT(5, 0, 0)
#      include <deal.II/sundials/sunlinsol_newempty.h>
#    endif
#  endif

#  include <iostream>

#  if DEAL_II_SUNDIALS_VERSION_LT(4, 0, 0)
// Make sure we know how to call sundials own ARKode() function
const auto &SundialsARKode = ARKode;
#  endif

DEAL_II_NAMESPACE_OPEN

namespace SUNDIALS
{
  using namespace internal;

  namespace
  {
    template <typename VectorType>
    int
    t_arkode_explicit_function(realtype tt,
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
    t_arkode_implicit_function(realtype tt,
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



#  if DEAL_II_SUNDIALS_VERSION_LT(4, 0, 0)
    template <typename VectorType>
    int
    t_arkode_setup_jacobian(ARKodeMem    arkode_mem,
                            int          convfail,
                            N_Vector     ypred,
                            N_Vector     fpred,
                            booleantype *jcurPtr,
                            N_Vector,
                            N_Vector,
                            N_Vector)
    {
      Assert(arkode_mem->ark_user_data != nullptr, ExcInternalError());
      ARKode<VectorType> &solver =
        *static_cast<ARKode<VectorType> *>(arkode_mem->ark_user_data);

      auto *src_ypred = internal::unwrap_nvector_const<VectorType>(ypred);
      auto *src_fpred = internal::unwrap_nvector_const<VectorType>(fpred);
      // avoid reinterpret_cast
      bool jcurPtr_tmp = false;
      int  err         = solver.setup_jacobian(convfail,
                                      arkode_mem->ark_tn,
                                      arkode_mem->ark_gamma,
                                      *src_ypred,
                                      *src_fpred,
                                      jcurPtr_tmp);
#    if DEAL_II_SUNDIALS_VERSION_GTE(2, 0, 0)
      *jcurPtr = jcurPtr_tmp ? SUNTRUE : SUNFALSE;
#    else
      *jcurPtr = jcurPtr_tmp ? TRUE : FALSE;
#    endif

      return err;
    }



    template <typename VectorType>
    int
    t_arkode_solve_jacobian(ARKodeMem arkode_mem,
                            N_Vector  b,
#    if DEAL_II_SUNDIALS_VERSION_LT(3, 0, 0)
                            N_Vector,
#    endif
                            N_Vector ycur,
                            N_Vector fcur)
    {
      Assert(arkode_mem->ark_user_data != nullptr, ExcInternalError());
      ARKode<VectorType> &solver =
        *static_cast<ARKode<VectorType> *>(arkode_mem->ark_user_data);

      auto *dst      = internal::unwrap_nvector<VectorType>(b);
      auto *src_ycur = internal::unwrap_nvector_const<VectorType>(ycur);
      auto *src_fcur = internal::unwrap_nvector_const<VectorType>(fcur);

      // make a temporary copy to work on in the user call
      VectorType src = *dst;

      int err = solver.solve_jacobian_system(arkode_mem->ark_tn,
                                             arkode_mem->ark_gamma,
                                             *src_ycur,
                                             *src_fcur,
                                             src,
                                             *dst);


      return err;
    }



    template <typename VectorType>
    int
    t_arkode_setup_mass(ARKodeMem arkode_mem, N_Vector, N_Vector, N_Vector)
    {
      Assert(arkode_mem->ark_user_data != nullptr, ExcInternalError());
      ARKode<VectorType> &solver =
        *static_cast<ARKode<VectorType> *>(arkode_mem->ark_user_data);
      int err = solver.setup_mass(arkode_mem->ark_tn);
      return err;
    }



    template <typename VectorType>
    int
    t_arkode_solve_mass(ARKodeMem arkode_mem,
#    if DEAL_II_SUNDIALS_VERSION_LT(3, 0, 0)
                        N_Vector b,
                        N_Vector
#    else
                        N_Vector b
#    endif
    )
    {
      Assert(arkode_mem->ark_user_data != nullptr, ExcInternalError());
      ARKode<VectorType> &solver =
        *static_cast<ARKode<VectorType> *>(arkode_mem->ark_user_data);

      auto *dst = internal::unwrap_nvector<VectorType>(b);

      // make a temporary copy to work on in the user call
      VectorType src = *dst;

      int err = solver.solve_mass_system(src, *dst);

      return err;
    }
#  else

    template <typename VectorType>
    int
    t_arkode_jac_times_vec_function(N_Vector v,
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
      GrowingVectorMemory<VectorType> mem;

      auto *src_v  = internal::unwrap_nvector_const<VectorType>(v);
      auto *src_y  = internal::unwrap_nvector_const<VectorType>(y);
      auto *src_fy = internal::unwrap_nvector_const<VectorType>(fy);

      auto *dst_Jv = internal::unwrap_nvector<VectorType>(Jv);

      return solver.jacobian_times_vector(*src_v, *dst_Jv, t, *src_y, *src_fy);
    }



    template <typename VectorType>
    int
    t_arkode_jac_times_setup_function(realtype t,
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
    t_arkode_prec_solve_function(realtype t,
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
    t_arkode_prec_setup_function(realtype     t,
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
    t_arkode_mass_times_vec_function(N_Vector v,
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
    t_arkode_mass_times_setup_function(realtype t, void *mtimes_data)
    {
      Assert(mtimes_data != nullptr, ExcInternalError());
      ARKode<VectorType> &solver =
        *static_cast<ARKode<VectorType> *>(mtimes_data);

      return solver.mass_times_setup(t);
    }



    template <typename VectorType>
    int
    t_arkode_mass_prec_solve_function(realtype t,
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
    t_arkode_mass_prec_setup_function(realtype t, void *user_data)
    {
      Assert(user_data != nullptr, ExcInternalError());
      ARKode<VectorType> &solver =
        *static_cast<ARKode<VectorType> *>(user_data);

      return solver.mass_preconditioner_setup(t);
    }



    /**
     * storage for internal content of the linear solver wrapper
     */
    template <typename VectorType>
    struct LinearSolverContent
    {
      ATimesFn a_times_fn;
      PSetupFn preconditioner_setup;
      PSolveFn preconditioner_solve;

      LinearSolveFunction<VectorType> lsolve;

      //! solver reference to forward all solver related calls to user-specified
      //! functions
      ARKode<VectorType> *solver;
      void *              P_data;
      void *              A_data;
    };



    /**
     * Access our LinearSolverContent from the generic content of the
     * SUNLinearSolver @p ls.
     */
    template <typename VectorType>
    LinearSolverContent<VectorType> *
    access_content(SUNLinearSolver ls)
    {
      Assert(ls->content != nullptr, ExcInternalError());
      return static_cast<LinearSolverContent<VectorType> *>(ls->content);
    }



    SUNLinearSolver_Type arkode_linsol_get_type(SUNLinearSolver)
    {
      return SUNLINEARSOLVER_ITERATIVE;
    }



    template <typename VectorType>
    int
    arkode_linsol_solve(SUNLinearSolver LS,
                        SUNMatrix,
                        N_Vector x,
                        N_Vector b,
                        realtype tol)
    {
      auto                content = access_content<VectorType>(LS);
      ARKode<VectorType> &solver  = *(content->solver);

      auto *src_b = internal::unwrap_nvector_const<VectorType>(b);
      auto *dst_x = internal::unwrap_nvector<VectorType>(x);

      SundialsOperator<VectorType> op(solver,
                                      content->A_data,
                                      content->a_times_fn);

      SundialsPreconditioner<VectorType> preconditioner(
        solver, content->P_data, content->preconditioner_solve, tol);

      return content->lsolve(op, preconditioner, *dst_x, *src_b, tol);
    }



    template <typename VectorType>
    int
    arkode_linsol_setup(SUNLinearSolver LS, SUNMatrix)
    {
      auto content = access_content<VectorType>(LS);
      if (content->preconditioner_setup)
        return content->preconditioner_setup(content->P_data);
      return 0;
    }



    template <typename VectorType>
    int arkode_linsol_initialize(SUNLinearSolver)
    {
      // this method is currently only provided because SUNDIALS 4.0.0 requires
      // it - no user-set action is implemented so far
      return 0;
    }



    template <typename VectorType>
    int
    arkode_linsol_set_a_times(SUNLinearSolver LS, void *A_data, ATimesFn ATimes)
    {
      auto content        = access_content<VectorType>(LS);
      content->A_data     = A_data;
      content->a_times_fn = ATimes;
      return 0;
    }



    template <typename VectorType>
    int
    arkode_linsol_set_preconditioner(SUNLinearSolver LS,
                                     void *          P_data,
                                     PSetupFn        p_setup,
                                     PSolveFn        p_solve)
    {
      auto content                  = access_content<VectorType>(LS);
      content->P_data               = P_data;
      content->preconditioner_setup = p_setup;
      content->preconditioner_solve = p_solve;
      return 0;
    }
#  endif

  } // namespace

#  if DEAL_II_SUNDIALS_VERSION_GTE(4, 0, 0)

  /*!
   * Attach wrapper functions to SUNDIALS' linear solver interface. We pretend
   * that the user-supplied linear solver is matrix-free, even though it can
   * be matrix-based. This way SUNDIALS does not need to understand our matrix
   * types.
   */
  template <typename VectorType>
  class SundialsLinearSolverWrapper
  {
  public:
    SundialsLinearSolverWrapper(ARKode<VectorType> &            solver,
                                LinearSolveFunction<VectorType> lsolve)
    {
      sun_linear_solver                  = SUNLinSolNewEmpty();
      sun_linear_solver->ops->gettype    = arkode_linsol_get_type;
      sun_linear_solver->ops->solve      = arkode_linsol_solve<VectorType>;
      sun_linear_solver->ops->setup      = arkode_linsol_setup<VectorType>;
      sun_linear_solver->ops->initialize = arkode_linsol_initialize<VectorType>;
      sun_linear_solver->ops->setatimes = arkode_linsol_set_a_times<VectorType>;
      sun_linear_solver->ops->setpreconditioner =
        arkode_linsol_set_preconditioner<VectorType>;

      content.solver             = &solver;
      content.lsolve             = lsolve;
      sun_linear_solver->content = &content;
    }

    ~SundialsLinearSolverWrapper()
    {
      SUNLinSolFreeEmpty(sun_linear_solver);
    }

    SUNLinearSolver
    get_wrapped_solver()
    {
      return sun_linear_solver;
    }

  private:
    SUNLinearSolver                 sun_linear_solver;
    LinearSolverContent<VectorType> content;
  };

#  endif

  template <typename VectorType>
  ARKode<VectorType>::ARKode(const AdditionalData &data,
                             const MPI_Comm &      mpi_comm)
    : data(data)
    , arkode_mem(nullptr)
    , communicator(is_serial_vector<VectorType>::value ?
                     MPI_COMM_SELF :
                     Utilities::MPI::duplicate_communicator(mpi_comm))
  {
    set_functions_to_trigger_an_assert();
  }

  template <typename VectorType>
  ARKode<VectorType>::~ARKode()
  {
    if (arkode_mem)
#  if DEAL_II_SUNDIALS_VERSION_LT(4, 0, 0)
      ARKodeFree(&arkode_mem);
#  else
      ARKStepFree(&arkode_mem);
#  endif

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
  ARKode<VectorType>::solve_ode(VectorType &solution)
  {
    DiscreteTime time(data.initial_time,
                      data.final_time,
                      data.initial_step_size);
    reset(time.get_current_time(), time.get_next_step_size(), solution);

    if (output_step)
      output_step(time.get_current_time(), solution, time.get_step_number());

    while (time.is_at_end() == false)
      {
        time.set_desired_next_step_size(data.output_period);
        double actual_next_time;
#  if DEAL_II_SUNDIALS_VERSION_LT(4, 0, 0)
        const auto status = SundialsARKode(
          arkode_mem, time.get_next_time(), yy, &actual_next_time, ARK_NORMAL);
#  else
        const auto status = ARKStepEvolve(
          arkode_mem, time.get_next_time(), yy, &actual_next_time, ARK_NORMAL);
#  endif
        (void)status;
        AssertARKode(status);

        time.set_next_step_size(actual_next_time - time.get_current_time());
        time.advance_time();

        while (solver_should_restart(time.get_current_time(), solution))
          reset(time.get_current_time(),
                time.get_previous_step_size(),
                solution);

        if (output_step)
          output_step(time.get_current_time(),
                      solution,
                      time.get_step_number());
      }

    return time.get_step_number();
  }

#  if DEAL_II_SUNDIALS_VERSION_LT(4, 0, 0)
  template <typename VectorType>
  void
  ARKode<VectorType>::reset(const double current_time,
                            const double current_time_step,
                            VectorType & solution)
  {
    if (arkode_mem)
      ARKodeFree(&arkode_mem);

    arkode_mem = ARKodeCreate();

    int status;
    (void)status;

    Assert(explicit_function || implicit_function,
           ExcFunctionNotProvided("explicit_function || implicit_function"));

    // just a view on the memory in solution, all write operations on yy by
    // ARKODE will automatically be mirrored to solution
    yy = internal::make_nvector_view(solution);

    status = ARKodeInit(
      arkode_mem,
      explicit_function ? &t_arkode_explicit_function<VectorType> : nullptr,
      implicit_function ? &t_arkode_implicit_function<VectorType> : nullptr,
      current_time,
      yy);
    AssertARKode(status);

    if (get_local_tolerances)
      {
        abs_tolls = make_nvector_view(get_local_tolerances());
        status =
          ARKodeSVtolerances(arkode_mem, data.relative_tolerance, abs_tolls);
        AssertARKode(status);
      }
    else
      {
        status = ARKodeSStolerances(arkode_mem,
                                    data.relative_tolerance,
                                    data.absolute_tolerance);
        AssertARKode(status);
      }

    status = ARKodeSetInitStep(arkode_mem, current_time_step);
    AssertARKode(status);

    status = ARKodeSetUserData(arkode_mem, this);
    AssertARKode(status);

    status = ARKodeSetStopTime(arkode_mem, data.final_time);
    AssertARKode(status);

    status =
      ARKodeSetMaxNonlinIters(arkode_mem, data.maximum_non_linear_iterations);
    AssertARKode(status);

    // Initialize solver
    auto ARKode_mem = static_cast<ARKodeMem>(arkode_mem);

    if (solve_jacobian_system)
      {
        status = ARKodeSetNewton(arkode_mem);
        AssertARKode(status);
        if (data.implicit_function_is_linear)
          {
            status = ARKodeSetLinear(
              arkode_mem, data.implicit_function_is_time_independent ? 0 : 1);
            AssertARKode(status);
          }


        ARKode_mem->ark_lsolve = t_arkode_solve_jacobian<VectorType>;
        if (setup_jacobian)
          {
            ARKode_mem->ark_lsetup = t_arkode_setup_jacobian<VectorType>;
#    if DEAL_II_SUNDIALS_VERSION_LT(3, 0, 0)
            ARKode_mem->ark_setupNonNull = true;
#    endif
          }
      }
    else
      {
        status =
          ARKodeSetFixedPoint(arkode_mem, data.maximum_non_linear_iterations);
        AssertARKode(status);
      }


    if (solve_mass_system)
      {
        ARKode_mem->ark_msolve = t_arkode_solve_mass<VectorType>;

        if (setup_mass)
          {
            ARKode_mem->ark_msetup = t_arkode_setup_mass<VectorType>;
#    if DEAL_II_SUNDIALS_VERSION_LT(3, 0, 0)
            ARKode_mem->ark_MassSetupNonNull = true;
#    endif
          }
      }

    status = ARKodeSetOrder(arkode_mem, data.maximum_order);
    AssertARKode(status);

    if (custom_setup)
      custom_setup(arkode_mem);
  }

#  else

  template <typename VectorType>
  void
  ARKode<VectorType>::reset(const double current_time,
                            const double current_time_step,
                            VectorType &solution)
  {
    if (arkode_mem)
      ARKStepFree(&arkode_mem);

    int status;
    (void)status;

    // just a view on the memory in solution, all write operations on yy by
    // ARKODE will automatically be mirrored to solution
    yy = internal::make_nvector_view(solution);

    Assert(explicit_function || implicit_function,
           ExcFunctionNotProvided("explicit_function || implicit_function"));

    arkode_mem = ARKStepCreate(
      explicit_function ? &t_arkode_explicit_function<VectorType> : nullptr,
      implicit_function ? &t_arkode_implicit_function<VectorType> : nullptr,
      current_time,
      yy);

    Assert(arkode_mem != nullptr, ExcInternalError());

    if (get_local_tolerances)
      {
        abs_tolls = internal::make_nvector_view(get_local_tolerances());
        status =
          ARKStepSVtolerances(arkode_mem, data.relative_tolerance, abs_tolls);
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

    setup_mass_solver();

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
              std::make_unique<SundialsLinearSolverWrapper<VectorType>>(
                *this, solve_linearized_system);
            sun_linear_solver = linear_solver->get_wrapped_solver();
          }
        else
          {
            // use default solver from SUNDIALS
            // TODO give user options
            sun_linear_solver =
              SUNLinSol_SPGMR(yy,
                              PREC_NONE,
                              0 /*krylov subvectors, 0 uses default*/);
          }
        status = ARKStepSetLinearSolver(arkode_mem, sun_linear_solver, nullptr);
        AssertARKode(status);
        status =
          ARKStepSetJacTimes(arkode_mem,
                             jacobian_times_setup ?
                               t_arkode_jac_times_setup_function<VectorType> :
                               nullptr,
                             t_arkode_jac_times_vec_function<VectorType>);
        AssertARKode(status);
        if (jacobian_preconditioner_solve)
          {
            status = ARKStepSetPreconditioner(
              arkode_mem,
              jacobian_preconditioner_setup ?
                t_arkode_prec_setup_function<VectorType> :
                nullptr,
              t_arkode_prec_solve_function<VectorType>);
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
        auto y_template = internal::make_nvector_view(solution);

        SUNNonlinearSolver fixed_point_solver =
          SUNNonlinSol_FixedPoint(y_template,
                                  data.anderson_acceleration_subspace);

        status = ARKStepSetNonlinearSolver(arkode_mem, fixed_point_solver);
        AssertARKode(status);
      }

    status =
      ARKStepSetMaxNonlinIters(arkode_mem, data.maximum_non_linear_iterations);
    AssertARKode(status);
  }



  template <typename VectorType>
  void
  ARKode<VectorType>::setup_mass_solver()
  {
    int status;
    (void)status;

    if (mass_times_vector)
      {
        SUNLinearSolver sun_mass_linear_solver;
        if (solve_mass)
          {
            mass_solver =
              std::make_unique<SundialsLinearSolverWrapper<VectorType>>(
                *this, solve_mass);
            sun_mass_linear_solver = mass_solver->get_wrapped_solver();
          }
        else
          {
            sun_mass_linear_solver =
              SUNLinSol_SPGMR(yy,
                              PREC_NONE,
                              0 /*krylov subvectors, 0 uses default*/);
          }
        booleantype mass_time_dependent =
          data.mass_is_time_independent ? SUNFALSE : SUNTRUE;
        status = ARKStepSetMassLinearSolver(arkode_mem,
                                            sun_mass_linear_solver,
                                            nullptr,
                                            mass_time_dependent);
        AssertARKode(status);
        status =
          ARKStepSetMassTimes(arkode_mem,
                              mass_times_setup ?
                                t_arkode_mass_times_setup_function<VectorType> :
                                nullptr,
                              t_arkode_mass_times_vec_function<VectorType>,
                              this);
        AssertARKode(status);

        if (mass_preconditioner_solve)
          {
            status = ARKStepSetMassPreconditioner(
              arkode_mem,
              mass_preconditioner_setup ?
                t_arkode_mass_prec_setup_function<VectorType> :
                nullptr,
              t_arkode_mass_prec_solve_function<VectorType>);
            AssertARKode(status);
          }
      }
  }
#  endif



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

#  if DEAL_II_SUNDIALS_VERSION_GTE(4, 0, 0)

  template <typename VectorType>
  SundialsOperator<VectorType>::SundialsOperator(ARKode<VectorType> &solver,
                                                 void *              A_data,
                                                 ATimesFn            a_times_fn)
    : solver(solver)
    , A_data(A_data)
    , a_times_fn(a_times_fn)

  {
    Assert(a_times_fn != nullptr, ExcInternalError());
  }



  template <typename VectorType>
  void
  SundialsOperator<VectorType>::vmult(VectorType &      dst,
                                      const VectorType &src) const
  {
    auto sun_dst = internal::make_nvector_view(dst);
    auto sun_src = internal::make_nvector_view(src);
    int  status  = a_times_fn(A_data, sun_src, sun_dst);
    (void)status;
    AssertARKode(status);
  }



  template <typename VectorType>
  SundialsPreconditioner<VectorType>::SundialsPreconditioner(
    ARKode<VectorType> &solver,
    void *              P_data,
    PSolveFn            p_solve_fn,
    double              tol)
    : solver(solver)
    , P_data(P_data)
    , p_solve_fn(p_solve_fn)
    , tol(tol)
  {}



  template <typename VectorType>
  void
  SundialsPreconditioner<VectorType>::vmult(VectorType &      dst,
                                            const VectorType &src) const
  {
    // apply identity preconditioner if nothing else specified
    if (!p_solve_fn)
      {
        dst = src;
        return;
      }

    auto sun_dst = internal::make_nvector_view(dst);
    auto sun_src = internal::make_nvector_view(src);
    // for custom preconditioners no distinction between left and right
    // preconditioning is made
    int status =
      p_solve_fn(P_data, sun_src, sun_dst, tol, 0 /*precondition_type*/);
    (void)status;
    AssertARKode(status);
  }
#  endif

  template class ARKode<Vector<double>>;
  template class ARKode<BlockVector<double>>;

#  if DEAL_II_SUNDIALS_VERSION_GTE(4, 0, 0)
  template struct SundialsOperator<Vector<double>>;
  template struct SundialsOperator<BlockVector<double>>;

  template struct SundialsPreconditioner<Vector<double>>;
  template struct SundialsPreconditioner<BlockVector<double>>;
#  endif

#  ifdef DEAL_II_WITH_MPI

#    ifdef DEAL_II_WITH_TRILINOS
  template class ARKode<TrilinosWrappers::MPI::Vector>;
  template class ARKode<TrilinosWrappers::MPI::BlockVector>;

#      if DEAL_II_SUNDIALS_VERSION_GTE(4, 0, 0)
  template struct SundialsOperator<TrilinosWrappers::MPI::Vector>;
  template struct SundialsOperator<TrilinosWrappers::MPI::BlockVector>;

  template struct SundialsPreconditioner<TrilinosWrappers::MPI::Vector>;
  template struct SundialsPreconditioner<TrilinosWrappers::MPI::BlockVector>;
#      endif
#    endif // DEAL_II_WITH_TRILINOS

#    ifdef DEAL_II_WITH_PETSC
#      ifndef PETSC_USE_COMPLEX
  template class ARKode<PETScWrappers::MPI::Vector>;
  template class ARKode<PETScWrappers::MPI::BlockVector>;

#        if DEAL_II_SUNDIALS_VERSION_GTE(4, 0, 0)
  template class SundialsOperator<PETScWrappers::MPI::Vector>;
  template class SundialsOperator<PETScWrappers::MPI::BlockVector>;

  template class SundialsPreconditioner<PETScWrappers::MPI::Vector>;
  template class SundialsPreconditioner<PETScWrappers::MPI::BlockVector>;
#        endif
#      endif // PETSC_USE_COMPLEX
#    endif   // DEAL_II_WITH_PETSC

#  endif // DEAL_II_WITH_MPI

} // namespace SUNDIALS

DEAL_II_NAMESPACE_CLOSE

#endif
