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

#include <deal.II/lac/vector_operation.h>

#include <deal.II/sundials/ida.h>

#ifdef DEAL_II_WITH_SUNDIALS
#  include <deal.II/base/utilities.h>

#  include <deal.II/lac/block_vector.h>

#  include <deal.II/sundials/n_vector.h>
#  include <deal.II/sundials/sunlinsol_wrapper.h>
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

#  include <idas/idas.h>
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
#    include <sundials/sundials_context.h>
#  endif

#  include <iomanip>
#  include <iostream>

DEAL_II_NAMESPACE_OPEN

namespace SUNDIALS
{
  template <typename VectorType>
  IDA<VectorType>::IDA(const AdditionalData &data)
    : IDA(data, MPI_COMM_SELF)
  {}



  template <typename VectorType>
  IDA<VectorType>::IDA(const AdditionalData &data, const MPI_Comm mpi_comm)
    : data(data)
    , ida_mem(nullptr)
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    , ida_ctx(nullptr)
#  endif
    , mpi_communicator(mpi_comm)
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
                        &ida_ctx);
    (void)status;
    AssertIDA(status);
#  elif DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    const int status =
      SUNContext_Create(mpi_communicator == MPI_COMM_SELF ? nullptr :
                                                            &mpi_communicator,
                        &ida_ctx);
    (void)status;
    AssertIDA(status);
#  endif
  }



  template <typename VectorType>
  IDA<VectorType>::~IDA()
  {
    IDAFree(&ida_mem);
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    const int status = SUNContext_Free(&ida_ctx);
    (void)status;
    AssertIDA(status);
#  endif

    Assert(pending_exception == nullptr, ExcInternalError());
  }



  template <typename VectorType>
  unsigned int
  IDA<VectorType>::solve_dae(VectorType &solution, VectorType &solution_dot)
  {
    double       t           = data.initial_time;
    double       h           = data.initial_step_size;
    unsigned int step_number = 0;

    int status;
    (void)status;

    reset(data.initial_time, data.initial_step_size, solution, solution_dot);

    // The solution is stored in
    // solution. Here we take only a
    // view of it.

    double next_time = data.initial_time;

    output_step(0, solution, solution_dot, 0);

    while (t < data.final_time)
      {
        next_time += data.output_period;

        auto yy = internal::make_nvector_view(solution
#  if !DEAL_II_SUNDIALS_VERSION_LT(6, 0, 0)
                                              ,
                                              ida_ctx
#  endif
        );
        auto yp = internal::make_nvector_view(solution_dot
#  if !DEAL_II_SUNDIALS_VERSION_LT(6, 0, 0)
                                              ,
                                              ida_ctx
#  endif
        );

        // Execute time steps. If we ended up with a pending exception,
        // see if it was recoverable; if it was, and if IDA recovered,
        // just continue on. If IDA did not recover, rethrow the exception.
        // Do the same if the exception was not recoverable.
        status = IDASolve(ida_mem, next_time, &t, yy, yp, IDA_NORMAL);
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
                  /* just eat the exception and continue */;
                else
                  throw;
              }
            catch (...)
              {
                pending_exception = nullptr;
                throw;
              }
          }
        AssertIDA(status);

        status = IDAGetLastStep(ida_mem, &h);
        AssertIDA(status);

        while (solver_should_restart(t, solution, solution_dot))
          reset(t, h, solution, solution_dot);

        ++step_number;

        output_step(t, solution, solution_dot, step_number);
      }
    long int n_steps;
    status = IDAGetNumSteps(ida_mem, &n_steps);
    AssertIDA(status);
    return n_steps;
  }



  template <typename VectorType>
  void
  IDA<VectorType>::reset(const double current_time,
                         const double current_time_step,
                         VectorType  &solution,
                         VectorType  &solution_dot)
  {
    bool first_step = (current_time == data.initial_time);

    int status;
    (void)status;

#  if DEAL_II_SUNDIALS_VERSION_GTE(7, 0, 0)
    status = SUNContext_Free(&ida_ctx);
    AssertIDA(status);

    // Same comment applies as in class constructor:
    status =
      SUNContext_Create(mpi_communicator == MPI_COMM_SELF ? SUN_COMM_NULL :
                                                            mpi_communicator,
                        &ida_ctx);
    AssertIDA(status);
#  elif DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    status = SUNContext_Free(&ida_ctx);
    AssertIDA(status);

    // Same comment applies as in class constructor:
    status =
      SUNContext_Create(mpi_communicator == MPI_COMM_SELF ? nullptr :
                                                            &mpi_communicator,
                        &ida_ctx);
    AssertIDA(status);
#  endif

    if (ida_mem)
      {
        IDAFree(&ida_mem);
        // Initialization is version-dependent: do that in a moment
      }

#  if DEAL_II_SUNDIALS_VERSION_LT(6, 0, 0)
    ida_mem = IDACreate();
#  else
    ida_mem = IDACreate(ida_ctx);
#  endif

    auto yy = internal::make_nvector_view(solution
#  if !DEAL_II_SUNDIALS_VERSION_LT(6, 0, 0)
                                          ,
                                          ida_ctx
#  endif
    );
    auto yp = internal::make_nvector_view(solution_dot
#  if !DEAL_II_SUNDIALS_VERSION_LT(6, 0, 0)
                                          ,
                                          ida_ctx
#  endif
    );

    status = IDAInit(
      ida_mem,
      [](SUNDIALS::realtype tt,
         N_Vector           yy,
         N_Vector           yp,
         N_Vector           rr,
         void              *user_data) -> int {
        IDA<VectorType> &solver = *static_cast<IDA<VectorType> *>(user_data);

        auto *src_yy   = internal::unwrap_nvector_const<VectorType>(yy);
        auto *src_yp   = internal::unwrap_nvector_const<VectorType>(yp);
        auto *residual = internal::unwrap_nvector<VectorType>(rr);

        return Utilities::call_and_possibly_capture_exception(
          solver.residual,
          solver.pending_exception,
          tt,
          *src_yy,
          *src_yp,
          *residual);
      },
      current_time,
      yy,
      yp);
    AssertIDA(status);
    if (get_local_tolerances)
      {
        const auto abs_tols = internal::make_nvector_view(get_local_tolerances()
#  if !DEAL_II_SUNDIALS_VERSION_LT(6, 0, 0)
                                                            ,
                                                          ida_ctx
#  endif
        );
        status = IDASVtolerances(ida_mem, data.relative_tolerance, abs_tols);
        AssertIDA(status);
      }
    else
      {
        status = IDASStolerances(ida_mem,
                                 data.relative_tolerance,
                                 data.absolute_tolerance);
        AssertIDA(status);
      }

    status = IDASetInitStep(ida_mem, current_time_step);
    AssertIDA(status);

    status = IDASetUserData(ida_mem, this);
    AssertIDA(status);

    if (data.ic_type == AdditionalData::use_y_diff ||
        data.reset_type == AdditionalData::use_y_diff ||
        data.ignore_algebraic_terms_for_errors)
      {
        VectorType diff_comp_vector(solution);
        diff_comp_vector = 0.0;
        for (const auto &component : differential_components())
          diff_comp_vector[component] = 1.0;
        diff_comp_vector.compress(VectorOperation::insert);

        const auto diff_id = internal::make_nvector_view(diff_comp_vector
#  if !DEAL_II_SUNDIALS_VERSION_LT(6, 0, 0)
                                                         ,
                                                         ida_ctx
#  endif
        );
        status = IDASetId(ida_mem, diff_id);
        AssertIDA(status);
      }

    status = IDASetSuppressAlg(ida_mem, data.ignore_algebraic_terms_for_errors);
    AssertIDA(status);

    //  status = IDASetMaxNumSteps(ida_mem, max_steps);
    status = IDASetStopTime(ida_mem, data.final_time);
    AssertIDA(status);

    status = IDASetMaxNonlinIters(ida_mem, data.maximum_non_linear_iterations);
    AssertIDA(status);

    // Initialize solver
    SUNMatrix       J  = nullptr;
    SUNLinearSolver LS = nullptr;

    // and attach it to the SUNLinSol object. The functions that will get
    // called do not actually receive the IDAMEM object, just the LS
    // object, so we have to store a pointer to the current
    // object in the LS object
#  if DEAL_II_SUNDIALS_VERSION_LT(6, 0, 0)
    LS = SUNLinSolNewEmpty();
#  else
    LS      = SUNLinSolNewEmpty(ida_ctx);
#  endif

    LS->content = this;

    LS->ops->gettype = [](SUNLinearSolver /*ignored*/) -> SUNLinearSolver_Type {
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

    AssertThrow(solve_with_jacobian,
                ExcFunctionNotProvided("solve_with_jacobian"));
    LS->ops->solve = [](SUNLinearSolver LS,
                        SUNMatrix /*ignored*/,
                        N_Vector           x,
                        N_Vector           b,
                        SUNDIALS::realtype tol) -> int {
      IDA<VectorType> &solver = *static_cast<IDA<VectorType> *>(LS->content);

      auto *src_b = internal::unwrap_nvector_const<VectorType>(b);
      auto *dst_x = internal::unwrap_nvector<VectorType>(x);
      return Utilities::call_and_possibly_capture_exception(
        solver.solve_with_jacobian,
        solver.pending_exception,
        *src_b,
        *dst_x,
        tol);
    };

    // When we set an iterative solver IDA requires that resid is provided. From
    // SUNDIALS docs If an iterative method computes the preconditioned initial
    // residual and returns with a successful solve without performing any
    // iterations (i.e., either the initial guess or the preconditioner is
    // sufficiently accurate), then this optional routine may be called by the
    // SUNDIALS package. This routine should return the N_Vector containing the
    // preconditioned initial residual vector.
    LS->ops->resid = [](SUNLinearSolver /*ignored*/) -> N_Vector {
      return nullptr;
    };
    // When we set an iterative solver IDA requires that last number of
    // iteration is provided. Since we can't know what kind of solver the user
    // has provided we set 1. This is clearly suboptimal.
    LS->ops->numiters = [](SUNLinearSolver /*ignored*/) -> int { return 1; };
    // Even though we don't use it, IDA still wants us to set some
    // kind of matrix object for the nonlinear solver. This is because
    // if we don't set it, it won't call the functions that set up
    // the matrix object (i.e., the argument to the 'IDASetJacFn'
    // function below).
#  if DEAL_II_SUNDIALS_VERSION_LT(6, 0, 0)
    J = SUNMatNewEmpty();
#  else
    J       = SUNMatNewEmpty(ida_ctx);
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
    status = IDASetLinearSolver(ida_mem, LS, J);
    AssertIDA(status);

    status = IDASetLSNormFactor(ida_mem, data.ls_norm_factor);
    AssertIDA(status);
    // Finally tell IDA about
    // it as well. The manual says that this must happen *after*
    // calling IDASetLinearSolver
    status = IDASetJacFn(
      ida_mem,
      [](SUNDIALS::realtype tt,
         SUNDIALS::realtype cj,
         N_Vector           yy,
         N_Vector           yp,
         N_Vector /* residual */,
         SUNMatrix /* ignored */,
         void *user_data,
         N_Vector /* tmp1 */,
         N_Vector /* tmp2 */,
         N_Vector /* tmp3 */) -> int {
        Assert(user_data != nullptr, ExcInternalError());
        IDA<VectorType> &solver = *static_cast<IDA<VectorType> *>(user_data);

        auto *src_yy = internal::unwrap_nvector_const<VectorType>(yy);
        auto *src_yp = internal::unwrap_nvector_const<VectorType>(yp);

        return Utilities::call_and_possibly_capture_exception(
          solver.setup_jacobian,
          solver.pending_exception,
          tt,
          *src_yy,
          *src_yp,
          cj);
      });
    AssertIDA(status);
    status = IDASetMaxOrd(ida_mem, data.maximum_order);
    AssertIDA(status);

    typename AdditionalData::InitialConditionCorrection type;
    if (first_step)
      type = data.ic_type;
    else
      type = data.reset_type;

    status =
      IDASetMaxNumItersIC(ida_mem, data.maximum_non_linear_iterations_ic);
    AssertIDA(status);

    if (type == AdditionalData::use_y_dot)
      {
        // (re)initialization of the vectors
        status =
          IDACalcIC(ida_mem, IDA_Y_INIT, current_time + current_time_step);
        AssertIDA(status);

        status = IDAGetConsistentIC(ida_mem, yy, yp);
        AssertIDA(status);
      }
    else if (type == AdditionalData::use_y_diff)
      {
        status =
          IDACalcIC(ida_mem, IDA_YA_YDP_INIT, current_time + current_time_step);
        AssertIDA(status);

        status = IDAGetConsistentIC(ida_mem, yy, yp);
        AssertIDA(status);
      }
  }

  template <typename VectorType>
  void
  IDA<VectorType>::set_functions_to_trigger_an_assert()
  {
    reinit_vector = [](VectorType &) {
      AssertThrow(false, ExcFunctionNotProvided("reinit_vector"));
    };

    residual = [](const double,
                  const VectorType &,
                  const VectorType &,
                  VectorType &) -> int {
      int ret = 0;
      AssertThrow(false, ExcFunctionNotProvided("residual"));
      return ret;
    };


    output_step = [](const double,
                     const VectorType &,
                     const VectorType &,
                     const unsigned int) { return; };

    solver_should_restart =
      [](const double, VectorType &, VectorType &) -> bool { return false; };

    differential_components = [&]() -> IndexSet {
      GrowingVectorMemory<VectorType>            mem;
      typename VectorMemory<VectorType>::Pointer v(mem);
      reinit_vector(*v);
      return v->locally_owned_elements();
    };
  }

  template class IDA<Vector<double>>;
  template class IDA<BlockVector<double>>;

#  ifdef DEAL_II_WITH_MPI

#    ifdef DEAL_II_WITH_TRILINOS
  template class IDA<TrilinosWrappers::MPI::Vector>;
  template class IDA<TrilinosWrappers::MPI::BlockVector>;
#    endif // DEAL_II_WITH_TRILINOS

#    ifdef DEAL_II_WITH_PETSC
#      ifndef PETSC_USE_COMPLEX
  template class IDA<PETScWrappers::MPI::Vector>;
  template class IDA<PETScWrappers::MPI::BlockVector>;
#      endif // PETSC_USE_COMPLEX
#    endif   // DEAL_II_WITH_PETSC

#  endif // DEAL_II_WITH_MPI

} // namespace SUNDIALS

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SUNDIALS
