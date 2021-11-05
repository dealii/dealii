//-----------------------------------------------------------
//
//    Copyright (C) 2015 - 2020 by the deal.II authors
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

#  include <sundials/sundials_config.h>
#  if DEAL_II_SUNDIALS_VERSION_LT(4, 0, 0)
#    ifdef DEAL_II_SUNDIALS_WITH_IDAS
#      include <idas/idas_impl.h>
#    else
#      include <ida/ida_impl.h>
#    endif
#  endif
#  if DEAL_II_SUNDIALS_VERSION_LT(5, 0, 0)
#    include <deal.II/sundials/sunlinsol_newempty.h>
#  endif
#  include <iomanip>
#  include <iostream>

DEAL_II_NAMESPACE_OPEN

namespace SUNDIALS
{
  using namespace internal;

  namespace
  {
    template <typename VectorType>
    int
    t_dae_residual(realtype tt,
                   N_Vector yy,
                   N_Vector yp,
                   N_Vector rr,
                   void *   user_data)
    {
      IDA<VectorType> &solver = *static_cast<IDA<VectorType> *>(user_data);

      auto *src_yy   = internal::unwrap_nvector_const<VectorType>(yy);
      auto *src_yp   = internal::unwrap_nvector_const<VectorType>(yp);
      auto *residual = internal::unwrap_nvector<VectorType>(rr);

      int err = solver.residual(tt, *src_yy, *src_yp, *residual);

      return err;
    }



#  if DEAL_II_SUNDIALS_VERSION_LT(4, 0, 0)
    template <typename VectorType>
    int
    t_dae_lsetup(IDAMem   IDA_mem,
                 N_Vector yy,
                 N_Vector yp,
                 N_Vector resp,
                 N_Vector tmp1,
                 N_Vector tmp2,
                 N_Vector tmp3)
    {
      (void)tmp1;
      (void)tmp2;
      (void)tmp3;
      (void)resp;
      IDA<VectorType> &solver =
        *static_cast<IDA<VectorType> *>(IDA_mem->ida_user_data);

      auto *src_yy = internal::unwrap_nvector_const<VectorType>(yy);
      auto *src_yp = internal::unwrap_nvector_const<VectorType>(yp);

      int err = solver.setup_jacobian(IDA_mem->ida_tn,
                                      *src_yy,
                                      *src_yp,
                                      IDA_mem->ida_cj);

      return err;
    }



    template <typename VectorType>
    int
    t_dae_solve(IDAMem   IDA_mem,
                N_Vector b,
                N_Vector weight,
                N_Vector yy,
                N_Vector yp,
                N_Vector resp)
    {
      (void)weight;
      (void)yy;
      (void)yp;
      (void)resp;
      IDA<VectorType> &solver =
        *static_cast<IDA<VectorType> *>(IDA_mem->ida_user_data);
      GrowingVectorMemory<VectorType> mem;

      typename VectorMemory<VectorType>::Pointer dst(mem);
      solver.reinit_vector(*dst);

      auto *src = internal::unwrap_nvector<VectorType>(b);

      int err = solver.solve_jacobian_system(*src, *dst);
      *src    = *dst;

      return err;
    }



#  else
    template <typename VectorType>
    int
    t_dae_jacobian_setup(realtype tt,
                         realtype cj,
                         N_Vector yy,
                         N_Vector yp,
                         N_Vector /* residual */,
                         SUNMatrix /* ignored */,
                         void *user_data,
                         N_Vector /* tmp1 */,
                         N_Vector /* tmp2 */,
                         N_Vector /* tmp3 */)
    {
      Assert(user_data != nullptr, ExcInternalError());
      IDA<VectorType> &solver = *static_cast<IDA<VectorType> *>(user_data);

      auto *src_yy = internal::unwrap_nvector_const<VectorType>(yy);
      auto *src_yp = internal::unwrap_nvector_const<VectorType>(yp);

      int err = solver.setup_jacobian(tt, *src_yy, *src_yp, cj);

      return err;
    }



    template <typename VectorType>
    int
    t_dae_solve_jacobian_system(SUNLinearSolver LS,
                                SUNMatrix /*ignored*/,
                                N_Vector x,
                                N_Vector b,
                                realtype /*tol*/)
    {
      const IDA<VectorType> &solver =
        *static_cast<const IDA<VectorType> *>(LS->content);

      auto *src_b = internal::unwrap_nvector_const<VectorType>(b);
      auto *dst_x = internal::unwrap_nvector<VectorType>(x);

      const int err = solver.solve_jacobian_system(*src_b, *dst_x);

      return err;
    }


    template <typename VectorType>
    int
    t_dae_solve_with_jacobian(SUNLinearSolver LS,
                              SUNMatrix /*ignored*/,
                              N_Vector x,
                              N_Vector b,
                              realtype tol)
    {
      IDA<VectorType> &solver = *static_cast<IDA<VectorType> *>(LS->content);

      auto *src_b = internal::unwrap_nvector_const<VectorType>(b);
      auto *dst_x = internal::unwrap_nvector<VectorType>(x);

      int       n_iter;
      const int err = solver.solve_with_jacobian(*src_b, *dst_x, n_iter, tol);
      solver.set_n_iter(n_iter > 0 ? n_iter : 1);

      return err;
    }
#  endif
  } // namespace



  template <typename VectorType>
  IDA<VectorType>::IDA(const AdditionalData &data, const MPI_Comm mpi_comm)
    : data(data)
    , ida_mem(nullptr)
    , communicator(is_serial_vector<VectorType>::value ?
                     MPI_COMM_SELF :
                     Utilities::MPI::duplicate_communicator(mpi_comm))
  {
    set_functions_to_trigger_an_assert();
  }

  template <typename VectorType>
  IDA<VectorType>::~IDA()
  {
    if (ida_mem)
      IDAFree(&ida_mem);
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
  IDA<VectorType>::solve_dae(VectorType &solution, VectorType &solution_dot)
  {
    double       t           = data.initial_time;
    double       h           = data.initial_step_size;
    unsigned int step_number = 0;

    this->n_iter = 1;
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

        auto yy = internal::make_nvector_view(solution);
        auto yp = internal::make_nvector_view(solution_dot);

        status = IDASolve(ida_mem, next_time, &t, yy, yp, IDA_NORMAL);
        AssertIDA(status);

        status = IDAGetLastStep(ida_mem, &h);
        AssertIDA(status);

        while (solver_should_restart(t, solution, solution_dot))
          reset(t, h, solution, solution_dot);

        step_number++;

        output_step(t, solution, solution_dot, step_number);
      }

    return step_number;
  }



  template <typename VectorType>
  void
  IDA<VectorType>::reset(const double current_time,
                         const double current_time_step,
                         VectorType & solution,
                         VectorType & solution_dot)
  {
    bool first_step = (current_time == data.initial_time);

    if (ida_mem)
      IDAFree(&ida_mem);

    ida_mem = IDACreate();

    int status;
    (void)status;

    auto yy = internal::make_nvector_view(solution);
    auto yp = internal::make_nvector_view(solution_dot);

    status = IDAInit(ida_mem, t_dae_residual<VectorType>, current_time, yy, yp);
    AssertIDA(status);
    if (get_local_tolerances)
      {
        const auto abs_tols =
          internal::make_nvector_view(get_local_tolerances());
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
        auto dc          = differential_components();
        for (auto i = dc.begin(); i != dc.end(); ++i)
          diff_comp_vector[*i] = 1.0;
        diff_comp_vector.compress(VectorOperation::insert);

        const auto diff_id = internal::make_nvector_view(diff_comp_vector);
        status             = IDASetId(ida_mem, diff_id);
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
#  if DEAL_II_SUNDIALS_VERSION_LT(4, 0, 0)
    auto IDA_mem = static_cast<IDAMem>(ida_mem);

    IDA_mem->ida_lsetup = t_dae_lsetup<VectorType>;

    if (solve_jacobian_system)
      IDA_mem->ida_lsolve = t_dae_solve<VectorType>;
    else
      AssertThrow(false, ExcFunctionNotProvided("solve_jacobian_system"));
#    if DEAL_II_SUNDIALS_VERSION_LT(3, 0, 0)
    IDA_mem->ida_setupNonNull = true;
#    endif
#  else
    SUNMatrix       J  = nullptr;
    SUNLinearSolver LS = nullptr;

    // and attach it to the SUNLinSol object. The functions that will get
    // called do not actually receive the IDAMEM object, just the LS
    // object, so we have to store a pointer to the current
    // object in the LS object
    LS          = SUNLinSolNewEmpty();
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
          free(LS->ops);
          LS->ops = nullptr;
        }
      free(LS);
      LS = nullptr;
      return 0;
    };

    if (solve_with_jacobian)
      {
        LS->ops->solve = t_dae_solve_with_jacobian<VectorType>;
      }
    else if (solve_jacobian_system)
      {
        LS->ops->solve = t_dae_solve_jacobian_system<VectorType>;
      }
    else
      {
        AssertThrow(false, ExcFunctionNotProvided("solve_with_jacobian"));
      }
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
    LS->ops->numiters = [](SUNLinearSolver LS) -> int {
      IDA<VectorType> &solver = *static_cast<IDA<VectorType> *>(LS->content);
      return solver.get_n_iter();
    };
    // Even though we don't use it, IDA still wants us to set some
    // kind of matrix object for the nonlinear solver. This is because
    // if we don't set it, it won't call the functions that set up
    // the matrix object (i.e., the argument to the 'IDASetJacFn'
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
    status = IDASetLinearSolver(ida_mem, LS, J);
    AssertIDA(status);

    status = IDASetLSNormFactor(ida_mem, data.ls_norm_factor);
    AssertIDA(status);
    // Finally tell IDA about
    // it as well. The manual says that this must happen *after*
    // calling IDASetLinearSolver
    status = IDASetJacFn(ida_mem, &t_dae_jacobian_setup<VectorType>);
    AssertIDA(status);
#  endif
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
      const unsigned int size = v->size();
      return complete_index_set(size);
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
