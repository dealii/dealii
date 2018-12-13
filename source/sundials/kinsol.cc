//-----------------------------------------------------------
//
//    Copyright (C) 2017 - 2018 by the deal.II authors
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
#  include <deal.II/base/utilities.h>

#  include <deal.II/sundials/copy.h>

#  include <sundials/sundials_config.h>
#  if DEAL_II_SUNDIALS_VERSION_GTE(3, 0, 0)
#    include <kinsol/kinsol_direct.h>
#    include <kinsol/kinsol_spils.h>
#    include <sunlinsol/sunlinsol_dense.h>
#    include <sunlinsol/sunlinsol_spfgmr.h>
#    include <sunlinsol/sunlinsol_spgmr.h>
#    include <sunmatrix/sunmatrix_dense.h>
#  else
#    include <kinsol/kinsol_dense.h>
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
    t_kinsol_function(N_Vector yy, N_Vector FF, void *user_data)
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



    template <typename VectorType>
    int
    t_kinsol_setup_jacobian(KINMem kinsol_mem)
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
    t_kinsol_solve_jacobian(KINMem    kinsol_mem,
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

    template <class VectorType>
    int
    t_kinsol_jacobian_vmult(
      N_Vector     v,     /* vector that must be multiplied by the jacobian */
      N_Vector     Jv,    /* where to store the result of the vmult */
      N_Vector     u,     /* current value of the dependent variable vector */
      booleantype *new_u, /* u has changed since last call to this function */
      void *       user_data)
    {
      KINSOL<VectorType> &solver =
        *static_cast<KINSOL<VectorType> *>(user_data);
      GrowingVectorMemory<VectorType> mem;

      typename VectorMemory<VectorType>::Pointer src_v(mem);
      solver.reinit_vector(*src_v);

      typename VectorMemory<VectorType>::Pointer src_u(mem);
      solver.reinit_vector(*src_u);

      typename VectorMemory<VectorType>::Pointer dst_jv(mem);
      solver.reinit_vector(*dst_jv);

      copy(*src_v, v);
      copy(*src_u, u);

      solver.jacobian_vmult(*src_v, *src_u, *dst_jv);

      copy(Jv, *dst_jv);
      *new_u = SUNFALSE;
      return 0;
    }

    template <class VectorType>
    int
    t_kinsol_setup_preconditioner(
      N_Vector u,      /* current unscaled value of the iterate            */
      N_Vector uscale, /* diagonal elements of the scaling matrix for u    */
      N_Vector fval,   /* F(u)                                             */
      N_Vector fscale, /* diagonal elements of the scaling matrix for fval */
      void *   user_data)
    {
      KINSOL<VectorType> &solver =
        *static_cast<KINSOL<VectorType> *>(user_data);
      GrowingVectorMemory<VectorType>            mem;
      typename VectorMemory<VectorType>::Pointer src_u(mem);
      solver.reinit_vector(*src_u);

      typename VectorMemory<VectorType>::Pointer src_f(mem);
      solver.reinit_vector(*src_f);

      copy(*src_u, u);
      copy(*src_f, fval);

      // if the user provided scaling vectors, we apply them
      // otherwise uscale and fscale are vectors of 1
      if (solver.get_solution_scaling)
        {
          typename VectorMemory<VectorType>::Pointer scale(mem);
          solver.reinit_vector(*scale);
          copy(*scale, uscale);
          for (auto i = 0llu; i < scale->size(); ++i)
            (*src_u)[i] *= (*scale)[i];
        }

      if (solver.get_function_scaling)
        {
          typename VectorMemory<VectorType>::Pointer scale(mem);
          solver.reinit_vector(*scale);
          copy(*scale, fscale);
          for (auto i = 0llu; i < scale->size(); ++i)
            (*src_f)[i] *= (*scale)[i];
        }
      solver.setup_preconditioner(*src_u, *src_f);

      return 0;
    }

    template <class VectorType, bool setup_free>
    int
    t_kinsol_solve_preconditioner(
      N_Vector u,      /* current unscaled value of the iterate            */
      N_Vector uscale, /* diagonal elements of the scaling matrix for u    */
      N_Vector fval,   /* F(u)                                             */
      N_Vector fscale, /* diagonal elements of the scaling matrix for fval */
      N_Vector v,      /* on input v=rhs, on output v = solution           */
      void *   user_data)
    {
      KINSOL<VectorType> &solver =
        *static_cast<KINSOL<VectorType> *>(user_data);
      GrowingVectorMemory<VectorType> mem;

      typename VectorMemory<VectorType>::Pointer dst_v(mem);
      solver.reinit_vector(*dst_v);

      typename VectorMemory<VectorType>::Pointer src_v(mem);
      solver.reinit_vector(*src_v);

      copy(*src_v, v);

      if (setup_free)
        {
          typename VectorMemory<VectorType>::Pointer src_u(mem);
          solver.reinit_vector(*src_u);

          typename VectorMemory<VectorType>::Pointer src_f(mem);
          solver.reinit_vector(*src_f);

          copy(*src_u, u);
          copy(*src_f, fval);

          // if the user provided scaling vectors, we apply them
          // otherwise uscale and fscale are vectors of 1
          if (solver.get_solution_scaling)
            {
              typename VectorMemory<VectorType>::Pointer scale(mem);
              solver.reinit_vector(*scale);
              copy(*scale, uscale);
              for (auto i = 0llu; i < scale->size(); ++i)
                (*src_u)[i] *= (*scale)[i];
            }

          if (solver.get_function_scaling)
            {
              typename VectorMemory<VectorType>::Pointer scale(mem);
              solver.reinit_vector(*scale);
              copy(*scale, fscale);
              for (auto i = 0llu; i < scale->size(); ++i)
                (*src_f)[i] *= (*scale)[i];
            }
          solver.solve_preconditioner_setup_free(*src_u,
                                                 *src_f,
                                                 *src_v,
                                                 *dst_v);
        }
      else
        solver.solve_preconditioner(*src_v, *dst_v);

      copy(v, *dst_v);
      return 0;
    }
  } // namespace

  template <typename VectorType>
  KINSOL<VectorType>::KINSOL(const AdditionalData &data,
                             const MPI_Comm        mpi_comm)
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

    int status = KINInit(kinsol_mem, t_kinsol_function<VectorType>, solution);
    static_cast<void>(status);
    AssertKINSOL(status);

    status = KINSetUserData(kinsol_mem, static_cast<void *>(this));
    AssertKINSOL(status);

    status = KINSetPrintLevel(kinsol_mem, data.verbosity);
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

#  if DEAL_II_SUNDIALS_VERSION_GTE(3, 0, 0)
    SUNMatrix       J  = nullptr;
    SUNLinearSolver LS = nullptr;
#  endif

    if (solve_jacobian_system)
      {
        auto KIN_mem        = static_cast<KINMem>(kinsol_mem);
        KIN_mem->kin_lsolve = t_kinsol_solve_jacobian<VectorType>;
        if (setup_jacobian)
          {
            KIN_mem->kin_lsetup = t_kinsol_setup_jacobian<VectorType>;
#  if DEAL_II_SUNDIALS_VERSION_LT(3, 0, 0)
            KIN_mem->kin_setupNonNull = true;
#  endif
          }
      }
    else
      {
        if (data.linear_solver == KINSOL::AdditionalData::LinearSolver::dense)
          {
#  if DEAL_II_SUNDIALS_VERSION_GTE(3, 0, 0)
            J      = SUNDenseMatrix(system_size, system_size);
            LS     = SUNDenseLinearSolver(u_scale, J);
            status = KINDlsSetLinearSolver(kinsol_mem, LS, J);
#  else
            status = KINDense(kinsol_mem, system_size);
#  endif
            AssertKINSOL(status);
          }
        else
          {
#  if DEAL_II_SUNDIALS_VERSION_GTE(3, 0, 0)
            switch (data.linear_solver)
              {
                case KINSOL::AdditionalData::LinearSolver::gmres:
                  if (solve_preconditioner || solve_preconditioner_setup_free)
                    {
                      // KINSOL wants a RIGHT preconditioner
                      // as per documentation, other type of preconditioning
                      // (i.e., left, both) are not well supported
                      LS = SUNSPGMR(u_scale, PREC_RIGHT, -1);
                    }
                  else
                    LS = SUNSPGMR(u_scale, PREC_NONE, -1);
                  break;


                case KINSOL::AdditionalData::LinearSolver::fgmres:
                  if (solve_preconditioner || solve_preconditioner_setup_free)
                    {
                      // KINSOL Flexible GMRES works only with a right
                      // preconditioner
                      LS = SUNSPFGMR(u_scale, PREC_RIGHT, -1);
                    }
                  else
                    LS = SUNSPFGMR(u_scale, PREC_NONE, -1);
                  break;

                default:
                  AssertThrow(
                    false,
                    ExcMessage(
                      "Supported since version 3.0 of Sundials. Please upgrade your library."));
                  break;
              }

            status = KINSpilsSetLinearSolver(kinsol_mem, LS);
            AssertKINSOL(status);
            status =
              KINSpilsSetJacTimesVecFn(kinsol_mem,
                                       t_kinsol_jacobian_vmult<VectorType>);
            AssertKINSOL(status);
            if (solve_preconditioner_setup_free)
              {
                status = KINSpilsSetPreconditioner(
                  kinsol_mem,
                  nullptr,
                  t_kinsol_solve_preconditioner<VectorType, true>);
                AssertKINSOL(status);
              }
            else if (solve_preconditioner and setup_preconditioner)
              {
                status = KINSpilsSetPreconditioner(
                  kinsol_mem,
                  t_kinsol_setup_preconditioner<VectorType>,
                  t_kinsol_solve_preconditioner<VectorType, false>);
                AssertKINSOL(status);
              }

#  else
            AssertThrow(false, ExcNotImplemented());
#  endif
          }
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

#  if DEAL_II_SUNDIALS_VERSION_GTE(3, 0, 0)
    SUNMatDestroy(J);
    SUNLinSolFree(LS);
#  endif
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
    jacobian_vmult = [](const VectorType &, const VectorType &, VectorType &) {
      AssertThrow(false, ExcFunctionNotProvided("jacobian_vmult"));
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
