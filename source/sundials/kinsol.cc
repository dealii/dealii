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

#  include <sundials/sundials_config.h>
#  if DEAL_II_SUNDIALS_VERSION_GTE(3, 0, 0)
#    include <kinsol/kinsol_direct.h>
#    include <sunlinsol/sunlinsol_dense.h>
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
#  if DEAL_II_SUNDIALS_VERSION_GTE(3, 0, 0)
        J      = SUNDenseMatrix(system_size, system_size);
        LS     = SUNDenseLinearSolver(u_scale, J);
        status = KINDlsSetLinearSolver(kinsol_mem, LS, J);
#  else
        status = KINDense(kinsol_mem, system_size);
#  endif
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
