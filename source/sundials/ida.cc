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

#include <deal.II/sundials/ida.h>

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

#  ifdef DEAL_II_SUNDIALS_WITH_IDAS
#    include <idas/idas_impl.h>
#  else
#    include <ida/ida_impl.h>
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
      GrowingVectorMemory<VectorType> mem;

      typename VectorMemory<VectorType>::Pointer src_yy(mem);
      solver.reinit_vector(*src_yy);

      typename VectorMemory<VectorType>::Pointer src_yp(mem);
      solver.reinit_vector(*src_yp);

      typename VectorMemory<VectorType>::Pointer residual(mem);
      solver.reinit_vector(*residual);

      copy(*src_yy, yy);
      copy(*src_yp, yp);

      int err = solver.residual(tt, *src_yy, *src_yp, *residual);

      copy(rr, *residual);

      return err;
    }



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
      GrowingVectorMemory<VectorType> mem;

      typename VectorMemory<VectorType>::Pointer src_yy(mem);
      solver.reinit_vector(*src_yy);

      typename VectorMemory<VectorType>::Pointer src_yp(mem);
      solver.reinit_vector(*src_yp);

      copy(*src_yy, yy);
      copy(*src_yp, yp);

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

      typename VectorMemory<VectorType>::Pointer src(mem);
      solver.reinit_vector(*src);

      typename VectorMemory<VectorType>::Pointer dst(mem);
      solver.reinit_vector(*dst);

      copy(*src, b);

      int err = solver.solve_jacobian_system(*src, *dst);
      copy(b, *dst);

      return err;
    }

  } // namespace

  template <typename VectorType>
  IDA<VectorType>::IDA(const AdditionalData &data, const MPI_Comm mpi_comm)
    : data(data)
    , ida_mem(nullptr)
    , yy(nullptr)
    , yp(nullptr)
    , abs_tolls(nullptr)
    , diff_id(nullptr)
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
    unsigned int system_size = solution.size();

    double       t           = data.initial_time;
    double       h           = data.initial_step_size;
    unsigned int step_number = 0;

    int status;
    (void)status;

    // The solution is stored in
    // solution. Here we take only a
    // view of it.
#  ifdef DEAL_II_WITH_MPI
    if (is_serial_vector<VectorType>::value == false)
      {
        const IndexSet    is                = solution.locally_owned_elements();
        const std::size_t local_system_size = is.n_elements();

        yy = N_VNew_Parallel(communicator, local_system_size, system_size);

        yp = N_VNew_Parallel(communicator, local_system_size, system_size);

        diff_id = N_VNew_Parallel(communicator, local_system_size, system_size);

        abs_tolls =
          N_VNew_Parallel(communicator, local_system_size, system_size);
      }
    else
#  endif
      {
        Assert(is_serial_vector<VectorType>::value,
               ExcInternalError(
                 "Trying to use a serial code with a parallel vector."));
        yy        = N_VNew_Serial(system_size);
        yp        = N_VNew_Serial(system_size);
        diff_id   = N_VNew_Serial(system_size);
        abs_tolls = N_VNew_Serial(system_size);
      }
    reset(data.initial_time, data.initial_step_size, solution, solution_dot);

    double next_time = data.initial_time;

    output_step(0, solution, solution_dot, 0);

    while (t < data.final_time)
      {
        next_time += data.output_period;

        status = IDASolve(ida_mem, next_time, &t, yy, yp, IDA_NORMAL);
        AssertIDA(status);

        status = IDAGetLastStep(ida_mem, &h);
        AssertIDA(status);

        copy(solution, yy);
        copy(solution_dot, yp);

        while (solver_should_restart(t, solution, solution_dot))
          reset(t, h, solution, solution_dot);

        step_number++;

        output_step(t, solution, solution_dot, step_number);
      }

      // Free the vectors which are no longer used.
#  ifdef DEAL_II_WITH_MPI
    if (is_serial_vector<VectorType>::value == false)
      {
        N_VDestroy_Parallel(yy);
        N_VDestroy_Parallel(yp);
        N_VDestroy_Parallel(abs_tolls);
        N_VDestroy_Parallel(diff_id);
      }
    else
#  endif
      {
        N_VDestroy_Serial(yy);
        N_VDestroy_Serial(yp);
        N_VDestroy_Serial(abs_tolls);
        N_VDestroy_Serial(diff_id);
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
    unsigned int system_size;
    bool         first_step = (current_time == data.initial_time);

    if (ida_mem)
      IDAFree(&ida_mem);

    ida_mem = IDACreate();


    // Free the vectors which are no longer used.
    if (yy)
      {
#  ifdef DEAL_II_WITH_MPI
        if (is_serial_vector<VectorType>::value == false)
          {
            N_VDestroy_Parallel(yy);
            N_VDestroy_Parallel(yp);
            N_VDestroy_Parallel(abs_tolls);
            N_VDestroy_Parallel(diff_id);
          }
        else
#  endif
          {
            N_VDestroy_Serial(yy);
            N_VDestroy_Serial(yp);
            N_VDestroy_Serial(abs_tolls);
            N_VDestroy_Serial(diff_id);
          }
      }

    int status;
    (void)status;
    system_size = solution.size();
#  ifdef DEAL_II_WITH_MPI
    if (is_serial_vector<VectorType>::value == false)
      {
        const IndexSet    is                = solution.locally_owned_elements();
        const std::size_t local_system_size = is.n_elements();

        yy = N_VNew_Parallel(communicator, local_system_size, system_size);

        yp = N_VNew_Parallel(communicator, local_system_size, system_size);

        diff_id = N_VNew_Parallel(communicator, local_system_size, system_size);

        abs_tolls =
          N_VNew_Parallel(communicator, local_system_size, system_size);
      }
    else
#  endif
      {
        yy        = N_VNew_Serial(system_size);
        yp        = N_VNew_Serial(system_size);
        diff_id   = N_VNew_Serial(system_size);
        abs_tolls = N_VNew_Serial(system_size);
      }

    copy(yy, solution);
    copy(yp, solution_dot);

    status = IDAInit(ida_mem, t_dae_residual<VectorType>, current_time, yy, yp);
    AssertIDA(status);

    if (get_local_tolerances)
      {
        copy(abs_tolls, get_local_tolerances());
        status = IDASVtolerances(ida_mem, data.relative_tolerance, abs_tolls);
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

        copy(diff_id, diff_comp_vector);
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
    auto IDA_mem = static_cast<IDAMem>(ida_mem);

    IDA_mem->ida_lsetup = t_dae_lsetup<VectorType>;
    IDA_mem->ida_lsolve = t_dae_solve<VectorType>;
#  if DEAL_II_SUNDIALS_VERSION_LT(3, 0, 0)
    IDA_mem->ida_setupNonNull = true;
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

        copy(solution, yy);
        copy(solution_dot, yp);
      }
    else if (type == AdditionalData::use_y_diff)
      {
        status =
          IDACalcIC(ida_mem, IDA_YA_YDP_INIT, current_time + current_time_step);
        AssertIDA(status);

        status = IDAGetConsistentIC(ida_mem, yy, yp);
        AssertIDA(status);

        copy(solution, yy);
        copy(solution_dot, yp);
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

    setup_jacobian = [](const double,
                        const VectorType &,
                        const VectorType &,
                        const double) -> int {
      int ret = 0;
      AssertThrow(false, ExcFunctionNotProvided("setup_jacobian"));
      return ret;
    };

    solve_jacobian_system = [](const VectorType &, VectorType &) -> int {
      int ret = 0;
      AssertThrow(false, ExcFunctionNotProvided("solve_jacobian_system"));
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
