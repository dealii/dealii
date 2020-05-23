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

#  include <arkode/arkode_impl.h>
#  include <sundials/sundials_config.h>

#  include <iomanip>
#  include <iostream>

// Make sure we know how to call sundials own ARKode() function
const auto &SundialsARKode = ARKode;

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
      ARKode<VectorType> &solver =
        *static_cast<ARKode<VectorType> *>(user_data);
      GrowingVectorMemory<VectorType> mem;

      typename VectorMemory<VectorType>::Pointer src_yy(mem);
      solver.reinit_vector(*src_yy);

      typename VectorMemory<VectorType>::Pointer dst_yp(mem);
      solver.reinit_vector(*dst_yp);

      copy(*src_yy, yy);

      int err = solver.explicit_function(tt, *src_yy, *dst_yp);

      copy(yp, *dst_yp);

      return err;
    }



    template <typename VectorType>
    int
    t_arkode_implicit_function(realtype tt,
                               N_Vector yy,
                               N_Vector yp,
                               void *   user_data)
    {
      ARKode<VectorType> &solver =
        *static_cast<ARKode<VectorType> *>(user_data);
      GrowingVectorMemory<VectorType> mem;

      typename VectorMemory<VectorType>::Pointer src_yy(mem);
      solver.reinit_vector(*src_yy);

      typename VectorMemory<VectorType>::Pointer dst_yp(mem);
      solver.reinit_vector(*dst_yp);

      copy(*src_yy, yy);

      int err = solver.implicit_function(tt, *src_yy, *dst_yp);

      copy(yp, *dst_yp);

      return err;
    }



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
      ARKode<VectorType> &solver =
        *static_cast<ARKode<VectorType> *>(arkode_mem->ark_user_data);
      GrowingVectorMemory<VectorType> mem;

      typename VectorMemory<VectorType>::Pointer src_ypred(mem);
      solver.reinit_vector(*src_ypred);

      typename VectorMemory<VectorType>::Pointer src_fpred(mem);
      solver.reinit_vector(*src_fpred);

      copy(*src_ypred, ypred);
      copy(*src_fpred, fpred);

      // avoid reinterpret_cast
      bool jcurPtr_tmp = false;
      int  err         = solver.setup_jacobian(convfail,
                                      arkode_mem->ark_tn,
                                      arkode_mem->ark_gamma,
                                      *src_ypred,
                                      *src_fpred,
                                      jcurPtr_tmp);
#  if DEAL_II_SUNDIALS_VERSION_GTE(2, 0, 0)
      *jcurPtr = jcurPtr_tmp ? SUNTRUE : SUNFALSE;
#  else
      *jcurPtr = jcurPtr_tmp ? TRUE : FALSE;
#  endif

      return err;
    }



    template <typename VectorType>
    int
    t_arkode_solve_jacobian(ARKodeMem arkode_mem,
                            N_Vector  b,
#  if DEAL_II_SUNDIALS_VERSION_LT(3, 0, 0)
                            N_Vector,
#  endif
                            N_Vector ycur,
                            N_Vector fcur)
    {
      ARKode<VectorType> &solver =
        *static_cast<ARKode<VectorType> *>(arkode_mem->ark_user_data);
      GrowingVectorMemory<VectorType> mem;

      typename VectorMemory<VectorType>::Pointer src(mem);
      solver.reinit_vector(*src);

      typename VectorMemory<VectorType>::Pointer src_ycur(mem);
      solver.reinit_vector(*src_ycur);

      typename VectorMemory<VectorType>::Pointer src_fcur(mem);
      solver.reinit_vector(*src_fcur);

      typename VectorMemory<VectorType>::Pointer dst(mem);
      solver.reinit_vector(*dst);

      copy(*src, b);
      copy(*src_ycur, ycur);
      copy(*src_fcur, fcur);

      int err = solver.solve_jacobian_system(arkode_mem->ark_tn,
                                             arkode_mem->ark_gamma,
                                             *src_ycur,
                                             *src_fcur,
                                             *src,
                                             *dst);
      copy(b, *dst);

      return err;
    }



    template <typename VectorType>
    int
    t_arkode_setup_mass(ARKodeMem arkode_mem, N_Vector, N_Vector, N_Vector)
    {
      ARKode<VectorType> &solver =
        *static_cast<ARKode<VectorType> *>(arkode_mem->ark_user_data);
      int err = solver.setup_mass(arkode_mem->ark_tn);
      return err;
    }



    template <typename VectorType>
    int
    t_arkode_solve_mass(ARKodeMem arkode_mem,
#  if DEAL_II_SUNDIALS_VERSION_LT(3, 0, 0)
                        N_Vector b,
                        N_Vector
#  else
                        N_Vector b
#  endif
    )
    {
      ARKode<VectorType> &solver =
        *static_cast<ARKode<VectorType> *>(arkode_mem->ark_user_data);
      GrowingVectorMemory<VectorType> mem;

      typename VectorMemory<VectorType>::Pointer src(mem);
      solver.reinit_vector(*src);

      typename VectorMemory<VectorType>::Pointer dst(mem);
      solver.reinit_vector(*dst);

      copy(*src, b);

      int err = solver.solve_mass_system(*src, *dst);
      copy(b, *dst);

      return err;
    }
  } // namespace

  template <typename VectorType>
  ARKode<VectorType>::ARKode(const AdditionalData &data,
                             const MPI_Comm        mpi_comm)
    : data(data)
    , arkode_mem(nullptr)
    , yy(nullptr)
    , abs_tolls(nullptr)
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
      ARKodeFree(&arkode_mem);
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
        abs_tolls = N_VNew_Serial(system_size);
      }
    reset(data.initial_time, data.initial_step_size, solution);

    double next_time = data.initial_time;

    if (output_step)
      output_step(0, solution, 0);

    while (t < data.final_time)
      {
        next_time += data.output_period;

        status = SundialsARKode(arkode_mem, next_time, yy, &t, ARK_NORMAL);

        AssertARKode(status);

        status = ARKodeGetLastStep(arkode_mem, &h);
        AssertARKode(status);

        copy(solution, yy);

        while (solver_should_restart(t, solution))
          reset(t, h, solution);

        step_number++;

        if (output_step)
          output_step(t, solution, step_number);
      }

      // Free the vectors which are no longer used.
#  ifdef DEAL_II_WITH_MPI
    if (is_serial_vector<VectorType>::value == false)
      {
        N_VDestroy_Parallel(yy);
        N_VDestroy_Parallel(abs_tolls);
      }
    else
#  endif
      {
        N_VDestroy_Serial(yy);
        N_VDestroy_Serial(abs_tolls);
      }

    return step_number;
  }

  template <typename VectorType>
  void
  ARKode<VectorType>::reset(const double      current_time,
                            const double      current_time_step,
                            const VectorType &solution)
  {
    unsigned int system_size;

    if (arkode_mem)
      ARKodeFree(&arkode_mem);

    arkode_mem = ARKodeCreate();

    // Free the vectors which are no longer used.
    if (yy)
      {
#  ifdef DEAL_II_WITH_MPI
        if (is_serial_vector<VectorType>::value == false)
          {
            N_VDestroy_Parallel(yy);
            N_VDestroy_Parallel(abs_tolls);
          }
        else
#  endif
          {
            N_VDestroy_Serial(yy);
            N_VDestroy_Serial(abs_tolls);
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

        abs_tolls =
          N_VNew_Parallel(communicator, local_system_size, system_size);
      }
    else
#  endif
      {
        yy        = N_VNew_Serial(system_size);
        abs_tolls = N_VNew_Serial(system_size);
      }

    copy(yy, solution);

    Assert(explicit_function || implicit_function,
           ExcFunctionNotProvided("explicit_function || implicit_function"));

    status = ARKodeInit(
      arkode_mem,
      explicit_function ? &t_arkode_explicit_function<VectorType> : nullptr,
      implicit_function ? &t_arkode_implicit_function<VectorType> : nullptr,
      current_time,
      yy);
    AssertARKode(status);

    if (get_local_tolerances)
      {
        copy(abs_tolls, get_local_tolerances());
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
#  if DEAL_II_SUNDIALS_VERSION_LT(3, 0, 0)
            ARKode_mem->ark_setupNonNull = true;
#  endif
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
#  if DEAL_II_SUNDIALS_VERSION_LT(3, 0, 0)
            ARKode_mem->ark_MassSetupNonNull = true;
#  endif
          }
      }

    status = ARKodeSetOrder(arkode_mem, data.maximum_order);
    AssertARKode(status);
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

  template class ARKode<Vector<double>>;
  template class ARKode<BlockVector<double>>;

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
