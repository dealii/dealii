//-----------------------------------------------------------
//
//    Copyright (C) 2015 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE at
//    the top level of the deal.II distribution.
//
//-----------------------------------------------------------


#include <deal.II/sundials/ida_interface.h>
#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SUNDIALS

#include <deal.II/base/utilities.h>
#include <deal.II/lac/block_vector.h>
#ifdef DEAL_II_WITH_TRILINOS
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#endif
#ifdef DEAL_II_WITH_PETSC
#include <deal.II/lac/petsc_parallel_block_vector.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#endif
#include <deal.II/base/utilities.h>
#include <deal.II/sundials/copy.h>

#include <iostream>
#include <iomanip>
#include <ida/ida_impl.h>

DEAL_II_NAMESPACE_OPEN

namespace SUNDIALS
{
  using namespace internal;

  namespace
  {
    template<typename VectorType>
    int t_dae_residual(realtype tt, N_Vector yy, N_Vector yp,
                       N_Vector rr, void *user_data)
    {
      IDAInterface<VectorType> &solver = *static_cast<IDAInterface<VectorType> *>(user_data);

      std::shared_ptr<VectorType> src_yy = solver.create_new_vector();
      std::shared_ptr<VectorType> src_yp = solver.create_new_vector();
      std::shared_ptr<VectorType> residual = solver.create_new_vector();

      copy(*src_yy, yy);
      copy(*src_yp, yp);

      int err = solver.residual(tt, *src_yy, *src_yp, *residual);

      copy(rr, *residual);

      return err;
    }



    template<typename VectorType>
    int t_dae_lsetup(IDAMem IDA_mem,
                     N_Vector yy,
                     N_Vector yp,
                     N_Vector resp,
                     N_Vector tmp1,
                     N_Vector tmp2,
                     N_Vector tmp3)
    {
      (void) tmp1;
      (void) tmp2;
      (void) tmp3;
      (void) resp;
      IDAInterface<VectorType> &solver = *static_cast<IDAInterface<VectorType> *>(IDA_mem->ida_user_data);

      std::shared_ptr<VectorType> src_yy = solver.create_new_vector();
      std::shared_ptr<VectorType> src_yp = solver.create_new_vector();

      copy(*src_yy, yy);
      copy(*src_yp, yp);

      int err = solver.setup_jacobian(IDA_mem->ida_tn,
                                      *src_yy,
                                      *src_yp,
                                      IDA_mem->ida_cj);
      return err;
    }


    template<typename VectorType>
    int t_dae_solve(IDAMem IDA_mem,
                    N_Vector b,
                    N_Vector weight,
                    N_Vector yy,
                    N_Vector yp,
                    N_Vector resp)
    {
      (void) weight;
      (void) yy;
      (void) yp;
      (void) resp;
      IDAInterface<VectorType> &solver = *static_cast<IDAInterface<VectorType> *>(IDA_mem->ida_user_data);

      std::shared_ptr<VectorType> dst = solver.create_new_vector();
      std::shared_ptr<VectorType> src = solver.create_new_vector();

      copy(*src, b);

      int err = solver.solve_jacobian_system(*src,*dst);
      copy(b, *dst);
      return err;
    }

  }

  template <typename VectorType>
  IDAInterface<VectorType>::IDAInterface(const MPI_Comm mpi_comm,
                                         const double &initial_time,
                                         const double &final_time,
                                         const double &initial_step_size,
                                         const double &min_step_size,
                                         const double &abs_tol,
                                         const double &rel_tol,
                                         const unsigned int &max_order,
                                         const double &output_period,
                                         const bool &ignore_algebraic_terms_for_errors,
                                         const std::string &ic_type,
                                         const std::string &reset_type,
                                         const double &ic_alpha,
                                         const unsigned int &ic_max_iter,
                                         const unsigned int &max_non_linear_iterations,
                                         const bool &verbose,
                                         const bool &use_local_tolerances) :
    initial_time(initial_time),
    final_time(final_time),
    initial_step_size(initial_step_size),
    min_step_size(min_step_size),
    abs_tol(abs_tol),
    rel_tol(rel_tol),
    max_order(max_order),
    output_period(output_period),
    ignore_algebraic_terms_for_errors(ignore_algebraic_terms_for_errors),
    ic_type(ic_type),
    reset_type(reset_type),
    ic_alpha(ic_alpha),
    ic_max_iter(ic_max_iter),
    max_non_linear_iterations(max_non_linear_iterations),
    verbose(verbose),
    use_local_tolerances(use_local_tolerances),
    ida_mem(nullptr),
    communicator(Utilities::MPI::duplicate_communicator(mpi_comm)),
    pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_comm)==0)
  {
    set_functions_to_trigger_an_assert();
  }

  template <typename VectorType>
  IDAInterface<VectorType>::~IDAInterface()
  {
    if (ida_mem)
      IDAFree(&ida_mem);
    MPI_Comm_free(&communicator);
  }

  template <typename VectorType>
  void IDAInterface<VectorType>::add_parameters(ParameterHandler &prm)
  {
    prm.add_parameter("Initial step size",initial_step_size);

    prm.add_parameter("Minimum step size", min_step_size);

    prm.add_parameter("Absolute error tolerance", abs_tol);

    prm.add_parameter("Relative error tolerance", rel_tol);

    prm.add_parameter("Initial time", initial_time);

    prm.add_parameter("Final time", final_time);

    prm.add_parameter("Time interval between each output", output_period);

    prm.add_parameter("Maximum order of BDF", max_order);

    prm.add_parameter("Maximum number of nonlinear iterations", max_non_linear_iterations);

    prm.add_parameter("Ignore algebraic terms for error computations", ignore_algebraic_terms_for_errors,
                      "Indicate whether or not to suppress algebraic variables "
                      "in the local error test.");

    prm.add_parameter("Initial condition type", ic_type,
                      "This is one of the following three options for the "
                      "initial condition calculation. \n"
                      " none: do not try to make initial conditions consistent. \n"
                      " use_y_diff: compute the algebraic components of y and differential\n"
                      "    components of y_dot, given the differential components of y. \n"
                      "    This option requires that the user specifies differential and \n"
                      "    algebraic components in the function get_differential_components.\n"
                      " use_y_dot: compute all components of y, given y_dot.",
                      Patterns::Selection("none|use_y_diff|use_y_dot"));

    prm.add_parameter("Initial condition type after restart", reset_type,
                      "This is one of the following three options for the "
                      "initial condition calculation. \n"
                      " none: do not try to make initial conditions consistent. \n"
                      " use_y_diff: compute the algebraic components of y and differential\n"
                      "    components of y_dot, given the differential components of y. \n"
                      "    This option requires that the user specifies differential and \n"
                      "    algebraic components in the function get_differential_components.\n"
                      " use_y_dot: compute all components of y, given y_dot.",
                      Patterns::Selection("none|use_y_diff|use_y_dot"));

    prm.add_parameter("Initial condition Newton parameter", ic_alpha);

    prm.add_parameter("Initial condition Newton max iterations", ic_max_iter);

    prm.add_parameter("Use local tolerances", use_local_tolerances);

    prm.add_parameter("Show output of time steps", verbose);
  }


  template <typename VectorType>
  unsigned int IDAInterface<VectorType>::solve_dae(VectorType &solution,
                                                   VectorType &solution_dot)
  {

    unsigned int system_size = solution.size();
    unsigned int local_system_size = system_size;

    double t = initial_time;
    double h = initial_step_size;
    unsigned int step_number = 0;

    int status;

    // The solution is stored in
    // solution. Here we take only a
    // view of it.
#ifdef DEAL_II_WITH_MPI
    if (is_serial_vector<VectorType>::value == false)
      {
        IndexSet is = solution.locally_owned_elements();
        local_system_size = is.n_elements();

        yy        = N_VNew_Parallel(communicator,
                                    local_system_size,
                                    system_size);

        yp        = N_VNew_Parallel(communicator,
                                    local_system_size,
                                    system_size);

        diff_id   = N_VNew_Parallel(communicator,
                                    local_system_size,
                                    system_size);

        abs_tolls = N_VNew_Parallel(communicator,
                                    local_system_size,
                                    system_size);
      }
    else
#endif
      {
        Assert(is_serial_vector<VectorType>::value,
               ExcInternalError("Trying to use a serial code with a parallel vector."));
        yy        = N_VNew_Serial(system_size);
        yp        = N_VNew_Serial(system_size);
        diff_id   = N_VNew_Serial(system_size);
        abs_tolls = N_VNew_Serial(system_size);
      }
    reset(initial_time,
          initial_step_size,
          solution,
          solution_dot);

    double next_time = initial_time;

    output_step( 0, solution, solution_dot, 0);

    while (t<final_time)
      {

        next_time += output_period;
        if (verbose)
          {
            pcout << " "//"\r"
                  << std::setw(5) << t << " ----> "
                  << std::setw(5) << next_time
                  << std::endl;
          }
        status = IDASolve(ida_mem, next_time, &t, yy, yp, IDA_NORMAL);

        status = IDAGetLastStep(ida_mem, &h);
        AssertThrow(status == 0, ExcMessage("Error in IDA Solver"));

        copy(solution, yy);
        copy(solution_dot, yp);

        while (solver_should_restart(t, solution, solution_dot))
          reset(t, h, solution, solution_dot);

        step_number++;

        output_step(t, solution, solution_dot, step_number);
      }

    pcout << std::endl;
    // Free the vectors which are no longer used.
#ifdef DEAL_II_WITH_MPI
    if (is_serial_vector<VectorType>::value == false)
      {
        N_VDestroy_Parallel(yy);
        N_VDestroy_Parallel(yp);
        N_VDestroy_Parallel(abs_tolls);
        N_VDestroy_Parallel(diff_id);
      }
    else
#endif
      {
        N_VDestroy_Serial(yy);
        N_VDestroy_Serial(yp);
        N_VDestroy_Serial(abs_tolls);
        N_VDestroy_Serial(diff_id);
      }

    return step_number;
  }

  template <typename VectorType>
  void IDAInterface<VectorType>::reset(const double &current_time,
                                       const double &current_time_step,
                                       VectorType &solution,
                                       VectorType &solution_dot)
  {

    unsigned int system_size;
    unsigned int local_system_size;
    bool first_step = (current_time == initial_time);

    if (ida_mem)
      IDAFree(&ida_mem);

    ida_mem = IDACreate();


    // Free the vectors which are no longer used.
    if (yy)
      {
#ifdef DEAL_II_WITH_MPI
        if (is_serial_vector<VectorType>::value == false)
          {
            N_VDestroy_Parallel(yy);
            N_VDestroy_Parallel(yp);
            N_VDestroy_Parallel(abs_tolls);
            N_VDestroy_Parallel(diff_id);
          }
        else
#endif
          {
            N_VDestroy_Serial(yy);
            N_VDestroy_Serial(yp);
            N_VDestroy_Serial(abs_tolls);
            N_VDestroy_Serial(diff_id);
          }
      }

    int status;
    system_size = solution.size();
#ifdef DEAL_II_WITH_MPI
    if (is_serial_vector<VectorType>::value == false)
      {
        IndexSet is = solution.locally_owned_elements();
        local_system_size = is.n_elements();

        yy        = N_VNew_Parallel(communicator,
                                    local_system_size,
                                    system_size);

        yp        = N_VNew_Parallel(communicator,
                                    local_system_size,
                                    system_size);

        diff_id   = N_VNew_Parallel(communicator,
                                    local_system_size,
                                    system_size);

        abs_tolls = N_VNew_Parallel(communicator,
                                    local_system_size,
                                    system_size);

      }
    else
#endif
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

    if (use_local_tolerances)
      {
        copy(abs_tolls, get_local_tolerances());
        status += IDASVtolerances(ida_mem, rel_tol, abs_tolls);
      }
    else
      {
        status += IDASStolerances(ida_mem, rel_tol, abs_tol);
      }

    status = IDASetInitStep(ida_mem, current_time_step);
    AssertIDA(status);

    status = IDASetUserData(ida_mem, (void *) this);
    AssertIDA(status);

    if (ic_type == "use_y_diff" || reset_type == "use_y_diff" || ignore_algebraic_terms_for_errors)
      {
        VectorType diff_comp_vector(solution);
        diff_comp_vector = 0.0;
        auto dc = differential_components();
        for (auto i = dc.begin(); i != dc.end(); ++i)
          diff_comp_vector[*i] = 1.0;

        copy(diff_id, diff_comp_vector);
        status = IDASetId(ida_mem, diff_id);
        AssertIDA(status);
      }

    status = IDASetSuppressAlg(ida_mem, ignore_algebraic_terms_for_errors);
    AssertIDA(status);

//  status = IDASetMaxNumSteps(ida_mem, max_steps);
    status = IDASetStopTime(ida_mem, final_time);
    AssertIDA(status);

    status = IDASetMaxNonlinIters(ida_mem, max_non_linear_iterations);
    AssertIDA(status);

    // Initialize solver
    IDAMem IDA_mem;
    IDA_mem = (IDAMem) ida_mem;

    IDA_mem->ida_lsetup = t_dae_lsetup<VectorType>;
    IDA_mem->ida_lsolve = t_dae_solve<VectorType>;
    IDA_mem->ida_setupNonNull = true;

    status = IDASetMaxOrd(ida_mem, max_order);
    AssertIDA(status);

    std::string type;
    if (first_step)
      type = ic_type;
    else
      type = reset_type;

    if (verbose)
      {
        pcout << "computing consistent initial conditions with the option "
              << type
              << " please be patient."
              << std::endl;
      }

    if (type == "use_y_dot")
      {
        // (re)initialization of the vectors
        status = IDACalcIC(ida_mem, IDA_Y_INIT, current_time+current_time_step);
        AssertIDA(status);

        status = IDAGetConsistentIC(ida_mem, yy, yp);
        AssertIDA(status);

        copy(solution, yy);
        copy(solution_dot, yp);
      }
    else if (type == "use_y_diff")
      {
        status = IDACalcIC(ida_mem, IDA_YA_YDP_INIT, current_time+current_time_step);
        AssertIDA(status);

        status = IDAGetConsistentIC(ida_mem, yy, yp);
        AssertIDA(status);

        copy(solution, yy);
        copy(solution_dot, yp);
      }

    if (verbose)
      {
        pcout << "compute initial conditions: done."
              << std::endl;
      }
  }

  template<typename VectorType>
  void IDAInterface<VectorType>::set_initial_time(const double &t)
  {
    initial_time = t;
  }

  template<typename VectorType>
  void IDAInterface<VectorType>::set_functions_to_trigger_an_assert()
  {

    create_new_vector = []() ->std::unique_ptr<VectorType>
    {
      std::unique_ptr<VectorType> p;
      AssertThrow(false, ExcFunctionNotProvided("create_new_vector"));
      return p;
    };

    residual = [](const double,
                  const VectorType &,
                  const VectorType &,
                  VectorType &) ->int
    {
      int ret=0;
      AssertThrow(false,  ExcFunctionNotProvided("residual"));
      return ret;
    };

    setup_jacobian = [](const double,
                        const VectorType &,
                        const VectorType &,
                        const double) ->int
    {
      int ret=0;
      AssertThrow(false, ExcFunctionNotProvided("setup_jacobian"));
      return ret;
    };

    solve_jacobian_system = [](const VectorType &,
                               VectorType &) ->int
    {
      int ret=0;
      AssertThrow(false, ExcFunctionNotProvided("solve_jacobian_system"));
      return ret;
    };

    output_step = [](const double,
                     const VectorType &,
                     const VectorType &,
                     const unsigned int)
    {
      return;
    };

    solver_should_restart = [](const double,
                               VectorType &,
                               VectorType &) ->bool
    {
      return false;
    };

    differential_components = []() ->IndexSet
    {
      IndexSet i;
      AssertThrow(false, ExcFunctionNotProvided("differential_components"));
      return i;
    };

    get_local_tolerances = []() ->VectorType &
    {
      std::shared_ptr<VectorType> y;
      AssertThrow(false, ExcFunctionNotProvided("get_local_tolerances"));
      return *y;
    };
  }

  template class IDAInterface<Vector<double> >;
  template class IDAInterface<BlockVector<double> >;

#ifdef DEAL_II_WITH_MPI

#ifdef DEAL_II_WITH_TRILINOS
  template class IDAInterface<TrilinosWrappers::MPI::Vector>;
  template class IDAInterface<TrilinosWrappers::MPI::BlockVector>;
#endif

#ifdef DEAL_II_WITH_PETSC
  template class IDAInterface<PETScWrappers::MPI::Vector>;
  template class IDAInterface<PETScWrappers::MPI::BlockVector>;
#endif

#endif

}

DEAL_II_NAMESPACE_CLOSE

#endif
