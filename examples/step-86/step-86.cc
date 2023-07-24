/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2013 - 2021 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Wolfgang Bangerth, Colorado State University, 2023
 */


#include <deal.II/base/utilities.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/timer.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/petsc_ts.h>

#include <fstream>
#include <iostream>


namespace Step86
{
  using namespace dealii;

  // We derive from ParameterAcceptor to use the ParameterAcceptor facilities in
  // order to parse and read options from a parameter file (see also step-70).
  template <int dim>
  class HeatEquation : public ParameterAcceptor
  {
  public:
    HeatEquation(const MPI_Comm mpi_communicator = MPI_COMM_WORLD);
    void run();

  private:
    // This is the usual function we need to setup dofs, constraints, and linear
    // algebra objects.
    void setup_system();

    // Here is where we start to deviate from the "manual" way of solving time
    // dependent problems. High level packages for the solution of ODEs and IDAs
    // usually expect the user to provide a function that computes the residual
    // w.r.t. the solution and the time derivative of the solution.
    //
    // This allows those packages to abstract away the details of how the time
    // derivative is defined (e.g., backward Euler, Crank-Nicolson, etc.) and
    // allow the user to provide a uniform interface, irrespective of the time
    // stepper used. PETSc and SUNDIALS are two examples of such packages, and
    // they require the user to provide C style call-backs to compute residuals,
    // Jacobians, etc. In deal.II, we wrap these C stype libraries with our own
    // wrappers, that use a style which is closer to c++, and that allows us to
    // simply define the callbacks via lambda functions.
    //
    // To make it clear what we are doing here, we start by defining functions
    // with the same name of the interface that we will use later on. These are
    // the function that will be called by the PETSc time stepper, and consist
    // in the residual function:
    void implicit_function(const double                      time,
                           const PETScWrappers::MPI::Vector &solution,
                           const PETScWrappers::MPI::Vector &solution_dot,
                           PETScWrappers::MPI::Vector &      residual) const;

    // the Jacobian function
    void
    assemble_implicit_jacobian(const double                      time,
                               const PETScWrappers::MPI::Vector &solution,
                               const PETScWrappers::MPI::Vector &solution_dot,
                               const double                      shift);

    // the function that actually solves the linear system
    void solve_with_jacobian(const PETScWrappers::MPI::Vector &src,
                             PETScWrappers::MPI::Vector &      residual) const;

    // and the function that writes the solution to file at each time step
    void output_results(const double                      time,
                        const PETScWrappers::MPI::Vector &solution,
                        const unsigned int timestep_number) const;

    // in addition to the above functions, we also choose to implement an extra
    // call back, which is used to distribute affine constraints to the solution
    // vector after each solve function. We could inserted a call to
    // `distribute` inside the solve_with_jacobian function, but by doing things
    // this way, we can also replace all-together our implementation of the
    // solve_with_jacobian function with a call to a PETSc solver (from the
    // command line), and we would still have the constraints correctly applied
    // to the solution vector.
    void distribute(const double                time,
                    PETScWrappers::MPI::Vector &residual) const;

    // An extremely important part of solving a time dependent PDE system is
    // related to how we treat boundary conditions and hanging nodes
    // constraints. In particular, these constraints are not part of the ODE
    // system, and are technically not true degrees of freedom. We need to
    // instruct the time stepper to ignore these components when computing the
    // error in the residuals, and we do so by providing a function that returns
    // an IndexSet with the indices of the algebraic components of the solution
    // vector. This is the function that does that.
    IndexSet algebraic_components() const;

    // In this tutorial program, similar to what we did in step-26, we want to
    // adapt the mesh at regular time intervals. However, if we are using an
    // external time stepper, we need to make sure that the time stepper is
    // aware of the mesh changes. In particular, the time stepper must support
    // the fact that the number of degrees of freedom may change, and it must
    // support the fact all internal vectors (e.g., the solution vector and all
    // intermediate stages used to compute the current time derivative) will
    // need to be transferred to the new mesh. Deal.II does not allow to do this
    // on each vector separately, so the time stepper need to be made aware of
    // when this is happening. This function will be called by the time stepper
    // whenever it is time to adapt the mesh,
    void
    prepare_for_coarsening_and_refinement(const PETScWrappers::MPI::Vector &y);

    // and as soon as the mesh has been adapted, we will need to interpolate the
    // solution and all intermediate stages to the new mesh. This is the
    // function that does that.
    void interpolate(const std::vector<PETScWrappers::MPI::Vector> &all_in,
                     std::vector<PETScWrappers::MPI::Vector> &      all_out);

    // For a time dependent problem, boundary conditions may generally depend on
    // time. For distributed solutions, we usually use AffineConstraint objects
    // for such scope, and we therefore need to update the constraints whenever
    // time changes.
    void update_constraints(const double time) const;

    // These are the standard member functions for a parallel computation (see,
    // for example, step-40).
    const MPI_Comm mpi_communicator;

    parallel::distributed::Triangulation<dim> triangulation;
    FE_Q<dim>                                 fe;
    DoFHandler<dim>                           dof_handler;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    // The only major difference between this program and step-26 is that we now
    // apply boundary conditions with AffineConstraints, and we allow time
    // dependent boundary conditions. During the computation of the residual, we
    // will need to impose two different types of boundary conditions:
    // homogeneous boundary conditions for the time derivative of the solution,
    // and non-homogeneous boundary conditions for the solution itself. Since
    // the computation of hanging node constraints may be expensive, especially
    // for high order and distributed problems, we actually create three
    // different objects: one for hanging nodes, one for homogeneous boundary
    // and one for non-homogeneous boundary conditions. The first two are
    // created only once, while the third is updated at each time step, starting
    // from a copy of the first one.
    AffineConstraints<double>         hanging_nodes_constraints;
    AffineConstraints<double>         homogeneous_constraints;
    mutable AffineConstraints<double> constraints;

    // Instead of assembling two matrices once, and then modifying them at each
    // time step, we let the time-stepper decide when to assemble the Jacobian
    // matrix, and provide a function that does so. This is the matrix that we
    // will use to store the Jacobian matrix.
    PETScWrappers::MPI::SparseMatrix jacobian_matrix;

    // Unlike in step-26, we do not need to store ourselves the solution vector
    // at previous time steps. Instead, we will need to access the solution and
    // the time derivative of the solution at the current time step, and we
    // therefore need to provide a way to access locally relevant degrees of
    // freedom during the assembly and during the output of the solution. We do
    // so by creating two vectors that will be used to store the locally
    // relevant solution and the locally relevant time derivative of the
    // solution.
    PETScWrappers::MPI::Vector         solution;
    mutable PETScWrappers::MPI::Vector locally_relevant_solution;
    mutable PETScWrappers::MPI::Vector locally_relevant_solution_dot;

    // Finally, this is the object that we will use to store the parameters of
    // the time stepper:
    PETScWrappers::TimeStepperData time_stepper_data;

    // Since we already derive our class from ParameterAcceptor, we exploit its
    // facilities to parse also the parameters of the functions that define the
    // initial value, the right hand side, and the boundary values. We therefore
    // create three ParameterAcceptorProxy objects, which wrap the actual
    // ParsedFunction class into objects that ParameterAcceptor can handle.
    ParameterAcceptorProxy<Functions::ParsedFunction<dim>>
      initial_value_function;

    mutable ParameterAcceptorProxy<Functions::ParsedFunction<dim>>
      right_hand_side_function;

    mutable ParameterAcceptorProxy<Functions::ParsedFunction<dim>>
      boundary_values_function;

    // Finally, we also store the parameters that control the mesh refinement,
    // and the adaptive mesh refinement
    unsigned int initial_global_refinement;
    unsigned int max_delta_refinement_level;
    unsigned int adaption_frequency;

    // As in step-40, we keep track of the time spent in each function, and we
    // make sure that output is only written when we are on processor zero.
    ConditionalOStream  pcout;
    mutable TimerOutput computing_timer;
  };


  template <int dim>
  HeatEquation<dim>::HeatEquation(const MPI_Comm mpi_communicator)
    : ParameterAcceptor("/Heat Equation/")
    , mpi_communicator(mpi_communicator)
    , triangulation(mpi_communicator,
                    typename Triangulation<dim>::MeshSmoothing(
                      Triangulation<dim>::smoothing_on_refinement |
                      Triangulation<dim>::smoothing_on_coarsening))
    , fe(1)
    , dof_handler(triangulation)
    , time_stepper_data("", "beuler", 0.0, 1.0, 0.025)
    , initial_value_function("/Heat Equation/Initial value")
    , right_hand_side_function("/Heat Equation/Right hand side")
    , boundary_values_function("/Heat Equation/Boundary values")
    , initial_global_refinement(5)
    , max_delta_refinement_level(2)
    , adaption_frequency(0)
    , pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
    , computing_timer(mpi_communicator,
                      pcout,
                      TimerOutput::summary,
                      TimerOutput::wall_times)
  {
    enter_subsection("Time stepper");
    enter_my_subsection(this->prm);
    time_stepper_data.add_parameters(this->prm);
    leave_my_subsection(this->prm);
    leave_subsection();

    add_parameter("Initial global refinement",
                  initial_global_refinement,
                  "Number of times the mesh is refined globally before "
                  "starting the time stepping.");
    add_parameter("Maximum delta refinement level",
                  max_delta_refinement_level,
                  "Maximum number of local refinement levels.");
    add_parameter("Adaption frequency",
                  adaption_frequency,
                  "When to adapt the mesh.");
  }



  template <int dim>
  void HeatEquation<dim>::setup_system()
  {
    TimerOutput::Scope t(computing_timer, "setup_system");
    dof_handler.distribute_dofs(fe);

    pcout << std::endl
          << "===========================================" << std::endl
          << "Number of active cells: " << triangulation.n_active_cells()
          << std::endl
          << "Number of degrees of freedom: " << dof_handler.n_dofs()
          << std::endl
          << std::endl;

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);

    // We first create hanging node constraints. These affine constraints won't
    // be used as standalone constraints, but as a starting point for the
    // non-homogeneous constraints.
    hanging_nodes_constraints.clear();
    hanging_nodes_constraints.reinit(locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler,
                                            hanging_nodes_constraints);
    hanging_nodes_constraints.make_consistent_in_parallel(locally_owned_dofs,
                                                          locally_relevant_dofs,
                                                          mpi_communicator);
    hanging_nodes_constraints.close();

    // We then create an object used for homogeneous constraints. In time
    // dependent problems, these are needed since we impose the constraints
    // through algebraic equations. While technically it would be possible to
    // use the time derivative of the boundary function as a boundary conditions
    // for the time derivative of the solution, this is not done here. Instead,
    // we impose the boundary conditions through algebraic equations, and
    // therefore the time derivative of the boundary conditions is not part of
    // the algebraic system, and we need zero boundary conditions on the time
    // derivative of the solution when computing the residual. We use the
    // homogeneous_constraints object for this purpose.
    homogeneous_constraints.clear();
    homogeneous_constraints.reinit(locally_relevant_dofs);
    homogeneous_constraints.merge(hanging_nodes_constraints);
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(),
                                             homogeneous_constraints);
    homogeneous_constraints.make_consistent_in_parallel(locally_owned_dofs,
                                                        locally_relevant_dofs,
                                                        mpi_communicator);
    homogeneous_constraints.close();

    // We then create the actual non-homogeneous constraints. These are used
    // during the assembly and during the residual evaluation. Since this is
    // something we need to do at every time step, we created a member function
    // for it.
    update_constraints(boundary_values_function.get_time());

    DynamicSparsityPattern dsp(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    homogeneous_constraints,
                                    false);
    SparsityTools::distribute_sparsity_pattern(dsp,
                                               locally_owned_dofs,
                                               mpi_communicator,
                                               locally_relevant_dofs);

    // directly initialize from dsp, no need for the regular sparsity pattern:
    jacobian_matrix.reinit(locally_owned_dofs,
                           locally_owned_dofs,
                           dsp,
                           mpi_communicator);

    // finally, we initialize the solution and the two locally_relevant vectors
    solution.reinit(locally_owned_dofs, mpi_communicator);
    locally_relevant_solution.reinit(locally_owned_dofs,
                                     locally_relevant_dofs,
                                     mpi_communicator);
    locally_relevant_solution_dot.reinit(locally_owned_dofs,
                                         locally_relevant_dofs,
                                         mpi_communicator);
  }


  // Since regenerating the constraints at each time step may be expensive, we
  // make sure that we only do so when the time changes. We track time change by
  // checking if the time of the boundary_values_function has changed, with
  // respect to the time of the last call to this function. This will work most
  // of the times, but not the very first time we call this function, since the
  // time then may be zero and the time of the boundary_values_function is zero
  // at construction time. We therefore also check if the number of constraints,
  // and if these are empty, we regenerate the constraints regardless of the
  // time variable.
  template <int dim>
  void HeatEquation<dim>::update_constraints(const double time) const
  {
    if (constraints.n_constraints() == 0 ||
        time != boundary_values_function.get_time())
      {
        TimerOutput::Scope t(computing_timer, "update_constraints");
        boundary_values_function.set_time(time);
        constraints.clear();
        constraints.reinit(locally_relevant_dofs);
        constraints.merge(hanging_nodes_constraints);
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 0,
                                                 boundary_values_function,
                                                 constraints);
        constraints.make_consistent_in_parallel(locally_owned_dofs,
                                                locally_relevant_dofs,
                                                mpi_communicator);
        constraints.close();
      }
  }


  template <int dim>
  void HeatEquation<dim>::implicit_function(
    const double                      time,
    const PETScWrappers::MPI::Vector &y,
    const PETScWrappers::MPI::Vector &y_dot,
    PETScWrappers::MPI::Vector &      residual) const
  {
    TimerOutput::Scope t(computing_timer, "implicit_function");
    right_hand_side_function.set_time(time);

    // We face two difficulties here: the first is that the y and y_dot vectors
    // are read only, and we need to make sure they satisfy the correct boundary
    // conditions, and the second is that we need to compute the residual, and
    // therefore in general we need to evaluate their values and gradients
    // inside locally owned cells, and we need access to degrees of freedom
    // which may be owned by neighboring processors. We therefore need to create
    // a (non-ghosted) copy of the vectors, apply boundary conditions and
    // hanging node constraints, and then copy to ghosted vectors before we can
    // do anything sensible with them.
    PETScWrappers::MPI::Vector tmp_solution(y);
    PETScWrappers::MPI::Vector tmp_solution_dot(y_dot);

    // Fix boundary conditions
    update_constraints(time);
    constraints.distribute(tmp_solution);
    homogeneous_constraints.distribute(tmp_solution_dot);

    // Copy to ghosted vectors
    locally_relevant_solution     = tmp_solution;
    locally_relevant_solution_dot = tmp_solution_dot;

    // Anything that follows here is standard, and has been done several times
    // already, from step-3, step-4, step-6, and step-40
    QGauss<dim>   quadrature_formula(fe.degree + 1);
    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    std::vector<Tensor<1, dim>> solution_gradients(n_q_points);
    std::vector<double>         solution_dot_values(n_q_points);

    Vector<double> cell_residual(dofs_per_cell);

    residual = 0;
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);

          // as we said before, we need to evaluate the solution and its
          // gradient inside locally owned cells, and for this, we need also
          // access to degrees of freedom that may be owned by neighboring
          // processors. We therefore use the locally_relevant_solution and and
          // locally_relevant_solution_dot vectors.
          fe_values.get_function_gradients(locally_relevant_solution,
                                           solution_gradients);
          fe_values.get_function_values(locally_relevant_solution_dot,
                                        solution_dot_values);

          cell->get_dof_indices(local_dof_indices);

          cell_residual = 0;
          for (const unsigned int q : fe_values.quadrature_point_indices())
            for (const unsigned int i : fe_values.dof_indices())
              {
                cell_residual(i) +=
                  (fe_values.shape_value(i, q) * solution_dot_values[q] +
                   fe_values.shape_grad(i, q) * solution_gradients[q] -
                   fe_values.shape_value(i, q) *
                     right_hand_side_function.value(
                       fe_values.quadrature_point(q))) *
                  fe_values.JxW(q);
              }
          constraints.distribute_local_to_global(cell_residual,
                                                 local_dof_indices,
                                                 residual);
        }
    residual.compress(VectorOperation::add);

    // Now we correct the entries corresponding to constrained degrees of
    // freedom. From the point of view of the residual vector, if the input y
    // vector does not contain the correct values on constrained degrees of
    // freedom (hanging nodes or boundary conditions), we need to communicate
    // this to the time stepper, and we do so by setting the residual to the
    // actual difference between the input y vector and the our local copy of
    // it, in which we have applied the constraints. Since we have made a copy
    // of the input vector for this purpose, we use it to compute the residual
    // value. However, there is a difference between hanging nodes constraints
    // and boundary conditions: we do not want to make hanging node constraints
    // actually depend on their dependent degrees of freedom, since this would
    // imply that we are actually solving for the dependent degrees of freedom.
    // This is not what we are actually doing, however, since hanging nodes are
    // not actually solved for. They are eliminated from the system by the call
    // to AffineConstraints::distribute_local_to_global() above. From the point
    // of view of the Jacobian matrix, we are effectively setting hanging nodes
    // to an artificial value (usually zero), and we simply want to make sure
    // that we solve for those degrees of freedom a posteriori, by calling the
    // function AffineConstraints::distribute().
    //
    // Here we therefore check that the residual is equal to the input value on
    // the constrained dofs corresponding to hanging nodes (i.e., those for
    // which the lines of the constraints contain at least one other entry), and
    // to the difference between the input vector and the actual solution on
    // those constraints that correspond to boundary conditions. 
    for (const auto &c : constraints.get_lines())
      if (locally_owned_dofs.is_element(c.index))
        {
          if (c.entries.empty())
            {
              // This is an essential boundary condition
              residual[c.index] = y[c.index] - tmp_solution[c.index];
            }
          else
            {
              // This is a hanging node.
              residual[c.index] = y[c.index];
            }
        }
    residual.compress(VectorOperation::insert);
  }


  template <int dim>
  void HeatEquation<dim>::assemble_implicit_jacobian(
    const double,
    const PETScWrappers::MPI::Vector &,
    const PETScWrappers::MPI::Vector &,
    const double shift)
  {
    TimerOutput::Scope t(computing_timer, "assemble_implicit_jacobian");
    QGauss<dim>        quadrature_formula(fe.degree + 1);
    FEValues<dim>      fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

    jacobian_matrix = 0;
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);

          cell->get_dof_indices(local_dof_indices);

          cell_matrix = 0;
          for (const unsigned int q : fe_values.quadrature_point_indices())
            for (const unsigned int i : fe_values.dof_indices())
              for (const unsigned int j : fe_values.dof_indices())
                {
                  cell_matrix(i, j) +=
                    (shift * fe_values.shape_value(i, q) *
                       fe_values.shape_value(j, q) +
                     fe_values.shape_grad(i, q) * fe_values.shape_grad(j, q)) *
                    fe_values.JxW(q);
                }
          constraints.distribute_local_to_global(cell_matrix,
                                                 local_dof_indices,
                                                 jacobian_matrix);
        }
    jacobian_matrix.compress(VectorOperation::add);
    // For any constrained degree of freedom, we set the diagonal of the
    // jacobian to be one. This makes the jacobian matrix invertible, consistent
    // with what the time stepper expects, and it also makes sure that if we did
    // not make a mistake in the residual and/or in the jacbian matrix, then
    // asking the time stepper to check the jacobian with a finite difference
    // method will produce the correct result. This can be activated at run time
    // via passing the `-snes_test_jacobian` option on the command line.
    for (const auto &c : constraints.get_lines())
      jacobian_matrix.set(c.index, c.index, 1.0);
    jacobian_matrix.compress(VectorOperation::insert);
  }

  // Whenever the time stepper solves a linear system, it calls this function to
  // make sure that hanging nodes and boundary conditions are correctly applied.
  // We could have done this inside the solve_with_jacobian function, but by not
  // doing so, we can also replace the solve_with_jacobian function with a call
  // to a PETSc solver, and we would still have the constraints correctly
  // applied to the solution vector.
  template <int dim>
  void HeatEquation<dim>::distribute(const double                time,
                                     PETScWrappers::MPI::Vector &residual) const
  {
    TimerOutput::Scope t(computing_timer, "distribute");
    update_constraints(time);
    constraints.distribute(residual);
  }

  // This is the function that actually solves the linear system. We could in
  // principle not provide this function to the time stepper, and instead select
  // a specific solver on the command line by using the `-ksp_*` options of
  // petsc. However, by providing this function, we can use a specific solver
  // and preconditioner for the linear system, and still have the possibility to
  // change them on the command line.
  //
  // Providing a specific solver is more inline with the way we usually do
  // things in other deal.II examples, while letting PETSc choose a generic
  // solver, and changing it on the command line via `-ksp_type` is more inline
  // with the way PETSc is usually used. Both options are available here, since
  // we can still change both the solver and the preconditioner on the command
  // line via `-user_ksp_type` and `-user_pc_type` options.
  template <int dim>
  void
  HeatEquation<dim>::solve_with_jacobian(const PETScWrappers::MPI::Vector &src,
                                         PETScWrappers::MPI::Vector &dst) const
  {
    TimerOutput::Scope      t(computing_timer, "solve_with_jacobian");
    SolverControl           solver_control(1000, 1e-8 * src.l2_norm());
    PETScWrappers::SolverCG cg(solver_control);
    cg.set_prefix("user_");

#if defined(PETSC_HAVE_HYPRE)
    PETScWrappers::PreconditionBoomerAMG preconditioner;
    preconditioner.initialize(jacobian_matrix);
#else
    PETScWrappers::PreconditionSSOR preconditioner;
    preconditioner.initialize(
      jacobian_matrix, PETScWrappers::PreconditionSSOR::AdditionalData(1.0));
#endif
    cg.solve(jacobian_matrix, dst, src, preconditioner);

    pcout << "     " << solver_control.last_step() << " linear iterations."
          << std::endl;
  }


  // This function is called by the time stepper whenever it deems appropriate
  // to monitor the output solution. We use it to write the solution to a file,
  // and provide graphical output through paraview or visit. We also write a pvd
  // file, which groups all meta informations about the vtu files into a single
  // that can be used to load the full time dependent solution in paraview.
  template <int dim>
  void
  HeatEquation<dim>::output_results(const double                      time,
                                    const PETScWrappers::MPI::Vector &y,
                                    const unsigned int timestep_number) const
  {
    TimerOutput::Scope t(computing_timer, "output_results");
    DataOut<dim>       data_out;

    locally_relevant_solution = y;

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(locally_relevant_solution, "U");

    data_out.build_patches();

    data_out.set_flags(DataOutBase::VtkFlags(time, timestep_number));

    const std::string filename =
      "solution-" + Utilities::int_to_string(timestep_number, 3) + ".vtu";
    data_out.write_vtu_in_parallel(filename, mpi_communicator);

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        static std::vector<std::pair<double, std::string>> times_and_names;
        times_and_names.emplace_back(time, filename);
        std::ofstream pvd_output("solution.pvd");
        DataOutBase::write_pvd_record(pvd_output, times_and_names);
      }
  }


  // This function is essentially identical to step-26, with the only difference
  // that it uses the parallel::distributed::GridRefinement namespace function
  // instead of the serial one. Once again, we make sure that we never fall
  // below the minimum refinement level, and above the maximum one, that we can
  // select from the parameter file.
  template <int dim>
  void HeatEquation<dim>::prepare_for_coarsening_and_refinement(
    const PETScWrappers::MPI::Vector &sol)
  {
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

    locally_relevant_solution = sol;
    KellyErrorEstimator<dim>::estimate(dof_handler,
                                       QGauss<dim - 1>(fe.degree + 1),
                                       {},
                                       locally_relevant_solution,
                                       estimated_error_per_cell);

    parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction(
      triangulation, estimated_error_per_cell, 0.6, 0.4);

    const auto max_grid_level =
      initial_global_refinement + max_delta_refinement_level;
    const auto min_grid_level = initial_global_refinement;

    if (triangulation.n_levels() > max_grid_level)
      for (const auto &cell :
           triangulation.active_cell_iterators_on_level(max_grid_level))
        cell->clear_refine_flag();
    for (const auto &cell :
         triangulation.active_cell_iterators_on_level(min_grid_level))
      cell->clear_coarsen_flag();
  }


  // This function is called by the time stepper whenever it requires to
  // transfer the solution and any intermediate stage vectors to a new mesh. We
  // must make sure that all input vectors are transformed into ghosted vectors
  // before the actual transfer is executed, and that we distribute the hanging
  // node constraints on the output vectors as soon as we have interpolated the
  // vectors to the new mesh.
  //
  // We have no way to enforce boundary conditions at this stage, since every
  // type of time advancing scheme may have computed the intermediate stage at
  // different times, and we have no way to know what boundary conditions to use
  // at what stage.
  //
  // While this could be a problem if we used the values of the solution and of
  // the intermediate stages on the constrained degrees of freedom to compute
  // the errors, we do not do so. Instead, we compute the errors on the
  // differential equation, and ignore algebraic constraints, therefore we do no
  // need to guarantee that the boundary conditions are satisfied also in the
  // intermediate stages.
  //
  // We have at our disposal the hanging node constraints alone, though, and
  // therefore we enforce them on the output vectors, even if this is not really
  // needed.
  template <int dim>
  void HeatEquation<dim>::interpolate(
    const std::vector<PETScWrappers::MPI::Vector> &all_in,
    std::vector<PETScWrappers::MPI::Vector> &      all_out)
  {
    parallel::distributed::SolutionTransfer<dim, PETScWrappers::MPI::Vector>
      solution_trans(dof_handler);

    std::vector<PETScWrappers::MPI::Vector> all_in_ghosted(all_in.size());
    std::vector<const PETScWrappers::MPI::Vector *> all_in_ghosted_ptr(
      all_in.size());
    std::vector<PETScWrappers::MPI::Vector *> all_out_ptr(all_in.size());
    for (unsigned int i = 0; i < all_in.size(); ++i)
      {
        all_in_ghosted[i].reinit(locally_owned_dofs,
                                 locally_relevant_dofs,
                                 mpi_communicator);
        all_in_ghosted[i]     = all_in[i];
        all_in_ghosted_ptr[i] = &all_in_ghosted[i];
      }

    triangulation.prepare_coarsening_and_refinement();
    solution_trans.prepare_for_coarsening_and_refinement(all_in_ghosted_ptr);
    triangulation.execute_coarsening_and_refinement();

    setup_system();

    all_out.resize(all_in.size());
    for (unsigned int i = 0; i < all_in.size(); ++i)
      {
        all_out[i].reinit(locally_owned_dofs, mpi_communicator);
        all_out_ptr[i] = &all_out[i];
      }
    solution_trans.interpolate(all_out_ptr);

    // Ignore boundary conditions, since we have no idea at what time the
    // functions that have been passed here were generated. This should not be
    // an issue, since Dirichlet boundary conditions are treated separately in
    // the residual
    for (auto &v : all_out)
      hanging_nodes_constraints.distribute(v);
  }


  // Both boundary conditions and hanging node constraints are indicated as
  // algebraic constraints to the time stepper. This function returns a set of
  // indices that indicate which degrees of freedom are algebraically
  // constrained. This is used by the time stepper to determine which degrees of
  // freedom should be ignored when computing the errors.
  template <int dim>
  IndexSet HeatEquation<dim>::algebraic_components() const
  {
    auto algebraic_set = DoFTools::extract_hanging_node_dofs(dof_handler);
    algebraic_set.add_indices(DoFTools::extract_boundary_dofs(dof_handler));
    return algebraic_set;
  }


  // This is the main function of the class. It sets up the system, instruct the
  // time stepper which callbacks to use, and then calls the time stepper to
  // solve the problem.
  template <int dim>
  void HeatEquation<dim>::run()
  {
    GridGenerator::hyper_L(triangulation);
    triangulation.refine_global(initial_global_refinement);

    setup_system();

    VectorTools::interpolate(dof_handler, initial_value_function, solution);

    PETScWrappers::TimeStepper<PETScWrappers::MPI::Vector,
                               PETScWrappers::MPI::SparseMatrix>
      petsc_ts(time_stepper_data);

    petsc_ts.set_matrices(jacobian_matrix, jacobian_matrix);

    petsc_ts.implicit_function =
      [&](const auto t, const auto &y, const auto &y_dot, auto &res) {
        this->implicit_function(t, y, y_dot, res);
      };

    petsc_ts.setup_jacobian =
      [&](const auto t, const auto &y, const auto &y_dot, const auto alpha) {
        this->assemble_implicit_jacobian(t, y, y_dot, alpha);
      };

    petsc_ts.solve_with_jacobian = [&](const auto &src, auto &dst) {
      this->solve_with_jacobian(src, dst);
    };

    petsc_ts.distribute = [&](const auto t, auto &y) {
      this->distribute(t, y);
    };

    petsc_ts.algebraic_components = [&]() {
      return this->algebraic_components();
    };

    petsc_ts.monitor =
      [&](const auto t, const auto &y, const auto step_number) {
        pcout << "Time step " << step_number << " at t=" << t << std::endl;
        this->output_results(t, y, step_number);
      };

    petsc_ts.prepare_for_coarsening_and_refinement =
      [&](const auto, const auto it, const auto &y, auto &resize) {
        if (it > 0 && this->adaption_frequency > 0 &&
            it % this->adaption_frequency == 0)
          {
            this->prepare_for_coarsening_and_refinement(y);
            resize = true;
          }
      };

    petsc_ts.interpolate = [&](const auto &all_in, auto &all_out) {
      this->interpolate(all_in, all_out);
    };

    petsc_ts.solve(solution);
  }
} // namespace Step86



int main(int argc, char **argv)
{
  try
    {
      using namespace Step86;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      HeatEquation<2>                  heat_equation_solver;
      ParameterAcceptor::initialize("heat_equation.prm",
                                    "heat_equation_used.prm");
      heat_equation_solver.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
