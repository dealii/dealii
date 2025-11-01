/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2000 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 *
 * Authors:
 *   Wolfgang Bangerth, Colorado State University, 2024
 *   Stefano Zampini, King Abdullah University of Science and Technology, 2024
 */


// The program starts with the usual prerequisites, all of which you will know
// by now from either step-26 (the heat equation solver) or step-40 (the
// parallel Laplace solver):
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

// The only new include file relevant here is the one that provides us with the
// PETScWrappers::TimerStepper class:
#include <deal.II/lac/petsc_ts.h>

#include <fstream>
#include <iostream>


namespace Step86
{
  using namespace dealii;

  // @sect3{The HeatEquation class}
  //
  // At its core, this program's principal structure can be understood quite
  // easily if you know step-26 (for the heat equation) and step-40 (for how
  // a parallel solver looks like). It has many of the usual member functions
  // and member variables that for convenience of documentation we will list
  // towards the top of the class so that we can document the remainder
  // separately below.
  //
  // We derive the main class from ParameterAcceptor to make dealing with run
  // time parameters and reading them a parameter file easier. step-60 and
  // step-70 have already explained how this works.
  template <int dim>
  class HeatEquation : public ParameterAcceptor
  {
  public:
    HeatEquation(const MPI_Comm mpi_communicator);
    void run();

  private:
    const MPI_Comm mpi_communicator;

    ConditionalOStream pcout;
    TimerOutput        computing_timer;

    parallel::distributed::Triangulation<dim> triangulation;
    FE_Q<dim>                                 fe;
    DoFHandler<dim>                           dof_handler;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;


    void setup_system(const double time);

    void output_results(const double                      time,
                        const unsigned int                timestep_number,
                        const PETScWrappers::MPI::Vector &solution);

    // At this point, we start to deviate from the "manual" way of solving time
    // dependent problems. High level packages for the solution of ODEs
    // usually expect the user to provide a function that computes the residual
    // of the equation and the time derivative of the solution.
    //
    // This allows those packages to abstract away the details of how the time
    // derivative is defined (e.g., backward Euler, Crank-Nicolson, etc.) and
    // allow the user to provide a uniform interface, irrespective of the time
    // stepper used. PETSc TS and SUNDIALS are two examples of such packages,
    // and they require the user to provide C style call-backs to compute
    // residuals, Jacobians, etc. In deal.II, we wrap these C style libraries
    // with our own wrappers that use a style closer to c++, and that allows us
    // to simply define the callbacks via lambda functions. Several of these
    // lambda functions will simply call member functions to do the actual work.
    //
    // To make it clear what we are doing here, we start by defining functions
    // with the same name of the interface that we will use later on (and
    // documented both in the introduction of this program as well in the
    // documentation of the PETScWrappers::TimeStepper class). These are
    // the functions that compute the right hand side, the "Jacobian matrix"
    // (see the introduction for how that is defined), and that can solve
    // a linear system with the Jacobian matrix. At the bottom of the following
    // block, we also declare the matrix object that will store this Jacobian.
    //
    // Note that all of these functions receive the current solution (and,
    // where necessary, its time derivative) as inputs. This is because *the*
    // solution vector -- i.e., the variable that stores the current state
    // of the solution -- is kept inside the time integrator object. This
    // is useful: If we kept a copy of the solution vector as a member variable
    // of the current class, one would continuously have to wonder whether it
    // is still in sync with the version of the solution vector the time
    // integrator object internally believes is the currently correct
    // version of this vector. As a consequence, we do not store such a
    // copy here: Whenever a function requires access to the current value
    // of the solution, it receives it from the time integrator as a const
    // argument. The same observation can be made about the variable that stores
    // the current time, or the current length of the time step, or the number
    // of time steps performed so far: They are all kept inside the time
    // stepping object.
    void implicit_function(const double                      time,
                           const PETScWrappers::MPI::Vector &solution,
                           const PETScWrappers::MPI::Vector &solution_dot,
                           PETScWrappers::MPI::Vector       &residual);

    void
    assemble_implicit_jacobian(const double                      time,
                               const PETScWrappers::MPI::Vector &solution,
                               const PETScWrappers::MPI::Vector &solution_dot,
                               const double                      shift);

    void solve_with_jacobian(const PETScWrappers::MPI::Vector &src,
                             PETScWrappers::MPI::Vector       &residual);

    PETScWrappers::MPI::SparseMatrix jacobian_matrix;


    // In this tutorial program, similar to what we did in step-26, we
    // want to adapt the mesh at regular time intervals. However, if
    // we are using an external time stepper, we need to make sure
    // that the time stepper is aware of the mesh changes -- see the
    // discussion on this topic in the introduction. In particular,
    // the time stepper must support the fact that the number of
    // degrees of freedom may change, and it must support the fact
    // that all internal vectors (e.g., the solution vector and all
    // intermediate stages used to compute the current time
    // derivative) will need to be transferred to the new mesh. For
    // reasons that will be discussed in the implementation below, we
    // split the mesh refinement operation into two functions: The
    // first will mark which cells to refine, based on a given
    // solution vector from which we can compute error indicators; the
    // second one will do the actual mesh refinement and transfer a
    // set of vectors from the old to the new mesh.
    void prepare_for_coarsening_and_refinement(
      const PETScWrappers::MPI::Vector &solution);

    void transfer_solution_vectors_to_new_mesh(
      const double                                   time,
      const std::vector<PETScWrappers::MPI::Vector> &all_in,
      std::vector<PETScWrappers::MPI::Vector>       &all_out);

    // As also discussed in the introduction, we also have to deal with
    // "algebraic" solution components, i.e., degrees of freedom that are not
    // really free but instead have values determined either by boundary
    // conditions or by hanging node constraints. While the values of the
    // boundary conditions can change from time step to time step (because the
    // function $g(\mathbf x,t)$ may indeed depend on time), hanging node
    // constraints remain the same as long as the mesh remains the same. As a
    // consequence, we will keep an AffineConstraints object that stores the
    // hanging node constraints and that is only updated when the mesh changes,
    // and then an AffineConstraints object `current_constraints` that we will
    // initialize with the hanging node constraints and then add the constraints
    // due to Dirichlet boundary values at a specified time. This evaluation of
    // boundary values and combining of constraints happens in the
    // `update_current_constraints()` function.
    //
    // At one place in the program, we will also need an object that constrains
    // the same degrees of freedom, but with zero values even if the boundary
    // values for the solution are non-zero; we will keep this modified set of
    // constraints in `homogeneous_constraints`.
    AffineConstraints<double> hanging_node_constraints;
    AffineConstraints<double> current_constraints;
    AffineConstraints<double> homogeneous_constraints;

    void update_current_constraints(const double time);

    // The remainder of the class is simply consumed with objects that
    // describe either the functioning of the time stepper, when and how to
    // do mesh refinement, and the objects that describe right hand side,
    // initial conditions, and boundary conditions for the PDE.
    // Since we already derive our class from ParameterAcceptor, we exploit its
    // facilities to parse also the parameters of the functions that define the
    // initial value, the right hand side, and the boundary values. We therefore
    // create three ParameterAcceptorProxy objects, which wrap the actual
    // ParsedFunction class into objects that ParameterAcceptor can handle.
    PETScWrappers::TimeStepperData time_stepper_data;

    unsigned int initial_global_refinement;
    unsigned int max_delta_refinement_level;
    unsigned int mesh_adaptation_frequency;

    ParameterAcceptorProxy<Functions::ParsedFunction<dim>>
      right_hand_side_function;
    ParameterAcceptorProxy<Functions::ParsedFunction<dim>>
      initial_value_function;
    ParameterAcceptorProxy<Functions::ParsedFunction<dim>>
      boundary_values_function;
  };


  // @sect4{The HeatEquation constructor}
  //
  // The constructor is responsible for initializing all member variables of
  // the class. This is relatively straightforward, and includes setting up
  // parameters to be set upon reading from an input file via the
  // ParameterAcceptor mechanism previously detailed in step-60 and step-70.
  template <int dim>
  HeatEquation<dim>::HeatEquation(const MPI_Comm mpi_communicator)
    : ParameterAcceptor("/Heat Equation/")
    , mpi_communicator(mpi_communicator)
    , pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
    , computing_timer(mpi_communicator,
                      pcout,
                      TimerOutput::summary,
                      TimerOutput::wall_times)
    , triangulation(mpi_communicator,
                    typename Triangulation<dim>::MeshSmoothing(
                      Triangulation<dim>::smoothing_on_refinement |
                      Triangulation<dim>::smoothing_on_coarsening))
    , fe(1)
    , dof_handler(triangulation)
    , time_stepper_data("",
                        "beuler",
                        /* start time */ 0.0,
                        /* end time */ 1.0,
                        /* initial time step */ 0.025)
    , initial_global_refinement(5)
    , max_delta_refinement_level(2)
    , mesh_adaptation_frequency(0)
    , right_hand_side_function("/Heat Equation/Right hand side")
    , initial_value_function("/Heat Equation/Initial value")
    , boundary_values_function("/Heat Equation/Boundary values")
  {
    enter_subsection("Time stepper");
    {
      enter_my_subsection(this->prm);
      {
        time_stepper_data.add_parameters(this->prm);
      }
      leave_my_subsection(this->prm);
    }
    leave_subsection();

    add_parameter("Initial global refinement",
                  initial_global_refinement,
                  "Number of times the mesh is refined globally before "
                  "starting the time stepping.");
    add_parameter("Maximum delta refinement level",
                  max_delta_refinement_level,
                  "Maximum number of local refinement levels.");
    add_parameter("Mesh adaptation frequency",
                  mesh_adaptation_frequency,
                  "When to adapt the mesh.");
  }



  // @sect4{The HeatEquation::setup_system() function}
  //
  // This function is not very different from what we do in many other programs
  // (including step-26). We enumerate degrees of freedom, output some
  // information about then, build constraint objects (recalling that we
  // put hanging node constraints into their separate object), and then
  // also build an AffineConstraint object that contains both the hanging
  // node constraints as well as constraints corresponding to zero Dirichlet
  // boundary conditions. This last object is needed since we impose the
  // constraints through algebraic equations. While technically it would be
  // possible to use the time derivative of the boundary function as a boundary
  // conditions for the time derivative of the solution, this is not done here.
  // Instead, we impose the boundary conditions through algebraic equations, and
  // therefore the time derivative of the boundary conditions is not part of
  // the algebraic system, and we need zero boundary conditions on the time
  // derivative of the solution when computing the residual. We use the
  // `homogeneous_constraints` object for this purpose.
  //
  // Note one detail here: The function
  // AffineConstraints::make_consistent_in_parallel() might expand the
  // underlying IndexSet for locally stored lines of the
  // `hanging_node_constraints` object depending on how constraints are set up
  // in parallel. Thus, at the point where we want to create a second affine
  // constraints object that gets the information from the hanging node
  // constraints, we need to make sure to use the local lines stored in the
  // hanging node constraints, not the `locally_relevant_dofs` that object was
  // originally initialized to. If this is not respected, errors might appear
  // during the course of the simulation.
  //
  // Finally, we create the actual non-homogeneous `current_constraints` by
  // calling `update_current_constraints). These are also used during the
  // assembly and during the residual evaluation.
  template <int dim>
  void HeatEquation<dim>::setup_system(const double time)
  {
    TimerOutput::Scope t(computing_timer, "setup system");

    dof_handler.distribute_dofs(fe);
    pcout << std::endl
          << "Number of active cells: " << triangulation.n_active_cells()
          << std::endl
          << "Number of degrees of freedom: " << dof_handler.n_dofs()
          << std::endl
          << std::endl;

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);


    hanging_node_constraints.clear();
    hanging_node_constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler,
                                            hanging_node_constraints);
    hanging_node_constraints.make_consistent_in_parallel(locally_owned_dofs,
                                                         locally_relevant_dofs,
                                                         mpi_communicator);
    hanging_node_constraints.close();


    homogeneous_constraints.clear();
    homogeneous_constraints.reinit(locally_owned_dofs,
                                   hanging_node_constraints.get_local_lines());
    homogeneous_constraints.merge(hanging_node_constraints);
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(),
                                             homogeneous_constraints);
    homogeneous_constraints.make_consistent_in_parallel(
      locally_owned_dofs,
      hanging_node_constraints.get_local_lines(),
      mpi_communicator);
    homogeneous_constraints.close();


    update_current_constraints(time);


    // The final block of code resets and initializes the matrix object with
    // the appropriate sparsity pattern. Recall that we do not store solution
    // vectors in this class (the time integrator object does that internally)
    // and so do not have to resize and initialize them either.
    DynamicSparsityPattern dsp(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    homogeneous_constraints,
                                    false);
    SparsityTools::distribute_sparsity_pattern(dsp,
                                               locally_owned_dofs,
                                               mpi_communicator,
                                               locally_relevant_dofs);

    jacobian_matrix.reinit(locally_owned_dofs,
                           locally_owned_dofs,
                           dsp,
                           mpi_communicator);
  }


  // @sect4{The HeatEquation::output_results() function}
  //
  // This function is called from "monitor" function that is called in turns
  // by the time stepper in each time step. We use it to write the solution to a
  // file, and provide graphical output through paraview or visit. We also write
  // a pvd file, which groups all metadata about the `.vtu` files into a single
  // file that can be used to load the full time dependent solution in paraview.
  template <int dim>
  void HeatEquation<dim>::output_results(const double       time,
                                         const unsigned int timestep_number,
                                         const PETScWrappers::MPI::Vector &y)
  {
    TimerOutput::Scope t(computing_timer, "output results");

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(y, "U");
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



  // @sect4{The HeatEquation::implicit_function() function}
  //
  // As discussed in the introduction, we describe the ODE system to the time
  // stepper via its residual,
  // @f[
  //   R(t,U,\dot U) = M \frac{\partial U(t)}{\partial t} + AU(t) - F(t).
  // @f]
  // The following function computes it, given vectors for $U,\dot U$ that
  // we will denote by `y` and `y_dot` because that's how they are called
  // in the documentation of the PETScWrappers::TimeStepper class.
  //
  // At the top of the function, we do the usual set up when computing
  // integrals. We face two minor difficulties here: the first is that
  // the `y` and `y_dot` vectors we get as input are read only, but we
  // need to make sure they satisfy the correct boundary conditions and so
  // have to set elements in these vectors. The second is that we need to
  // compute the residual, and therefore in general we need to evaluate solution
  // values and gradients inside locally owned cells, and for this need access
  // to degrees of freedom which may be owned by neighboring processors. To
  // address these issues, we create (non-ghosted) writable copies of the input
  // vectors, apply boundary conditions and hanging node current_constraints;
  // and then copy these vectors to ghosted vectors before we can do anything
  // sensible with them.
  template <int dim>
  void
  HeatEquation<dim>::implicit_function(const double                      time,
                                       const PETScWrappers::MPI::Vector &y,
                                       const PETScWrappers::MPI::Vector &y_dot,
                                       PETScWrappers::MPI::Vector &residual)
  {
    TimerOutput::Scope t(computing_timer, "implicit function");

    PETScWrappers::MPI::Vector tmp_solution(locally_owned_dofs,
                                            mpi_communicator);
    PETScWrappers::MPI::Vector tmp_solution_dot(locally_owned_dofs,
                                                mpi_communicator);
    tmp_solution     = y;
    tmp_solution_dot = y_dot;

    update_current_constraints(time);
    current_constraints.distribute(tmp_solution);
    homogeneous_constraints.distribute(tmp_solution_dot);

    PETScWrappers::MPI::Vector locally_relevant_solution(locally_owned_dofs,
                                                         locally_relevant_dofs,
                                                         mpi_communicator);
    PETScWrappers::MPI::Vector locally_relevant_solution_dot(
      locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    locally_relevant_solution     = tmp_solution;
    locally_relevant_solution_dot = tmp_solution_dot;


    const QGauss<dim> quadrature_formula(fe.degree + 1);
    FEValues<dim>     fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    std::vector<Tensor<1, dim>> solution_gradients(n_q_points);
    std::vector<double>         solution_dot_values(n_q_points);

    Vector<double> cell_residual(dofs_per_cell);

    right_hand_side_function.set_time(time);

    // Now for computing the actual residual. Recall that we wan to compute the
    // vector
    // @f[
    //   R(t,U,\dot U) = M \frac{\partial U(t)}{\partial t} + AU(t) - F(t).
    // @f]
    // We could do that by actually forming the matrices $M$ and $A$, but this
    // is not efficient. Instead, recall (by writing out how the elements of
    // $M$ and $A$ are defined, and exchanging integrals and sums) that the
    // $i$th element of the residual vector is given by
    // @f{align*}{
    //   R(t,U,\dot U)_i
    //     &= \sum_j \int_\Omega \varphi_i(\mathbf x, t) \varphi_j(\mathbf x, t)
    //     {\partial U_j(t)}{\partial t}
    //       + \sum_j \int_\Omega \nabla \varphi_i(\mathbf x, t) \cdot \nabla
    //       \varphi_j(\mathbf x, t) U_j(t)
    //       - \int_\Omega \varphi_i f(\mathbf x, t)
    //     \\ &=
    //     \int_\Omega \varphi_i(\mathbf x, t) u_h(\mathbf x, t)
    //       + \int_\Omega \nabla \varphi_i(\mathbf x, t) \cdot \nabla
    //         u_h(\mathbf x, t)
    //       - \int_\Omega \varphi_i f(\mathbf x, t).
    // @f}
    // We can compute these integrals efficiently by breaking them up into
    // a sum over all cells and then applying quadrature. For the integrand,
    // we need to evaluate the solution and its gradient at the quadrature
    // points within each locally owned cell, and for this, we need also
    // access to degrees of freedom that may be owned by neighboring
    // processors. We therefore use the locally_relevant_solution and and
    // locally_relevant_solution_dot vectors.
    residual = 0;
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);

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
                  (fe_values.shape_value(i, q) *       // [phi_i(x_q) *
                     solution_dot_values[q]            //  u(x_q)
                   +                                   //  +
                   fe_values.shape_grad(i, q) *        //  grad phi_i(x_q) *
                     solution_gradients[q]             //  grad u(x_q)
                   -                                   //  -
                   fe_values.shape_value(i, q) *       //  phi_i(x_q) *
                     right_hand_side_function.value(   //
                       fe_values.quadrature_point(q))) //  f(x_q)]
                  * fe_values.JxW(q);                  // * dx
              }
          current_constraints.distribute_local_to_global(cell_residual,
                                                         local_dof_indices,
                                                         residual);
        }
    residual.compress(VectorOperation::add);

    // The end result of the operations above is a vector that contains the
    // residual vector, having taken into account the constraints due to
    // hanging nodes and Dirichlet boundary conditions (by virtue of having
    // used `current_constraints.distribute_local_to_global()` to add the
    // local contributions to the global vector. At the end of the day, the
    // residual vector $r$ will be used in the solution of linear systems
    // of the form $J z = r$ with the "Jacobian" matrix that we define
    // below. We want to achieve that for algebraic components, the algebraic
    // components of $z$ have predictable values that achieve the purposes
    // discussed in the following. We do this by ensuring that the entries
    // corresponding to algebraic components in the residual $r$ have specific
    // values, and then we will do the same in the next function for the
    // matrix; for this, you will have to know that the rows and columns
    // of the matrix corresponding to constrained entries are zero with the
    // exception of the diagonal entries. We will manually set that diagonal
    // entry to one, and so $z_i=r_i$ for algebraic components.
    //
    // From the point of view of the residual vector, if the input `y`
    // vector does not contain the correct values on constrained degrees of
    // freedom (hanging nodes or boundary conditions), we need to communicate
    // this to the time stepper, and we do so by setting the residual to the
    // actual difference between the input `y` vector and the our local copy of
    // it, in which we have applied the constraints (see the top of the
    // function where we called `current_constraints.distribute(tmp_solution)`
    // and a similar operation on the time derivative). Since we have made a
    // copy of the input vector for this purpose, we use it to compute the
    // residual value. However, there is a difference between hanging nodes
    // constraints and boundary conditions: we do not want to make hanging node
    // constraints actually depend on their dependent degrees of freedom, since
    // this would imply that we are actually solving for the dependent degrees
    // of freedom. This is not what we are actually doing, however, since
    // hanging nodes are not actually solved for. They are eliminated from the
    // system by the call to AffineConstraints::distribute_local_to_global()
    // above. From the point of view of the Jacobian matrix, we are effectively
    // setting hanging nodes to an artificial value (usually zero), and we
    // simply want to make sure that we solve for those degrees of freedom a
    // posteriori, by calling the function AffineConstraints::distribute().
    //
    // Here we therefore check that the residual is equal to the input value on
    // the constrained dofs corresponding to hanging nodes (i.e., those for
    // which the lines of the `current_constraints` contain at least one other
    // entry), and to the difference between the input vector and the actual
    // solution on those constraints that correspond to boundary conditions.
    for (const auto &c : current_constraints.get_lines())
      if (locally_owned_dofs.is_element(c.index))
        {
          if (c.entries.empty()) /* no dependencies -> a Dirichlet node */
            residual[c.index] = y[c.index] - tmp_solution[c.index];
          else /* has dependencies -> a hanging node */
            residual[c.index] = y[c.index];
        }
    residual.compress(VectorOperation::insert);
  }


  // @sect4{The HeatEquation::assemble_implicit_jacobian() function}
  //
  // The next operation is to compute the "Jacobian", which PETSc TS defines
  // as the matrix
  // @f[
  //   J_\alpha = \dfrac{\partial R}{\partial y} + \alpha \dfrac{\partial
  //   R}{\partial \dot y}
  // @f]
  // which, for the current linear problem, is simply
  // @f[
  //   J_\alpha = A + \alpha M
  // @f]
  // and which is in particular independent of time and the current solution
  // vectors $y$ and $\dot y$.
  //
  // Having seen the assembly of matrices before, there is little that should
  // surprise you in the actual assembly here:
  template <int dim>
  void HeatEquation<dim>::assemble_implicit_jacobian(
    const double /* time */,
    const PETScWrappers::MPI::Vector & /* y */,
    const PETScWrappers::MPI::Vector & /* y_dot */,
    const double alpha)
  {
    TimerOutput::Scope t(computing_timer, "assemble implicit Jacobian");

    const QGauss<dim> quadrature_formula(fe.degree + 1);
    FEValues<dim>     fe_values(fe,
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
                    (fe_values.shape_grad(i, q) *      // grad phi_i(x_q) *
                       fe_values.shape_grad(j, q)      // grad phi_j(x_q)
                     + alpha *                         //
                         fe_values.shape_value(i, q) * // phi_i(x_q) *
                         fe_values.shape_value(j, q)   // phi_j(x_q)
                     ) *
                    fe_values.JxW(q); // * dx
                }
          current_constraints.distribute_local_to_global(cell_matrix,
                                                         local_dof_indices,
                                                         jacobian_matrix);
        }
    jacobian_matrix.compress(VectorOperation::add);

    // The only interesting part is the following. Recall that we modified the
    // residual vector's entries corresponding to the algebraic components of
    // the solution in the previous function. The outcome of calling
    // `current_constraints.distribute_local_to_global()` a few lines
    // above is that the global matrix has zero rows and columns for the
    // algebraic (constrained) components of the solution; the function
    // puts a value on the diagonal that is nonzero and has about the
    // same size as the remaining diagonal entries of the matrix. What
    // this diagonal value is is unknown to us -- in other cases where
    // we call `current_constraints.distribute_local_to_global()` on
    // both the left side matrix and the right side vector, as in most
    // other tutorial programs, the matrix diagonal entries and right
    // hand side values are chosen in such a way that the result of
    // solving a linear system is what we want it to be, but the scaling
    // is done automatically.
    //
    // This is not good enough for us here, because we are building the
    // right hand side independently from the matrix in different functions.
    // Thus, for any constrained degree of freedom, we set the diagonal of the
    // Jacobian to be one. This leaves the Jacobian matrix invertible,
    // consistent with what the time stepper expects, and it also makes sure
    // that if we did not make a mistake in the residual and/or in the Jacbian
    // matrix, then asking the time stepper to check the Jacobian with a finite
    // difference method will produce the correct result. This can be activated
    // at run time via passing the `-snes_test_jacobian` option on the command
    // line.
    for (const auto &c : current_constraints.get_lines())
      jacobian_matrix.set(c.index, c.index, 1.0);
    jacobian_matrix.compress(VectorOperation::insert);
  }


  // @sect4{The HeatEquation::solve_with_jacobian() function}
  //
  // This is the function that actually solves the linear system with the
  // Jacobian matrix we have previously built (in a call to the previous
  // function during the current time step or another earlier one -- time
  // steppers are quite sophisticated in determining internally whether it is
  // necessary to update the Jacobian matrix, or whether one can reuse it for
  // another time step without rebuilding it; this is similar to how one can
  // re-use the Newton matrix for several Newton steps, see for example the
  // discussion in step-77). We could in principle not provide this function to
  // the time stepper, and instead select a specific solver on the command line
  // by using the `-ksp_*` options of PETSc. However, by providing this
  // function, we can use a specific solver and preconditioner for the linear
  // system, and still have the possibility to change them on the command line.
  //
  // Providing a specific solver is more in line with the way we usually do
  // things in other deal.II examples, while letting PETSc choose a generic
  // solver, and changing it on the command line via `-ksp_type` is more in line
  // with the way PETSc is usually used, and it can be a convenient approach
  // when we are experimenting to find an optimal solver for our problem. Both
  // options are available here, since we can still change both the solver and
  // the preconditioner on the command line via `-user_ksp_type` and
  // `-user_pc_type` options.
  //
  // In any case, recall that the Jacobian we built in the previous function is
  // always of the form
  // @f[
  //   J_\alpha = \alpha M + A
  // @f]
  // where $M$ is a mass matrix and $A$ a Laplace matrix. $M$ is symmetric and
  // positive definite; $A$ is symmetric and at least positive semidefinite;
  // $\alpha> 0$. As a consequence, the Jacobian matrix is a symmetric and
  // positive definite matrix, which we can efficiently solve with the Conjugate
  // Gradient method, along with either SSOR or (if available) the algebraic
  // multigrid implementation provided by PETSc (via the Hypre package) as
  // preconditioner. In practice, if you wanted to solve "real" problems, one
  // would spend some time finding which preconditioner is optimal, perhaps
  // using PETSc's ability to read solver and preconditioner choices from the
  // command line. But this is not the focus of this tutorial program, and so
  // we just go with the following:
  template <int dim>
  void
  HeatEquation<dim>::solve_with_jacobian(const PETScWrappers::MPI::Vector &src,
                                         PETScWrappers::MPI::Vector       &dst)
  {
    TimerOutput::Scope t(computing_timer, "solve with Jacobian");

#if defined(PETSC_HAVE_HYPRE)
    PETScWrappers::PreconditionBoomerAMG preconditioner;
    preconditioner.initialize(jacobian_matrix);
#else
    PETScWrappers::PreconditionSSOR preconditioner;
    preconditioner.initialize(
      jacobian_matrix, PETScWrappers::PreconditionSSOR::AdditionalData(1.0));
#endif

    SolverControl           solver_control(1000, 1e-8 * src.l2_norm());
    PETScWrappers::SolverCG cg(solver_control);
    cg.set_prefix("user_");

    cg.solve(jacobian_matrix, dst, src, preconditioner);

    pcout << "     " << solver_control.last_step() << " linear iterations."
          << std::endl;
  }


  // @sect4{The HeatEquation::prepare_for_coarsening_and_refinement() function}
  //
  // The next block of functions deals with mesh refinement. We split this
  // process up into a "decide whether and what you want to refine" and a
  // "please transfer these vectors from old to new mesh" phase, where the first
  // also deals with marking cells for refinement. (The decision whether or
  // not to refine is done in the lambda function that calls the current
  // function.)
  //
  // Breaking things into a "mark cells" function and into a "execute mesh
  // adaptation and transfer solution vectors" function is awkward, though
  // conceptually not difficult to understand. These two pieces of code
  // should really be part of the same function, as they are in step-26.
  // The issue is with what PETScWrappers::TimeStepper provides us with in these
  // callbacks. Specifically, the "decide whether and what you want to refine"
  // callback has access to the current solution, and so can evaluate
  // (spatial) error estimators to decide which cells to refine. The
  // second callback that transfers vectors from old to new mesh gets a bunch
  // of vectors, but without the semantic information on which of these is
  // the current solution vector. As a consequence, it cannot do the marking
  // of cells for refinement or coarsening, and we have to do that from the
  // first callback.
  //
  // In practice, however, the problem is minor. The first of these
  // two functions is essentially identical to the
  // first half of the corresponding function in step-26, with the only
  // difference that it uses the parallel::distributed::GridRefinement namespace
  // function instead of the serial one. Once again, we make sure that we never
  // fall below the minimum refinement level, and above the maximum one, that we
  // can select from the parameter file.
  template <int dim>
  void HeatEquation<dim>::prepare_for_coarsening_and_refinement(
    const PETScWrappers::MPI::Vector &y)
  {
    PETScWrappers::MPI::Vector locally_relevant_solution(locally_owned_dofs,
                                                         locally_relevant_dofs,
                                                         mpi_communicator);
    locally_relevant_solution = y;

    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
    KellyErrorEstimator<dim>::estimate(dof_handler,
                                       QGauss<dim - 1>(fe.degree + 1),
                                       {},
                                       locally_relevant_solution,
                                       estimated_error_per_cell);

    parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction(
      triangulation, estimated_error_per_cell, 0.6, 0.4);

    const unsigned int max_grid_level =
      initial_global_refinement + max_delta_refinement_level;
    const unsigned int min_grid_level = initial_global_refinement;

    if (triangulation.n_levels() > max_grid_level)
      for (const auto &cell :
           triangulation.active_cell_iterators_on_level(max_grid_level))
        cell->clear_refine_flag();
    for (const auto &cell :
         triangulation.active_cell_iterators_on_level(min_grid_level))
      cell->clear_coarsen_flag();
  }


  // @sect4{The HeatEquation::transfer_solution_vectors_to_new_mesh() function}
  //
  // The following function then is the second half of the correspond function
  // in step-26. It is called by the time stepper whenever it requires to
  // transfer the solution and any intermediate stage vectors to a new mesh. We
  // must make sure that all input vectors are transformed into ghosted vectors
  // before the actual transfer is executed, and that we distribute the hanging
  // node constraints on the output vectors as soon as we have interpolated the
  // vectors to the new mesh -- i.e., that all constraints are satisfied on
  // the vectors we transfer.
  //
  // We have no way to enforce boundary conditions at this stage, however. This
  // is because the various vectors may correspond to solutions at previous time
  // steps if the method used here is a multistep time integrator, and so may
  // correspond to different time points that we are not privy to.
  //
  // While this could be a problem if we used the values of the solution and of
  // the intermediate stages on the constrained degrees of freedom to compute
  // the errors, we do not do so. Instead, we compute the errors on the
  // differential equation, and ignore algebraic constraints; therefore we do no
  // need to guarantee that the boundary conditions are satisfied also in the
  // intermediate stages.
  //
  // We have at our disposal the hanging node current_constraints alone, though,
  // and therefore we enforce them on the output vectors, even if this is not
  // really needed.
  template <int dim>
  void HeatEquation<dim>::transfer_solution_vectors_to_new_mesh(
    const double                                   time,
    const std::vector<PETScWrappers::MPI::Vector> &all_in,
    std::vector<PETScWrappers::MPI::Vector>       &all_out)
  {
    SolutionTransfer<dim, PETScWrappers::MPI::Vector> solution_trans(
      dof_handler);

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

    setup_system(time);

    all_out.resize(all_in.size());
    for (unsigned int i = 0; i < all_in.size(); ++i)
      {
        all_out[i].reinit(locally_owned_dofs, mpi_communicator);
        all_out_ptr[i] = &all_out[i];
      }
    solution_trans.interpolate(all_out_ptr);

    for (PETScWrappers::MPI::Vector &v : all_out)
      hanging_node_constraints.distribute(v);
  }


  // @sect4{The HeatEquation::update_current_constraints() function}
  //
  // Since regenerating the constraints at each time step may be expensive, we
  // make sure that we only do so when the time changes. We track time change by
  // checking if the time of the boundary_values_function has changed, with
  // respect to the time of the last call to this function. This will work most
  // of the times, but not the very first time we call this function, since the
  // time then may be zero and the time of the `boundary_values_function` is
  // zero at construction time. We therefore also check if the number
  // constraints in `current_constraints`, and if these are empty, we regenerate
  // the constraints regardless of the time variable.
  template <int dim>
  void HeatEquation<dim>::update_current_constraints(const double time)
  {
    if (Utilities::MPI::sum(current_constraints.n_constraints(),
                            mpi_communicator) == 0 ||
        time != boundary_values_function.get_time())
      {
        TimerOutput::Scope t(computing_timer, "update current constraints");

        boundary_values_function.set_time(time);
        current_constraints.clear();
        current_constraints.reinit(locally_owned_dofs,
                                   hanging_node_constraints.get_local_lines());
        current_constraints.merge(hanging_node_constraints);
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 0,
                                                 boundary_values_function,
                                                 current_constraints);
        current_constraints.make_consistent_in_parallel(locally_owned_dofs,
                                                        locally_relevant_dofs,
                                                        mpi_communicator);
        current_constraints.close();
      }
  }


  // @sect4{The HeatEquation::run() function}
  //
  // We have finally arrived at the main function of the class. At the top, it
  // creates the mesh and sets up the variables that make up the linear system:
  template <int dim>
  void HeatEquation<dim>::run()
  {
    GridGenerator::hyper_L(triangulation);
    triangulation.refine_global(initial_global_refinement);

    setup_system(/* time */ 0);

    // We then set up the time stepping object and associate the matrix we will
    // build whenever requested for both the Jacobian matrix (see the definition
    // above of what the "Jacobian" actually refers to) and for the matrix
    // that will be used as a preconditioner for the Jacobian.
    PETScWrappers::TimeStepper<PETScWrappers::MPI::Vector,
                               PETScWrappers::MPI::SparseMatrix>
      petsc_ts(time_stepper_data);

    petsc_ts.set_matrices(jacobian_matrix, jacobian_matrix);


    // The real work setting up the time stepping object starts here. As
    // discussed in the introduction, the way the PETScWrappers::TimeStepper
    // class is used is by inverting control: At the end of this function, we
    // will call PETScWrappers::TimeStepper::solve() which internally will
    // run the loop over time steps, and at the appropriate places *call back*
    // into user code for whatever functionality is required. What we need to
    // do is to hook up these callback functions by assigning appropriate
    // lambda functions to member variables of the `petsc_ts` object.
    //
    // We start by creating lambda functions that provide information about
    // the "implicit function" (i.e., that part of the right hand side of the
    // ODE system that we want to treat implicitly -- which in our case is
    // the entire right hand side), a function that assembles the Jacobian
    // matrix, and a function that solves a linear system with the Jacobian.
    petsc_ts.implicit_function = [&](const double                      time,
                                     const PETScWrappers::MPI::Vector &y,
                                     const PETScWrappers::MPI::Vector &y_dot,
                                     PETScWrappers::MPI::Vector       &res) {
      this->implicit_function(time, y, y_dot, res);
    };

    petsc_ts.setup_jacobian = [&](const double                      time,
                                  const PETScWrappers::MPI::Vector &y,
                                  const PETScWrappers::MPI::Vector &y_dot,
                                  const double                      alpha) {
      this->assemble_implicit_jacobian(time, y, y_dot, alpha);
    };

    petsc_ts.solve_with_jacobian = [&](const PETScWrappers::MPI::Vector &src,
                                       PETScWrappers::MPI::Vector       &dst) {
      this->solve_with_jacobian(src, dst);
    };

    // The next two callbacks deal with identifying and setting variables
    // that are considered "algebraic" (rather than "differential"), i.e., for
    // which we know what values they are supposed to have rather than letting
    // their values be determined by the differential equation. We need to
    // instruct the time stepper to ignore these components when computing the
    // error in the residuals, and we do so by first providing a function that
    // returns an IndexSet with the indices of these algebraic components of the
    // solution vector (or rather, that subset of the locally-owned part of the
    // vector that is algebraic, in case we are running in parallel). This first
    // of the following two functions does that. Specifically, both nodes at
    // which Dirichlet boundary conditions are applied, and hanging nodes are
    // algebraically constrained. This function then returns a set of
    // indices that is initially empty (but knows about the size of the index
    // space) and which we then construct as the union of boundary and hanging
    // node indices.
    //
    // Following this, we then also need a function that, given a solution
    // vector `y` and the current time, sets the algebraic components of that
    // vector to their correct value. This comes down to ensuring that we have
    // up to date constraints in the `constraints` variable, and then applying
    // these constraints to the solution vector via
    // AffineConstraints::distribute(). (It is perhaps worth noting that we
    // *could* have achieved the same in `solve_with_jacobian()`. Whenever the
    // time stepper solves a linear system, it follows up the call to the solver
    // by calling the callback to set algebraic components correct. We could
    // also have put the calls to `update_current_constraints()` and
    // `distribute()` into the `solve_with_jacobian` function, but by not doing
    // so, we can also replace the `solve_with_jacobian` function with a call to
    // a PETSc solver, and we would still have the current_constraints correctly
    // applied to the solution vector.)
    petsc_ts.algebraic_components = [&]() {
      IndexSet algebraic_set(dof_handler.n_dofs());
      algebraic_set.add_indices(DoFTools::extract_boundary_dofs(dof_handler));
      algebraic_set.add_indices(
        DoFTools::extract_hanging_node_dofs(dof_handler));
      return algebraic_set;
    };

    petsc_ts.update_constrained_components =
      [&](const double time, PETScWrappers::MPI::Vector &y) {
        TimerOutput::Scope t(computing_timer, "set algebraic components");
        update_current_constraints(time);
        current_constraints.distribute(y);
      };


    // The next two callbacks relate to mesh refinement. As discussed in the
    // introduction, PETScWrappers::TimeStepper knows how to deal with the
    // situation where we want to change the mesh. All we have to provide
    // is a callback that returns `true` if we are at a point where we want
    // to refine the mesh (and `false` otherwise) and that if we want to
    // do mesh refinement does some prep work for that in the form of
    // calling the `prepare_for_coarsening_and_refinement` function.
    //
    // If the first callback below returns `true`, then PETSc TS will
    // do some clean-up operations, and call the second of the
    // callback functions
    // (`petsc_ts.transfer_solution_vectors_to_new_mesh`) with a
    // collection of vectors that need to be interpolated from the old
    // to the new mesh. This may include the current solution, perhaps
    // the current time derivative of the solution, and in the case of
    // [multistep time
    // integrators](https://en.wikipedia.org/wiki/Linear_multistep_method) also
    // the solutions of some previous time steps. We hand all of these over to
    // the `interpolate()` member function of this class.
    petsc_ts.decide_and_prepare_for_remeshing =
      [&](const double /* time */,
          const unsigned int                step_number,
          const PETScWrappers::MPI::Vector &y) -> bool {
      if (step_number > 0 && this->mesh_adaptation_frequency > 0 &&
          step_number % this->mesh_adaptation_frequency == 0)
        {
          pcout << std::endl << "Adapting the mesh..." << std::endl;
          this->prepare_for_coarsening_and_refinement(y);
          return true;
        }
      else
        return false;
    };

    petsc_ts.transfer_solution_vectors_to_new_mesh =
      [&](const double                                   time,
          const std::vector<PETScWrappers::MPI::Vector> &all_in,
          std::vector<PETScWrappers::MPI::Vector>       &all_out) {
        this->transfer_solution_vectors_to_new_mesh(time, all_in, all_out);
      };

    // The final callback is a "monitor" that is called in each
    // time step. Here we use it to create a graphical output. Perhaps a better
    // scheme would output the solution at fixed time intervals, rather
    // than in every time step, but this is not the main point of
    // this program and so we go with the easy approach:
    petsc_ts.monitor = [&](const double                      time,
                           const PETScWrappers::MPI::Vector &y,
                           const unsigned int                step_number) {
      pcout << "Time step " << step_number << " at t=" << time << std::endl;
      this->output_results(time, step_number, y);
    };


    // With all of this out of the way, the rest of the function is
    // anticlimactic: We just have to initialize the solution vector with
    // the initial conditions and call the function that does the time
    // stepping, and everything else will happen automatically:
    PETScWrappers::MPI::Vector solution(locally_owned_dofs, mpi_communicator);
    VectorTools::interpolate(dof_handler, initial_value_function, solution);

    petsc_ts.solve(solution);
  }
} // namespace Step86


// @sect3{The main() function}
//
// The rest of the program is as it always looks. We read run-time parameters
// from an input file via the ParameterAcceptor class in the same way as
// we showed in step-60 and step-70.
int main(int argc, char **argv)
{
  try
    {
      using namespace Step86;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      HeatEquation<2>                  heat_equation_solver(MPI_COMM_WORLD);

      const std::string input_filename =
        (argc > 1 ? argv[1] : "heat_equation.prm");
      ParameterAcceptor::initialize(input_filename, "heat_equation_used.prm");
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
