/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2011 - 2013 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Authors: Joerg Frohne, Texas A&M University and
 *                        University of Siegen, 2011, 2012
 *          Wolfgang Bangerth, Texas A&M University, 2012
 */


// @sect3{Include files}

// As usual, at the beginning we include all the header files we need in
// here. With the exception of the various files that provide interfaces to
// the Trilinos library, there are no surprises:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_precondition.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>
#include <list>


namespace Step41
{
  using namespace dealii;

  // @sect3{The <code>ObstacleProblem</code> class template}

  // This class supplies all function and variables needed to describe the
  // obstacle problem. It is close to what we had to do in step-4, and so
  // relatively simple. The only real new components are the
  // update_solution_and_constraints function that computes the active set and
  // a number of variables that are necessary to describe the original
  // (unconstrained) form of the linear system
  // (<code>complete_system_matrix</code> and
  // <code>complete_system_rhs</code>) as well as the active set itself and
  // the diagonal of the mass matrix $B$ used in scaling Lagrange multipliers
  // in the active set formulation. The rest is as in step-4:
  template <int dim>
  class ObstacleProblem
  {
  public:
    ObstacleProblem ();
    void run ();

  private:
    void make_grid ();
    void setup_system();
    void assemble_system ();
    void assemble_mass_matrix_diagonal (TrilinosWrappers::SparseMatrix &mass_matrix);
    void update_solution_and_constraints ();
    void solve ();
    void output_results (const unsigned int iteration) const;

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;
    ConstraintMatrix     constraints;
    IndexSet             active_set;

    TrilinosWrappers::SparseMatrix system_matrix;
    TrilinosWrappers::SparseMatrix complete_system_matrix;

    TrilinosWrappers::Vector       solution;
    TrilinosWrappers::Vector       system_rhs;
    TrilinosWrappers::Vector       complete_system_rhs;
    TrilinosWrappers::Vector       diagonal_of_mass_matrix;
    TrilinosWrappers::Vector       contact_force;
  };


  // @sect3{Right hand side, boundary values, and the obstacle}

  // In the following, we define classes that describe the right hand side
  // function, the Dirichlet boundary values, and the height of the obstacle
  // as a function of $\mathbf x$. In all three cases, we derive these classes
  // from Function@<dim@>, although in the case of <code>RightHandSide</code>
  // and <code>Obstacle</code> this is more out of convention than necessity
  // since we never pass such objects to the library. In any case, the
  // definition of the right hand side and boundary values classes is obvious
  // given our choice of $f=-10$, $u|_{\partial\Omega}=0$:
  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    RightHandSide () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
  };

  template <int dim>
  double RightHandSide<dim>::value (const Point<dim> &p,
                                    const unsigned int component) const
  {
    Assert (component == 0, ExcNotImplemented());

    return -10;
  }



  template <int dim>
  class BoundaryValues : public Function<dim>
  {
  public:
    BoundaryValues () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
  };

  template <int dim>
  double BoundaryValues<dim>::value (const Point<dim> &p,
                                     const unsigned int component) const
  {
    Assert (component == 0, ExcNotImplemented());

    return 0;
  }



  // We describe the obstacle function by a cascaded barrier (think: stair
  // steps):
  template <int dim>
  class Obstacle : public Function<dim>
  {
  public:
    Obstacle () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
  };

  template <int dim>
  double Obstacle<dim>::value (const Point<dim> &p,
                               const unsigned int component) const
  {
    Assert (component == 0, ExcNotImplemented());

    if (p (0) < -0.5)
      return -0.2;
    else if (p (0) >= -0.5 && p (0) < 0.0)
      return -0.4;
    else if (p (0) >= 0.0 && p (0) < 0.5)
      return -0.6;
    else
      return -0.8;
  }



  // @sect3{Implementation of the <code>ObstacleProblem</code> class}


  // @sect4{ObstacleProblem::ObstacleProblem}

  // To everyone who has taken a look at the first few tutorial programs, the
  // constructor is completely obvious:
  template <int dim>
  ObstacleProblem<dim>::ObstacleProblem ()
    :
    fe (1),
    dof_handler (triangulation)
  {}


  // @sect4{ObstacleProblem::make_grid}

  // We solve our obstacle problem on the square $[-1,1]\times [-1,1]$ in
  // 2D. This function therefore just sets up one of the simplest possible
  // meshes.
  template <int dim>
  void ObstacleProblem<dim>::make_grid ()
  {
    GridGenerator::hyper_cube (triangulation, -1, 1);
    triangulation.refine_global (7);

    std::cout << "Number of active cells: "
              << triangulation.n_active_cells()
              << std::endl
              << "Total number of cells: "
              << triangulation.n_cells()
              << std::endl;
  }


  // @sect4{ObstacleProblem::setup_system}

  // In this first function of note, we set up the degrees of freedom handler,
  // resize vectors and matrices, and deal with the constraints. Initially,
  // the constraints are, of course, only given by boundary values, so we
  // interpolate them towards the top of the function.
  template <int dim>
  void ObstacleProblem<dim>::setup_system ()
  {
    dof_handler.distribute_dofs (fe);
    active_set.set_size (dof_handler.n_dofs());

    std::cout << "Number of degrees of freedom: "
              << dof_handler.n_dofs()
              << std::endl
              << std::endl;

    VectorTools::interpolate_boundary_values (dof_handler,
                                              0,
                                              BoundaryValues<dim>(),
                                              constraints);
    constraints.close ();

    CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler,
                                     c_sparsity,
                                     constraints,
                                     false);

    system_matrix.reinit (c_sparsity);
    complete_system_matrix.reinit (c_sparsity);

    solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());
    complete_system_rhs.reinit (dof_handler.n_dofs());
    contact_force.reinit (dof_handler.n_dofs());

    // The only other thing to do here is to compute the factors in the $B$
    // matrix which is used to scale the residual. As discussed in the
    // introduction, we'll use a little trick to make this mass matrix
    // diagonal, and in the following then first compute all of this as a
    // matrix and then extract the diagonal elements for later use:
    TrilinosWrappers::SparseMatrix mass_matrix;
    mass_matrix.reinit (c_sparsity);
    assemble_mass_matrix_diagonal (mass_matrix);
    diagonal_of_mass_matrix.reinit (dof_handler.n_dofs());
    for (unsigned int j=0; j<solution.size (); j++)
      diagonal_of_mass_matrix (j) = mass_matrix.diag_element (j);
  }


  // @sect4{ObstacleProblem::assemble_system}

  // This function at once assembles the system matrix and right-hand-side and
  // applied the constraints (both due to the active set as well as from
  // boundary values) to our system. Otherwise, it is functionally equivalent
  // to the corresponding function in, for example, step-4.
  template <int dim>
  void ObstacleProblem<dim>::assemble_system ()
  {
    std::cout << "   Assembling system..." << std::endl;

    system_matrix = 0;
    system_rhs    = 0;

    const QGauss<dim>         quadrature_formula(fe.degree+1);
    const RightHandSide<dim>  right_hand_side;

    FEValues<dim>             fe_values (fe, quadrature_formula,
                                         update_values   | update_gradients |
                                         update_quadrature_points |
                                         update_JxW_values);

    const unsigned int        dofs_per_cell = fe.dofs_per_cell;
    const unsigned int        n_q_points    = quadrature_formula.size();

    FullMatrix<double>        cell_matrix (dofs_per_cell, dofs_per_cell);
    TrilinosWrappers::Vector  cell_rhs (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    for (; cell!=endc; ++cell)
      {
        fe_values.reinit (cell);
        cell_matrix = 0;
        cell_rhs = 0;

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              for (unsigned int j=0; j<dofs_per_cell; ++j)
                cell_matrix(i,j) += (fe_values.shape_grad (i, q_point) *
                                     fe_values.shape_grad (j, q_point) *
                                     fe_values.JxW (q_point));

              cell_rhs(i) += (fe_values.shape_value (i, q_point) *
                              right_hand_side.value (fe_values.quadrature_point (q_point)) *
                              fe_values.JxW (q_point));
            }

        cell->get_dof_indices (local_dof_indices);

        constraints.distribute_local_to_global (cell_matrix,
                                                cell_rhs,
                                                local_dof_indices,
                                                system_matrix,
                                                system_rhs,
                                                true);
      }
  }



  // @sect4{ObstacleProblem::assemble_mass_matrix_diagonal}

  // The next function is used in the computation of the diagonal mass matrix
  // $B$ used to scale variables in the active set method. As discussed in the
  // introduction, we get the mass matrix to be diagonal by choosing the
  // trapezoidal rule for quadrature. Doing so we don't really need the triple
  // loop over quadrature points, indices $i$ and indices $j$ any more and
  // can, instead, just use a double loop. The rest of the function is obvious
  // given what we have discussed in many of the previous tutorial programs.
  //
  // Note that at the time this function is called, the constraints object
  // only contains boundary value constraints; we therefore do not have to pay
  // attention in the last copy-local-to-global step to preserve the values of
  // matrix entries that may later on be constrained by the active set.
  //
  // Note also that the trick with the trapezoidal rule only works if we have
  // in fact $Q_1$ elements. For higher order elements, one would need to use
  // a quadrature formula that has quadrature points at all the support points
  // of the finite element. Constructing such a quadrature formula isn't
  // really difficult, but not the point here, and so we simply assert at the
  // top of the function that our implicit assumption about the finite element
  // is in fact satisfied.
  template <int dim>
  void
  ObstacleProblem<dim>::
  assemble_mass_matrix_diagonal (TrilinosWrappers::SparseMatrix &mass_matrix)
  {
    Assert (fe.degree == 1, ExcNotImplemented());

    const QTrapez<dim>        quadrature_formula;
    FEValues<dim>             fe_values (fe,
                                         quadrature_formula,
                                         update_values   |
                                         update_JxW_values);

    const unsigned int        dofs_per_cell = fe.dofs_per_cell;
    const unsigned int        n_q_points    = quadrature_formula.size();

    FullMatrix<double>        cell_matrix (dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    for (; cell!=endc; ++cell)
      {
        fe_values.reinit (cell);
        cell_matrix = 0;

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            cell_matrix(i,i) += (fe_values.shape_value (i, q_point) *
                                 fe_values.shape_value (i, q_point) *
                                 fe_values.JxW (q_point));

        cell->get_dof_indices (local_dof_indices);

        constraints.distribute_local_to_global (cell_matrix,
                                                local_dof_indices,
                                                mass_matrix);
      }
  }


  // @sect4{ObstacleProblem::update_solution_and_constraints}

  // In a sense, this is the central function of this program.  It updates the
  // active set of constrained degrees of freedom as discussed in the
  // introduction and computes a ConstraintMatrix object from it that can then
  // be used to eliminate constrained degrees of freedom from the solution of
  // the next iteration. At the same time we set the constrained degrees of
  // freedom of the solution to the correct value, namely the height of the
  // obstacle.
  //
  // Fundamentally, the function is rather simple: We have to loop over all
  // degrees of freedom and check the sign of the function $\Lambda^k_i +
  // c([BU^k]_i - G_i) = \Lambda^k_i + cB_i(U^k_i - [g_h]_i)$ because in our
  // case $G_i = B_i[g_h]_i$. To this end, we use the formula given in the
  // introduction by which we can compute the Lagrange multiplier as the
  // residual of the original linear system (given via the variables
  // <code>complete_system_matrix</code> and <code>complete_system_rhs</code>.
  // At the top of this function, we compute this residual using a function
  // that is part of the matrix classes.
  template <int dim>
  void
  ObstacleProblem<dim>::update_solution_and_constraints ()
  {
    std::cout << "   Updating active set..." << std::endl;

    const double penalty_parameter = 100.0;

    TrilinosWrappers::Vector lambda (dof_handler.n_dofs());
    complete_system_matrix.residual (lambda,
                                     solution, complete_system_rhs);
    contact_force.ratio (lambda, diagonal_of_mass_matrix);
    contact_force *= -1;

    // The next step is to reset the active set and constraints objects and to
    // start the loop over all degrees of freedom. This is made slightly more
    // complicated by the fact that we can't just loop over all elements of
    // the solution vector since there is no way for us then to find out what
    // location a DoF is associated with; however, we need this location to
    // test whether the displacement of a DoF is larger or smaller than the
    // height of the obstacle at this location.
    //
    // We work around this by looping over all cells and DoFs defined on each
    // of these cells. We use here that the displacement is described using a
    // $Q_1$ function for which degrees of freedom are always located on the
    // vertices of the cell; thus, we can get the index of each degree of
    // freedom and its location by asking the vertex for this information. On
    // the other hand, this clearly wouldn't work for higher order elements,
    // and so we add an assertion that makes sure that we only deal with
    // elements for which all degrees of freedom are located in vertices to
    // avoid tripping ourselves with non-functional code in case someone wants
    // to play with increasing the polynomial degree of the solution.
    //
    // The price to pay for having to loop over cells rather than DoFs is that
    // we may encounter some degrees of freedom more than once, namely each
    // time we visit one of the cells adjacent to a given vertex. We will
    // therefore have to keep track which vertices we have already touched and
    // which we haven't so far. We do so by using an array of flags
    // <code>dof_touched</code>:
    constraints.clear();
    active_set.clear ();

    const Obstacle<dim> obstacle;
    std::vector<bool>   dof_touched (dof_handler.n_dofs(), false);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
        {
          Assert (dof_handler.get_fe().dofs_per_cell ==
                  GeometryInfo<dim>::vertices_per_cell,
                  ExcNotImplemented());

          const unsigned int dof_index = cell->vertex_dof_index (v,0);

          if (dof_touched[dof_index] == false)
            dof_touched[dof_index] = true;
          else
            continue;

          // Now that we know that we haven't touched this DoF yet, let's get
          // the value of the displacement function there as well as the value
          // of the obstacle function and use this to decide whether the
          // current DoF belongs to the active set. For that we use the
          // function given above and in the introduction.
          //
          // If we decide that the DoF should be part of the active set, we
          // add its index to the active set, introduce an inhomogeneous
          // equality constraint in the ConstraintMatrix object, and reset the
          // solution value to the height of the obstacle. Finally, the
          // residual of the non-contact part of the system serves as an
          // additional control (the residual equals the remaining,
          // unaccounted forces, and should be zero outside the contact zone),
          // so we zero out the components of the residual vector (i.e., the
          // Lagrange multiplier lambda) that correspond to the area where the
          // body is in contact; at the end of the loop over all cells, the
          // residual will therefore only consist of the residual in the
          // non-contact zone. We output the norm of this residual along with
          // the size of the active set after the loop.
          const double obstacle_value = obstacle.value (cell->vertex(v));
          const double solution_value = solution (dof_index);

          if (lambda (dof_index) +
              penalty_parameter *
              diagonal_of_mass_matrix(dof_index) *
              (solution_value - obstacle_value)
              <
              0)
            {
              active_set.add_index (dof_index);
              constraints.add_line (dof_index);
              constraints.set_inhomogeneity (dof_index, obstacle_value);

              solution (dof_index) = obstacle_value;

              lambda (dof_index) = 0;
            }
        }
    std::cout << "      Size of active set: " << active_set.n_elements()
              << std::endl;

    std::cout << "   Residual of the non-contact part of the system: "
              << lambda.l2_norm()
              << std::endl;

    // In a final step, we add to the set of constraints on DoFs we have so
    // far from the active set those that result from Dirichlet boundary
    // values, and close the constraints object:
    VectorTools::interpolate_boundary_values (dof_handler,
                                              0,
                                              BoundaryValues<dim>(),
                                              constraints);
    constraints.close ();
  }

  // @sect4{ObstacleProblem::solve}

  // There is nothing to say really about the solve function. In the context
  // of a Newton method, we are not typically interested in very high accuracy
  // (why ask for a highly accurate solution of a linear problem that we know
  // only gives us an approximation of the solution of the nonlinear problem),
  // and so we use the ReductionControl class that stops iterations when
  // either an absolute tolerance is reached (for which we choose $10^{-12}$)
  // or when the residual is reduced by a certain factor (here, $10^{-3}$).
  template <int dim>
  void ObstacleProblem<dim>::solve ()
  {
    std::cout << "   Solving system..." << std::endl;

    ReductionControl                    reduction_control (100, 1e-12, 1e-3);
    SolverCG<TrilinosWrappers::Vector>  solver (reduction_control);
    TrilinosWrappers::PreconditionAMG   precondition;
    precondition.initialize (system_matrix);

    solver.solve (system_matrix, solution, system_rhs, precondition);
    constraints.distribute (solution);

    std::cout << "      Error: " << reduction_control.initial_value()
              << " -> " << reduction_control.last_value()
              << " in "
              <<  reduction_control.last_step()
              << " CG iterations."
              << std::endl;
  }


  // @sect4{ObstacleProblem::output_results}

  // We use the vtk-format for the output.  The file contains the displacement
  // and a numerical representation of the active set. The function looks
  // standard but note that we can add an IndexSet object to the DataOut
  // object in exactly the same way as a regular solution vector: it is simply
  // interpreted as a function that is either zero (when a degree of freedom
  // is not part of the IndexSet) or one (if it is).
  template <int dim>
  void ObstacleProblem<dim>::output_results (const unsigned int iteration) const
  {
    std::cout << "   Writing graphical output..." << std::endl;

    DataOut<dim> data_out;

    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, "displacement");
    data_out.add_data_vector (active_set, "active_set");
    data_out.add_data_vector (contact_force, "lambda");

    data_out.build_patches ();

    std::ofstream output_vtk ((std::string("output_") +
                               Utilities::int_to_string (iteration, 3) +
                               ".vtk").c_str ());
    data_out.write_vtk (output_vtk);
  }



  // @sect4{ObstacleProblem::run}

  // This is the function which has the top-level control over everything.  It
  // is not very long, and in fact rather straightforward: in every iteration
  // of the active set method, we assemble the linear system, solve it, update
  // the active set and project the solution back to the feasible set, and
  // then output the results. The iteration is terminated whenever the active
  // set has not changed in the previous iteration.
  //
  // The only trickier part is that we have to save the linear system (i.e.,
  // the matrix and right hand side) after assembling it in the first
  // iteration. The reason is that this is the only step where we can access
  // the linear system as built without any of the contact constraints
  // active. We need this to compute the residual of the solution at other
  // iterations, but in other iterations that linear system we form has the
  // rows and columns that correspond to constrained degrees of freedom
  // eliminated, and so we can no longer access the full residual of the
  // original equation.
  template <int dim>
  void ObstacleProblem<dim>::run ()
  {
    make_grid();
    setup_system ();

    IndexSet active_set_old (active_set);
    for (unsigned int iteration=0; iteration<=solution.size (); ++iteration)
      {
        std::cout << "Newton iteration " << iteration << std::endl;

        assemble_system ();

        if (iteration == 0)
          {
            complete_system_matrix.copy_from (system_matrix);
            complete_system_rhs = system_rhs;
          }

        solve ();
        update_solution_and_constraints ();
        output_results (iteration);

        if (active_set == active_set_old)
          break;

        active_set_old = active_set;

        std::cout << std::endl;
      }
  }
}


// @sect3{The <code>main</code> function}

// And this is the main function. It follows the pattern of all other main
// functions. The call to initialize MPI exists because the Trilinos library
// upon which we build our linear solvers in this program requires it.
int main (int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      using namespace Step41;

      deallog.depth_console (0);

      Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);

      ObstacleProblem<2> obstacle_problem;
      obstacle_problem.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
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
      std::cerr << std::endl << std::endl
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
