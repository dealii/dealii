/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2012 - 2024 by the deal.II authors
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
 * Author: Sven Wetterauer, University of Heidelberg, 2012
 */


// @sect3{Include files}

// The first few files have already been covered in previous examples and will
// thus not be further commented on.
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>


#include <fstream>
#include <iostream>

// We will use adaptive mesh refinement between Newton iterations. To do so,
// we need to be able to work with a solution on the new mesh, although it was
// computed on the old one. The SolutionTransfer class transfers the solution
// from the old to the new mesh:

#include <deal.II/numerics/solution_transfer.h>

// We then open a namespace for this program and import everything from the
// dealii namespace into it, as in previous programs:
namespace Step15
{
  using namespace dealii;


  // @sect3{The <code>MinimalSurfaceProblem</code> class template}

  // The class template is basically the same as in step-6.  Three additions
  // are made:
  // - There are two solution vectors, one for the Newton update
  //   $\delta u^n$, and one for the current iterate $u^n$.
  // - The single AffineConstraints<> object in step-6 that is used to store
  //   boundary conditions and hanging node constraints, is replaced by two
  //   different objects of the same type: `zero_constraints` and
  //   `nonzero_constraints`. The former contains homogeneous boundary
  //   conditions to be used for the residual and solution updates, while the
  //   latter contains the correct boundary conditions for the solution. Both
  //   objects also contain the hanging nodes constraints.
  // - The <code>setup_system</code> function takes an argument that denotes
  //   whether this is the first time it is called or not. The difference is
  //   that the first time around we need to distribute the degrees of freedom
  //   and set the solution vector for $u^n$ to the correct size. The following
  //   times, the function is called after we have already done these steps as
  //   part of refining the mesh in <code>refine_mesh</code>.
  // - We then also need a few new functions:
  //   <code>compute_residual()</code> is a function that computes
  //   the norm of the nonlinear (discrete) residual. We use this function to
  //   monitor convergence of the Newton iteration. The function takes a step
  //   length $\alpha^n$ as argument to compute the residual of $u^n + \alpha^n
  //   \; \delta u^n$. This is something one typically needs for step length
  //   control, although we will not use this feature here. Finally,
  //   <code>determine_step_length()</code> computes the step length $\alpha^n$
  //   in each Newton iteration. As discussed in the introduction, we here use a
  //   fixed step length and leave implementing a better strategy as an
  //   exercise. (step-77 does this differently: It simply uses an
  //   external package for the whole solution process, and a good
  //   line search strategy is part of what that package provides.)

  template <int dim>
  class MinimalSurfaceProblem
  {
  public:
    MinimalSurfaceProblem();
    void run();

  private:
    void   setup_system();
    void   assemble_system();
    void   solve();
    void   refine_mesh();
    double compute_residual(const double alpha) const;
    double determine_step_length() const;
    void   output_results(const unsigned int refinement_cycle) const;

    Triangulation<dim> triangulation;

    DoFHandler<dim> dof_handler;
    const FE_Q<dim> fe;

    AffineConstraints<double> zero_constraints;
    AffineConstraints<double> nonzero_constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> current_solution;
    Vector<double> newton_update;
    Vector<double> system_rhs;
  };

  // @sect3{Boundary condition}

  // The boundary condition is implemented just like in step-4.  It is chosen
  // as $g(x,y)=\sin(2 \pi (x+y))$:

  template <int dim>
  class BoundaryValues : public Function<dim>
  {
  public:
    virtual double value(const Point<dim>  &p,
                         const unsigned int component = 0) const override;
  };


  template <int dim>
  double BoundaryValues<dim>::value(const Point<dim> &p,
                                    const unsigned int /*component*/) const
  {
    return std::sin(2 * numbers::PI * (p[0] + p[1]));
  }

  // @sect3{The <code>MinimalSurfaceProblem</code> class implementation}

  // @sect4{MinimalSurfaceProblem::MinimalSurfaceProblem}

  // The constructor and destructor of the class are the same as in the first
  // few tutorials.

  template <int dim>
  MinimalSurfaceProblem<dim>::MinimalSurfaceProblem()
    : dof_handler(triangulation)
    , fe(2)
  {}


  // @sect4{MinimalSurfaceProblem::setup_system}

  // As always in the setup-system function, we set up the variables of the
  // finite element method. There are some differences to step-6, because
  // we need to construct two AffineConstraint<> objects.
  template <int dim>
  void MinimalSurfaceProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);
    current_solution.reinit(dof_handler.n_dofs());

    zero_constraints.clear();
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(),
                                             zero_constraints);
    DoFTools::make_hanging_node_constraints(dof_handler, zero_constraints);
    zero_constraints.close();

    nonzero_constraints.clear();
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             BoundaryValues<dim>(),
                                             nonzero_constraints);

    DoFTools::make_hanging_node_constraints(dof_handler, nonzero_constraints);
    nonzero_constraints.close();

    newton_update.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp, zero_constraints);

    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);
  }

  // @sect4{MinimalSurfaceProblem::assemble_system}

  // This function does the same as in the previous tutorials except that now,
  // of course, the matrix and right hand side functions depend on the
  // previous iteration's solution. As discussed in the introduction, we need
  // to use zero boundary values for the Newton updates; this is done by using
  // the `zero_constraints` object when assembling into the global matrix and
  // vector.
  //
  // The top of the function contains the usual boilerplate code, setting up
  // the objects that allow us to evaluate shape functions at quadrature
  // points and temporary storage locations for the local matrices and
  // vectors, as well as for the gradients of the previous solution at the
  // quadrature points. We then start the loop over all cells:
  template <int dim>
  void MinimalSurfaceProblem<dim>::assemble_system()
  {
    const QGauss<dim> quadrature_formula(fe.degree + 1);

    system_matrix = 0;
    system_rhs    = 0;

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_gradients | update_quadrature_points |
                              update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    std::vector<Tensor<1, dim>> old_solution_gradients(n_q_points);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_matrix = 0;
        cell_rhs    = 0;

        fe_values.reinit(cell);

        // For the assembly of the linear system, we have to obtain the values
        // of the previous solution's gradients at the quadrature
        // points. There is a standard way of doing this: the
        // FEValues::get_function_gradients function takes a vector that
        // represents a finite element field defined on a DoFHandler, and
        // evaluates the gradients of this field at the quadrature points of the
        // cell with which the FEValues object has last been reinitialized.
        // The values of the gradients at all quadrature points are then written
        // into the second argument:
        fe_values.get_function_gradients(current_solution,
                                         old_solution_gradients);

        // With this, we can then do the integration loop over all quadrature
        // points and shape functions.  Having just computed the gradients of
        // the old solution in the quadrature points, we are able to compute
        // the coefficients $a_{n}$ in these points.  The assembly of the
        // system itself then looks similar to what we always do with the
        // exception of the nonlinear terms, as does copying the results from
        // the local objects into the global ones:
        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            const double coeff =
              1.0 / std::sqrt(1 + old_solution_gradients[q] *
                                    old_solution_gradients[q]);

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  cell_matrix(i, j) +=
                    (((fe_values.shape_grad(i, q)      // ((\nabla \phi_i
                       * coeff                         //   * a_n
                       * fe_values.shape_grad(j, q))   //   * \nabla \phi_j)
                      -                                //  -
                      (fe_values.shape_grad(i, q)      //  (\nabla \phi_i
                       * coeff * coeff * coeff         //   * a_n^3
                       * (fe_values.shape_grad(j, q)   //   * (\nabla \phi_j
                          * old_solution_gradients[q]) //      * \nabla u_n)
                       * old_solution_gradients[q]))   //   * \nabla u_n)))
                     * fe_values.JxW(q));              // * dx

                cell_rhs(i) -= (fe_values.shape_grad(i, q)  // \nabla \phi_i
                                * coeff                     // * a_n
                                * old_solution_gradients[q] // * \nabla u_n
                                * fe_values.JxW(q));        // * dx
              }
          }

        cell->get_dof_indices(local_dof_indices);
        zero_constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }
  }



  // @sect4{MinimalSurfaceProblem::solve}

  // The solve function is the same as always. At the end of the solution
  // process we update the current solution by setting
  // $u^{n+1}=u^n+\alpha^n\;\delta u^n$.
  template <int dim>
  void MinimalSurfaceProblem<dim>::solve()
  {
    SolverControl            solver_control(system_rhs.size(),
                                 system_rhs.l2_norm() * 1e-6);
    SolverCG<Vector<double>> solver(solver_control);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    solver.solve(system_matrix, newton_update, system_rhs, preconditioner);

    zero_constraints.distribute(newton_update);

    const double alpha = determine_step_length();
    current_solution.add(alpha, newton_update);
  }


  // @sect4{MinimalSurfaceProblem::refine_mesh}

  // The first part of this function is the same as in step-6... However,
  // after refining the mesh we have to transfer the old solution to the new
  // one which we do with the help of the SolutionTransfer class. The process
  // is slightly convoluted, so let us describe it in detail:
  template <int dim>
  void MinimalSurfaceProblem<dim>::refine_mesh()
  {
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

    KellyErrorEstimator<dim>::estimate(
      dof_handler,
      QGauss<dim - 1>(fe.degree + 1),
      std::map<types::boundary_id, const Function<dim> *>(),
      current_solution,
      estimated_error_per_cell);

    GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                    estimated_error_per_cell,
                                                    0.3,
                                                    0.03);

    // Then we need an additional step: if, for example, you flag a cell that
    // is once more refined than its neighbor, and that neighbor is not
    // flagged for refinement, we would end up with a jump of two refinement
    // levels across a cell interface.  To avoid these situations, the library
    // will silently also have to refine the neighbor cell once. It does so by
    // calling the Triangulation::prepare_coarsening_and_refinement function
    // before actually doing the refinement and coarsening.  This function
    // flags a set of additional cells for refinement or coarsening, to
    // enforce rules like the one-hanging-node rule.  The cells that are
    // flagged for refinement and coarsening after calling this function are
    // exactly the ones that will actually be refined or coarsened. Usually,
    // you don't have to do this by hand
    // (Triangulation::execute_coarsening_and_refinement does this for
    // you). However, we need to initialize the SolutionTransfer class and it
    // needs to know the final set of cells that will be coarsened or refined
    // in order to store the data from the old mesh and transfer to the new
    // one. Thus, we call the function by hand:
    triangulation.prepare_coarsening_and_refinement();

    // With this out of the way, we initialize a SolutionTransfer object with
    // the present DoFHandler. We make a copy of the solution vector and attach
    // it to the SolutionTransfer. Now we can actually execute the refinement
    // and create the new matrices and vectors including the vector
    // `current_solution`, that will hold the current solution on the new mesh
    // after calling SolutionTransfer::interpolate():
    SolutionTransfer<dim> solution_transfer(dof_handler);
    const Vector<double>  coarse_solution = current_solution;
    solution_transfer.prepare_for_coarsening_and_refinement(coarse_solution);

    triangulation.execute_coarsening_and_refinement();

    setup_system();

    solution_transfer.interpolate(current_solution);

    // On the new mesh, there are different hanging nodes, computed in
    // `setup_system()` above. To be on the safe side, we should  make sure that
    // the current solution's vector entries satisfy the hanging node
    // constraints (see the discussion in the documentation of the
    // SolutionTransfer class for why this is necessary) and boundary values. As
    // explained at the end of the introduction, the interpolated solution does
    // not automatically satisfy the boundary values even if the solution before
    // refinement had the correct boundary values.
    nonzero_constraints.distribute(current_solution);
  }



  // @sect4{MinimalSurfaceProblem::compute_residual}

  // In order to monitor convergence, we need a way to compute the norm of the
  // (discrete) residual, i.e., the norm of the vector
  // $\left<F(u^n),\varphi_i\right>$ with $F(u)=-\nabla \cdot \left(
  // \frac{1}{\sqrt{1+|\nabla u|^{2}}}\nabla u \right)$ as discussed in the
  // introduction. It turns out that (although we don't use this feature in
  // the current version of the program) one needs to compute the residual
  // $\left<F(u^n+\alpha^n\;\delta u^n),\varphi_i\right>$ when determining
  // optimal step lengths, and so this is what we implement here: the function
  // takes the step length $\alpha^n$ as an argument. The original
  // functionality is of course obtained by passing a zero as argument.
  //
  // In the function below, we first set up a vector for the residual, and
  // then a vector for the evaluation point $u^n+\alpha^n\;\delta u^n$. This
  // is followed by the same boilerplate code we use for all integration
  // operations:
  template <int dim>
  double MinimalSurfaceProblem<dim>::compute_residual(const double alpha) const
  {
    Vector<double> residual(dof_handler.n_dofs());

    Vector<double> evaluation_point(dof_handler.n_dofs());
    evaluation_point = current_solution;
    evaluation_point.add(alpha, newton_update);

    const QGauss<dim> quadrature_formula(fe.degree + 1);
    FEValues<dim>     fe_values(fe,
                            quadrature_formula,
                            update_gradients | update_quadrature_points |
                              update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    Vector<double>              cell_residual(dofs_per_cell);
    std::vector<Tensor<1, dim>> gradients(n_q_points);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_residual = 0;
        fe_values.reinit(cell);

        // The actual computation is much as in
        // <code>assemble_system()</code>. We first evaluate the gradients of
        // $u^n+\alpha^n\,\delta u^n$ at the quadrature points, then compute
        // the coefficient $a_n$, and then plug it all into the formula for
        // the residual:
        fe_values.get_function_gradients(evaluation_point, gradients);


        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            const double coeff =
              1. / std::sqrt(1 + gradients[q] * gradients[q]);

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              cell_residual(i) -= (fe_values.shape_grad(i, q) // \nabla \phi_i
                                   * coeff                    // * a_n
                                   * gradients[q]             // * \nabla u_n
                                   * fe_values.JxW(q));       // * dx
          }

        cell->get_dof_indices(local_dof_indices);
        zero_constraints.distribute_local_to_global(cell_residual,
                                                    local_dof_indices,
                                                    residual);
      }

    return residual.l2_norm();
  }



  // @sect4{MinimalSurfaceProblem::determine_step_length}

  // As discussed in the introduction, Newton's method frequently does not
  // converge if we always take full steps, i.e., compute $u^{n+1}=u^n+\delta
  // u^n$. Rather, one needs a damping parameter (step length) $\alpha^n$ and
  // set $u^{n+1}=u^n+\alpha^n\delta u^n$. This function is the one called
  // to compute $\alpha^n$.
  //
  // Here, we simply always return 0.1. This is of course a sub-optimal
  // choice: ideally, what one wants is that the step size goes to one as we
  // get closer to the solution, so that we get to enjoy the rapid quadratic
  // convergence of Newton's method. We will discuss better strategies below
  // in the results section, and step-77 also covers this aspect.
  template <int dim>
  double MinimalSurfaceProblem<dim>::determine_step_length() const
  {
    return 0.1;
  }



  // @sect4{MinimalSurfaceProblem::output_results}

  // This last function to be called from `run()` outputs the current solution
  // (and the Newton update) in graphical form as a VTU file. It is entirely the
  // same as what has been used in previous tutorials.
  template <int dim>
  void MinimalSurfaceProblem<dim>::output_results(
    const unsigned int refinement_cycle) const
  {
    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(current_solution, "solution");
    data_out.add_data_vector(newton_update, "update");
    data_out.build_patches();

    const std::string filename =
      "solution-" + Utilities::int_to_string(refinement_cycle, 2) + ".vtu";
    std::ofstream output(filename);
    data_out.write_vtu(output);
  }


  // @sect4{MinimalSurfaceProblem::run}

  // In the run function, we build the first grid and then have the top-level
  // logic for the Newton iteration.
  //
  // As described in the introduction, the domain is the unit disk around
  // the origin, created in the same way as shown in step-6. The mesh is
  // globally refined twice followed later on by several adaptive cycles.
  //
  // Before starting the Newton loop, we also need to do
  // ensure that the first Newton iterate already has the correct
  // boundary values, as discussed in the introduction.
  template <int dim>
  void MinimalSurfaceProblem<dim>::run()
  {
    GridGenerator::hyper_ball(triangulation);
    triangulation.refine_global(2);

    setup_system();
    nonzero_constraints.distribute(current_solution);

    // The Newton iteration starts next. We iterate until the (norm of the)
    // residual computed at the end of the previous iteration is less than
    // $10^{-3}$, as checked at the end of the `do { ... } while` loop that
    // starts here. Because we don't have a reasonable value to initialize
    // the variable, we just use the largest value that can be represented
    // as a `double`.
    double       last_residual_norm = std::numeric_limits<double>::max();
    unsigned int refinement_cycle   = 0;
    do
      {
        std::cout << "Mesh refinement step " << refinement_cycle << std::endl;

        if (refinement_cycle != 0)
          refine_mesh();

        // On every mesh we do exactly five Newton steps. We print the initial
        // residual here and then start the iterations on this mesh.
        //
        // In every Newton step the system matrix and the right hand side have
        // to be computed first, after which we store the norm of the right
        // hand side as the residual to check against when deciding whether to
        // stop the iterations. We then solve the linear system (the function
        // also updates $u^{n+1}=u^n+\alpha^n\;\delta u^n$) and output the
        // norm of the residual at the end of this Newton step.
        //
        // After the end of this loop, we then also output the solution on the
        // current mesh in graphical form and increment the counter for the
        // mesh refinement cycle.
        std::cout << "  Initial residual: " << compute_residual(0) << std::endl;

        for (unsigned int inner_iteration = 0; inner_iteration < 5;
             ++inner_iteration)
          {
            assemble_system();
            last_residual_norm = system_rhs.l2_norm();

            solve();

            std::cout << "  Residual: " << compute_residual(0) << std::endl;
          }

        output_results(refinement_cycle);

        ++refinement_cycle;
        std::cout << std::endl;
      }
    while (last_residual_norm > 1e-2);
  }
} // namespace Step15

// @sect4{The main function}

// Finally the main function. This follows the scheme of all other main
// functions:
int main()
{
  try
    {
      using namespace Step15;

      MinimalSurfaceProblem<2> problem;
      problem.run();
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
