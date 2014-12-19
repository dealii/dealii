/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2012 - 2014 by the deal.II authors
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
 * Author: Sven Wetterauer, University of Heidelberg, 2012
 */


// @sect3{Include files}

// The first few files have already been covered in previous examples and will
// thus not be further commented on.
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
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
  // - The <code>setup_system</code> function takes an argument that denotes whether
  //   this is the first time it is called or not. The difference is that the
  //   first time around we need to distribute the degrees of freedom and set the
  //   solution vector for $u^n$ to the correct size. The following times, the
  //   function is called after we have already done these steps as part of
  //   refining the mesh in <code>refine_mesh</code>.
  // - We then also need new functions: <code>set_boundary_values()</code>
  //   takes care of setting the boundary values on the solution vector
  //   correctly, as discussed at the end of the
  //   introduction. <code>compute_residual()</code> is a function that computes
  //   the norm of the nonlinear (discrete) residual. We use this function to
  //   monitor convergence of the Newton iteration. The function takes a step
  //   length $\alpha^n$ as argument to compute the residual of $u^n + \alpha^n
  //   \; \delta u^n$. This is something one typically needs for step length
  //   control, although we will not use this feature here. Finally,
  //   <code>determine_step_length()</code> computes the step length $\alpha^n$
  //   in each Newton iteration. As discussed in the introduction, we here use a
  //   fixed step length and leave implementing a better strategy as an
  //   exercise.

  template <int dim>
  class MinimalSurfaceProblem
  {
  public:
    MinimalSurfaceProblem ();
    ~MinimalSurfaceProblem ();

    void run ();

  private:
    void setup_system (const bool initial_step);
    void assemble_system ();
    void solve ();
    void refine_mesh ();
    void set_boundary_values ();
    double compute_residual (const double alpha) const;
    double determine_step_length () const;

    Triangulation<dim>   triangulation;

    DoFHandler<dim>      dof_handler;
    FE_Q<dim>            fe;

    ConstraintMatrix     hanging_node_constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       present_solution;
    Vector<double>       newton_update;
    Vector<double>       system_rhs;
  };

  // @sect3{Boundary condition}

  // The boundary condition is implemented just like in step-4.  It is chosen
  // as $g(x,y)=\sin(2 \pi (x+y))$:

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
                                     const unsigned int /*component*/) const
  {
    return std::sin(2 * numbers::PI * (p[0]+p[1]));
  }

  // @sect3{The <code>MinimalSurfaceProblem</code> class implementation}

  // @sect4{MinimalSurfaceProblem::MinimalSurfaceProblem}

  // The constructor and destructor of the class are the same as in the first
  // few tutorials.

  template <int dim>
  MinimalSurfaceProblem<dim>::MinimalSurfaceProblem ()
    :
    dof_handler (triangulation),
    fe (2)
  {}



  template <int dim>
  MinimalSurfaceProblem<dim>::~MinimalSurfaceProblem ()
  {
    dof_handler.clear ();
  }

  // @sect4{MinimalSurfaceProblem::setup_system}

  // As always in the setup-system function, we setup the variables of the
  // finite element method. There are same differences to step-6, because
  // there we start solving the PDE from scratch in every refinement cycle
  // whereas here we need to take the solution from the previous mesh onto the
  // current mesh. Consequently, we can't just reset solution vectors. The
  // argument passed to this function thus indicates whether we can
  // distributed degrees of freedom (plus compute constraints) and set the
  // solution vector to zero or whether this has happened elsewhere already
  // (specifically, in <code>refine_mesh()</code>).

  template <int dim>
  void MinimalSurfaceProblem<dim>::setup_system (const bool initial_step)
  {
    if (initial_step)
      {
        dof_handler.distribute_dofs (fe);
        present_solution.reinit (dof_handler.n_dofs());

        hanging_node_constraints.clear ();
        DoFTools::make_hanging_node_constraints (dof_handler,
                                                 hanging_node_constraints);
        hanging_node_constraints.close ();
      }


    // The remaining parts of the function are the same as in step-6.

    newton_update.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());

    CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);

    hanging_node_constraints.condense (c_sparsity);

    sparsity_pattern.copy_from(c_sparsity);
    system_matrix.reinit (sparsity_pattern);
  }

  // @sect4{MinimalSurfaceProblem::assemble_system}

  // This function does the same as in the previous tutorials except that now,
  // of course, the matrix and right hand side functions depend on the
  // previous iteration's solution. As discussed in the introduction, we need
  // to use zero boundary values for the Newton updates; we compute them at
  // the end of this function.
  //
  // The top of the function contains the usual boilerplate code, setting up
  // the objects that allow us to evaluate shape functions at quadrature
  // points and temporary storage locations for the local matrices and
  // vectors, as well as for the gradients of the previous solution at the
  // quadrature points. We then start the loop over all cells:
  template <int dim>
  void MinimalSurfaceProblem<dim>::assemble_system ()
  {
    const QGauss<dim>  quadrature_formula(3);

    system_matrix = 0;
    system_rhs = 0;

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_gradients         |
                             update_quadrature_points |
                             update_JxW_values);

    const unsigned int           dofs_per_cell = fe.dofs_per_cell;
    const unsigned int           n_q_points    = quadrature_formula.size();

    FullMatrix<double>           cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>               cell_rhs (dofs_per_cell);

    std::vector<Tensor<1, dim> > old_solution_gradients(n_q_points);

    std::vector<types::global_dof_index>    local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        cell_matrix = 0;
        cell_rhs = 0;

        fe_values.reinit (cell);

        // For the assembly of the linear system, we have to obtain the values
        // of the previous solution's gradients at the quadrature
        // points. There is a standard way of doing this: the
        // FEValues::get_function_gradients function takes a vector that
        // represents a finite element field defined on a DoFHandler, and
        // evaluates the gradients of this field at the quadrature points of the
        // cell with which the FEValues object has last been reinitialized.
        // The values of the gradients at all quadrature points are then written
        // into the second argument:
        fe_values.get_function_gradients(present_solution,
                                         old_solution_gradients);

        // With this, we can then do the integration loop over all quadrature
        // points and shape functions.  Having just computed the gradients of
        // the old solution in the quadrature points, we are able to compute
        // the coefficients $a_{n}$ in these points.  The assembly of the
        // system itself then looks similar to what we always do with the
        // exception of the nonlinear terms, as does copying the results from
        // the local objects into the global ones:
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          {
            const double coeff
              = 1.0 / std::sqrt(1 +
                                old_solution_gradients[q_point] *
                                old_solution_gradients[q_point]);

            for (unsigned int i=0; i<dofs_per_cell; ++i)
              {
                for (unsigned int j=0; j<dofs_per_cell; ++j)
                  {
                    cell_matrix(i, j) += (fe_values.shape_grad(i, q_point)
                                          * coeff
                                          * (fe_values.shape_grad(j, q_point)
                                             -
                                             coeff * coeff
                                             * (fe_values.shape_grad(j, q_point)
                                                *
                                                old_solution_gradients[q_point])
                                             * old_solution_gradients[q_point]
                                            )
                                          * fe_values.JxW(q_point));
                  }

                cell_rhs(i) -= (fe_values.shape_grad(i, q_point)
                                * coeff
                                * old_solution_gradients[q_point]
                                * fe_values.JxW(q_point));
              }
          }

        cell->get_dof_indices (local_dof_indices);
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              system_matrix.add (local_dof_indices[i],
                                 local_dof_indices[j],
                                 cell_matrix(i,j));

            system_rhs(local_dof_indices[i]) += cell_rhs(i);
          }
      }

    // Finally, we remove hanging nodes from the system and apply zero
    // boundary values to the linear system that defines the Newton updates
    // $\delta u^n$:
    hanging_node_constraints.condense (system_matrix);
    hanging_node_constraints.condense (system_rhs);

    std::map<types::global_dof_index,double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler,
                                              0,
                                              ZeroFunction<dim>(),
                                              boundary_values);
    MatrixTools::apply_boundary_values (boundary_values,
                                        system_matrix,
                                        newton_update,
                                        system_rhs);
  }



  // @sect4{MinimalSurfaceProblem::solve}

  // The solve function is the same as always. At the end of the solution
  // process we update the current solution by setting
  // $u^{n+1}=u^n+\alpha^n\;\delta u^n$.
  template <int dim>
  void MinimalSurfaceProblem<dim>::solve ()
  {
    SolverControl solver_control (system_rhs.size(),
                                  system_rhs.l2_norm()*1e-6);
    SolverCG<>    solver (solver_control);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    solver.solve (system_matrix, newton_update, system_rhs,
                  preconditioner);

    hanging_node_constraints.distribute (newton_update);

    const double alpha = determine_step_length();
    present_solution.add (alpha, newton_update);
  }


  // @sect4{MinimalSurfaceProblem::refine_mesh}

  // The first part of this function is the same as in step-6... However,
  // after refining the mesh we have to transfer the old solution to the new
  // one which we do with the help of the SolutionTransfer class. The process
  // is slightly convoluted, so let us describe it in detail:
  template <int dim>
  void MinimalSurfaceProblem<dim>::refine_mesh ()
  {
    Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

    KellyErrorEstimator<dim>::estimate (dof_handler,
                                        QGauss<dim-1>(3),
                                        typename FunctionMap<dim>::type(),
                                        present_solution,
                                        estimated_error_per_cell);

    GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                     estimated_error_per_cell,
                                                     0.3, 0.03);

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
    triangulation.prepare_coarsening_and_refinement ();

    // With this out of the way, we initialize a SolutionTransfer object with
    // the present DoFHandler and attach the solution vector to it, followed
    // by doing the actual refinement and distribution of degrees of freedom
    // on the new mesh
    SolutionTransfer<dim> solution_transfer(dof_handler);
    solution_transfer.prepare_for_coarsening_and_refinement(present_solution);

    triangulation.execute_coarsening_and_refinement();

    dof_handler.distribute_dofs(fe);

    // Finally, we retrieve the old solution interpolated to the new
    // mesh. Since the SolutionTransfer function does not actually store the
    // values of the old solution, but rather indices, we need to preserve the
    // old solution vector until we have gotten the new interpolated
    // values. Thus, we have the new values written into a temporary vector,
    // and only afterwards write them into the solution vector object. Once we
    // have this solution we have to make sure that the $u^n$ we now have
    // actually has the correct boundary values. As explained at the end of
    // the introduction, this is not automatically the case even if the
    // solution before refinement had the correct boundary values, and so we
    // have to explicitly make sure that it now has:
    Vector<double> tmp(dof_handler.n_dofs());
    solution_transfer.interpolate(present_solution, tmp);
    present_solution = tmp;

    set_boundary_values ();

    // On the new mesh, there are different hanging nodes, which we have to
    // compute again. To ensure there are no hanging nodes of the old mesh in
    // the object, it's first cleared.  To be on the safe side, we then also
    // make sure that the current solution's vector entries satisfy the
    // hanging node constraints (see the discussion in the documentation of
    // the SolutionTransfer class for why this is necessary):
    hanging_node_constraints.clear();

    DoFTools::make_hanging_node_constraints(dof_handler,
                                            hanging_node_constraints);
    hanging_node_constraints.close();

    hanging_node_constraints.distribute (present_solution);

    // We end the function by updating all the remaining data structures,
    // indicating to <code>setup_dofs()</code> that this is not the first
    // go-around and that it needs to preserve the content of the solution
    // vector:
    setup_system (false);
  }



  // @sect4{MinimalSurfaceProblem::set_boundary_values}

  // The next function ensures that the solution vector's entries respect the
  // boundary values for our problem.  Having refined the mesh (or just
  // started computations), there might be new nodal points on the
  // boundary. These have values that are simply interpolated from the
  // previous mesh (or are just zero), instead of the correct boundary
  // values. This is fixed up by setting all boundary nodes explicit to the
  // right value:
  template <int dim>
  void MinimalSurfaceProblem<dim>::set_boundary_values ()
  {
    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler,
                                              0,
                                              BoundaryValues<dim>(),
                                              boundary_values);
    for (std::map<types::global_dof_index, double>::const_iterator
         p = boundary_values.begin();
         p != boundary_values.end(); ++p)
      present_solution(p->first) = p->second;
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
  double MinimalSurfaceProblem<dim>::compute_residual (const double alpha) const
  {
    Vector<double> residual (dof_handler.n_dofs());

    Vector<double> evaluation_point (dof_handler.n_dofs());
    evaluation_point = present_solution;
    evaluation_point.add (alpha, newton_update);

    const QGauss<dim>  quadrature_formula(3);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_gradients         |
                             update_quadrature_points |
                             update_JxW_values);

    const unsigned int           dofs_per_cell = fe.dofs_per_cell;
    const unsigned int           n_q_points    = quadrature_formula.size();

    Vector<double>               cell_residual (dofs_per_cell);
    std::vector<Tensor<1, dim> > gradients(n_q_points);

    std::vector<types::global_dof_index>    local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        cell_residual = 0;
        fe_values.reinit (cell);

        // The actual computation is much as in
        // <code>assemble_system()</code>. We first evaluate the gradients of
        // $u^n+\alpha^n\,\delta u^n$ at the quadrature points, then compute
        // the coefficient $a_n$, and then plug it all into the formula for
        // the residual:
        fe_values.get_function_gradients (evaluation_point,
                                          gradients);


        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          {
            const double coeff = 1/std::sqrt(1 +
                                             gradients[q_point] *
                                             gradients[q_point]);

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              cell_residual(i) -= (fe_values.shape_grad(i, q_point)
                                   * coeff
                                   * gradients[q_point]
                                   * fe_values.JxW(q_point));
          }

        cell->get_dof_indices (local_dof_indices);
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          residual(local_dof_indices[i]) += cell_residual(i);
      }

    // At the end of this function we also have to deal with the hanging node
    // constraints and with the issue of boundary values. With regard to the
    // latter, we have to set to zero the elements of the residual vector for
    // all entries that correspond to degrees of freedom that sit at the
    // boundary. The reason is that because the value of the solution there is
    // fixed, they are of course no "real" degrees of freedom and so, strictly
    // speaking, we shouldn't have assembled entries in the residual vector
    // for them. However, as we always do, we want to do exactly the same
    // thing on every cell and so we didn't not want to deal with the question
    // of whether a particular degree of freedom sits at the boundary in the
    // integration above. Rather, we will simply set to zero these entries
    // after the fact. To this end, we first need to determine which degrees
    // of freedom do in fact belong to the boundary and then loop over all of
    // those and set the residual entry to zero. This happens in the following
    // lines which we have already seen used in step-11:
    hanging_node_constraints.condense (residual);

    std::vector<bool> boundary_dofs (dof_handler.n_dofs());
    DoFTools::extract_boundary_dofs (dof_handler,
                                     ComponentMask(),
                                     boundary_dofs);
    for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
      if (boundary_dofs[i] == true)
        residual(i) = 0;

    // At the end of the function, we return the norm of the residual:
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
  // in the results section.
  template <int dim>
  double MinimalSurfaceProblem<dim>::determine_step_length() const
  {
    return 0.1;
  }



  // @sect4{MinimalSurfaceProblem::run}

  // In the run function, we build the first grid and then have the top-level
  // logic for the Newton iteration. The function has two variables, one that
  // indicates whether this is the first time we solve for a Newton update and
  // one that indicates the refinement level of the mesh:
  template <int dim>
  void MinimalSurfaceProblem<dim>::run ()
  {
    unsigned int refinement = 0;
    bool         first_step = true;

    // As described in the introduction, the domain is the unit disk around
    // the origin, created in the same way as shown in step-6. The mesh is
    // globally refined twice followed later on by several adaptive cycles:
    GridGenerator::hyper_ball (triangulation);
    static const HyperBallBoundary<dim> boundary;
    triangulation.set_boundary (0, boundary);
    triangulation.refine_global(2);

    // The Newton iteration starts next. During the first step we do not have
    // information about the residual prior to this step and so we continue
    // the Newton iteration until we have reached at least one iteration and
    // until residual is less than $10^{-3}$.
    //
    // At the beginning of the loop, we do a bit of setup work. In the first
    // go around, we compute the solution on the twice globally refined mesh
    // after setting up the basic data structures. In all following mesh
    // refinement loops, the mesh will be refined adaptively.
    double previous_res = 0;
    while (first_step || (previous_res>1e-3))
      {
        if (first_step == true)
          {
            std::cout << "******** Initial mesh "
                      << " ********"
                      << std::endl;

            setup_system (true);
            set_boundary_values ();
          }
        else
          {
            ++refinement;
            std::cout << "******** Refined mesh " << refinement
                      << " ********"
                      << std::endl;

            refine_mesh();
          }

        // On every mesh we do exactly five Newton steps. We print the initial
        // residual here and then start the iterations on this mesh.
        //
        // In every Newton step the system matrix and the right hand side have
        // to be computed first, after which we store the norm of the right
        // hand side as the residual to check against when deciding whether to
        // stop the iterations. We then solve the linear system (the function
        // also updates $u^{n+1}=u^n+\alpha^n\;\delta u^n$) and output the
        // residual at the end of this Newton step:
        std::cout << "  Initial residual: "
                  << compute_residual(0)
                  << std::endl;

        for (unsigned int inner_iteration=0; inner_iteration<5; ++inner_iteration)
          {
            assemble_system ();
            previous_res = system_rhs.l2_norm();

            solve ();

            first_step = false;
            std::cout << "  Residual: "
                      << compute_residual(0)
                      << std::endl;
          }

        // Every fifth iteration, i.e., just before we refine the mesh again,
        // we output the solution as well as the Newton update. This happens
        // as in all programs before:
        DataOut<dim> data_out;

        data_out.attach_dof_handler (dof_handler);
        data_out.add_data_vector (present_solution, "solution");
        data_out.add_data_vector (newton_update, "update");
        data_out.build_patches ();
        const std::string filename = "solution-" +
                                     Utilities::int_to_string (refinement, 2) +
                                     ".vtk";
        std::ofstream output (filename.c_str());
        data_out.write_vtk (output);

      }
  }
}

// @sect4{The main function}

// Finally the main function. This follows the scheme of all other main
// functions:
int main ()
{
  try
    {
      using namespace dealii;
      using namespace Step15;

      deallog.depth_console (0);

      MinimalSurfaceProblem<2> laplace_problem_2d;
      laplace_problem_2d.run ();
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
