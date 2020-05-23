/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2000 - 2020 by the deal.II authors
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
 * Author: Wolfgang Bangerth, University of Heidelberg, 2000
 */


// @sect3{Include files}

// As usual, the first few include files are already known, so we will not
// comment on them further.
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

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
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

// In this example, we need vector-valued finite elements. The support for
// these can be found in the following include file:
#include <deal.II/fe/fe_system.h>
// We will compose the vector-valued finite elements from regular Q1 elements
// which can be found here, as usual:
#include <deal.II/fe/fe_q.h>

// This again is C++:
#include <fstream>
#include <iostream>

// The last step is as in previous programs. In particular, just like in
// step-7, we pack everything that's specific to this program into a namespace
// of its own.
namespace Step8
{
  using namespace dealii;

  // @sect3{The <code>ElasticProblem</code> class template}

  // The main class is, except for its name, almost unchanged with respect to
  // the step-6 example.
  //
  // The only change is the use of a different class for the <code>fe</code>
  // variable: Instead of a concrete finite element class such as FE_Q, we now
  // use a more generic one, FESystem. In fact, FESystem is not really a
  // finite element itself in that it does not implement shape functions of
  // its own. Rather, it is a class that can be used to stack several other
  // elements together to form one vector-valued finite element. In our case,
  // we will compose the vector-valued element of <code>FE_Q(1)</code>
  // objects, as shown below in the constructor of this class.
  template <int dim>
  class ElasticProblem
  {
  public:
    ElasticProblem();
    void run();

  private:
    void setup_system();
    void assemble_system();
    void solve();
    void refine_grid();
    void output_results(const unsigned int cycle) const;

    Triangulation<dim> triangulation;
    DoFHandler<dim>    dof_handler;

    FESystem<dim> fe;

    AffineConstraints<double> constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> system_rhs;
  };


  // @sect3{Right hand side values}

  // Before going over to the implementation of the main class, we declare and
  // define the function which describes the right hand side. This time, the
  // right hand side is vector-valued, as is the solution, so we will describe
  // the changes required for this in some more detail.
  //
  // To prevent cases where the return vector has not previously been set to
  // the right size we test for this case and otherwise throw an exception at
  // the beginning of the function. Note that enforcing that output arguments
  // already have the correct size is a convention in deal.II, and enforced
  // almost everywhere. The reason is that we would otherwise have to check at
  // the beginning of the function and possibly change the size of the output
  // vector. This is expensive, and would almost always be unnecessary (the
  // first call to the function would set the vector to the right size, and
  // subsequent calls would only have to do redundant checks). In addition,
  // checking and possibly resizing the vector is an operation that can not be
  // removed if we can't rely on the assumption that the vector already has
  // the correct size; this is in contract to the Assert call that is
  // completely removed if the program is compiled in optimized mode.
  //
  // Likewise, if by some accident someone tried to compile and run the
  // program in only one space dimension (in which the elastic equations do
  // not make much sense since they reduce to the ordinary Laplace equation),
  // we terminate the program in the second assertion. The program will work
  // just fine in 3d, however.
  template <int dim>
  void right_hand_side(const std::vector<Point<dim>> &points,
                       std::vector<Tensor<1, dim>> &  values)
  {
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));
    Assert(dim >= 2, ExcNotImplemented());

    // The rest of the function implements computing force values. We will use
    // a constant (unit) force in x-direction located in two little circles
    // (or spheres, in 3d) around points (0.5,0) and (-0.5,0), and y-force in
    // an area around the origin; in 3d, the z-component of these centers is
    // zero as well.
    //
    // For this, let us first define two objects that denote the centers of
    // these areas. Note that upon construction of the Point objects, all
    // components are set to zero.
    Point<dim> point_1, point_2;
    point_1(0) = 0.5;
    point_2(0) = -0.5;

    for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
      {
        // If <code>points[point_n]</code> is in a circle (sphere) of radius
        // 0.2 around one of these points, then set the force in x-direction
        // to one, otherwise to zero:
        if (((points[point_n] - point_1).norm_square() < 0.2 * 0.2) ||
            ((points[point_n] - point_2).norm_square() < 0.2 * 0.2))
          values[point_n][0] = 1.0;
        else
          values[point_n][0] = 0.0;

        // Likewise, if <code>points[point_n]</code> is in the vicinity of the
        // origin, then set the y-force to one, otherwise to zero:
        if (points[point_n].norm_square() < 0.2 * 0.2)
          values[point_n][1] = 1.0;
        else
          values[point_n][1] = 0.0;
      }
  }



  // @sect3{The <code>ElasticProblem</code> class implementation}

  // @sect4{ElasticProblem::ElasticProblem constructor}

  // Following is the constructor of the main class. As said before, we would
  // like to construct a vector-valued finite element that is composed of
  // several scalar finite elements (i.e., we want to build the vector-valued
  // element so that each of its vector components consists of the shape
  // functions of a scalar element). Of course, the number of scalar finite
  // elements we would like to stack together equals the number of components
  // the solution function has, which is <code>dim</code> since we consider
  // displacement in each space direction. The FESystem class can handle this:
  // we pass it the finite element of which we would like to compose the
  // system of, and how often it shall be repeated:

  template <int dim>
  ElasticProblem<dim>::ElasticProblem()
    : dof_handler(triangulation)
    , fe(FE_Q<dim>(1), dim)
  {}
  // In fact, the FESystem class has several more constructors which can
  // perform more complex operations than just stacking together several
  // scalar finite elements of the same type into one; we will get to know
  // these possibilities in later examples.


  // @sect4{ElasticProblem::setup_system}

  // Setting up the system of equations is identical to the function used in
  // the step-6 example. The DoFHandler class and all other classes used here
  // are fully aware that the finite element we want to use is vector-valued,
  // and take care of the vector-valuedness of the finite element
  // themselves. (In fact, they do not, but this does not need to bother you:
  // since they only need to know how many degrees of freedom there are per
  // vertex, line and cell, and they do not ask what they represent,
  // i.e. whether the finite element under consideration is vector-valued or
  // whether it is, for example, a scalar Hermite element with several degrees
  // of freedom on each vertex).
  template <int dim>
  void ElasticProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);
    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(dim),
                                             constraints);
    constraints.close();

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    /*keep_constrained_dofs = */ false);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);
  }


  // @sect4{ElasticProblem::assemble_system}

  // The big changes in this program are in the creation of matrix and right
  // hand side, since they are problem-dependent. We will go through that
  // process \step-by-step, since it is a bit more complicated than in previous
  // examples.
  //
  // The first parts of this function are the same as before, however: setting
  // up a suitable quadrature formula, initializing an FEValues object for the
  // (vector-valued) finite element we use as well as the quadrature object,
  // and declaring a number of auxiliary arrays. In addition, we declare the
  // ever same two abbreviations: <code>n_q_points</code> and
  // <code>dofs_per_cell</code>. The number of degrees of freedom per cell we
  // now obviously ask from the composed finite element rather than from the
  // underlying scalar Q1 element. Here, it is <code>dim</code> times the
  // number of degrees of freedom per cell of the Q1 element, though this is
  // not explicit knowledge we need to care about:
  template <int dim>
  void ElasticProblem<dim>::assemble_system()
  {
    QGauss<dim> quadrature_formula(fe.degree + 1);

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // As was shown in previous examples as well, we need a place where to
    // store the values of the coefficients at all the quadrature points on a
    // cell. In the present situation, we have two coefficients, lambda and
    // mu.
    std::vector<double> lambda_values(n_q_points);
    std::vector<double> mu_values(n_q_points);

    // Well, we could as well have omitted the above two arrays since we will
    // use constant coefficients for both lambda and mu, which can be declared
    // like this. They both represent functions always returning the constant
    // value 1.0. Although we could omit the respective factors in the
    // assemblage of the matrix, we use them here for purpose of
    // demonstration.
    Functions::ConstantFunction<dim> lambda(1.), mu(1.);

    // Like the two constant functions above, we will call the function
    // right_hand_side just once per cell to make things simpler.
    std::vector<Tensor<1, dim>> rhs_values(n_q_points);

    // Now we can begin with the loop over all cells:
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_matrix = 0;
        cell_rhs    = 0;

        fe_values.reinit(cell);

        // Next we get the values of the coefficients at the quadrature
        // points. Likewise for the right hand side:
        lambda.value_list(fe_values.get_quadrature_points(), lambda_values);
        mu.value_list(fe_values.get_quadrature_points(), mu_values);
        right_hand_side(fe_values.get_quadrature_points(), rhs_values);

        // Then assemble the entries of the local stiffness matrix and right
        // hand side vector. This follows almost one-to-one the pattern
        // described in the introduction of this example.  One of the few
        // comments in place is that we can compute the number
        // <code>comp(i)</code>, i.e. the index of the only nonzero vector
        // component of shape function <code>i</code> using the
        // <code>fe.system_to_component_index(i).first</code> function call
        // below.
        //
        // (By accessing the <code>first</code> variable of the return value
        // of the <code>system_to_component_index</code> function, you might
        // already have guessed that there is more in it. In fact, the
        // function returns a <code>std::pair@<unsigned int, unsigned
        // int@></code>, of which the first element is <code>comp(i)</code>
        // and the second is the value <code>base(i)</code> also noted in the
        // introduction, i.e.  the index of this shape function within all the
        // shape functions that are nonzero in this component,
        // i.e. <code>base(i)</code> in the diction of the introduction. This
        // is not a number that we are usually interested in, however.)
        //
        // With this knowledge, we can assemble the local matrix
        // contributions:
        for (const unsigned int i : fe_values.dof_indices())
          {
            const unsigned int component_i =
              fe.system_to_component_index(i).first;

            for (const unsigned int j : fe_values.dof_indices())
              {
                const unsigned int component_j =
                  fe.system_to_component_index(j).first;

                for (const unsigned int q_point :
                     fe_values.quadrature_point_indices())
                  {
                    cell_matrix(i, j) +=
                      // The first term is $\lambda \partial_i u_i, \partial_j
                      // v_j) + (\mu \partial_i u_j, \partial_j v_i)$. Note
                      // that <code>shape_grad(i,q_point)</code> returns the
                      // gradient of the only nonzero component of the i-th
                      // shape function at quadrature point q_point. The
                      // component <code>comp(i)</code> of the gradient, which
                      // is the derivative of this only nonzero vector
                      // component of the i-th shape function with respect to
                      // the comp(i)th coordinate is accessed by the appended
                      // brackets.
                      (                                                  //
                        (fe_values.shape_grad(i, q_point)[component_i] * //
                         fe_values.shape_grad(j, q_point)[component_j] * //
                         lambda_values[q_point])                         //
                        +                                                //
                        (fe_values.shape_grad(i, q_point)[component_j] * //
                         fe_values.shape_grad(j, q_point)[component_i] * //
                         mu_values[q_point])                             //
                        +                                                //
                        // The second term is $(\mu \nabla u_i, \nabla
                        // v_j)$. We need not access a specific component of
                        // the gradient, since we only have to compute the
                        // scalar product of the two gradients, of which an
                        // overloaded version of <tt>operator*</tt> takes
                        // care, as in previous examples.
                        //
                        // Note that by using the <tt>?:</tt> operator, we only
                        // do this if <tt>component_i</tt> equals
                        // <tt>component_j</tt>, otherwise a zero is added
                        // (which will be optimized away by the compiler).
                        ((component_i == component_j) ?        //
                           (fe_values.shape_grad(i, q_point) * //
                            fe_values.shape_grad(j, q_point) * //
                            mu_values[q_point]) :              //
                           0)                                  //
                        ) *                                    //
                      fe_values.JxW(q_point);                  //
                  }
              }
          }

        // Assembling the right hand side is also just as discussed in the
        // introduction:
        for (const unsigned int i : fe_values.dof_indices())
          {
            const unsigned int component_i =
              fe.system_to_component_index(i).first;

            for (const unsigned int q_point :
                 fe_values.quadrature_point_indices())
              cell_rhs(i) += fe_values.shape_value(i, q_point) *
                             rhs_values[q_point][component_i] *
                             fe_values.JxW(q_point);
          }

        // The transfer from local degrees of freedom into the global matrix
        // and right hand side vector does not depend on the equation under
        // consideration, and is thus the same as in all previous
        // examples.
        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }
  }



  // @sect4{ElasticProblem::solve}

  // The solver does not care about where the system of equations comes, as
  // long as it stays positive definite and symmetric (which are the
  // requirements for the use of the CG solver), which the system indeed
  // is. Therefore, we need not change anything.
  template <int dim>
  void ElasticProblem<dim>::solve()
  {
    SolverControl            solver_control(1000, 1e-12);
    SolverCG<Vector<double>> cg(solver_control);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve(system_matrix, solution, system_rhs, preconditioner);

    constraints.distribute(solution);
  }


  // @sect4{ElasticProblem::refine_grid}

  // The function that does the refinement of the grid is the same as in the
  // step-6 example. The quadrature formula is adapted to the linear elements
  // again. Note that the error estimator by default adds up the estimated
  // obtained from all components of the finite element solution, i.e., it
  // uses the displacement in all directions with the same weight. If we would
  // like the grid to be adapted to the x-displacement only, we could pass the
  // function an additional parameter which tells it to do so and do not
  // consider the displacements in all other directions for the error
  // indicators. However, for the current problem, it seems appropriate to
  // consider all displacement components with equal weight.
  template <int dim>
  void ElasticProblem<dim>::refine_grid()
  {
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

    KellyErrorEstimator<dim>::estimate(dof_handler,
                                       QGauss<dim - 1>(fe.degree + 1),
                                       {},
                                       solution,
                                       estimated_error_per_cell);

    GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                    estimated_error_per_cell,
                                                    0.3,
                                                    0.03);

    triangulation.execute_coarsening_and_refinement();
  }


  // @sect4{ElasticProblem::output_results}

  // The output happens mostly as has been shown in previous examples
  // already. The only difference is that the solution function is vector
  // valued. The DataOut class takes care of this automatically, but we have
  // to give each component of the solution vector a different name.
  //
  // To do this, the DataOut::add_vector() function wants a vector of
  // strings. Since the number of components is the same as the number
  // of dimensions we are working in, we use the <code>switch</code>
  // statement below.
  //
  // We note that some graphics programs have restriction on what
  // characters are allowed in the names of variables. deal.II therefore
  // supports only the minimal subset of these characters that is supported
  // by all programs. Basically, these are letters, numbers, underscores,
  // and some other characters, but in particular no whitespace and
  // minus/hyphen. The library will throw an exception otherwise, at least
  // if in debug mode.
  //
  // After listing the 1d, 2d, and 3d case, it is good style to let the
  // program die if we run upon a case which we did not consider. Remember
  // that the Assert macro generates an exception if the condition in the
  // first parameter is not satisfied. Of course, the condition
  // <code>false</code> can never be satisfied, so the program will always
  // abort whenever it gets to the default statement:
  template <int dim>
  void ElasticProblem<dim>::output_results(const unsigned int cycle) const
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);

    std::vector<std::string> solution_names;
    switch (dim)
      {
        case 1:
          solution_names.emplace_back("displacement");
          break;
        case 2:
          solution_names.emplace_back("x_displacement");
          solution_names.emplace_back("y_displacement");
          break;
        case 3:
          solution_names.emplace_back("x_displacement");
          solution_names.emplace_back("y_displacement");
          solution_names.emplace_back("z_displacement");
          break;
        default:
          Assert(false, ExcNotImplemented());
      }

    // After setting up the names for the different components of the
    // solution vector, we can add the solution vector to the list of
    // data vectors scheduled for output. Note that the following
    // function takes a vector of strings as second argument, whereas
    // the one which we have used in all previous examples accepted a
    // string there. (In fact, the function we had used before would
    // convert the single string into a vector with only one element
    // and forwards that to the other function.)
    data_out.add_data_vector(solution, solution_names);
    data_out.build_patches();

    std::ofstream output("solution-" + std::to_string(cycle) + ".vtk");
    data_out.write_vtk(output);
  }



  // @sect4{ElasticProblem::run}

  // The <code>run</code> function does the same things as in step-6, for
  // example. This time, we use the square [-1,1]^d as domain, and we refine
  // it globally four times before starting the first iteration.
  //
  // The reason for refining is a bit accidental: we use the QGauss
  // quadrature formula with two points in each direction for integration of the
  // right hand side; that means that there are four quadrature points on each
  // cell (in 2D). If we only refine the initial grid once globally, then there
  // will be only four quadrature points in each direction on the
  // domain. However, the right hand side function was chosen to be rather
  // localized and in that case, by pure chance, it happens that all quadrature
  // points lie at points where the right hand side function is zero (in
  // mathematical terms, the quadrature points happen to be at points outside
  // the <i>support</i> of the right hand side function). The right hand side
  // vector computed with quadrature will then contain only zeroes (even though
  // it would of course be nonzero if we had computed the right hand side vector
  // exactly using the integral) and the solution of the system of
  // equations is the zero vector, i.e., a finite element function that is zero
  // everywhere. In a sense, we
  // should not be surprised that this is happening since we have chosen
  // an initial grid that is totally unsuitable for the problem at hand.
  //
  // The unfortunate thing is that if the discrete solution is constant, then
  // the error indicators computed by the KellyErrorEstimator class are zero
  // for each cell as well, and the call to
  // Triangulation::refine_and_coarsen_fixed_number() will not flag any cells
  // for refinement (why should it if the indicated error is zero for each
  // cell?). The grid in the next iteration will therefore consist of four
  // cells only as well, and the same problem occurs again.
  //
  // The conclusion needs to be: while of course we will not choose the
  // initial grid to be well-suited for the accurate solution of the problem,
  // we must at least choose it such that it has the chance to capture the
  // important features of the solution. In this case, it needs to be able to
  // see the right hand side. Thus, we refine globally four times. (Any larger
  // number of global refinement steps would of course also work.)
  template <int dim>
  void ElasticProblem<dim>::run()
  {
    for (unsigned int cycle = 0; cycle < 8; ++cycle)
      {
        std::cout << "Cycle " << cycle << ':' << std::endl;

        if (cycle == 0)
          {
            GridGenerator::hyper_cube(triangulation, -1, 1);
            triangulation.refine_global(4);
          }
        else
          refine_grid();

        std::cout << "   Number of active cells:       "
                  << triangulation.n_active_cells() << std::endl;

        setup_system();

        std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
                  << std::endl;

        assemble_system();
        solve();
        output_results(cycle);
      }
  }
} // namespace Step8

// @sect3{The <code>main</code> function}

// After closing the <code>Step8</code> namespace in the last line above, the
// following is the main function of the program and is again exactly like in
// step-6 (apart from the changed class names, of course).
int main()
{
  try
    {
      Step8::ElasticProblem<2> elastic_problem_2d;
      elastic_problem_2d.run();
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
