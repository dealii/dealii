/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2000 - 2013 by the deal.II authors
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
 * Author: Wolfgang Bangerth, University of Heidelberg, 2000
 */


// @sect3{Include files}

// As usual, the first few include files are already known, so we will not
// comment on them further.
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
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
  // variable: Instead of a concrete finite element class such as
  // <code>FE_Q</code>, we now use a more generic one,
  // <code>FESystem</code>. In fact, <code>FESystem</code> is not really a
  // finite element itself in that it does not implement shape functions of
  // its own.  Rather, it is a class that can be used to stack several other
  // elements together to form one vector-valued finite element. In our case,
  // we will compose the vector-valued element of <code>FE_Q(1)</code>
  // objects, as shown below in the constructor of this class.
  template <int dim>
  class ElasticProblem
  {
  public:
    ElasticProblem ();
    ~ElasticProblem ();
    void run ();

  private:
    void setup_system ();
    void assemble_system ();
    void solve ();
    void refine_grid ();
    void output_results (const unsigned int cycle) const;

    Triangulation<dim>   triangulation;
    DoFHandler<dim>      dof_handler;

    FESystem<dim>        fe;

    ConstraintMatrix     hanging_node_constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       system_rhs;
  };


  // @sect3{Right hand side values}

  // Before going over to the implementation of the main class, we declare and
  // define the class which describes the right hand side. This time, the
  // right hand side is vector-valued, as is the solution, so we will describe
  // the changes required for this in some more detail.
  //
  // The first thing is that vector-valued functions have to have a
  // constructor, since they need to pass down to the base class of how many
  // components the function consists. The default value in the constructor of
  // the base class is one (i.e.: a scalar function), which is why we did not
  // need not define a constructor for the scalar function used in previous
  // programs.
  template <int dim>
  class RightHandSide :  public Function<dim>
  {
  public:
    RightHandSide ();

    // The next change is that we want a replacement for the
    // <code>value</code> function of the previous examples. There, a second
    // parameter <code>component</code> was given, which denoted which
    // component was requested. Here, we implement a function that returns the
    // whole vector of values at the given place at once, in the second
    // argument of the function. The obvious name for such a replacement
    // function is <code>vector_value</code>.
    //
    // Secondly, in analogy to the <code>value_list</code> function, there is
    // a function <code>vector_value_list</code>, which returns the values of
    // the vector-valued function at several points at once:
    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &values) const;

    virtual void vector_value_list (const std::vector<Point<dim> > &points,
                                    std::vector<Vector<double> >   &value_list) const;
  };


  // This is the constructor of the right hand side class. As said above, it
  // only passes down to the base class the number of components, which is
  // <code>dim</code> in the present case (one force component in each of the
  // <code>dim</code> space directions).
  //
  // Some people would have moved the definition of such a short function
  // right into the class declaration. We do not do that, as a matter of
  // style: the deal.II style guides require that class declarations contain
  // only declarations, and that definitions are always to be found
  // outside. This is, obviously, as much as matter of taste as indentation,
  // but we try to be consistent in this direction.
  template <int dim>
  RightHandSide<dim>::RightHandSide ()
    :
    Function<dim> (dim)
  {}


  // Next the function that returns the whole vector of values at the point
  // <code>p</code> at once.
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
  // the correct size; this is in contract to the <code>Assert</code> call
  // that is completely removed if the program is compiled in optimized mode.
  //
  // Likewise, if by some accident someone tried to compile and run the
  // program in only one space dimension (in which the elastic equations do
  // not make much sense since they reduce to the ordinary Laplace equation),
  // we terminate the program in the second assertion. The program will work
  // just fine in 3d, however.
  template <int dim>
  inline
  void RightHandSide<dim>::vector_value (const Point<dim> &p,
                                         Vector<double>   &values) const
  {
    Assert (values.size() == dim,
            ExcDimensionMismatch (values.size(), dim));
    Assert (dim >= 2, ExcNotImplemented());

    // The rest of the function implements computing force values. We will use
    // a constant (unit) force in x-direction located in two little circles
    // (or spheres, in 3d) around points (0.5,0) and (-0.5,0), and y-force in
    // an area around the origin; in 3d, the z-component of these centers is
    // zero as well.
    //
    // For this, let us first define two objects that denote the centers of
    // these areas. Note that upon construction of the <code>Point</code>
    // objects, all components are set to zero.
    Point<dim> point_1, point_2;
    point_1(0) = 0.5;
    point_2(0) = -0.5;

    // If now the point <code>p</code> is in a circle (sphere) of radius 0.2
    // around one of these points, then set the force in x-direction to one,
    // otherwise to zero:
    if (((p-point_1).square() < 0.2*0.2) ||
        ((p-point_2).square() < 0.2*0.2))
      values(0) = 1;
    else
      values(0) = 0;

    // Likewise, if <code>p</code> is in the vicinity of the origin, then set
    // the y-force to 1, otherwise to zero:
    if (p.square() < 0.2*0.2)
      values(1) = 1;
    else
      values(1) = 0;
  }



  // Now, this is the function of the right hand side class that returns the
  // values at several points at once. The function starts out with checking
  // that the number of input and output arguments is equal (the sizes of the
  // individual output vectors will be checked in the function that we call
  // further down below). Next, we define an abbreviation for the number of
  // points which we shall work on, to make some things simpler below.
  template <int dim>
  void RightHandSide<dim>::vector_value_list (const std::vector<Point<dim> > &points,
                                              std::vector<Vector<double> >   &value_list) const
  {
    Assert (value_list.size() == points.size(),
            ExcDimensionMismatch (value_list.size(), points.size()));

    const unsigned int n_points = points.size();

    // Finally we treat each of the points. In one of the previous examples,
    // we have explained why the
    // <code>value_list</code>/<code>vector_value_list</code> function had
    // been introduced: to prevent us from calling virtual functions too
    // frequently. On the other hand, we now need to implement the same
    // function twice, which can lead to confusion if one function is changed
    // but the other is not.
    //
    // We can prevent this situation by calling
    // <code>RightHandSide::vector_value</code> on each point in the input
    // list. Note that by giving the full name of the function, including the
    // class name, we instruct the compiler to explicitly call this function,
    // and not to use the virtual function call mechanism that would be used
    // if we had just called <code>vector_value</code>. This is important,
    // since the compiler generally can't make any assumptions which function
    // is called when using virtual functions, and it therefore can't inline
    // the called function into the site of the call. On the contrary, here we
    // give the fully qualified name, which bypasses the virtual function
    // call, and consequently the compiler knows exactly which function is
    // called and will inline above function into the present location. (Note
    // that we have declared the <code>vector_value</code> function above
    // <code>inline</code>, though modern compilers are also able to inline
    // functions even if they have not been declared as inline).
    //
    // It is worth noting why we go to such length explaining what we
    // do. Using this construct, we manage to avoid any inconsistency: if we
    // want to change the right hand side function, it would be difficult to
    // always remember that we always have to change two functions in the same
    // way. Using this forwarding mechanism, we only have to change a single
    // place (the <code>vector_value</code> function), and the second place
    // (the <code>vector_value_list</code> function) will always be consistent
    // with it. At the same time, using virtual function call bypassing, the
    // code is no less efficient than if we had written it twice in the first
    // place:
    for (unsigned int p=0; p<n_points; ++p)
      RightHandSide<dim>::vector_value (points[p],
                                        value_list[p]);
  }



  // @sect3{The <code>ElasticProblem</code> class implementation}

  // @sect4{ElasticProblem::ElasticProblem}

  // Following is the constructor of the main class. As said before, we would
  // like to construct a vector-valued finite element that is composed of
  // several scalar finite elements (i.e., we want to build the vector-valued
  // element so that each of its vector components consists of the shape
  // functions of a scalar element). Of course, the number of scalar finite
  // elements we would like to stack together equals the number of components
  // the solution function has, which is <code>dim</code> since we consider
  // displacement in each space direction. The <code>FESystem</code> class can
  // handle this: we pass it the finite element of which we would like to
  // compose the system of, and how often it shall be repeated:

  template <int dim>
  ElasticProblem<dim>::ElasticProblem ()
    :
    dof_handler (triangulation),
    fe (FE_Q<dim>(1), dim)
  {}
  // In fact, the <code>FESystem</code> class has several more constructors
  // which can perform more complex operations than just stacking together
  // several scalar finite elements of the same type into one; we will get to
  // know these possibilities in later examples.



  // @sect4{ElasticProblem::~ElasticProblem}

  // The destructor, on the other hand, is exactly as in step-6:
  template <int dim>
  ElasticProblem<dim>::~ElasticProblem ()
  {
    dof_handler.clear ();
  }


  // @sect4{ElasticProblem::setup_system}

  // Setting up the system of equations is identical to the function used in
  // the step-6 example. The <code>DoFHandler</code> class and all other
  // classes used here are fully aware that the finite element we want to use
  // is vector-valued, and take care of the vector-valuedness of the finite
  // element themselves. (In fact, they do not, but this does not need to
  // bother you: since they only need to know how many degrees of freedom
  // there are per vertex, line and cell, and they do not ask what they
  // represent, i.e. whether the finite element under consideration is
  // vector-valued or whether it is, for example, a scalar Hermite element
  // with several degrees of freedom on each vertex).
  template <int dim>
  void ElasticProblem<dim>::setup_system ()
  {
    dof_handler.distribute_dofs (fe);
    hanging_node_constraints.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler,
                                             hanging_node_constraints);
    hanging_node_constraints.close ();
    sparsity_pattern.reinit (dof_handler.n_dofs(),
                             dof_handler.n_dofs(),
                             dof_handler.max_couplings_between_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);

    hanging_node_constraints.condense (sparsity_pattern);

    sparsity_pattern.compress();

    system_matrix.reinit (sparsity_pattern);

    solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());
  }


  // @sect4{ElasticProblem::assemble_system}

  // The big changes in this program are in the creation of matrix and right
  // hand side, since they are problem-dependent. We will go through that
  // process step-by-step, since it is a bit more complicated than in previous
  // examples.
  //
  // The first parts of this function are the same as before, however: setting
  // up a suitable quadrature formula, initializing an <code>FEValues</code>
  // object for the (vector-valued) finite element we use as well as the
  // quadrature object, and declaring a number of auxiliary arrays. In
  // addition, we declare the ever same two abbreviations:
  // <code>n_q_points</code> and <code>dofs_per_cell</code>. The number of
  // degrees of freedom per cell we now obviously ask from the composed finite
  // element rather than from the underlying scalar Q1 element. Here, it is
  // <code>dim</code> times the number of degrees of freedom per cell of the
  // Q1 element, though this is not explicit knowledge we need to care about:
  template <int dim>
  void ElasticProblem<dim>::assemble_system ()
  {
    QGauss<dim>  quadrature_formula(2);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values   | update_gradients |
                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       cell_rhs (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    // As was shown in previous examples as well, we need a place where to
    // store the values of the coefficients at all the quadrature points on a
    // cell. In the present situation, we have two coefficients, lambda and
    // mu.
    std::vector<double>     lambda_values (n_q_points);
    std::vector<double>     mu_values (n_q_points);

    // Well, we could as well have omitted the above two arrays since we will
    // use constant coefficients for both lambda and mu, which can be declared
    // like this. They both represent functions always returning the constant
    // value 1.0. Although we could omit the respective factors in the
    // assemblage of the matrix, we use them here for purpose of
    // demonstration.
    ConstantFunction<dim> lambda(1.), mu(1.);

    // Then again, we need to have the same for the right hand side. This is
    // exactly as before in previous examples. However, we now have a
    // vector-valued right hand side, which is why the data type of the
    // <code>rhs_values</code> array is changed. We initialize it by
    // <code>n_q_points</code> elements, each of which is a
    // <code>Vector@<double@></code> with <code>dim</code> elements.
    RightHandSide<dim>      right_hand_side;
    std::vector<Vector<double> > rhs_values (n_q_points,
                                             Vector<double>(dim));


    // Now we can begin with the loop over all cells:
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        cell_matrix = 0;
        cell_rhs = 0;

        fe_values.reinit (cell);

        // Next we get the values of the coefficients at the quadrature
        // points. Likewise for the right hand side:
        lambda.value_list (fe_values.get_quadrature_points(), lambda_values);
        mu.value_list     (fe_values.get_quadrature_points(), mu_values);

        right_hand_side.vector_value_list (fe_values.get_quadrature_points(),
                                           rhs_values);

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
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const unsigned int
            component_i = fe.system_to_component_index(i).first;

            for (unsigned int j=0; j<dofs_per_cell; ++j)
              {
                const unsigned int
                component_j = fe.system_to_component_index(j).first;

                for (unsigned int q_point=0; q_point<n_q_points;
                     ++q_point)
                  {
                    cell_matrix(i,j)
                    +=
                      // The first term is (lambda d_i u_i, d_j v_j) + (mu d_i
                      // u_j, d_j v_i).  Note that
                      // <code>shape_grad(i,q_point)</code> returns the
                      // gradient of the only nonzero component of the i-th
                      // shape function at quadrature point q_point. The
                      // component <code>comp(i)</code> of the gradient, which
                      // is the derivative of this only nonzero vector
                      // component of the i-th shape function with respect to
                      // the comp(i)th coordinate is accessed by the appended
                      // brackets.
                      (
                        (fe_values.shape_grad(i,q_point)[component_i] *
                         fe_values.shape_grad(j,q_point)[component_j] *
                         lambda_values[q_point])
                        +
                        (fe_values.shape_grad(i,q_point)[component_j] *
                         fe_values.shape_grad(j,q_point)[component_i] *
                         mu_values[q_point])
                        +
                        // The second term is (mu nabla u_i, nabla v_j).  We
                        // need not access a specific component of the
                        // gradient, since we only have to compute the scalar
                        // product of the two gradients, of which an
                        // overloaded version of the operator* takes care, as
                        // in previous examples.
                        //
                        // Note that by using the ?: operator, we only do this
                        // if comp(i) equals comp(j), otherwise a zero is
                        // added (which will be optimized away by the
                        // compiler).
                        ((component_i == component_j) ?
                         (fe_values.shape_grad(i,q_point) *
                          fe_values.shape_grad(j,q_point) *
                          mu_values[q_point])  :
                         0)
                      )
                      *
                      fe_values.JxW(q_point);
                  }
              }
          }

        // Assembling the right hand side is also just as discussed in the
        // introduction:
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const unsigned int
            component_i = fe.system_to_component_index(i).first;

            for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
              cell_rhs(i) += fe_values.shape_value(i,q_point) *
                             rhs_values[q_point](component_i) *
                             fe_values.JxW(q_point);
          }

        // The transfer from local degrees of freedom into the global matrix
        // and right hand side vector does not depend on the equation under
        // consideration, and is thus the same as in all previous
        // examples. The same holds for the elimination of hanging nodes from
        // the matrix and right hand side, once we are done with assembling
        // the entire linear system:
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

    hanging_node_constraints.condense (system_matrix);
    hanging_node_constraints.condense (system_rhs);

    // The interpolation of the boundary values needs a small modification:
    // since the solution function is vector-valued, so need to be the
    // boundary values. The <code>ZeroFunction</code> constructor accepts a
    // parameter that tells it that it shall represent a vector valued,
    // constant zero function with that many components. By default, this
    // parameter is equal to one, in which case the <code>ZeroFunction</code>
    // object would represent a scalar function. Since the solution vector has
    // <code>dim</code> components, we need to pass <code>dim</code> as number
    // of components to the zero function as well.
    std::map<types::global_dof_index,double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler,
                                              0,
                                              ZeroFunction<dim>(dim),
                                              boundary_values);
    MatrixTools::apply_boundary_values (boundary_values,
                                        system_matrix,
                                        solution,
                                        system_rhs);
  }



  // @sect4{ElasticProblem::solve}

  // The solver does not care about where the system of equations comes, as
  // long as it stays positive definite and symmetric (which are the
  // requirements for the use of the CG solver), which the system indeed
  // is. Therefore, we need not change anything.
  template <int dim>
  void ElasticProblem<dim>::solve ()
  {
    SolverControl           solver_control (1000, 1e-12);
    SolverCG<>              cg (solver_control);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve (system_matrix, solution, system_rhs,
              preconditioner);

    hanging_node_constraints.distribute (solution);
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
  void ElasticProblem<dim>::refine_grid ()
  {
    Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

    KellyErrorEstimator<dim>::estimate (dof_handler,
                                        QGauss<dim-1>(2),
                                        typename FunctionMap<dim>::type(),
                                        solution,
                                        estimated_error_per_cell);

    GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                     estimated_error_per_cell,
                                                     0.3, 0.03);

    triangulation.execute_coarsening_and_refinement ();
  }


  // @sect4{ElasticProblem::output_results}

  // The output happens mostly as has been shown in previous examples
  // already. The only difference is that the solution function is vector
  // valued. The <code>DataOut</code> class takes care of this automatically,
  // but we have to give each component of the solution vector a different
  // name.
  template <int dim>
  void ElasticProblem<dim>::output_results (const unsigned int cycle) const
  {
    std::string filename = "solution-";
    filename += ('0' + cycle);
    Assert (cycle < 10, ExcInternalError());

    filename += ".vtk";
    std::ofstream output (filename.c_str());

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);



    // As said above, we need a different name for each component of the
    // solution function. To pass one name for each component, a vector of
    // strings is used. Since the number of components is the same as the
    // number of dimensions we are working in, the following
    // <code>switch</code> statement is used.
    //
    // We note that some graphics programs have restriction as to what
    // characters are allowed in the names of variables. The library therefore
    // supports only the minimal subset of these characters that is supported
    // by all programs. Basically, these are letters, numbers, underscores,
    // and some other characters, but in particular no whitespace and
    // minus/hyphen. The library will throw an exception otherwise, at least
    // if in debug mode.
    //
    // After listing the 1d, 2d, and 3d case, it is good style to let the
    // program die if we run upon a case which we did not consider. Remember
    // that the <code>Assert</code> macro generates an exception if the
    // condition in the first parameter is not satisfied. Of course, the
    // condition <code>false</code> can never be satisfied, so the program
    // will always abort whenever it gets to the default statement:
    std::vector<std::string> solution_names;
    switch (dim)
      {
      case 1:
        solution_names.push_back ("displacement");
        break;
      case 2:
        solution_names.push_back ("x_displacement");
        solution_names.push_back ("y_displacement");
        break;
      case 3:
        solution_names.push_back ("x_displacement");
        solution_names.push_back ("y_displacement");
        solution_names.push_back ("z_displacement");
        break;
      default:
        Assert (false, ExcNotImplemented());
      }

    // After setting up the names for the different components of the solution
    // vector, we can add the solution vector to the list of data vectors
    // scheduled for output. Note that the following function takes a vector
    // of strings as second argument, whereas the one which we have used in
    // all previous examples accepted a string there. In fact, the latter
    // function is only a shortcut for the function which we call here: it
    // puts the single string that is passed to it into a vector of strings
    // with only one element and forwards that to the other function.
    data_out.add_data_vector (solution, solution_names);
    data_out.build_patches ();
    data_out.write_vtk (output);
  }



  // @sect4{ElasticProblem::run}

  // The <code>run</code> function does the same things as in step-6, for
  // example. This time, we use the square [-1,1]^d as domain, and we refine
  // it twice globally before starting the first iteration.
  //
  // The reason is the following: we use the <code>Gauss</code> quadrature
  // formula with two points in each direction for integration of the right
  // hand side; that means that there are four quadrature points on each cell
  // (in 2D). If we only refine the initial grid once globally, then there
  // will be only four quadrature points in each direction on the
  // domain. However, the right hand side function was chosen to be rather
  // localized and in that case all quadrature points lie outside the support
  // of the right hand side function. The right hand side vector will then
  // contain only zeroes and the solution of the system of equations is the
  // zero vector, i.e. a finite element function that it zero everywhere. We
  // should not be surprised about such things happening, since we have chosen
  // an initial grid that is totally unsuitable for the problem at hand.
  //
  // The unfortunate thing is that if the discrete solution is constant, then
  // the error indicators computed by the <code>KellyErrorEstimator</code>
  // class are zero for each cell as well, and the call to
  // <code>refine_and_coarsen_fixed_number</code> on the
  // <code>triangulation</code> object will not flag any cells for refinement
  // (why should it if the indicated error is zero for each cell?). The grid
  // in the next iteration will therefore consist of four cells only as well,
  // and the same problem occurs again.
  //
  // The conclusion needs to be: while of course we will not choose the
  // initial grid to be well-suited for the accurate solution of the problem,
  // we must at least choose it such that it has the chance to capture the
  // most striking features of the solution. In this case, it needs to be able
  // to see the right hand side. Thus, we refine twice globally. (Note that
  // the <code>refine_global</code> function is not part of the
  // <code>GridRefinement</code> class in which
  // <code>refine_and_coarsen_fixed_number</code> is declared, for
  // example. The reason is first that it is not an algorithm that computed
  // refinement flags from indicators, but more importantly that it actually
  // performs the refinement, in contrast to the functions in
  // <code>GridRefinement</code> that only flag cells without actually
  // refining the grid.)
  template <int dim>
  void ElasticProblem<dim>::run ()
  {
    for (unsigned int cycle=0; cycle<8; ++cycle)
      {
        std::cout << "Cycle " << cycle << ':' << std::endl;

        if (cycle == 0)
          {
            GridGenerator::hyper_cube (triangulation, -1, 1);
            triangulation.refine_global (2);
          }
        else
          refine_grid ();

        std::cout << "   Number of active cells:       "
                  << triangulation.n_active_cells()
                  << std::endl;

        setup_system ();

        std::cout << "   Number of degrees of freedom: "
                  << dof_handler.n_dofs()
                  << std::endl;

        assemble_system ();
        solve ();
        output_results (cycle);
      }
  }
}

// @sect3{The <code>main</code> function}

// After closing the <code>Step8</code> namespace in the last line above, the
// following is the main function of the program and is again exactly like in
// step-6 (apart from the changed class names, of course).
int main ()
{
  try
    {
      dealii::deallog.depth_console (0);

      Step8::ElasticProblem<2> elastic_problem_2d;
      elastic_problem_2d.run ();
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
