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

// The first few files have already been covered in previous examples and will
// thus not be further commented on.
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>

// From the following include file we will import the declaration of
// H1-conforming finite element shape functions. This family of finite
// elements is called <code>FE_Q</code>, and was used in all examples before
// already to define the usual bi- or tri-linear elements, but we will now use
// it for bi-quadratic elements:
#include <deal.II/fe/fe_q.h>
// We will not read the grid from a file as in the previous example, but
// generate it using a function of the library. However, we will want to write
// out the locally refined grids (just the grid, not the solution) in each
// step, so we need the following include file instead of
// <code>grid_in.h</code>:
#include <deal.II/grid/grid_out.h>


// When using locally refined grids, we will get so-called <code>hanging
// nodes</code>. However, the standard finite element methods assumes that the
// discrete solution spaces be continuous, so we need to make sure that the
// degrees of freedom on hanging nodes conform to some constraints such that
// the global solution is continuous. We are also going to store the boundary
// conditions in this object. The following file contains a class which is
// used to handle these constraints:
#include <deal.II/lac/affine_constraints.h>

// In order to refine our grids locally, we need a function from the library
// that decides which cells to flag for refinement or coarsening based on the
// error indicators we have computed. This function is defined here:
#include <deal.II/grid/grid_refinement.h>

// Finally, we need a simple way to actually compute the refinement indicators
// based on some error estimate. While in general, adaptivity is very
// problem-specific, the error indicator in the following file often yields
// quite nicely adapted grids for a wide class of problems.
#include <deal.II/numerics/error_estimator.h>

// Finally, this is as in previous programs:
using namespace dealii;


// @sect3{The <code>Step6</code> class template}

// The main class is again almost unchanged. Two additions, however, are made:
// we have added the <code>refine_grid</code> function, which is used to
// adaptively refine the grid (instead of the global refinement in the
// previous examples), and a variable which will hold the constraints.
template <int dim>
class Step6
{
public:
  Step6();

  void run();

private:
  void setup_system();
  void assemble_system();
  void solve();
  void refine_grid();
  void output_results(const unsigned int cycle) const;

  Triangulation<dim> triangulation;

  FE_Q<dim>       fe;
  DoFHandler<dim> dof_handler;


  // This is the new variable in the main class. We need an object which holds
  // a list of constraints to hold the hanging nodes and the boundary
  // conditions.
  AffineConstraints<double> constraints;

  SparseMatrix<double> system_matrix;
  SparsityPattern      sparsity_pattern;

  Vector<double> solution;
  Vector<double> system_rhs;
};


// @sect3{Nonconstant coefficients}

// The implementation of nonconstant coefficients is copied verbatim from
// step-5:
template <int dim>
double coefficient(const Point<dim> &p)
{
  if (p.square() < 0.5 * 0.5)
    return 20;
  else
    return 1;
}



// @sect3{The <code>Step6</code> class implementation}

// @sect4{Step6::Step6}

// The constructor of this class is mostly the same as before, but this time
// we want to use the quadratic element. To do so, we only have to replace the
// constructor argument (which was <code>1</code> in all previous examples) by
// the desired polynomial degree (here <code>2</code>):
template <int dim>
Step6<dim>::Step6()
  : fe(2)
  , dof_handler(triangulation)
{}



// @sect4{Step6::setup_system}

// The next function sets up all the variables that describe the linear
// finite element problem, such as the DoFHandler, matrices, and
// vectors. The difference to what we did in step-5 is only that we now also
// have to take care of hanging node constraints. These constraints are
// handled almost exclusively by the library, i.e. you only need to know
// that they exist and how to get them, but you do not have to know how they
// are formed or what exactly is done with them.
//
// At the beginning of the function, you find all the things that are the same
// as in step-5: setting up the degrees of freedom (this time we have
// quadratic elements, but there is no difference from a user code perspective
// to the linear -- or any other degree, for that matter -- case), generating
// the sparsity pattern, and initializing the solution and right hand side
// vectors. Note that the sparsity pattern will have significantly more
// entries per row now, since there are now 9 degrees of freedom per cell
// (rather than only four), that can couple with each other.
template <int dim>
void Step6<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());

  // We may now populate the AffineConstraints object with the hanging node
  // constraints. Since we will call this function in a loop we first clear
  // the current set of constraints from the last system and then compute new
  // ones:
  constraints.clear();
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);


  // Now we are ready to interpolate the boundary values with indicator 0 (the
  // whole boundary) and store the resulting constraints in our
  // <code>constraints</code> object. Note that we do not to apply the
  // boundary conditions after assembly, like we did in earlier steps: instead
  // we put all constraints on our function space in the AffineConstraints
  // object. We can add constraints to the AffineConstraints object in either
  // order: if two constraints conflict then the constraint matrix either abort
  // or throw an exception via the Assert macro.
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           Functions::ZeroFunction<dim>(),
                                           constraints);

  // After all constraints have been added, they need to be sorted and
  // rearranged to perform some actions more efficiently. This postprocessing
  // is done using the <code>close()</code> function, after which no further
  // constraints may be added any more:
  constraints.close();

  // Now we first build our compressed sparsity pattern like we did in the
  // previous examples. Nevertheless, we do not copy it to the final sparsity
  // pattern immediately.  Note that we call a variant of
  // make_sparsity_pattern that takes the AffineConstraints object as the third
  // argument. We are letting the routine know that we will never write into
  // the locations given by <code>constraints</code> by setting the argument
  // <code>keep_constrained_dofs</code> to false (in other words, that we will
  // never write into entries of the matrix that correspond to constrained
  // degrees of freedom). If we were to condense the
  // constraints after assembling, we would have to pass <code>true</code>
  // instead because then we would first write into these locations only to
  // later set them to zero again during condensation.
  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler,
                                  dsp,
                                  constraints,
                                  /*keep_constrained_dofs = */ false);

  // Now all non-zero entries of the matrix are known (i.e. those from
  // regularly assembling the matrix and those that were introduced by
  // eliminating constraints). We may copy our intermediate object to the
  // sparsity pattern:
  sparsity_pattern.copy_from(dsp);

  // We may now, finally, initialize the sparse matrix:
  system_matrix.reinit(sparsity_pattern);
}


// @sect4{Step6::assemble_system}

// Next, we have to assemble the matrix. However, to copy the local matrix and
// vector on each cell into the global system, we are no longer using a
// hand-written loop. Instead, we use
// AffineConstraints::distribute_local_to_global() that internally executes
// this loop while performing Gaussian elimination on rows and columns
// corresponding to constrained degrees on freedom.
//
// The rest of the code that forms the local contributions remains
// unchanged. It is worth noting, however, that under the hood several things
// are different than before. First, the variable <code>dofs_per_cell</code>
// and return value of <code>quadrature_formula.size()</code> now are 9 each,
// where they were 4 before. Introducing such variables as abbreviations is a
// good strategy to make code work with different elements without having to
// change too much code. Secondly, the <code>fe_values</code> object of course
// needs to do other things as well, since the shape functions are now
// quadratic, rather than linear, in each coordinate variable. Again, however,
// this is something that is completely handled by the library.
template <int dim>
void Step6<dim>::assemble_system()
{
  const QGauss<dim> quadrature_formula(fe.degree + 1);

  FEValues<dim> fe_values(fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      cell_matrix = 0;
      cell_rhs    = 0;

      fe_values.reinit(cell);

      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
          const double current_coefficient =
            coefficient<dim>(fe_values.quadrature_point(q_index));
          for (const unsigned int i : fe_values.dof_indices())
            {
              for (const unsigned int j : fe_values.dof_indices())
                cell_matrix(i, j) +=
                  (current_coefficient *              // a(x_q)
                   fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                   fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                   fe_values.JxW(q_index));           // dx

              cell_rhs(i) += (1.0 *                               // f(x)
                              fe_values.shape_value(i, q_index) * // phi_i(x_q)
                              fe_values.JxW(q_index));            // dx
            }
        }

      // Finally, transfer the contributions from @p cell_matrix and
      // @p cell_rhs into the global objects.
      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(
        cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
    }
  // Now we are done assembling the linear system. The constraint matrix took
  // care of applying the boundary conditions and also eliminated hanging node
  // constraints. The constrained nodes are still in the linear system (there
  // is a nonzero entry, chosen in a way that the matrix is well conditioned,
  // on the diagonal of the matrix and all other entries for this line are set
  // to zero) but the computed values are invalid (i.e., the corresponding
  // entries in <code>system_rhs</code> are currently meaningless). We compute
  // the correct values for these nodes at the end of the <code>solve</code>
  // function.
}


// @sect4{Step6::solve}

// We continue with gradual improvements. The function that solves the linear
// system again uses the SSOR preconditioner, and is again unchanged except
// that we have to incorporate hanging node constraints. As mentioned above,
// the degrees of freedom from the AffineConstraints object corresponding to
// hanging node constraints and boundary values have been removed from the
// linear system by giving the rows and columns of the matrix a special
// treatment. This way, the values for these degrees of freedom have wrong,
// but well-defined values after solving the linear system. What we then have
// to do is to use the constraints to assign to them the values that they
// should have. This process, called <code>distributing</code> constraints,
// computes the values of constrained nodes from the values of the
// unconstrained ones, and requires only a single additional function call
// that you find at the end of this function:

template <int dim>
void Step6<dim>::solve()
{
  SolverControl            solver_control(1000, 1e-12);
  SolverCG<Vector<double>> solver(solver_control);

  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  solver.solve(system_matrix, solution, system_rhs, preconditioner);

  constraints.distribute(solution);
}


// @sect4{Step6::refine_grid}

// We use a sophisticated error estimation scheme to refine the mesh instead
// of global refinement. We will use the KellyErrorEstimator class which
// implements an error estimator for the Laplace equation; it can in principle
// handle variable coefficients, but we will not use these advanced features,
// but rather use its most simple form since we are not interested in
// quantitative results but only in a quick way to generate locally refined
// grids.
//
// Although the error estimator derived by Kelly et al. was originally
// developed for the Laplace equation, we have found that it is also well
// suited to quickly generate locally refined grids for a wide class of
// problems. This error estimator uses the solution gradient's jump at
// cell faces (which is a measure for the second derivatives) and
// scales it by the size of the cell. It is therefore a measure for the local
// smoothness of the solution at the place of each cell and it is thus
// understandable that it yields reasonable grids also for hyperbolic
// transport problems or the wave equation as well, although these grids are
// certainly suboptimal compared to approaches specially tailored to the
// problem. This error estimator may therefore be understood as a quick way to
// test an adaptive program.
//
// The way the estimator works is to take a <code>DoFHandler</code> object
// describing the degrees of freedom and a vector of values for each degree of
// freedom as input and compute a single indicator value for each active cell
// of the triangulation (i.e. one value for each of the active cells). To do
// so, it needs two additional pieces of information: a face quadrature formula,
// i.e., a quadrature formula on <code>dim-1</code> dimensional objects. We use
// a 3-point Gauss rule again, a choice that is consistent and appropriate with
// the bi-quadratic finite element shape functions in this program.
// (What constitutes a suitable quadrature rule here of course depends on
// knowledge of the way the error estimator evaluates the solution field. As
// said above, the jump of the gradient is integrated over each face, which
// would be a quadratic function on each face for the quadratic elements in
// use in this example. In fact, however, it is the square of the jump of the
// gradient, as explained in the documentation of that class, and that is a
// quartic function, for which a 3 point Gauss formula is sufficient since it
// integrates polynomials up to order 5 exactly.)
//
// Secondly, the function wants a list of boundary indicators for those
// boundaries where we have imposed Neumann values of the kind
// $\partial_n u(\mathbf x) = h(\mathbf x)$, along with a function $h(\mathbf
// x)$ for each such boundary. This information is represented by a map from
// boundary indicators to function objects describing the Neumann boundary
// values. In the present example program, we do not use Neumann boundary
// values, so this map is empty, and in fact constructed using the default
// constructor of the map in the place where the function call expects the
// respective function argument.
//
// The output is a vector of values for all active cells. While it may
// make sense to compute the <b>value</b> of a solution degree of freedom
// very accurately, it is usually not necessary to compute the <b>error
// indicator</b> corresponding to the solution on a cell particularly
// accurately. We therefore typically use a vector of floats instead of a vector
// of doubles to represent error indicators.
template <int dim>
void Step6<dim>::refine_grid()
{
  Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

  KellyErrorEstimator<dim>::estimate(dof_handler,
                                     QGauss<dim - 1>(fe.degree + 1),
                                     {},
                                     solution,
                                     estimated_error_per_cell);

  // The above function returned one error indicator value for each cell in
  // the <code>estimated_error_per_cell</code> array. Refinement is now done
  // as follows: refine those 30 per cent of the cells with the highest error
  // values, and coarsen the 3 per cent of cells with the lowest values.
  //
  // One can easily verify that if the second number were zero, this would
  // approximately result in a doubling of cells in each step in two space
  // dimensions, since for each of the 30 per cent of cells, four new would be
  // replaced, while the remaining 70 per cent of cells remain untouched. In
  // practice, some more cells are usually produced since it is disallowed
  // that a cell is refined twice while the neighbor cell is not refined; in
  // that case, the neighbor cell would be refined as well.
  //
  // In many applications, the number of cells to be coarsened would be set to
  // something larger than only three per cent. A non-zero value is useful
  // especially if for some reason the initial (coarse) grid is already rather
  // refined. In that case, it might be necessary to refine it in some
  // regions, while coarsening in some other regions is useful. In our case
  // here, the initial grid is very coarse, so coarsening is only necessary in
  // a few regions where over-refinement may have taken place. Thus a small,
  // non-zero value is appropriate here.
  //
  // The following function now takes these refinement indicators and flags
  // some cells of the triangulation for refinement or coarsening using the
  // method described above. It is from a class that implements several
  // different algorithms to refine a triangulation based on cell-wise error
  // indicators.
  GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                  estimated_error_per_cell,
                                                  0.3,
                                                  0.03);

  // After the previous function has exited, some cells are flagged for
  // refinement, and some other for coarsening. The refinement or coarsening
  // itself is not performed by now, however, since there are cases where
  // further modifications of these flags is useful. Here, we don't want to do
  // any such thing, so we can tell the triangulation to perform the actions
  // for which the cells are flagged:
  triangulation.execute_coarsening_and_refinement();
}


// @sect4{Step6::output_results}

// At the end of computations on each grid, and just before we continue the
// next cycle with mesh refinement, we want to output the results from this
// cycle.
//
// We have already seen in step-1 how this can be achieved for the
// mesh itself. Here, we change a few things:
// <ol>
//   <li>We use two different formats: gnuplot and VTU.</li>
//   <li>We embed the cycle number in the output file name.</li>
//   <li>For gnuplot output, we set up a GridOutFlags::Gnuplot object to
//   provide a few extra visualization arguments so that edges appear
//   curved. This is explained in further detail in step-10.</li>
// </ol>
template <int dim>
void Step6<dim>::output_results(const unsigned int cycle) const
{
  {
    GridOut               grid_out;
    std::ofstream         output("grid-" + std::to_string(cycle) + ".gnuplot");
    GridOutFlags::Gnuplot gnuplot_flags(false, 5);
    grid_out.set_flags(gnuplot_flags);
    MappingQGeneric<dim> mapping(3);
    grid_out.write_gnuplot(triangulation, output, &mapping);
  }

  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "solution");
    data_out.build_patches();

    std::ofstream output("solution-" + std::to_string(cycle) + ".vtu");
    data_out.write_vtu(output);
  }
}


// @sect4{Step6::run}

// The final function before <code>main()</code> is again the main driver of
// the class, <code>run()</code>. It is similar to the one of step-5, except
// that we generate a file in the program again instead of reading it from
// disk, in that we adaptively instead of globally refine the mesh, and that
// we output the solution on the final mesh in the present function.
//
// The first block in the main loop of the function deals with mesh generation.
// If this is the first cycle of the program, instead of reading the grid from
// a file on disk as in the previous example, we now again create it using a
// library function. The domain is again a circle with center at the origin and
// a radius of one (these are the two hidden arguments to the function, which
// have default values).
//
// You will notice by looking at the coarse grid that it is of inferior
// quality than the one which we read from the file in the previous example:
// the cells are less equally formed. However, using the library function this
// program works in any space dimension, which was not the case before.
//
// In case we find that this is not the first cycle, we want to refine the
// grid. Unlike the global refinement employed in the last example program, we
// now use the adaptive procedure described above.
//
// The rest of the loop looks as before:
template <int dim>
void Step6<dim>::run()
{
  for (unsigned int cycle = 0; cycle < 8; ++cycle)
    {
      std::cout << "Cycle " << cycle << ':' << std::endl;

      if (cycle == 0)
        {
          GridGenerator::hyper_ball(triangulation);
          triangulation.refine_global(1);
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


// @sect3{The <code>main</code> function}

// The main function is unaltered in its functionality from the previous
// example, but we have taken a step of additional caution. Sometimes,
// something goes wrong (such as insufficient disk space upon writing an
// output file, not enough memory when trying to allocate a vector or a
// matrix, or if we can't read from or write to a file for whatever reason),
// and in these cases the library will throw exceptions. Since these are
// run-time problems, not programming errors that can be fixed once and for
// all, this kind of exceptions is not switched off in optimized mode, in
// contrast to the <code>Assert</code> macro which we have used to test
// against programming errors. If uncaught, these exceptions propagate the
// call tree up to the <code>main</code> function, and if they are not caught
// there either, the program is aborted. In many cases, like if there is not
// enough memory or disk space, we can't do anything but we can at least print
// some text trying to explain the reason why the program failed. A way to do
// so is shown in the following. It is certainly useful to write any larger
// program in this way, and you can do so by more or less copying this
// function except for the <code>try</code> block that actually encodes the
// functionality particular to the present application.
int main()
{
  // The general idea behind the layout of this function is as follows: let's
  // try to run the program as we did before...
  try
    {
      Step6<2> laplace_problem_2d;
      laplace_problem_2d.run();
    }
  // ...and if this should fail, try to gather as much information as
  // possible. Specifically, if the exception that was thrown is an object of
  // a class that is derived from the C++ standard class
  // <code>exception</code>, then we can use the <code>what</code> member
  // function to get a string which describes the reason why the exception was
  // thrown.
  //
  // The deal.II exception classes are all derived from the standard class,
  // and in particular, the <code>exc.what()</code> function will return
  // approximately the same string as would be generated if the exception was
  // thrown using the <code>Assert</code> macro. You have seen the output of
  // such an exception in the previous example, and you then know that it
  // contains the file and line number of where the exception occurred, and
  // some other information. This is also what the following statements would
  // print.
  //
  // Apart from this, there isn't much that we can do except exiting the
  // program with an error code (this is what the <code>return 1;</code>
  // does):
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
  // If the exception that was thrown somewhere was not an object of a class
  // derived from the standard <code>exception</code> class, then we can't do
  // anything at all. We then simply print an error message and exit.
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

  // If we got to this point, there was no exception which propagated up to
  // the main function (there may have been exceptions, but they were caught
  // somewhere in the program or the library). Therefore, the program
  // performed as was expected and we can return without error.
  return 0;
}
