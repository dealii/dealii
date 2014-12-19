/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2013 by the deal.II authors
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
 * Author: Wolfgang Bangerth, University of Heidelberg, 1999
 */


// @sect3{Include files}

// Again, the first few include files are already known, so we won't comment
// on them:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

// This one is new. We want to read a triangulation from disk, and the class
// which does this is declared in the following file:
#include <deal.II/grid/grid_in.h>

// We will use a circular domain, and the object describing the boundary of it
// comes from this file:
#include <deal.II/grid/tria_boundary_lib.h>

// This is C++ ...
#include <fstream>
// ... and this is too: We will convert integers to strings using the C++
// stringstream class <code>ostringstream</code>:
#include <sstream>

// Finally, this has been discussed in previous tutorial programs before:
using namespace dealii;


// @sect3{The <code>Step5</code> class template}

// The main class is mostly as in the previous example. The most visible
// change is that the function <code>make_grid_and_dofs</code> has been
// removed, since creating the grid is now done in the <code>run</code>
// function and the rest of its functionality is now in
// <code>setup_system</code>. Apart from this, everything is as before.
template <int dim>
class Step5
{
public:
  Step5 ();
  void run ();

private:
  void setup_system ();
  void assemble_system ();
  void solve ();
  void output_results (const unsigned int cycle) const;

  Triangulation<dim>   triangulation;
  FE_Q<dim>            fe;
  DoFHandler<dim>      dof_handler;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double>       solution;
  Vector<double>       system_rhs;
};


// @sect3{Nonconstant coefficients, using <code>Assert</code>}

// In step-4, we showed how to use non-constant boundary values and right hand
// side.  In this example, we want to use a variable coefficient in the
// elliptic operator instead. Of course, the suitable object is a
// <code>Function</code>, as we have used for the right hand side and boundary
// values in the last example. We will use it again, but we implement another
// function <code>value_list</code> which takes a list of points and returns
// the values of the function at these points as a list. The reason why such a
// function is reasonable although we can get all the information from the
// <code>value</code> function as well will be explained below when assembling
// the matrix.
//
// The need to declare a seemingly useless default constructor exists here
// just as in the previous example.
template <int dim>
class Coefficient : public Function<dim>
{
public:
  Coefficient ()  : Function<dim>() {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;

  virtual void value_list (const std::vector<Point<dim> > &points,
                           std::vector<double>            &values,
                           const unsigned int              component = 0) const;
};



// This is the implementation of the coefficient function for a single
// point. We let it return 20 if the distance to the origin is less than 0.5,
// and 1 otherwise. As in the previous example, we simply ignore the second
// parameter of the function that is used to denote different components of
// vector-valued functions (we deal only with a scalar function here, after
// all):
template <int dim>
double Coefficient<dim>::value (const Point<dim> &p,
                                const unsigned int /*component*/) const
{
  if (p.square() < 0.5*0.5)
    return 20;
  else
    return 1;
}



// And this is the function that returns the value of the coefficient at a
// whole list of points at once. Of course, we need to make sure that the
// values are the same as if we would ask the <code>value</code> function for
// each point individually.
//
// This method takes three parameters: a list of points at which to evaluate
// the function, a list that will hold the values at these points, and the
// vector component that should be zero here since we only have a single
// scalar function.  Now, of course the size of the output array
// (<code>values</code>) must be the same as that of the input array
// (<code>points</code>), and we could simply assume that. However, in
// practice, it turns out that more than 90 per cent of programming errors are
// invalid function parameters such as invalid array sizes, etc, so we should
// try to make sure that the parameters are valid. For this, the
// <code>Assert</code> macro is a good means, since it makes sure that the
// condition which is given as first argument is valid, and if not throws an
// exception (its second argument) which will usually terminate the program
// giving information where the error occurred and what the reason was. This
// generally reduces the time to find programming errors dramatically and we
// have found assertions an invaluable means to program fast.
//
// On the other hand, all these checks (there are more than 4200 of them in
// the library at present) should not slow down the program too much if you
// want to do large computations. To this end, the <code>Assert</code> macro
// is only used in debug mode and expands to nothing if in optimized
// mode. Therefore, while you test your program on small problems and debug
// it, the assertions will tell you where the problems are.  Once your program
// is stable, you can switch off debugging and the program will run your real
// computations without the assertions and at maximum speed. (In fact, it
// turns out the switching off all the checks in the library that prevent you
// from calling functions with the wrong arguments by switching to optimized
// mode, makes most programs run faster by about a factor of four. This
// should, however, not try to induce you to always run in optimized mode:
// Most people who have tried that soon realize that they introduce lots of
// errors that would have easily been caught had they run the program in debug
// mode while developing.) For those who want to try: The way to switch from
// debug mode to optimized mode is to go edit the Makefile in this
// directory. It should have a line <code>debug-mode = on</code>; simply
// replace it by <code>debug-mode = off</code> and recompile your program. The
// output of the <code>make</code> program should already indicate to you that
// the program is now compiled in optimized mode, and it will later also be
// linked to libraries that have been compiled for optimized mode.
//
// Here, as has been said above, we would like to make sure that the size of
// the two arrays is equal, and if not throw an exception. Comparing the sizes
// of two arrays is one of the most frequent checks, which is why there is
// already an exception class <code>ExcDimensionMismatch</code> that takes the
// sizes of two vectors and prints some output in case the condition is
// violated:

template <int dim>
void Coefficient<dim>::value_list (const std::vector<Point<dim> > &points,
                                   std::vector<double>            &values,
                                   const unsigned int              component) const
{
  Assert (values.size() == points.size(),
          ExcDimensionMismatch (values.size(), points.size()));
  // Since examples are not very good if they do not demonstrate their point,
  // we will show how to trigger this exception at the end of the main
  // program, and what output results from this (see the <code>Results</code>
  // section of this example program). You will certainly notice that the
  // output is quite well suited to quickly find what the problem is and what
  // parameters are expected. An additional plus is that if the program is run
  // inside a debugger, it will stop at the point where the exception is
  // triggered, so you can go up the call stack to immediately find the place
  // where the the array with the wrong size was set up.

  // While we're at it, we can do another check: the coefficient is a scalar,
  // but the <code>Function</code> class also represents vector-valued
  // function. A scalar function must therefore be considered as a
  // vector-valued function with only one component, so the only valid
  // component for which a user might ask is zero (we always count from
  // zero). The following assertion checks this. If the condition in the
  // <code>Assert</code> call is violated, an exception of type
  // <code>ExcRange</code> will be triggered; that class takes the violating
  // index as first argument, and the second and third arguments denote a
  // range that includes the left point but is open at the right, i.e. here
  // the interval [0,1). For integer arguments, this means that the only value
  // in the range is the zero, of course. (The interval is half open since we
  // also want to write exceptions like <code>ExcRange(i,0,v.size())</code>,
  // where an index must be between zero but less than the size of an
  // array. To save us the effort of writing <code>v.size()-1</code> in many
  // places, the range is defined as half-open.)
  Assert (component == 0,
          ExcIndexRange (component, 0, 1));

  // The rest of the function is uneventful: we define <code>n_q_points</code>
  // as an abbreviation for the number of points for which function values are
  // requested, and then simply fill the output value:
  const unsigned int n_points = points.size();

  for (unsigned int i=0; i<n_points; ++i)
    {
      if (points[i].square() < 0.5*0.5)
        values[i] = 20;
      else
        values[i] = 1;
    }
}


// @sect3{The <code>Step5</code> class implementation}

// @sect4{Step5::Step5}

// This function is as before.
template <int dim>
Step5<dim>::Step5 () :
  fe (1),
  dof_handler (triangulation)
{}



// @sect4{Step5::setup_system}

// This is the function <code>make_grid_and_dofs</code> from the previous
// example, minus the generation of the grid. Everything else is unchanged:
template <int dim>
void Step5<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);

  std::cout << "   Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << std::endl;

  CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);
  sparsity_pattern.copy_from(c_sparsity);

  system_matrix.reinit (sparsity_pattern);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
}



// @sect4{Step5::assemble_system}

// As in the previous examples, this function is not changed much with regard
// to its functionality, but there are still some optimizations which we will
// show. For this, it is important to note that if efficient solvers are used
// (such as the preconditions CG method), assembling the matrix and right hand
// side can take a comparable time, and you should think about using one or
// two optimizations at some places.
//
// What we will show here is how we can avoid calls to the shape_value,
// shape_grad, and quadrature_point functions of the FEValues object, and in
// particular optimize away most of the virtual function calls of the Function
// object. The way to do so will be explained in the following, while those
// parts of this function that are not changed with respect to the previous
// example are not commented on.
//
// The first parts of the function are completely unchanged from before:
template <int dim>
void Step5<dim>::assemble_system ()
{
  QGauss<dim>  quadrature_formula(2);

  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values    |  update_gradients |
                           update_quadrature_points  |  update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  // Here is one difference: for this program, we will again use a constant
  // right hand side function and zero boundary values, but a variable
  // coefficient. We have already declared the class that represents this
  // coefficient above, so we only have to declare a corresponding object
  // here.
  //
  // Then, below, we will ask the <code>coefficient</code> function object to
  // compute the values of the coefficient at all quadrature points on one
  // cell at once. The reason for this is that, if you look back at how we did
  // this in step-4, you will realize that we called the function computing
  // the right hand side value inside nested loops over all degrees of freedom
  // and over all quadrature points, i.e. dofs_per_cell*n_q_points times. For
  // the coefficient that is used inside the matrix, this would actually be
  // dofs_per_cell*dofs_per_cell*n_q_points. On the other hand, the function
  // will of course return the same value every time it is called with the
  // same quadrature point, independently of what shape function we presently
  // treat; secondly, these are virtual function calls, so are rather
  // expensive. Obviously, there are only n_q_point different values, and we
  // shouldn't call the function more often than that. Or, even better than
  // this, compute all of these values at once, and get away with a single
  // function call per cell.
  //
  // This is exactly what we are going to do. For this, we need some space to
  // store the values in. We therefore also have to declare an array to hold
  // these values:
  const Coefficient<dim> coefficient;
  std::vector<double>    coefficient_values (n_q_points);

  // Next is the typical loop over all cells to compute local contributions
  // and then to transfer them into the global matrix and vector.
  //
  // The only two things in which this loop differs from step-4 is that we
  // want to compute the value of the coefficient in all quadrature points on
  // the present cell at the beginning, and then use it in the computation of
  // the local contributions. This is what we do in the call to
  // <code>coefficient.value_list</code> in the fourth line of the loop.
  //
  // The second change is how we make use of this coefficient in computing the
  // cell matrix contributions. This is in the obvious way, and not worth more
  // comments. For the right hand side, we use a constant value again.
  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      cell_matrix = 0;
      cell_rhs = 0;

      fe_values.reinit (cell);

      coefficient.value_list (fe_values.get_quadrature_points(),
                              coefficient_values);

      for (unsigned int q_index=0; q_index<n_q_points; ++q_index)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += (coefficient_values[q_index] *
                                   fe_values.shape_grad(i,q_index) *
                                   fe_values.shape_grad(j,q_index) *
                                   fe_values.JxW(q_index));

            cell_rhs(i) += (fe_values.shape_value(i,q_index) *
                            1.0 *
                            fe_values.JxW(q_index));
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

  // With the matrix so built, we use zero boundary values again:
  std::map<types::global_dof_index,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
                                            0,
                                            ZeroFunction<dim>(),
                                            boundary_values);
  MatrixTools::apply_boundary_values (boundary_values,
                                      system_matrix,
                                      solution,
                                      system_rhs);
}


// @sect4{Step5::solve}

// The solution process again looks mostly like in the previous
// examples. However, we will now use a preconditioned conjugate gradient
// algorithm. It is not very difficult to make this change. In fact, the only
// thing we have to alter is that we need an object which will act as a
// preconditioner. We will use SSOR (symmetric successive overrelaxation),
// with a relaxation factor of 1.2. For this purpose, the
// <code>SparseMatrix</code> class has a function which does one SSOR step,
// and we need to package the address of this function together with the
// matrix on which it should act (which is the matrix to be inverted) and the
// relaxation factor into one object. The <code>PreconditionSSOR</code> class
// does this for us. (<code>PreconditionSSOR</code> class takes a template
// argument denoting the matrix type it is supposed to work on. The default
// value is <code>SparseMatrix@<double@></code>, which is exactly what we need
// here, so we simply stick with the default and do not specify anything in
// the angle brackets.)
//
// Note that for the present case, SSOR doesn't really perform much better
// than most other preconditioners (though better than no preconditioning at
// all). A brief comparison of different preconditioners is presented in the
// Results section of the next tutorial program, step-6.
//
// With this, the rest of the function is trivial: instead of the
// <code>PreconditionIdentity</code> object we have created before, we now use
// the preconditioner we have declared, and the CG solver will do the rest for
// us:
template <int dim>
void Step5<dim>::solve ()
{
  SolverControl           solver_control (1000, 1e-12);
  SolverCG<>              solver (solver_control);

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  solver.solve (system_matrix, solution, system_rhs,
                preconditioner);

  std::cout << "   " << solver_control.last_step()
            << " CG iterations needed to obtain convergence."
            << std::endl;
}


// @sect4{Step5::output_results and setting output flags}

// Writing output to a file is mostly the same as for the previous example,
// but here we will show how to modify some output options and how to
// construct a different filename for each refinement cycle.
template <int dim>
void Step5<dim>::output_results (const unsigned int cycle) const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");

  data_out.build_patches ();

  // For this example, we would like to write the output directly to a file in
  // Encapsulated Postscript (EPS) format. The library supports this, but
  // things may be a bit more difficult sometimes, since EPS is a printing
  // format, unlike most other supported formats which serve as input for
  // graphical tools. Therefore, you can't scale or rotate the image after it
  // has been written to disk, and you have to decide about the viewpoint or
  // the scaling in advance.
  //
  // The defaults in the library are usually quite reasonable, and regarding
  // viewpoint and scaling they coincide with the defaults of
  // Gnuplot. However, since this is a tutorial, we will demonstrate how to
  // change them. For this, we first have to generate an object describing the
  // flags for EPS output (similar flag classes exist for all supported output
  // formats):
  DataOutBase::EpsFlags eps_flags;
  // They are initialized with the default values, so we only have to change
  // those that we don't like. For example, we would like to scale the z-axis
  // differently (stretch each data point in z-direction by a factor of four):
  eps_flags.z_scaling = 4;
  // Then we would also like to alter the viewpoint from which we look at the
  // solution surface. The default is at an angle of 60 degrees down from the
  // vertical axis, and 30 degrees rotated against it in mathematical positive
  // sense. We raise our viewpoint a bit and look more along the y-axis:
  eps_flags.azimut_angle = 40;
  eps_flags.turn_angle   = 10;
  // That shall suffice. There are more flags, for example whether to draw the
  // mesh lines, which data vectors to use for colorization of the interior of
  // the cells, and so on. You may want to take a look at the documentation of
  // the EpsFlags structure to get an overview of what is possible.
  //
  // The only thing still to be done, is to tell the output object to use
  // these flags:
  data_out.set_flags (eps_flags);
  // The above way to modify flags requires recompilation each time we would
  // like to use different flags. This is inconvenient, and we will see more
  // advanced ways in step-19 where the output flags are determined at run
  // time using an input file (step-19 doesn't show many other things; you
  // should feel free to read over it even if you haven't done step-6 to
  // step-18 yet).

  // Finally, we need the filename to which the results are to be written. We
  // would like to have it of the form <code>solution-N.eps</code>, where N is
  // the number of the refinement cycle. Thus, we have to convert an integer
  // to a part of a string; this can be done using the <code>sprintf</code>
  // function, but in C++ there is a more elegant way: write everything into a
  // special stream (just like writing into a file or to the screen) and
  // retrieve what you wrote as a string. This applies the usual conversions
  // from integer to strings, and one could as well use stream modifiers such
  // as <code>setw</code>, <code>setprecision</code>, and so on. In C++, you
  // can do this by using the so-called stringstream classes:
  std::ostringstream filename;

  // In order to now actually generate a filename, we fill the stringstream
  // variable with the base of the filename, then the number part, and finally
  // the suffix indicating the file type:
  filename << "solution-"
           << cycle
           << ".eps";

  // We can get whatever we wrote to the stream using the <code>str()</code>
  // function. The result is a string which we have to convert to a char*
  // using the <code>c_str()</code> function. Use that as filename for the
  // output stream and then write the data to the file:
  std::ofstream output (filename.str().c_str());

  data_out.write_eps (output);
}



// @sect4{Step5::run}

// The second to last thing in this program is the definition of the
// <code>run()</code> function. In contrast to the previous programs, we will
// compute on a sequence of meshes that after each iteration is globally
// refined. The function therefore consists of a loop over 6 cycles. In each
// cycle, we first print the cycle number, and then have to decide what to do
// with the mesh. If this is not the first cycle, we simply refine the
// existing mesh once globally. Before running through these cycles, however,
// we have to generate a mesh:

// In previous examples, we have already used some of the functions from the
// <code>GridGenerator</code> class. Here we would like to read a grid from a
// file where the cells are stored and which may originate from someone else,
// or may be the product of a mesh generator tool.
//
// In order to read a grid from a file, we generate an object of data type
// GridIn and associate the triangulation to it (i.e. we tell it to fill our
// triangulation object when we ask it to read the file). Then we open the
// respective file and initialize the triangulation with the data in the file:
template <int dim>
void Step5<dim>::run ()
{
  GridIn<dim> grid_in;
  grid_in.attach_triangulation (triangulation);
  std::ifstream input_file("circle-grid.inp");
  // We would now like to read the file. However, the input file is only for a
  // two-dimensional triangulation, while this function is a template for
  // arbitrary dimension. Since this is only a demonstration program, we will
  // not use different input files for the different dimensions, but rather
  // kill the whole program if we are not in 2D:
  Assert (dim==2, ExcInternalError());
  // ExcInternalError is a globally defined exception, which may be thrown
  // whenever something is terribly wrong. Usually, one would like to use more
  // specific exceptions, and particular in this case one would of course try
  // to do something else if <code>dim</code> is not equal to two, e.g. create
  // a grid using library functions. Aborting a program is usually not a good
  // idea and assertions should really only be used for exceptional cases
  // which should not occur, but might due to stupidity of the programmer,
  // user, or someone else. The situation above is not a very clever use of
  // Assert, but again: this is a tutorial and it might be worth to show what
  // not to do, after all.

  // So if we got past the assertion, we know that dim==2, and we can now
  // actually read the grid. It is in UCD (unstructured cell data) format (but
  // the ending of the <code>UCD</code>-file is <code>inp</code>), as
  // supported as input format by the AVS Explorer (a visualization program),
  // for example:
  grid_in.read_ucd (input_file);
  // If you like to use another input format, you have to use an other
  // <code>grid_in.read_xxx</code> function. (See the documentation of the
  // <code>GridIn</code> class to find out what input formats are presently
  // supported.)

  // The grid in the file describes a circle. Therefore we have to use a
  // boundary object which tells the triangulation where to put new points on
  // the boundary when the grid is refined. This works in the same way as in
  // the first example. Note that the HyperBallBoundary constructor takes two
  // parameters, the center of the ball and the radius, but that their default
  // (the origin and 1.0) are the ones which we would like to use here.
  static const HyperBallBoundary<dim> boundary;
  triangulation.set_boundary (0, boundary);

  for (unsigned int cycle=0; cycle<6; ++cycle)
    {
      std::cout << "Cycle " << cycle << ':' << std::endl;

      if (cycle != 0)
        triangulation.refine_global (1);

      // Now that we have a mesh for sure, we write some output and do all the
      // things that we have already seen in the previous examples.
      std::cout << "   Number of active cells: "
                << triangulation.n_active_cells()
                << std::endl
                << "   Total number of cells: "
                << triangulation.n_cells()
                << std::endl;

      setup_system ();
      assemble_system ();
      solve ();
      output_results (cycle);
    }
}


// @sect3{The <code>main</code> function}

// The main function looks mostly like the one in the previous example, so we
// won't comment on it further:
int main ()
{
  deallog.depth_console (0);

  Step5<2> laplace_problem_2d;
  laplace_problem_2d.run ();

  // Finally, we have promised to trigger an exception in the
  // <code>Coefficient</code> class through the <code>Assert</code> macro we
  // have introduced there. For this, we have to call its
  // <code>value_list</code> function with two arrays of different size (the
  // number in parentheses behind the declaration of the object). We have
  // commented out these lines in order to allow the program to exit
  // gracefully in normal situations (we use the program in day-to-day testing
  // of changes to the library as well), so you will only get the exception by
  // un-commenting the following lines. Take a look at the Results section of
  // the program to see what happens when the code is actually run:
  /*
    Coefficient<2>    coefficient;
    std::vector<Point<2> > points (2);
    std::vector<double>    coefficient_values (1);
    coefficient.value_list (points, coefficient_values);
  */

  return 0;
}
