/* $Id$ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 1999 */

/*    $Id$       */
/*    Version: $Name$                                          */
/*                                                                */
/*    Copyright (C) 1999, 2000, 2001, 2002 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

				 // Again, the first few include files
				 // are already known, so we won't
				 // comment on them:
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/vector_memory.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>

				 // This one is new. We want to read a
				 // triangulation from disk, and the
				 // class which does this is declared
				 // in the following file:
#include <grid/grid_in.h>

				 // We will use a circular domain, and
				 // the object describing the boundary
				 // of it comes from this file:
#include <grid/tria_boundary_lib.h>

				 // This is C++ ...
#include <fstream>
				 // ... and this is too: We will
				 // convert integers to strings using
				 // the C++ stringstream class
				 // ``ostringstream''. One annoying
				 // complication arises here in that
				 // the classes ``std::istringstream''
				 // and ``std::ostringstream'' (with
				 // these names) have not been part of
				 // standard libraries of C++
				 // compilers for long. They have only
				 // been part of C++ compilers since
				 // around the time the C++ standard
				 // was made in 1999. For example, the
				 // gcc compiler up to and including
				 // version 2.95.2 did not have them,
				 // but instead provided classes
				 // ``istrstream'' and ``ostrstream''
				 // with a similar, but nevertheless
				 // slightly different
				 // interface. Furthermore, they were
				 // declared in the include file
				 // ``<strstream>'', while the new
				 // standards conforming classes are
				 // declared in ``<sstream>''. Many
				 // other compilers followed the gcc
				 // scheme, so whenever we want to
				 // support versions of compilers that
				 // appeared before approximately
				 // 2000/2001, we have to support
				 // these old classes.
				 //
				 // Since we do want to support these
				 // compilers, the ``./configure''
				 // script you run as the very first
				 // step of installing the library
				 // determines whether the compiler
				 // you want to use supports the new
				 // classes, or whether we have to
				 // fall back on the old ones. If the
				 // new classes are supported, then
				 // the preprocessor variable
				 // ``HAVE_STD_STRINGSTREAM'' is set
				 // in the ``base/config.h'' include
				 // file, that all include files in
				 // the library also include. Since we
				 // have included quite a number of
				 // files from the library at this
				 // point, the definition or
				 // non-definition of this
				 // preprocessor variable can now be
				 // used to decide whether old or new
				 // header names have to be used to
				 // import string stream classes:
#ifdef HAVE_STD_STRINGSTREAM
#  include <sstream>
#else
#  include <strstream>
#endif


				 // The main class is mostly as in the
				 // previous example. The most visible
				 // change is that the function
				 // ``make_grid_and_dofs'' has been
				 // removed, since making of the grid
				 // is now done in the ``run''
				 // function and the rest of its
				 // functionality now is in
				 // ``setup_system''. Apart from this,
				 // everything is as before.
template <int dim>
class LaplaceProblem 
{
  public:
    LaplaceProblem ();
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



				 // In this example, we want to use a
				 // variable coefficient in the
				 // elliptic operator. Of course, the
				 // suitable object is a Function, as
				 // we have used it for the right hand
				 // side and boundary values in the
				 // last example. We will use it
				 // again, but we implement another
				 // function ``value_list'' which
				 // takes a list of points and returns
				 // the values of the function at
				 // these points as a list. The reason
				 // why such a function is reasonable
				 // although we can get all the
				 // information from the ``value''
				 // function as well will be explained
				 // below when assembling the matrix.
				 //
				 // The need to declare a seemingly
				 // useless default constructor exists
				 // here just as in the previous
				 // example.
template <int dim>
class Coefficient : public Function<dim> 
{
  public:
    Coefficient ()  : Function<dim>() {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
    
    virtual void value_list (const std::vector<Point<dim> > &points,
			     std::vector<double>            &values,
			     const unsigned int              component = 0) const;
};



				 // This is the implementation of the
				 // coefficient function for a single
				 // point. We let it return 20 if the
				 // distance to the point of origin is
				 // less than 0.5, and 1 otherwise:
template <int dim>
double Coefficient<dim>::value (const Point<dim> &p,
				const unsigned int) const 
{
  if (p.square() < 0.5*0.5)
    return 20;
  else
    return 1;
};



				 // And this is the function that
				 // returns the value of the
				 // coefficient at a whole list of
				 // points at once. Of course, the
				 // values are the same as if we would
				 // ask the ``value'' function.
template <int dim>
void Coefficient<dim>::value_list (const std::vector<Point<dim> > &points,
				   std::vector<double>            &values,
				   const unsigned int              component) const 
{
				   // Use n_q_points as an
				   // abbreviation for the number of
				   // points for which function values
				   // are requested:
  const unsigned int n_points = points.size();

				   // Now, of course the size of the
				   // output array (``values'') must
				   // be the same as that of the input
				   // array (``points''), and we could
				   // simply assume that. However, in
				   // practice more than 90 per cent
				   // of programming errors are
				   // invalid function parameters such
				   // as invalid array sizes, etc, so
				   // we should try to make sure that
				   // the parameters are valid. For
				   // this, the Assert macro is a good
				   // means, since it asserts that the
				   // condition which is given as
				   // first argument is valid, and if
				   // not throws an exception (its
				   // second argument) which will
				   // usually terminate the program
				   // giving information where the
				   // error occured and what the
				   // reason was. This generally
				   // reduces the time to find
				   // programming errors dramatically
				   // and we have found assertions an
				   // invaluable means to program
				   // fast.
				   //
				   // On the other hand, all these
				   // checks (there are more than 2000
				   // of them in the library) should
				   // not slow down the program too
				   // much, which is why the Assert
				   // macro is only used in debug mode
				   // and expands to nothing if in
				   // optimized mode. Therefore, while
				   // you test your program and debug
				   // it, the assertions will tell you
				   // where the problems are, and once
				   // your program is stable you can
				   // switch off debugging and the
				   // program will run without the
				   // assertions and at maximum speed.
				   //
				   // Here, as has been said above, we
				   // would like to make sure that the
				   // size of the two arrays is equal,
				   // and if not throw an
				   // exception. Since the following
				   // test is rather frequent for the
				   // classes derived from
				   // ``Function'', that class
				   // declares an exception
				   // ``ExcDimensionMismatch'' which
				   // takes the sizes of two vectors
				   // and prints some output in case
				   // the condition is violated:
  Assert (values.size() == n_points, 
	  ExcDimensionMismatch (values.size(), n_points));
				   // Since examples are not very good
				   // if they do not demonstrate their
				   // point, we will show how to
				   // trigger this exception at the
				   // end of the main program, and
				   // what output results from this
				   // (see the ``Results'' section of
				   // this example program). You will
				   // certainly notice that the output
				   // is quite well suited to quickly
				   // find what the problem is and
				   // what parameters are expected. An
				   // additional plus is that if the
				   // program is run inside a
				   // debugger, it will stop at the
				   // point where the exception is
				   // triggered, so you can go up the
				   // call stack to immediately find
				   // the place where the the array
				   // with the wrong size was set up.
  
				   // While we're at it, we can do
				   // another check: the coefficient
				   // is a scalar, but the Function
				   // class also represents
				   // vector-valued function. A scalar
				   // function must therefore be
				   // considered as a vector-valued
				   // function with only one
				   // component, so the only valid
				   // component for which a user might
				   // ask is zero (we always count
				   // from zero). The following
				   // assertion checks this. (The
				   // ``1'' is denotes the number of
				   // components that this function
				   // has.)
  Assert (component == 0, 
	  ExcIndexRange (component, 0, 1));
  
  for (unsigned int i=0; i<n_points; ++i)
    {
      if (points[i].square() < 0.5*0.5)
	values[i] = 20;
      else
	values[i] = 1;
    };
};


				 // This function is as before.
template <int dim>
LaplaceProblem<dim>::LaplaceProblem () :
                fe (1),
		dof_handler (triangulation)
{};



				 // This is the function
				 // ``make_grid_and_dofs'' from the
				 // previous example, minus the
				 // generation of the grid. Everything
				 // else is unchanged.
template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);

  std::cout << "   Number of degrees of freedom: "
	    << dof_handler.n_dofs()
	    << std::endl;

  sparsity_pattern.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
};



				 // As in the previous examples, this
				 // function is not changed much with
				 // regard to its functionality, but
				 // there are still some optimizations
				 // which we will show. For this, it
				 // is important to note that if
				 // efficient solvers are used (such
				 // as the preconditions CG method),
				 // assembling the matrix and right
				 // hand side can take a comparable
				 // time, and you should think about
				 // using one or two optimizations at
				 // some places.
				 //
				 // What we will show here is how we
				 // can avoid calls to the
				 // shape_value, shape_grad, and
				 // quadrature_point functions of the
				 // FEValues object, and in particular
				 // optimize away most of the virtual
				 // function calls of the Function
				 // object. The way to do so will be
				 // explained in the following, while
				 // those parts of this function that
				 // are not changed with respect to
				 // the previous example are not
				 // commented on.
template <int dim>
void LaplaceProblem<dim>::assemble_system () 
{  
				   // This time, we will again use a
				   // constant right hand side
				   // function, but a variable
				   // coefficient. The following
				   // object will be used for this:
  const Coefficient<dim> coefficient;

  QGauss2<dim>  quadrature_formula;

  FEValues<dim> fe_values (fe, quadrature_formula, 
			   UpdateFlags(update_values    |
				       update_gradients |
				       update_q_points  |
				       update_JxW_values));

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.n_quadrature_points;

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

				   // Below, we will ask the
				   // Coefficient class to compute the
				   // values of the coefficient at all
				   // quadrature points on one cell at
				   // once. For this, we need some
				   // space to store the values in,
				   // which we use the following
				   // variable for:
  std::vector<double>     coefficient_values (n_q_points);

  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
						 endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      cell_matrix.clear ();
      cell_rhs.clear ();

				       // As before, we want the
				       // FEValues object to compute
				       // the quantities which we told
				       // him to compute in the
				       // constructor using the update
				       // flags.
      fe_values.reinit (cell);
				       // Now, these quantities are
				       // stored in arrays in the
				       // FEValues object. Usually,
				       // the internals of how and
				       // where they are stored is not
				       // something that the outside
				       // world should know, but since
				       // this is a time critical
				       // function we decided to
				       // publicize these arrays a
				       // little bit, and provide
				       // facilities to export the
				       // address where this data is
				       // stored.
				       //
				       // For example, the values of
				       // shape function j at
				       // quadrature point q is stored
				       // in a matrix, of which we can
				       // get the address as follows
				       // (note that this is a
				       // reference to the matrix,
				       // symbolized by the ampersand ``&'',
				       // and that it must be a
				       // constant reference, since
				       // only read-only access is
				       // granted):
      const FullMatrix<double> 
	& shape_values = fe_values.get_shape_values();
				       // Instead of writing
				       // fe_values.shape_value(j,q)
				       // we can now write
				       // shape_values[j][q], i.e. the
				       // function call needed
				       // previously for each access
				       // will be optimized away.
				       //
				       // There are alike functions
				       // for almost all data elements
				       // in the FEValues class. The
				       // gradient are accessed as
				       // follows:
      const std::vector<std::vector<Tensor<1,dim> > >
	& shape_grads  = fe_values.get_shape_grads();
				       // The data type looks a bit
				       // unwieldy, since each entry
				       // in the matrix (j,q) now
				       // needs to be the gradient of
				       // the shape function, which is
				       // a tensor rather than a
				       // scalar.
				       //
				       // Similarly, access to the
				       // place where quadrature
				       // points and the determinants
				       // of the Jacobian matrices
				       // times the weights of the
				       // respective quadrature points
				       // are stored, can be obtained
				       // like this:
      const std::vector<double>
	& JxW_values   = fe_values.get_JxW_values();
      const std::vector<Point<dim> >
	& q_points     = fe_values.get_quadrature_points();
				       // Admittedly, the declarations
				       // above are not easily
				       // readable, but they can save
				       // many function calls in the
				       // inner loops and can thus
				       // make assemblage faster.
				       //
				       // An additional advantage is
				       // that the inner loops are
				       // simpler to read, since the
				       // fe_values object is no more
				       // explicitely needed to access
				       // the different fields (see
				       // below).

				       // There is one more thing: in
				       // this example, we want to use
				       // a non-constant
				       // coefficient. In the previous
				       // example, we have called the
				       // ``value'' function of the
				       // right hand side object for
				       // each quadrature
				       // point. Unfortunately, that
				       // is a virtual function, so
				       // calling it is relatively
				       // expensive. Therefore, we use
				       // a function of the Function
				       // class which returns the
				       // values at all quadrature
				       // points at once; that
				       // function is still virtual,
				       // but it needs to be computed
				       // once per cell only, not once
				       // in the inner loop:
      coefficient.value_list (q_points, coefficient_values);
				       // It should be noted that the
				       // creation of the
				       // coefficient_values object is
				       // done outside the loop over
				       // all cells to avoid memory
				       // allocation each time we
				       // visit a new cell. Contrary
				       // to this, the other variables
				       // above were created inside
				       // the loop, but they were only
				       // references to memory that
				       // has already been allocated
				       // (i.e. they are pointers to
				       // that memory) and therefore,
				       // no new memory needs to be
				       // allocated; in particular, by
				       // declaring the pointers as
				       // close to their use as
				       // possible, we give the
				       // compiler a better choice to
				       // optimize them away
				       // altogether, something which
				       // it definitely can't do with
				       // the coefficient_values
				       // object since it is too
				       // complicated, but mostly
				       // because it's address is
				       // passed to a virtual function
				       // which is not knows at
				       // compile time.
      
				       // Using the various
				       // abbreviations, the loops
				       // then look like this (the
				       // parentheses around the
				       // product of the two gradients
				       // are needed to indicate the
				       // dot product; we have to
				       // overrule associativity of
				       // the operator* here, since
				       // the compiler would otherwise
				       // complain about an undefined
				       // product of double*gradient
				       // since it parses
				       // left-to-right):
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  {
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      cell_matrix(i,j) += (coefficient_values[q_point] *
				   (shape_grads[i][q_point]    *
				    shape_grads[j][q_point])   *
				   JxW_values[q_point]);

					     // For the right hand
					     // side, a constant value
					     // is used again:
	    cell_rhs(i) += (shape_values[i][q_point] *
			    1.0 *
			    JxW_values[q_point]);
	  };


      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    system_matrix.add (local_dof_indices[i],
			       local_dof_indices[j],
			       cell_matrix(i,j));
	  
	  system_rhs(local_dof_indices[i]) += cell_rhs(i);
	};
    };

				   // Again use zero boundary values:
  std::map<unsigned int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
					    0,
					    ZeroFunction<dim>(),
					    boundary_values);
  MatrixTools::apply_boundary_values (boundary_values,
				      system_matrix,
				      solution,
				      system_rhs);
};



				 // The solution process again looks
				 // mostly like in the previous
				 // examples. However, we will now use
				 // a preconditioned conjugate
				 // gradient algorithm. It is not very
				 // difficult to make this change:
template <int dim>
void LaplaceProblem<dim>::solve () 
{
  SolverControl           solver_control (1000, 1e-12);
  PrimitiveVectorMemory<> vector_memory;
  SolverCG<>              cg (solver_control, vector_memory);

				   // The only thing we have to alter
				   // is that we need an object which
				   // will act as a preconditioner. We
				   // will use SSOR (symmetric
				   // successive overrelaxation), with
				   // a relaxation factor of 1.2. For
				   // this purpose, the SparseMatrix
				   // class has a function which does
				   // one SSOR step, and we need to
				   // package the address of this
				   // function together with the
				   // matrix on which it should act
				   // (which is the matrix to be
				   // inverted) and the relaxation
				   // factor into one object. This can
				   // be done like this:
  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);
				   // (Note that we did not have to
				   // explicitely pass the address of
				   // the SSOR function of the matrix
				   // to this objects, rather it is
				   // hardcoded into the object, thus
				   // the name.)
				   //
				   // The default template parameters
				   // of the ``PreconditionRelaxation''
				   // class is the matrix type, which
				   // defaults to the types used in
				   // this program.

				   // Calling the solver now looks
				   // mostly like in the example
				   // before, but where there was an
				   // object of type
				   // PreconditionIdentity before,
				   // there now is the newly generated
				   // preconditioner object.
  cg.solve (system_matrix, solution, system_rhs,
	    preconditioner);

  std::cout << "   " << solver_control.last_step()
	    << " CG iterations needed to obtain convergence."
	    << std::endl;
};



				 // Writing output to a file is mostly
				 // the same as for the previous
				 // example, but here we will show how
				 // to modify some output options and
				 // how to construct a different
				 // filename for each refinement
				 // cycle.
template <int dim>
void LaplaceProblem<dim>::output_results (const unsigned int cycle) const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");

  data_out.build_patches ();

				   // For this example, we would like
				   // to write the output directly to
				   // a file in Encapsulated
				   // Postscript (EPS) format. The
				   // library supports this, but
				   // things may be a bit more
				   // difficult sometimes, since EPS
				   // is a printing format, unlike
				   // most other supported formats
				   // which serve as input for
				   // graphical tools. Therefore, you
				   // can't scale or rotate the image
				   // after it has been written to
				   // disk, and you have to decide
				   // about the viewpoint or the
				   // scaling in advance.
				   //
				   // The defaults in the library are
				   // usually quite reasonable, and
				   // regarding viewpoint and scaling
				   // they coincide with the defaults
				   // of Gnuplot. However, since this
				   // is a tutorial, we will
				   // demonstrate how to change
				   // them. For this, we first have to
				   // generate an object describing
				   // the flags for EPS output:
  DataOutBase::EpsFlags eps_flags;
				   // They are initialized with the
				   // default values, so we only have
				   // to change those that we don't
				   // like. For example, we would like
				   // to scale the z-axis differently
				   // (stretch each data point in
				   // z-direction by a factor of four):
  eps_flags.z_scaling = 4;
				   // Then we would also like to alter
				   // the viewpoint from which we look
				   // at the solution surface. The
				   // default is at an angle of 60
				   // degrees down from the vertical
				   // axis, and 30 degrees rotated
				   // against it in mathematical
				   // positive sense. We raise our
				   // viewpoint a bit and look more
				   // along the y-axis:
  eps_flags.azimut_angle = 40;
  eps_flags.turn_angle   = 10;
				   // That shall suffice. There are
				   // more flags, for example whether
				   // to draw the mesh lines, which
				   // data vectors to use for
				   // colorization of the interior of
				   // the cells, and so on. You may
				   // want to take a look at the
				   // documentation of the EpsFlags
				   // structure to get an overview of
				   // what is possible.
				   //
				   // The only thing still to be done,
				   // is to tell the output object to
				   // use these flags:
  data_out.set_flags (eps_flags);
				   // The above way to modify flags
				   // requires recompilation each time
				   // we would like to use different
				   // flags. This is inconvenient, and
				   // we will see more advanced ways
				   // in following examples where the
				   // output flags are determined at
				   // run time using an input file.

				   // Finally, we need the filename to
				   // which the results are to be
				   // written. We would like to have
				   // it of the form
				   // ``solution-N.eps'', where N is
				   // the number of the refinement
				   // cycle. Thus, we have to convert
				   // an integer to a part of a
				   // string; this can be done using
				   // the ``sprintf'' function, but in
				   // C++ there is a more elegant way:
				   // write everything into a special
				   // stream (just like writing into a
				   // file or to the screen) and
				   // retrieve what you wrote as a
				   // string. This applies the usual
				   // conversions from integer to
				   // strings, and one could as well
				   // give stream modifiers such as
				   // ``setw'', ``setprecision'', and
				   // so on.
				   //
				   // In C++, you can do this by using
				   // the so-called stringstream
				   // classes. As already discussed at
				   // the point of inclusion of the
				   // respective header file above,
				   // there is some historical
				   // confusion we have to work around
				   // here, since the class we'd like
				   // to use used to be called
				   // ``ostrstream'', but now is named
				   // ``ostringstream''. In the same
				   // way as done above in deciding
				   // which file to include, we here
				   // decide which class name to use:
#ifdef HAVE_STD_STRINGSTREAM
  std::ostringstream filename;
#else
  std::ostrstream filename;
#endif
				   // Fortunately, the interface of
				   // the two classes which we might
				   // now be using, depending on which
				   // one is available, is close
				   // enough that we need to take care
				   // about the differences only once
				   // below, so we can use them in a
				   // rather straightforward way, even
				   // if they are not identical.

				   // In order to now actually
				   // generate a filename, we fill the
				   // stringstream variable with the
				   // base of the filename, then the
				   // number part, and finally the
				   // suffix indicating the file type:
  filename << "solution-"
	   << cycle
	   << ".eps";
  
				   // For the old string stream
				   // classes, we have to append the
				   // final '\0' that appears at the
				   // end of ``char *''
				   // variables. This is done by the
				   // following construct:
#ifndef HAVE_STD_STRINGSTREAM
  filename << std::ends;
#endif
				   // We can get whatever we wrote to
				   // the stream using the ``str()''
				   // function. If the new
				   // stringstream classes are used,
				   // then the result is a string
				   // which we have to convert to a
				   // char* using the ``c_str()''
				   // function, otherwise the result
				   // is a char* right away. Use that
				   // as filename for the output
				   // stream:
#ifdef HAVE_STD_STRINGSTREAM
  std::ofstream output (filename.str().c_str());
#else
  std::ofstream output (filename.str());
#endif
				   // And then write the data to the
				   // file.
  data_out.write_eps (output);
};



template <int dim>
void LaplaceProblem<dim>::run () 
{
  for (unsigned int cycle=0; cycle<6; ++cycle)
    {
      std::cout << "Cycle " << cycle << ':' << std::endl;

				       // If this is the first round,
				       // then we have no grid yet,
				       // and we will create it
				       // here. In previous examples,
				       // we have already used some of
				       // the functions from the
				       // GridGenerator class. Here we
				       // would like to read a grid
				       // from a file where the cells
				       // are stored and which may
				       // originate from someone else,
				       // or may be the product of a
				       // mesh generator tool.
				       //
				       // In order to read a grid from
				       // a file, we generate an
				       // object of data type GridIn
				       // and associate the
				       // triangulation to it (i.e. we
				       // tell it to fill our
				       // triangulation object when we
				       // ask it to read the
				       // file). Then we open the
				       // respective file and
				       // initialize the triangulation
				       // with the data in the file:
      if (cycle == 0)
	{
	  GridIn<dim> grid_in;
	  grid_in.attach_triangulation (triangulation);
	  std::ifstream input_file("circle-grid.inp");
					   // We would now like to
					   // read the file. However,
					   // the input file is only
					   // for a two-dimensional
					   // triangulation, while
					   // this function is a
					   // template for arbitrary
					   // dimension. Since this is
					   // only a demonstration
					   // program, we will not use
					   // different input files
					   // for the different
					   // dimensions, but rather
					   // kill the whole program
					   // if we are not in 2D:
	  Assert (dim==2, ExcInternalError());
					   // ExcInternalError is a
					   // globally defined
					   // exception, which may be
					   // thrown whenever
					   // something is terribly
					   // wrong. Usually, one
					   // would like to use more
					   // specific exceptions, and
					   // particular in this case
					   // one would of course try
					   // to do something else if
					   // ``dim'' is not equal to
					   // two, e.g. create a grid
					   // using library
					   // functions. Aborting a
					   // program is usually not a
					   // good idea and assertions
					   // should really only be
					   // used for exceptional
					   // cases which should not
					   // occur, but might due to
					   // stupidity of the
					   // programmer, user, or
					   // someone else. The
					   // situation above is not a
					   // very clever use of
					   // Assert, but again: this
					   // is a tutorial and it
					   // might be worth to show
					   // what not to do, after
					   // all.
	  
					   // We can now actually read
					   // the grid. It is in UCD
					   // (unstructured cell data)
					   // format (but the ending
					   // of the ``UCD''-file is
					   // ``inp''), as supported
					   // as input format by the
					   // AVS Explorer (a
					   // visualization program),
					   // for example:
	  grid_in.read_ucd (input_file);
                                           // If you like to use
                                           // another input format,
                                           // you have to use an other
                                           // ``grid_in.read_xxx''
                                           // function. (See the
                                           // documentation of the
                                           // ``GridIn'' class to find
                                           // out what input formats
                                           // are presently
                                           // supported.)

					   // The grid in the file
					   // describes a
					   // circle. Therefore we
					   // have to use a boundary
					   // object which tells the
					   // triangulation where to
					   // put new points on the
					   // boundary when the grid
					   // is refined. This works
					   // in the same way as in
					   // the first example. Note
					   // that the
					   // HyperBallBoundary
					   // constructor takes two
					   // parameters, the center
					   // of the ball and the
					   // radius, but that their
					   // default (the origin and
					   // 1.0) are the ones which
					   // we would like to use
					   // here.
	  static const HyperBallBoundary<dim> boundary;
	  triangulation.set_boundary (0, boundary);
	}
				       // If this is not the first
				       // cycle, then simply refine
				       // the grid once globally.
      else
	triangulation.refine_global (1);

				       // Write some output and do all
				       // the things that we have
				       // already seen in the previous
				       // examples.
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
    };
};

    

				 // The main function looks mostly
				 // like the one in the previous
				 // example, so we won't comment on it
				 // further.
int main () 
{
  deallog.depth_console (0);

  LaplaceProblem<2> laplace_problem_2d;
  laplace_problem_2d.run ();

				   // Finally, we have promised to
				   // trigger an exception in the
				   // Coefficient class. For this, we
				   // have to call its ``value_list''
				   // function with two arrays of
				   // different size (the number in
				   // parentheses behind the name of
				   // the object). We have commented
				   // out these lines in order to
				   // allow the program to exit
				   // gracefully in normal situations
				   // (we use the program in
				   // day-to-day testing of changes to
				   // the library as well), so you
				   // will only get the exception by
				   // un-commenting the following
				   // lines.
/*  
  Coefficient<2>    coefficient;
  std::vector<Point<2> > points (2);
  std::vector<double>    coefficient_values (1);
  coefficient.value_list (points, coefficient_values);
*/
  
  return 0;
};
