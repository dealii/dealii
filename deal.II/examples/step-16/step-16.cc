/* $Id$ */
/* Author: Guido Kanschat, University of Heidelberg, 2003 */

/*    $Id$       */
/*    Version: $Name$                                          */
/*                                                                */
/*    Copyright (C) 2003 by the deal.II authors */
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
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
				 // These are the new include files
				 // required for multi-level methods.
				 // First, the file defining the
				 // multigrid method itself.
#include <multigrid/multigrid.h>
				 // The DoFHandler is replaced by an
				 // MGDoFHandler which is defined
				 // here.
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>

				 // Then, we need some pre-made
				 // transfer routines between grids.
#include <multigrid/mg_transfer.h>

				 // This is C++ ... see step 5 for
				 // further comments.
#include <fstream>
#ifdef HAVE_STD_STRINGSTREAM
#  include <sstream>
#else
#  include <strstream>
#endif


				 // This class is based on the same
				 // class in step 5. Remark that we
				 // replaced the DoFHandler by
				 // MGDoFHandler. since this inherits
				 // from DoFHandler, the new object
				 // incorporates the old functionality
				 // plus the new functions for degrees
				 // of freedom on different
				 // levels. Furthermore, we added
				 // MGLevelObjects for sparsity
				 // patterns and matrices.
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
    MGDoFHandler<dim>      mg_dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    MGLevelObject<SparsityPattern> mg_sparsity;
    MGLevelObject<SparseMatrix<float> > mg_matrices;
    
    Vector<double>       solution;
    Vector<double>       system_rhs;
};


				 // This function is as before.
template <int dim>
LaplaceProblem<dim>::LaplaceProblem () :
                fe (1),
		mg_dof_handler (triangulation)
{}



				 // This is the function of step 5
				 // augmented by the setup of the
				 // multi-level matrices in the end.
template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
  mg_dof_handler.distribute_dofs (fe);

  std::cout << "   Number of degrees of freedom: "
	    << mg_dof_handler.n_dofs()
	    << std::endl;

  sparsity_pattern.reinit (mg_dof_handler.n_dofs(),
			   mg_dof_handler.n_dofs(),
			   mg_dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (mg_dof_handler, sparsity_pattern);
  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);

  solution.reinit (mg_dof_handler.n_dofs());
  system_rhs.reinit (mg_dof_handler.n_dofs());

				   // The multi-level objects are
				   // resized to hold matrices for
				   // every level. The coarse level is
				   // zero (this is mandatory right
				   // now but may change in a future
				   // revision). Remark, that the
				   // finest level is nlevels-1.
  const unsigned int nlevels = triangulation.n_levels();
  mg_sparsity.resize(0, nlevels-1);
  mg_matrices.resize(0, nlevels-1);
  
  for (unsigned int level=0;level<nlevels;++level)
    {
    }
}



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

  typename DoFHandler<dim>::active_cell_iterator cell = mg_dof_handler.begin_active(),
						 endc = mg_dof_handler.end();
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
				       // It should be noted that the
				       // creation of the
				       // coefficient_values object is
				       // done outside the loop over
				       // all cells to avoid memory
				       // allocation each time we
				       // visit a new cell.
      
				       // With all this, the loops
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
	      cell_matrix(i,j) += (fe_values.shape_grad(i,q_point)
				   * fe_values.shape_grad(j,q_point)
				   * fe_values.JxW(q_point));

					     // For the right hand
					     // side, a constant value
					     // is used again:
	    cell_rhs(i) += (fe_values.shape_value(i,q_point)
			    * 1.0 * fe_values.JxW(q_point));
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
  VectorTools::interpolate_boundary_values (mg_dof_handler,
					    0,
					    ZeroFunction<dim>(),
					    boundary_values);
  MatrixTools::apply_boundary_values (boundary_values,
				      system_matrix,
				      solution,
				      system_rhs);
}



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
}



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

  data_out.attach_dof_handler (mg_dof_handler);
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
}



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
					   // Generate grid here!
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
}

    

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
}
