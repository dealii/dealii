/* $Id$ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 2000 */

				 // The first few files have already
				 // been covered in previous examples
				 // and will thus not be further
				 // commented on.
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
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_lib.lagrange.h>
#include <fe/fe_values.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>

				 // We will not read the grid from a
				 // file as in the previous example,
				 // but generate it using a function
				 // of the library. However, we will
				 // want to write out the locally
				 // refined grids in each step, so we
				 // need the following include file
				 // instead of ``grid_in.h'':
#include <grid/grid_out.h>

				 // When using locally refined grids,
				 // we will get so-called ``hanging
				 // nodes''. However, the standard
				 // finite element methods assumes
				 // that the discrete solution spaces
				 // be continuous, so we need to make
				 // sure that the degrees of freedom
				 // on hanging nodes conform to some
				 // constraints such that the global
				 // solution is continuous. The
				 // following file contains a class
				 // which is used to handle these
				 // constraints:
#include <dofs/dof_constraints.h>

				 // Finally, we would like to use a
				 // simple way to adaptively refine
				 // the grid. While in general,
				 // adaptivity is very
				 // problem-specific, the error
				 // indicator in the following file
				 // often yields quite nicely adapted
				 // grids for a wide class of
				 // problems.
#include <numerics/error_estimator.h>

#include <fstream>


				 // The main class is again almost
				 // unchanged. Two additions, however,
				 // are made: we have added the
				 // ``refine'' function, which is used
				 // to adaptively refine the grid
				 // (instead of the global refinement
				 // in the previous examples), and a
				 // variable which will hold the
				 // constraints associated to the
				 // hanging nodes.
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
    void refine_grid ();
    void output_results (const unsigned int cycle) const;

    Triangulation<dim>   triangulation;
    FEQ1<dim>            fe;
    DoFHandler<dim>      dof_handler;

				     // This is the new variable in
				     // the main class. We need an
				     // object which holds a list of
				     // the constraints from the
				     // hanging nodes:
    ConstraintMatrix     hanging_node_constraints;

    SparseMatrixStruct   sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       system_rhs;
};



template <int dim>
class Coefficient : public Function<dim> 
{
  public:
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
    
    virtual void value_list (const vector<Point<dim> > &points,
			     vector<double>            &values,
			     const unsigned int         component = 0) const;
};



template <int dim>
double Coefficient<dim>::value (const Point<dim> &p,
				const unsigned int) const 
{
  if (p.square() < 0.5*0.5)
    return 20;
  else
    return 1;
};



template <int dim>
void Coefficient<dim>::value_list (const vector<Point<dim> > &points,
				   vector<double>            &values,
				   const unsigned int component) const 
{
  const unsigned int n_points = points.size();

  Assert (values.size() == n_points, 
	  ExcVectorHasWrongSize (values.size(), n_points));
  
  Assert (component == 0, 
	  ExcWrongComponent (component, 1));
  
  for (unsigned int i=0; i<n_points; ++i)
    {
      if (points[i].square() < 0.5*0.5)
	values[i] = 20;
      else
	values[i] = 1;
    };
};


template <int dim>
LaplaceProblem<dim>::LaplaceProblem () :
		dof_handler (triangulation)
{};



template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);

				   // After setting up all the degrees
				   // of freedoms, we can make up the
				   // list of constraints associated
				   // with the hanging nodes. This is
				   // done using the following
				   // function calls (the first clears
				   // the contents of the object,
				   // which is still there from the
				   // previous cycle, i.e. before the
				   // grid was refined):
  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
					   hanging_node_constraints);
				   // In principle, the
				   // ConstraintMatrix class can hold
				   // other constraints as well,
				   // i.e. constraints that do not
				   // stem from hanging
				   // nodes. Sometimes, it is useful
				   // to use such constraints, in
				   // which case they may be added to
				   // the ConstraintMatrix object
				   // after the hanging node
				   // constraints were computed. After
				   // all constraints have been added,
				   // they need to be sorted and
				   // rearranged to perform some
				   // actions more efficiently. This
				   // postprocessing is done using the
				   // ``close'' function, after which
				   // no further constraints may be
				   // added any more.
  hanging_node_constraints.close ();

  sparsity_pattern.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);

				   // The constrained hanging nodes
				   // will later be eliminated from
				   // the linear system of
				   // equations. When doing so, some
				   // additional entries in the global
				   // matrix will be set to non-zero
				   // values, so we have to reserve
				   // some space for them here. Since
				   // the process of elimination of
				   // these constrained nodes is
				   // called ``condensation'', the
				   // functions that eliminate them
				   // are called ``condense'' for both
				   // the system matrix and right hand
				   // side, as well as for teh
				   // sparsity pattern.
  hanging_node_constraints.condense (sparsity_pattern);

				   // Now all non-zero entries of the
				   // matrix are known (i.e. those
				   // from regularly assembling the
				   // matrix and those that were
				   // introduced by eliminating
				   // constraints). We can thus close
				   // the sparsity pattern and remove
				   // unneeded space:
  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
};



template <int dim>
void LaplaceProblem<dim>::assemble_system () 
{  
  const Coefficient<dim> coefficient;

  QGauss3<dim>  quadrature_formula;

  FEValues<dim> fe_values (fe, quadrature_formula, 
			   UpdateFlags(update_values    |
				       update_gradients |
				       update_q_points  |
				       update_JxW_values));

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.n_quadrature_points;

  FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs (dofs_per_cell);

  vector<int>        local_dof_indices (dofs_per_cell);

  vector<double>     coefficient_values (n_q_points);

  DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
					endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      cell_matrix.clear ();
      cell_rhs.clear ();

      fe_values.reinit (cell);
      const FullMatrix<double> 
	& shape_values = fe_values.get_shape_values();
      const vector<vector<Tensor<1,dim> > >
	& shape_grads  = fe_values.get_shape_grads();
      const vector<double>
	& JxW_values   = fe_values.get_JxW_values();
      const vector<Point<dim> >
	& q_points     = fe_values.get_quadrature_points();

      coefficient.value_list (q_points, coefficient_values);
      
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  {
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      cell_matrix(i,j) += (coefficient_values[q_point] *
				   (shape_grads[i][q_point]    *
				    shape_grads[j][q_point])   *
				   JxW_values[q_point]);

	    cell_rhs(i) += (shape_values (i,q_point) *
			    1.0 *
			    fe_values.JxW (q_point));
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

  map<int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
					    0,
					    ZeroFunction<dim>(),
					    boundary_values);
  MatrixTools<dim>::apply_boundary_values (boundary_values,
					   system_matrix,
					   solution,
					   system_rhs);

				   // After the system of equations
				   // has been assembled just as for
				   // the previous examples, we still
				   // have to eliminate the
				   // constraints due to hanging
				   // nodes. This is done using the
				   // following two function calls:
  hanging_node_constraints.condense (system_matrix);
  hanging_node_constraints.condense (system_rhs);
				   // Using them, degrees of freedom
				   // associated to hanging nodes have
				   // been removed from the linear
				   // system and the independent
				   // variables are only regular
				   // nodes. The constrained nodes are
				   // still in the linear system
				   // (there is a one on the diagonal
				   // of the matrix and all other
				   // entries for this line are set to
				   // zero) but the computed values
				   // are invalid. They are set to
				   // reasonable values in the
				   // ``solve'' function.
};



template <int dim>
void LaplaceProblem<dim>::solve () 
{
  SolverControl           solver_control (1000, 1e-12);
  PrimitiveVectorMemory<> vector_memory;
  SolverCG<>              cg (solver_control, vector_memory);

  PreconditionRelaxation<>
    preconditioner(system_matrix,
		   &SparseMatrix<double>::template precondition_SSOR<double>,
		   1.2);

  cg.solve (system_matrix, solution, system_rhs,
	    preconditioner);

				   // To set the constrained nodes to
				   // resonable values, you have to
				   // use the following function. It
				   // computes the values of these
				   // nodes from the values of the
				   // unconstrained nodes, which are
				   // the solutions of the linear
				   // system just solved.
  hanging_node_constraints.distribute (solution);
};


				 // Instead of global refinement, we
				 // now use a slightly more elaborate
				 // scheme. We will use the
				 // ``KellyErrorEstimator'' class
				 // which implements an error
				 // estimator for the Laplace
				 // equation; it can in principle
				 // handle variable coefficients, but
				 // we will not use these advanced
				 // features, but rather use its most
				 // simple form since we are not
				 // interested in quantitative results
				 // but only in a quick way to
				 // generate locally refined grids.
				 //
				 // Although the error estimator
				 // derived by Kelly et al. was
				 // originally developed for Laplace's
				 // equation, we have found that it is
				 // also well suited to quickly
				 // generate locally refined grids for
				 // a wide class of
				 // problems. Basically, it looks at
				 // the jumps of the gradients of the
				 // solution over the faces of cells
				 // (which is a measure for the second
				 // derivatives) and scales it by the
				 // size of the cell. It is therefore
				 // a measure for the local smoothness
				 // of the solution at the place of
				 // each cell and it is thus
				 // understandable that it yields
				 // reasonable grids also for
				 // hyperbolic transport problems or
				 // the wave equation as well,
				 // although these grids are certainly
				 // suboptimal compared to approaches
				 // specially tailored to the
				 // problem. This error estimator may
				 // therefore be understood as a quick
				 // way to test an adaptive program.
template <int dim>
void LaplaceProblem<dim>::refine_grid ()
{
				   // The output of the error
				   // estimator class is an error
				   // indicator for each cell. We
				   // therefore need a vector with as
				   // many elements as there are
				   // active cells. Since accuracy is
				   // not that important here, the
				   // data type for the error values
				   // on each cell is ``float''
				   // instead of ``double''.
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

				   // Next, the error estimator can
				   // handle Neumann boundary
				   // conditions. For this, it needs
				   // to know which parts of the
				   // boundary have Neumann boundary
				   // conditions and teh respective
				   // boundary values there. This
				   // information is mediated by a map
				   // in which the keys are the
				   // boundary part numbers and the
				   // values are pointers to the
				   // boundary value functions. We
				   // create such a map, but since we
				   // do not use Neumann boundary
				   // conditions, the map will not
				   // contain entries.
  KellyErrorEstimator<dim>::FunctionMap neumann_boundary;

				   // Now we call the error
				   // estimator. The parameters should
				   // be clear apart from the
				   // quadrature formula: as said
				   // above, the jump of the gradients
				   // of the solution across the faces
				   // of a cell are considered. They
				   // are integrated along the face,
				   // but as usual in finite element
				   // programs the integration is done
				   // using quadrature. Since the
				   // error estimator class can't know
				   // itself which quadrature formula
				   // might be appropriate, we have to
				   // pass one to the function. Note
				   // that since the quadrature has to
				   // take place along faces, the
				   // dimension of the quadrature
				   // formula is ``dim-1'' rather then
				   // ``dim''.
  KellyErrorEstimator<dim>::estimate (dof_handler,
				      QGauss3<dim-1>(),
				      neumann_boundary,
				      solution,
				      estimated_error_per_cell);

				   // The above function returned one
				   // error indicator value for each
				   // cell in the
				   // ``estimated_error_per_cell''
				   // array. Refinement is now done as
				   // follows: refine those 30 per
				   // cent of the cells with the
				   // highest error values, and
				   // coarsen the 3 per cent of cells
				   // with the lowest values.
				   //
				   // One can easily verify that if
				   // the second number were zero,
				   // this would approximately result
				   // in a doubling of cells in each
				   // step in two space dimensions,
				   // since for each of the 30 per
				   // cent of cells four new would be
				   // replaced. In practice, some more
				   // cells are usually produced since
				   // it is disallowed that a cell is
				   // refined twice while the neighbor
				   // cell is not refined; in that
				   // case, the neighbor cell would be
				   // refined as well.
				   //
				   // In many applications, the number
				   // of cells to be coarsened would
				   // be set to something larger than
				   // only three per cent. A non-zero
				   // value is useful especially if
				   // for some reason the initial
				   // (coarse) grid is already rather
				   // refined. In that case, it might
				   // be necessary to refine it in
				   // some regions, while coarsening
				   // in some other regions is
				   // useful. In our case here, the
				   // initial grid is very coarse, so
				   // coarsening is only necessary in
				   // a few regions where
				   // over-refinement may have taken
				   // place. Thus a small, non-zero
				   // value is appropriate here.
  triangulation.refine_and_coarsen_fixed_number (estimated_error_per_cell,
						 0.3, 0.03);

				   // After the previous function has
				   // exited, some cells are flagged
				   // for refinement, and some other
				   // for coarsening. The refinement
				   // or coarsening itself is not
				   // performed by now, however, since
				   // there are many cases where
				   // further modifications of these
				   // flags is useful. Here, we don't
				   // want to do any such thing, so we
				   // can tell the triangulation to
				   // perform the actions for which
				   // the cells are flagged.
  triangulation.execute_coarsening_and_refinement ();
};



template <int dim>
void LaplaceProblem<dim>::output_results (const unsigned int cycle) const
{
				   // We want to write the grid in
				   // each cycle. Here is another way
				   // to quickly produce a filename
				   // based on the cycle number. It
				   // assumes that the numbers `0'
				   // through `9' are represented
				   // consecutively in the character
				   // set (which is the case in all
				   // known character sets). However,
				   // this will only work if the cycle
				   // number is less than ten, which
				   // we check by an assertion.
  string filename = "grid-";
  filename += ('0' + cycle);
  Assert (cycle < 10, ExcInternalError());
  
  filename += ".eps";
  ofstream output (filename.c_str());

				   // Using this filename, we write
				   // each grid as a postscript file.
  GridOut grid_out;
  grid_out.write_eps (triangulation, output);
};


template <int dim>
void LaplaceProblem<dim>::run () 
{
  for (unsigned int cycle=0; cycle<8; ++cycle)
    {
      cout << "Cycle " << cycle << ':' << endl;

      if (cycle == 0)
	{
					   // Instead of reading the
					   // grid from a file on disk
					   // as in the previous
					   // example, we now again
					   // create it using a
					   // library function. The
					   // domain is again a
					   // circle, which is why we
					   // have to provide a
					   // suitable boundary object
					   // as well.
					   //
					   // You will notice by
					   // looking at the coarse
					   // grid that it is of
					   // inferior quality than
					   // the one which we read
					   // from the file in the
					   // previous example: the
					   // cells are less equally
					   // formed. However, using
					   // the library function
					   // this program works in
					   // any space dimension,
					   // which was not the case
					   // before.
	  GridGenerator::hyper_ball (triangulation);

	  static const HyperBallBoundary<dim> boundary;
	  triangulation.set_boundary (0, boundary);

	  triangulation.refine_global (1);
	}
      else
					 // In case this is not the
					 // first cycle, we want to
					 // refine the grid. Unlike
					 // the global refinement
					 // employed in the last
					 // example, we now use the
					 // adaptive procedure
					 // described in the function
					 // which we now call:
	{
	  refine_grid ();
	};
      

      cout << "   Number of active cells: "
	   << triangulation.n_active_cells()
	   << endl;

      setup_system ();
      assemble_system ();
      solve ();
      output_results (cycle);
    };

				   // The solution on the final grid
				   // is now written to a file. As
				   // already done in one of the
				   // previous examples, we use the
				   // EPS format for output, and to
				   // obtain a reasonable view on the
				   // solution, we rescale the z-axis
				   // by a factor of four.
  DataOut<dim>::EpsFlags eps_flags;
  eps_flags.z_scaling = 4;
  
  DataOut<dim> data_out;
  data_out.set_flags (eps_flags);

  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");
  data_out.build_patches ();
  
  ofstream output ("final-solution.eps");
  data_out.write_eps (output);
};

    
				 // The main function is unaltered in
				 // its functionality against the
				 // previous example, but we have
				 // taken a step of additional
				 // caution. Sometimes, something goes
				 // wrong (such as insufficient disk
				 // space upon writing an output file,
				 // not enough memory when trying to
				 // allocate a vector or a matrix, or
				 // if we can't read from or write to
				 // a file for whatever reason), and
				 // in these cases the library will
				 // throw exceptions. Since they do
				 // not constitute programming errors,
				 // these exceptions also are not
				 // switched off in optimized mode, in
				 // constrast to the ``Assert'' macro
				 // which we have used to test against
				 // programming errors. If uncought,
				 // these exceptions propagate the
				 // call tree up to the ``main''
				 // function, and if they are not
				 // caught there either, the program
				 // is aborted. In many cases, like if
				 // there is not enough memory or disk
				 // space, we can't do anything but we
				 // can at least print some text
				 // trying to explain the reason why
				 // the program failed. A way to do so
				 // is shown in the following. It is
				 // certainly useful to write any
				 // larger program in this way, and
				 // you can do so by more or less
				 // copying this function apart from
				 // the ``try'' block which contains
				 // the code that constitutes the
				 // actual functionality.
int main () 
{

				   // The general idea behind the
				   // layout of this function is as
				   // follows: let's try to run the
				   // program as we did before...
  try
    {
      deallog.depth_console (0);

      LaplaceProblem<2> laplace_problem_2d;
      laplace_problem_2d.run ();
    }
				   // ...and if this should fail, try
				   // to gather as much information as
				   // possible. Specifically, if the
				   // exception that was thrown is an
				   // object of a class that is
				   // derived from the C++ standard
				   // class ``exception'', then we can
				   // use the ``what'' member function
				   // to get a string which describes
				   // the reason why the exception was
				   // thrown. 
				   //
				   // The deal.II exception classes
				   // are all derived from the
				   // standard class, and in
				   // particular, the ``exc.what()''
				   // function will return
				   // approximately the same string as
				   // would be generated if the
				   // exception was thrown using the
				   // ``Assert'' macro. You have seen
				   // the output of such an exception
				   // in the previous example, and you
				   // then know that it contains the
				   // file and line number of where
				   // the exception occured, and some
				   // other information. This is also
				   // what would be printed in the
				   // following.
  catch (exception &exc)
    {
      cerr << endl << endl
	   << "----------------------------------------------------"
	   << endl;
      cerr << "Exception on processing: " << exc.what() << endl
	   << "Aborting!" << endl
	   << "----------------------------------------------------"
	   << endl;
				       // We can't do much more than
				       // printing as much information
				       // as we can get to, so abort
				       // with error:
      return 1;
    }
				   // If the exception that was thrown
				   // somewhere was not an object of a
				   // class derived from the standard
				   // ``exception'' class, then we
				   // can't do anything at all. We
				   // then simply print an error
				   // message and exit.
  catch (...) 
    {
      cerr << endl << endl
	   << "----------------------------------------------------"
	   << endl;
      cerr << "Unknown exception!" << endl
	   << "Aborting!" << endl
	   << "----------------------------------------------------"
	   << endl;
      return 1;
    };

				   // If we got to this point, there
				   // was no exception which
				   // propagated up to the main
				   // functino (maybe there were some,
				   // but they were caught somewhere
				   // in the program or the
				   // library). Therefore, the program
				   // performed as was expected and we
				   // can return without error.
  return 0;
};
