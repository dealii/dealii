/* $Id$ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 1999 */

				 // The first few (many?) include
				 // files have already been used in
				 // the previous example, so we will
				 // not explain their meaning here
				 // again.
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_q.h>
#include <dofs/dof_tools.h>
#include <fe/fe_values.h>
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/vector_memory.h>
#include <lac/precondition.h>

#include <numerics/data_out.h>
#include <fstream>

				 // This is new, however: in the
				 // previous example we got some
				 // unwanted output from the linear
				 // solvers. If we want to suppress
				 // it, we have to include this file
				 // and add a line somewhere to the
				 // program; in this program, it was
				 // added to the main function.
#include <base/logstream.h>



				 // This is again the same
				 // LaplaceProblem class as in the
				 // previous example. The only
				 // difference is that we have now
				 // declared it as a class with a
				 // template parameter, and the
				 // template parameter is of course
				 // the spatial dimension in which we
				 // would like to solve the Laplace
				 // equation. Of course, several of
				 // the member variables depend on
				 // this dimension as well, in
				 // particular the Triangulation
				 // class, which has to represent
				 // quadrilaterals or hexahedra,
				 // respectively. Apart from this,
				 // everything is as before.
template <int dim>
class LaplaceProblem 
{
  public:
    LaplaceProblem ();
    void run ();
    
  private:
    void make_grid_and_dofs ();
    void assemble_system ();
    void solve ();
    void output_results () const;

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       system_rhs;
};


				 // In the following, we declare two
				 // more classes, which will represent
				 // the functions of the
				 // dim-dimensional space denoting the
				 // right hand side and the
				 // non-homogeneous Dirichlet boundary
				 // values.
				 //
				 // Each of these classes is derived
				 // from a common, abstract base class
				 // Function, which declares the
				 // common interface which all
				 // functions have to follow. In
				 // particular, concrete classes have
				 // to overload the `value' function,
				 // which takes a point in
				 // dim-dimensional space as
				 // parameters and shall return the
				 // value at that point as a `double'
				 // variable.
				 //
				 // The `value' function takes a
				 // second argument, which we have
				 // here named `component': This is
				 // only meant for vector valued
				 // functions, where you may want to
				 // access a certain component of the
				 // vector at the point `p'. However,
				 // our functions are scalar, so we
				 // need not worry about this
				 // parameter and we will not use it
				 // in the implementation of the
				 // functions. Note that in the base
				 // class (Function), the declaration
				 // of the `value' function has a
				 // default value of zero for the
				 // component, so we will access the
				 // `value' function of the right hand
				 // side with only one parameter,
				 // namely the point where we want to
				 // evaluate the function.
template <int dim>
class RightHandSide : public Function<dim> 
{
  public:
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};



template <int dim>
class BoundaryValues : public Function<dim> 
{
  public:
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};




				 // We wanted the right hand side
				 // function to be 4*(x**4+y**4) in
				 // 2D, or 4*(x**4+y**4+z**4) in
				 // 3D. Unfortunately, this is not as
				 // elegantly feasible dimension
				 // independently as much of the rest
				 // of this program, so we have to do
				 // it using a small
				 // loop. Fortunately, the compiler
				 // knows the size of the loop at
				 // compile time, i.e. the number of
				 // times the body will be executed,
				 // so it can optimize away the
				 // overhead needed for the loop and
				 // the result will be as fast as if
				 // we had used the formulas above
				 // right away.
				 //
				 // Note that the different
				 // coordinates (i.e. `x', `y', ...)
				 // of the point are accessed using
				 // the () operator.
template <int dim>
double RightHandSide<dim>::value (const Point<dim> &p,
				  const unsigned int) const 
{
  double return_value = 0;
  for (unsigned int i=0; i<dim; ++i)
    return_value += 4*std::pow(p(i), 4);

  return return_value;
};


				 // The boundary values were to be
				 // chosen to be x*x+y*y in 2D, and
				 // x*x+y*y+z*z in 3D. This happens to
				 // be equal to the square of the
				 // vector from the origin to the
				 // point at which we would like to
				 // evaluate the function,
				 // irrespective of the dimension. So
				 // that is what we return:
template <int dim>
double BoundaryValues<dim>::value (const Point<dim> &p,
				   const unsigned int) const 
{
  return p.square();
};




				 // This is the constructor of the
				 // LaplaceProblem class. It specifies
				 // the desired polynomial degree of
				 // the finite elements and associates
				 // the DoFHandler to the
				 // triangulation just as in the
				 // previous example.
template <int dim>
LaplaceProblem<dim>::LaplaceProblem () :
                fe (1),
		dof_handler (triangulation)
{};



				 // Grid creation is something
				 // inherently dimension
				 // dependent. However, as long as the
				 // domains are sufficiently similar
				 // in 2D or 3D, the library can
				 // abstract for you. In our case, we
				 // would like to again solve on the
				 // square [-1,1]x[-1,1] in 2D, or on
				 // the cube [-1,1]x[-1,1]x[-1,1] in
				 // 3D; both can be termed
				 // ``hyper_cube'', so we may use the
				 // same function in whatever
				 // dimension we are. Of course, the
				 // functions that create a hypercube
				 // in two and three dimensions are
				 // very much different, but that is
				 // something you need not care
				 // about. Let the library handle the
				 // difficult things.
				 //
				 // Likewise, associating a degree of
				 // freedom with each vertex is
				 // something which certainly looks
				 // different in 2D and 3D, but that
				 // does not need to bother you. This
				 // function therefore looks exactly
				 // like in the previous example,
				 // although it performs actions that
				 // in their details are quite
				 // different. The only significant
				 // difference is the number of cells
				 // resulting, which is much higher in
				 // three than in two space
				 // dimensions!
template <int dim>
void LaplaceProblem<dim>::make_grid_and_dofs ()
{
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (4);
  
  std::cout << "   Number of active cells: "
	    << triangulation.n_active_cells()
	    << std::endl
	    << "   Total number of cells: "
	    << triangulation.n_cells()
	    << std::endl;

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



				 // Unlike in the previous example, we
				 // would now like to use a
				 // non-constant right hand side
				 // function and non-zero boundary
				 // values. Both are tasks that are
				 // readily achieved with a only a few
				 // new lines of code in the
				 // assemblage of the matrix and right
				 // hand side.
				 //
				 // More interesting, though, is the
				 // way we assemble matrix and right
				 // hand side vector dimension
				 // independently: there is simply no
				 // difference to the pure
				 // two-dimensional case. Since the
				 // important objects used in this
				 // function (quadrature formula,
				 // FEValues) depend on the dimension
				 // by way of a template parameter as
				 // well, they can take care of
				 // setting up properly everything for
				 // the dimension for which this
				 // function is compiled. By declaring
				 // all classes which might depend on
				 // the dimension using a template
				 // parameter, the library can make
				 // nearly all work for you and you
				 // don't have to care about most
				 // things.
template <int dim>
void LaplaceProblem<dim>::assemble_system () 
{  
  QGauss2<dim>  quadrature_formula;

				   // We wanted to have a non-constant
				   // right hand side, so we use an
				   // object of the class declared
				   // above to generate the necessary
				   // data. Since this right hand side
				   // object is only used in this
				   // function, we only declare it
				   // here, rather than as a member
				   // variable of the LaplaceProblem
				   // class, or somewhere else.
  const RightHandSide<dim> right_hand_side;

				   // Compared to the previous
				   // example, in order to evaluate
				   // the non-constant right hand side
				   // function we now also need the
				   // quadrature points on the cell we
				   // are presently on (previously,
				   // they were only needed on the
				   // unit cell, in order to compute
				   // the values and gradients of the
				   // shape function, which are
				   // defined on the unit cell
				   // however). We can tell the
				   // FEValues object to do for us by
				   // giving it the update_q_points
				   // flag:
  FEValues<dim> fe_values (fe, quadrature_formula, 
			   UpdateFlags(update_values    |
				       update_gradients |
				       update_q_points  |
				       update_JxW_values));

				   // Note that the following numbers
				   // depend on the dimension which we
				   // are presently using. However,
				   // the FE and Quadrature classes do
				   // all the necessary work for you
				   // and you don't have to care about
				   // the dimension dependent parts:
  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.n_quadrature_points;

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

				   // Note here, that a cell is a
				   // quadrilateral in two space
				   // dimensions, but a hexahedron in
				   // 3D. In fact, the
				   // active_cell_iterator data type
				   // is something different,
				   // depending on the dimension we
				   // are in, but to the outside world
				   // they look alike and you will
				   // probably never see a difference
				   // although they are totally
				   // unrelated.
  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
						 endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      cell_matrix.clear ();
      cell_rhs.clear ();

				       // Now we have to assemble the
				       // local matrix and right hand
				       // side. This is done exactly
				       // like in the previous
				       // example, but now we revert
				       // the order of the loops
				       // (which we can safely do
				       // since they are independent
				       // of each other) and merge the
				       // loops for the local matrix
				       // and the local vector as far
				       // as possible; this makes
				       // things a bit faster.
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  {
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      cell_matrix(i,j) += (fe_values.shape_grad (i, q_point) *
				   fe_values.shape_grad (j, q_point) *
				   fe_values.JxW (q_point));

					     // Here is about the only
					     // difference to the
					     // previous example:
					     // instead of using a
					     // constant right hand
					     // side, we use the
					     // respective object and
					     // evaluate it at the
					     // quadrature points.
	    cell_rhs(i) += (fe_values.shape_value (i, q_point) *
			    right_hand_side.value (fe_values.quadrature_point (q_point)) *
			    fe_values.JxW (q_point));
	  };
      
				       // The transfer into the global
				       // matrix and right hand side
				       // is done exactly as before,
				       // but here we have again
				       // merged some loops for
				       // efficiency:
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

  
				   // We wanted to have
				   // non-homogeneous boundary values
				   // in this example, contrary to the
				   // one before. This is a simple
				   // task, we only have to replace
				   // the ZeroFunction used there by
				   // an object of the class which
				   // describes the boundary values we
				   // would like to use (i.e. the
				   // BoundaryValues class declared
				   // above):
  std::map<unsigned int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
					    0,
					    BoundaryValues<dim>(),
					    boundary_values);
  MatrixTools<dim>::apply_boundary_values (boundary_values,
					   system_matrix,
					   solution,
					   system_rhs);
};


				 // Solving the linear system of
				 // equation is something that looks
				 // almost identical in most
				 // programs. In particular, it is
				 // dimension independent, so this
				 // function is mostly copied from the
				 // previous example.
template <int dim>
void LaplaceProblem<dim>::solve () 
{
  SolverControl           solver_control (1000, 1e-12);
  PrimitiveVectorMemory<> vector_memory;
  SolverCG<>              cg (solver_control, vector_memory);
  cg.solve (system_matrix, solution, system_rhs,
	    PreconditionIdentity());

				   // We have made one addition,
				   // though: since we suppress output
				   // from the linear solvers, we have
				   // to print the number of
				   // iterations by hand.
  std::cout << "   " << solver_control.last_step()
	    << " CG iterations needed to obtain convergence."
	    << std::endl;
};



				 // This function also does what the
				 // respective one did in the previous
				 // example. No changes here for
				 // dimension independence either.
template <int dim>
void LaplaceProblem<dim>::output_results () const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");

  data_out.build_patches ();

				   // Only difference to the previous
				   // example: write output in GMV
				   // format, rather than for
				   // gnuplot. We use the dimension in
				   // the filename to generate
				   // distinct filenames for each run
				   // (in a better program, one would
				   // check whether `dim' can have
				   // other values than 2 or 3, but we
				   // neglect this here for the sake
				   // of brevity).
  std::ofstream output (dim == 2 ?
			"solution-2d.gmv" :
			"solution-3d.gmv");
  data_out.write_gmv (output);
};



				 // This is the function which has the
				 // top-level control over
				 // everything. Apart from one line of
				 // additional output, it is the same
				 // as for the previous example.
template <int dim>
void LaplaceProblem<dim>::run () 
{
  std::cout << "Solving problem in " << dim << " space dimensions." << std::endl;
  
  make_grid_and_dofs();
  assemble_system ();
  solve ();
  output_results ();
};

    

				 // And this is the main function. It
				 // also looks mostly like in the
				 // previous example:
int main () 
{
				   // In the previous example, we had
				   // the output from the linear
				   // solvers about the starting
				   // residual and the number of the
				   // iteration where convergence was
				   // detected. This can be suppressed
				   // like this:
  deallog.depth_console (0);
				   // The rationale here is the
				   // following: the deallog
				   // (i.e. deal-log, not de-allog)
				   // variable represents a stream to
				   // which some parts of the library
				   // write output. It redirects this
				   // output to the console and if
				   // required to a file. The output
				   // is nested in a way that each
				   // function can use a prefix string
				   // (separated by colons) for each
				   // line of output; if it calls
				   // another function, that may also
				   // use its prefix which is then
				   // printed after the one of the
				   // calling function. Since output
				   // from functions which are nested
				   // deep below is usually not as
				   // important as top-level output,
				   // you can give the deallog
				   // variable a maximal depth of
				   // nested output for output to
				   // console and file. The depth zero
				   // which we gave here means that no
				   // output is written.

				   // After having done this
				   // administrative stuff, we can go
				   // on just as before: define one of
				   // these top-level objects and
				   // transfer control to
				   // it. Actually, now is the point
				   // where we have to tell the
				   // compiler which dimension we
				   // would like to use; all functions
				   // up to now including the classes
				   // were only templates and nothing
				   // has been compiled by now, but by
				   // declaring the following objects,
				   // the compiler will start to
				   // compile all the functions at the
				   // top using the template parameter
				   // replaced with a concrete value.
				   //
				   // For demonstration, we will first
				   // let the whole thing run in 2D
				   // and then in 3D:
  LaplaceProblem<2> laplace_problem_2d;
  laplace_problem_2d.run ();

  LaplaceProblem<3> laplace_problem_3d;
  laplace_problem_3d.run ();
  
  return 0;
};
