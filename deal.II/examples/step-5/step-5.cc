/* $Id$ */

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
#include <fe/fe_lib.lagrange.h>
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

#include <base/logstream.h>


template <int dim>
class LaplaceProblem 
{
  public:
    LaplaceProblem ();
    void run ();
    
  private:
    void make_grid_and_dofs (const unsigned int refinement);
    void assemble_system ();
    void solve ();
    void output_results ();
    void clear ();

    Triangulation<dim>   triangulation;
    FEQ1<dim>            fe;
    DoFHandler<dim>      dof_handler;

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
    return 10;
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
    if (points[i].square() < 0.5*0.5)
      values[i] = 10;
    else
      values[i] = 1;
};



template <int dim>
LaplaceProblem<dim>::LaplaceProblem () :
		dof_handler (triangulation)
{};



template <int dim>
void LaplaceProblem<dim>::make_grid_and_dofs (const unsigned int refinement)
{
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (refinement);
  
  cout << "   Number of active cells: "
       << triangulation.n_active_cells()
       << endl
       << "   Total number of cells: "
       << triangulation.n_cells()
       << endl;

  dof_handler.distribute_dofs (fe);

  cout << "   Number of degrees of freedom: "
       << dof_handler.n_dofs()
       << endl;

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
				 // time, and it is worth the effort
				 // to use one or two optimizations at
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

				   // ...
  vector<double>     coefficient_values (n_q_points);

  DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
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
				       // symbolized by the ampersand,
				       // and that it must be a
				       // constant reference, since
				       // only read-only access is
				       // granted):
      const FullMatrix<double> 
	& shape_values = fe_values.get_shape_values();
				       // Instead of writing
				       // fe_values.shape_value(j,q)
				       // we can now write
				       // shape_values(j,q), i.e. the
				       // function call needed
				       // previously for each access
				       // has been otimized away.
				       //
				       // There are alike functions
				       // for almost all data elements
				       // in the FEValues class. The
				       // gradient are accessed as
				       // follows:
      const vector<vector<Tensor<1,dim> > >
	& shape_grads  = fe_values.get_shape_grads();
				       // The data type looks a bit
				       // unwieldy, since each entry
				       // in the matrix (j,q) now
				       // needs to be the gradient of
				       // the shape function, which is
				       // a vector.
				       //
				       // Similarly, access to the
				       // place where quadrature
				       // points and the determinants
				       // of the Jacobian matrices
				       // times the weights of the
				       // respective quadrature points
				       // are stored, can be obtained
				       // like this:
      const vector<double>
	& JxW_values   = fe_values.get_JxW_values();
      const vector<Point<dim> >
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
				       // below). Unfortunately,
				       // things became a bit
				       // inconsistent, since the
				       // shape values are accessed
				       // via the FullMatrix operator
				       // (), i.e. using parentheses,
				       // while all the other fields
				       // are accessed through vector
				       // operator [], i.e. using
				       // brackets. This is due to
				       // historical reasons and
				       // frequently leads to a bit of
				       // confusion, but since the
				       // places where this happens
				       // are few in well-written
				       // programs, this is not too
				       // big a problem.

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

				   // Again use zero boundary values:
  map<int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
					    0,
					    ZeroFunction<dim>(),
					    boundary_values);
  MatrixTools<dim>::apply_boundary_values (boundary_values,
					   system_matrix,
					   solution,
					   system_rhs);
};



template <int dim>
void LaplaceProblem<dim>::solve () 
{
  SolverControl           solver_control (1000, 1e-12);
  PrimitiveVectorMemory<> vector_memory;
  SolverCG<>              cg (solver_control, vector_memory);

				   // ...
  PreconditionRelaxation<>
    preconditioner(system_matrix,
		   &SparseMatrix<double>::template precondition_SSOR<double>,
		   1.2);

  cg.solve (system_matrix, solution, system_rhs,
	    preconditioner);

  cout << "   " << solver_control.last_step()
       << " CG iterations needed to obtain convergence."
       << endl;
};



template <int dim>
void LaplaceProblem<dim>::output_results () 
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");

  data_out.build_patches ();

  ofstream output (dim == 2 ?
		   "solution-2d.gmv" :
		   "solution-3d.gmv");
				   // ...
  data_out.write_gnuplot (output);
};



template <int dim>
void LaplaceProblem<dim>::clear () 
{
  system_rhs.reinit (0);
  solution.reinit (0);
  system_matrix.reinit ();
  sparsity_pattern.reinit (0, 0, 0);
  dof_handler.clear ();
  triangulation.clear ();
};



template <int dim>
void LaplaceProblem<dim>::run () 
{
  cout << "Solving problem in " << dim << " space dimensions." << endl;
  
  for (unsigned int refinement=0; refinement<7; ++refinement)
    {
      cout << "Refinement step: " << refinement << endl;
      
      make_grid_and_dofs(refinement);
      assemble_system ();
      solve ();
      output_results ();

      clear ();
    };
};

    

int main () 
{
  deallog.depth_console (0);

  LaplaceProblem<2> laplace_problem_2d;
  laplace_problem_2d.run ();
  
  return 0;
};
