/* $Id$ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 2000 */


				 // still unfinished
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
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_lib.lagrange.h>
#include <fe/fe_values.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>

				 // out statt in
#include <grid/grid_out.h>

#include <grid/tria_boundary_lib.h>

				 //...
#include <dofs/dof_constraints.h>
#include <numerics/error_estimator.h>

#include <fstream>



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

				     // This is the only addition to
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
				 // scheme.
				 // ...
template <int dim>
void LaplaceProblem<dim>::refine_grid ()
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

  KellyErrorEstimator<dim>::FunctionMap neumann_boundary;
  KellyErrorEstimator<dim>::estimate (dof_handler,
				      QGauss3<dim-1>(),
				      neumann_boundary,
				      solution,
				      estimated_error_per_cell);

  triangulation.refine_and_coarsen_fixed_number (estimated_error_per_cell,
						 0.3, 0.03);
  triangulation.execute_coarsening_and_refinement ();
};

template <int dim>
void LaplaceProblem<dim>::output_results (const unsigned int cycle) const
{
				   // ...
  string filename = "grid-";
  filename += ('0' + cycle);
  Assert (cycle < 10, ExcInternalError());
  
  filename += ".eps";
  ofstream output (filename.c_str());

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
					   //...
	  GridGenerator::hyper_ball (triangulation);

	  static const HyperBallBoundary<dim> boundary;
	  triangulation.set_boundary (0, boundary);

	  triangulation.refine_global (1);
	}
      else
					 // ...
	refine_grid ();

      cout << "   Number of active cells: "
	   << triangulation.n_active_cells()
	   << endl;

      setup_system ();
      assemble_system ();
      solve ();
      output_results (cycle);
    };

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

    

int main () 
{
  deallog.depth_console (0);

  LaplaceProblem<2> laplace_problem_2d;
  laplace_problem_2d.run ();
  
  return 0;
};
