/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */



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
#include <grid/grid_refinement.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_values.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <fe/fe_lib.lagrange.h>
#include <grid/grid_out.h>
#include <dofs/dof_constraints.h>
#include <numerics/error_estimator.h>

#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>
#include <multigrid/mg_dof_tools.h>
#include <multigrid/mg_base.h>
#include <multigrid/mg_smoother.h>
#include <multigrid/multigrid.h>

#include <lac/solver_richardson.h>

#include <fstream>



template <int dim>
class LaplaceProblem 
{
  public:
    LaplaceProblem ();
    ~LaplaceProblem ();
    void run ();
    
  private:
    void setup_system ();
    void assemble_system ();
    void solve ();
    void refine_grid ();
    void output_results (const unsigned int cycle) const;

    Triangulation<dim>   triangulation;
    MGDoFHandler<dim>    mg_dof_handler;

    FEQ1<dim>            fe;

    ConstraintMatrix     hanging_node_constraints;

    SparsityPattern      global_sparsity_pattern;
    SparseMatrix<double> global_system_matrix;

    MGLevelObject<SparsityPattern>       level_sparsity_patterns;
    MGLevelObject<SparseMatrix<double> > level_system_matrices;
    
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
	  ExcDimensionMismatch (values.size(), n_points));
  
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




class MGSmootherLAC : public MGSmootherBase
{
  private:
    SmartPointer<MGLevelObject<SparseMatrix<double> > >matrices;
  public:
    MGSmootherLAC(MGLevelObject<SparseMatrix<double> >&);
    
    virtual void smooth (const unsigned int level,
			 Vector<double> &u,
			 const Vector<double> &rhs) const;    
};


MGSmootherLAC::MGSmootherLAC(MGLevelObject<SparseMatrix<double> >& matrix)
		:
		matrices(&matrix)
{}


void
MGSmootherLAC::smooth (const unsigned int level,
		       Vector<double> &u,
		       const Vector<double> &rhs) const
{
  SolverControl control(2,1.e-300,false,false);
  PrimitiveVectorMemory<> mem;
  SolverRichardson<> rich(control, mem);
  PreconditionSSOR<> prec;
  prec.initialize((*matrices)[level], 1.);

  rich.solve((*matrices)[level], u, rhs, prec);
}



template <int dim>
LaplaceProblem<dim>::LaplaceProblem () :
		mg_dof_handler (triangulation)
{};



template <int dim>
LaplaceProblem<dim>::~LaplaceProblem () 
{
  mg_dof_handler.clear ();
};



template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
  mg_dof_handler.distribute_dofs (fe);

  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (mg_dof_handler,
					   hanging_node_constraints);
  hanging_node_constraints.close ();
  global_sparsity_pattern.reinit (mg_dof_handler.DoFHandler<dim>::n_dofs(),
				  mg_dof_handler.DoFHandler<dim>::n_dofs(),
				  mg_dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (mg_dof_handler, global_sparsity_pattern);
  hanging_node_constraints.condense (global_sparsity_pattern);
  global_sparsity_pattern.compress();

  global_system_matrix.reinit (global_sparsity_pattern);

  solution.reinit (mg_dof_handler.DoFHandler<dim>::n_dofs());
  system_rhs.reinit (mg_dof_handler.DoFHandler<dim>::n_dofs());


  const unsigned int n_levels = triangulation.n_levels();
  level_system_matrices.resize (0, n_levels);
  level_sparsity_patterns.resize (0, n_levels);
  
  for (unsigned int level=0; level<n_levels; ++level) 
    {
      level_sparsity_patterns[level].reinit (mg_dof_handler.n_dofs(level),
					     mg_dof_handler.n_dofs(level),
					     mg_dof_handler.max_couplings_between_dofs()); //xxx
      MGDoFTools::make_sparsity_pattern (mg_dof_handler,
					 level_sparsity_patterns[level],
					 level);
      level_sparsity_patterns[level].compress();

      level_system_matrices[level].reinit (level_sparsity_patterns[level]);
    };
};



template <int dim>
void LaplaceProblem<dim>::assemble_system () 
{  
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

  vector<unsigned int> local_dof_indices (dofs_per_cell);

				   // FIX
  vector<double>       coefficient_values (n_q_points, 1.0);

				   // not only active cells
  MGDoFHandler<dim>::cell_iterator cell = mg_dof_handler.begin(),
				   endc = mg_dof_handler.end();
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

				       // FIX
//      coefficient.value_list (q_points, coefficient_values);
      
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  {
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      cell_matrix(i,j) += (coefficient_values[q_point] *
				   (shape_grads[i][q_point]    *
				    shape_grads[j][q_point] +
				    shape_values(i,q_point)    *
				    shape_values(j,q_point)  )   *
				   JxW_values[q_point]);

	    cell_rhs(i) += (shape_values (i,q_point) *
			    sin(4*sqrt(q_points[q_point].square())) *
			    fe_values.JxW (q_point));
	  };


      cell->get_mg_dof_indices (local_dof_indices);
      const unsigned int level = cell->level();
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  level_system_matrices[level].add (local_dof_indices[i],
					    local_dof_indices[j],
					    cell_matrix(i,j));
      
				       // if active, then also into
				       // global matrix
      if (cell->active())
	{
	  cell->get_dof_indices (local_dof_indices);
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    {
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
		global_system_matrix.add (local_dof_indices[i],
					  local_dof_indices[j],
					  cell_matrix(i,j));
	      
	      system_rhs(local_dof_indices[i]) += cell_rhs(i);
	    };
	};
    };

  hanging_node_constraints.condense (global_system_matrix);
  hanging_node_constraints.condense (system_rhs);

//    map<unsigned int,double> boundary_values;
//    VectorTools::interpolate_boundary_values (mg_dof_handler,
//  					    0,
//  					    ZeroFunction<dim>(),
//  					    boundary_values);
//    MatrixTools<dim>::apply_boundary_values (boundary_values,
//  					   global_system_matrix,
//  					   solution,
//  					   system_rhs);
};



template <int dim>
void LaplaceProblem<dim>::solve () 
{

    {
      SolverControl           solver_control (1000, 1e-12);
      PrimitiveVectorMemory<> vector_memory;
      SolverCG<>              cg (solver_control, vector_memory);

      SolverControl           coarse_grid_solver_control (1000, 1e-12);
      PrimitiveVectorMemory<> coarse_grid_vector_memory;
      
      SolverCG<>              coarse_grid_cg (coarse_grid_solver_control,
					      coarse_grid_vector_memory);
      
//        PreconditionRelaxation<>
//  	coarse_grid_solver_preconditioner(level_system_matrices[level_system_matrices.get_minlevel()],
//  					  &SparseMatrix<double>::template precondition_SSOR<double>,
//  					  1.2);
      PreconditionIdentity coarse_grid_solver_preconditioner;
      
      MGCoarseGridLACIteration<SolverCG<>, SparseMatrix<double>, PreconditionIdentity>
	coarse_grid_solver (coarse_grid_cg,
			    level_system_matrices[level_system_matrices.get_minlevel()],
			    coarse_grid_solver_preconditioner);
      
      MGSmootherLAC      smoother (level_system_matrices);
      MGTransferPrebuilt grid_transfer;
      grid_transfer.build_matrices (mg_dof_handler);
      
      Multigrid<2> multigrid (mg_dof_handler,
			      hanging_node_constraints,
			      level_sparsity_patterns,
			      level_system_matrices,
			      grid_transfer);
      
      PreconditionMG<Multigrid<2> >
	mg_precondition (multigrid, smoother, smoother, coarse_grid_solver);

      solution.clear ();
      cg.solve (global_system_matrix, solution, system_rhs,
		mg_precondition);

      cout << "   MG Outer iterations:       " << solver_control.last_step()
	   << endl;

      cout << "   MG Total inner iterations: " << coarse_grid_solver_control.last_step()
	   << endl;
    };
  
    {
      SolverControl           solver_control (1000, 1e-12);
      PrimitiveVectorMemory<> vector_memory;
      SolverCG<>              cg (solver_control, vector_memory);

      PreconditionSSOR<> preconditioner;
      preconditioner.initialize(global_system_matrix, 1.2);
      
      solution.clear ();
      cg.solve (global_system_matrix, solution, system_rhs,
		preconditioner);
      
      cout << "   CG Outer iterations:       " << solver_control.last_step()
	   << endl;
    };

  hanging_node_constraints.distribute (solution);
};


template <int dim>
void LaplaceProblem<dim>::refine_grid ()
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

  KellyErrorEstimator<dim>::FunctionMap neumann_boundary;
  KellyErrorEstimator<dim>::estimate (mg_dof_handler,
				      QGauss3<dim-1>(),
				      neumann_boundary,
				      solution,
				      estimated_error_per_cell);

  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
						   estimated_error_per_cell,
						   0.3, 0.03);
  triangulation.execute_coarsening_and_refinement ();
};



template <int dim>
void LaplaceProblem<dim>::output_results (const unsigned int cycle) const
{
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
	  GridGenerator::hyper_cube (triangulation);
	  triangulation.refine_global (1);
	}
      else
	{
	  refine_grid ();
	};
      

      cout << "   Number of active cells:       "
	   << triangulation.n_active_cells()
	   << endl;

      setup_system ();

      cout << "   Number of degrees of freedom: "
	   << mg_dof_handler.DoFHandler<dim>::n_dofs()
	   << endl;
      
      assemble_system ();
      solve ();
      output_results (cycle);

  typename DataOut<dim>::EpsFlags eps_flags;
  eps_flags.z_scaling = 4;
  
  DataOut<dim> data_out;
  data_out.set_flags (eps_flags);

  data_out.attach_dof_handler (mg_dof_handler);
  data_out.add_data_vector (solution, "solution");
  data_out.build_patches ();
  
  ofstream output ("final-solution.eps");
  data_out.write_eps (output);
    };
};


    
int main () 
{
  try
    {
      deallog.depth_console (0);

      LaplaceProblem<2> laplace_problem_2d;
      laplace_problem_2d.run ();
    }
  catch (exception &exc)
    {
      cerr << endl << endl
	   << "----------------------------------------------------"
	   << endl;
      cerr << "Exception on processing: " << endl
	   << exc.what() << endl
	   << "Aborting!" << endl
	   << "----------------------------------------------------"
	   << endl;
      return 1;
    }
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

  return 0;
};
