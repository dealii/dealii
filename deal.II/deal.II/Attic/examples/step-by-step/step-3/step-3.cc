#include <base/function.h>
#include <base/quadrature_lib.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_lib.lagrange.h>
#include <fe/fe_values.h>
#include <dofs/dof_tools.h>
#include <numerics/data_out.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>

#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/vector_memory.h>
#include <lac/precondition.h>

#include <fstream>


class LaplaceProblem 
{
  public:
    LaplaceProblem ();
    
    void make_grid_and_dofs ();
    void assemble_system ();
    void solve ();
    void output_results ();

    void run ();
    
  private:
    Triangulation<2> triangulation;
    FEQ1<2>          fe;
    DoFHandler<2>    dof_handler;
    
    SparseMatrixStruct   sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       system_rhs;

    Vector<double>       solution;
};


LaplaceProblem::LaplaceProblem () :
		dof_handler (triangulation)
{};



void LaplaceProblem::make_grid_and_dofs ()
{
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (5);

  dof_handler.distribute_dofs (fe);

  sparsity_pattern.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
};



void LaplaceProblem::assemble_system () 
{
  QGauss3<2>  quadrature_formula;
  FEValues<2> fe_values (fe, quadrature_formula, 
			 UpdateFlags(update_gradients |
				     update_JxW_values));

  const unsigned int n_q_points    = quadrature_formula.n_quadrature_points;
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  
  FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs (dofs_per_cell);
  
  vector<int>        local_dof_indices (dofs_per_cell);

  DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(),
				      endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      
      cell_matrix.clear ();
      cell_rhs.clear ();
      
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	    cell_matrix(i,j) += (fe_values.shape_grad (i, q_point) *
				 fe_values.shape_grad (j, q_point) *
				 fe_values.JxW (q_point));

      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	  cell_rhs(i) += (fe_values.shape_value (i, q_point) *
			  1 *
			  fe_values.JxW (q_point));

      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  system_matrix.add (local_dof_indices[i],
			     local_dof_indices[j],
			     cell_matrix(i,j));
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	system_rhs(local_dof_indices[i]) += cell_rhs(i);
    };


  map<int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
					       0,
					       ZeroFunction<2>(),
					       boundary_values);
  MatrixTools<2>::apply_boundary_values (boundary_values,
					 system_matrix,
					 solution,
					 system_rhs);
};



void LaplaceProblem::solve () 
{
  SolverControl           solver_control (1000, 1e-12);
  PrimitiveVectorMemory<> vector_memory;
  SolverCG<>              cg (solver_control, vector_memory);
  
  cg.solve (system_matrix, solution, system_rhs,
	    PreconditionIdentity());
};


void LaplaceProblem::output_results () 
{
  DataOut<2> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");
  data_out.build_patches ();
  
  ofstream output ("solution.gpl");
  data_out.write_gnuplot (output);
};


void LaplaceProblem::run () 
{
  make_grid_and_dofs();
  assemble_system ();
  solve ();
  output_results ();
};

    

int main () 
{
  LaplaceProblem laplace_problem;
  laplace_problem.run ();
  return 0;
};
