/* $Id$ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 1999 */
/* Program is based on /examples/step-3 
   Purpose: compare the results when using a normal matrix and
   a block matrix
*/

#include <base/logstream.h>
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
#include <lac/block_vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/vector_memory.h>
#include <lac/precondition.h>

#include <numerics/data_out.h>
#include <fstream>


template <class Vector, class Matrix, class Sparsity>
class LaplaceProblem 
{
  public:
    LaplaceProblem ();

    void run ();
    void reinit_sparsity ();
    void reinit_vectors ();
    
    Vector       solution;

  private:
    void make_grid_and_dofs ();
    void assemble_system ();
    void solve ();

    Triangulation<2>     triangulation;
    FEQ1<2>              fe;
    DoFHandler<2>        dof_handler;

    Sparsity      sparsity_pattern;
    Matrix system_matrix;

    Vector       system_rhs;
};


template <class Vector, class Matrix, class Sparsity>
LaplaceProblem<Vector,Matrix,Sparsity>::LaplaceProblem () :
		dof_handler (triangulation)
{};


template <class Vector, class Matrix, class Sparsity>
void LaplaceProblem<Vector,Matrix,Sparsity>::make_grid_and_dofs ()
{
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (5);
  deallog << "Number of active cells: "
       << triangulation.n_active_cells()
       << endl;
  deallog << "Total number of cells: "
       << triangulation.n_cells()
       << endl;
  
  dof_handler.distribute_dofs (fe);

  deallog << "Number of degrees of freedom: "
       << dof_handler.n_dofs()
       << endl;

  reinit_sparsity ();
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);
  reinit_vectors ();
};


template <>
void LaplaceProblem<Vector<double>,SparseMatrix<double>,SparsityPattern>::reinit_sparsity () 
{
  sparsity_pattern.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
};



template <>
void LaplaceProblem<Vector<double>,SparseMatrix<double>,SparsityPattern>::reinit_vectors () 
{
  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
};



template <>
void LaplaceProblem<Vector<float>,SparseMatrix<float>,SparsityPattern>::reinit_sparsity () 
{
  sparsity_pattern.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
};



template <>
void LaplaceProblem<Vector<float>,SparseMatrix<float>,SparsityPattern>::reinit_vectors () 
{
  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
};



template <>
void LaplaceProblem<BlockVector<2,double>,BlockSparseMatrix<double,2,2>,BlockSparsityPattern<2,2> >::reinit_sparsity () 
{
  const unsigned int n_dofs = dof_handler.n_dofs();
  const unsigned int block_size[2] = { n_dofs/3, n_dofs - n_dofs/3 };

  for (unsigned int i=0; i<2; ++i)
    for (unsigned int j=0; j<2; ++j)
      sparsity_pattern.block(i,j).reinit (block_size[i], block_size[j],
					  dof_handler.max_couplings_between_dofs());
  sparsity_pattern.collect_sizes ();
};



template <>
void LaplaceProblem<BlockVector<2,double>,BlockSparseMatrix<double,2,2>,BlockSparsityPattern<2,2> >::reinit_vectors () 
{
  const unsigned int n_dofs = dof_handler.n_dofs();
  const unsigned int block_size_[2] = { n_dofs/3, n_dofs - n_dofs/3 };
  const vector<unsigned int> block_size (&block_size_[0],
					 &block_size_[2]);
					 
  solution.reinit (block_size);
  system_rhs.reinit (block_size);
};



template <>
void LaplaceProblem<BlockVector<3,double>,BlockSparseMatrix<double,3,3>,BlockSparsityPattern<3,3> >::reinit_sparsity () 
{
  const unsigned int n_dofs = dof_handler.n_dofs();
  const unsigned int block_size[3] = { n_dofs/5, n_dofs/7, n_dofs - n_dofs/5 - n_dofs/7 };

  for (unsigned int i=0; i<3; ++i)
    for (unsigned int j=0; j<3; ++j)
      sparsity_pattern.block(i,j).reinit (block_size[i], block_size[j],
					  dof_handler.max_couplings_between_dofs());
  sparsity_pattern.collect_sizes ();
};



template <>
void LaplaceProblem<BlockVector<3,double>,BlockSparseMatrix<double,3,3>,BlockSparsityPattern<3,3> >::reinit_vectors () 
{
  const unsigned int n_dofs = dof_handler.n_dofs();
  const unsigned int block_size_[3] = { n_dofs/5, n_dofs/7, n_dofs - n_dofs/5 - n_dofs/7 };
  const vector<unsigned int> block_size (&block_size_[0],
					 &block_size_[3]);
					 
  solution.reinit (block_size);
  system_rhs.reinit (block_size);
};



template <class Vector, class Matrix, class Sparsity>
void LaplaceProblem<Vector,Matrix,Sparsity>::assemble_system () 
{
  QGauss2<2>  quadrature_formula;
  FEValues<2> fe_values (fe, quadrature_formula, 
			 UpdateFlags(update_values    |
				     update_gradients |
				     update_JxW_values));

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.n_quadrature_points;

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  ::Vector<double>     cell_rhs (dofs_per_cell);

  vector<unsigned int> local_dof_indices (dofs_per_cell);

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


  map<unsigned int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
					    0,
					    ZeroFunction<2>(),
					    boundary_values);
  MatrixTools<2>::apply_boundary_values (boundary_values,
					 system_matrix,
					 solution,
					 system_rhs);
};


template <class Vector, class Matrix, class Sparsity>
void LaplaceProblem<Vector,Matrix,Sparsity>::solve () 
{
  SolverControl           solver_control (1000, 1e-12, false, false);
  PrimitiveVectorMemory<Vector> vector_memory;
  SolverCG<Vector>        cg (solver_control, vector_memory);

  PreconditionRelaxation<Matrix,Vector> preconditioner
    (system_matrix, &Matrix::precondition_Jacobi, 0.8);
  
  cg.solve (system_matrix, solution, system_rhs,
	    preconditioner);
};


template <class Vector, class Matrix, class Sparsity>
void LaplaceProblem<Vector,Matrix,Sparsity>::run () 
{
  make_grid_and_dofs ();
  assemble_system ();
  solve ();

  for (unsigned int i=0; i<solution.size(); ++i)
    deallog << typeid(Vector).name ()
	    << ' '
	    << typeid(Matrix).name ()
	    << '-'
	    << i << ' ' << solution(i) << endl;
};

    

int main () 
{
  ofstream logfile("block_matrices.output");
  logfile.precision(4);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  

				   // vector of solution vectors
  vector<vector<double> > solutions;

  if (true)
    {
      LaplaceProblem<Vector<double>,SparseMatrix<double>,SparsityPattern>
	laplace_problem;  
      laplace_problem.run ();
      
      solutions.push_back (vector<double>());
      solutions.back().resize (laplace_problem.solution.size());
      for (unsigned int i=0; i<laplace_problem.solution.size(); ++i)
	solutions.back()[i] = laplace_problem.solution(i);
    };
  
  if (true)
    {
      LaplaceProblem<Vector<float>,SparseMatrix<float>,SparsityPattern>
	laplace_problem;  
      laplace_problem.run ();
      
      solutions.push_back (vector<double>());
      solutions.back().resize (laplace_problem.solution.size());
      for (unsigned int i=0; i<laplace_problem.solution.size(); ++i)
	solutions.back()[i] = laplace_problem.solution(i);
    };
  
  if (true)
    {
      LaplaceProblem<BlockVector<2,double>,BlockSparseMatrix<double,2,2>,BlockSparsityPattern<2,2> >
	laplace_problem;  
      laplace_problem.run ();
      
      solutions.push_back (vector<double>());
      solutions.back().resize (laplace_problem.solution.size());
      for (unsigned int i=0; i<laplace_problem.solution.size(); ++i)
	solutions.back()[i] = laplace_problem.solution(i);
    };

  if (true)
    {
      LaplaceProblem<BlockVector<3,double>,BlockSparseMatrix<double,3,3>,BlockSparsityPattern<3,3> >
	laplace_problem;  
      laplace_problem.run ();
      
      solutions.push_back (vector<double>());
      solutions.back().resize (laplace_problem.solution.size());
      for (unsigned int i=0; i<laplace_problem.solution.size(); ++i)
	solutions.back()[i] = laplace_problem.solution(i);
    };

  const unsigned int n_datasets = solutions.size();
  deallog << "Checking " << n_datasets << " data sets." << endl;
  
  for (unsigned int i=1; i<n_datasets; ++i)
    Assert (solutions[i].size() == solutions[i].size(),
	    ExcInternalError());
  
  logfile.precision(16);
  for (unsigned int i=1; i<n_datasets; ++i)
    {
				       // relative accuracy. data set
				       // 1 is computed using floats
				       // instead of doubles, so lower
				       // our requirements
      const double accuracy = (i==1 ? 1e-6 : 1e-12);
      
      for (unsigned int j=0; j<solutions[0].size(); ++j)
	if ( fabs(solutions[i][j] - solutions[0][j]) >
	     accuracy*fabs(solutions[i][j] + solutions[0][j]))
	  {
	    deallog << "Discrepancy: i=" << i << ", j=" << j
		    << ", sol[i][j]=" << solutions[i][j]
		    << ", sol[0][j]=" << solutions[0][j]
		    << endl;
	    deallog << flush;
	    Assert (false, ExcInternalError());
	  };
    };
  
    
  return 0;
};
