/* $Id$ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 2000 */

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
#include <fe/fe_values.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <fe/fe_lib.lagrange.h>
#include <dofs/dof_constraints.h>
#include <numerics/error_estimator.h>

#include <numerics/dof_renumbering.h>
#include <base/smartpointer.h>

#include <fstream>


template <int dim>
class LaplaceProblem 
{
  public:
    enum RefinementMode {
	  global_refinement, adaptive_refinement
    };
    
    LaplaceProblem (const FiniteElement<dim> &fe,
		    const RefinementMode      refinement_mode);
    ~LaplaceProblem ();
    void run ();
    
  private:
    void setup_system ();
    void assemble_system ();
    void solve ();
    void refine_grid ();
    void process_solution (const unsigned int cycle) const;

    Triangulation<dim>                      triangulation;
    DoFHandler<dim>                         dof_handler;
				     //...
    SmartPointer<const FiniteElement<dim> > fe;
    ConstraintMatrix                        hanging_node_constraints;

    SparsityPattern                         sparsity_pattern;
    SparseMatrix<double>                    system_matrix;

    Vector<double>                          solution;
    Vector<double>                          system_rhs;

    RefinementMode                          refinement_mode;
};



template <int dim>
class SolutionBase 
{
  protected:
    static const unsigned int n_source_centers = 3;    
    static const Point<dim>   source_centers[n_source_centers];
    static const double       width;
};


template <int dim>
class Solution : public Function<dim>,
		 protected SolutionBase<dim>
{
  public:
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
				    const unsigned int  component = 0) const;
};



template <int dim>
class RightHandSide : public Function<dim>,
		      protected SolutionBase<dim>
{
  public:
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};



template <>
const Point<1>
SolutionBase<1>::source_centers[SolutionBase<1>::n_source_centers]
= { Point<1>(-1.0 / 3.0), 
    Point<1>(0.0), 
    Point<1>(+1.0 / 3.0)   };

template <>
const Point<2>
SolutionBase<2>::source_centers[SolutionBase<2>::n_source_centers]
= { Point<2>(-0.5, +0.5), 
    Point<2>(-0.5, -0.5), 
    Point<2>(+0.5, -0.5)   };

template <int dim>
const double SolutionBase<dim>::width = 0.15;



template <int dim>
double Solution<dim>::value (const Point<dim>   &p,
			     const unsigned int) const
{
  double return_value = 0;
  for (unsigned int i=0; i<n_source_centers; ++i)
    {
      const Point<dim> shifted_point = p-source_centers[i];
      
      return_value += exp(-shifted_point.square() / (width*width));
    };
  
  return return_value;
};



template <int dim>
Tensor<1,dim> Solution<dim>::gradient (const Point<dim>   &p,
				       const unsigned int) const
{
  Tensor<1,dim> return_value;
  for (unsigned int i=0; i<n_source_centers; ++i)
    {
      const Point<dim> shifted_point = p-source_centers[i];
      
      return_value += (-2 / (width*width) *
		       exp(-shifted_point.square() / (width*width)) *
		       shifted_point);
    };
  
  return return_value;
};



template <int dim>
double RightHandSide<dim>::value (const Point<dim>   &p,
				  const unsigned int) const
{
  double return_value = 0;
  for (unsigned int i=0; i<n_source_centers; ++i)
    {
      const Point<dim> shifted_point = p-source_centers[i];
      
      return_value += ((2*dim - 4*shifted_point.square()/(width*width)) / (width*width) *
		       exp(-shifted_point.square() / (width*width)));
    };
  
  return return_value;
};



template <int dim>
LaplaceProblem<dim>::LaplaceProblem (const FiniteElement<dim> &fe,
				     const RefinementMode refinement_mode) :
		dof_handler (triangulation),
		fe (&fe),
		refinement_mode (refinement_mode)
{};



template <int dim>
LaplaceProblem<dim>::~LaplaceProblem () 
{
  dof_handler.clear ();
};



template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
  dof_handler.distribute_dofs (*fe);
				   // Renumber the degrees of freedom...
  DoFRenumbering::Cuthill_McKee (dof_handler);

  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
					   hanging_node_constraints);
  hanging_node_constraints.close ();

  sparsity_pattern.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  hanging_node_constraints.condense (sparsity_pattern);
  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
};



template <int dim>
void LaplaceProblem<dim>::assemble_system () 
{  
  QGauss3<dim>  quadrature_formula;
  FEValues<dim> fe_values (*fe, quadrature_formula, 
			   UpdateFlags(update_values    |
				       update_gradients |
				       update_q_points  |
				       update_JxW_values));

  const unsigned int   dofs_per_cell = fe->dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.n_quadrature_points;

  RightHandSide<dim>   right_hand_side;
  vector<double>       rhs_values (n_q_points);

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  vector<unsigned int> local_dof_indices (dofs_per_cell);

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

      right_hand_side.value_list (q_points, rhs_values);
      
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  {
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      cell_matrix(i,j) += (shape_grads[i][q_point] *
				   shape_grads[j][q_point] *
				   JxW_values[q_point]);

	    cell_rhs(i) += (shape_values (i,q_point) *
			    rhs_values [q_point] *
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

  hanging_node_constraints.condense (system_matrix);
  hanging_node_constraints.condense (system_rhs);

  map<unsigned int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
					    0,
					    Solution<dim>(),
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

  PreconditionRelaxation<>
    preconditioner(system_matrix,
		   &SparseMatrix<double>::template precondition_SSOR<double>,
		   1.2);

  cg.solve (system_matrix, solution, system_rhs,
	    preconditioner);

  hanging_node_constraints.distribute (solution);
};



template <int dim>
void LaplaceProblem<dim>::refine_grid ()
{
  switch (refinement_mode) 
    {
      case global_refinement:
      {
	triangulation.refine_global (1);
	break;
      };
       
      case adaptive_refinement:
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

	break;
      };
    };
};


#include <numerics/data_out.h>

template <int dim>
void LaplaceProblem<dim>::process_solution (const unsigned int cycle) const
{
  Vector<float> difference_per_cell (triangulation.n_active_cells());
  
  VectorTools::integrate_difference (dof_handler,
				     solution,
				     Solution<dim>(),
				     difference_per_cell,
				     QGauss3<dim>(),
				     L2_norm);
  const double L2_error = difference_per_cell.l2_norm();

  VectorTools::integrate_difference (dof_handler,
				     solution,
				     Solution<dim>(),
				     difference_per_cell,
				     QGauss3<dim>(),
				     H1_seminorm);
  const double H1_error = difference_per_cell.l2_norm();

  VectorTools::integrate_difference (dof_handler,
				     solution,
				     Solution<dim>(),
				     difference_per_cell,
				     QGauss3<dim>(),
				     Linfty_norm);
  const double Linfty_error = difference_per_cell.linfty_norm();
  
  cout << "Cycle " << cycle << ':' 
       << endl
       << "   Number of active cells:       "
       << triangulation.n_active_cells()
       << endl
       << "   Number of degrees of freedom: "
       << dof_handler.n_dofs()
       << endl;

  cout << "   L2     error: " << L2_error      << endl
       << "   H1     error: " << H1_error      << endl
       << "   Linfty error: " << Linfty_error  << endl;
};



template <int dim>
void LaplaceProblem<dim>::run () 
{
  for (unsigned int cycle=0; cycle<12; ++cycle)
    {
      if (cycle == 0)
	{
	  GridGenerator::hyper_cube (triangulation, -1, 1);
	  triangulation.refine_global (1);
	}
      else
	refine_grid ();      

      setup_system ();
      
      assemble_system ();
      solve ();
      process_solution (cycle);
    };
  
  string filename;
  switch (refinement_mode)
    {
      case global_refinement:
	    filename = "solution-global";
	    break;
      case adaptive_refinement:
	    filename = "solution-adaptive";
	    break;
      default:
	    Assert (false, ExcInternalError());
    };
  filename += ".gmv";
	    
  ofstream output (filename.c_str());


  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");
  data_out.build_patches ();
  data_out.write_gmv (output);
};



int main () 
{
  try
    {
      deallog.depth_console (0);

      FEQ1<2> fe;
      LaplaceProblem<2> laplace_problem_2d (fe, LaplaceProblem<2>::adaptive_refinement);
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
