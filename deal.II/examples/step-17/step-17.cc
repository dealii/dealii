/* $Id$ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 2000 */

/*    $Id$       */
/*    Version: $Name$                                          */
/*                                                                */
/*    Copyright (C) 2004 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>

                                 // xxx
#include <lac/petsc_vector.h>
#include <lac/petsc_sparse_matrix.h>
#include <lac/petsc_solver.h>
#include <lac/petsc_precondition.h>

#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_refinement.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_values.h>
#include <fe/fe_system.h>
#include <fe/fe_q.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <dofs/dof_constraints.h>
#include <numerics/error_estimator.h>

                                 // xxx
#include <grid/grid_tools.h>

#include <fstream>
#include <iostream>


template <int dim>
class ElasticProblem 
{
  public:
    ElasticProblem ();
    ~ElasticProblem ();
    void run ();
    
  private:
    void setup_system ();
    void assemble_system ();
    void solve ();
    void refine_grid ();
    void output_results (const unsigned int cycle) const;

    Triangulation<dim>   triangulation;
    DoFHandler<dim>      dof_handler;

    FESystem<dim>        fe;

    ConstraintMatrix     hanging_node_constraints;

                                     // xxx no sparsity
    PETScWrappers::SparseMatrix system_matrix;

    PETScWrappers::Vector       solution;
    PETScWrappers::Vector       system_rhs;

                                     // xxx
    const unsigned int n_partitions;
    const unsigned int this_partition;
};


template <int dim>
class RightHandSide :  public Function<dim> 
{
  public:
    RightHandSide ();
    
    virtual void vector_value (const Point<dim> &p,
			       Vector<double>   &values) const;

    virtual void vector_value_list (const std::vector<Point<dim> > &points,
				    std::vector<Vector<double> >   &value_list) const;
};


template <int dim>
RightHandSide<dim>::RightHandSide () :
		Function<dim> (dim)
{}


template <int dim>
inline
void RightHandSide<dim>::vector_value (const Point<dim> &p,
				       Vector<double>   &values) const 
{
  Assert (values.size() == dim, 
	  ExcDimensionMismatch (values.size(), dim));
  Assert (dim >= 2, ExcInternalError());
  
  Point<dim> point_1, point_2;
  point_1(0) = 0.5;
  point_2(0) = -0.5;
  
  if (((p-point_1).square() < 0.2*0.2) ||
      ((p-point_2).square() < 0.2*0.2))
    values(0) = 1;
  else
    values(0) = 0;
  
  if (p.square() < 0.2*0.2)
    values(1) = 1;
  else
    values(1) = 0;    
}



template <int dim>
void RightHandSide<dim>::vector_value_list (const std::vector<Point<dim> > &points,
					    std::vector<Vector<double> >   &value_list) const 
{
  const unsigned int n_points = points.size();

  Assert (value_list.size() == n_points, 
	  ExcDimensionMismatch (value_list.size(), n_points));

  for (unsigned int p=0; p<n_points; ++p)
    RightHandSide<dim>::vector_value (points[p],
				      value_list[p]);
}




template <int dim>
ElasticProblem<dim>::ElasticProblem ()
                :
		dof_handler (triangulation),
		fe (FE_Q<dim>(1), dim),
                                                 // xxx
                n_partitions (8),
                this_partition (0)
{}



template <int dim>
ElasticProblem<dim>::~ElasticProblem () 
{
  dof_handler.clear ();
}


template <int dim>
void ElasticProblem<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);
  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
					   hanging_node_constraints);
  hanging_node_constraints.close ();

                                   // no sparsity pattern
  system_matrix.reinit (dof_handler.n_dofs(),
                        dof_handler.n_dofs(),
                        dof_handler.max_couplings_between_dofs());

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
}


template <int dim>
void ElasticProblem<dim>::assemble_system () 
{
                                   // move to front
  std::map<unsigned int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
					    0,
					    ZeroFunction<dim>(dim),
					    boundary_values);

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

  std::vector<double>     lambda_values (n_q_points);
  std::vector<double>     mu_values (n_q_points);

  ConstantFunction<dim> lambda(1.), mu(1.);

  RightHandSide<dim>      right_hand_side;
  std::vector<Vector<double> > rhs_values (n_q_points,
					   Vector<double>(dim));


  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
						 endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      cell_matrix.clear ();
      cell_rhs.clear ();

      fe_values.reinit (cell);
      
      lambda.value_list (fe_values.get_quadrature_points(), lambda_values);
      mu.value_list     (fe_values.get_quadrature_points(), mu_values);

      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  const unsigned int 
	    component_i = fe.system_to_component_index(i).first;
	  
	  for (unsigned int j=0; j<dofs_per_cell; ++j) 
	    {
	      const unsigned int 
		component_j = fe.system_to_component_index(j).first;
	      
	      for (unsigned int q_point=0; q_point<n_q_points;
		   ++q_point)
		{
		  cell_matrix(i,j) 
		    += 
		    (
		      (fe_values.shape_grad(i,q_point)[component_i] *
		       fe_values.shape_grad(j,q_point)[component_j] *
		       lambda_values[q_point])
		      +
		      (fe_values.shape_grad(i,q_point)[component_j] *
		       fe_values.shape_grad(j,q_point)[component_i] *
		       mu_values[q_point])
		      +
		      ((component_i == component_j) ?
		       (fe_values.shape_grad(i,q_point) *
			fe_values.shape_grad(j,q_point) *
			mu_values[q_point])  :
		       0)
		    )
		    *
		    fe_values.JxW(q_point);
		};
	    };
	};

      right_hand_side.vector_value_list (fe_values.get_quadrature_points(),
					 rhs_values);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  const unsigned int 
	    component_i = fe.system_to_component_index(i).first;
	  
	  for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	    cell_rhs(i) += fe_values.shape_value(i,q_point) *
			   rhs_values[q_point](component_i) *
			   fe_values.JxW(q_point);
	};

      cell->get_dof_indices (local_dof_indices);

                                       //xxx
      MatrixTools::local_apply_boundary_values (boundary_values,
                                                local_dof_indices,
                                                cell_matrix,
                                                cell_rhs,
                                                false);

                                       // xxx
      hanging_node_constraints
        .distribute_local_to_global (cell_matrix,
                                     local_dof_indices,
                                     system_matrix);

      hanging_node_constraints
        .distribute_local_to_global (cell_rhs,
                                     local_dof_indices,
                                     system_rhs);
    }

                                   //xxx no condense necessary, no apply_b_v
                                   //either


                                   // xxx
  system_matrix.compress ();
  system_rhs.compress ();
}



template <int dim>
void ElasticProblem<dim>::solve () 
{
                                   // xxx
  SolverControl           solver_control (1000, 1e-10);
  PETScWrappers::SolverCG cg (solver_control);

  PETScWrappers::PreconditionSSOR preconditioner(system_matrix, 1.2);

  cg.solve (system_matrix, solution, system_rhs,
	    preconditioner);

  hanging_node_constraints.distribute (solution);
}



template <int dim>
void ElasticProblem<dim>::refine_grid ()
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

  typename FunctionMap<dim>::type neumann_boundary;
  KellyErrorEstimator<dim>::estimate (dof_handler,
				      QGauss2<dim-1>(),
				      neumann_boundary,
				      solution,
				      estimated_error_per_cell);

  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
						   estimated_error_per_cell,
						   0.3, 0.03);

  triangulation.execute_coarsening_and_refinement ();

                                   // xxx
  GridTools::partition_triangulation (n_partitions, triangulation);
}


template <int dim>
void ElasticProblem<dim>::output_results (const unsigned int cycle) const
{
  std::string filename = "solution-";
  filename += ('0' + cycle);
  Assert (cycle < 10, ExcInternalError());
  
  filename += ".gmv";
  std::ofstream output (filename.c_str());

  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);

 

  std::vector<std::string> solution_names;
  switch (dim)
    {
      case 1:
	    solution_names.push_back ("displacement");
	    break;
      case 2:
	    solution_names.push_back ("x_displacement");	    
	    solution_names.push_back ("y_displacement");
	    break;
      case 3:
	    solution_names.push_back ("x_displacement");	    
	    solution_names.push_back ("y_displacement");
	    solution_names.push_back ("z_displacement");
	    break;
      default:
	    Assert (false, ExcInternalError());
    };
	     
  data_out.add_data_vector (solution, solution_names);
  data_out.build_patches ();
  data_out.write_gmv (output);
}



template <int dim>
void ElasticProblem<dim>::run () 
{
  for (unsigned int cycle=0; cycle<12; ++cycle)
    {
      std::cout << "Cycle " << cycle << ':' << std::endl;

      if (cycle == 0)
	{
	  GridGenerator::hyper_cube (triangulation, -1, 1);
	  triangulation.refine_global (2);

                                           // xxx
          GridTools::partition_triangulation (n_partitions, triangulation);
	}
      else
	refine_grid ();

      std::cout << "   Number of active cells:       "
		<< triangulation.n_active_cells()
		<< std::endl;

      setup_system ();

      std::cout << "   Number of degrees of freedom: "
		<< dof_handler.n_dofs()
		<< std::endl;
      
      assemble_system ();
      solve ();
      output_results (cycle);
    };
}


int main (int argc, char **argv) 
{
  try
    {
                                       // xxx
      PetscInitialize(&argc,&argv,0,0);
      deallog.depth_console (0);

                                       // xxx localize scope
      {
        ElasticProblem<2> elastic_problem_2d;
        elastic_problem_2d.run ();
      }

                                       // xxx
      PetscFinalize();      
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      
      return 1;
    }
  catch (...) 
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    };

  return 0;
}
