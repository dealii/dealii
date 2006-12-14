/* $Id$ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 2000 */

/*    $Id$       */
/*    Version: $Name$                                          */
/*                                                                */
/*    Copyright (C) 2000, 2001, 2002, 2003, 2004, 2006 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

                                 // @sect3{Include files}

				 // The first few files have already
				 // been covered in previous examples
				 // and will thus not be further
				 // commented on.
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/logstream.h>
#include <base/utilities.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <dofs/hp_dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/hp_fe_values.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>

#include <fstream>
#include <iostream>

#include <fe/fe_q.h>
#include <grid/grid_out.h>
#include <dofs/dof_constraints.h>
#include <grid/grid_refinement.h>
#include <numerics/error_estimator.h>

				 // Finally, this is as in previous
				 // programs:
using namespace dealii;

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

    hp::DoFHandler<dim>      dof_handler;
    hp::FECollection<dim>    fe_collection;
    hp::QCollection<dim>     quadrature_collection;

    ConstraintMatrix     hanging_node_constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       system_rhs;
};


template <int dim>
LaplaceProblem<dim>::LaplaceProblem () :
		dof_handler (triangulation)
{
  for (unsigned int degree=1; degree<5; ++degree)
    {
      fe_collection.push_back (FE_Q<dim>(degree));
      quadrature_collection.push_back (QGauss<dim>(degree+2));
    }
}


template <int dim>
LaplaceProblem<dim>::~LaplaceProblem () 
{
  dof_handler.clear ();
}

template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe_collection);

  sparsity_pattern.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());

  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
					   hanging_node_constraints);

  hanging_node_constraints.close ();

  hanging_node_constraints.condense (sparsity_pattern);
  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);
}

template <int dim>
void LaplaceProblem<dim>::assemble_system () 
{
  hp::FEValues<dim> hp_fe_values (fe_collection,
				  quadrature_collection, 
				  update_values    |  update_gradients |
				  update_q_points  |  update_JxW_values);

  typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      const unsigned int   dofs_per_cell = cell->get_fe().dofs_per_cell;
      FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
      Vector<double>       cell_rhs (dofs_per_cell);

      std::vector<unsigned int> local_dof_indices (dofs_per_cell);

      cell_matrix = 0;
      cell_rhs = 0;

      hp_fe_values.reinit (cell);

      const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
      
      for (unsigned int q_point=0; q_point<fe_values.n_quadrature_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  {
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      cell_matrix(i,j) += (fe_values.shape_grad(i,q_point) *
				   fe_values.shape_grad(j,q_point) *
				   fe_values.JxW(q_point));

	    cell_rhs(i) += (fe_values.shape_value(i,q_point) *
			    1.0 *
			    fe_values.JxW(q_point));
	  }

      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    system_matrix.add (local_dof_indices[i],
			       local_dof_indices[j],
			       cell_matrix(i,j));
	  
	  system_rhs(local_dof_indices[i]) += cell_rhs(i);
	}
    }

  hanging_node_constraints.condense (system_matrix);
  hanging_node_constraints.condense (system_rhs);
  std::map<unsigned int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
					    0,
					    ZeroFunction<dim>(),
					    boundary_values);
  MatrixTools::apply_boundary_values (boundary_values,
				      system_matrix,
				      solution,
				      system_rhs);
}

template <int dim>
void LaplaceProblem<dim>::solve () 
{
  SolverControl           solver_control (1000, 1e-12);
  SolverCG<>              cg (solver_control);

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  cg.solve (system_matrix, solution, system_rhs,
	    preconditioner);

  hanging_node_constraints.distribute (solution);
}

template <int dim>
void LaplaceProblem<dim>::refine_grid ()
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

  KellyErrorEstimator<dim>::estimate (dof_handler,
				      QGauss<dim-1>(3),
				      typename FunctionMap<dim>::type(),
				      solution,
				      estimated_error_per_cell);

  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
						   estimated_error_per_cell,
						   0.3, 0.03);

  triangulation.execute_coarsening_and_refinement ();
}

template <int dim>
void LaplaceProblem<dim>::output_results (const unsigned int cycle) const
{
  Assert (cycle < 10, ExcNotImplemented());
  
  {
    const std::string filename = "grid-" +
				 Utilities::int_to_string (cycle, 2) +
				 ".eps";
    std::ofstream output (filename.c_str());
    
    GridOut grid_out;
    grid_out.write_eps (triangulation, output);
  }
  
  {
    const std::string filename = "solution-" +
				 Utilities::int_to_string (cycle, 2) +
				 ".gnuplot";
    DataOut<dim,hp::DoFHandler<dim> > data_out;

    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, "solution");
    data_out.build_patches ();
  
    std::ofstream output (filename.c_str());
    data_out.write_gnuplot (output);
  }
}


void
create_coarse_grid (Triangulation<2> &coarse_grid)
{
  const unsigned int dim = 2;
  static const Point<2> vertices_1[]
    = {  Point<2> (-1.,   -1.),
         Point<2> (-1./2, -1.),
         Point<2> (0.,    -1.),
         Point<2> (+1./2, -1.),
         Point<2> (+1,    -1.),
	     
         Point<2> (-1.,   -1./2.),
         Point<2> (-1./2, -1./2.),
         Point<2> (0.,    -1./2.),
         Point<2> (+1./2, -1./2.),
         Point<2> (+1,    -1./2.),
	     
         Point<2> (-1.,   0.),
         Point<2> (-1./2, 0.),
         Point<2> (+1./2, 0.),
         Point<2> (+1,    0.),
	     
         Point<2> (-1.,   1./2.),
         Point<2> (-1./2, 1./2.),
         Point<2> (0.,    1./2.),
         Point<2> (+1./2, 1./2.),
         Point<2> (+1,    1./2.),
	     
         Point<2> (-1.,   1.),
         Point<2> (-1./2, 1.),
         Point<2> (0.,    1.),			  
         Point<2> (+1./2, 1.),
         Point<2> (+1,    1.)    };
  const unsigned int
    n_vertices = sizeof(vertices_1) / sizeof(vertices_1[0]);
  const std::vector<Point<dim> > vertices (&vertices_1[0],
                                           &vertices_1[n_vertices]);
  static const int cell_vertices[][GeometryInfo<dim>::vertices_per_cell]
    = {{0, 1, 5, 6},
       {1, 2, 6, 7},
       {2, 3, 7, 8},
       {3, 4, 8, 9},
       {5, 6, 10, 11},
       {8, 9, 12, 13},
       {10, 11, 14, 15},
       {12, 13, 17, 18},
       {14, 15, 19, 20},
       {15, 16, 20, 21},
       {16, 17, 21, 22},
       {17, 18, 22, 23}};
  const unsigned int
    n_cells = sizeof(cell_vertices) / sizeof(cell_vertices[0]);

  std::vector<CellData<dim> > cells (n_cells, CellData<dim>());
  for (unsigned int i=0; i<n_cells; ++i) 
    {
      for (unsigned int j=0;
           j<GeometryInfo<dim>::vertices_per_cell;
           ++j)
        cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    }

  coarse_grid.create_triangulation (vertices,
                                    cells,
                                    SubCellData());
  coarse_grid.refine_global (1);
}



template <int dim>
void LaplaceProblem<dim>::run () 
{
  for (unsigned int cycle=0; cycle<5; ++cycle)
    {
      std::cout << "Cycle " << cycle << ':' << std::endl;

      if (cycle == 0)
	create_coarse_grid (triangulation);
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
    }
}

int main () 
{
  try
    {
      deallog.depth_console (0);

      LaplaceProblem<2> laplace_problem_2d;
      laplace_problem_2d.run ();
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
    }

  return 0;
}
