/* $Id$ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 2000 */

/*    $Id$       */
/*                                                                */
/*    Copyright (C) 2011 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vectors.h>
#include <deal.II/numerics/matrices.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>

#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_out.h>


#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/grid/grid_refinement.h>

#include <deal.II/numerics/error_estimator.h>

using namespace dealii;



template <int dim>
class LaplaceProblem
{
  public:
    LaplaceProblem ();
    ~LaplaceProblem ();

    void run ();

  private:
    bool interface_intersects_cell (const typename Triangulation<dim>::cell_iterator &cell) const;

    void setup_system ();
    void assemble_system ();
    void solve ();
    void refine_grid ();
    void output_results (const unsigned int cycle) const;

    Triangulation<dim>    triangulation;

    hp::DoFHandler<dim>   dof_handler;
    hp::FECollection<dim> fe_collection;

    ConstraintMatrix      constraints;

    SparsityPattern       sparsity_pattern;
    SparseMatrix<double>  system_matrix;

    Vector<double>        solution;
    Vector<double>        system_rhs;
};




template <int dim>
class Coefficient : public Function<dim>
{
  public:
    Coefficient () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;

    virtual void value_list (const std::vector<Point<dim> > &points,
			     std::vector<double>            &values,
			     const unsigned int              component = 0) const;
};



template <int dim>
double Coefficient<dim>::value (const Point<dim> &p,
				const unsigned int) const
{
  if (p.square() < 0.5*0.5)
    return 20;
  else
    return 1;
}



template <int dim>
void Coefficient<dim>::value_list (const std::vector<Point<dim> > &points,
				   std::vector<double>            &values,
				   const unsigned int              component) const
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
    }
}




template <int dim>
LaplaceProblem<dim>::LaplaceProblem ()
		:
		dof_handler (triangulation)
{
  fe_collection.push_back (FESystem<dim> (FE_Q<dim>(1), 1,
					  FE_Nothing<dim>(), 1));
  fe_collection.push_back (FESystem<dim> (FE_Q<dim>(1), 1,
					  FE_Q<dim>(1), 1));
}



template <int dim>
LaplaceProblem<dim>::~LaplaceProblem ()
{
  dof_handler.clear ();
}



template <int dim>
bool
LaplaceProblem<dim>::
interface_intersects_cell (const typename Triangulation<dim>::cell_iterator &cell) const
{
  return false;
}



template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
				   // decide which element to use
				   // where. to do so, we need to know
				   // which elements are intersected
				   // by the interface, or are at
				   // least adjacent. to this end, in
				   // a first step, loop over all
				   // cells and record which vertices
				   // are on cells that are
				   // intersected
  std::vector<bool> vertex_is_on_intersected_cell (triangulation.n_vertices(),
						   false);
  for (typename Triangulation<dim>::cell_iterator cell
	 = triangulation.begin_active();
       cell != triangulation.end(); ++cell)
    if (interface_intersects_cell(cell))
      for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
	vertex_is_on_intersected_cell[cell->vertex_index(v)] = true;

				   // now loop over all cells
				   // again. if one of the vertices of
				   // a given cell is part of a cell
				   // that is intersected, then we
				   // need to use the enriched space
				   // there. otherwise, use the normal
				   // space
  for (typename hp::DoFHandler<dim>::cell_iterator cell
	 = dof_handler.begin_active();
       cell != dof_handler.end(); ++cell)
    {
      bool use_enriched_space = false;
      for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
	if (vertex_is_on_intersected_cell[cell->vertex_index(v)] == true)
	  {
	    use_enriched_space = true;
	    break;
	  }

      if (use_enriched_space == false)
	cell->set_active_fe_index(0);
      else
	cell->set_active_fe_index(1);
    }

  dof_handler.distribute_dofs (fe_collection);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());


  constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
					   constraints);
  constraints.close ();

				   // now constrain those enriched
				   // DoFs that are on cells that are
				   // not intersected but that are
				   // adjacent to cells that are
  for (typename hp::DoFHandler<dim>::cell_iterator cell
	 = dof_handler.begin_active();
       cell != dof_handler.end(); ++cell)
    if ((cell->active_fe_index() == 1)
	&&
	(interface_intersects_cell(cell) == false))
				       // we are on an enriched cell
				       // but it isn't
				       // intersected. see which
				       // vertices are not part of
				       // intersected cells and
				       // constrain these DoFs
      for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
	if (vertex_is_on_intersected_cell[cell->vertex_index(v)] == false)
	  constraints.add_line (cell->vertex_dof_index(v,1));
  constraints.close();


  CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);

  constraints.condense (c_sparsity);

  sparsity_pattern.copy_from(c_sparsity);

  system_matrix.reinit (sparsity_pattern);
}


template <int dim>
void LaplaceProblem<dim>::assemble_system ()
{
  const QGauss<dim>  quadrature_formula(3);

  FEValues<dim> plain_fe_values (fe_collection[0], quadrature_formula,
				 update_values    |  update_gradients |
				 update_quadrature_points  |  update_JxW_values);

  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix;
  Vector<double>       cell_rhs;

  std::vector<unsigned int> local_dof_indices;

  const Coefficient<dim> coefficient;
  std::vector<double>    coefficient_values (n_q_points);

  typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
      cell_matrix.reinit (dofs_per_cell, dofs_per_cell);
      cell_rhs.reinit (dofs_per_cell);

      cell_matrix = 0;
      cell_rhs = 0;

      plain_fe_values.reinit (cell);

      coefficient.value_list (plain_fe_values.get_quadrature_points(),
			      coefficient_values);

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  {
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      cell_matrix(i,j) += (coefficient_values[q_point] *
				   plain_fe_values.shape_grad(i,q_point) *
				   plain_fe_values.shape_grad(j,q_point) *
				   plain_fe_values.JxW(q_point));

	    cell_rhs(i) += (plain_fe_values.shape_value(i,q_point) *
			    1.0 *
			    plain_fe_values.JxW(q_point));
	  }

      local_dof_indices.resize (dofs_per_cell);
      cell->get_dof_indices (local_dof_indices);
      constraints.distribute_local_to_global (cell_matrix, cell_rhs,
					      local_dof_indices,
					      system_matrix, system_rhs);
    }

  std::map<unsigned int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
					    0,
					    ZeroFunction<dim>(2),
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
  SolverCG<>              solver (solver_control);

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  solver.solve (system_matrix, solution, system_rhs,
		preconditioner);

  constraints.distribute (solution);
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

  std::string filename = "solution-";
  filename += ('0' + cycle);
  filename += ".vtk";

  std::ofstream output (filename.c_str());

  DataOut<dim,hp::DoFHandler<dim> > data_out;

  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");
  data_out.build_patches ();

  data_out.write_vtk (output);
}




template <int dim>
void LaplaceProblem<dim>::run ()
{
  for (unsigned int cycle=0; cycle<8; ++cycle)
    {
      std::cout << "Cycle " << cycle << ':' << std::endl;

      if (cycle == 0)
	{
	  GridGenerator::hyper_ball (triangulation);

	  static const HyperBallBoundary<dim> boundary;
	  triangulation.set_boundary (0, boundary);

	  triangulation.refine_global (1);
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
