/* $Id$ */
/* Author: Wolfgang Bangerth, Texas A&M University, 2011 */

/*    $Id$       */
/*                                                                */
/*    Copyright (C) 2011 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */



#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <base/function.h>
#include <base/utilities.h>

#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/sparse_direct.h>
#include <lac/constraint_matrix.h>

#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_tools.h>
#include <grid/grid_refinement.h>

#include <dofs/dof_handler.h>
#include <dofs/dof_renumbering.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>

#include <fe/fe_q.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>
#include <fe/mapping_q1.h>

#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <numerics/error_estimator.h>

#include <fstream>
#include <sstream>

using namespace dealii;


template <int dim>
class StokesProblem
{
  public:
    StokesProblem (const unsigned int degree);
    void run ();

  private:
    void setup_dofs ();
    void assemble_system ();
    void solve ();
    void output_results (const unsigned int refinement_cycle) const;
    void refine_mesh ();

    const unsigned int   degree;

    Triangulation<dim>   triangulation;
    FESystem<dim>        fe;
    DoFHandler<dim>      dof_handler;

    ConstraintMatrix     constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> system_rhs;
};



template <int dim>
class BoundaryValues : public Function<dim>
{
  public:
    BoundaryValues () : Function<dim>(dim+1) {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &value) const;
};


template <int dim>
double
BoundaryValues<dim>::value (const Point<dim>  &p,
			    const unsigned int component) const
{
  Assert (component < this->n_components,
	  ExcIndexRange (component, 0, this->n_components));

  if (component == 0)
    return (p[0] < 0 ? -1 : (p[0] > 0 ? 1 : 0));
  return 0;
}


template <int dim>
void
BoundaryValues<dim>::vector_value (const Point<dim> &p,
				   Vector<double>   &values) const
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values(c) = BoundaryValues<dim>::value (p, c);
}



template <int dim>
class RightHandSide : public Function<dim>
{
  public:
    RightHandSide () : Function<dim>(dim+1) {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &value) const;

};


template <int dim>
double
RightHandSide<dim>::value (const Point<dim>  &/*p*/,
                           const unsigned int /*component*/) const
{
  return 0;
}


template <int dim>
void
RightHandSide<dim>::vector_value (const Point<dim> &p,
                                  Vector<double>   &values) const
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values(c) = RightHandSide<dim>::value (p, c);
}







template <int dim>
StokesProblem<dim>::StokesProblem (const unsigned int degree)
                :
                degree (degree),
                triangulation (Triangulation<dim>::maximum_smoothing),
                fe (FE_Q<dim>(degree+1), dim,
                    FE_Q<dim>(degree), 1),
                dof_handler (triangulation)
{}




template <int dim>
void StokesProblem<dim>::setup_dofs ()
{
  system_matrix.clear ();

  dof_handler.distribute_dofs (fe);
  DoFRenumbering::Cuthill_McKee (dof_handler);

  std::vector<unsigned int> block_component (dim+1,0);
  block_component[dim] = 1;
  DoFRenumbering::component_wise (dof_handler, block_component);

  {
    constraints.clear ();
    std::vector<bool> component_mask (dim+1, true);
    component_mask[dim] = false;
    DoFTools::make_hanging_node_constraints (dof_handler,
					     constraints);
    VectorTools::interpolate_boundary_values (dof_handler,
					      1,
					      BoundaryValues<dim>(),
					      constraints,
					      component_mask);
  }

  constraints.close ();

  std::vector<unsigned int> dofs_per_block (2);
  DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, block_component);
  const unsigned int n_u = dofs_per_block[0],
                     n_p = dofs_per_block[1];

  std::cout << "   Number of active cells: "
            << triangulation.n_active_cells()
            << std::endl
            << "   Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << " (" << n_u << '+' << n_p << ')'
            << std::endl;

  {
    CompressedSimpleSparsityPattern csp (dof_handler.n_dofs(),
					 dof_handler.n_dofs());

    DoFTools::make_sparsity_pattern (dof_handler, csp, constraints, false);
    sparsity_pattern.copy_from (csp);
  }

  system_matrix.reinit (sparsity_pattern);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
}



template <int dim>
void StokesProblem<dim>::assemble_system ()
{
  system_matrix=0;
  system_rhs=0;

  QGauss<dim>   quadrature_formula(degree+2);

  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values    |
                           update_quadrature_points  |
                           update_JxW_values |
                           update_gradients);

  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;

  const unsigned int   n_q_points      = quadrature_formula.size();

  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       local_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  const RightHandSide<dim>          right_hand_side;
  std::vector<Vector<double> >      rhs_values (n_q_points,
                                                Vector<double>(dim+1));

  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);

  std::vector<SymmetricTensor<2,dim> > phi_grads_u (dofs_per_cell);
  std::vector<double>                  div_phi_u   (dofs_per_cell);
  std::vector<double>                  phi_p       (dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      local_matrix = 0;
      local_rhs = 0;

      right_hand_side.vector_value_list(fe_values.get_quadrature_points(),
                                        rhs_values);

      for (unsigned int q=0; q<n_q_points; ++q)
	{
	  for (unsigned int k=0; k<dofs_per_cell; ++k)
	    {
	      phi_grads_u[k] = fe_values[velocities].symmetric_gradient (k, q);
	      div_phi_u[k]   = fe_values[velocities].divergence (k, q);
	      phi_p[k]       = fe_values[pressure].value (k, q);
	    }

	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    {
	      for (unsigned int j=0; j<=i; ++j)
		{
		  local_matrix(i,j) += (phi_grads_u[i] * phi_grads_u[j]
					- div_phi_u[i] * phi_p[j]
					- phi_p[i] * div_phi_u[j]
					+ phi_p[i] * phi_p[j])
				       * fe_values.JxW(q);

		}

	      const unsigned int component_i =
		fe.system_to_component_index(i).first;
	      local_rhs(i) += fe_values.shape_value(i,q) *
			      rhs_values[q](component_i) *
			      fe_values.JxW(q);
	    }
	}


      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=i+1; j<dofs_per_cell; ++j)
	  local_matrix(i,j) = local_matrix(j,i);

      cell->get_dof_indices (local_dof_indices);
      constraints.distribute_local_to_global (local_matrix, local_rhs,
					      local_dof_indices,
					      system_matrix, system_rhs);
    }
}




template <int dim>
void StokesProblem<dim>::solve ()
{
  SparseDirectUMFPACK direct_solver;
  direct_solver.initialize (system_matrix);
  direct_solver.vmult (solution, system_rhs);
}




template <int dim>
void
StokesProblem<dim>::output_results (const unsigned int refinement_cycle)  const
{
  std::vector<std::string> solution_names (dim, "velocity");
  solution_names.push_back ("pressure");

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation
    (dim, DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation
    .push_back (DataComponentInterpretation::component_is_scalar);

  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, solution_names,
			    DataOut<dim>::type_dof_data,
			    data_component_interpretation);
  data_out.build_patches ();

  std::ostringstream filename;
  filename << "solution-"
           << Utilities::int_to_string (refinement_cycle, 2)
           << ".vtk";

  std::ofstream output (filename.str().c_str());
  data_out.write_vtk (output);
}



template <int dim>
void
StokesProblem<dim>::refine_mesh ()
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

  std::vector<bool> component_mask (dim+1, false);
  component_mask[dim] = true;
  KellyErrorEstimator<dim>::estimate (dof_handler,
                                      QGauss<dim-1>(degree+1),
                                      typename FunctionMap<dim>::type(),
                                      solution,
                                      estimated_error_per_cell,
                                      component_mask);

  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                   estimated_error_per_cell,
                                                   0.3, 0.0);
  triangulation.execute_coarsening_and_refinement ();
}



template <int dim>
void StokesProblem<dim>::run ()
{
  {
    std::vector<unsigned int> subdivisions (dim, 1);
    subdivisions[0] = 4;

    const Point<dim> bottom_left = (dim == 2 ?
				    Point<dim>(-2,-1) :
				    Point<dim>(-2,0,-1));
    const Point<dim> top_right   = (dim == 2 ?
				    Point<dim>(2,0) :
				    Point<dim>(2,1,0));

    GridGenerator::subdivided_hyper_rectangle (triangulation,
					       subdivisions,
					       bottom_left,
					       top_right);
  }

  for (typename Triangulation<dim>::active_cell_iterator
	 cell = triangulation.begin_active();
       cell != triangulation.end(); ++cell)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      if (cell->face(f)->center()[dim-1] == 0)
	cell->face(f)->set_all_boundary_indicators(1);


  triangulation.refine_global (4-dim);

  for (unsigned int refinement_cycle = 0; refinement_cycle<6;
       ++refinement_cycle)
    {
      std::cout << "Refinement cycle " << refinement_cycle << std::endl;

      if (refinement_cycle > 0)
        refine_mesh ();

      setup_dofs ();

      std::cout << "   Assembling..." << std::endl << std::flush;
      assemble_system ();

      std::cout << "   Solving..." << std::flush;
      solve ();

      output_results (refinement_cycle);

      std::cout << std::endl;
    }
}



int main ()
{
  try
    {
      deallog.depth_console (0);

      StokesProblem<2> flow_problem(1);
      flow_problem.run ();
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
