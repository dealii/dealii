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

#include <dofs/dof_renumbering.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_accessor.h>

#include <fe/fe_q.h>
#include <fe/fe_nothing.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>

#include <hp/dof_handler.h>
#include <hp/fe_collection.h>
#include <hp/fe_values.h>

#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <numerics/error_estimator.h>

#include <fstream>
#include <sstream>

using namespace dealii;


template <int dim>
class FluidStructureProblem
{
  public:
    FluidStructureProblem (const unsigned int stokes_degree,
			   const unsigned int elasticity_degree);
    void run ();

  private:
    void setup_dofs ();
    void assemble_system ();
    void solve ();
    void output_results (const unsigned int refinement_cycle) const;
    void refine_mesh ();

    const unsigned int    stokes_degree;
    const unsigned int    elasticity_degree;

    Triangulation<dim>    triangulation;
    FESystem<dim>         stokes_fe;
    FESystem<dim>         elasticity_fe;
    hp::FECollection<dim> fe_collection;
    hp::DoFHandler<dim>   dof_handler;

    ConstraintMatrix      constraints;

    SparsityPattern       sparsity_pattern;
    SparseMatrix<double>  system_matrix;

    Vector<double> solution;
    Vector<double> system_rhs;

    const double lambda;
    const double mu;
};



template <int dim>
class BoundaryValues : public Function<dim>
{
  public:
    BoundaryValues () : Function<dim>(dim+1+dim) {}

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

  if (component == 1)
    return std::sin(numbers::PI*p[0]);

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
FluidStructureProblem<dim>::FluidStructureProblem (const unsigned int stokes_degree,
				   const unsigned int elasticity_degree)
                :
                stokes_degree (stokes_degree),
		elasticity_degree (elasticity_degree),
                triangulation (Triangulation<dim>::maximum_smoothing),
                stokes_fe (FE_Q<dim>(stokes_degree+1), dim,
			   FE_Q<dim>(stokes_degree), 1,
			   FE_Nothing<dim>(), dim),
                elasticity_fe (FE_Nothing<dim>(), dim,
			       FE_Nothing<dim>(), 1,
			       FE_Q<dim>(elasticity_degree), dim),
                dof_handler (triangulation),
		lambda (1),
		mu (1)
{
  fe_collection.push_back (stokes_fe);
  fe_collection.push_back (elasticity_fe);
}




template <int dim>
void FluidStructureProblem<dim>::setup_dofs ()
{
  system_matrix.clear ();

  for (typename hp::DoFHandler<dim>::active_cell_iterator
	 cell = dof_handler.begin_active();
       cell != dof_handler.end(); ++cell)
    if (((std::fabs(cell->center()[0]) < 0.25)
	 &&
	 (cell->center()[1] > 0.5))
	||
	((std::fabs(cell->center()[0]) >= 0.25)
	 &&
	 (cell->center()[1] > -0.5)))
      cell->set_active_fe_index (0);
    else
      cell->set_active_fe_index (1);

  dof_handler.distribute_dofs (fe_collection);
  DoFRenumbering::Cuthill_McKee (dof_handler);

  {
    constraints.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler,
					     constraints);

    std::vector<bool> velocity_mask (dim+1+dim, false);
    for (unsigned int d=0; d<dim; ++d)
      velocity_mask[d] = true;
    VectorTools::interpolate_boundary_values (dof_handler,
					      1,
					      BoundaryValues<dim>(),
					      constraints,
					      velocity_mask);
    std::vector<bool> elasticity_mask (dim+1+dim, false);
    for (unsigned int d=dim+1; d<dim+1+dim; ++d)
      elasticity_mask[d] = true;
    VectorTools::interpolate_boundary_values (dof_handler,
					      0,
					      ZeroFunction<dim>(dim+1+dim),
					      constraints,
					      elasticity_mask);
  }

				   // make sure velocity is zero at
				   // the interface
  {
    std::vector<unsigned int> local_face_dof_indices (stokes_fe.dofs_per_face);
    for (typename hp::DoFHandler<dim>::active_cell_iterator
	   cell = dof_handler.begin_active();
	 cell != dof_handler.end(); ++cell)
      if (cell->active_fe_index() == 0)
	for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
	  if (not cell->at_boundary(f))
	    {
	      bool face_is_on_interface = false;

	      if ((cell->neighbor(f)->has_children() == false)
		  &&
		  (cell->neighbor(f)->active_fe_index() == 1))
		face_is_on_interface = true;
	      else if (cell->neighbor(f)->has_children() == true)
		{
						   // neighbor does
						   // have
						   // children. see if
						   // any of the cells
						   // on the other
						   // side are elastic
		  for (unsigned int sf=0; sf<cell->face(f)->n_children(); ++sf)
		    if (cell->neighbor_child_on_subface(f, sf)->active_fe_index() == 1)
		      {
			face_is_on_interface = true;
			break;
		      }
		}

	      if (face_is_on_interface)
		{
		  cell->face(f)->get_dof_indices (local_face_dof_indices, 0);
		  for (unsigned int i=0; i<local_face_dof_indices.size(); ++i)
		    constraints.add_line (local_face_dof_indices[i]);
		}
	    }
  }


  constraints.close ();

  std::cout << "   Number of active cells: "
            << triangulation.n_active_cells()
            << std::endl
            << "   Number of degrees of freedom: "
            << dof_handler.n_dofs()
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
void FluidStructureProblem<dim>::assemble_system ()
{
  system_matrix=0;
  system_rhs=0;

  const double viscosity = 2;

  const QGauss<dim> stokes_quadrature(stokes_degree+2);
  const QGauss<dim> elasticity_quadrature(elasticity_degree+2);

  hp::QCollection<dim>  q_collection;
  q_collection.push_back (stokes_quadrature);
  q_collection.push_back (elasticity_quadrature);

  hp::FEValues<dim> hp_fe_values (fe_collection, q_collection,
				  update_values    |
				  update_quadrature_points  |
				  update_JxW_values |
				  update_gradients);

  const QGauss<dim-1> face_quadrature(std::max (stokes_degree+2,
						elasticity_degree+2));

  FEFaceValues<dim> stokes_fe_face_values (stokes_fe, face_quadrature,
					   update_JxW_values |
					   update_normal_vectors |
					   update_gradients);
  FEFaceValues<dim> elasticity_fe_face_values (elasticity_fe, face_quadrature,
					       update_values);

  const unsigned int   stokes_dofs_per_cell     = stokes_fe.dofs_per_cell;
  const unsigned int   elasticity_dofs_per_cell = elasticity_fe.dofs_per_cell;

  FullMatrix<double>   local_matrix;
  Vector<double>       local_rhs;

  std::vector<unsigned int> local_dof_indices;

  const RightHandSide<dim>          right_hand_side;

  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);
  const FEValuesExtractors::Vector displacements (dim+1);

  std::vector<SymmetricTensor<2,dim> > stokes_phi_grads_u (stokes_dofs_per_cell);
  std::vector<double>                  stokes_div_phi_u   (stokes_dofs_per_cell);
  std::vector<double>                  stokes_phi_p       (stokes_dofs_per_cell);

  std::vector<Tensor<2,dim> >          elasticity_phi_grad (elasticity_dofs_per_cell);
  std::vector<double>                  elasticity_phi_div  (elasticity_dofs_per_cell);
  std::vector<Tensor<1,dim> >          elasticity_phi  (elasticity_dofs_per_cell);

  typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      hp_fe_values.reinit (cell);

      const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

      local_matrix.reinit (cell->get_fe().dofs_per_cell,
			   cell->get_fe().dofs_per_cell);
      local_rhs.reinit (cell->get_fe().dofs_per_cell);

      if (cell->active_fe_index() == 0)
	{
	  const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
	  Assert (dofs_per_cell == stokes_dofs_per_cell,
		  ExcInternalError());

	  for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q)
	    {
	      for (unsigned int k=0; k<dofs_per_cell; ++k)
		{
		  stokes_phi_grads_u[k] = fe_values[velocities].symmetric_gradient (k, q);
		  stokes_div_phi_u[k]   = fe_values[velocities].divergence (k, q);
		  stokes_phi_p[k]       = fe_values[pressure].value (k, q);
		}

	      for (unsigned int i=0; i<dofs_per_cell; ++i)
		for (unsigned int j=0; j<dofs_per_cell; ++j)
		  local_matrix(i,j) += (2 * viscosity * stokes_phi_grads_u[i] * stokes_phi_grads_u[j]
					- stokes_div_phi_u[i] * stokes_phi_p[j]
					- stokes_phi_p[i] * stokes_div_phi_u[j])
				       * fe_values.JxW(q);
	    }
	}
      else
	{
	  const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
	  Assert (dofs_per_cell == elasticity_dofs_per_cell,
		  ExcInternalError());

	  for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q)
	    {
	      for (unsigned int k=0; k<dofs_per_cell; ++k)
		{
		  elasticity_phi_grad[k] = fe_values[displacements].gradient (k, q);
		  elasticity_phi_div[k]  = fe_values[displacements].divergence (k, q);
		}

	    for (unsigned int i=0; i<dofs_per_cell; ++i)
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
		{
		  local_matrix(i,j)
		    +=  (lambda *
			 elasticity_phi_div[i] * elasticity_phi_div[j]
			 +
			 mu *
			 scalar_product(elasticity_phi_grad[i], elasticity_phi_grad[j])
			 +
			 mu *
			 scalar_product(elasticity_phi_grad[i], transpose(elasticity_phi_grad[j]))
		      )
		      *
		      fe_values.JxW(q);
		  }
	      }
	}

      const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
      local_dof_indices.resize (dofs_per_cell);
      cell->get_dof_indices (local_dof_indices);

				       // local_rhs==0, but need to do
				       // this here because of
				       // boundary values
      constraints.distribute_local_to_global (local_matrix, local_rhs,
					      local_dof_indices,
					      system_matrix, system_rhs);

				       // see about face terms
      if (cell->active_fe_index() == 0)
					 // we are on a fluid cell
	for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
	  if (cell->at_boundary(f) == false)
	    {
	      if ((cell->neighbor(f)->level() == cell->level())
		  &&
		  (cell->neighbor(f)->has_children() == false))
		{
						   // cells are on same level
		  if (cell->neighbor(f)->active_fe_index() == 0)
						     // neighbor is also a fluid cell
		    continue;

						   // same size
						   // neighbors;
						   // neighbor is
						   // solids cell
		  stokes_fe_face_values.reinit (cell, f);
		  elasticity_fe_face_values.reinit (cell->neighbor(f),
						    cell->neighbor_of_neighbor(f));

		  local_matrix.reinit (stokes_dofs_per_cell,
				       elasticity_dofs_per_cell);
		  for (unsigned int q=0; q<face_quadrature.size(); ++q)
		    {
		      const Tensor<1,dim> normal_vector = stokes_fe_face_values.normal_vector(q);

		      for (unsigned int k=0; k<stokes_dofs_per_cell; ++k)
			stokes_phi_grads_u[k] = stokes_fe_face_values[velocities].symmetric_gradient (k, q);
		      for (unsigned int k=0; k<elasticity_dofs_per_cell; ++k)
			elasticity_phi[k] = elasticity_fe_face_values[displacements].value (k,q);

		      for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
			for (unsigned int j=0; j<elasticity_dofs_per_cell; ++j)
			  local_matrix(i,j) += (2 * viscosity *
						(stokes_phi_grads_u[i] *
						 normal_vector) *
						elasticity_phi[j] *
						stokes_fe_face_values.JxW(q));
		    }
		}
	      else if ((cell->neighbor(f)->level() == cell->level())
		       &&
		       (cell->neighbor(f)->has_children() == true))
		{
						   // neighbor has children
		}
	      else
		{
						   // neighbor is coarser
		  if (cell->active_fe_index() == cell->neighbor(f)->active_fe_index())
		    continue;
		}
	    }
    }
}




template <int dim>
void FluidStructureProblem<dim>::solve ()
{
  SparseDirectUMFPACK direct_solver;
  direct_solver.initialize (system_matrix);
  direct_solver.vmult (solution, system_rhs);

  constraints.distribute (solution);
}




template <int dim>
void
FluidStructureProblem<dim>::output_results (const unsigned int refinement_cycle)  const
{
  std::vector<std::string> solution_names (dim, "velocity");
  solution_names.push_back ("pressure");
  for (unsigned int d=0; d<dim; ++d)
    solution_names.push_back ("displacement");

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation
    (dim, DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation
    .push_back (DataComponentInterpretation::component_is_scalar);
  for (unsigned int d=0; d<dim; ++d)
    data_component_interpretation
      .push_back (DataComponentInterpretation::component_is_part_of_vector);

  DataOut<dim,hp::DoFHandler<dim> > data_out;
  data_out.attach_dof_handler (dof_handler);

  data_out.add_data_vector (solution, solution_names,
			    DataOut<dim,hp::DoFHandler<dim> >::type_dof_data,
			    data_component_interpretation);
  data_out.build_patches (stokes_fe.degree);

  std::ostringstream filename;
  filename << "solution-"
           << Utilities::int_to_string (refinement_cycle, 2)
           << ".vtk";

  std::ofstream output (filename.str().c_str());
  data_out.write_vtk (output);
}



template <int dim>
void
FluidStructureProblem<dim>::refine_mesh ()
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

  const QGauss<dim-1> stokes_quadrature(stokes_degree+2);
  const QGauss<dim-1> elasticity_quadrature(elasticity_degree+2);

  hp::QCollection<dim-1> face_q_collection;
  face_q_collection.push_back (stokes_quadrature);
  face_q_collection.push_back (elasticity_quadrature);

  std::vector<bool> component_mask (dim+1+dim, false);
  component_mask[dim] = true;
  KellyErrorEstimator<dim>::estimate (dof_handler,
                                      face_q_collection,
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
void FluidStructureProblem<dim>::run ()
{
  GridGenerator::hyper_cube (triangulation, -1, 1);
  for (typename Triangulation<dim>::active_cell_iterator
	 cell = triangulation.begin_active();
       cell != triangulation.end(); ++cell)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      if (cell->face(f)->at_boundary()
	  &&
	  (cell->face(f)->center()[dim-1] == 1))
	cell->face(f)->set_all_boundary_indicators(1);
  triangulation.refine_global (3);

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

      FluidStructureProblem<2> flow_problem(1, 1);
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
