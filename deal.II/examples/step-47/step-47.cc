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
    unsigned int compute_quadrature(const Quadrature<dim> &plain_quadrature, const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const std::vector<double> &level_set_values);
		void append_quadrature(Quadrature<dim> plain_quadrature, std::vector<Point<dim> > v);

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
  const QGauss<dim>  quadrature_formula(2);


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

	std::vector<double> level_set_values;
	level_set_values.push_back(-1./2.);
	level_set_values.push_back(1);
	level_set_values.push_back(-1./2.);
	level_set_values.push_back(1.);
  for (; cell!=endc; ++cell)
    {
			compute_quadrature(quadrature_formula, cell, level_set_values);
			assert(0);

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

// To integrate the enriched elements we have to find the geometrical decomposition
// of the original element in subelements. The subelements are used to integrate
// the elements on both sides of the discontinuity. The disontinuity line is approximated
// by a piece-wise linear interpolation between the intersection of the discontinuity
// with the edges of the elements. The vector level_set_values has the values of
// the level set function at the vertices of the elements. From these values can be found
// by linear interpolation the intersections. There are three kind of decomposition that
// are considered.
// Type 1: there is not cut. Type 2: a corner of the element is cut. Type 3: two corners are cut.

	template <int dim>
	std::pair<unsigned int, Quadrature<dim> >
	LaplaceProblem<dim>::compute_quadrature (const Quadrature<dim> &plain_quadrature,
		const typename hp::DoFHandler<dim>::active_cell_iterator &cell,
		const std::vector<double> &level_set_values                    )
{

	unsigned int type = 0;

	// find the type of cut
	int sign_ls[GeometryInfo<dim>::vertices_per_cell];
	for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
	{
		if (level_set_values[v] > 0) sign_ls[v] = 1;
		else if (level_set_values[v] < 0) sign_ls[v] = -1;
		else sign_ls[v] = 0;
	}

	// the sign of the level set function at the 4 nodes of the elements can be positive + or negative -
	// depending on the sign of the level set function we have the folloing three classes of decomposition
  // type 1: ++++, ----
	// type 2: -+++, +-++, ++-+, +++-, +---, -+--, --+-, ---+
	// type 3: +--+, ++--, +-+-, -++-, --++, -+-+

	if ( sign_ls[0]==sign_ls[1] & sign_ls[0]==sign_ls[2] & sign_ls[0]==sign_ls[3] ) type =1;
	else if ( sign_ls[0]*sign_ls[1]*sign_ls[2]*sign_ls[3] < 0 ) type = 2;
	else type = 3;

	unsigned int Pos = 100;
	Point<dim> A(0,0);
	Point<dim> B(0,0);
	Point<dim> C(0,0);
	Point<dim> D(0,0);
	Point<dim> E(0,0);
	Point<dim> F(0,0);

	if (type == 1)
	  return std::pair<unsigned int, Quadrature<dim> >(1, plain_quadrature);

	if (type==2)
	{
		const unsigned int   n_q_points    = plain_quadrature.size();

		// loop over all subelements for integration
		// in type 2 there are 5 subelements

		Quadrature<dim> xfem_quadrature(5*n_q_points);

		std::vector<Point<dim> > v(GeometryInfo<dim>::vertices_per_cell);

		if (sign_ls[0]!=sign_ls[1] && sign_ls[0]!=sign_ls[2] && sign_ls[0]!=sign_ls[3]) Pos = 0;
    else if (sign_ls[1]!=sign_ls[0] && sign_ls[1]!=sign_ls[2] && sign_ls[1]!=sign_ls[3]) Pos = 1;
    else if (sign_ls[2]!=sign_ls[0] && sign_ls[2]!=sign_ls[1] && sign_ls[2]!=sign_ls[3]) Pos = 2;
    else if (sign_ls[3]!=sign_ls[0] && sign_ls[3]!=sign_ls[1] && sign_ls[3]!=sign_ls[2]) Pos = 3;
		else assert(0); // error message

		std::cout << "Pos " << Pos << std::endl;

		// Find cut coordinates

		// deal.ii local coordinates

		//    2-------3
		//    |       |
		//		|       |
		//		|       |
		//		0-------1

		if (Pos == 0)
		{
			A[0] = 1. - level_set_values[1]/(level_set_values[1]-level_set_values[0]);
			B[1] = 1. - level_set_values[2]/(level_set_values[2]-level_set_values[0]);
			A(1) = 0.;
			B(0) = 0.;
			C(0) = 0.5*( A(0) + B(0) );
			C(1) = 0.5*( A(1) + B(1) );
			D(0) = 2./3. * C(0);
			D(1) = 2./3. * C(1);
			E(0) = 0.5*A(0);
			E(1) = 0.;
			F(0) = 0.;
			F(1) = 0.5*B(1);
		}
		else if (Pos == 1)
		{
			A[0] = level_set_values[0]/(level_set_values[0]-level_set_values[1]);
			B[1] = 1 - level_set_values[3]/(level_set_values[3]-level_set_values[1]);
			A(1) = 0.;
			B(0) = 1.;
			C(0) = 0.5*( A(0) + B(0) );
			C(1) = 0.5*( A(1) + B(1) );
			D(0) = 2./3. * C(0);
			D(1) = 2./3. * C(1);
			E(0) = 0.5*A(0);
			E(1) = 0.;
			F(0) = 1.;
			F(1) = 0.5*B(1);
			}
		else if (Pos == 2)
		{
			A[0] = 1 - level_set_values[3]/(level_set_values[3]-level_set_values[2]);
			B[1] = level_set_values[0]/(level_set_values[0]-level_set_values[2]);
			A(1) = 1.;
			B(0) = 0.;
			C(0) = 0.5*( A(0) + B(0) );
			C(1) = 0.5*( A(1) + B(1) );
			D(0) = 2./3. * C(0);
			D(1) = 1./3. + 2./3. * C(1);
			E(0) = 0.5* A(0);
			E(1) = 0.;
			F(0) = 0.;
			F(1) = 0.5*( 1. + B(1) );
			}
		else if (Pos == 3)
		{
			A[0] = level_set_values[2]/(level_set_values[2]-level_set_values[3]);
			B[1] = level_set_values[1]/(level_set_values[1]-level_set_values[3]);
			A(1) = 1.;
			B(0) = 1.;
			C(0) = 0.5*( A(0) + B(0) );
			C(1) = 0.5*( A(1) + B(1) );
			D(0) = 1./3. + 2./3. * C(0);
			D(1) = 1./3. + 2./3. * C(1);
			E(0) = 0.5*( 1. + A(0) );
			E(1) = 1.;
			F(0) = 1.;
			F(1) = 0.5*( 1. + B(1) );
			}

    Point<dim> v0(0,0);
    Point<dim> v1(1,0);
    Point<dim> v2(0,1);
    Point<dim> v3(1,1);

		std::cout << A << std::endl;
		std::cout << B << std::endl;
		std::cout << C << std::endl;
		std::cout << D << std::endl;
		std::cout << E << std::endl;
		std::cout << F << std::endl;

		Point<dim> subcell_vertices[10];
		subcell_vertices[0] = v0;
		subcell_vertices[1] = v1;
		subcell_vertices[2] = v2;
		subcell_vertices[3] = v3;
		subcell_vertices[4] = A;
		subcell_vertices[5] = B;
		subcell_vertices[6] = C;
		subcell_vertices[7] = D;
		subcell_vertices[8] = E;
		subcell_vertices[9] = F;

		// lookup table for the decomposition

		if (dim==2)
		{
			 unsigned int subcell_v_indices[4][5][4] = {
				 {{0,8,9,7}, {9,7,5,6}, {8,4,7,6}, {5,6,2,3}, {4,1,6,3}},
				 {{8,1,7,9}, {4,8,6,7}, {7,9,6,8}, {0,4,2,6}, {2,6,3,5}},
				 {{7,9,8,3}, {4,6,8,7}, {6,5,7,9}, {0,6,2,4}, {0,1,6,5}},
				 {{9,7,2,8}, {5,6,9,7}, {6,4,7,8}, {0,1,5,6}, {6,1,4,3}}
			 };

			 for (unsigned int subcell = 0; subcell<5; subcell++)
			 {
				 std::vector<Point<dim> > vertices;
				 for (unsigned int i=0; i<4; i++)
				 {
					 vertices.push_back( subcell_vertices[subcell_v_indices[Pos][subcell][i]] );
					 //std::cout << "Pos       : " << Pos << std::endl;
					 //std::cout << "subcell   : " << subcell << std::endl;
					 //std::cout << "i         : " << i << std::endl;
					 //std::cout << "subcell v : " << subcell_v_indices[Pos][subcell][i] << std::endl;
					 //std::cout << vertices[i](0) << "  " << vertices[i](1) << std::endl;
				 }
				 std::cout << std::endl;
				 // create quadrature rule
				 append_quadrature( xfem_quadrature,
							                vertices     );
			 }

		}

		return std::pair<unsigned int, Quadrature<dim> >(2, xfem_quadrature);
	}

	// Type three decomposition
	// (+--+, ++--, +-+-, -++-, --++, -+-+)

	if (type==3)
	{
		const unsigned int   n_q_points    = plain_quadrature.size();

		// loop over all subelements for integration
		// in type 2 there are 5 subelements

		Quadrature<dim> xfem_quadrature(5*n_q_points);

		std::vector<Point<dim> > v(GeometryInfo<dim>::vertices_per_cell);

		if ( sign_ls[0]==sign_ls[1] && sign_ls[2]==sign_ls[3] )
		{
			Pos = 0;
			A(1) = level_set_values[0]/((level_set_values[0]-level_set_values[2]));
			B(1) = level_set_values[1]/((level_set_values[1]-level_set_values[3]));
		}
		else if ( sign_ls[0]==sign_ls[2] && sign_ls[1]==sign_ls[3] )
		{
			Pos = 1;
			A(0) = level_set_values[0]/((level_set_values[0]-level_set_values[1]));
			B(0) = level_set_values[2]/((level_set_values[2]-level_set_values[3]));
		}
		else if ( sign_ls[0]==sign_ls[3] && sign_ls[1]==sign_ls[2] )
		{
			std::cout << "Error: the element has two cut lines and this is not allowed" << std::endl;
			assert(0);
		}
		else
		{
			std::cout << "Error: the level set function has not the right values" << std::endl;
			assert(0);
		}

		std::cout << "Pos " << Pos << std::endl;
		std::cout << A << std::endl;
		std::cout << B << std::endl;

//TODO: fill xfem_quadrature

		return std::pair<unsigned int, Quadrature<dim> >(3, xfem_quadrature);
	}



	return 0;

}

	template <int dim>
void LaplaceProblem<dim>::append_quadrature ( const Quadrature<dim> &plain_quadrature,
					      const std::vector<Point<dim> > &v,
					      std::vector<Point<dim> > &xfem_points,
					      std::vector<double>      &xfem_weights)

{
	// Project integration points into sub-elements.
	// Map F1.
	// The map F1 maps quadrature points from a reference element to a subelement of a reference element.
	// To implement the action of this map the coordinates of the subelements have been calculated (A(0)...F(0),A(1)...F(1))
	// the coordinates of the quadrature points are given by the bi-linear map defined by the form functions
	// $x^\prime_i = \sum_j v^\prime \phi_j(x^hat_i)$, where the $\phi_j$ are the shape functions of the FEQ.

	unsigned int n_v = GeometryInfo<dim>::vertices_per_cell;

	std::vector<Point<dim> > q_points = plain_quadrature.get_points();
	std::vector<Point<dim> > q_transf(q_points.size());
	std::vector<double> W = plain_quadrature.get_weights();
	std::vector<double> phi(n_v);
	std::vector<Tensor<1,dim> > grad_phi(n_v);

	const unsigned int   n_q_points    = plain_quadrature.size();

	std::vector<double> JxW(n_q_points);

	for ( unsigned int i = 0; i < n_q_points; i++)
	{
	  switch (dim)
	    {
	      case 2:
	      {
		double xi  = q_points[i](0);
		double eta = q_points[i](1);

						 // Define shape functions on reference element
						 // we consider a bi-linear mapping
		phi[0] = (1. - xi) * (1. - eta);
		phi[1] = xi * (1. - eta);
		phi[2] = (1. - xi) * eta;
		phi[3] = xi * eta;

		grad_phi[0][0] = (-1. + eta);
		grad_phi[1][0] = (1. - eta);
		grad_phi[2][0] = -eta;
		grad_phi[3][0] = eta;

		grad_phi[0][1] = (-1. + xi);
		grad_phi[1][1] = -xi;
		grad_phi[2][1] = -xi;
		grad_phi[3][1] = xi;

		break;
	      }

	      default:
		    Assert (false, ExcNotImplemented());
	    }


	  Tensor<2,dim> jacobian;

					   // Calculate Jacobian of transformation
	  for (unsigned int d=0; d<dim; ++d)
	    for (unsigned int e=0; e<dim; ++e)
	      for (unsigned int j = 0; j<GeometryInfo<dim>::vertices_per_cell; j++)
		jacobian[d][e] += grad_phi[j][d] * v[j](e);

/*

	  for (unsigned int j = 0; j<GeometryInfo<dim>::vertices_per_cell; j++)
	    {
	      dx_dxi  += dphi_dxi[j]  * v[j](0);
	      dx_deta += dphi_deta[j] * v[j](0);
	      dy_dxi  += dphi_dxi[j]  * v[j](1);
	      dy_deta += dphi_deta[j] * v[j](1);
	    }
*/
	  double detJ = determinant(jacobian);
	  xfem_weights.push_back (W[i] * detJ);

		// Map integration points from reference element to subcell of reference elemment
		Point<dim> q_prime;
		for (unsigned int d=0; d<dim; ++d)
		  for (unsigned int j = 0; j<GeometryInfo<dim>::vertices_per_cell; j++)
		    q_prime[d] += v[j](d) * phi[j];
		xfem_points.push_back(q_prime);
	}

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
