/*      $Id$                 */


#include <fe/fe.h>
#include <fe/fe_values.h>
#include <fe/fe_update_flags.h>
#include <fe/quadrature.h>
#include <numerics/error_estimator.h>
#include <grid/dof.h>
#include <grid/tria_iterator.h>
#include <grid/dof_accessor.h>
#include <lac/dvector.h>

#include <numeric>
#include <algorithm>



inline double sqr (const double x) {
  return x*x;
};



void KellyErrorEstimator<1>::estimate_error (const DoFHandler<1> &,
					     const Quadrature<0> &,
					     const FiniteElement<1> &,
					     const Boundary<1> &,
					     const FunctionMap &,
					     const dVector &,
					     dVector &) const {
  Assert(false, ExcNotImplemented());
};



template <int dim>
void KellyErrorEstimator<dim>::estimate_error (const DoFHandler<dim>    &dof,
					       const Quadrature<dim-1>  &quadrature,
					       const FiniteElement<dim> &fe,
					       const Boundary<dim>      &boundary,
					       const FunctionMap        &neumann_bc,
					       const dVector            &solution,
					       dVector                  &error) const {
  Assert (neumann_bc.find(255) == neumann_bc.end(),
	  ExcInvalidBoundaryIndicator());

				   // reserve one slot for each cell and set
				   // it to zero
  error.reinit (dof.get_tria().n_active_cells());

				   // number of integration points per face
  const unsigned int n_q_points = quadrature.n_quadrature_points;
				   // number of dofs per cell
  const unsigned int n_dofs = fe.total_dofs;
  
				   // make up a fe face values object for the
				   // restriction of the finite element function
				   // to a face, for the present cell and its
				   // neighbor.
  FEFaceValues<dim> fe_face_values_cell (fe, quadrature,
					 UpdateFlags(update_gradients | update_JxW_values |
						     update_jacobians | update_q_points |
						     update_normal_vectors)); 
  FEFaceValues<dim> fe_face_values_neighbor (fe, quadrature,
					     UpdateFlags(update_gradients)); 
  DoFHandler<dim>::active_cell_iterator cell = dof.begin_active(),
					endc = dof.end();
  
				   // loop over all cells
  for (unsigned int present_cell=0; cell!=endc; ++cell, ++present_cell)
				     // loop over all faces of this cell
    for (unsigned int face_no=0; face_no<2*dim; ++face_no)
      {
	const unsigned char boundary_indicator = cell->face(face_no)->boundary_indicator();
	if ((boundary_indicator != 255) &&
	    neumann_bc.find(boundary_indicator)==neumann_bc.end())
					   // this face is part of the boundary
					   // but not of the neumann boundary
					   // -> nothing to do
	  continue;

	
					 // initialize data of the restriction
					 // of this cell to the present face
	fe_face_values_cell.reinit (cell, face_no, fe, boundary);

					 // set up a vector of the gradients
					 // of the finite element function
					 // on this cell at the quadrature
					 // points
					 //
					 // let psi be a short name for
					 // [grad u_h]
	vector<Point<dim> > psi(n_q_points);
	
					 // get a list of the values of
					 // the solution for the ansatz
					 // functions on this cell
	vector<double>   dof_values;
	cell->get_dof_values (solution, dof_values);

					 // get a list of the gradients of
					 // the ansatz functions on this
					 // cell at the quadrature points
	const vector<vector<Point<dim> > > &shape_grads(fe_face_values_cell.get_shape_grads());

					 // compute the gradients of the solution
					 // at the quadrature points by summing
					 // over the ansatz functions.
	for (unsigned int j=0; j<n_q_points; ++j) 
	  for (unsigned int i=0; i<n_dofs; ++i)
	    psi[j] += dof_values[i]*shape_grads[i][j];


	
					 // now compute over the other side of
					 // the face
	if (boundary_indicator == 255)
					   // internal face; integrate jump
					   // of gradient across this face
	  {
	    Assert (cell->neighbor(face_no).state() == valid,
		    ExcInternalError());
	    unsigned int neighbor_neighbor;
	    DoFHandler<dim>::active_cell_iterator neighbor = cell->neighbor(face_no);

					     // find which number the current
					     // face has relative to the neighboring
					     // cell
	    for (neighbor_neighbor=0; neighbor_neighbor<2*dim; ++neighbor_neighbor)
	      if (neighbor->neighbor(neighbor_neighbor) == cell)
		break;

	    Assert (neighbor_neighbor<dim*2, ExcInternalError());

					     // get restriction of finite element
					     // function of #neighbor# to the
					     // common face.
	    fe_face_values_neighbor.reinit (neighbor, neighbor_neighbor, fe, boundary);

					     // get a list of the values of
					     // the solution for the ansatz
					     // functions on this neighbor
	    neighbor->get_dof_values (solution, dof_values);
					     // get a list of the gradients of the
					     // 
	    const vector<vector<Point<dim> > > &neighbor_grads (fe_face_values_cell.
								get_shape_grads());
					     // subtract the gradients of the
					     // solution on the neigbor cell
					     // at the quadrature points from
					     // those of the present cell
	    for (unsigned int j=0; j<n_q_points; ++j) 
	      for (unsigned int i=0; i<n_dofs; ++i)
		psi[j] -= dof_values[i]*neighbor_grads[i][j];	    
	  };



	
					 // now psi contains the following:
					 // - for an internal face, psi=[grad u]
					 // - for a neumann boundary face,
					 //   psi=grad u
					 // each component being the
					 // mentioned value at one of the
					 // quadrature points

					 // next we have to multiply this with
					 // the normal vector. Since we have
					 // taken the difference of gradients
					 // for internal faces, we may chose
					 // the normal vector of one cell,
					 // taking that of the neighbor
					 // would only change the sign. We take
					 // the outward normal.
					 //
					 // let phi be the name of the integrand
	vector<double> phi(n_q_points,0);
	const vector<Point<dim> > &normal_vectors(fe_face_values_cell.get_normal_vectors());
	for (unsigned int point=0; point<n_q_points; ++point)
	  phi[point] = psi[point]*normal_vectors[point];
	
	
	if (boundary_indicator != 255)
					   // neumann boundary face. compute
					   // difference between normal
					   // derivative and boundary function
	  {
					     // get the values of the boundary
					     // function at the quadrature
					     // points
	    vector<double> g(n_q_points);
	    neumann_bc.find(boundary_indicator)->second
	      ->value_list (fe_face_values_cell.get_quadrature_points(),
			    g);

	    for (unsigned int point=0; point<n_q_points; ++point)
	      phi[point] -= g[point];
	  };

	
					 // now phi contains the following:
					 // - for an internal face, phi=[du/dn]
					 // - for a neumann boundary face,
					 //   phi=du/dn-g
					 // each component being the
					 // mentioned value at one of the
					 // quadrature points

					 // take the square of the phi[i]
					 // for integration
	transform (phi.begin(), phi.end(),
		   phi.begin(), ptr_fun(sqr));

					 // perform integration by multiplication
					 // with weights and summation. Add the
					 // contribution of this face to the
					 // estimator of this cell
	
	error(present_cell)
	  += sqrt(inner_product (phi.begin(), phi.end(),
				 fe_face_values_cell.get_JxW_values().begin(),
				 0.0));
      };
};





// explicit instantiations

template class KellyErrorEstimator<1>;
template class KellyErrorEstimator<2>;
