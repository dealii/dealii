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

				   // create a map of integrals indexed by
				   // the corresponding face. In this map
				   // we store the integrated jump of the
				   // gradient for each face. By doing so,
				   // we can check whether we have already
				   // integrated along this face by testing
				   // whether the respective face is already
				   // a key in this map.
				   // At the end of the function, we again
				   // loop over the cells and collect the
				   // conrtibutions of the different faces
				   // of the cell.
  map<DoFHandler<dim>::face_iterator, double> face_integrals;
  
				   // number of integration points per face
  const unsigned int n_q_points = quadrature.n_quadrature_points;
  
				   // make up a fe face values object for the
				   // restriction of the finite element function
				   // to a face, for the present cell and its
				   // neighbor.
  FEFaceValues<dim> fe_face_values_cell (fe, quadrature,
					 UpdateFlags(update_gradients  |
						     update_JxW_values |
						     update_jacobians  |
						     update_q_points   |
						     update_normal_vectors)); 
  FEFaceValues<dim> fe_face_values_neighbor (fe, quadrature,
					     UpdateFlags(update_gradients |
							 update_jacobians));

				   // loop variables
  DoFHandler<dim>::active_cell_iterator cell = dof.begin_active(),
					endc = dof.end();
				   // loop over all cells
  for (; cell!=endc; ++cell)
				     // loop over all faces of this cell
    for (unsigned int face_no=0; face_no<2*dim; ++face_no)
      {
					 // if we already visited this
					 // face: do nothing
	if (face_integrals.find(cell->face(face_no)) !=
	    face_integrals.end())
	  continue;

					 // if the neighboring cell is less
					 // refined than the present one, then
					 // do nothing since we integrate
					 // over the subfaces when we visit
					 // the coarse cells.
	if (cell->at_boundary(face_no) == false)
	  if (cell->neighbor(face_no)->level() < cell->level())
	    continue;
	
					 // if this face is part of the boundary
					 // but not of the neumann boundary
					 // -> nothing to do. However, to make
					 // things easier when summing up the
					 // contributions of the faces of cells,
					 // we enter this face into the list
					 // of faces with contribution zero.
	const unsigned char boundary_indicator
	  = cell->face(face_no)->boundary_indicator();
	if ((boundary_indicator != 255) &&
	    neumann_bc.find(boundary_indicator)==neumann_bc.end()) 
	  {
	    face_integrals[cell->face(face_no)] = 0;
	    continue;
	  };


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
	fe_face_values_cell.get_function_grads (solution, psi);
	
	
					 // now compute over the other side of
					 // the face
	if (boundary_indicator == 255)
					   // internal face; integrate jump
					   // of gradient across this face
	  {
	    Assert (cell->neighbor(face_no).state() == valid,
		    ExcInternalError());
	    unsigned int neighbor_neighbor;
	    DoFHandler<dim>::active_cell_iterator neighbor
	      = cell->neighbor(face_no);

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
	    fe_face_values_neighbor.reinit (neighbor, neighbor_neighbor,
					    fe, boundary);

					     // get a list of the gradients of
					     // the finite element solution
					     // restricted to the neighbor cell
	    vector<Point<dim> > neighbor_psi (n_q_points);
	    fe_face_values_neighbor.get_function_grads (solution, neighbor_psi);

					     // compute the jump in the gradients
	    transform (psi.begin(), psi.end(),
		       neighbor_psi.begin(),
		       psi.begin(),
		       minus<Point<dim> >());
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
	const vector<Point<dim> > &normal_vectors(fe_face_values_cell.
						  get_normal_vectors());
	
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
					 // with weights and summation.
	face_integrals[cell->face(face_no)]
	  = (inner_product (phi.begin(), phi.end(),
			    fe_face_values_cell.get_JxW_values().begin(),
			    0.0) *
	     cell->diameter() / 24);
      };


				   // finally add up the contributions of the
				   // faces for each cell
  
				   // reserve one slot for each cell and set
				   // it to zero
  error.reinit (dof.get_tria().n_active_cells());

  unsigned int present_cell=0;
  for (cell=dof.begin_active(); cell!=endc; ++cell, ++present_cell)
    {
				       // loop over all faces of this cell
      for (unsigned int face_no=0; face_no<2*dim; ++face_no) 
	{
	  Assert (face_integrals.find(cell->face(face_no)) !=
		  face_integrals.end(),
		  ExcInternalError());
	  error(present_cell) += face_integrals[cell->face(face_no)];
	};
	   
      error(present_cell) = sqrt(error(present_cell));
    };
};





// explicit instantiations

template class KellyErrorEstimator<1>;
template class KellyErrorEstimator<2>;
