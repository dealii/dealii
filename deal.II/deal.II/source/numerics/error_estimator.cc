/*      $Id$                 */


#include <fe/fe.h>
#include <fe/fe_values.h>
#include <fe/fe_update_flags.h>
#include <base/quadrature.h>
#include <base/quadrature_lib.h>
#include <numerics/error_estimator.h>
#include <grid/dof.h>
#include <grid/tria_iterator.h>
#include <grid/dof_accessor.h>
#include <grid/geometry_info.h>
#include <lac/vector.h>
#include <lac/vector.h>

#include <numeric>
#include <algorithm>
#include <cmath>



inline static double sqr (const double x) {
  return x*x;
};



#if deal_II_dimension == 1

template <>
void KellyErrorEstimator<1>::estimate (const DoFHandler<1>  &dof,
				       const Quadrature<0>  &,
				       const FunctionMap    &neumann_bc,
				       const Vector<double> &solution,
				       Vector<float>        &error,
				       const Function<1>    *coefficient,
				       const unsigned int    selected_component)
{
  Assert (selected_component < dof.get_fe().n_components,
	  ExcInvalidComponent (selected_component, dof.get_fe().n_components));
  Assert (coefficient->n_components == 1,
	  ExcInternalError());
  
  const unsigned int dim=1;

				   // reserve one slot for each cell and set
				   // it to zero
  error.reinit (dof.get_tria().n_active_cells());

				   // loop over all cells. note that the
				   // error indicator is only a sum over
				   // the two contributions from the two
				   // vertices of each cell.
  QTrapez<1> quadrature;
  FEValues<dim> fe_values (dof.get_fe(), quadrature, update_gradients);
  DoFHandler<dim>::active_cell_iterator cell = dof.begin_active();
  for (unsigned int cell_index=0; cell != dof.end(); ++cell, ++cell_index)
    {
				       // loop over te two points bounding
				       // this line. n==0 is left point,
				       // n==1 is right point
      for (unsigned int n=0; n<2; ++n)
	{
					   // find right active neighbor
	  DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(n);
	  if (neighbor.state() == valid)
	    while (neighbor->has_children())
	      neighbor = neighbor->child(n==0 ? 1 : 0);
      
					   // now get the gradients on the
					   // both sides of the point
	  vector<vector<Tensor<1,dim> > >
	    gradients (2, vector<Tensor<1,1> >(dof.get_fe().n_components));
	  
	  fe_values.reinit (cell);
	  fe_values.get_function_grads (solution, gradients);
	  const double grad_here = gradients[n][selected_component][0];

	  double grad_neighbor;
	  if (neighbor.state() == valid)
	    {
	      fe_values.reinit (neighbor);
	      fe_values.get_function_grads (solution, gradients);
	      grad_neighbor = gradients[n==0 ? 1 : 0][selected_component][0];
	    }
	  else
	    if (neumann_bc.find(n) != neumann_bc.end())
	      grad_neighbor = neumann_bc.find(n)->second->value(cell->vertex(0));
	    else
	      grad_neighbor = 0;
	    
	  const double jump = (grad_here - grad_neighbor) *
			      (coefficient != 0 ?
			       coefficient->value(cell->vertex(n)) :
			       1);
	  error(cell_index) += jump*jump * cell->diameter();
	};
      error(cell_index) = sqrt(error(cell_index));
    };
};
	
#endif


template <int dim>
void KellyErrorEstimator<dim>::estimate (const DoFHandler<dim>    &dof,
					 const Quadrature<dim-1>  &quadrature,
					 const FunctionMap        &neumann_bc,
					 const Vector<double>     &solution,
					 Vector<float>            &error,
					 const Function<dim>      *coefficient,
					 const unsigned int        selected_component)
{
  Assert (neumann_bc.find(255) == neumann_bc.end(),
	  ExcInvalidBoundaryIndicator());
  Assert (selected_component < dof.get_fe().n_components,
	  ExcInvalidComponent (selected_component, dof.get_fe().n_components));

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
  FaceIntegrals face_integrals;
  
				   // number of integration points per face
  const unsigned int n_q_points = quadrature.n_quadrature_points;
  
				   // make up a fe face values object for the
				   // restriction of the finite element function
				   // to a face, for the present cell and its
				   // neighbor. In principle we would only need
				   // one at a time, but this way we can
				   // have more fine grained access to what
				   // values really need to be computed (we
				   // need not compute all values on the
				   // neighbor cells, so using two objects
				   // gives us a performance gain).
  FEFaceValues<dim> fe_face_values_cell (dof.get_fe(), quadrature,
					 UpdateFlags(update_gradients      |
						     update_JxW_values     |
						     ((!neumann_bc.empty() ||
						       (coefficient != 0))  ?
						      update_q_points : 0) |
						     update_normal_vectors)); 
  FEFaceValues<dim> fe_face_values_neighbor (dof.get_fe(), quadrature, update_gradients);
  FESubfaceValues<dim> fe_subface_values (dof.get_fe(), quadrature, update_gradients);
  
				   // loop over all cells
  const active_cell_iterator endc = dof.end();
  for (active_cell_iterator cell = dof.begin_active(); cell!=endc; ++cell)
				     // loop over all faces of this cell
    for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
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

	
	if (cell->face(face_no)->has_children() == false)
					   // if the face is a regular one, i.e.
					   // either on the other side there is
					   // nirvana (face is at boundary), or
					   // the other side's refinement level
					   // is the same as that of this side,
					   // then handle the integration of
					   // these both cases together
	  integrate_over_regular_face (cell, face_no, neumann_bc,
				       n_q_points,
				       fe_face_values_cell,
				       fe_face_values_neighbor,
				       face_integrals,
				       solution,
				       dof.get_fe().n_components,
				       selected_component,
				       coefficient);
	else
					   // otherwise we need to do some
					   // special computations which do
					   // not fit into the framework of
					   // the above function
	  integrate_over_irregular_face (cell, face_no,
					 n_q_points,
					 fe_face_values_cell,
					 fe_subface_values,
					 face_integrals, solution,
					 dof.get_fe().n_components,
					 selected_component,
					 coefficient);
      };


				   // finally add up the contributions of the
				   // faces for each cell
  
				   // reserve one slot for each cell and set
				   // it to zero
  error.reinit (dof.get_tria().n_active_cells());

  unsigned int present_cell=0;
  for (active_cell_iterator cell=dof.begin_active(); cell!=endc; ++cell, ++present_cell)
    {
				       // loop over all faces of this cell
      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no) 
	{
	  Assert (face_integrals.find(cell->face(face_no)) !=
		  face_integrals.end(),
		  ExcInternalError());
	  error(present_cell) += (face_integrals[cell->face(face_no)] *
				  cell->diameter() / 24);
	};
	   
      error(present_cell) = sqrt(error(present_cell));
    };
};



#if deal_II_dimension == 1

template <>
void KellyErrorEstimator<1>::integrate_over_regular_face (const active_cell_iterator &,
							  const unsigned int      ,
							  const FunctionMap      &,
							  const unsigned int      ,
							  FEFaceValues<1>        &,
							  FEFaceValues<1>        &,
							  FaceIntegrals          &,
							  const Vector<double>   &,
							  const unsigned int      ,
							  const unsigned int      ,
							  const Function<1>      *) {
  Assert (false, ExcInternalError());
};



template <>
void KellyErrorEstimator<1>::
integrate_over_irregular_face (const active_cell_iterator &,
			       const unsigned int          ,
			       const unsigned int          ,
			       FEFaceValues<1>            &,
			       FESubfaceValues<1>         &,
			       FaceIntegrals              &,
			       const Vector<double>       &,
			       const unsigned int          ,
			       const unsigned int          ,
			       const Function<1>          *) {
  Assert (false, ExcInternalError());
};

#endif



template <int dim>
void KellyErrorEstimator<dim>::
integrate_over_regular_face (const active_cell_iterator &cell,
			     const unsigned int          face_no,
			     const FunctionMap          &neumann_bc,
			     const unsigned int          n_q_points,
			     FEFaceValues<dim>          &fe_face_values_cell,
			     FEFaceValues<dim>          &fe_face_values_neighbor,
			     FaceIntegrals              &face_integrals,
			     const Vector<double>       &solution,
			     const unsigned int          n_components,
			     const unsigned int          selected_component,
			     const Function<dim>        *coefficient) {
  const DoFHandler<dim>::face_iterator face = cell->face(face_no);
  
				   // initialize data of the restriction
				   // of this cell to the present face
  fe_face_values_cell.reinit (cell, face_no);
  
				   // set up a vector of the gradients
				   // of the finite element function
				   // on this cell at the quadrature
				   // points
				   //
				   // let psi be a short name for
				   // [a grad u_h]
  vector<Tensor<1,dim> > psi(n_q_points);
  if (n_components == 1)
    fe_face_values_cell.get_function_grads (solution, psi);
  else
    {
      vector<vector<Tensor<1,dim> > > tmp (n_q_points,
					   vector<Tensor<1,dim> >(n_components));
      fe_face_values_cell.get_function_grads (solution, tmp);
      for (unsigned int i=0; i<n_q_points; ++i)
	psi[i] = tmp[i][selected_component];
    };
  
  
  
				   // now compute over the other side of
				   // the face
  if (face->at_boundary() == false)
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
      for (neighbor_neighbor=0; neighbor_neighbor<GeometryInfo<dim>::faces_per_cell;
	   ++neighbor_neighbor)
	if (neighbor->neighbor(neighbor_neighbor) == cell)
	  break;
      
      Assert (neighbor_neighbor<GeometryInfo<dim>::faces_per_cell, ExcInternalError());
      
				       // get restriction of finite element
				       // function of #neighbor# to the
				       // common face.
      fe_face_values_neighbor.reinit (neighbor, neighbor_neighbor);

				       // get a list of the gradients of
				       // the finite element solution
				       // restricted to the neighbor cell
      vector<Tensor<1,dim> > neighbor_psi (n_q_points);
      if (n_components == 1)
	fe_face_values_neighbor.get_function_grads (solution, neighbor_psi);
      else
	{
	  vector<vector<Tensor<1,dim> > > tmp (n_q_points,
					       vector<Tensor<1,dim> >(n_components));
	  fe_face_values_neighbor.get_function_grads (solution, tmp);
	  for (unsigned int i=0; i<n_q_points; ++i)
	    psi[i] = tmp[i][selected_component];
	};

      
				       // compute the jump in the gradients
      transform (psi.begin(), psi.end(),
		 neighbor_psi.begin(),
		 psi.begin(),
		 minus<Tensor<1,dim> >());
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

				   // if a coefficient was given: use that
				   // to scale the jump in the gradient
  if (coefficient != 0)
    {
      vector<double>      coefficient_values (n_q_points);
      coefficient->value_list (fe_face_values_cell.get_quadrature_points(),
			       coefficient_values);
      for (unsigned int point=0; point<n_q_points; ++point)
	phi[point] *= coefficient_values[point];
    };
      
  
  if (face->at_boundary() == true)
				     // neumann boundary face. compute
				     // difference between normal
				     // derivative and boundary function
    {
      const unsigned char boundary_indicator = face->boundary_indicator();

      Assert (neumann_bc.find(boundary_indicator) != neumann_bc.end(),
	      ExcInternalError ());
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
				   // - for an internal face, phi=[a du/dn]
				   // - for a neumann boundary face,
				   //   phi=a du/dn-g
				   // each component being the
				   // mentioned value at one of the
				   // quadrature points
  
				   // take the square of the phi[i]
				   // for integration
  transform (phi.begin(), phi.end(),
	     phi.begin(), ptr_fun(sqr));
  
				   // perform integration by multiplication
				   // with weights and summation.
  face_integrals[face] = inner_product (phi.begin(), phi.end(),
					fe_face_values_cell.get_JxW_values().begin(),
					0.0);
};




template <int dim>
void KellyErrorEstimator<dim>::
integrate_over_irregular_face (const active_cell_iterator &cell,
			       const unsigned int          face_no,
			       const unsigned int          n_q_points,
			       FEFaceValues<dim>          &fe_face_values,
			       FESubfaceValues<dim>       &fe_subface_values,
			       FaceIntegrals              &face_integrals,
			       const Vector<double>       &solution,
			       const unsigned int          n_components,
			       const unsigned int          selected_component,
			       const Function<dim>        *coefficient) {
  const DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face_no);
  Assert (neighbor.state() == valid, ExcInternalError());
  Assert (neighbor->has_children(), ExcInternalError());
				   // set up a vector of the gradients
				   // of the finite element function
				   // on this cell at the quadrature
				   // points
				   //
				   // let psi be a short name for
				   // [a grad u_h]
  vector<Tensor<1,dim> > psi(n_q_points);

				   // store which number #cell# has in the
				   // list of neighbors of #neighbor#
  unsigned int neighbor_neighbor;
  for (neighbor_neighbor=0; neighbor_neighbor<GeometryInfo<dim>::faces_per_cell;
       ++neighbor_neighbor)
    if (neighbor->neighbor(neighbor_neighbor) == cell)
      break;
  Assert (neighbor_neighbor<GeometryInfo<dim>::faces_per_cell, ExcInternalError());
  
				   // loop over all subfaces
  for (unsigned int subface_no=0; subface_no<GeometryInfo<dim>::subfaces_per_face;
       ++subface_no)
    {
				       // get an iterator pointing to the
				       // cell behind the present subface
      const active_cell_iterator neighbor_child
	= neighbor->child(GeometryInfo<dim>::
			  child_cell_on_face(neighbor_neighbor,subface_no));
      Assert (neighbor_child->face(neighbor_neighbor) ==
	      cell->face(face_no)->child(subface_no),
	      ExcInternalError());
      Assert (!neighbor->child(GeometryInfo<dim>::
			       child_cell_on_face(neighbor_neighbor,subface_no))->has_children(),
	      ExcInternalError());
            
				       // restrict the finite element on the
				       // present cell to the subface and
				       // store the gradient of the solution
				       // in psi
      fe_subface_values.reinit (cell, face_no, subface_no);
      if (n_components == 1)
	fe_subface_values.get_function_grads (solution, psi);
      else
	{
	  vector<vector<Tensor<1,dim> > > tmp (n_q_points,
					       vector<Tensor<1,dim> >(n_components));
	  fe_subface_values.get_function_grads (solution, tmp);
	  for (unsigned int i=0; i<n_q_points; ++i)
	    psi[i] = tmp[i][selected_component];
	};

				       // restrict the finite element on the
				       // neighbor cell to the common #subface#.
				       // store the gradient in #neighbor_psi#
      vector<Tensor<1,dim> > neighbor_psi (n_q_points);
      fe_face_values.reinit (neighbor_child, neighbor_neighbor);
      if (n_components == 1)
	fe_face_values.get_function_grads (solution, neighbor_psi);
      else
	{
	  vector<vector<Tensor<1,dim> > > tmp (n_q_points,
					       vector<Tensor<1,dim> >(n_components));
	  fe_face_values.get_function_grads (solution, tmp);
	  for (unsigned int i=0; i<n_q_points; ++i)
	    psi[i] = tmp[i][selected_component];
	};
      
				       // compute the jump in the gradients
      transform (psi.begin(), psi.end(),
		 neighbor_psi.begin(),
		 psi.begin(),
		 minus<Tensor<1,dim> >());

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
      const vector<Point<dim> > &normal_vectors(fe_face_values.
						get_normal_vectors());
  
      for (unsigned int point=0; point<n_q_points; ++point)
	phi[point] = psi[point]*normal_vectors[point];
      
				       // if a coefficient was given: use that
				       // to scale the jump in the gradient
      if (coefficient != 0)
	{
	  vector<double>      coefficient_values (n_q_points);
	  coefficient->value_list (fe_face_values.get_quadrature_points(),
				   coefficient_values);
	  for (unsigned int point=0; point<n_q_points; ++point)
	    phi[point] *= coefficient_values[point];
	};

				   // take the square of the phi[i]
				       // for integration
      transform (phi.begin(), phi.end(),
		 phi.begin(), ptr_fun(sqr));

				       // perform integration by multiplication
				       // with weights and summation.
      face_integrals[neighbor_child->face(neighbor_neighbor)]
	= inner_product (phi.begin(), phi.end(),
			 fe_face_values.get_JxW_values().
			 begin(),
			 0.0);
    };

  
  				   // finally loop over all subfaces to
				   // collect the contributions of the
				   // subfaces and store them with the
				   // mother face
  double sum=0;
  DoFHandler<dim>::face_iterator face = cell->face(face_no);
  for (unsigned int subface_no=0; subface_no<GeometryInfo<dim>::subfaces_per_face;
       ++subface_no) 
    {
      Assert (face_integrals.find(face->child(subface_no)) !=
	      face_integrals.end(),
	      ExcInternalError());
      sum += face_integrals[face->child(subface_no)];
    };
  face_integrals[face] = sum;
};



// explicit instantiations

template class KellyErrorEstimator<deal_II_dimension>;
