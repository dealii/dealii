//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2004, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <base/quadrature.h>
#include <base/quadrature_lib.h>
#include <base/qprojector.h>
#include <base/table.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe.h>
#include <fe/mapping.h>
#include <fe/fe_raviart_thomas.h>
#include <fe/fe_values.h>
#include <fe/fe_tools.h>

#include <sstream>
#include <iostream>


template <int dim>
FE_RaviartThomas<dim>::FE_RaviartThomas (const unsigned int deg)
		:
		FE_PolyTensor<PolynomialsRaviartThomas<dim>, dim> (
		  deg,
		  FiniteElementData<dim>(get_dpo_vector(deg),
					 dim, deg+1, FiniteElementData<dim>::Hdiv, 1),
		  std::vector<bool>(PolynomialsRaviartThomas<dim>::compute_n_pols(deg), true),
		  std::vector<std::vector<bool> >(PolynomialsRaviartThomas<dim>::compute_n_pols(deg),
						  std::vector<bool>(dim,true))),
		rt_order(deg)
{
  Assert (dim >= 2, ExcImpossibleInDim(dim));
  const unsigned int n_dofs = this->dofs_per_cell;
  
  this->mapping_type = this->contravariant;
				   // First, initialize the
				   // generalized support points and
				   // quadrature weights, since they
				   // are required for interpolation.
  initialize_support_points(deg);
				   // Now compute the inverse node
				   //matrix, generating the correct
				   //basis functions from the raw
				   //ones.
  FullMatrix<double> M(n_dofs, n_dofs);
  FETools::compute_node_matrix(M, *this);
  this->inverse_node_matrix.reinit(n_dofs, n_dofs);
  this->inverse_node_matrix.invert(M);
				   // From now on, the shape functions
				   // will be the correct ones, not
				   // the raw shape functions anymore.
  

				   // initialize the various matrices
  for (unsigned int i=0; i<GeometryInfo<dim>::children_per_cell; ++i)
    {
      this->prolongation[i].reinit (n_dofs, n_dofs);
      this->restriction[i].reinit (n_dofs, n_dofs);
    }
  
  FETools::compute_embedding_matrices (*this, &this->prolongation[0]);
  initialize_restriction();
  
  std::vector<FullMatrix<double> >
    face_embeddings(1<<(dim-1), FullMatrix<double>(this->dofs_per_face,
						   this->dofs_per_face));
  FETools::compute_face_embedding_matrices(*this, &face_embeddings[0], 0, 0);
  this->interface_constraints.reinit((1<<(dim-1)) * this->dofs_per_face,
				     this->dofs_per_face);
  unsigned int target_row=0;
  for (unsigned int d=0;d<face_embeddings.size();++d)
    for (unsigned int i=0;i<face_embeddings[d].m();++i)
      {
	for (unsigned int j=0;j<face_embeddings[d].n();++j)
	  this->interface_constraints(target_row,j) = face_embeddings[d](i,j);
	++target_row;
      }
}



template <int dim>
std::string
FE_RaviartThomas<dim>::get_name () const
{
				   // note that the
				   // FETools::get_fe_from_name
				   // function depends on the
				   // particular format of the string
				   // this function returns, so they
				   // have to be kept in synch

  std::ostringstream namebuf;  
  namebuf << "FE_RaviartThomas<" << dim << ">(" << rt_order << ")";

  return namebuf.str();
}



template <int dim>
FiniteElement<dim> *
FE_RaviartThomas<dim>::clone() const
{
  return new FE_RaviartThomas<dim>(rt_order);
}


//---------------------------------------------------------------------------
// Auxiliary and internal functions
//---------------------------------------------------------------------------


#if deal_II_dimension == 1

template <int dim>
void
FE_RaviartThomas<dim>::initialize_support_points (const unsigned int deg)
{
  return;
  
  Assert (false, ExcNotImplemented());
  
  QGauss<dim> cell_quadrature(deg+1);
  const unsigned int n_interior_points
    = (deg>0) ? cell_quadrature.n_quadrature_points : 0;
  
  this->generalized_support_points.resize (2 + n_interior_points);
  
				   // Number of the point being entered
  unsigned int current = 0;

  
  if (deg==0) return;

  interior_weights.reinit(TableIndices<3>(2+n_interior_points, 0, dim));
  
  for (unsigned int k=0;k<cell_quadrature.n_quadrature_points;++k)
    this->generalized_support_points[current++] = cell_quadrature.point(k);
  
  Assert (current == this->generalized_support_points.size(),
	  ExcInternalError());
}

#else

// Version for 2d and higher. See above for 1d version
template <int dim>
void
FE_RaviartThomas<dim>::initialize_support_points (const unsigned int deg)
{
  QGauss<dim> cell_quadrature(deg+1);
  const unsigned int n_interior_points
    = (deg>0) ? cell_quadrature.n_quadrature_points : 0;
  
  unsigned int n_face_points = (dim>1) ? 1 : 0;
				   // compute (deg+1)^(dim-1)
  for (unsigned int d=1;d<dim;++d)
    n_face_points *= deg+1;

  
  this->generalized_support_points.resize (GeometryInfo<dim>::faces_per_cell*n_face_points
					   + n_interior_points);
  this->generalized_face_support_points.resize (n_face_points);
  
				   // Number of the point being entered
  unsigned int current = 0;

  if (dim>1)
    {
      QGauss<dim-1> face_points (deg+1);
      TensorProductPolynomials<dim-1> legendre
	= Polynomials::Legendre::generate_complete_basis(deg);

      boundary_weights.reinit(n_face_points, legendre.n());
      
//       Assert (face_points.n_quadrature_points == this->dofs_per_face,
// 	      ExcInternalError());
      
      for (unsigned int k=0;k<n_face_points;++k)
	{
	  this->generalized_face_support_points[k] = face_points.point(k);
					   // Compute its quadrature
					   // contribution for each
					   // moment.
	  for (unsigned int i=0;i<legendre.n();++i)
	    {
	      boundary_weights(k, i)
		= face_points.weight(k)
		* legendre.compute_value(i, face_points.point(k));
	    }
	}

      Quadrature<dim> faces = QProjector<dim>::project_to_all_faces(face_points);
      for (;current<GeometryInfo<dim>::faces_per_cell*n_face_points;
	   ++current)
	{
					   // Enter the support point
					   // into the vector
	  this->generalized_support_points[current] = faces.point(current);
	}
    }
  
  if (deg==0) return;
  
				   // Create Legendre basis for the
				   // space D_xi Q_k
  std::vector<AnisotropicPolynomials<dim>* > polynomials(dim);
  for (unsigned int dd=0;dd<dim;++dd)
    {
      std::vector<std::vector<Polynomials::Polynomial<double> > > poly(dim);
      for (unsigned int d=0;d<dim;++d)
	poly[d] = Polynomials::Legendre::generate_complete_basis(deg);
      poly[dd] = Polynomials::Legendre::generate_complete_basis(deg-1);

      polynomials[dd] = new AnisotropicPolynomials<dim>(poly);
    }
  
  interior_weights.reinit(TableIndices<3>(n_interior_points, polynomials[0]->n(), dim));
  
  for (unsigned int k=0;k<cell_quadrature.n_quadrature_points;++k)
    {
      this->generalized_support_points[current++] = cell_quadrature.point(k);
      for (unsigned int i=0;i<polynomials[0]->n();++i)
	for (unsigned int d=0;d<dim;++d)
	  interior_weights(k,i,d) = cell_quadrature.weight(k)
				    * polynomials[d]->compute_value(i,cell_quadrature.point(k));
    }

  for (unsigned int d=0;d<dim;++d)
    delete polynomials[d];
  
  Assert (current == this->generalized_support_points.size(),
	  ExcInternalError());
}

#endif


#if deal_II_dimension == 1

template <int dim>
void
FE_RaviartThomas<dim>::initialize_restriction()
{
  for (unsigned int i=0;i<GeometryInfo<dim>::children_per_cell;++i)
    this->restriction[i].reinit(0,0);
}

#else

// This function is the same Raviart-Thomas interpolation performed by
// interpolate. Still, we cannot use interpolate, since it was written
// for smooth functions. Thefunctions interpolated here are not
// smooth, maybe even not continuous. Therefore, we must double the
// number of quadrature points in each direction in order to integrate
// only smooth functions.

// Then again, the interpolated function is chosen such that the
// moments coincide with the function to be interpolated.

template <int dim>
void
FE_RaviartThomas<dim>::initialize_restriction()
{
  QGauss<dim-1> q_base (rt_order+1);
  const unsigned int n_face_points = q_base.n_quadrature_points;
				   // First, compute interpolation on
				   // subfaces
  for (unsigned int face=0;face<GeometryInfo<dim>::faces_per_cell;++face)
    {
				       // The shape functions of the
				       // child cell are evaluated
				       // in the quadrature points
				       // of a full face.
      Quadrature<dim> q_face
	= QProjector<dim>::project_to_face(q_base, face);
				       // Store shape values, since the
				       // evaluation suffers if not
				       // ordered by point
      Table<2,double> cached_values(this->dofs_per_cell, q_face.n_quadrature_points);
      for (unsigned int k=0;k<q_face.n_quadrature_points;++k)
	for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
	  cached_values(i,k)
	    = this->shape_value_component(i, q_face.point(k),
					  GeometryInfo<dim>::unit_normal_direction[face]);

      for (unsigned int sub=0;sub<GeometryInfo<dim>::subfaces_per_face;++sub)
	{
					   // The weight fuctions for
					   // the coarse face are
					   // evaluated on the subface
					   // only.
	  Quadrature<dim> q_sub
	    = QProjector<dim>::project_to_subface(q_base, face, sub);
	  const unsigned int child
	    = GeometryInfo<dim>::child_cell_on_face(face, sub);
	  
					   // On a certain face, we must
					   // compute the moments of ALL
					   // fine level functions with
					   // the coarse level weight
					   // functions belonging to
					   // that face. Due to the
					   // orthogonalization process
					   // when building the shape
					   // functions, these weights
					   // are equal to the
					   // corresponding shpe
					   // functions.
	  for (unsigned int k=0;k<n_face_points;++k)
	    for (unsigned int i_child = 0; i_child < this->dofs_per_cell; ++i_child)
	      for (unsigned int i_face = 0; i_face < this->dofs_per_face; ++i_face)
		{
						   // The quadrature
						   // weights on the
						   // subcell are NOT
						   // transformed, so we
						   // have to do it here.
		  this->restriction[child](face*this->dofs_per_face+i_face,
				     i_child)
		    += std::pow(.5, dim-1.) * q_sub.weight(k)
		    * cached_values(i_child, k)
		    * this->shape_value_component(face*this->dofs_per_face+i_face,
						  q_sub.point(k),
						  GeometryInfo<dim>::unit_normal_direction[face]);
		}
	}
    }
  
  if (rt_order==0) return;
  
				   // Create Legendre basis for the
				   // space D_xi Q_k. Here, we cannot
				   // use the shape functions
  std::vector<AnisotropicPolynomials<dim>* > polynomials(dim);
  for (unsigned int dd=0;dd<dim;++dd)
    {
      std::vector<std::vector<Polynomials::Polynomial<double> > > poly(dim);
      for (unsigned int d=0;d<dim;++d)
	poly[d] = Polynomials::Legendre::generate_complete_basis(rt_order);
      poly[dd] = Polynomials::Legendre::generate_complete_basis(rt_order-1);

      polynomials[dd] = new AnisotropicPolynomials<dim>(poly);
    }
  
  QGauss<dim> q_cell(rt_order+1);
  const unsigned int start_cell_dofs
    = GeometryInfo<dim>::faces_per_cell*this->dofs_per_face;

				   // Store shape values, since the
				   // evaluation suffers if not
				   // ordered by point
  Table<3,double> cached_values(this->dofs_per_cell, q_cell.n_quadrature_points, dim);
  for (unsigned int k=0;k<q_cell.n_quadrature_points;++k)
    for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
      for (unsigned int d=0;d<dim;++d)
	cached_values(i,k,d) = this->shape_value_component(i, q_cell.point(k), d);
  
  for (unsigned int child=0;child<GeometryInfo<dim>::children_per_cell;++child)
    {
      Quadrature<dim> q_sub = QProjector<dim>::project_to_child(q_cell, child);
      
      for (unsigned int k=0;k<q_sub.n_quadrature_points;++k)
	for (unsigned int i_child = 0; i_child < this->dofs_per_cell; ++i_child)
	  for (unsigned int d=0;d<dim;++d)
	    for (unsigned int i_weight=0;i_weight<polynomials[d]->n();++i_weight)
	      {
		this->restriction[child](start_cell_dofs+i_weight*dim+d,
				   i_child)
		  += q_sub.weight(k)
		  * cached_values(i_child, k, d)
		  * polynomials[d]->compute_value(i_weight, q_sub.point(k));
	      }
    }
  
  for (unsigned int d=0;d<dim;++d)
    delete polynomials[d];
}

#endif

#if deal_II_dimension == 1

template <>
std::vector<unsigned int>
FE_RaviartThomas<1>::get_dpo_vector (const unsigned int)
{
  Assert (false, ExcImpossibleInDim(1));
  return std::vector<unsigned int>();
}

#endif


template <int dim>
std::vector<unsigned int>
FE_RaviartThomas<dim>::get_dpo_vector (const unsigned int rt_order)
{
                                   // the element is face-based (not
                                   // to be confused with George
                                   // W. Bush's Faith Based
                                   // Initiative...), and we have
                                   // (rt_order+1)^(dim-1) DoFs per face
  unsigned int dofs_per_face = 1;
  for (unsigned int d=0; d<dim-1; ++d)
    dofs_per_face *= rt_order+1;

                                   // and then there are interior dofs
  const unsigned int
    interior_dofs = dim*rt_order*dofs_per_face;
  
  std::vector<unsigned int> dpo(dim+1);
  dpo[dim-1] = dofs_per_face;
  dpo[dim]   = interior_dofs;
  
  return dpo;
}



template <int dim>
UpdateFlags
FE_RaviartThomas<dim>::update_once (const UpdateFlags) const
{
				   // even the values have to be
				   // computed on the real cell, so
				   // nothing can be done in advance
  return update_default;
}



template <int dim>
UpdateFlags
FE_RaviartThomas<dim>::update_each (const UpdateFlags flags) const
{
  UpdateFlags out = update_default;

  if (flags & update_values)
    out |= update_values             | update_covariant_transformation
                                     | update_contravariant_transformation 
                                     | update_JxW_values;
  if (flags & update_gradients)
    out |= update_gradients          | update_covariant_transformation 
                                     | update_contravariant_transformation
                                     | update_JxW_values;
  //TODO: Set update flags appropriately and figure out, how the second
  // derivatives for the RT elements can be computed correctly.
  if (flags & update_second_derivatives)
    out |= update_second_derivatives | update_contravariant_transformation;

  return out;
}

//---------------------------------------------------------------------------
// Data field initialization
//---------------------------------------------------------------------------




template <int dim>
unsigned int
FE_RaviartThomas<dim>::n_base_elements () const
{
  return 1;
}



template <int dim>
const FiniteElement<dim> &
FE_RaviartThomas<dim>::base_element (const unsigned int index) const
{
  Assert (index==0, ExcIndexRange(index, 0, 1));
  return *this;
}



template <int dim>
unsigned int
FE_RaviartThomas<dim>::element_multiplicity (const unsigned int index) const
{
  Assert (index==0, ExcIndexRange(index, 0, 1));
  return 1;
}



template <int dim>
bool
FE_RaviartThomas<dim>::has_support_on_face (const unsigned int shape_index,
                                            const unsigned int face_index) const
{
  Assert (shape_index < this->dofs_per_cell,
	  ExcIndexRange (shape_index, 0, this->dofs_per_cell));
  Assert (face_index < GeometryInfo<dim>::faces_per_cell,
	  ExcIndexRange (face_index, 0, GeometryInfo<dim>::faces_per_cell));

				   // Return computed values if we
				   // know them easily. Otherwise, it
				   // is always safe to return true.
  switch (rt_order)
    {
      case 0:
      {
        switch (dim)
          {
            case 2:
            {
                                               // only on the one
                                               // non-adjacent face
                                               // are the values
                                               // actually zero. list
                                               // these in a table
              return (face_index != GeometryInfo<dim>::opposite_face[shape_index]);
            }
            
            default:
	      return true;
          };
      };
      
      default:  // other rt_order
	return true;
    };
  
  return true;
}



template <int dim>
void
FE_RaviartThomas<dim>::interpolate(
  std::vector<double>&,
  const std::vector<double>&) const
{
  Assert(false, ExcNotImplemented());
}



template <int dim>
void
FE_RaviartThomas<dim>::interpolate(
  std::vector<double>&    local_dofs,
  const std::vector<Vector<double> >& values,
  unsigned int offset) const
{
  Assert (values.size() == this->generalized_support_points.size(),
	  ExcDimensionMismatch(values.size(), this->generalized_support_points.size()));
  Assert (local_dofs.size() == this->dofs_per_cell,
	  ExcDimensionMismatch(local_dofs.size(),this->dofs_per_cell));
  Assert (values[0].size() >= offset+this->n_components(),
	  ExcDimensionMismatch(values[0].size(),offset+this->n_components()));

  std::fill(local_dofs.begin(), local_dofs.end(), 0.);

  const unsigned int n_face_points = boundary_weights.size(0);
  for (unsigned int face=0;face<GeometryInfo<dim>::faces_per_cell;++face)
    for (unsigned int k=0;k<n_face_points;++k)
      for (unsigned int i=0;i<boundary_weights.size(1);++i)
      {
	local_dofs[i+face*this->dofs_per_face] += boundary_weights(k,i)
			 * values[face*n_face_points+k](GeometryInfo<dim>::unit_normal_direction[face]+offset);
      }
  
  const unsigned start_cell_dofs = GeometryInfo<dim>::faces_per_cell*this->dofs_per_face;
  const unsigned start_cell_points = GeometryInfo<dim>::faces_per_cell*n_face_points;
  
  for (unsigned int k=0;k<interior_weights.size(0);++k)
    for (unsigned int i=0;i<interior_weights.size(1);++i)
      for (unsigned int d=0;d<dim;++d)
	local_dofs[start_cell_dofs+i*dim+d] += interior_weights(k,i,d) * values[k+start_cell_points](d+offset);
}


template <int dim>
void
FE_RaviartThomas<dim>::interpolate(
  std::vector<double>& local_dofs,
  const VectorSlice<const std::vector<std::vector<double> > >& values) const
{
  Assert (values.size() == this->n_components(),
	  ExcDimensionMismatch(values.size(), this->n_components()));
  Assert (values[0].size() == this->generalized_support_points.size(),
	  ExcDimensionMismatch(values[0].size(), this->generalized_support_points.size()));
  Assert (local_dofs.size() == this->dofs_per_cell,
	  ExcDimensionMismatch(local_dofs.size(),this->dofs_per_cell));

  std::fill(local_dofs.begin(), local_dofs.end(), 0.);

  const unsigned int n_face_points = boundary_weights.size(0);
  for (unsigned int face=0;face<GeometryInfo<dim>::faces_per_cell;++face)
    for (unsigned int k=0;k<n_face_points;++k)
      for (unsigned int i=0;i<boundary_weights.size(1);++i)
      {
	local_dofs[i+face*this->dofs_per_face] += boundary_weights(k,i)
			 * values[GeometryInfo<dim>::unit_normal_direction[face]][face*n_face_points+k];
      }
  
  const unsigned start_cell_dofs = GeometryInfo<dim>::faces_per_cell*this->dofs_per_face;
  const unsigned start_cell_points = GeometryInfo<dim>::faces_per_cell*n_face_points;
  
  for (unsigned int k=0;k<interior_weights.size(0);++k)
    for (unsigned int i=0;i<interior_weights.size(1);++i)
      for (unsigned int d=0;d<dim;++d)
	local_dofs[start_cell_dofs+i*dim+d] += interior_weights(k,i,d) * values[d][k+start_cell_points];
}


template <int dim>
unsigned int
FE_RaviartThomas<dim>::memory_consumption () const
{
  Assert (false, ExcNotImplemented ());
  return 0;
}


template class FE_RaviartThomas<deal_II_dimension>;
