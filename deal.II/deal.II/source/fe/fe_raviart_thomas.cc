//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2004, 2005 by the deal.II authors
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

#ifdef HAVE_STD_STRINGSTREAM
#  include <sstream>
#else
#  include <strstream>
#endif

#include <iostream>
using namespace std;


template <int dim>
FE_RaviartThomas<dim>::FE_RaviartThomas (const unsigned int deg)
		:
		FE_PolyTensor<PolynomialsRaviartThomas<dim>, dim> (
		  deg,
		  FiniteElementData<dim>(get_dpo_vector(deg),
					 dim, deg+1, FiniteElementData<dim>::Hdiv),
		  get_ria_vector (deg),
		  std::vector<std::vector<bool> >(
		    FiniteElementData<dim>(get_dpo_vector(deg),
					   dim,deg+1).dofs_per_cell,
		    std::vector<bool>(dim,true))),
		rt_order(deg)
{
  Assert (dim >= 2, ExcImpossibleInDim(dim));
  const unsigned int n_dofs = this->dofs_per_cell;
  
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
    this->prolongation[i].reinit (n_dofs,
				  n_dofs);
  FETools::compute_embedding_matrices (*this, &this->prolongation[0]);
  
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
//TODO:[WB] What is this?
                                   // then make
                                   // system_to_component_table
                                   // invalid, since this has no
                                   // meaning for the present element
  std::vector<std::pair<unsigned,unsigned> > tmp1, tmp2;
  this->system_to_component_table.swap (tmp1);
  this->face_system_to_component_table.swap (tmp2);
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

#ifdef HAVE_STD_STRINGSTREAM
  std::ostringstream namebuf;
#else
  std::ostrstream namebuf;
#endif
  
  namebuf << "FE_RaviartThomas<" << dim << ">(" << rt_order << ")";

#ifndef HAVE_STD_STRINGSTREAM
  namebuf << std::ends;
#endif
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



#if deal_II_dimension == 1

template <>
std::vector<bool>
FE_RaviartThomas<1>::get_ria_vector (const unsigned int)
{
  Assert (false, ExcImpossibleInDim(1));
  return std::vector<bool>();
}

#endif


template <int dim>
std::vector<bool>
FE_RaviartThomas<dim>::get_ria_vector (const unsigned int rt_order)
{
  unsigned int dofs_per_cell, dofs_per_face;
  switch (dim)
    {
      case 2:
	    dofs_per_face = rt_order+1;
	    dofs_per_cell = 2*(rt_order+1)*(rt_order+2);
	    break;
      case 3:
	    dofs_per_face = (rt_order+1)*(rt_order+1);
	    dofs_per_cell = 3*(rt_order+1)*(rt_order+1)*(rt_order+2);
	    break;
      default:
	    Assert (false, ExcNotImplemented());
    }
  Assert (FiniteElementData<dim>(get_dpo_vector(rt_order),dim).dofs_per_cell ==
	  dofs_per_cell,
	  ExcInternalError());
  Assert (FiniteElementData<dim>(get_dpo_vector(rt_order),dim).dofs_per_face ==
	  dofs_per_face,
	  ExcInternalError());
  
				   // all face dofs need to be
				   // non-additive, since they have
				   // continuity requirements.
				   // however, the interior dofs are
				   // made additive
  std::vector<bool> ret_val(dofs_per_cell,false);
  for (unsigned int i=GeometryInfo<dim>::faces_per_cell*dofs_per_face;
       i < dofs_per_cell; ++i)
    ret_val[i] = true;

  return ret_val;
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
    out |= update_values             | update_covariant_transformation;
  if (flags & update_gradients)
    out |= update_gradients          | update_covariant_transformation;
  if (flags & update_second_derivatives)
    out |= update_second_derivatives | update_covariant_transformation;

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
              const unsigned int
                opposite_faces[GeometryInfo<2>::faces_per_cell]
                = { 2, 3, 0, 1};
              
              return (face_index != opposite_faces[shape_index]);
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
