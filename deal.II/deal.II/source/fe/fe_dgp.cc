//----------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------

#include <base/quadrature.h>
#include <base/polynomial.h>
#include <base/polynomial_space.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe.h>
#include <fe/mapping.h>
#include <fe/mapping_q1.h>
#include <fe/fe_dgp.h>
#include <fe/fe_values.h>



template <int dim>
FE_DGP<dim>::FE_DGP (unsigned int degree)
		:
		FiniteElement<dim> (FiniteElementData<dim>(get_dpo_vector(degree),1),
				    std::vector<bool>(1,true)),
		degree(degree),
		poly(0)
{
				   // create array of Legendre polynomials
  std::vector<Legendre<double> > v;
  for (unsigned int i=0;i<=degree;++i)
    v.push_back(Legendre<double>(i));
  poly = new PolynomialSpace<dim> (v);

				   // if defined, copy over matrices
				   // from precomputed arrays
//    if ((degree < Matrices::n_embedding_matrices) &&
//        (Matrices::embedding[degree] != 0))
//      {
//        prolongation[0].fill (Matrices::embedding[degree]);
//      }
//    else
//  				     // matrix undefined, set size to zero
//      for (unsigned int i=0;i<GeometryInfo<dim>::children_per_cell;++i)
//        prolongation[i].reinit(0);

//  				   // same as above: copy over matrix
//  				   // from predefined values and
//  				   // generate all others by rotation
//    if ((degree < Matrices::n_projection_matrices) &&
//        (Matrices::projection_matrices[degree] != 0))
//      {
//        restriction[0].fill (Matrices::projection_matrices[degree]);
//      }
//    else
//  				     // matrix undefined, set size to zero
//      for (unsigned int i=0;i<GeometryInfo<dim>::children_per_cell;++i)
//        restriction[i].reinit(0);  
};



template <int dim>
FE_DGP<dim>::~FE_DGP ()
{
				   // delete poly member and set it to
				   // zero to prevent accidental use
  delete poly;
  poly = 0;
}



template <int dim>
FiniteElement<dim> *
FE_DGP<dim>::clone() const
{
  return new FE_DGP<dim>(degree);
}



template <int dim>
double
FE_DGP<dim>::shape_value (const unsigned int i,
			  const Point<dim> &p) const
{
  return poly->compute_value(i, p);
}



template <int dim>
Tensor<1,dim>
FE_DGP<dim>::shape_grad (const unsigned int i,
			 const Point<dim> &p) const
{
  return poly->compute_grad(i, p);
}



template <int dim>
Tensor<2,dim>
FE_DGP<dim>::shape_grad_grad (const unsigned int i,
			      const Point<dim> &p) const
{
  return poly->compute_grad_grad(i, p);
}


//----------------------------------------------------------------------
// Auxiliary functions
//----------------------------------------------------------------------


template <int dim>
std::vector<unsigned int>
FE_DGP<dim>::get_dpo_vector(unsigned int deg)
{
  std::vector<unsigned int> dpo(dim+1, static_cast<unsigned int>(0));
  dpo[dim] = ++deg;
  for (unsigned int i=1;i<dim;++i)
    {
      dpo[dim] *= deg+i;
      dpo[dim] /= i+1;
    }
  return dpo;
}


template <int dim>
UpdateFlags
FE_DGP<dim>::update_once (const UpdateFlags flags) const
{
  UpdateFlags out = update_default;

  if (flags & update_values)
    out |= update_values;

  return out;
}


template <int dim>
UpdateFlags
FE_DGP<dim>::update_each (const UpdateFlags flags) const
{
  UpdateFlags out = update_default;

  if (flags & update_gradients)
    out |= update_gradients | update_covariant_transformation;

  if (flags & update_second_derivatives)
    out |= update_second_derivatives | update_covariant_transformation;

  return out;
}


//----------------------------------------------------------------------
// Data field initialization
//----------------------------------------------------------------------

template <int dim>
typename Mapping<dim>::InternalDataBase *
FE_DGP<dim>::get_data (const UpdateFlags      update_flags,
		       const Mapping<dim>    &mapping,
		       const Quadrature<dim> &quadrature) const
{
  InternalData* data = new InternalData;
  std::vector<double> values(0);
  std::vector<Tensor<1,dim> > grads(0);
  std::vector<Tensor<2,dim> > grad_grads(0);
  
				   // check what needs to be
				   // initialized only once and what
				   // on every cell/face/subface we
				   // visit
  data->update_once = update_once(update_flags);
  data->update_each = update_each(update_flags);
  data->update_flags = data->update_once | data->update_each;

  const UpdateFlags flags(data->update_flags);
  const unsigned int n_q_points = quadrature.n_quadrature_points;
  
  if (flags & update_values)
    {
      values.resize (dofs_per_cell);
      data->shape_values.resize(dofs_per_cell,
				std::vector<double>(n_q_points));
    }

  if (flags & update_gradients)
    {
      grads.resize (dofs_per_cell);
      data->shape_gradients.resize(dofs_per_cell,
				   std::vector<Tensor<1,dim> >(n_q_points));
    }

				   // if second derivatives through
				   // finite differencing is required,
				   // then initialize some objects for
				   // that
  if (flags & update_second_derivatives)
    data->initialize_2nd (this, mapping, quadrature);
  
  
  if (flags & (update_values | update_gradients))
    for (unsigned int i=0; i<n_q_points; ++i)
      {
	poly->compute(quadrature.point(i), values, grads, grad_grads);
	for (unsigned int k=0; k<dofs_per_cell; ++k)
	  {
	    if (flags & update_values)
	      data->shape_values[k][i] = values[k];
	    if (flags & update_gradients)
	      data->shape_gradients[k][i] = grads[k];
	  }
      }
  return data;
}



//----------------------------------------------------------------------
// Fill data of FEValues
//----------------------------------------------------------------------

template <int dim>
void
FE_DGP<dim>::fill_fe_values (const Mapping<dim>                   &mapping,
			     const typename DoFHandler<dim>::cell_iterator &cell,
			     const Quadrature<dim>                &quadrature,
			     typename Mapping<dim>::InternalDataBase       &mapping_data,
			     typename Mapping<dim>::InternalDataBase       &fedata,
			     FEValuesData<dim>                    &data) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  InternalData &fe_data = dynamic_cast<InternalData &> (fedata);
  
  const UpdateFlags flags(fe_data.current_update_flags());

  for (unsigned int k=0; k<dofs_per_cell; ++k)
    {
      for (unsigned int i=0;i<quadrature.n_quadrature_points;++i)
	if (flags & update_values)
	  data.shape_values(k,i) = fe_data.shape_values[k][i];
      
      if (flags & update_gradients)
	mapping.transform_covariant(data.shape_gradients[k],
				    fe_data.shape_gradients[k],
				    mapping_data, 0);
    }
  
  if (flags & update_second_derivatives)
    compute_2nd (mapping, cell, 0, mapping_data, fe_data, data);
  
  fe_data.first_cell = false;
}



template <int dim>
void
FE_DGP<dim>::fill_fe_face_values (const Mapping<dim>                   &mapping,
				  const typename DoFHandler<dim>::cell_iterator &cell,
				  const unsigned int                    face,
				  const Quadrature<dim-1>              &quadrature,
				  typename Mapping<dim>::InternalDataBase       &mapping_data,
				  typename Mapping<dim>::InternalDataBase       &fedata,
				  FEValuesData<dim>                    &data) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  InternalData &fe_data = dynamic_cast<InternalData &> (fedata);

				   // offset determines which data set
				   // to take (all data sets for all
				   // faces are stored contiguously)
  const unsigned int offset = face * quadrature.n_quadrature_points;
  
  const UpdateFlags flags(fe_data.update_once | fe_data.update_each);

  for (unsigned int k=0; k<dofs_per_cell; ++k)
    {
      for (unsigned int i=0;i<quadrature.n_quadrature_points;++i)
	if (flags & update_values)
	  data.shape_values(k,i) = fe_data.shape_values[k][i+offset];
      
      if (flags & update_gradients)
	mapping.transform_covariant(data.shape_gradients[k],
				    fe_data.shape_gradients[k],
				    mapping_data, offset);
    }

  if (flags & update_second_derivatives)
    compute_2nd (mapping, cell, offset, mapping_data, fe_data, data);
  
  fe_data.first_cell = false;
}



template <int dim>
void
FE_DGP<dim>::fill_fe_subface_values (const Mapping<dim>                   &mapping,
				     const typename DoFHandler<dim>::cell_iterator &cell,
				     const unsigned int                    face,
				     const unsigned int                    subface,
				     const Quadrature<dim-1>              &quadrature,
				     typename Mapping<dim>::InternalDataBase       &mapping_data,
				     typename Mapping<dim>::InternalDataBase       &fedata,
				     FEValuesData<dim>                    &data) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  InternalData &fe_data = dynamic_cast<InternalData &> (fedata);
  
				   // offset determines which data set
				   // to take (all data sets for all
				   // sub-faces are stored contiguously)
  const unsigned int offset = (face * GeometryInfo<dim>::subfaces_per_face + subface)
			      * quadrature.n_quadrature_points;

  const UpdateFlags flags(fe_data.update_once | fe_data.update_each);

  for (unsigned int k=0; k<dofs_per_cell; ++k)
    {
      for (unsigned int i=0;i<quadrature.n_quadrature_points;++i)
	if (flags & update_values)
	  data.shape_values(k,i) = fe_data.shape_values[k][i+offset];
      
      if (flags & update_gradients)
	mapping.transform_covariant(data.shape_gradients[k],
				    fe_data.shape_gradients[k],
				    mapping_data, offset);
    }
  
  if (flags & update_second_derivatives)
    compute_2nd (mapping, cell, offset, mapping_data, fe_data, data);
  
  fe_data.first_cell = false;
}



template <int dim>
unsigned int
FE_DGP<dim>::n_base_elements () const
{
  return 1;
};



template <int dim>
const FiniteElement<dim> &
FE_DGP<dim>::base_element (const unsigned int index) const
{
  Assert (index==0, ExcIndexRange(index, 0, 1));
  return *this;
};



template <int dim>
bool
FE_DGP<dim>::has_support_on_face (const unsigned int,
				  const unsigned int) const
{
  return true;
}



template <int dim>
unsigned int
FE_DGP<dim>::memory_consumption () const
{
  Assert (false, ExcNotImplemented ());
  return 0;
}



template <int dim>
unsigned int
FE_DGP<dim>::get_degree () const
{
  return degree;
};




template class FE_DGP<deal_II_dimension>;
