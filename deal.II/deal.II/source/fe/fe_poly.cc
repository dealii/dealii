//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <base/qprojector.h>
#include <base/tensor_product_polynomials.h>
#include <base/polynomials_p.h>
#include <fe/fe_poly.h>
#include <fe/fe_values.h>


template <class POLY, int dim>
FE_Poly<POLY,dim>::FE_Poly (const POLY& poly_space,
			    const FiniteElementData<dim> &fe_data,
			    const std::vector<bool> &restriction_is_additive_flags,
			    const std::vector<std::vector<bool> > &nonzero_components):
		FiniteElement<dim> (fe_data,
				    restriction_is_additive_flags,
				    nonzero_components),
                poly_space(poly_space)
{}


template <class POLY, int dim>
unsigned int
FE_Poly<POLY,dim>::get_degree () const
{
  return this->degree;
}


template <class POLY, int dim>
double
FE_Poly<POLY,dim>::shape_value (const unsigned int i,
				const Point<dim> &p) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  return poly_space.compute_value(i, p);
}


template <class POLY, int dim>
double
FE_Poly<POLY,dim>::shape_value_component (const unsigned int i,
					   const Point<dim> &p,
					   const unsigned int component) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  Assert (component == 0, ExcIndexRange (component, 0, 1));
  return poly_space.compute_value(i, p);
}



template <class POLY, int dim>
Tensor<1,dim>
FE_Poly<POLY,dim>::shape_grad (const unsigned int i,
				const Point<dim> &p) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  return poly_space.compute_grad(i, p);
}



template <class POLY, int dim>
Tensor<1,dim>
FE_Poly<POLY,dim>::shape_grad_component (const unsigned int i,
					  const Point<dim> &p,
					  const unsigned int component) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  Assert (component == 0, ExcIndexRange (component, 0, 1));
  return poly_space.compute_grad(i, p);
}



template <class POLY, int dim>
Tensor<2,dim>
FE_Poly<POLY,dim>::shape_grad_grad (const unsigned int i,
				     const Point<dim> &p) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  return poly_space.compute_grad_grad(i, p);
}



template <class POLY, int dim>
Tensor<2,dim>
FE_Poly<POLY,dim>::shape_grad_grad_component (const unsigned int i,
					       const Point<dim> &p,
					       const unsigned int component) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  Assert (component == 0, ExcIndexRange (component, 0, 1));
  return poly_space.compute_grad_grad(i, p);
}



//---------------------------------------------------------------------------
// Auxiliary functions
//---------------------------------------------------------------------------




template <class POLY, int dim>
UpdateFlags
FE_Poly<POLY,dim>::update_once (const UpdateFlags flags) const
{
				   // for this kind of elements, only
				   // the values can be precomputed
				   // once and for all. set this flag
				   // if the values are requested at
				   // all
  return (update_default | (flags & update_values));
}



template <class POLY, int dim>
UpdateFlags
FE_Poly<POLY,dim>::update_each (const UpdateFlags flags) const
{
  UpdateFlags out = update_default;

  if (flags & update_gradients)
    out |= update_gradients | update_covariant_transformation;
  if (flags & update_second_derivatives)
    out |= update_second_derivatives | update_covariant_transformation;

  return out;
}



//---------------------------------------------------------------------------
// Data field initialization
//---------------------------------------------------------------------------

template <class POLY, int dim>
typename Mapping<dim>::InternalDataBase *
FE_Poly<POLY,dim>::get_data (const UpdateFlags      update_flags,
			     const Mapping<dim>    &mapping,
			     const Quadrature<dim> &quadrature) const
{
				   // generate a new data object and
				   // initialize some fields
  InternalData* data = new InternalData;

				   // check what needs to be
				   // initialized only once and what
				   // on every cell/face/subface we
				   // visit
  data->update_once = update_once(update_flags);
  data->update_each = update_each(update_flags);
  data->update_flags = data->update_once | data->update_each;

  const UpdateFlags flags(data->update_flags);
  const unsigned int n_q_points = quadrature.n_quadrature_points;

				   // some scratch arrays
  std::vector<double> values(0);
  std::vector<Tensor<1,dim> > grads(0);
  std::vector<Tensor<2,dim> > grad_grads(0);

				   // initialize fields only if really
				   // necessary. otherwise, don't
				   // allocate memory
  if (flags & update_values)
    {
      values.resize (this->dofs_per_cell);
      data->shape_values.resize (this->dofs_per_cell,
				 std::vector<double> (n_q_points));
    }

  if (flags & update_gradients)
    {
      grads.resize (this->dofs_per_cell);
      data->shape_gradients.resize (this->dofs_per_cell,
				    std::vector<Tensor<1,dim> > (n_q_points));
    }

				   // if second derivatives through
				   // finite differencing is required,
				   // then initialize some objects for
				   // that
  if (flags & update_second_derivatives)
    data->initialize_2nd (this, mapping, quadrature);

				   // next already fill those fields
				   // of which we have information by
				   // now. note that the shape
				   // gradients are only those on the
				   // unit cell, and need to be
				   // transformed when visiting an
				   // actual cell
  if (flags & (update_values | update_gradients))
    for (unsigned int i=0; i<n_q_points; ++i)
      {
	poly_space.compute(quadrature.point(i),
			   values, grads, grad_grads);
	
	if (flags & update_values)
	  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
	    data->shape_values[k][i] = values[k];
	
	if (flags & update_gradients)
	  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
	    data->shape_gradients[k][i] = grads[k];
      }
  return data;
}




//---------------------------------------------------------------------------
// Fill data of FEValues
//---------------------------------------------------------------------------

template <class POLY, int dim>
void
FE_Poly<POLY,dim>::fill_fe_values (const Mapping<dim>                   &mapping,
				   const typename Triangulation<dim>::cell_iterator &cell,
				   const Quadrature<dim>                &quadrature,
				   typename Mapping<dim>::InternalDataBase &mapping_data,
				   typename Mapping<dim>::InternalDataBase &fedata,
				   FEValuesData<dim>                    &data) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  InternalData &fe_data = dynamic_cast<InternalData &> (fedata);
  
  const UpdateFlags flags(fe_data.current_update_flags());

  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
    {
      if (flags & update_values)
	for (unsigned int i=0; i<quadrature.n_quadrature_points; ++i)
	  data.shape_values(k,i) = fe_data.shape_values[k][i];
      
      if (flags & update_gradients)
	{
	  Assert (data.shape_gradients[k].size() <=
		  fe_data.shape_gradients[k].size(),
		  ExcInternalError());	  
	  mapping.transform_covariant(fe_data.shape_gradients[k], 0,
                                      data.shape_gradients[k],
				      mapping_data);
	}
    }

  const typename QProjector<dim>::DataSetDescriptor dsd;
  if (flags & update_second_derivatives)
    this->compute_2nd (mapping, cell, dsd.cell(),
		       mapping_data, fe_data, data);
}



template <class POLY, int dim>
void
FE_Poly<POLY,dim>::fill_fe_face_values (const Mapping<dim>                   &mapping,
				const typename Triangulation<dim>::cell_iterator &cell,
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
  
  const typename QProjector<dim>::DataSetDescriptor dsd;
  const typename QProjector<dim>::DataSetDescriptor offset
    = dsd.face (face, cell->face_orientation(face),
		quadrature.n_quadrature_points);
  
  const UpdateFlags flags(fe_data.update_once | fe_data.update_each);

  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
    {
      if (flags & update_values)
        for (unsigned int i=0; i<quadrature.n_quadrature_points; ++i)
	  data.shape_values(k,i) = fe_data.shape_values[k][i+offset];
      
      if (flags & update_gradients)
	{
	  Assert (data.shape_gradients[k].size() + offset <=
		  fe_data.shape_gradients[k].size(),
		  ExcInternalError());
	  mapping.transform_covariant(fe_data.shape_gradients[k], offset,
                                      data.shape_gradients[k],
				      mapping_data);
	}
    }

  if (flags & update_second_derivatives)
    this->compute_2nd (mapping, cell, offset, mapping_data, fe_data, data);
}



template <class POLY, int dim>
void
FE_Poly<POLY,dim>::fill_fe_subface_values (const Mapping<dim>                   &mapping,
				   const typename Triangulation<dim>::cell_iterator &cell,
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

  const typename QProjector<dim>::DataSetDescriptor dsd;
  const typename QProjector<dim>::DataSetDescriptor offset
    = dsd.sub_face (face, subface, cell->face_orientation(face),
		    quadrature.n_quadrature_points);

  const UpdateFlags flags(fe_data.update_once | fe_data.update_each);

  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
    {
      if (flags & update_values)
        for (unsigned int i=0; i<quadrature.n_quadrature_points; ++i)
	  data.shape_values(k,i) = fe_data.shape_values[k][i+offset];
      
      if (flags & update_gradients)
	{
	  Assert (data.shape_gradients[k].size() + offset <=
		  fe_data.shape_gradients[k].size(),
		  ExcInternalError());
	  mapping.transform_covariant(fe_data.shape_gradients[k], offset,
                                      data.shape_gradients[k],
				      mapping_data);
	}
    }
  
  if (flags & update_second_derivatives)
    this->compute_2nd (mapping, cell, offset, mapping_data, fe_data, data);
}



template <class POLY, int dim>
unsigned int
FE_Poly<POLY,dim>::n_base_elements () const
{
  return 1;
}



template <class POLY, int dim>
const FiniteElement<dim> &
FE_Poly<POLY,dim>::base_element (const unsigned int index) const
{
  Assert (index==0, ExcIndexRange(index, 0, 1));
  return *this;
}



template <class POLY, int dim>
unsigned int
FE_Poly<POLY,dim>::element_multiplicity (const unsigned int index) const
{
  Assert (index==0, ExcIndexRange(index, 0, 1));
  return 1;
}



template class FE_Poly<TensorProductPolynomials<deal_II_dimension>, deal_II_dimension>;
template class FE_Poly<PolynomialSpace<deal_II_dimension>, deal_II_dimension>;
template class FE_Poly<PolynomialsP<deal_II_dimension>, deal_II_dimension>;
