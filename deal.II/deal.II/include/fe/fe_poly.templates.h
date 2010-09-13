//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2006, 2007, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <base/qprojector.h>
#include <base/polynomial_space.h>
#include <fe/fe_values.h>
#include <fe/fe_poly.h>


DEAL_II_NAMESPACE_OPEN

template <class POLY, int dim, int spacedim>
FE_Poly<POLY,dim,spacedim>::FE_Poly (const POLY& poly_space,
			    const FiniteElementData<dim> &fe_data,
			    const std::vector<bool> &restriction_is_additive_flags,
			    const std::vector<std::vector<bool> > &nonzero_components):
		FiniteElement<dim,spacedim> (fe_data,
				    restriction_is_additive_flags,
				    nonzero_components),
                poly_space(poly_space)
{
  AssertDimension(dim, POLY::dimension);
}


template <class POLY, int dim, int spacedim>
unsigned int
FE_Poly<POLY,dim,spacedim>::get_degree () const
{
  return this->degree;
}


template <class POLY, int dim, int spacedim>
double
FE_Poly<POLY,dim,spacedim>::shape_value (const unsigned int i,
				const Point<dim> &p) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  return poly_space.compute_value(i, p);
}


template <class POLY, int dim, int spacedim>
double
FE_Poly<POLY,dim,spacedim>::shape_value_component (const unsigned int i,
					   const Point<dim> &p,
					   const unsigned int component) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  Assert (component == 0, ExcIndexRange (component, 0, 1));
  return poly_space.compute_value(i, p);
}



template <class POLY, int dim, int spacedim>
Tensor<1,dim>
FE_Poly<POLY,dim,spacedim>::shape_grad (const unsigned int i,
				const Point<dim> &p) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  return poly_space.compute_grad(i, p);
}



template <class POLY, int dim, int spacedim>
Tensor<1,dim>
FE_Poly<POLY,dim,spacedim>::shape_grad_component (const unsigned int i,
					  const Point<dim> &p,
					  const unsigned int component) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  Assert (component == 0, ExcIndexRange (component, 0, 1));
  return poly_space.compute_grad(i, p);
}



template <class POLY, int dim, int spacedim>
Tensor<2,dim>
FE_Poly<POLY,dim,spacedim>::shape_grad_grad (const unsigned int i,
				     const Point<dim> &p) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  return poly_space.compute_grad_grad(i, p);
}



template <class POLY, int dim, int spacedim>
Tensor<2,dim>
FE_Poly<POLY,dim,spacedim>::shape_grad_grad_component (const unsigned int i,
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




template <class POLY, int dim, int spacedim>
UpdateFlags
FE_Poly<POLY,dim,spacedim>::update_once (const UpdateFlags flags) const
{
				   // for this kind of elements, only
				   // the values can be precomputed
				   // once and for all. set this flag
				   // if the values are requested at
				   // all
  return (update_default | (flags & update_values));
}



template <class POLY, int dim, int spacedim>
UpdateFlags
FE_Poly<POLY,dim,spacedim>::update_each (const UpdateFlags flags) const
{
  UpdateFlags out = update_default;

  if (flags & update_gradients)
    out |= update_gradients | update_covariant_transformation;
  if (flags & update_hessians)
    out |= update_hessians | update_covariant_transformation;
  if (flags & update_cell_normal_vectors)
    out |= update_cell_normal_vectors | update_JxW_values;
  
  return out;
}



//---------------------------------------------------------------------------
// Data field initialization
//---------------------------------------------------------------------------

template <class POLY, int dim, int spacedim>
typename Mapping<dim,spacedim>::InternalDataBase *
FE_Poly<POLY,dim,spacedim>::get_data (const UpdateFlags      update_flags,
				      const Mapping<dim,spacedim>    &mapping,
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
  const unsigned int n_q_points = quadrature.size();

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
  if (flags & update_hessians)
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


template <>
inline
void
FE_Poly<TensorProductPolynomials<1>,1,2>::fill_fe_values 
  (const Mapping<1,2>                        &mapping,
   const Triangulation<1,2>::cell_iterator   &/*cell*/,
   const Quadrature<1>                       &quadrature,
   Mapping<1,2>::InternalDataBase            &mapping_data,
   Mapping<1,2>::InternalDataBase            &fedata,
   FEValuesData<1,2>                         &data,
   CellSimilarity::Similarity           &cell_similarity) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  Assert (dynamic_cast<InternalData *> (&fedata) != 0, ExcInternalError());
  InternalData &fe_data = static_cast<InternalData &> (fedata);
  
  const UpdateFlags flags(fe_data.current_update_flags());
  
  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
    {
      if (flags & update_values)
	for (unsigned int i=0; i<quadrature.size(); ++i)
	  data.shape_values(k,i) = fe_data.shape_values[k][i];
      
      if (flags & update_gradients && cell_similarity != CellSimilarity::translation)
	mapping.transform(fe_data.shape_gradients[k], data.shape_gradients[k], 
			  mapping_data, mapping_covariant);
    }
  
  //  const typename QProjector<1>::DataSetDescriptor dsd;
  if (flags & update_hessians && cell_similarity != CellSimilarity::translation)
    {
      AssertThrow(false, ExcNotImplemented());
/*       this->compute_2nd (mapping, cell, dsd.cell(), */
/* 			 mapping_data, fe_data, data); */
    }
}



template <>
inline
void
FE_Poly<TensorProductPolynomials<2>,2,3>::fill_fe_values 
  (const Mapping<2,3>                      &mapping,
   const Triangulation<2,3>::cell_iterator &/* cell */,
   const Quadrature<2>                     &quadrature,
   Mapping<2,3>::InternalDataBase          &mapping_data,
   Mapping<2,3>::InternalDataBase          &fedata,
   FEValuesData<2,3>                       &data,
   CellSimilarity::Similarity         &cell_similarity) const
{

				   // assert that the following dynamics
				   // cast is really well-defined.
  Assert (dynamic_cast<InternalData *> (&fedata) != 0, ExcInternalError());
  InternalData &fe_data = static_cast<InternalData &> (fedata);
  
  const UpdateFlags flags(fe_data.current_update_flags());
  
  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
    {
      if (flags & update_values)
	for (unsigned int i=0; i<quadrature.size(); ++i)
	  data.shape_values(k,i) = fe_data.shape_values[k][i];     

       if (flags & update_gradients && cell_similarity != CellSimilarity::translation)
	 mapping.transform(fe_data.shape_gradients[k], data.shape_gradients[k],
			   mapping_data, mapping_covariant);
    }
  
/*   const QProjector<2>::DataSetDescriptor dsd; */
  if (flags & update_hessians && cell_similarity != CellSimilarity::translation)
    {
      AssertThrow(false, ExcNotImplemented());
/*       this->compute_2nd (mapping, cell, dsd.cell(), */
/* 			 mapping_data, fe_data, data); */
    }

}


template <>
inline
void
FE_Poly<PolynomialSpace<1>,1,2>::fill_fe_values (
  const Mapping<1,2>                      &mapping,
  const Triangulation<1,2>::cell_iterator & /*cell*/,
  const Quadrature<1>                     &quadrature,
  Mapping<1,2>::InternalDataBase          &mapping_data,
  Mapping<1,2>::InternalDataBase          &fedata,
  FEValuesData<1,2>                       &data,
  CellSimilarity::Similarity         &cell_similarity) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  	 
  Assert (dynamic_cast<InternalData *> (&fedata) != 0, ExcInternalError());
  InternalData &fe_data = static_cast<InternalData &> (fedata);
  
  const UpdateFlags flags(fe_data.current_update_flags());
  
  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
    {
      if (flags & update_values)
	for (unsigned int i=0; i<quadrature.size(); ++i)
	  data.shape_values(k,i) = fe_data.shape_values[k][i];
      
// TODO: I would think this should work. Guido
      
      if (flags & update_gradients && cell_similarity != CellSimilarity::translation)
	mapping.transform(fe_data.shape_gradients[k], data.shape_gradients[k],
			  mapping_data, mapping_covariant);
    }
  
				   //  const typename QProjector<1>::DataSetDescriptor dsd;
  if (flags & update_hessians && cell_similarity != CellSimilarity::translation)
    {
      AssertThrow(false, ExcNotImplemented());
/*       this->compute_2nd (mapping, cell, dsd.cell(), */
/* 			 mapping_data, fe_data, data); */
    }
}


template <>
inline
void
FE_Poly<PolynomialSpace<2>,2,3>::fill_fe_values 
  (const Mapping<2,3>                      &mapping,
   const Triangulation<2,3>::cell_iterator &/* cell */,
   const Quadrature<2>                     &quadrature,
   Mapping<2,3>::InternalDataBase          &mapping_data,
   Mapping<2,3>::InternalDataBase          &fedata,
   FEValuesData<2,3>                       &data,
   CellSimilarity::Similarity         &cell_similarity) const
{
  Assert (dynamic_cast<InternalData *> (&fedata) != 0, ExcInternalError());
  InternalData &fe_data = static_cast<InternalData &> (fedata);
  
  const UpdateFlags flags(fe_data.current_update_flags());
  
  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
    {
      if (flags & update_values)
	for (unsigned int i=0; i<quadrature.size(); ++i)
	  data.shape_values(k,i) = fe_data.shape_values[k][i];
      

       if (flags & update_gradients && cell_similarity != CellSimilarity::translation)
	 mapping.transform(fe_data.shape_gradients[k], data.shape_gradients[k],
  				      mapping_data, mapping_covariant);
    }
  
/*   const QProjector<2>::DataSetDescriptor dsd; */
  if (flags & update_hessians && cell_similarity != CellSimilarity::translation)
    {
      AssertThrow(false, ExcNotImplemented());
/*       this->compute_2nd (mapping, cell, dsd.cell(), */
/* 			 mapping_data, fe_data, data); */
    }

}


template <class POLY, int dim, int spacedim>
void
FE_Poly<POLY,dim,spacedim>::fill_fe_values 
  (const Mapping<dim,spacedim>                      &mapping,
   const typename Triangulation<dim,spacedim>::cell_iterator &cell,
   const Quadrature<dim>                            &quadrature,
   typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
   typename Mapping<dim,spacedim>::InternalDataBase &fedata,
   FEValuesData<dim,spacedim>                       &data,
   CellSimilarity::Similarity                  &cell_similarity) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible	 
  Assert (dynamic_cast<InternalData *> (&fedata) != 0, ExcInternalError());
  InternalData &fe_data = static_cast<InternalData &> (fedata);
  
  const UpdateFlags flags(fe_data.current_update_flags());
  
  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
    {
      if (flags & update_values)
	for (unsigned int i=0; i<quadrature.size(); ++i)
	  data.shape_values(k,i) = fe_data.shape_values[k][i];
      
      if (flags & update_gradients && cell_similarity != CellSimilarity::translation)
	mapping.transform(fe_data.shape_gradients[k], data.shape_gradients[k],
			  mapping_data, mapping_covariant);
    }
  
  const typename QProjector<dim>::DataSetDescriptor dsd;
  if (flags & update_hessians && cell_similarity != CellSimilarity::translation)
    this->compute_2nd (mapping, cell, dsd.cell(),
		       mapping_data, fe_data, data);
}



template <>
inline
void
FE_Poly<TensorProductPolynomials<1>,1,2>::fill_fe_face_values 
  (const Mapping<1,2> &,
   const Triangulation<1,2>::cell_iterator &,
   const unsigned int,
   const Quadrature<0> &, 
   Mapping<1,2>::InternalDataBase &,
   Mapping<1,2>::InternalDataBase &,
   FEValuesData<1,2> &) const
{
  AssertThrow(false, ExcNotImplemented());
}  


 
template <>
inline
void
FE_Poly<TensorProductPolynomials<2>,2,3>::fill_fe_face_values (const Mapping<2,3> &,
							       const Triangulation<2,3>::cell_iterator &,
							       const unsigned int,
							       const Quadrature<1> &, 
							       Mapping<2,3>::InternalDataBase &,
							       Mapping<2,3>::InternalDataBase &,
							       FEValuesData<2,3> &) const
{
  AssertThrow(false, ExcNotImplemented());
}  



template <>
inline
void
FE_Poly<PolynomialSpace<1>,1,2>::fill_fe_face_values (const Mapping<1,2> &,
						      const Triangulation<1,2>::cell_iterator &,
						      const unsigned int,
						      const Quadrature<0> &, 
						      Mapping<1,2>::InternalDataBase &,
						      Mapping<1,2>::InternalDataBase &,
						      FEValuesData<1,2> &) const
{
  AssertThrow(false, ExcNotImplemented());
}  


 
template <>
inline
void
FE_Poly<PolynomialSpace<2>,2,3>::fill_fe_face_values (const Mapping<2,3> &,
						      const Triangulation<2,3>::cell_iterator &,
						      const unsigned int,
						      const Quadrature<1> &, 
						      Mapping<2,3>::InternalDataBase &,
						      Mapping<2,3>::InternalDataBase &,
						      FEValuesData<2,3> &) const
{
  AssertThrow(false, ExcNotImplemented());
}  


template <class POLY, int dim, int spacedim>
void
FE_Poly<POLY,dim,spacedim>::
fill_fe_face_values (const Mapping<dim,spacedim>                   &mapping,
		     const typename Triangulation<dim,spacedim>::cell_iterator &cell,
		     const unsigned int                    face,
		     const Quadrature<dim-1>              &quadrature,
		     typename Mapping<dim,spacedim>::InternalDataBase       &mapping_data,
		     typename Mapping<dim,spacedim>::InternalDataBase       &fedata,
		     FEValuesData<dim,spacedim>                    &data) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  Assert (dynamic_cast<InternalData *> (&fedata) != 0, ExcInternalError());
  InternalData &fe_data = static_cast<InternalData &> (fedata);

				   // offset determines which data set
				   // to take (all data sets for all
				   // faces are stored contiguously)
  
  const typename QProjector<dim>::DataSetDescriptor dsd;
  const typename QProjector<dim>::DataSetDescriptor offset
    = dsd.face (face,
		cell->face_orientation(face),
		cell->face_flip(face),
		cell->face_rotation(face),
		quadrature.size());
  
  const UpdateFlags flags(fe_data.update_once | fe_data.update_each);

  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
    {
      if (flags & update_values)
        for (unsigned int i=0; i<quadrature.size(); ++i)
	  data.shape_values(k,i) = fe_data.shape_values[k][i+offset];
      
      if (flags & update_gradients)
	mapping.transform(make_slice(fe_data.shape_gradients[k], offset, quadrature.size()),
			  data.shape_gradients[k],
			  mapping_data, mapping_covariant);
    }

  if (flags & update_hessians)
    this->compute_2nd (mapping, cell, offset, mapping_data, fe_data, data);
}


//codimension 1
template <>
inline
void
FE_Poly<TensorProductPolynomials<1>,1,2>::fill_fe_subface_values (const Mapping<1,2> &,
							       const Triangulation<1,2>::cell_iterator &,
							       const unsigned int,
							       const unsigned int,	  
							       const Quadrature<0> &, 
							       Mapping<1,2>::InternalDataBase &,
							       Mapping<1,2>::InternalDataBase &,
							       FEValuesData<1,2> &) const
{
  AssertThrow(false, ExcNotImplemented());
}  


template <>
inline
void
FE_Poly<TensorProductPolynomials<2>,2,3>::fill_fe_subface_values (const Mapping<2,3> &,
							       const Triangulation<2,3>::cell_iterator &,
							       const unsigned int,
							       const unsigned int,	  
							       const Quadrature<1> &, 
							       Mapping<2,3>::InternalDataBase &,
							       Mapping<2,3>::InternalDataBase &,
							       FEValuesData<2,3> &) const
{
  AssertThrow(false, ExcNotImplemented());
}  


template <>
inline
void
FE_Poly<PolynomialSpace<1>,1,2>::fill_fe_subface_values (const Mapping<1,2> &,
							       const Triangulation<1,2>::cell_iterator &,
							       const unsigned int,
							       const unsigned int,	  
							       const Quadrature<0> &, 
							       Mapping<1,2>::InternalDataBase &,
							       Mapping<1,2>::InternalDataBase &,
							       FEValuesData<1,2> &) const
{
  AssertThrow(false, ExcNotImplemented());
}  


template <>
inline
void
FE_Poly<PolynomialSpace<2>,2,3>::fill_fe_subface_values (const Mapping<2,3> &,
							       const Triangulation<2,3>::cell_iterator &,
							       const unsigned int,
							       const unsigned int,	  
							       const Quadrature<1> &, 
							       Mapping<2,3>::InternalDataBase &,
							       Mapping<2,3>::InternalDataBase &,
							       FEValuesData<2,3> &) const
{
  AssertThrow(false, ExcNotImplemented());
}  




template <class POLY, int dim, int spacedim>
void
FE_Poly<POLY,dim,spacedim>::fill_fe_subface_values (const Mapping<dim,spacedim>                   &mapping,
				   const typename Triangulation<dim,spacedim>::cell_iterator &cell,
				   const unsigned int                    face,
				   const unsigned int                    subface,
				   const Quadrature<dim-1>              &quadrature,
				   typename Mapping<dim,spacedim>::InternalDataBase       &mapping_data,
				   typename Mapping<dim,spacedim>::InternalDataBase       &fedata,
				   FEValuesData<dim,spacedim>                    &data) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  Assert (dynamic_cast<InternalData *> (&fedata) != 0, ExcInternalError());
  InternalData &fe_data = static_cast<InternalData &> (fedata);

				   // offset determines which data set
				   // to take (all data sets for all
				   // sub-faces are stored contiguously)

  const typename QProjector<dim>::DataSetDescriptor dsd;
  const typename QProjector<dim>::DataSetDescriptor offset
    = dsd.subface (face, subface,
		   cell->face_orientation(face),
		   cell->face_flip(face),
		   cell->face_rotation(face),
		   quadrature.size(),
		   cell->subface_case(face));

  const UpdateFlags flags(fe_data.update_once | fe_data.update_each);

  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
    {
      if (flags & update_values)
        for (unsigned int i=0; i<quadrature.size(); ++i)
	  data.shape_values(k,i) = fe_data.shape_values[k][i+offset];
      
      if (flags & update_gradients)
	mapping.transform(make_slice(fe_data.shape_gradients[k], offset, quadrature.size()),
			  data.shape_gradients[k],
			  mapping_data, mapping_covariant);
    }
  
  if (flags & update_hessians)
    this->compute_2nd (mapping, cell, offset, mapping_data, fe_data, data);
}



template <class POLY, int dim, int spacedim>
unsigned int
FE_Poly<POLY,dim,spacedim>::n_base_elements () const
{
  return 1;
}



template <class POLY, int dim, int spacedim>
const FiniteElement<dim,spacedim> &
FE_Poly<POLY,dim,spacedim>::base_element (const unsigned int index) const
{
  Assert (index==0, ExcIndexRange(index, 0, 1));
  return *this;
}



template <class POLY, int dim, int spacedim>
unsigned int
FE_Poly<POLY,dim,spacedim>::element_multiplicity (const unsigned int index) const
{
  Assert (index==0, ExcIndexRange(index, 0, 1));
  return 1;
}


DEAL_II_NAMESPACE_CLOSE
