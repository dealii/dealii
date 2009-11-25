//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009 by the deal.II authors
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
#include <fe/fe_poly_face.h>


DEAL_II_NAMESPACE_OPEN

template <class POLY, int dim, int spacedim>
FE_PolyFace<POLY,dim,spacedim>::FE_PolyFace (
  const POLY& poly_space,
  const FiniteElementData<dim> &fe_data,
  const std::vector<bool> &restriction_is_additive_flags):
		FiniteElement<dim,spacedim> (fe_data,
					     restriction_is_additive_flags,
					     std::vector<std::vector<bool> > (1, std::vector<bool>(1,true))),
		poly_space(poly_space)
{
  AssertDimension(dim, POLY::dimension+1);
}


template <class POLY, int dim, int spacedim>
unsigned int
FE_PolyFace<POLY,dim,spacedim>::get_degree () const
{
  return this->degree;
}


//---------------------------------------------------------------------------
// Auxiliary functions
//---------------------------------------------------------------------------




template <class POLY, int dim, int spacedim>
UpdateFlags
FE_PolyFace<POLY,dim,spacedim>::update_once (const UpdateFlags flags) const
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
FE_PolyFace<POLY,dim,spacedim>::update_each (const UpdateFlags flags) const
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
FE_PolyFace<POLY,dim,spacedim>::get_data (
  const UpdateFlags,
  const Mapping<dim,spacedim>&,
  const Quadrature<dim>&) const
{
  Assert(false, ExcNotImplemented());
  return 0;
}


template <class POLY, int dim, int spacedim>
typename Mapping<dim,spacedim>::InternalDataBase *
FE_PolyFace<POLY,dim,spacedim>::get_face_data (
  const UpdateFlags update_flags,
  const Mapping<dim,spacedim>&,
  const Quadrature<dim-1>& quadrature) const
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
  std::vector<Tensor<1,dim-1> > grads(0);
  std::vector<Tensor<2,dim-1> > grad_grads(0);

				   // initialize fields only if really
				   // necessary. otherwise, don't
				   // allocate memory
  if (flags & update_values)
    {
      values.resize (this->dofs_per_cell);
      data->shape_values.resize (this->dofs_per_cell,
				 std::vector<double> (n_q_points));
      for (unsigned int i=0; i<n_q_points; ++i)
	{
	  poly_space.compute(quadrature.point(i),
			     values, grads, grad_grads);
	  
	  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
	    data->shape_values[k][i] = values[k];
	}
    }
				   // No derivatives of this element
				   // are implemented.
  if (flags & update_gradients || flags & update_hessians)
    {
      Assert(false, ExcNotImplemented());
    }
  
  return data;
}



template <class POLY, int dim, int spacedim>
typename Mapping<dim,spacedim>::InternalDataBase *
FE_PolyFace<POLY,dim,spacedim>::get_subface_data (
  const UpdateFlags flags,
  const Mapping<dim,spacedim>& mapping,
  const Quadrature<dim-1>& quadrature) const
{
  return get_face_data (flags, mapping,
			QProjector<dim-1>::project_to_all_children(quadrature));
}





//---------------------------------------------------------------------------
// Fill data of FEValues
//---------------------------------------------------------------------------
template <class POLY, int dim, int spacedim>
void
FE_PolyFace<POLY,dim,spacedim>::fill_fe_values 
  (const Mapping<dim,spacedim>&,
   const typename Triangulation<dim,spacedim>::cell_iterator&,
   const Quadrature<dim>&,
   typename Mapping<dim,spacedim>::InternalDataBase&,
   typename Mapping<dim,spacedim>::InternalDataBase&,
   FEValuesData<dim,spacedim>&,
   CellSimilarity::Similarity&) const
{
  Assert(false, ExcNotImplemented());
}



template <class POLY, int dim, int spacedim>
void
FE_PolyFace<POLY,dim,spacedim>::fill_fe_face_values (
  const Mapping<dim,spacedim>&,
  const typename Triangulation<dim,spacedim>::cell_iterator&,
  const unsigned int,
  const Quadrature<dim-1>& quadrature,
  typename Mapping<dim,spacedim>::InternalDataBase&,
  typename Mapping<dim,spacedim>::InternalDataBase& fedata,
  FEValuesData<dim,spacedim>& data) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  Assert (dynamic_cast<InternalData *> (&fedata) != 0, ExcInternalError());
  InternalData &fe_data = static_cast<InternalData &> (fedata);
  
  const UpdateFlags flags(fe_data.update_once | fe_data.update_each);
  
  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
    {
      if (flags & update_values)
        for (unsigned int i=0; i<quadrature.size(); ++i)
	  data.shape_values(k,i) = fe_data.shape_values[k][i];
    }
}


template <class POLY, int dim, int spacedim>
void
FE_PolyFace<POLY,dim,spacedim>::fill_fe_subface_values (
  const Mapping<dim,spacedim>&,
  const typename Triangulation<dim,spacedim>::cell_iterator&,
  const unsigned int,
  const unsigned int subface,
  const Quadrature<dim-1>& quadrature,
  typename Mapping<dim,spacedim>::InternalDataBase&,
  typename Mapping<dim,spacedim>::InternalDataBase& fedata,
  FEValuesData<dim,spacedim>& data) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  Assert (dynamic_cast<InternalData *> (&fedata) != 0, ExcInternalError());
  InternalData &fe_data = static_cast<InternalData &> (fedata);

  const UpdateFlags flags(fe_data.update_once | fe_data.update_each);

  unsigned int offset = subface*quadrature.size();
  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
    {
      if (flags & update_values)
        for (unsigned int i=0; i<quadrature.size(); ++i)
	  data.shape_values(k,i) = fe_data.shape_values[k][i+offset];
      
    }
  Assert (!(flags & update_gradients), ExcNotImplemented());
  Assert (!(flags & update_hessians), ExcNotImplemented());
}



template <class POLY, int dim, int spacedim>
unsigned int
FE_PolyFace<POLY,dim,spacedim>::n_base_elements () const
{
  return 1;
}



template <class POLY, int dim, int spacedim>
const FiniteElement<dim,spacedim> &
FE_PolyFace<POLY,dim,spacedim>::base_element (const unsigned int index) const
{
  Assert (index==0, ExcIndexRange(index, 0, 1));
  return *this;
}



template <class POLY, int dim, int spacedim>
unsigned int
FE_PolyFace<POLY,dim,spacedim>::element_multiplicity (const unsigned int index) const
{
  Assert (index==0, ExcIndexRange(index, 0, 1));
  return 1;
}


DEAL_II_NAMESPACE_CLOSE
