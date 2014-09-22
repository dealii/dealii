// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include <deal.II/base/qprojector.h>
#include <deal.II/base/polynomial_space.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_poly_face.h>


DEAL_II_NAMESPACE_OPEN

template <class POLY, int dim, int spacedim>
FE_PolyFace<POLY,dim,spacedim>::FE_PolyFace (
  const POLY &poly_space,
  const FiniteElementData<dim> &fe_data,
  const std::vector<bool> &restriction_is_additive_flags):
  FiniteElement<dim,spacedim> (fe_data,
                               restriction_is_additive_flags,
                               std::vector<ComponentMask> (1, ComponentMask(1,true))),
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
FE_PolyFace<POLY,dim,spacedim>::update_once (const UpdateFlags) const
{
  // for this kind of elements, only the values can be precomputed once and
  // for all. set this flag if the values are requested at all
  return update_default;
}



template <class POLY, int dim, int spacedim>
UpdateFlags
FE_PolyFace<POLY,dim,spacedim>::update_each (const UpdateFlags flags) const
{
  UpdateFlags out = flags & update_values;
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
  const Mapping<dim,spacedim> &,
  const Quadrature<dim> &) const
{
  InternalData *data = new InternalData;
  return data;
}


template <class POLY, int dim, int spacedim>
typename Mapping<dim,spacedim>::InternalDataBase *
FE_PolyFace<POLY,dim,spacedim>::get_face_data (
  const UpdateFlags update_flags,
  const Mapping<dim,spacedim> &,
  const Quadrature<dim-1>& quadrature) const
{
  // generate a new data object and
  // initialize some fields
  InternalData *data = new InternalData;

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
      values.resize (poly_space.n());
      data->shape_values.resize (poly_space.n(),
                                 std::vector<double> (n_q_points));
      for (unsigned int i=0; i<n_q_points; ++i)
        {
          poly_space.compute(quadrature.point(i),
                             values, grads, grad_grads);

          for (unsigned int k=0; k<poly_space.n(); ++k)
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
  const Mapping<dim,spacedim> &mapping,
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
(const Mapping<dim,spacedim> &,
 const typename Triangulation<dim,spacedim>::cell_iterator &,
 const Quadrature<dim> &,
 typename Mapping<dim,spacedim>::InternalDataBase &,
 typename Mapping<dim,spacedim>::InternalDataBase &,
 FEValuesData<dim,spacedim> &,
 CellSimilarity::Similarity &) const
{
  // Do nothing, since we do not have
  // values in the interior
}



template <class POLY, int dim, int spacedim>
void
FE_PolyFace<POLY,dim,spacedim>::fill_fe_face_values (
  const Mapping<dim,spacedim> &,
  const typename Triangulation<dim,spacedim>::cell_iterator &,
  const unsigned int face,
  const Quadrature<dim-1>& quadrature,
  typename Mapping<dim,spacedim>::InternalDataBase &,
  typename Mapping<dim,spacedim>::InternalDataBase &fedata,
  FEValuesData<dim,spacedim> &data) const
{
  // convert data object to internal
  // data for this class. fails with
  // an exception if that is not
  // possible
  Assert (dynamic_cast<InternalData *> (&fedata) != 0, ExcInternalError());
  InternalData &fe_data = static_cast<InternalData &> (fedata);

  const UpdateFlags flags(fe_data.update_once | fe_data.update_each);

  if (flags & update_values)
    for (unsigned int i=0; i<quadrature.size(); ++i)
      {
        for (unsigned int k=0; k<this->dofs_per_cell; ++k)
          data.shape_values(k,i) = 0.;
        switch (dim)
          {
          case 3:
          {
            // Fill data for quad shape functions
            if (this->dofs_per_quad !=0)
              {
                const unsigned int foffset = this->first_quad_index + this->dofs_per_quad * face;
                for (unsigned int k=0; k<this->dofs_per_quad; ++k)
                  data.shape_values(foffset+k,i) = fe_data.shape_values[k+this->first_face_quad_index][i];
              }
          }
          case 2:
          {
            // Fill data for line shape functions
            if (this->dofs_per_line != 0)
              {
                const unsigned int foffset = this->first_line_index;
                for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_face; ++line)
                  {
                    for (unsigned int k=0; k<this->dofs_per_line; ++k)
                      data.shape_values(foffset+GeometryInfo<dim>::face_to_cell_lines(face, line)*this->dofs_per_line+k,i) = fe_data.shape_values[k+(line*this->dofs_per_line)+this->first_face_line_index][i];
                  }
              }
          }
          case 1:
          {
            // Fill data for vertex shape functions
            if (this->dofs_per_vertex != 0)
              for (unsigned int lvertex=0; lvertex<GeometryInfo<dim>::vertices_per_face; ++lvertex)
                data.shape_values(GeometryInfo<dim>::face_to_cell_vertices(face, lvertex),i) = fe_data.shape_values[lvertex][i];
            break;
          }
          }
      }
}


template <class POLY, int dim, int spacedim>
void
FE_PolyFace<POLY,dim,spacedim>::fill_fe_subface_values (
  const Mapping<dim,spacedim> &,
  const typename Triangulation<dim,spacedim>::cell_iterator &,
  const unsigned int face,
  const unsigned int subface,
  const Quadrature<dim-1>& quadrature,
  typename Mapping<dim,spacedim>::InternalDataBase &,
  typename Mapping<dim,spacedim>::InternalDataBase &fedata,
  FEValuesData<dim,spacedim> &data) const
{
  // convert data object to internal
  // data for this class. fails with
  // an exception if that is not
  // possible
  Assert (dynamic_cast<InternalData *> (&fedata) != 0, ExcInternalError());
  InternalData &fe_data = static_cast<InternalData &> (fedata);

  const UpdateFlags flags(fe_data.update_once | fe_data.update_each);

  const unsigned int foffset = fe_data.shape_values.size() * face;
  const unsigned int offset = subface*quadrature.size();

  if (flags & update_values)
    for (unsigned int i=0; i<quadrature.size(); ++i)
      {
        for (unsigned int k=0; k<this->dofs_per_cell; ++k)
          data.shape_values(k,i) = 0.;
        for (unsigned int k=0; k<fe_data.shape_values.size(); ++k)
          data.shape_values(foffset+k,i) = fe_data.shape_values[k][i+offset];
      }

  Assert (!(flags & update_gradients), ExcNotImplemented());
  Assert (!(flags & update_hessians), ExcNotImplemented());
}

DEAL_II_NAMESPACE_CLOSE
