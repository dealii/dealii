// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2016 by the deal.II authors
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

#ifndef dealii__fe_poly_face_templates_h
#define dealii__fe_poly_face_templates_h


#include <deal.II/base/qprojector.h>
#include <deal.II/base/polynomial_space.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_poly_face.h>


DEAL_II_NAMESPACE_OPEN

template <class PolynomialType, int dim, int spacedim>
FE_PolyFace<PolynomialType,dim,spacedim>::FE_PolyFace (
  const PolynomialType &poly_space,
  const FiniteElementData<dim> &fe_data,
  const std::vector<bool> &restriction_is_additive_flags):
  FiniteElement<dim,spacedim> (fe_data,
                               restriction_is_additive_flags,
                               std::vector<ComponentMask> (1, ComponentMask(1,true))),
  poly_space(poly_space)
{
  AssertDimension(dim, PolynomialType::dimension+1);
}


template <class PolynomialType, int dim, int spacedim>
unsigned int
FE_PolyFace<PolynomialType,dim,spacedim>::get_degree () const
{
  return this->degree;
}


//---------------------------------------------------------------------------
// Auxiliary functions
//---------------------------------------------------------------------------


template <class PolynomialType, int dim, int spacedim>
UpdateFlags
FE_PolyFace<PolynomialType,dim,spacedim>::requires_update_flags (const UpdateFlags flags) const
{
  UpdateFlags out = flags & update_values;
  if (flags & update_gradients)
    out |= update_gradients | update_covariant_transformation;
  if (flags & update_hessians)
    out |= update_hessians | update_covariant_transformation;
  if (flags & update_normal_vectors)
    out |= update_normal_vectors | update_JxW_values;

  return out;
}


//---------------------------------------------------------------------------
// Fill data of FEValues
//---------------------------------------------------------------------------
template <class PolynomialType, int dim, int spacedim>
void
FE_PolyFace<PolynomialType,dim,spacedim>::
fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator &,
                const CellSimilarity::Similarity                                     ,
                const Quadrature<dim> &,
                const Mapping<dim,spacedim> &,
                const typename Mapping<dim,spacedim>::InternalDataBase &,
                const dealii::internal::FEValues::MappingRelatedData<dim, spacedim> &,
                const typename FiniteElement<dim,spacedim>::InternalDataBase &,
                dealii::internal::FEValues::FiniteElementRelatedData<dim, spacedim> &) const
{
  // Do nothing, since we do not have
  // values in the interior
}



template <class PolynomialType, int dim, int spacedim>
void
FE_PolyFace<PolynomialType,dim,spacedim>::
fill_fe_face_values (const typename Triangulation<dim,spacedim>::cell_iterator &,
                     const unsigned int                                                   face_no,
                     const Quadrature<dim-1>                                             &quadrature,
                     const Mapping<dim,spacedim> &,
                     const typename Mapping<dim,spacedim>::InternalDataBase &,
                     const dealii::internal::FEValues::MappingRelatedData<dim, spacedim> &,
                     const typename FiniteElement<dim,spacedim>::InternalDataBase        &fe_internal,
                     dealii::internal::FEValues::FiniteElementRelatedData<dim, spacedim> &output_data) const
{
  // convert data object to internal
  // data for this class. fails with
  // an exception if that is not
  // possible
  Assert (dynamic_cast<const InternalData *> (&fe_internal) != nullptr, ExcInternalError());
  const InternalData &fe_data = static_cast<const InternalData &> (fe_internal);

  if (fe_data.update_each & update_values)
    for (unsigned int i=0; i<quadrature.size(); ++i)
      {
        for (unsigned int k=0; k<this->dofs_per_cell; ++k)
          output_data.shape_values(k,i) = 0.;
        switch (dim)
          {
          case 3:
          {
            // Fill data for quad shape functions
            if (this->dofs_per_quad !=0)
              {
                const unsigned int foffset = this->first_quad_index + this->dofs_per_quad * face_no;
                for (unsigned int k=0; k<this->dofs_per_quad; ++k)
                  output_data.shape_values(foffset+k,i) = fe_data.shape_values[k+this->first_face_quad_index][i];
              }
          }
          DEAL_II_FALLTHROUGH;

          case 2:
          {
            // Fill data for line shape functions
            if (this->dofs_per_line != 0)
              {
                const unsigned int foffset = this->first_line_index;
                for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_face; ++line)
                  {
                    for (unsigned int k=0; k<this->dofs_per_line; ++k)
                      output_data.shape_values(foffset+GeometryInfo<dim>::face_to_cell_lines(face_no, line)*this->dofs_per_line+k,i)
                        = fe_data.shape_values[k+(line*this->dofs_per_line)+this->first_face_line_index][i];
                  }
              }
          }
          DEAL_II_FALLTHROUGH;

          case 1:
          {
            // Fill data for vertex shape functions
            if (this->dofs_per_vertex != 0)
              for (unsigned int lvertex=0; lvertex<GeometryInfo<dim>::vertices_per_face; ++lvertex)
                output_data.shape_values(GeometryInfo<dim>::face_to_cell_vertices(face_no, lvertex),i)
                  = fe_data.shape_values[lvertex][i];
            break;
          }
          }
      }
}


template <class PolynomialType, int dim, int spacedim>
void
FE_PolyFace<PolynomialType,dim,spacedim>::
fill_fe_subface_values (const typename Triangulation<dim,spacedim>::cell_iterator &,
                        const unsigned int                                                   face_no,
                        const unsigned int                                                   sub_no,
                        const Quadrature<dim-1>                                             &quadrature,
                        const Mapping<dim,spacedim> &,
                        const typename Mapping<dim,spacedim>::InternalDataBase &,
                        const dealii::internal::FEValues::MappingRelatedData<dim, spacedim> &,
                        const typename FiniteElement<dim,spacedim>::InternalDataBase        &fe_internal,
                        dealii::internal::FEValues::FiniteElementRelatedData<dim, spacedim> &output_data) const
{
  // convert data object to internal
  // data for this class. fails with
  // an exception if that is not
  // possible
  Assert (dynamic_cast<const InternalData *> (&fe_internal) != nullptr, ExcInternalError());
  const InternalData &fe_data = static_cast<const InternalData &> (fe_internal);

  const unsigned int foffset = fe_data.shape_values.size() * face_no;
  const unsigned int offset = sub_no*quadrature.size();

  if (fe_data.update_each & update_values)
    {
      for (unsigned int k=0; k<this->dofs_per_cell; ++k)
        for (unsigned int i=0; i<quadrature.size(); ++i)
          output_data.shape_values(k,i) = 0.;
      for (unsigned int k=0; k<fe_data.shape_values.size(); ++k)
        for (unsigned int i=0; i<quadrature.size(); ++i)
          output_data.shape_values(foffset+k,i) = fe_data.shape_values[k][i+offset];
    }

  Assert (!(fe_data.update_each & update_gradients), ExcNotImplemented());
  Assert (!(fe_data.update_each & update_hessians), ExcNotImplemented());
}

DEAL_II_NAMESPACE_CLOSE

#endif
