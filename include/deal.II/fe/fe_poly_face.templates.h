// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_poly_face_templates_h
#define dealii_fe_poly_face_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/polynomial_space.h>
#include <deal.II/base/qprojector.h>

#include <deal.II/fe/fe_poly_face.h>
#include <deal.II/fe/fe_values.h>


DEAL_II_NAMESPACE_OPEN

template <typename PolynomialType, int dim, int spacedim>
FE_PolyFace<PolynomialType, dim, spacedim>::FE_PolyFace(
  const PolynomialType         &poly_space,
  const FiniteElementData<dim> &fe_data,
  const std::vector<bool>      &restriction_is_additive_flags)
  : FiniteElement<dim, spacedim>(
      fe_data,
      restriction_is_additive_flags,
      std::vector<ComponentMask>(1, ComponentMask(1, true)))
  , poly_space(poly_space)
{
  AssertDimension(dim, PolynomialType::dimension + 1);
}


template <typename PolynomialType, int dim, int spacedim>
unsigned int
FE_PolyFace<PolynomialType, dim, spacedim>::get_degree() const
{
  return this->degree;
}


//---------------------------------------------------------------------------
// Auxiliary functions
//---------------------------------------------------------------------------


template <typename PolynomialType, int dim, int spacedim>
UpdateFlags
FE_PolyFace<PolynomialType, dim, spacedim>::requires_update_flags(
  const UpdateFlags flags) const
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
template <typename PolynomialType, int dim, int spacedim>
void
FE_PolyFace<PolynomialType, dim, spacedim>::fill_fe_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &,
  const CellSimilarity::Similarity,
  const Quadrature<dim> &,
  const Mapping<dim, spacedim> &,
  const typename Mapping<dim, spacedim>::InternalDataBase &,
  const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim> &,
  const typename FiniteElement<dim, spacedim>::InternalDataBase &,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    &) const
{
  // Do nothing, since we do not have values in the interior. Since
  // FEValues initializes the output variables for this function
  // with invalid values, this means that we simply leave them at
  // the invalid values -- typically, signaling NaNs. This means
  // that when you later look at those components of shape
  // functions or solution vectors that correspond to the
  // face element, you will see signaling_NaNs. This simply means
  // that you should not use them -- the shape functions and the
  // solution do not actually live inside the cell, and so any
  // attempt at evaluating it there *should* yield an invalid
  // result.
}



template <typename PolynomialType, int dim, int spacedim>
void
FE_PolyFace<PolynomialType, dim, spacedim>::fill_fe_face_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &,
  const unsigned int              face_no,
  const hp::QCollection<dim - 1> &quadrature,
  const Mapping<dim, spacedim> &,
  const typename Mapping<dim, spacedim>::InternalDataBase &,
  const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim> &,
  const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    &output_data) const
{
  // convert data object to internal
  // data for this class. fails with
  // an exception if that is not
  // possible
  Assert(dynamic_cast<const InternalData *>(&fe_internal) != nullptr,
         ExcInternalError());
  const InternalData &fe_data = static_cast<const InternalData &>(fe_internal);

  AssertDimension(quadrature.size(), 1);

  if (fe_data.update_each & update_values)
    for (unsigned int i = 0; i < quadrature[0].size(); ++i)
      {
        for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
          output_data.shape_values(k, i) = 0.;
        switch (dim)
          {
            case 3:
              {
                // Fill data for quad shape functions
                if (this->n_dofs_per_quad(face_no) != 0)
                  {
                    const unsigned int foffset =
                      this->get_first_quad_index(face_no);
                    for (unsigned int k = 0; k < this->n_dofs_per_quad(face_no);
                         ++k)
                      output_data.shape_values(foffset + k, i) =
                        fe_data
                          .shape_values[k + this->get_first_face_quad_index(
                                              face_no)][i];
                  }
              }
              DEAL_II_FALLTHROUGH;

            case 2:
              {
                // Fill data for line shape functions
                if (this->n_dofs_per_line() != 0)
                  {
                    const unsigned int foffset = this->get_first_line_index();
                    for (unsigned int line = 0;
                         line < GeometryInfo<dim>::lines_per_face;
                         ++line)
                      {
                        for (unsigned int k = 0; k < this->n_dofs_per_line();
                             ++k)
                          output_data.shape_values(
                            foffset +
                              GeometryInfo<dim>::face_to_cell_lines(face_no,
                                                                    line) *
                                this->n_dofs_per_line() +
                              k,
                            i) =
                            fe_data.shape_values
                              [k + (line * this->n_dofs_per_line()) +
                               this->get_first_face_line_index(face_no)][i];
                      }
                  }
              }
              DEAL_II_FALLTHROUGH;

            case 1:
              {
                // Fill data for vertex shape functions
                if (this->n_dofs_per_vertex() != 0)
                  for (unsigned int lvertex = 0;
                       lvertex < GeometryInfo<dim>::vertices_per_face;
                       ++lvertex)
                    output_data.shape_values(
                      GeometryInfo<dim>::face_to_cell_vertices(face_no,
                                                               lvertex),
                      i) = fe_data.shape_values[lvertex][i];
                break;
              }
          }
      }
}


template <typename PolynomialType, int dim, int spacedim>
void
FE_PolyFace<PolynomialType, dim, spacedim>::fill_fe_subface_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &,
  const unsigned int         face_no,
  const unsigned int         sub_no,
  const Quadrature<dim - 1> &quadrature,
  const Mapping<dim, spacedim> &,
  const typename Mapping<dim, spacedim>::InternalDataBase &,
  const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim> &,
  const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    &output_data) const
{
  // convert data object to internal
  // data for this class. fails with
  // an exception if that is not
  // possible
  Assert(dynamic_cast<const InternalData *>(&fe_internal) != nullptr,
         ExcInternalError());
  const InternalData &fe_data = static_cast<const InternalData &>(fe_internal);

  const unsigned int foffset = fe_data.shape_values.size() * face_no;
  const unsigned int offset  = sub_no * quadrature.size();

  if (fe_data.update_each & update_values)
    {
      for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
        for (unsigned int i = 0; i < quadrature.size(); ++i)
          output_data.shape_values(k, i) = 0.;
      for (unsigned int k = 0; k < fe_data.shape_values.size(); ++k)
        for (unsigned int i = 0; i < quadrature.size(); ++i)
          output_data.shape_values(foffset + k, i) =
            fe_data.shape_values[k][i + offset];
    }

  Assert(!(fe_data.update_each & update_gradients), ExcNotImplemented());
  Assert(!(fe_data.update_each & update_hessians), ExcNotImplemented());
}

DEAL_II_NAMESPACE_CLOSE

#endif
