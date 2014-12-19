// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
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
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/tensor_product_polynomials_const.h>
#include <deal.II/base/polynomials_p.h>
#include <deal.II/base/polynomials_piecewise.h>
#include <deal.II/fe/fe_poly.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/fe/fe_poly.templates.h>

DEAL_II_NAMESPACE_OPEN


template <>
void
FE_Poly<TensorProductPolynomials<1>,1,2>::fill_fe_values
(const Mapping<1,2>                        &mapping,
 const Triangulation<1,2>::cell_iterator   &cell,
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

  if (flags & update_hessians && cell_similarity != CellSimilarity::translation)
    this->compute_2nd (mapping, cell, QProjector<1>::DataSetDescriptor::cell(),
                       mapping_data, fe_data, data);
}



template <>
void
FE_Poly<TensorProductPolynomials<2>,2,3>::fill_fe_values
(const Mapping<2,3>                      &mapping,
 const Triangulation<2,3>::cell_iterator &cell,
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

  if (flags & update_hessians && cell_similarity != CellSimilarity::translation)
    this->compute_2nd (mapping, cell, QProjector<2>::DataSetDescriptor::cell(),
                       mapping_data, fe_data, data);
}


template <>
void
FE_Poly<PolynomialSpace<1>,1,2>::fill_fe_values (
  const Mapping<1,2>                      &mapping,
  const Triangulation<1,2>::cell_iterator &cell,
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

      if (flags & update_gradients && cell_similarity != CellSimilarity::translation)
        mapping.transform(fe_data.shape_gradients[k], data.shape_gradients[k],
                          mapping_data, mapping_covariant);
    }

  if (flags & update_hessians && cell_similarity != CellSimilarity::translation)
    this->compute_2nd (mapping, cell, QProjector<1>::DataSetDescriptor::cell(),
                       mapping_data, fe_data, data);
}


template <>
void
FE_Poly<PolynomialSpace<2>,2,3>::fill_fe_values
(const Mapping<2,3>                      &mapping,
 const Triangulation<2,3>::cell_iterator &cell,
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

  if (flags & update_hessians && cell_similarity != CellSimilarity::translation)
    this->compute_2nd (mapping, cell, QProjector<2>::DataSetDescriptor::cell(),
                       mapping_data, fe_data, data);
}


#include "fe_poly.inst"

DEAL_II_NAMESPACE_CLOSE
