// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2015 by the deal.II authors
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
FE_Poly<TensorProductPolynomials<1>,1,2>::
fill_fe_values (const Mapping<1,2>                                &mapping,
                const Triangulation<1,2>::cell_iterator &,
                const Quadrature<1>                               &quadrature,
                const Mapping<1,2>::InternalDataBase              &mapping_internal,
                const FiniteElement<1,2>::InternalDataBase        &fedata,
                const internal::FEValues::MappingRelatedData<1,2> &mapping_data,
                internal::FEValues::FiniteElementRelatedData<1,2> &output_data,
                const CellSimilarity::Similarity                  cell_similarity) const
{
  // convert data object to internal
  // data for this class. fails with
  // an exception if that is not
  // possible
  Assert (dynamic_cast<const InternalData *> (&fedata) != 0, ExcInternalError());
  const InternalData &fe_data = static_cast<const InternalData &> (fedata);

  const UpdateFlags flags(fe_data.current_update_flags());

  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
    {
      if (flags & update_values)
        for (unsigned int i=0; i<quadrature.size(); ++i)
          output_data.shape_values(k,i) = fe_data.shape_values[k][i];

      if (flags & update_gradients && cell_similarity != CellSimilarity::translation)
        mapping.transform(fe_data.shape_gradients[k],
                          output_data.shape_gradients[k],
                          mapping_internal, mapping_covariant);

      if (flags & update_hessians && cell_similarity != CellSimilarity::translation)
        {
          // compute the hessians in the unit cell (accounting for the Jacobian gradiant)
          for (unsigned int i=0; i<quadrature.size(); ++i)
            {
              fe_data.untransformed_shape_hessians[i] = fe_data.shape_hessians[k][i];
            }

          correct_untransformed_hessians (fe_data.untransformed_shape_hessians,
                                          mapping_data, output_data, quadrature.size(), k);

          mapping.transform(fe_data.untransformed_shape_hessians,
                            output_data.shape_hessians[k],
                            mapping_internal, mapping_covariant_gradient);
        }
    }
}



template <>
void
FE_Poly<TensorProductPolynomials<2>,2,3>::
fill_fe_values (const Mapping<2,3>                                &mapping,
                const Triangulation<2,3>::cell_iterator &,
                const Quadrature<2>                               &quadrature,
                const Mapping<2,3>::InternalDataBase              &mapping_internal,
                const FiniteElement<2,3>::InternalDataBase        &fedata,
                const internal::FEValues::MappingRelatedData<2,3> &mapping_data,
                internal::FEValues::FiniteElementRelatedData<2,3> &output_data,
                const CellSimilarity::Similarity                  cell_similarity) const
{

  // assert that the following dynamics
  // cast is really well-defined.
  Assert (dynamic_cast<const InternalData *> (&fedata) != 0, ExcInternalError());
  const InternalData &fe_data = static_cast<const InternalData &> (fedata);

  const UpdateFlags flags(fe_data.current_update_flags());

  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
    {
      if (flags & update_values)
        for (unsigned int i=0; i<quadrature.size(); ++i)
          output_data.shape_values(k,i) = fe_data.shape_values[k][i];

      if (flags & update_gradients && cell_similarity != CellSimilarity::translation)
        mapping.transform(fe_data.shape_gradients[k],
                          output_data.shape_gradients[k],
                          mapping_internal, mapping_covariant);

      if (flags & update_hessians && cell_similarity != CellSimilarity::translation)
        {
          // compute the hessians in the unit cell (accounting for the Jacobian gradiant)
          for (unsigned int i=0; i<quadrature.size(); ++i)
            {
              fe_data.untransformed_shape_hessians[i] = fe_data.shape_hessians[k][i];
            }

          correct_untransformed_hessians (fe_data.untransformed_shape_hessians,
                                          mapping_data, output_data, quadrature.size(), k);

          mapping.transform(fe_data.untransformed_shape_hessians,
                            output_data.shape_hessians[k],
                            mapping_internal, mapping_covariant_gradient);
        }
    }
}


template <>
void
FE_Poly<PolynomialSpace<1>,1,2>::
fill_fe_values (const Mapping<1,2>                                &mapping,
                const Triangulation<1,2>::cell_iterator &,
                const Quadrature<1>                               &quadrature,
                const Mapping<1,2>::InternalDataBase              &mapping_internal,
                const FiniteElement<1,2>::InternalDataBase        &fedata,
                const internal::FEValues::MappingRelatedData<1,2> &mapping_data,
                internal::FEValues::FiniteElementRelatedData<1,2> &output_data,
                const CellSimilarity::Similarity                  cell_similarity) const
{
  // convert data object to internal
  // data for this class. fails with
  // an exception if that is not
  // possible

  Assert (dynamic_cast<const InternalData *> (&fedata) != 0, ExcInternalError());
  const InternalData &fe_data = static_cast<const InternalData &> (fedata);

  const UpdateFlags flags(fe_data.current_update_flags());

  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
    {
      if (flags & update_values)
        for (unsigned int i=0; i<quadrature.size(); ++i)
          output_data.shape_values(k,i) = fe_data.shape_values[k][i];

      if (flags & update_gradients && cell_similarity != CellSimilarity::translation)
        mapping.transform(fe_data.shape_gradients[k],
                          output_data.shape_gradients[k],
                          mapping_internal, mapping_covariant);

      if (flags & update_hessians && cell_similarity != CellSimilarity::translation)
        {
          // compute the hessians in the unit cell (accounting for the Jacobian gradiant)
          for (unsigned int i=0; i<quadrature.size(); ++i)
            {
              fe_data.untransformed_shape_hessians[i] = fe_data.shape_hessians[k][i];
            }

          correct_untransformed_hessians (fe_data.untransformed_shape_hessians,
                                          mapping_data, output_data, quadrature.size(), k);

          mapping.transform(fe_data.untransformed_shape_hessians,
                            output_data.shape_hessians[k],
                            mapping_internal, mapping_covariant_gradient);
        }
    }
}


template <>
void
FE_Poly<PolynomialSpace<2>,2,3>::
fill_fe_values (const Mapping<2,3>                                &mapping,
                const Triangulation<2,3>::cell_iterator &,
                const Quadrature<2>                               &quadrature,
                const Mapping<2,3>::InternalDataBase              &mapping_internal,
                const FiniteElement<2,3>::InternalDataBase        &fedata,
                const internal::FEValues::MappingRelatedData<2,3> &mapping_data,
                internal::FEValues::FiniteElementRelatedData<2,3> &output_data,
                const CellSimilarity::Similarity                  cell_similarity) const
{
  Assert (dynamic_cast<const InternalData *> (&fedata) != 0, ExcInternalError());
  const InternalData &fe_data = static_cast<const InternalData &> (fedata);

  const UpdateFlags flags(fe_data.current_update_flags());

  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
    {
      if (flags & update_values)
        for (unsigned int i=0; i<quadrature.size(); ++i)
          output_data.shape_values(k,i) = fe_data.shape_values[k][i];


      if (flags & update_gradients && cell_similarity != CellSimilarity::translation)
        mapping.transform(fe_data.shape_gradients[k],
                          output_data.shape_gradients[k],
                          mapping_internal, mapping_covariant);

      if (flags & update_hessians && cell_similarity != CellSimilarity::translation)
        {
          // compute the hessians in the unit cell (accounting for the Jacobian gradiant)
          for (unsigned int i=0; i<quadrature.size(); ++i)
            {
              fe_data.untransformed_shape_hessians[i] = fe_data.shape_hessians[k][i];
            }

          correct_untransformed_hessians (fe_data.untransformed_shape_hessians,
                                          mapping_data, output_data, quadrature.size(), k);

          mapping.transform(fe_data.untransformed_shape_hessians,
                            output_data.shape_hessians[k],
                            mapping_internal, mapping_covariant_gradient);
        }
    }
}


#include "fe_poly.inst"

DEAL_II_NAMESPACE_CLOSE
