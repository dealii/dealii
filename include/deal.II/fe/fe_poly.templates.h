// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2015 by the deal.II authors
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
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/tensor_product_polynomials_const.h>
#include <deal.II/base/tensor_product_polynomials_bubbles.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_poly.h>


DEAL_II_NAMESPACE_OPEN

template <class POLY, int dim, int spacedim>
FE_Poly<POLY,dim,spacedim>::FE_Poly (const POLY &poly_space,
                                     const FiniteElementData<dim> &fe_data,
                                     const std::vector<bool> &restriction_is_additive_flags,
                                     const std::vector<ComponentMask> &nonzero_components):
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
  (void)component;
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
  return poly_space.template compute_derivative<1>(i, p);
}



template <class POLY, int dim, int spacedim>
Tensor<1,dim>
FE_Poly<POLY,dim,spacedim>::shape_grad_component (const unsigned int i,
                                                  const Point<dim> &p,
                                                  const unsigned int component) const
{
  (void)component;
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  Assert (component == 0, ExcIndexRange (component, 0, 1));
  return poly_space.template compute_derivative<1>(i, p);
}



template <class POLY, int dim, int spacedim>
Tensor<2,dim>
FE_Poly<POLY,dim,spacedim>::shape_grad_grad (const unsigned int i,
                                             const Point<dim> &p) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  return poly_space.template compute_derivative<2>(i, p);
}



template <class POLY, int dim, int spacedim>
Tensor<2,dim>
FE_Poly<POLY,dim,spacedim>::shape_grad_grad_component (const unsigned int i,
                                                       const Point<dim> &p,
                                                       const unsigned int component) const
{
  (void)component;
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  Assert (component == 0, ExcIndexRange (component, 0, 1));
  return poly_space.template compute_derivative<2>(i, p);
}



template <class POLY, int dim, int spacedim>
Tensor<3,dim>
FE_Poly<POLY,dim,spacedim>::shape_3rd_derivative (const unsigned int i,
                                                  const Point<dim> &p) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  return poly_space.template compute_derivative<3>(i, p);
}



template <class POLY, int dim, int spacedim>
Tensor<3,dim>
FE_Poly<POLY,dim,spacedim>::shape_3rd_derivative_component (const unsigned int i,
                                                            const Point<dim> &p,
                                                            const unsigned int component) const
{
  (void)component;
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  Assert (component == 0, ExcIndexRange (component, 0, 1));
  return poly_space.template compute_derivative<3>(i, p);
}



template <class POLY, int dim, int spacedim>
Tensor<4,dim>
FE_Poly<POLY,dim,spacedim>::shape_4th_derivative (const unsigned int i,
                                                  const Point<dim> &p) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  return poly_space.template compute_derivative<4>(i, p);
}



template <class POLY, int dim, int spacedim>
Tensor<4,dim>
FE_Poly<POLY,dim,spacedim>::shape_4th_derivative_component (const unsigned int i,
                                                            const Point<dim> &p,
                                                            const unsigned int component) const
{
  (void)component;
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  Assert (component == 0, ExcIndexRange (component, 0, 1));
  return poly_space.template compute_derivative<4>(i, p);
}



//---------------------------------------------------------------------------
// Auxiliary functions
//---------------------------------------------------------------------------


template <class POLY, int dim, int spacedim>
UpdateFlags
FE_Poly<POLY,dim,spacedim>::requires_update_flags (const UpdateFlags flags) const
{
  return update_once(flags) | update_each(flags);
}



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
    out |= update_hessians | update_covariant_transformation
           | update_gradients | update_jacobian_pushed_forward_grads;
  if (flags & update_3rd_derivatives)
    out |= update_3rd_derivatives | update_covariant_transformation
           | update_hessians | update_gradients
           | update_jacobian_pushed_forward_grads
           | update_jacobian_pushed_forward_2nd_derivatives;
  if (flags & update_cell_normal_vectors)
    out |= update_cell_normal_vectors | update_JxW_values;

  return out;
}



//---------------------------------------------------------------------------
// Fill data of FEValues
//---------------------------------------------------------------------------


template <class POLY, int dim, int spacedim>
void
FE_Poly<POLY,dim,spacedim>::
fill_fe_values (const Mapping<dim,spacedim>                                  &mapping,
                const typename Triangulation<dim,spacedim>::cell_iterator &,
                const Quadrature<dim>                                        &quadrature,
                const typename Mapping<dim,spacedim>::InternalDataBase       &mapping_internal,
                const typename FiniteElement<dim,spacedim>::InternalDataBase &fedata,
                const internal::FEValues::MappingRelatedData<dim,spacedim>   &mapping_data,
                internal::FEValues::FiniteElementRelatedData<dim,spacedim>   &output_data,
                const CellSimilarity::Similarity                              cell_similarity) const
{
  // convert data object to internal
  // data for this class. fails with
  // an exception if that is not
  // possible
  Assert (dynamic_cast<const InternalData *> (&fedata) != 0, ExcInternalError());
  const InternalData &fe_data = static_cast<const InternalData &> (fedata);

  const UpdateFlags flags(fe_data.update_each);

  if (flags & update_values)
    for (unsigned int k=0; k<this->dofs_per_cell; ++k)
      for (unsigned int i=0; i<quadrature.size(); ++i)
        output_data.shape_values(k,i) = fe_data.shape_values[k][i];

  if (flags & update_gradients && cell_similarity != CellSimilarity::translation)
    for (unsigned int k=0; k<this->dofs_per_cell; ++k)
      mapping.transform (fe_data.shape_gradients[k],
                         mapping_covariant,
                         mapping_internal,
                         output_data.shape_gradients[k]);

  if (flags & update_hessians && cell_similarity != CellSimilarity::translation)
    {
      for (unsigned int k=0; k<this->dofs_per_cell; ++k)
        mapping.transform (fe_data.shape_hessians[k],
                           mapping_covariant_gradient,
                           mapping_internal,
                           output_data.shape_hessians[k]);

      for (unsigned int k=0; k<this->dofs_per_cell; ++k)
        for (unsigned int i=0; i<quadrature.size(); ++i)
          for (unsigned int j=0; j<spacedim; ++j)
            output_data.shape_hessians[k][i] -=
              mapping_data.jacobian_pushed_forward_grads[i][j]
              * output_data.shape_gradients[k][i][j];
    }

  if (flags & update_3rd_derivatives && cell_similarity != CellSimilarity::translation)
    {
      for (unsigned int k=0; k<this->dofs_per_cell; ++k)
        mapping.transform (fe_data.shape_3rd_derivatives[k],
                           mapping_covariant_hessian,
                           mapping_internal,
                           output_data.shape_3rd_derivatives[k]);

      for (unsigned int k=0; k<this->dofs_per_cell; ++k)
        correct_third_derivatives(output_data, mapping_data, quadrature.size(), k);
    }
}



template <class POLY, int dim, int spacedim>
void
FE_Poly<POLY,dim,spacedim>::
fill_fe_face_values (const Mapping<dim,spacedim>                                  &mapping,
                     const typename Triangulation<dim,spacedim>::cell_iterator    &cell,
                     const unsigned int                                            face,
                     const Quadrature<dim-1>                                      &quadrature,
                     const typename Mapping<dim,spacedim>::InternalDataBase       &mapping_internal,
                     const typename FiniteElement<dim,spacedim>::InternalDataBase &fedata,
                     const internal::FEValues::MappingRelatedData<dim,spacedim>   &mapping_data,
                     internal::FEValues::FiniteElementRelatedData<dim,spacedim>   &output_data) const
{
  // convert data object to internal
  // data for this class. fails with
  // an exception if that is not
  // possible
  Assert (dynamic_cast<const InternalData *> (&fedata) != 0, ExcInternalError());
  const InternalData &fe_data = static_cast<const InternalData &> (fedata);

  // offset determines which data set
  // to take (all data sets for all
  // faces are stored contiguously)

  const typename QProjector<dim>::DataSetDescriptor offset
    = QProjector<dim>::DataSetDescriptor::face (face,
                                                cell->face_orientation(face),
                                                cell->face_flip(face),
                                                cell->face_rotation(face),
                                                quadrature.size());

  const UpdateFlags flags(fe_data.update_each);

  if (flags & update_values)
    for (unsigned int k=0; k<this->dofs_per_cell; ++k)
      for (unsigned int i=0; i<quadrature.size(); ++i)
        output_data.shape_values(k,i) = fe_data.shape_values[k][i+offset];

  if (flags & update_gradients)
    for (unsigned int k=0; k<this->dofs_per_cell; ++k)
      mapping.transform (make_slice(fe_data.shape_gradients[k], offset, quadrature.size()),
                         mapping_covariant,
                         mapping_internal,
                         output_data.shape_gradients[k]);

  if (flags & update_hessians)
    {
      for (unsigned int k=0; k<this->dofs_per_cell; ++k)
        mapping.transform (make_slice(fe_data.shape_hessians[k],
                                      offset,
                                      quadrature.size()),
                           mapping_covariant_gradient,
                           mapping_internal,
                           output_data.shape_hessians[k]);

      for (unsigned int k=0; k<this->dofs_per_cell; ++k)
        for (unsigned int i=0; i<quadrature.size(); ++i)
          for (unsigned int j=0; j<spacedim; ++j)
            output_data.shape_hessians[k][i] -=
              mapping_data.jacobian_pushed_forward_grads[i][j]
              * output_data.shape_gradients[k][i][j];
    }

  if (flags & update_3rd_derivatives)
    {
      for (unsigned int k=0; k<this->dofs_per_cell; ++k)
        mapping.transform (make_slice(fe_data.shape_3rd_derivatives[k],
                                      offset,
                                      quadrature.size()),
                           mapping_covariant_hessian,
                           mapping_internal,
                           output_data.shape_3rd_derivatives[k]);

      for (unsigned int k=0; k<this->dofs_per_cell; ++k)
        correct_third_derivatives(output_data, mapping_data, quadrature.size(), k);
    }
}

template <class POLY, int dim, int spacedim>
void
FE_Poly<POLY,dim,spacedim>::
fill_fe_subface_values (const Mapping<dim,spacedim>                                  &mapping,
                        const typename Triangulation<dim,spacedim>::cell_iterator    &cell,
                        const unsigned int                                            face,
                        const unsigned int                                            subface,
                        const Quadrature<dim-1>                                      &quadrature,
                        const typename Mapping<dim,spacedim>::InternalDataBase       &mapping_internal,
                        const typename FiniteElement<dim,spacedim>::InternalDataBase &fedata,
                        const internal::FEValues::MappingRelatedData<dim,spacedim>   &mapping_data,
                        internal::FEValues::FiniteElementRelatedData<dim,spacedim>   &output_data) const
{
  // convert data object to internal
  // data for this class. fails with
  // an exception if that is not
  // possible
  Assert (dynamic_cast<const InternalData *> (&fedata) != 0, ExcInternalError());
  const InternalData &fe_data = static_cast<const InternalData &> (fedata);

  // offset determines which data set
  // to take (all data sets for all
  // sub-faces are stored contiguously)

  const typename QProjector<dim>::DataSetDescriptor offset
    = QProjector<dim>::DataSetDescriptor::subface (face, subface,
                                                   cell->face_orientation(face),
                                                   cell->face_flip(face),
                                                   cell->face_rotation(face),
                                                   quadrature.size(),
                                                   cell->subface_case(face));

  const UpdateFlags flags(fe_data.update_each);

  if (flags & update_values)
    for (unsigned int k=0; k<this->dofs_per_cell; ++k)
      for (unsigned int i=0; i<quadrature.size(); ++i)
        output_data.shape_values(k,i) = fe_data.shape_values[k][i+offset];

  if (flags & update_gradients)
    for (unsigned int k=0; k<this->dofs_per_cell; ++k)
      mapping.transform (make_slice(fe_data.shape_gradients[k], offset, quadrature.size()),
                         mapping_covariant,
                         mapping_internal,
                         output_data.shape_gradients[k]);

  if (flags & update_hessians)
    {
      for (unsigned int k=0; k<this->dofs_per_cell; ++k)
        mapping.transform (make_slice(fe_data.shape_hessians[k],
                                      offset,
                                      quadrature.size()),
                           mapping_covariant_gradient,
                           mapping_internal,
                           output_data.shape_hessians[k]);

      for (unsigned int k=0; k<this->dofs_per_cell; ++k)
        for (unsigned int i=0; i<quadrature.size(); ++i)
          for (unsigned int j=0; j<spacedim; ++j)
            output_data.shape_hessians[k][i] -=
              mapping_data.jacobian_pushed_forward_grads[i][j]
              * output_data.shape_gradients[k][i][j];
    }

  if (flags & update_3rd_derivatives)
    {
      for (unsigned int k=0; k<this->dofs_per_cell; ++k)
        mapping.transform (make_slice(fe_data.shape_3rd_derivatives[k],
                                      offset,
                                      quadrature.size()),
                           mapping_covariant_hessian,
                           mapping_internal,
                           output_data.shape_3rd_derivatives[k]);

      for (unsigned int k=0; k<this->dofs_per_cell; ++k)
        correct_third_derivatives(output_data, mapping_data, quadrature.size(), k);
    }
}



template <class POLY, int dim, int spacedim>
inline void
FE_Poly<POLY,dim,spacedim>::
correct_third_derivatives (internal::FEValues::FiniteElementRelatedData<dim,spacedim>       &output_data,
                           const internal::FEValues::MappingRelatedData<dim,spacedim>       &mapping_data,
                           const unsigned int                                                n_q_points,
                           const unsigned int                                                dof) const
{
  for (unsigned int i=0; i<n_q_points; ++i)
    for (unsigned int j=0; j<spacedim; ++j)
      for (unsigned int k=0; k<spacedim; ++k)
        for (unsigned int l=0; l<spacedim; ++l)
          for (unsigned int m=0; m<spacedim; ++m)
            {
              output_data.shape_3rd_derivatives[dof][i][j][k][l] -=
                (mapping_data.jacobian_pushed_forward_grads[i][m][j][l] *
                 output_data.shape_hessians[dof][i][k][m])
                + (mapping_data.jacobian_pushed_forward_grads[i][m][k][l] *
                   output_data.shape_hessians[dof][i][j][m])
                + (mapping_data.jacobian_pushed_forward_grads[i][m][j][k] *
                   output_data.shape_hessians[dof][i][l][m])
                + (mapping_data.jacobian_pushed_forward_2nd_derivatives[i][m][j][k][l] *
                   output_data.shape_gradients[dof][i][m]);
            }

}

namespace internal
{
  template <class POLY>
  inline
  std::vector<unsigned int>
  get_poly_space_numbering (const POLY &)
  {
    Assert (false, ExcNotImplemented());
    return std::vector<unsigned int>();
  }

  template <class POLY>
  inline
  std::vector<unsigned int>
  get_poly_space_numbering_inverse (const POLY &)
  {
    Assert (false, ExcNotImplemented());
    return std::vector<unsigned int>();
  }

  template <int dim, typename POLY>
  inline
  std::vector<unsigned int>
  get_poly_space_numbering (const TensorProductPolynomials<dim,POLY> &poly)
  {
    return poly.get_numbering();
  }

  template <int dim, typename POLY>
  inline
  std::vector<unsigned int>
  get_poly_space_numbering_inverse (const TensorProductPolynomials<dim,POLY> &poly)
  {
    return poly.get_numbering_inverse();
  }

  template <int dim>
  inline
  std::vector<unsigned int>
  get_poly_space_numbering (const TensorProductPolynomialsConst<dim> &poly)
  {
    return poly.get_numbering();
  }

  template <int dim>
  inline
  std::vector<unsigned int>
  get_poly_space_numbering_inverse (const TensorProductPolynomialsConst<dim> &poly)
  {
    return poly.get_numbering_inverse();
  }
}



template <class POLY, int dim, int spacedim>
std::vector<unsigned int>
FE_Poly<POLY,dim,spacedim>::get_poly_space_numbering () const
{
  return internal::get_poly_space_numbering (poly_space);
}




template <class POLY, int dim, int spacedim>
std::vector<unsigned int>
FE_Poly<POLY,dim,spacedim>::get_poly_space_numbering_inverse () const
{
  return internal::get_poly_space_numbering_inverse (poly_space);
}



DEAL_II_NAMESPACE_CLOSE
