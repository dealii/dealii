// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/config.h>

#include <deal.II/base/polynomial_space.h>
#include <deal.II/base/polynomials_p.h>
#include <deal.II/base/polynomials_piecewise.h>
#include <deal.II/base/polynomials_rannacher_turek.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/tensor_product_polynomials_bubbles.h>
#include <deal.II/base/tensor_product_polynomials_const.h>

#include <deal.II/fe/fe_poly.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_cartesian.h>

DEAL_II_NAMESPACE_OPEN

#ifndef DOXYGEN

template <int dim, int spacedim>
FE_Poly<dim, spacedim>::FE_Poly(const FE_Poly &fe)
  : FiniteElement<dim, spacedim>(fe)
  , poly_space(fe.poly_space->clone())
{}

template <int dim, int spacedim>
FE_Poly<dim, spacedim>::FE_Poly(
  const ScalarPolynomialsBase<dim> &poly_space,
  const FiniteElementData<dim>     &fe_data,
  const std::vector<bool>          &restriction_is_additive_flags,
  const std::vector<ComponentMask> &nonzero_components)
  : FiniteElement<dim, spacedim>(fe_data,
                                 restriction_is_additive_flags,
                                 nonzero_components)
  , poly_space(poly_space.clone())
{}


template <int dim, int spacedim>
unsigned int
FE_Poly<dim, spacedim>::get_degree() const
{
  return this->degree;
}


template <int dim, int spacedim>
double
FE_Poly<dim, spacedim>::shape_value(const unsigned int i,
                                    const Point<dim>  &p) const
{
  AssertIndexRange(i, this->n_dofs_per_cell());
  return poly_space->compute_value(i, p);
}


template <int dim, int spacedim>
double
FE_Poly<dim, spacedim>::shape_value_component(
  const unsigned int i,
  const Point<dim>  &p,
  const unsigned int component) const
{
  (void)component;
  AssertIndexRange(i, this->n_dofs_per_cell());
  AssertIndexRange(component, 1);
  return poly_space->compute_value(i, p);
}



template <int dim, int spacedim>
Tensor<1, dim>
FE_Poly<dim, spacedim>::shape_grad(const unsigned int i,
                                   const Point<dim>  &p) const
{
  AssertIndexRange(i, this->n_dofs_per_cell());
  return poly_space->template compute_derivative<1>(i, p);
}



template <int dim, int spacedim>
Tensor<1, dim>
FE_Poly<dim, spacedim>::shape_grad_component(const unsigned int i,
                                             const Point<dim>  &p,
                                             const unsigned int component) const
{
  (void)component;
  AssertIndexRange(i, this->n_dofs_per_cell());
  AssertIndexRange(component, 1);
  return poly_space->template compute_derivative<1>(i, p);
}



template <int dim, int spacedim>
Tensor<2, dim>
FE_Poly<dim, spacedim>::shape_grad_grad(const unsigned int i,
                                        const Point<dim>  &p) const
{
  AssertIndexRange(i, this->n_dofs_per_cell());
  return poly_space->template compute_derivative<2>(i, p);
}



template <int dim, int spacedim>
Tensor<2, dim>
FE_Poly<dim, spacedim>::shape_grad_grad_component(
  const unsigned int i,
  const Point<dim>  &p,
  const unsigned int component) const
{
  (void)component;
  AssertIndexRange(i, this->n_dofs_per_cell());
  AssertIndexRange(component, 1);
  return poly_space->template compute_derivative<2>(i, p);
}



template <int dim, int spacedim>
Tensor<3, dim>
FE_Poly<dim, spacedim>::shape_3rd_derivative(const unsigned int i,
                                             const Point<dim>  &p) const
{
  AssertIndexRange(i, this->n_dofs_per_cell());
  return poly_space->template compute_derivative<3>(i, p);
}



template <int dim, int spacedim>
Tensor<3, dim>
FE_Poly<dim, spacedim>::shape_3rd_derivative_component(
  const unsigned int i,
  const Point<dim>  &p,
  const unsigned int component) const
{
  (void)component;
  AssertIndexRange(i, this->n_dofs_per_cell());
  AssertIndexRange(component, 1);
  return poly_space->template compute_derivative<3>(i, p);
}



template <int dim, int spacedim>
Tensor<4, dim>
FE_Poly<dim, spacedim>::shape_4th_derivative(const unsigned int i,
                                             const Point<dim>  &p) const
{
  AssertIndexRange(i, this->n_dofs_per_cell());
  return poly_space->template compute_derivative<4>(i, p);
}



template <int dim, int spacedim>
Tensor<4, dim>
FE_Poly<dim, spacedim>::shape_4th_derivative_component(
  const unsigned int i,
  const Point<dim>  &p,
  const unsigned int component) const
{
  (void)component;
  AssertIndexRange(i, this->n_dofs_per_cell());
  AssertIndexRange(component, 1);
  return poly_space->template compute_derivative<4>(i, p);
}



//---------------------------------------------------------------------------
// Auxiliary functions
//---------------------------------------------------------------------------


template <int dim, int spacedim>
UpdateFlags
FE_Poly<dim, spacedim>::requires_update_flags(const UpdateFlags flags) const
{
  UpdateFlags out = update_default;

  if (flags & update_values)
    out |= update_values;
  if (flags & update_gradients)
    out |= update_gradients | update_covariant_transformation;
  if (flags & update_hessians)
    out |= update_hessians | update_covariant_transformation |
           update_gradients | update_jacobian_pushed_forward_grads;
  if (flags & update_3rd_derivatives)
    out |= update_3rd_derivatives | update_covariant_transformation |
           update_hessians | update_gradients |
           update_jacobian_pushed_forward_grads |
           update_jacobian_pushed_forward_2nd_derivatives;
  if (flags & update_normal_vectors)
    out |= update_normal_vectors | update_JxW_values;

  return out;
}



//---------------------------------------------------------------------------
// Fill data of FEValues
//---------------------------------------------------------------------------



/**
 * Returns whether we need to correct the Hessians and third derivatives with
 * the derivatives of the Jacobian. This is determined by checking if
 * the jacobian_pushed_forward are zero.
 *
 * Especially for the third derivatives, the correction term is very expensive,
 * which is why we check if the derivatives are zero before computing the
 * correction.
 */
template <int dim, int spacedim>
bool
higher_derivatives_need_correcting(
  const Mapping<dim, spacedim> &mapping,
  const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
                    &mapping_data,
  const unsigned int n_q_points,
  const UpdateFlags  update_flags)
{
  // If higher derivatives weren't requested we don't need to correct them.
  const bool update_higher_derivatives =
    (update_flags & update_hessians) || (update_flags & update_3rd_derivatives);
  if (!update_higher_derivatives)
    return false;

  // If we have a Cartesian mapping, we know that jacoban_pushed_forward_grads
  // are identically zero.
  if (dynamic_cast<const MappingCartesian<dim> *>(&mapping))
    return false;

  // Here, we should check if jacobian_pushed_forward_grads are zero at the
  // quadrature points. This is yet to be implemented.
  (void)mapping_data;
  (void)n_q_points;

  return true;
}



template <int dim, int spacedim>
void
FE_Poly<dim, spacedim>::fill_fe_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &,
  const CellSimilarity::Similarity                         cell_similarity,
  const Quadrature<dim>                                   &quadrature,
  const Mapping<dim, spacedim>                            &mapping,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
  const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
                                                                &mapping_data,
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

  const UpdateFlags flags(fe_data.update_each);

  const bool need_to_correct_higher_derivatives =
    higher_derivatives_need_correcting(mapping,
                                       mapping_data,
                                       quadrature.size(),
                                       flags);

  // transform gradients and higher derivatives. there is nothing to do
  // for values since we already emplaced them into output_data when
  // we were in get_data()
  if ((flags & update_gradients) &&
      (cell_similarity != CellSimilarity::translation))
    for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
      mapping.transform(make_array_view(fe_data.shape_gradients, k),
                        mapping_covariant,
                        mapping_internal,
                        make_array_view(output_data.shape_gradients, k));

  if ((flags & update_hessians) &&
      (cell_similarity != CellSimilarity::translation))
    {
      for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
        mapping.transform(make_array_view(fe_data.shape_hessians, k),
                          mapping_covariant_gradient,
                          mapping_internal,
                          make_array_view(output_data.shape_hessians, k));

      if (need_to_correct_higher_derivatives)
        correct_hessians(output_data, mapping_data, quadrature.size());
    }

  if ((flags & update_3rd_derivatives) &&
      (cell_similarity != CellSimilarity::translation))
    {
      for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
        mapping.transform(make_array_view(fe_data.shape_3rd_derivatives, k),
                          mapping_covariant_hessian,
                          mapping_internal,
                          make_array_view(output_data.shape_3rd_derivatives,
                                          k));

      if (need_to_correct_higher_derivatives)
        correct_third_derivatives(output_data, mapping_data, quadrature.size());
    }
}



template <int dim, int spacedim>
void
FE_Poly<dim, spacedim>::fill_fe_face_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const hp::QCollection<dim - 1>                             &quadrature,
  const Mapping<dim, spacedim>                               &mapping,
  const typename Mapping<dim, spacedim>::InternalDataBase    &mapping_internal,
  const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
                                                                &mapping_data,
  const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    &output_data) const
{
  const unsigned int n_q_points =
    quadrature[quadrature.size() == 1 ? 0 : face_no].size();

  // convert data object to internal
  // data for this class. fails with
  // an exception if that is not
  // possible
  Assert(dynamic_cast<const InternalData *>(&fe_internal) != nullptr,
         ExcInternalError());
  const InternalData &fe_data = static_cast<const InternalData &>(fe_internal);

  // offset determines which data set
  // to take (all data sets for all
  // faces are stored contiguously)

  const auto offset =
    QProjector<dim>::DataSetDescriptor::face(this->reference_cell(),
                                             face_no,
                                             cell->combined_face_orientation(
                                               face_no),
                                             quadrature);

  const UpdateFlags flags(fe_data.update_each);

  const bool need_to_correct_higher_derivatives =
    higher_derivatives_need_correcting(mapping,
                                       mapping_data,
                                       n_q_points,
                                       flags);

  // transform gradients and higher derivatives. we also have to copy
  // the values (unlike in the case of fill_fe_values()) since
  // we need to take into account the offsets
  if (flags & update_values)
    for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
      for (unsigned int i = 0; i < n_q_points; ++i)
        output_data.shape_values(k, i) = fe_data.shape_values[k][i + offset];

  if (flags & update_gradients)
    for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
      mapping.transform(
        make_array_view(fe_data.shape_gradients, k, offset, n_q_points),
        mapping_covariant,
        mapping_internal,
        make_array_view(output_data.shape_gradients, k));

  if (flags & update_hessians)
    {
      for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
        mapping.transform(
          make_array_view(fe_data.shape_hessians, k, offset, n_q_points),
          mapping_covariant_gradient,
          mapping_internal,
          make_array_view(output_data.shape_hessians, k));

      if (need_to_correct_higher_derivatives)
        correct_hessians(output_data, mapping_data, n_q_points);
    }

  if (flags & update_3rd_derivatives)
    {
      for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
        mapping.transform(
          make_array_view(fe_data.shape_3rd_derivatives, k, offset, n_q_points),
          mapping_covariant_hessian,
          mapping_internal,
          make_array_view(output_data.shape_3rd_derivatives, k));

      if (need_to_correct_higher_derivatives)
        correct_third_derivatives(output_data, mapping_data, n_q_points);
    }
}



template <int dim, int spacedim>
void
FE_Poly<dim, spacedim>::fill_fe_subface_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const unsigned int                                          sub_no,
  const Quadrature<dim - 1>                                  &quadrature,
  const Mapping<dim, spacedim>                               &mapping,
  const typename Mapping<dim, spacedim>::InternalDataBase    &mapping_internal,
  const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
                                                                &mapping_data,
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

  // offset determines which data set
  // to take (all data sets for all
  // sub-faces are stored contiguously)

  const auto offset =
    QProjector<dim>::DataSetDescriptor::subface(this->reference_cell(),
                                                face_no,
                                                sub_no,
                                                cell->combined_face_orientation(
                                                  face_no),
                                                quadrature.size(),
                                                cell->subface_case(face_no));

  const UpdateFlags flags(fe_data.update_each);

  const bool need_to_correct_higher_derivatives =
    higher_derivatives_need_correcting(mapping,
                                       mapping_data,
                                       quadrature.size(),
                                       flags);

  // transform gradients and higher derivatives. we also have to copy
  // the values (unlike in the case of fill_fe_values()) since
  // we need to take into account the offsets
  if (flags & update_values)
    for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
      for (unsigned int i = 0; i < quadrature.size(); ++i)
        output_data.shape_values(k, i) = fe_data.shape_values[k][i + offset];

  if (flags & update_gradients)
    for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
      mapping.transform(
        make_array_view(fe_data.shape_gradients, k, offset, quadrature.size()),
        mapping_covariant,
        mapping_internal,
        make_array_view(output_data.shape_gradients, k));

  if (flags & update_hessians)
    {
      for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
        mapping.transform(
          make_array_view(fe_data.shape_hessians, k, offset, quadrature.size()),
          mapping_covariant_gradient,
          mapping_internal,
          make_array_view(output_data.shape_hessians, k));

      if (need_to_correct_higher_derivatives)
        correct_hessians(output_data, mapping_data, quadrature.size());
    }

  if (flags & update_3rd_derivatives)
    {
      for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
        mapping.transform(make_array_view(fe_data.shape_3rd_derivatives,
                                          k,
                                          offset,
                                          quadrature.size()),
                          mapping_covariant_hessian,
                          mapping_internal,
                          make_array_view(output_data.shape_3rd_derivatives,
                                          k));

      if (need_to_correct_higher_derivatives)
        correct_third_derivatives(output_data, mapping_data, quadrature.size());
    }
}



template <int dim, int spacedim>
inline void
FE_Poly<dim, spacedim>::correct_hessians(
  internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim>
    &output_data,
  const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
                    &mapping_data,
  const unsigned int n_q_points) const
{
  for (unsigned int dof = 0; dof < this->n_dofs_per_cell(); ++dof)
    for (unsigned int i = 0; i < n_q_points; ++i)
      for (unsigned int j = 0; j < spacedim; ++j)
        output_data.shape_hessians[dof][i] -=
          mapping_data.jacobian_pushed_forward_grads[i][j] *
          output_data.shape_gradients[dof][i][j];
}



template <int dim, int spacedim>
inline void
FE_Poly<dim, spacedim>::correct_third_derivatives(
  internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim>
    &output_data,
  const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
                    &mapping_data,
  const unsigned int n_q_points) const
{
  for (unsigned int dof = 0; dof < this->n_dofs_per_cell(); ++dof)
    for (unsigned int i = 0; i < n_q_points; ++i)
      for (unsigned int j = 0; j < spacedim; ++j)
        for (unsigned int k = 0; k < spacedim; ++k)
          for (unsigned int l = 0; l < spacedim; ++l)
            for (unsigned int m = 0; m < spacedim; ++m)
              {
                output_data.shape_3rd_derivatives[dof][i][j][k][l] -=
                  (mapping_data.jacobian_pushed_forward_grads[i][m][j][l] *
                   output_data.shape_hessians[dof][i][k][m]) +
                  (mapping_data.jacobian_pushed_forward_grads[i][m][k][l] *
                   output_data.shape_hessians[dof][i][j][m]) +
                  (mapping_data.jacobian_pushed_forward_grads[i][m][j][k] *
                   output_data.shape_hessians[dof][i][l][m]) +
                  (mapping_data
                     .jacobian_pushed_forward_2nd_derivatives[i][m][j][k][l] *
                   output_data.shape_gradients[dof][i][m]);
              }
}



template <int dim, int spacedim>
inline const ScalarPolynomialsBase<dim> &
FE_Poly<dim, spacedim>::get_poly_space() const
{
  return *poly_space;
}



template <int dim, int spacedim>
std::vector<unsigned int>
FE_Poly<dim, spacedim>::get_poly_space_numbering() const
{
  auto *const space_tensor_prod =
    dynamic_cast<TensorProductPolynomials<dim> *>(this->poly_space.get());
  if (space_tensor_prod != nullptr)
    return space_tensor_prod->get_numbering();

  auto *const space_tensor_prod_aniso =
    dynamic_cast<AnisotropicPolynomials<dim> *>(this->poly_space.get());
  if (space_tensor_prod_aniso != nullptr)
    return space_tensor_prod_aniso->get_numbering();

  auto *const space_tensor_prod_piecewise = dynamic_cast<
    TensorProductPolynomials<dim, Polynomials::PiecewisePolynomial<double>> *>(
    this->poly_space.get());
  if (space_tensor_prod_piecewise != nullptr)
    return space_tensor_prod_piecewise->get_numbering();

  auto *const space_tensor_prod_bubbles =
    dynamic_cast<TensorProductPolynomialsBubbles<dim> *>(
      this->poly_space.get());
  if (space_tensor_prod_bubbles != nullptr)
    return space_tensor_prod_bubbles->get_numbering();

  auto *const space_tensor_prod_const =
    dynamic_cast<TensorProductPolynomialsConst<dim> *>(this->poly_space.get());
  if (space_tensor_prod_const != nullptr)
    return space_tensor_prod_const->get_numbering();

  DEAL_II_NOT_IMPLEMENTED();
  return std::vector<unsigned int>();
}



template <int dim, int spacedim>
std::vector<unsigned int>
FE_Poly<dim, spacedim>::get_poly_space_numbering_inverse() const
{
  return Utilities::invert_permutation(get_poly_space_numbering());
}



template <int dim, int spacedim>
std::size_t
FE_Poly<dim, spacedim>::memory_consumption() const
{
  return FiniteElement<dim, spacedim>::memory_consumption() +
         poly_space->memory_consumption();
}

#endif

#include "fe/fe_poly.inst"

DEAL_II_NAMESPACE_CLOSE
