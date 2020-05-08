// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#include <deal.II/base/array_view.h>
#include <deal.II/base/derivative_form.h>
#include <deal.II/base/polynomials_abf.h>
#include <deal.II/base/polynomials_bdm.h>
#include <deal.II/base/polynomials_bernardi_raugel.h>
#include <deal.II/base/polynomials_nedelec.h>
#include <deal.II/base/polynomials_raviart_thomas.h>
#include <deal.II/base/polynomials_rt_bubbles.h>
#include <deal.II/base/qprojector.h>

#include <deal.II/fe/fe_poly_tensor.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_cartesian.h>

#include <deal.II/grid/tria.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace FE_PolyTensor
  {
    namespace
    {
      //---------------------------------------------------------------------------
      // Utility method, which is used to determine the change of sign for
      // the DoFs on the faces of the given cell.
      //---------------------------------------------------------------------------

      /**
       * On non-Cartesian grids, the sign of the DoFs associated with the faces
       * of the elements has to be changed in some cases.  This procedure
       * implements an algorithm that determines those DoFs that need this
       * sign change for a given cell.
       */
      template <int spacedim>
      void
      get_face_sign_change_rt(const dealii::Triangulation<1>::cell_iterator &,
                              const FiniteElement<1, spacedim> &,
                              const std::vector<MappingKind> &,
                              std::vector<double> &)
      {
        // nothing to do in 1d
      }



      //      template<int spacedim>
      void
      get_face_sign_change_rt(
        const dealii::Triangulation<2>::cell_iterator &cell,
        const FiniteElement<2, 2> &                    fe,
        const std::vector<MappingKind> &               mapping_kind,
        std::vector<double> &                          face_sign)
      {
        const unsigned int dim      = 2;
        const unsigned int spacedim = 2;

        for (unsigned int f = GeometryInfo<dim>::faces_per_cell / 2;
             f < GeometryInfo<dim>::faces_per_cell;
             ++f)
          {
            dealii::Triangulation<dim, spacedim>::face_iterator face =
              cell->face(f);
            if (!face->at_boundary())
              {
                const unsigned int nn = cell->neighbor_face_no(f);

                if (nn < GeometryInfo<dim>::faces_per_cell / 2)
                  for (unsigned int j = 0; j < fe.dofs_per_face; ++j)
                    {
                      const unsigned int cell_j = fe.face_to_cell_index(j, f);

                      Assert(f * fe.dofs_per_face + j < face_sign.size(),
                             ExcInternalError());
                      Assert(mapping_kind.size() == 1 ||
                               cell_j < mapping_kind.size(),
                             ExcInternalError());

                      // TODO: This is probably only going to work for those
                      // elements for which all dofs are face dofs
                      if ((mapping_kind.size() > 1 ?
                             mapping_kind[cell_j] :
                             mapping_kind[0]) == mapping_raviart_thomas)
                        face_sign[f * fe.dofs_per_face + j] = -1.0;
                    }
              }
          }
      }



      template <int spacedim>
      void
      get_face_sign_change_rt(
        const dealii::Triangulation<3>::cell_iterator & /*cell*/,
        const FiniteElement<3, spacedim> & /*fe*/,
        const std::vector<MappingKind> & /*mapping_kind*/,
        std::vector<double> & /*face_sign*/)
      {
        // TODO: think about what it would take here
      }



      template <int spacedim>
      void
      get_face_sign_change_nedelec(
        const dealii::Triangulation<1>::cell_iterator & /*cell*/,
        const FiniteElement<1, spacedim> & /*fe*/,
        const std::vector<MappingKind> & /*mapping_kind*/,
        std::vector<double> & /*face_sign*/)
      {
        // nothing to do in 1d
      }



      template <int spacedim>
      void
      get_face_sign_change_nedelec(
        const dealii::Triangulation<2>::cell_iterator & /*cell*/,
        const FiniteElement<2, spacedim> & /*fe*/,
        const std::vector<MappingKind> & /*mapping_kind*/,
        std::vector<double> & /*face_sign*/)
      {
        // TODO: think about what it would take here
      }


      template <int spacedim>
      void
      get_face_sign_change_nedelec(
        const dealii::Triangulation<3>::cell_iterator &cell,
        const FiniteElement<3, spacedim> & /*fe*/,
        const std::vector<MappingKind> &mapping_kind,
        std::vector<double> &           face_sign)
      {
        const unsigned int dim = 3;
        // TODO: This is probably only going to work for those elements for
        // which all dofs are face dofs
        for (unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell; ++l)
          if (!(cell->line_orientation(l)) &&
              mapping_kind[0] == mapping_nedelec)
            face_sign[l] = -1.0;
      }
    } // namespace
  }   // namespace FE_PolyTensor
} // namespace internal



template <int dim, int spacedim>
FE_PolyTensor<dim, spacedim>::FE_PolyTensor(
  const TensorPolynomialsBase<dim> &polynomials,
  const FiniteElementData<dim> &    fe_data,
  const std::vector<bool> &         restriction_is_additive_flags,
  const std::vector<ComponentMask> &nonzero_components)
  : FiniteElement<dim, spacedim>(fe_data,
                                 restriction_is_additive_flags,
                                 nonzero_components)
  , mapping_kind({MappingKind::mapping_none})
  , poly_space(polynomials.clone())
{
  cached_point(0) = -1;
  // Set up the table converting
  // components to base
  // components. Since we have only
  // one base element, everything
  // remains zero except the
  // component in the base, which is
  // the component itself
  for (unsigned int comp = 0; comp < this->n_components(); ++comp)
    this->component_to_base_table[comp].first.second = comp;
}



template <int dim, int spacedim>
FE_PolyTensor<dim, spacedim>::FE_PolyTensor(const FE_PolyTensor &fe)
  : FiniteElement<dim, spacedim>(fe)
  , mapping_kind(fe.mapping_kind)
  , poly_space(fe.poly_space->clone())
  , inverse_node_matrix(fe.inverse_node_matrix)
{}



template <int dim, int spacedim>
bool
FE_PolyTensor<dim, spacedim>::single_mapping_kind() const
{
  return mapping_kind.size() == 1;
}



template <int dim, int spacedim>
MappingKind
FE_PolyTensor<dim, spacedim>::get_mapping_kind(const unsigned int i) const
{
  if (single_mapping_kind())
    return mapping_kind[0];

  AssertIndexRange(i, mapping_kind.size());
  return mapping_kind[i];
}



template <int dim, int spacedim>
double
FE_PolyTensor<dim, spacedim>::shape_value(const unsigned int,
                                          const Point<dim> &) const

{
  Assert(false, (typename FiniteElement<dim, spacedim>::ExcFENotPrimitive()));
  return 0.;
}



template <int dim, int spacedim>
double
FE_PolyTensor<dim, spacedim>::shape_value_component(
  const unsigned int i,
  const Point<dim> & p,
  const unsigned int component) const
{
  AssertIndexRange(i, this->dofs_per_cell);
  AssertIndexRange(component, dim);

  std::lock_guard<std::mutex> lock(cache_mutex);

  if (cached_point != p || cached_values.size() == 0)
    {
      cached_point = p;
      cached_values.resize(poly_space->n());

      std::vector<Tensor<4, dim>> dummy1;
      std::vector<Tensor<5, dim>> dummy2;
      poly_space->evaluate(
        p, cached_values, cached_grads, cached_grad_grads, dummy1, dummy2);
    }

  double s = 0;
  if (inverse_node_matrix.n_cols() == 0)
    return cached_values[i][component];
  else
    for (unsigned int j = 0; j < inverse_node_matrix.n_cols(); ++j)
      s += inverse_node_matrix(j, i) * cached_values[j][component];
  return s;
}



template <int dim, int spacedim>
Tensor<1, dim>
FE_PolyTensor<dim, spacedim>::shape_grad(const unsigned int,
                                         const Point<dim> &) const
{
  Assert(false, (typename FiniteElement<dim, spacedim>::ExcFENotPrimitive()));
  return Tensor<1, dim>();
}



template <int dim, int spacedim>
Tensor<1, dim>
FE_PolyTensor<dim, spacedim>::shape_grad_component(
  const unsigned int i,
  const Point<dim> & p,
  const unsigned int component) const
{
  AssertIndexRange(i, this->dofs_per_cell);
  AssertIndexRange(component, dim);

  std::lock_guard<std::mutex> lock(cache_mutex);

  if (cached_point != p || cached_grads.size() == 0)
    {
      cached_point = p;
      cached_grads.resize(poly_space->n());

      std::vector<Tensor<4, dim>> dummy1;
      std::vector<Tensor<5, dim>> dummy2;
      poly_space->evaluate(
        p, cached_values, cached_grads, cached_grad_grads, dummy1, dummy2);
    }

  Tensor<1, dim> s;
  if (inverse_node_matrix.n_cols() == 0)
    return cached_grads[i][component];
  else
    for (unsigned int j = 0; j < inverse_node_matrix.n_cols(); ++j)
      s += inverse_node_matrix(j, i) * cached_grads[j][component];

  return s;
}



template <int dim, int spacedim>
Tensor<2, dim>
FE_PolyTensor<dim, spacedim>::shape_grad_grad(const unsigned int,
                                              const Point<dim> &) const
{
  Assert(false, (typename FiniteElement<dim, spacedim>::ExcFENotPrimitive()));
  return Tensor<2, dim>();
}



template <int dim, int spacedim>
Tensor<2, dim>
FE_PolyTensor<dim, spacedim>::shape_grad_grad_component(
  const unsigned int i,
  const Point<dim> & p,
  const unsigned int component) const
{
  AssertIndexRange(i, this->dofs_per_cell);
  AssertIndexRange(component, dim);

  std::lock_guard<std::mutex> lock(cache_mutex);

  if (cached_point != p || cached_grad_grads.size() == 0)
    {
      cached_point = p;
      cached_grad_grads.resize(poly_space->n());

      std::vector<Tensor<4, dim>> dummy1;
      std::vector<Tensor<5, dim>> dummy2;
      poly_space->evaluate(
        p, cached_values, cached_grads, cached_grad_grads, dummy1, dummy2);
    }

  Tensor<2, dim> s;
  if (inverse_node_matrix.n_cols() == 0)
    return cached_grad_grads[i][component];
  else
    for (unsigned int j = 0; j < inverse_node_matrix.n_cols(); ++j)
      s += inverse_node_matrix(i, j) * cached_grad_grads[j][component];

  return s;
}


//---------------------------------------------------------------------------
// Fill data of FEValues
//---------------------------------------------------------------------------

template <int dim, int spacedim>
void
FE_PolyTensor<dim, spacedim>::fill_fe_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const CellSimilarity::Similarity                            cell_similarity,
  const Quadrature<dim> &                                     quadrature,
  const Mapping<dim, spacedim> &                              mapping,
  const typename Mapping<dim, spacedim>::InternalDataBase &   mapping_internal,
  const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                     spacedim>
    &                                                            mapping_data,
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

  const unsigned int n_q_points = quadrature.size();

  Assert(!(fe_data.update_each & update_values) ||
           fe_data.shape_values.size()[0] == this->dofs_per_cell,
         ExcDimensionMismatch(fe_data.shape_values.size()[0],
                              this->dofs_per_cell));
  Assert(!(fe_data.update_each & update_values) ||
           fe_data.shape_values.size()[1] == n_q_points,
         ExcDimensionMismatch(fe_data.shape_values.size()[1], n_q_points));

  // Create table with sign changes, due to the special structure of the RT
  // elements.
  // TODO: Preliminary hack to demonstrate the overall principle!

  // Compute eventual sign changes depending on the neighborhood
  // between two faces.
  std::fill(fe_data.sign_change.begin(), fe_data.sign_change.end(), 1.0);

  internal::FE_PolyTensor::get_face_sign_change_rt(cell,
                                                   *this,
                                                   this->mapping_kind,
                                                   fe_data.sign_change);

  internal::FE_PolyTensor::get_face_sign_change_nedelec(cell,
                                                        *this,
                                                        this->mapping_kind,
                                                        fe_data.sign_change);


  for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
    {
      const MappingKind mapping_kind = get_mapping_kind(i);

      const unsigned int first =
        output_data.shape_function_to_row_table[i * this->n_components() +
                                                this->get_nonzero_components(i)
                                                  .first_selected_component()];

      // update the shape function values as necessary
      //
      // we only need to do this if the current cell is not a translation of
      // the previous one; or, even if it is a translation, if we use mappings
      // other than the standard mappings that require us to recompute values
      // and derivatives because of possible sign changes
      if (fe_data.update_each & update_values &&
          ((cell_similarity != CellSimilarity::translation) ||
           ((mapping_kind == mapping_piola) ||
            (mapping_kind == mapping_raviart_thomas) ||
            (mapping_kind == mapping_nedelec))))
        {
          switch (mapping_kind)
            {
              case mapping_none:
                {
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_values(first + d, k) =
                        fe_data.shape_values[i][k][d];
                  break;
                }

              case mapping_covariant:
              case mapping_contravariant:
                {
                  mapping.transform(make_array_view(fe_data.shape_values, i),
                                    mapping_kind,
                                    mapping_internal,
                                    make_array_view(
                                      fe_data.transformed_shape_values));

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_values(first + d, k) =
                        fe_data.transformed_shape_values[k][d];

                  break;
                }

              case mapping_raviart_thomas:
              case mapping_piola:
                {
                  mapping.transform(make_array_view(fe_data.shape_values, i),
                                    mapping_piola,
                                    mapping_internal,
                                    make_array_view(
                                      fe_data.transformed_shape_values));
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_values(first + d, k) =
                        fe_data.sign_change[i] *
                        fe_data.transformed_shape_values[k][d];
                  break;
                }

              case mapping_nedelec:
                {
                  mapping.transform(make_array_view(fe_data.shape_values, i),
                                    mapping_covariant,
                                    mapping_internal,
                                    make_array_view(
                                      fe_data.transformed_shape_values));

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_values(first + d, k) =
                        fe_data.sign_change[i] *
                        fe_data.transformed_shape_values[k][d];

                  break;
                }

              default:
                Assert(false, ExcNotImplemented());
            }
        }

      // update gradients. apply the same logic as above
      if (fe_data.update_each & update_gradients &&
          ((cell_similarity != CellSimilarity::translation) ||
           ((mapping_kind == mapping_piola) ||
            (mapping_kind == mapping_raviart_thomas) ||
            (mapping_kind == mapping_nedelec))))

        {
          switch (mapping_kind)
            {
              case mapping_none:
                {
                  mapping.transform(make_array_view(fe_data.shape_grads, i),
                                    mapping_covariant,
                                    mapping_internal,
                                    make_array_view(
                                      fe_data.transformed_shape_grads));
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_gradients[first + d][k] =
                        fe_data.transformed_shape_grads[k][d];
                  break;
                }
              case mapping_covariant:
                {
                  mapping.transform(make_array_view(fe_data.shape_grads, i),
                                    mapping_covariant_gradient,
                                    mapping_internal,
                                    make_array_view(
                                      fe_data.transformed_shape_grads));

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        fe_data.transformed_shape_grads[k][d] -=
                          output_data.shape_values(first + n, k) *
                          mapping_data.jacobian_pushed_forward_grads[k][n][d];

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_gradients[first + d][k] =
                        fe_data.transformed_shape_grads[k][d];

                  break;
                }
              case mapping_contravariant:
                {
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    fe_data.untransformed_shape_grads[k] =
                      fe_data.shape_grads[i][k];
                  mapping.transform(
                    make_array_view(fe_data.untransformed_shape_grads),
                    mapping_contravariant_gradient,
                    mapping_internal,
                    make_array_view(fe_data.transformed_shape_grads));

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        fe_data.transformed_shape_grads[k][d] +=
                          output_data.shape_values(first + n, k) *
                          mapping_data.jacobian_pushed_forward_grads[k][d][n];


                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_gradients[first + d][k] =
                        fe_data.transformed_shape_grads[k][d];

                  break;
                }
              case mapping_raviart_thomas:
              case mapping_piola:
                {
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    fe_data.untransformed_shape_grads[k] =
                      fe_data.shape_grads[i][k];
                  mapping.transform(
                    make_array_view(fe_data.untransformed_shape_grads),
                    mapping_piola_gradient,
                    mapping_internal,
                    make_array_view(fe_data.transformed_shape_grads));

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        fe_data.transformed_shape_grads[k][d] +=
                          (output_data.shape_values(first + n, k) *
                           mapping_data
                             .jacobian_pushed_forward_grads[k][d][n]) -
                          (output_data.shape_values(first + d, k) *
                           mapping_data.jacobian_pushed_forward_grads[k][n][n]);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_gradients[first + d][k] =
                        fe_data.sign_change[i] *
                        fe_data.transformed_shape_grads[k][d];

                  break;
                }

              case mapping_nedelec:
                {
                  // treat the gradients of
                  // this particular shape
                  // function at all
                  // q-points. if Dv is the
                  // gradient of the shape
                  // function on the unit
                  // cell, then
                  // (J^-T)Dv(J^-1) is the
                  // value we want to have on
                  // the real cell.
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    fe_data.untransformed_shape_grads[k] =
                      fe_data.shape_grads[i][k];

                  mapping.transform(
                    make_array_view(fe_data.untransformed_shape_grads),
                    mapping_covariant_gradient,
                    mapping_internal,
                    make_array_view(fe_data.transformed_shape_grads));

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        fe_data.transformed_shape_grads[k][d] -=
                          output_data.shape_values(first + n, k) *
                          mapping_data.jacobian_pushed_forward_grads[k][n][d];

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_gradients[first + d][k] =
                        fe_data.sign_change[i] *
                        fe_data.transformed_shape_grads[k][d];

                  break;
                }

              default:
                Assert(false, ExcNotImplemented());
            }
        }

      // update hessians. apply the same logic as above
      if (fe_data.update_each & update_hessians &&
          ((cell_similarity != CellSimilarity::translation) ||
           ((mapping_kind == mapping_piola) ||
            (mapping_kind == mapping_raviart_thomas) ||
            (mapping_kind == mapping_nedelec))))

        {
          switch (mapping_kind)
            {
              case mapping_none:
                {
                  mapping.transform(
                    make_array_view(fe_data.shape_grad_grads, i),
                    mapping_covariant_gradient,
                    mapping_internal,
                    make_array_view(fe_data.transformed_shape_hessians));

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        fe_data.transformed_shape_hessians[k][d] -=
                          output_data.shape_gradients[first + d][k][n] *
                          mapping_data.jacobian_pushed_forward_grads[k][n];

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_hessians[first + d][k] =
                        fe_data.transformed_shape_hessians[k][d];

                  break;
                }
              case mapping_covariant:
                {
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    fe_data.untransformed_shape_hessian_tensors[k] =
                      fe_data.shape_grad_grads[i][k];

                  mapping.transform(
                    make_array_view(
                      fe_data.untransformed_shape_hessian_tensors),
                    mapping_covariant_hessian,
                    mapping_internal,
                    make_array_view(fe_data.transformed_shape_hessians));

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        for (unsigned int i = 0; i < spacedim; ++i)
                          for (unsigned int j = 0; j < spacedim; ++j)
                            {
                              fe_data.transformed_shape_hessians[k][d][i][j] -=
                                (output_data.shape_values(first + n, k) *
                                 mapping_data
                                   .jacobian_pushed_forward_2nd_derivatives
                                     [k][n][d][i][j]) +
                                (output_data.shape_gradients[first + d][k][n] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][i][j]) +
                                (output_data.shape_gradients[first + n][k][i] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][d][j]) +
                                (output_data.shape_gradients[first + n][k][j] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][i][d]);
                            }

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_hessians[first + d][k] =
                        fe_data.transformed_shape_hessians[k][d];

                  break;
                }
              case mapping_contravariant:
                {
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    fe_data.untransformed_shape_hessian_tensors[k] =
                      fe_data.shape_grad_grads[i][k];

                  mapping.transform(
                    make_array_view(
                      fe_data.untransformed_shape_hessian_tensors),
                    mapping_contravariant_hessian,
                    mapping_internal,
                    make_array_view(fe_data.transformed_shape_hessians));

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        for (unsigned int i = 0; i < spacedim; ++i)
                          for (unsigned int j = 0; j < spacedim; ++j)
                            {
                              fe_data.transformed_shape_hessians[k][d][i][j] +=
                                (output_data.shape_values(first + n, k) *
                                 mapping_data
                                   .jacobian_pushed_forward_2nd_derivatives
                                     [k][d][n][i][j]) +
                                (output_data.shape_gradients[first + n][k][i] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][d][n][j]) +
                                (output_data.shape_gradients[first + n][k][j] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][d][i][n]) -
                                (output_data.shape_gradients[first + d][k][n] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][i][j]);
                              for (unsigned int m = 0; m < spacedim; ++m)
                                fe_data
                                  .transformed_shape_hessians[k][d][i][j] -=
                                  (mapping_data
                                     .jacobian_pushed_forward_grads[k][d][i]
                                                                   [m] *
                                   mapping_data
                                     .jacobian_pushed_forward_grads[k][m][n]
                                                                   [j] *
                                   output_data.shape_values(first + n, k)) +
                                  (mapping_data
                                     .jacobian_pushed_forward_grads[k][d][m]
                                                                   [j] *
                                   mapping_data
                                     .jacobian_pushed_forward_grads[k][m][i]
                                                                   [n] *
                                   output_data.shape_values(first + n, k));
                            }

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_hessians[first + d][k] =
                        fe_data.transformed_shape_hessians[k][d];

                  break;
                }
              case mapping_raviart_thomas:
              case mapping_piola:
                {
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    fe_data.untransformed_shape_hessian_tensors[k] =
                      fe_data.shape_grad_grads[i][k];

                  mapping.transform(
                    make_array_view(
                      fe_data.untransformed_shape_hessian_tensors),
                    mapping_piola_hessian,
                    mapping_internal,
                    make_array_view(fe_data.transformed_shape_hessians));

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        for (unsigned int i = 0; i < spacedim; ++i)
                          for (unsigned int j = 0; j < spacedim; ++j)
                            {
                              fe_data.transformed_shape_hessians[k][d][i][j] +=
                                (output_data.shape_values(first + n, k) *
                                 mapping_data
                                   .jacobian_pushed_forward_2nd_derivatives
                                     [k][d][n][i][j]) +
                                (output_data.shape_gradients[first + n][k][i] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][d][n][j]) +
                                (output_data.shape_gradients[first + n][k][j] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][d][i][n]) -
                                (output_data.shape_gradients[first + d][k][n] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][i][j]);

                              fe_data.transformed_shape_hessians[k][d][i][j] -=
                                (output_data.shape_values(first + d, k) *
                                 mapping_data
                                   .jacobian_pushed_forward_2nd_derivatives
                                     [k][n][n][i][j]) +
                                (output_data.shape_gradients[first + d][k][i] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][n][j]) +
                                (output_data.shape_gradients[first + d][k][j] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][n][i]);

                              for (unsigned int m = 0; m < spacedim; ++m)
                                {
                                  fe_data
                                    .transformed_shape_hessians[k][d][i][j] -=
                                    (mapping_data
                                       .jacobian_pushed_forward_grads[k][d][i]
                                                                     [m] *
                                     mapping_data
                                       .jacobian_pushed_forward_grads[k][m][n]
                                                                     [j] *
                                     output_data.shape_values(first + n, k)) +
                                    (mapping_data
                                       .jacobian_pushed_forward_grads[k][d][m]
                                                                     [j] *
                                     mapping_data
                                       .jacobian_pushed_forward_grads[k][m][i]
                                                                     [n] *
                                     output_data.shape_values(first + n, k));

                                  fe_data
                                    .transformed_shape_hessians[k][d][i][j] +=
                                    (mapping_data
                                       .jacobian_pushed_forward_grads[k][n][i]
                                                                     [m] *
                                     mapping_data
                                       .jacobian_pushed_forward_grads[k][m][n]
                                                                     [j] *
                                     output_data.shape_values(first + d, k)) +
                                    (mapping_data
                                       .jacobian_pushed_forward_grads[k][n][m]
                                                                     [j] *
                                     mapping_data
                                       .jacobian_pushed_forward_grads[k][m][i]
                                                                     [n] *
                                     output_data.shape_values(first + d, k));
                                }
                            }

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_hessians[first + d][k] =
                        fe_data.sign_change[i] *
                        fe_data.transformed_shape_hessians[k][d];

                  break;
                }

              case mapping_nedelec:
                {
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    fe_data.untransformed_shape_hessian_tensors[k] =
                      fe_data.shape_grad_grads[i][k];

                  mapping.transform(
                    make_array_view(
                      fe_data.untransformed_shape_hessian_tensors),
                    mapping_covariant_hessian,
                    mapping_internal,
                    make_array_view(fe_data.transformed_shape_hessians));

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        for (unsigned int i = 0; i < spacedim; ++i)
                          for (unsigned int j = 0; j < spacedim; ++j)
                            {
                              fe_data.transformed_shape_hessians[k][d][i][j] -=
                                (output_data.shape_values(first + n, k) *
                                 mapping_data
                                   .jacobian_pushed_forward_2nd_derivatives
                                     [k][n][d][i][j]) +
                                (output_data.shape_gradients[first + d][k][n] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][i][j]) +
                                (output_data.shape_gradients[first + n][k][i] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][d][j]) +
                                (output_data.shape_gradients[first + n][k][j] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][i][d]);
                            }

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_hessians[first + d][k] =
                        fe_data.sign_change[i] *
                        fe_data.transformed_shape_hessians[k][d];

                  break;
                }

              default:
                Assert(false, ExcNotImplemented());
            }
        }

      // third derivatives are not implemented
      if (fe_data.update_each & update_3rd_derivatives &&
          ((cell_similarity != CellSimilarity::translation) ||
           ((mapping_kind == mapping_piola) ||
            (mapping_kind == mapping_raviart_thomas) ||
            (mapping_kind == mapping_nedelec))))
        {
          Assert(false, ExcNotImplemented())
        }
    }
}



template <int dim, int spacedim>
void
FE_PolyTensor<dim, spacedim>::fill_fe_face_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const Quadrature<dim - 1> &                                 quadrature,
  const Mapping<dim, spacedim> &                              mapping,
  const typename Mapping<dim, spacedim>::InternalDataBase &   mapping_internal,
  const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                     spacedim>
    &                                                            mapping_data,
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

  const unsigned int n_q_points = quadrature.size();
  // offset determines which data set
  // to take (all data sets for all
  // faces are stored contiguously)

  const typename QProjector<dim>::DataSetDescriptor offset =
    QProjector<dim>::DataSetDescriptor::face(face_no,
                                             cell->face_orientation(face_no),
                                             cell->face_flip(face_no),
                                             cell->face_rotation(face_no),
                                             n_q_points);

  // TODO: Size assertions

  // Create table with sign changes, due to the special structure of the RT
  // elements.
  // TODO: Preliminary hack to demonstrate the overall prinicple!

  // Compute eventual sign changes depending
  // on the neighborhood between two faces.
  std::fill(fe_data.sign_change.begin(), fe_data.sign_change.end(), 1.0);

  internal::FE_PolyTensor::get_face_sign_change_rt(cell,
                                                   *this,
                                                   this->mapping_kind,
                                                   fe_data.sign_change);

  internal::FE_PolyTensor::get_face_sign_change_nedelec(cell,
                                                        *this,
                                                        this->mapping_kind,
                                                        fe_data.sign_change);

  for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
    {
      const MappingKind mapping_kind = get_mapping_kind(i);

      const unsigned int first =
        output_data.shape_function_to_row_table[i * this->n_components() +
                                                this->get_nonzero_components(i)
                                                  .first_selected_component()];

      if (fe_data.update_each & update_values)
        {
          switch (mapping_kind)
            {
              case mapping_none:
                {
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_values(first + d, k) =
                        fe_data.shape_values[i][k + offset][d];
                  break;
                }

              case mapping_covariant:
              case mapping_contravariant:
                {
                  const ArrayView<Tensor<1, spacedim>>
                    transformed_shape_values =
                      make_array_view(fe_data.transformed_shape_values,
                                      offset,
                                      n_q_points);
                  mapping.transform(make_array_view(fe_data.shape_values,
                                                    i,
                                                    offset,
                                                    n_q_points),
                                    mapping_kind,
                                    mapping_internal,
                                    transformed_shape_values);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_values(first + d, k) =
                        transformed_shape_values[k][d];

                  break;
                }
              case mapping_raviart_thomas:
              case mapping_piola:
                {
                  const ArrayView<Tensor<1, spacedim>>
                    transformed_shape_values =
                      make_array_view(fe_data.transformed_shape_values,
                                      offset,
                                      n_q_points);
                  mapping.transform(make_array_view(fe_data.shape_values,
                                                    i,
                                                    offset,
                                                    n_q_points),
                                    mapping_piola,
                                    mapping_internal,
                                    transformed_shape_values);
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_values(first + d, k) =
                        fe_data.sign_change[i] * transformed_shape_values[k][d];
                  break;
                }

              case mapping_nedelec:
                {
                  const ArrayView<Tensor<1, spacedim>>
                    transformed_shape_values =
                      make_array_view(fe_data.transformed_shape_values,
                                      offset,
                                      n_q_points);
                  mapping.transform(make_array_view(fe_data.shape_values,
                                                    i,
                                                    offset,
                                                    n_q_points),
                                    mapping_covariant,
                                    mapping_internal,
                                    transformed_shape_values);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_values(first + d, k) =
                        fe_data.sign_change[i] * transformed_shape_values[k][d];

                  break;
                }

              default:
                Assert(false, ExcNotImplemented());
            }
        }

      if (fe_data.update_each & update_gradients)
        {
          switch (mapping_kind)
            {
              case mapping_none:
                {
                  const ArrayView<Tensor<2, spacedim>> transformed_shape_grads =
                    make_array_view(fe_data.transformed_shape_grads,
                                    offset,
                                    n_q_points);
                  mapping.transform(
                    make_array_view(fe_data.shape_grads, i, offset, n_q_points),
                    mapping_covariant,
                    mapping_internal,
                    transformed_shape_grads);
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_gradients[first + d][k] =
                        transformed_shape_grads[k][d];
                  break;
                }

              case mapping_covariant:
                {
                  const ArrayView<Tensor<2, spacedim>> transformed_shape_grads =
                    make_array_view(fe_data.transformed_shape_grads,
                                    offset,
                                    n_q_points);
                  mapping.transform(
                    make_array_view(fe_data.shape_grads, i, offset, n_q_points),
                    mapping_covariant_gradient,
                    mapping_internal,
                    transformed_shape_grads);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        transformed_shape_grads[k][d] -=
                          output_data.shape_values(first + n, k) *
                          mapping_data.jacobian_pushed_forward_grads[k][n][d];

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_gradients[first + d][k] =
                        transformed_shape_grads[k][d];
                  break;
                }

              case mapping_contravariant:
                {
                  const ArrayView<Tensor<2, spacedim>> transformed_shape_grads =
                    make_array_view(fe_data.transformed_shape_grads,
                                    offset,
                                    n_q_points);
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    fe_data.untransformed_shape_grads[k + offset] =
                      fe_data.shape_grads[i][k + offset];
                  mapping.transform(
                    make_array_view(fe_data.untransformed_shape_grads,
                                    offset,
                                    n_q_points),
                    mapping_contravariant_gradient,
                    mapping_internal,
                    transformed_shape_grads);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        transformed_shape_grads[k][d] +=
                          output_data.shape_values(first + n, k) *
                          mapping_data.jacobian_pushed_forward_grads[k][d][n];

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_gradients[first + d][k] =
                        transformed_shape_grads[k][d];

                  break;
                }

              case mapping_raviart_thomas:
              case mapping_piola:
                {
                  const ArrayView<Tensor<2, spacedim>> transformed_shape_grads =
                    make_array_view(fe_data.transformed_shape_grads,
                                    offset,
                                    n_q_points);
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    fe_data.untransformed_shape_grads[k + offset] =
                      fe_data.shape_grads[i][k + offset];
                  mapping.transform(
                    make_array_view(fe_data.untransformed_shape_grads,
                                    offset,
                                    n_q_points),
                    mapping_piola_gradient,
                    mapping_internal,
                    transformed_shape_grads);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        transformed_shape_grads[k][d] +=
                          (output_data.shape_values(first + n, k) *
                           mapping_data
                             .jacobian_pushed_forward_grads[k][d][n]) -
                          (output_data.shape_values(first + d, k) *
                           mapping_data.jacobian_pushed_forward_grads[k][n][n]);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_gradients[first + d][k] =
                        fe_data.sign_change[i] * transformed_shape_grads[k][d];

                  break;
                }

              case mapping_nedelec:
                {
                  // treat the gradients of
                  // this particular shape
                  // function at all
                  // q-points. if Dv is the
                  // gradient of the shape
                  // function on the unit
                  // cell, then
                  // (J^-T)Dv(J^-1) is the
                  // value we want to have on
                  // the real cell.
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    fe_data.untransformed_shape_grads[k + offset] =
                      fe_data.shape_grads[i][k + offset];

                  const ArrayView<Tensor<2, spacedim>> transformed_shape_grads =
                    make_array_view(fe_data.transformed_shape_grads,
                                    offset,
                                    n_q_points);
                  mapping.transform(
                    make_array_view(fe_data.untransformed_shape_grads,
                                    offset,
                                    n_q_points),
                    mapping_covariant_gradient,
                    mapping_internal,
                    transformed_shape_grads);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        transformed_shape_grads[k][d] -=
                          output_data.shape_values(first + n, k) *
                          mapping_data.jacobian_pushed_forward_grads[k][n][d];

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_gradients[first + d][k] =
                        fe_data.sign_change[i] * transformed_shape_grads[k][d];

                  break;
                }

              default:
                Assert(false, ExcNotImplemented());
            }
        }

      if (fe_data.update_each & update_hessians)
        {
          switch (mapping_kind)
            {
              case mapping_none:
                {
                  const ArrayView<Tensor<3, spacedim>>
                    transformed_shape_hessians =
                      make_array_view(fe_data.transformed_shape_hessians,
                                      offset,
                                      n_q_points);
                  mapping.transform(make_array_view(fe_data.shape_grad_grads,
                                                    i,
                                                    offset,
                                                    n_q_points),
                                    mapping_covariant_gradient,
                                    mapping_internal,
                                    transformed_shape_hessians);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        transformed_shape_hessians[k][d] -=
                          output_data.shape_gradients[first + d][k][n] *
                          mapping_data.jacobian_pushed_forward_grads[k][n];

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_hessians[first + d][k] =
                        transformed_shape_hessians[k][d];

                  break;
                }
              case mapping_covariant:
                {
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    fe_data.untransformed_shape_hessian_tensors[k + offset] =
                      fe_data.shape_grad_grads[i][k + offset];

                  const ArrayView<Tensor<3, spacedim>>
                    transformed_shape_hessians =
                      make_array_view(fe_data.transformed_shape_hessians,
                                      offset,
                                      n_q_points);
                  mapping.transform(
                    make_array_view(fe_data.untransformed_shape_hessian_tensors,
                                    offset,
                                    n_q_points),
                    mapping_covariant_hessian,
                    mapping_internal,
                    transformed_shape_hessians);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        for (unsigned int i = 0; i < spacedim; ++i)
                          for (unsigned int j = 0; j < spacedim; ++j)
                            {
                              transformed_shape_hessians[k][d][i][j] -=
                                (output_data.shape_values(first + n, k) *
                                 mapping_data
                                   .jacobian_pushed_forward_2nd_derivatives
                                     [k][n][d][i][j]) +
                                (output_data.shape_gradients[first + d][k][n] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][i][j]) +
                                (output_data.shape_gradients[first + n][k][i] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][d][j]) +
                                (output_data.shape_gradients[first + n][k][j] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][i][d]);
                            }

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_hessians[first + d][k] =
                        transformed_shape_hessians[k][d];

                  break;
                }

              case mapping_contravariant:
                {
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    fe_data.untransformed_shape_hessian_tensors[k + offset] =
                      fe_data.shape_grad_grads[i][k + offset];

                  const ArrayView<Tensor<3, spacedim>>
                    transformed_shape_hessians =
                      make_array_view(fe_data.transformed_shape_hessians,
                                      offset,
                                      n_q_points);
                  mapping.transform(
                    make_array_view(fe_data.untransformed_shape_hessian_tensors,
                                    offset,
                                    n_q_points),
                    mapping_contravariant_hessian,
                    mapping_internal,
                    transformed_shape_hessians);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        for (unsigned int i = 0; i < spacedim; ++i)
                          for (unsigned int j = 0; j < spacedim; ++j)
                            {
                              transformed_shape_hessians[k][d][i][j] +=
                                (output_data.shape_values(first + n, k) *
                                 mapping_data
                                   .jacobian_pushed_forward_2nd_derivatives
                                     [k][d][n][i][j]) +
                                (output_data.shape_gradients[first + n][k][i] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][d][n][j]) +
                                (output_data.shape_gradients[first + n][k][j] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][d][i][n]) -
                                (output_data.shape_gradients[first + d][k][n] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][i][j]);
                              for (unsigned int m = 0; m < spacedim; ++m)
                                transformed_shape_hessians[k][d][i][j] -=
                                  (mapping_data
                                     .jacobian_pushed_forward_grads[k][d][i]
                                                                   [m] *
                                   mapping_data
                                     .jacobian_pushed_forward_grads[k][m][n]
                                                                   [j] *
                                   output_data.shape_values(first + n, k)) +
                                  (mapping_data
                                     .jacobian_pushed_forward_grads[k][d][m]
                                                                   [j] *
                                   mapping_data
                                     .jacobian_pushed_forward_grads[k][m][i]
                                                                   [n] *
                                   output_data.shape_values(first + n, k));
                            }

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_hessians[first + d][k] =
                        transformed_shape_hessians[k][d];

                  break;
                }

              case mapping_raviart_thomas:
              case mapping_piola:
                {
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    fe_data.untransformed_shape_hessian_tensors[k + offset] =
                      fe_data.shape_grad_grads[i][k + offset];

                  const ArrayView<Tensor<3, spacedim>>
                    transformed_shape_hessians =
                      make_array_view(fe_data.transformed_shape_hessians,
                                      offset,
                                      n_q_points);
                  mapping.transform(
                    make_array_view(fe_data.untransformed_shape_hessian_tensors,
                                    offset,
                                    n_q_points),
                    mapping_piola_hessian,
                    mapping_internal,
                    transformed_shape_hessians);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        for (unsigned int i = 0; i < spacedim; ++i)
                          for (unsigned int j = 0; j < spacedim; ++j)
                            {
                              transformed_shape_hessians[k][d][i][j] +=
                                (output_data.shape_values(first + n, k) *
                                 mapping_data
                                   .jacobian_pushed_forward_2nd_derivatives
                                     [k][d][n][i][j]) +
                                (output_data.shape_gradients[first + n][k][i] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][d][n][j]) +
                                (output_data.shape_gradients[first + n][k][j] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][d][i][n]) -
                                (output_data.shape_gradients[first + d][k][n] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][i][j]);

                              transformed_shape_hessians[k][d][i][j] -=
                                (output_data.shape_values(first + d, k) *
                                 mapping_data
                                   .jacobian_pushed_forward_2nd_derivatives
                                     [k][n][n][i][j]) +
                                (output_data.shape_gradients[first + d][k][i] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][n][j]) +
                                (output_data.shape_gradients[first + d][k][j] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][n][i]);

                              for (unsigned int m = 0; m < spacedim; ++m)
                                {
                                  transformed_shape_hessians[k][d][i][j] -=
                                    (mapping_data
                                       .jacobian_pushed_forward_grads[k][d][i]
                                                                     [m] *
                                     mapping_data
                                       .jacobian_pushed_forward_grads[k][m][n]
                                                                     [j] *
                                     output_data.shape_values(first + n, k)) +
                                    (mapping_data
                                       .jacobian_pushed_forward_grads[k][d][m]
                                                                     [j] *
                                     mapping_data
                                       .jacobian_pushed_forward_grads[k][m][i]
                                                                     [n] *
                                     output_data.shape_values(first + n, k));

                                  transformed_shape_hessians[k][d][i][j] +=
                                    (mapping_data
                                       .jacobian_pushed_forward_grads[k][n][i]
                                                                     [m] *
                                     mapping_data
                                       .jacobian_pushed_forward_grads[k][m][n]
                                                                     [j] *
                                     output_data.shape_values(first + d, k)) +
                                    (mapping_data
                                       .jacobian_pushed_forward_grads[k][n][m]
                                                                     [j] *
                                     mapping_data
                                       .jacobian_pushed_forward_grads[k][m][i]
                                                                     [n] *
                                     output_data.shape_values(first + d, k));
                                }
                            }

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_hessians[first + d][k] =
                        fe_data.sign_change[i] *
                        transformed_shape_hessians[k][d];

                  break;
                }

              case mapping_nedelec:
                {
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    fe_data.untransformed_shape_hessian_tensors[k + offset] =
                      fe_data.shape_grad_grads[i][k + offset];

                  const ArrayView<Tensor<3, spacedim>>
                    transformed_shape_hessians =
                      make_array_view(fe_data.transformed_shape_hessians,
                                      offset,
                                      n_q_points);
                  mapping.transform(
                    make_array_view(fe_data.untransformed_shape_hessian_tensors,
                                    offset,
                                    n_q_points),
                    mapping_covariant_hessian,
                    mapping_internal,
                    transformed_shape_hessians);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        for (unsigned int i = 0; i < spacedim; ++i)
                          for (unsigned int j = 0; j < spacedim; ++j)
                            {
                              transformed_shape_hessians[k][d][i][j] -=
                                (output_data.shape_values(first + n, k) *
                                 mapping_data
                                   .jacobian_pushed_forward_2nd_derivatives
                                     [k][n][d][i][j]) +
                                (output_data.shape_gradients[first + d][k][n] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][i][j]) +
                                (output_data.shape_gradients[first + n][k][i] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][d][j]) +
                                (output_data.shape_gradients[first + n][k][j] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][i][d]);
                            }

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_hessians[first + d][k] =
                        fe_data.sign_change[i] *
                        transformed_shape_hessians[k][d];

                  break;
                }

              default:
                Assert(false, ExcNotImplemented());
            }
        }

      // third derivatives are not implemented
      if (fe_data.update_each & update_3rd_derivatives)
        {
          Assert(false, ExcNotImplemented())
        }
    }
}



template <int dim, int spacedim>
void
FE_PolyTensor<dim, spacedim>::fill_fe_subface_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const unsigned int                                          sub_no,
  const Quadrature<dim - 1> &                                 quadrature,
  const Mapping<dim, spacedim> &                              mapping,
  const typename Mapping<dim, spacedim>::InternalDataBase &   mapping_internal,
  const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                     spacedim>
    &                                                            mapping_data,
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

  const unsigned int n_q_points = quadrature.size();

  // offset determines which data set
  // to take (all data sets for all
  // sub-faces are stored contiguously)
  const typename QProjector<dim>::DataSetDescriptor offset =
    QProjector<dim>::DataSetDescriptor::subface(face_no,
                                                sub_no,
                                                cell->face_orientation(face_no),
                                                cell->face_flip(face_no),
                                                cell->face_rotation(face_no),
                                                n_q_points,
                                                cell->subface_case(face_no));

  //   Assert(mapping_kind == independent
  //       || ( mapping_kind == independent_on_cartesian
  //            && dynamic_cast<const MappingCartesian<dim>*>(&mapping) != 0),
  //       ExcNotImplemented());
  // TODO: Size assertions

  // TODO: Sign change for the face DoFs!

  // Compute eventual sign changes depending
  // on the neighborhood between two faces.
  std::fill(fe_data.sign_change.begin(), fe_data.sign_change.end(), 1.0);

  internal::FE_PolyTensor::get_face_sign_change_rt(cell,
                                                   *this,
                                                   this->mapping_kind,
                                                   fe_data.sign_change);

  internal::FE_PolyTensor::get_face_sign_change_nedelec(cell,
                                                        *this,
                                                        this->mapping_kind,
                                                        fe_data.sign_change);

  for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
    {
      const MappingKind mapping_kind = get_mapping_kind(i);

      const unsigned int first =
        output_data.shape_function_to_row_table[i * this->n_components() +
                                                this->get_nonzero_components(i)
                                                  .first_selected_component()];

      if (fe_data.update_each & update_values)
        {
          switch (mapping_kind)
            {
              case mapping_none:
                {
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_values(first + d, k) =
                        fe_data.shape_values[i][k + offset][d];
                  break;
                }

              case mapping_covariant:
              case mapping_contravariant:
                {
                  const ArrayView<Tensor<1, spacedim>>
                    transformed_shape_values =
                      make_array_view(fe_data.transformed_shape_values,
                                      offset,
                                      n_q_points);
                  mapping.transform(make_array_view(fe_data.shape_values,
                                                    i,
                                                    offset,
                                                    n_q_points),
                                    mapping_kind,
                                    mapping_internal,
                                    transformed_shape_values);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_values(first + d, k) =
                        transformed_shape_values[k][d];

                  break;
                }

              case mapping_raviart_thomas:
              case mapping_piola:
                {
                  const ArrayView<Tensor<1, spacedim>>
                    transformed_shape_values =
                      make_array_view(fe_data.transformed_shape_values,
                                      offset,
                                      n_q_points);

                  mapping.transform(make_array_view(fe_data.shape_values,
                                                    i,
                                                    offset,
                                                    n_q_points),
                                    mapping_piola,
                                    mapping_internal,
                                    transformed_shape_values);
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_values(first + d, k) =
                        fe_data.sign_change[i] * transformed_shape_values[k][d];
                  break;
                }

              case mapping_nedelec:
                {
                  const ArrayView<Tensor<1, spacedim>>
                    transformed_shape_values =
                      make_array_view(fe_data.transformed_shape_values,
                                      offset,
                                      n_q_points);

                  mapping.transform(make_array_view(fe_data.shape_values,
                                                    i,
                                                    offset,
                                                    n_q_points),
                                    mapping_covariant,
                                    mapping_internal,
                                    transformed_shape_values);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_values(first + d, k) =
                        fe_data.sign_change[i] * transformed_shape_values[k][d];

                  break;
                }

              default:
                Assert(false, ExcNotImplemented());
            }
        }

      if (fe_data.update_each & update_gradients)
        {
          const ArrayView<Tensor<2, spacedim>> transformed_shape_grads =
            make_array_view(fe_data.transformed_shape_grads,
                            offset,
                            n_q_points);
          switch (mapping_kind)
            {
              case mapping_none:
                {
                  mapping.transform(
                    make_array_view(fe_data.shape_grads, i, offset, n_q_points),
                    mapping_covariant,
                    mapping_internal,
                    transformed_shape_grads);
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_gradients[first + d][k] =
                        transformed_shape_grads[k][d];
                  break;
                }

              case mapping_covariant:
                {
                  mapping.transform(
                    make_array_view(fe_data.shape_grads, i, offset, n_q_points),
                    mapping_covariant_gradient,
                    mapping_internal,
                    transformed_shape_grads);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        transformed_shape_grads[k][d] -=
                          output_data.shape_values(first + n, k) *
                          mapping_data.jacobian_pushed_forward_grads[k][n][d];

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_gradients[first + d][k] =
                        transformed_shape_grads[k][d];

                  break;
                }

              case mapping_contravariant:
                {
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    fe_data.untransformed_shape_grads[k + offset] =
                      fe_data.shape_grads[i][k + offset];

                  mapping.transform(
                    make_array_view(fe_data.untransformed_shape_grads,
                                    offset,
                                    n_q_points),
                    mapping_contravariant_gradient,
                    mapping_internal,
                    transformed_shape_grads);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        transformed_shape_grads[k][d] +=
                          output_data.shape_values(first + n, k) *
                          mapping_data.jacobian_pushed_forward_grads[k][d][n];

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_gradients[first + d][k] =
                        transformed_shape_grads[k][d];

                  break;
                }

              case mapping_raviart_thomas:
              case mapping_piola:
                {
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    fe_data.untransformed_shape_grads[k + offset] =
                      fe_data.shape_grads[i][k + offset];

                  mapping.transform(
                    make_array_view(fe_data.untransformed_shape_grads,
                                    offset,
                                    n_q_points),
                    mapping_piola_gradient,
                    mapping_internal,
                    transformed_shape_grads);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        transformed_shape_grads[k][d] +=
                          (output_data.shape_values(first + n, k) *
                           mapping_data
                             .jacobian_pushed_forward_grads[k][d][n]) -
                          (output_data.shape_values(first + d, k) *
                           mapping_data.jacobian_pushed_forward_grads[k][n][n]);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_gradients[first + d][k] =
                        fe_data.sign_change[i] * transformed_shape_grads[k][d];

                  break;
                }

              case mapping_nedelec:
                {
                  // this particular shape
                  // function at all
                  // q-points. if Dv is the
                  // gradient of the shape
                  // function on the unit
                  // cell, then
                  // (J^-T)Dv(J^-1) is the
                  // value we want to have on
                  // the real cell.
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    fe_data.untransformed_shape_grads[k + offset] =
                      fe_data.shape_grads[i][k + offset];

                  mapping.transform(
                    make_array_view(fe_data.untransformed_shape_grads,
                                    offset,
                                    n_q_points),
                    mapping_covariant_gradient,
                    mapping_internal,
                    transformed_shape_grads);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        transformed_shape_grads[k][d] -=
                          output_data.shape_values(first + n, k) *
                          mapping_data.jacobian_pushed_forward_grads[k][n][d];

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_gradients[first + d][k] =
                        fe_data.sign_change[i] * transformed_shape_grads[k][d];

                  break;
                }

              default:
                Assert(false, ExcNotImplemented());
            }
        }

      if (fe_data.update_each & update_hessians)
        {
          switch (mapping_kind)
            {
              case mapping_none:
                {
                  const ArrayView<Tensor<3, spacedim>>
                    transformed_shape_hessians =
                      make_array_view(fe_data.transformed_shape_hessians,
                                      offset,
                                      n_q_points);
                  mapping.transform(make_array_view(fe_data.shape_grad_grads,
                                                    i,
                                                    offset,
                                                    n_q_points),
                                    mapping_covariant_gradient,
                                    mapping_internal,
                                    transformed_shape_hessians);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        transformed_shape_hessians[k][d] -=
                          output_data.shape_gradients[first + d][k][n] *
                          mapping_data.jacobian_pushed_forward_grads[k][n];

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_hessians[first + d][k] =
                        transformed_shape_hessians[k][d];

                  break;
                }
              case mapping_covariant:
                {
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    fe_data.untransformed_shape_hessian_tensors[k + offset] =
                      fe_data.shape_grad_grads[i][k + offset];

                  const ArrayView<Tensor<3, spacedim>>
                    transformed_shape_hessians =
                      make_array_view(fe_data.transformed_shape_hessians,
                                      offset,
                                      n_q_points);
                  mapping.transform(
                    make_array_view(fe_data.untransformed_shape_hessian_tensors,
                                    offset,
                                    n_q_points),
                    mapping_covariant_hessian,
                    mapping_internal,
                    transformed_shape_hessians);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        for (unsigned int i = 0; i < spacedim; ++i)
                          for (unsigned int j = 0; j < spacedim; ++j)
                            {
                              transformed_shape_hessians[k][d][i][j] -=
                                (output_data.shape_values(first + n, k) *
                                 mapping_data
                                   .jacobian_pushed_forward_2nd_derivatives
                                     [k][n][d][i][j]) +
                                (output_data.shape_gradients[first + d][k][n] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][i][j]) +
                                (output_data.shape_gradients[first + n][k][i] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][d][j]) +
                                (output_data.shape_gradients[first + n][k][j] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][i][d]);
                            }

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_hessians[first + d][k] =
                        transformed_shape_hessians[k][d];

                  break;
                }

              case mapping_contravariant:
                {
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    fe_data.untransformed_shape_hessian_tensors[k + offset] =
                      fe_data.shape_grad_grads[i][k + offset];

                  const ArrayView<Tensor<3, spacedim>>
                    transformed_shape_hessians =
                      make_array_view(fe_data.transformed_shape_hessians,
                                      offset,
                                      n_q_points);
                  mapping.transform(
                    make_array_view(fe_data.untransformed_shape_hessian_tensors,
                                    offset,
                                    n_q_points),
                    mapping_contravariant_hessian,
                    mapping_internal,
                    transformed_shape_hessians);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        for (unsigned int i = 0; i < spacedim; ++i)
                          for (unsigned int j = 0; j < spacedim; ++j)
                            {
                              transformed_shape_hessians[k][d][i][j] +=
                                (output_data.shape_values(first + n, k) *
                                 mapping_data
                                   .jacobian_pushed_forward_2nd_derivatives
                                     [k][d][n][i][j]) +
                                (output_data.shape_gradients[first + n][k][i] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][d][n][j]) +
                                (output_data.shape_gradients[first + n][k][j] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][d][i][n]) -
                                (output_data.shape_gradients[first + d][k][n] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][i][j]);
                              for (unsigned int m = 0; m < spacedim; ++m)
                                transformed_shape_hessians[k][d][i][j] -=
                                  (mapping_data
                                     .jacobian_pushed_forward_grads[k][d][i]
                                                                   [m] *
                                   mapping_data
                                     .jacobian_pushed_forward_grads[k][m][n]
                                                                   [j] *
                                   output_data.shape_values(first + n, k)) +
                                  (mapping_data
                                     .jacobian_pushed_forward_grads[k][d][m]
                                                                   [j] *
                                   mapping_data
                                     .jacobian_pushed_forward_grads[k][m][i]
                                                                   [n] *
                                   output_data.shape_values(first + n, k));
                            }

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_hessians[first + d][k] =
                        transformed_shape_hessians[k][d];

                  break;
                }

              case mapping_raviart_thomas:
              case mapping_piola:
                {
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    fe_data.untransformed_shape_hessian_tensors[k + offset] =
                      fe_data.shape_grad_grads[i][k + offset];

                  const ArrayView<Tensor<3, spacedim>>
                    transformed_shape_hessians =
                      make_array_view(fe_data.transformed_shape_hessians,
                                      offset,
                                      n_q_points);
                  mapping.transform(
                    make_array_view(fe_data.untransformed_shape_hessian_tensors,
                                    offset,
                                    n_q_points),
                    mapping_piola_hessian,
                    mapping_internal,
                    transformed_shape_hessians);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        for (unsigned int i = 0; i < spacedim; ++i)
                          for (unsigned int j = 0; j < spacedim; ++j)
                            {
                              transformed_shape_hessians[k][d][i][j] +=
                                (output_data.shape_values(first + n, k) *
                                 mapping_data
                                   .jacobian_pushed_forward_2nd_derivatives
                                     [k][d][n][i][j]) +
                                (output_data.shape_gradients[first + n][k][i] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][d][n][j]) +
                                (output_data.shape_gradients[first + n][k][j] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][d][i][n]) -
                                (output_data.shape_gradients[first + d][k][n] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][i][j]);

                              transformed_shape_hessians[k][d][i][j] -=
                                (output_data.shape_values(first + d, k) *
                                 mapping_data
                                   .jacobian_pushed_forward_2nd_derivatives
                                     [k][n][n][i][j]) +
                                (output_data.shape_gradients[first + d][k][i] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][n][j]) +
                                (output_data.shape_gradients[first + d][k][j] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][n][i]);
                              for (unsigned int m = 0; m < spacedim; ++m)
                                {
                                  transformed_shape_hessians[k][d][i][j] -=
                                    (mapping_data
                                       .jacobian_pushed_forward_grads[k][d][i]
                                                                     [m] *
                                     mapping_data
                                       .jacobian_pushed_forward_grads[k][m][n]
                                                                     [j] *
                                     output_data.shape_values(first + n, k)) +
                                    (mapping_data
                                       .jacobian_pushed_forward_grads[k][d][m]
                                                                     [j] *
                                     mapping_data
                                       .jacobian_pushed_forward_grads[k][m][i]
                                                                     [n] *
                                     output_data.shape_values(first + n, k));

                                  transformed_shape_hessians[k][d][i][j] +=
                                    (mapping_data
                                       .jacobian_pushed_forward_grads[k][n][i]
                                                                     [m] *
                                     mapping_data
                                       .jacobian_pushed_forward_grads[k][m][n]
                                                                     [j] *
                                     output_data.shape_values(first + d, k)) +
                                    (mapping_data
                                       .jacobian_pushed_forward_grads[k][n][m]
                                                                     [j] *
                                     mapping_data
                                       .jacobian_pushed_forward_grads[k][m][i]
                                                                     [n] *
                                     output_data.shape_values(first + d, k));
                                }
                            }

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_hessians[first + d][k] =
                        fe_data.sign_change[i] *
                        transformed_shape_hessians[k][d];

                  break;
                }

              case mapping_nedelec:
                {
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    fe_data.untransformed_shape_hessian_tensors[k + offset] =
                      fe_data.shape_grad_grads[i][k + offset];

                  const ArrayView<Tensor<3, spacedim>>
                    transformed_shape_hessians =
                      make_array_view(fe_data.transformed_shape_hessians,
                                      offset,
                                      n_q_points);
                  mapping.transform(
                    make_array_view(fe_data.untransformed_shape_hessian_tensors,
                                    offset,
                                    n_q_points),
                    mapping_covariant_hessian,
                    mapping_internal,
                    transformed_shape_hessians);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < spacedim; ++d)
                      for (unsigned int n = 0; n < spacedim; ++n)
                        for (unsigned int i = 0; i < spacedim; ++i)
                          for (unsigned int j = 0; j < spacedim; ++j)
                            {
                              transformed_shape_hessians[k][d][i][j] -=
                                (output_data.shape_values(first + n, k) *
                                 mapping_data
                                   .jacobian_pushed_forward_2nd_derivatives
                                     [k][n][d][i][j]) +
                                (output_data.shape_gradients[first + d][k][n] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][i][j]) +
                                (output_data.shape_gradients[first + n][k][i] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][d][j]) +
                                (output_data.shape_gradients[first + n][k][j] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[k][n][i][d]);
                            }

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_hessians[first + d][k] =
                        fe_data.sign_change[i] *
                        transformed_shape_hessians[k][d];

                  break;
                }

              default:
                Assert(false, ExcNotImplemented());
            }
        }

      // third derivatives are not implemented
      if (fe_data.update_each & update_3rd_derivatives)
        {
          Assert(false, ExcNotImplemented())
        }
    }
}



template <int dim, int spacedim>
UpdateFlags
FE_PolyTensor<dim, spacedim>::requires_update_flags(
  const UpdateFlags flags) const
{
  UpdateFlags out = update_default;

  for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
    {
      const MappingKind mapping_kind = get_mapping_kind(i);

      switch (mapping_kind)
        {
          case mapping_none:
            {
              if (flags & update_values)
                out |= update_values;

              if (flags & update_gradients)
                out |= update_gradients | update_values |
                       update_jacobian_pushed_forward_grads;

              if (flags & update_hessians)
                out |= update_hessians | update_values | update_gradients |
                       update_jacobian_pushed_forward_grads |
                       update_jacobian_pushed_forward_2nd_derivatives;
              break;
            }
          case mapping_raviart_thomas:
          case mapping_piola:
            {
              if (flags & update_values)
                out |= update_values | update_piola;

              if (flags & update_gradients)
                out |= update_gradients | update_values | update_piola |
                       update_jacobian_pushed_forward_grads |
                       update_covariant_transformation |
                       update_contravariant_transformation;

              if (flags & update_hessians)
                out |= update_hessians | update_piola | update_values |
                       update_gradients | update_jacobian_pushed_forward_grads |
                       update_jacobian_pushed_forward_2nd_derivatives |
                       update_covariant_transformation;

              break;
            }


          case mapping_contravariant:
            {
              if (flags & update_values)
                out |= update_values | update_piola;

              if (flags & update_gradients)
                out |= update_gradients | update_values |
                       update_jacobian_pushed_forward_grads |
                       update_covariant_transformation |
                       update_contravariant_transformation;

              if (flags & update_hessians)
                out |= update_hessians | update_piola | update_values |
                       update_gradients | update_jacobian_pushed_forward_grads |
                       update_jacobian_pushed_forward_2nd_derivatives |
                       update_covariant_transformation;

              break;
            }

          case mapping_nedelec:
          case mapping_covariant:
            {
              if (flags & update_values)
                out |= update_values | update_covariant_transformation;

              if (flags & update_gradients)
                out |= update_gradients | update_values |
                       update_jacobian_pushed_forward_grads |
                       update_covariant_transformation;

              if (flags & update_hessians)
                out |= update_hessians | update_values | update_gradients |
                       update_jacobian_pushed_forward_grads |
                       update_jacobian_pushed_forward_2nd_derivatives |
                       update_covariant_transformation;

              break;
            }

          default:
            {
              Assert(false, ExcNotImplemented());
            }
        }
    }

  return out;
}


// explicit instantiations
#include "fe_poly_tensor.inst"


DEAL_II_NAMESPACE_CLOSE
