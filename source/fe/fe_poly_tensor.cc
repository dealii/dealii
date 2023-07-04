// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2023 by the deal.II authors
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

#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_orientation.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace FE_PolyTensor
  {
    namespace
    {
      template <int spacedim>
      void
      get_dof_sign_change_h_div(
        const typename dealii::Triangulation<1, spacedim>::cell_iterator &,
        const FiniteElement<1, spacedim> &,
        const std::vector<MappingKind> &,
        std::vector<double> &)
      {
        // Nothing to do in 1d.
      }



      // TODO: This function is not a consistent fix of the orientation issue
      // like in 3d. It is rather kept not to break legacy behavior in 2d but
      // should be replaced. See also the implementation of
      // FE_RaviartThomas<dim>::initialize_quad_dof_index_permutation_and_sign_change()
      // or other H(div) conforming elements such as FE_ABF<dim> and
      // FE_BDM<dim>.
      template <int spacedim>
      void
      get_dof_sign_change_h_div(
        const typename dealii::Triangulation<2, spacedim>::cell_iterator &cell,
        const FiniteElement<2, spacedim> &                                fe,
        const std::vector<MappingKind> &mapping_kind,
        std::vector<double> &           face_sign)
      {
        const unsigned int dim = 2;
        // const unsigned int spacedim = 2;

        const CellId this_cell_id = cell->id();

        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
          {
            typename dealii::Triangulation<dim, spacedim>::face_iterator face =
              cell->face(f);
            if (!face->at_boundary())
              {
                const unsigned int nn = cell->neighbor_face_no(f);
                const typename dealii::Triangulation<dim,
                                                     spacedim>::cell_iterator
                             neighbor_cell_at_face = cell->neighbor(f);
                const CellId neigbor_cell_id = neighbor_cell_at_face->id();

                // Only fix sign if the orientation is opposite and only do so
                // on the face dofs on the cell with smaller cell_id.
                if (((nn + f) % 2 == 0) && this_cell_id < neigbor_cell_id)
                  for (unsigned int j = 0; j < fe.n_dofs_per_face(f); ++j)
                    {
                      const unsigned int cell_j = fe.face_to_cell_index(j, f);

                      Assert(f * fe.n_dofs_per_face(f) + j < face_sign.size(),
                             ExcInternalError());
                      Assert(mapping_kind.size() == 1 ||
                               cell_j < mapping_kind.size(),
                             ExcInternalError());

                      // TODO: This is probably only going to work for those
                      // elements for which all dofs are face dofs
                      if ((mapping_kind.size() > 1 ?
                             mapping_kind[cell_j] :
                             mapping_kind[0]) == mapping_raviart_thomas)
                        face_sign[f * fe.n_dofs_per_face(f) + j] = -1.0;
                    }
              }
          }
      }



      template <int spacedim>
      void
      get_dof_sign_change_h_div(
        const typename dealii::Triangulation<3, spacedim>::cell_iterator
          & /*cell*/,
        const FiniteElement<3, spacedim> & /*fe*/,
        const std::vector<MappingKind> & /*mapping_kind*/,
        std::vector<double> & /*face_sign*/)
      {
        // Nothing to do. In 3d we take care of it through the
        // adjust_quad_dof_sign_for_face_orientation_table
      }

      template <int spacedim>
      void
      get_dof_sign_change_nedelec(
        const typename dealii::Triangulation<1, spacedim>::cell_iterator
          & /*cell*/,
        const FiniteElement<1, spacedim> & /*fe*/,
        const std::vector<MappingKind> & /*mapping_kind*/,
        std::vector<double> & /*line_dof_sign*/)
      {
        // nothing to do in 1d
      }



      template <int spacedim>
      void
      get_dof_sign_change_nedelec(
        const typename dealii::Triangulation<2, spacedim>::cell_iterator &cell,
        const FiniteElement<2, spacedim> & /*fe*/,
        const std::vector<MappingKind> &mapping_kind,
        std::vector<double> &           line_dof_sign)
      {
        const unsigned int dim = 2;
        // TODO: This fixes only lowest order
        for (unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell; ++l)
          if (!(cell->line_orientation(l)) &&
              mapping_kind[0] == mapping_nedelec)
            line_dof_sign[l] = -1.0;
      }


      template <int spacedim>
      void
      get_dof_sign_change_nedelec(
        const typename dealii::Triangulation<3, spacedim>::cell_iterator &cell,
        const FiniteElement<3, spacedim> & /*fe*/,
        const std::vector<MappingKind> &mapping_kind,
        std::vector<double> &           line_dof_sign)
      {
        const unsigned int dim = 3;
        // TODO: This is probably only going to work for those elements for
        // which all dofs are face dofs
        for (unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell; ++l)
          if (!(cell->line_orientation(l)) &&
              mapping_kind[0] == mapping_nedelec)
            line_dof_sign[l] = -1.0;
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

  if (dim == 3)
    {
      adjust_quad_dof_sign_for_face_orientation_table.resize(
        this->n_unique_2d_subobjects());

      for (unsigned int f = 0; f < this->n_unique_2d_subobjects(); ++f)
        {
          adjust_quad_dof_sign_for_face_orientation_table[f] =
            Table<2, bool>(this->n_dofs_per_quad(f),
                           this->reference_cell().n_face_orientations(f));
          adjust_quad_dof_sign_for_face_orientation_table[f].fill(false);
        }
    }
}



template <int dim, int spacedim>
FE_PolyTensor<dim, spacedim>::FE_PolyTensor(const FE_PolyTensor &fe)
  : FiniteElement<dim, spacedim>(fe)
  , mapping_kind(fe.mapping_kind)
  , adjust_quad_dof_sign_for_face_orientation_table(
      fe.adjust_quad_dof_sign_for_face_orientation_table)
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
bool
FE_PolyTensor<dim, spacedim>::adjust_quad_dof_sign_for_face_orientation(
  const unsigned int index,
  const unsigned int face,
  const bool         face_orientation,
  const bool         face_flip,
  const bool         face_rotation) const
{
  // do nothing in 1d and 2d
  if (dim < 3)
    return false;

  // The exception are discontinuous
  // elements for which there should be no
  // face dofs anyway (i.e. dofs_per_quad==0
  // in 3d), so we don't need the table, but
  // the function should also not have been
  // called
  AssertIndexRange(index, this->n_dofs_per_quad(face));
  Assert(adjust_quad_dof_sign_for_face_orientation_table
             [this->n_unique_2d_subobjects() == 1 ? 0 : face]
               .n_elements() ==
           this->reference_cell().n_face_orientations(face) *
             this->n_dofs_per_quad(face),
         ExcInternalError());

  return adjust_quad_dof_sign_for_face_orientation_table
    [this->n_unique_2d_subobjects() == 1 ? 0 : face](
      index,
      internal::combined_face_orientation(face_orientation,
                                          face_rotation,
                                          face_flip));
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
  AssertIndexRange(i, this->n_dofs_per_cell());
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
  AssertIndexRange(i, this->n_dofs_per_cell());
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
  AssertIndexRange(i, this->n_dofs_per_cell());
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
  const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
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
           fe_data.shape_values.size()[0] == this->n_dofs_per_cell(),
         ExcDimensionMismatch(fe_data.shape_values.size()[0],
                              this->n_dofs_per_cell()));
  Assert(!(fe_data.update_each & update_values) ||
           fe_data.shape_values.size()[1] == n_q_points,
         ExcDimensionMismatch(fe_data.shape_values.size()[1], n_q_points));

  // TODO: The dof_sign_change only affects Nedelec elements and is not the
  // correct thing on complicated meshes for higher order Nedelec elements.
  // Something similar to FE_Q should be done to permute dofs and to change the
  // dof signs. A static way using tables (as done in the RaviartThomas<dim>
  // class) is preferable.
  std::fill(fe_data.dof_sign_change.begin(),
            fe_data.dof_sign_change.end(),
            1.0);
  internal::FE_PolyTensor::get_dof_sign_change_nedelec(cell,
                                                       *this,
                                                       this->mapping_kind,
                                                       fe_data.dof_sign_change);

  // TODO: This, similarly to the Nedelec case, is just a legacy function in 2d
  // and affects only face_dofs of H(div) conformal FEs. It does nothing in 1d.
  // Also nothing in 3d since we take care of it by using the
  // adjust_quad_dof_sign_for_face_orientation_table.
  internal::FE_PolyTensor::get_dof_sign_change_h_div(cell,
                                                     *this,
                                                     this->mapping_kind,
                                                     fe_data.dof_sign_change);

  // What is the first dof_index on a quad?
  const unsigned int first_quad_index = this->get_first_quad_index();
  // How many dofs per quad and how many quad dofs do we have at all?
  const unsigned int n_dofs_per_quad = this->n_dofs_per_quad();
  const unsigned int n_quad_dofs =
    n_dofs_per_quad * GeometryInfo<dim>::faces_per_cell;

  for (unsigned int dof_index = 0; dof_index < this->n_dofs_per_cell();
       ++dof_index)
    {
      /*
       * This assumes that the dofs are ordered by first vertices, lines, quads
       * and volume dofs. Note that in 2d this always gives false.
       */
      const bool is_quad_dof =
        (dim == 2 ? false :
                    (first_quad_index <= dof_index) &&
                      (dof_index < first_quad_index + n_quad_dofs));

      // TODO: This hack is not pretty and it is only here to handle the 2d
      // case and the Nedelec legacy case. In 2d dof_sign of a face_dof is never
      // handled by the
      // >>if(is_quad_dof){...}<< but still a possible dof sign change must be
      // handled, also for line_dofs in 3d such as in Nedelec. In these cases
      // this is encoded in the array fe_data.dof_sign_change[dof_index]. In 3d
      // it is handles with a table. This array is allocated in
      // fe_poly_tensor.h.
      double dof_sign = 1.0;
      // under some circumstances fe_data.dof_sign_change is not allocated
      if (fe_data.update_each & update_values)
        dof_sign = fe_data.dof_sign_change[dof_index];

      if (is_quad_dof)
        {
          /*
           * Find the face belonging to this dof_index. This is integer
           * division.
           */
          const unsigned int face_index_from_dof_index =
            (dof_index - first_quad_index) / (n_dofs_per_quad);

          const unsigned int local_quad_dof_index = dof_index % n_dofs_per_quad;

          // Correct the dof_sign if necessary
          if (adjust_quad_dof_sign_for_face_orientation(
                local_quad_dof_index,
                face_index_from_dof_index,
                cell->face_orientation(face_index_from_dof_index),
                cell->face_flip(face_index_from_dof_index),
                cell->face_rotation(face_index_from_dof_index)))
            dof_sign = -1.0;
        }

      const MappingKind mapping_kind = get_mapping_kind(dof_index);

      const unsigned int first =
        output_data.shape_function_to_row_table
          [dof_index * this->n_components() +
           this->get_nonzero_components(dof_index).first_selected_component()];

      // update the shape function values as necessary
      //
      // we only need to do this if the current cell is not a translation of
      // the previous one; or, even if it is a translation, if we use
      // mappings other than the standard mappings that require us to
      // recompute values and derivatives because of possible sign changes
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
                        fe_data.shape_values[dof_index][k][d];
                  break;
                }

              case mapping_covariant:
              case mapping_contravariant:
                {
                  mapping.transform(
                    make_array_view(fe_data.shape_values, dof_index),
                    mapping_kind,
                    mapping_internal,
                    make_array_view(fe_data.transformed_shape_values));

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_values(first + d, k) =
                        fe_data.transformed_shape_values[k][d];

                  break;
                }

              case mapping_raviart_thomas:
              case mapping_piola:
                {
                  mapping.transform(
                    make_array_view(fe_data.shape_values, dof_index),
                    mapping_piola,
                    mapping_internal,
                    make_array_view(fe_data.transformed_shape_values));
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_values(first + d, k) =
                        dof_sign * fe_data.transformed_shape_values[k][d];
                  break;
                }

              case mapping_nedelec:
                {
                  mapping.transform(
                    make_array_view(fe_data.shape_values, dof_index),
                    mapping_covariant,
                    mapping_internal,
                    make_array_view(fe_data.transformed_shape_values));

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_values(first + d, k) =
                        dof_sign * fe_data.transformed_shape_values[k][d];

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
                  mapping.transform(
                    make_array_view(fe_data.shape_grads, dof_index),
                    mapping_covariant,
                    mapping_internal,
                    make_array_view(fe_data.transformed_shape_grads));
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_gradients[first + d][k] =
                        fe_data.transformed_shape_grads[k][d];
                  break;
                }
              case mapping_covariant:
                {
                  mapping.transform(
                    make_array_view(fe_data.shape_grads, dof_index),
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
                        fe_data.transformed_shape_grads[k][d];

                  break;
                }
              case mapping_contravariant:
                {
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    fe_data.untransformed_shape_grads[k] =
                      fe_data.shape_grads[dof_index][k];
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
                      fe_data.shape_grads[dof_index][k];
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
                        dof_sign * fe_data.transformed_shape_grads[k][d];

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
                      fe_data.shape_grads[dof_index][k];

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
                        dof_sign * fe_data.transformed_shape_grads[k][d];

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
                    make_array_view(fe_data.shape_grad_grads, dof_index),
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
                      fe_data.shape_grad_grads[dof_index][k];

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
                      fe_data.shape_grad_grads[dof_index][k];

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
                      fe_data.shape_grad_grads[dof_index][k];

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
                        dof_sign * fe_data.transformed_shape_hessians[k][d];

                  break;
                }

              case mapping_nedelec:
                {
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    fe_data.untransformed_shape_hessian_tensors[k] =
                      fe_data.shape_grad_grads[dof_index][k];

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
                        dof_sign * fe_data.transformed_shape_hessians[k][d];

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
  const hp::QCollection<dim - 1> &                            quadrature,
  const Mapping<dim, spacedim> &                              mapping,
  const typename Mapping<dim, spacedim>::InternalDataBase &   mapping_internal,
  const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &                                                            mapping_data,
  const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    &output_data) const
{
  AssertDimension(quadrature.size(), 1);

  // convert data object to internal
  // data for this class. fails with
  // an exception if that is not
  // possible
  Assert(dynamic_cast<const InternalData *>(&fe_internal) != nullptr,
         ExcInternalError());
  const InternalData &fe_data = static_cast<const InternalData &>(fe_internal);

  const unsigned int n_q_points = quadrature[0].size();
  // offset determines which data set
  // to take (all data sets for all
  // faces are stored contiguously)

  const auto offset =
    QProjector<dim>::DataSetDescriptor::face(this->reference_cell(),
                                             face_no,
                                             cell->face_orientation(face_no),
                                             cell->face_flip(face_no),
                                             cell->face_rotation(face_no),
                                             n_q_points);

  // TODO: Size assertions

  // TODO: The dof_sign_change only affects Nedelec elements and is not the
  // correct thing on complicated meshes for higher order Nedelec elements.
  // Something similar to FE_Q should be done to permute dofs and to change the
  // dof signs. A static way using tables (as done in the RaviartThomas<dim>
  // class) is preferable.
  std::fill(fe_data.dof_sign_change.begin(),
            fe_data.dof_sign_change.end(),
            1.0);
  internal::FE_PolyTensor::get_dof_sign_change_nedelec(cell,
                                                       *this,
                                                       this->mapping_kind,
                                                       fe_data.dof_sign_change);

  // TODO: This, similarly to the Nedelec case, is just a legacy function in 2d
  // and affects only face_dofs of H(div) conformal FEs. It does nothing in 1d.
  // Also nothing in 3d since we take care of it by using the
  // adjust_quad_dof_sign_for_face_orientation_table.
  internal::FE_PolyTensor::get_dof_sign_change_h_div(cell,
                                                     *this,
                                                     this->mapping_kind,
                                                     fe_data.dof_sign_change);

  // What is the first dof_index on a quad?
  const unsigned int first_quad_index = this->get_first_quad_index();
  // How many dofs per quad and how many quad dofs do we have at all?
  const unsigned int n_dofs_per_quad = this->n_dofs_per_quad();
  const unsigned int n_quad_dofs =
    n_dofs_per_quad * GeometryInfo<dim>::faces_per_cell;

  for (unsigned int dof_index = 0; dof_index < this->n_dofs_per_cell();
       ++dof_index)
    {
      /*
       * This assumes that the dofs are ordered by first vertices, lines, quads
       * and volume dofs. Note that in 2d this always gives false.
       */
      const bool is_quad_dof =
        (dim == 2 ? false :
                    (first_quad_index <= dof_index) &&
                      (dof_index < first_quad_index + n_quad_dofs));

      // TODO: This hack is not pretty and it is only here to handle the 2d
      // case and the Nedelec legacy case. In 2d dof_sign of a face_dof is never
      // handled by the
      // >>if(is_quad_dof){...}<< but still a possible dof sign change must be
      // handled, also for line_dofs in 3d such as in Nedelec. In these cases
      // this is encoded in the array fe_data.dof_sign_change[dof_index]. In 3d
      // it is handles with a table. This array is allocated in
      // fe_poly_tensor.h.
      double dof_sign = 1.0;
      // under some circumstances fe_data.dof_sign_change is not allocated
      if (fe_data.update_each & update_values)
        dof_sign = fe_data.dof_sign_change[dof_index];

      if (is_quad_dof)
        {
          /*
           * Find the face belonging to this dof_index. This is integer
           * division.
           */
          unsigned int face_index_from_dof_index =
            (dof_index - first_quad_index) / (n_dofs_per_quad);

          unsigned int local_quad_dof_index = dof_index % n_dofs_per_quad;

          // Correct the dof_sign if necessary
          if (adjust_quad_dof_sign_for_face_orientation(
                local_quad_dof_index,
                face_index_from_dof_index,
                cell->face_orientation(face_index_from_dof_index),
                cell->face_flip(face_index_from_dof_index),
                cell->face_rotation(face_index_from_dof_index)))
            dof_sign = -1.0;
        }

      const MappingKind mapping_kind = get_mapping_kind(dof_index);

      const unsigned int first =
        output_data.shape_function_to_row_table
          [dof_index * this->n_components() +
           this->get_nonzero_components(dof_index).first_selected_component()];

      if (fe_data.update_each & update_values)
        {
          switch (mapping_kind)
            {
              case mapping_none:
                {
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_values(first + d, k) =
                        fe_data.shape_values[dof_index][k + offset][d];
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
                                                    dof_index,
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
                                                    dof_index,
                                                    offset,
                                                    n_q_points),
                                    mapping_piola,
                                    mapping_internal,
                                    transformed_shape_values);
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_values(first + d, k) =
                        dof_sign * transformed_shape_values[k][d];
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
                                                    dof_index,
                                                    offset,
                                                    n_q_points),
                                    mapping_covariant,
                                    mapping_internal,
                                    transformed_shape_values);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_values(first + d, k) =
                        dof_sign * transformed_shape_values[k][d];

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
                  mapping.transform(make_array_view(fe_data.shape_grads,
                                                    dof_index,
                                                    offset,
                                                    n_q_points),
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
                  mapping.transform(make_array_view(fe_data.shape_grads,
                                                    dof_index,
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
                      fe_data.shape_grads[dof_index][k + offset];
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
                      fe_data.shape_grads[dof_index][k + offset];
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
                        dof_sign * transformed_shape_grads[k][d];

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
                      fe_data.shape_grads[dof_index][k + offset];

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
                        dof_sign * transformed_shape_grads[k][d];

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
                                                    dof_index,
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
                      fe_data.shape_grad_grads[dof_index][k + offset];

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
                      fe_data.shape_grad_grads[dof_index][k + offset];

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
                      fe_data.shape_grad_grads[dof_index][k + offset];

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
                        dof_sign * transformed_shape_hessians[k][d];

                  break;
                }

              case mapping_nedelec:
                {
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    fe_data.untransformed_shape_hessian_tensors[k + offset] =
                      fe_data.shape_grad_grads[dof_index][k + offset];

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
                        dof_sign * transformed_shape_hessians[k][d];

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
  const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
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
  const auto offset =
    QProjector<dim>::DataSetDescriptor::subface(this->reference_cell(),
                                                face_no,
                                                sub_no,
                                                cell->face_orientation(face_no),
                                                cell->face_flip(face_no),
                                                cell->face_rotation(face_no),
                                                n_q_points,
                                                cell->subface_case(face_no));

  // TODO: Size assertions

  // TODO: The dof_sign_change only affects Nedelec elements and is not the
  // correct thing on complicated meshes for higher order Nedelec elements.
  // Something similar to FE_Q should be done to permute dofs and to change the
  // dof signs. A static way using tables (as done in the RaviartThomas<dim>
  // class) is preferable.
  std::fill(fe_data.dof_sign_change.begin(),
            fe_data.dof_sign_change.end(),
            1.0);
  internal::FE_PolyTensor::get_dof_sign_change_nedelec(cell,
                                                       *this,
                                                       this->mapping_kind,
                                                       fe_data.dof_sign_change);

  // TODO: This, similarly to the Nedelec case, is just a legacy function in 2d
  // and affects only face_dofs of H(div) conformal FEs. It does nothing in 1d.
  // Also nothing in 3d since we take care of it by using the
  // adjust_quad_dof_sign_for_face_orientation_table.
  internal::FE_PolyTensor::get_dof_sign_change_h_div(cell,
                                                     *this,
                                                     this->mapping_kind,
                                                     fe_data.dof_sign_change);

  // What is the first dof_index on a quad?
  const unsigned int first_quad_index = this->get_first_quad_index();
  // How many dofs per quad and how many quad dofs do we have at all?
  const unsigned int n_dofs_per_quad = this->n_dofs_per_quad();
  const unsigned int n_quad_dofs =
    n_dofs_per_quad * GeometryInfo<dim>::faces_per_cell;

  for (unsigned int dof_index = 0; dof_index < this->n_dofs_per_cell();
       ++dof_index)
    {
      /*
       * This assumes that the dofs are ordered by first vertices, lines, quads
       * and volume dofs. Note that in 2d this always gives false.
       */
      const bool is_quad_dof =
        (dim == 2 ? false :
                    (first_quad_index <= dof_index) &&
                      (dof_index < first_quad_index + n_quad_dofs));

      // TODO: This hack is not pretty and it is only here to handle the 2d
      // case and the Nedelec legacy case. In 2d dof_sign of a face_dof is never
      // handled by the
      // >>if(is_quad_dof){...}<< but still a possible dof sign change must be
      // handled, also for line_dofs in 3d such as in Nedelec. In these cases
      // this is encoded in the array fe_data.dof_sign_change[dof_index]. In 3d
      // it is handles with a table. This array is allocated in
      // fe_poly_tensor.h.
      double dof_sign = 1.0;
      // under some circumstances fe_data.dof_sign_change is not allocated
      if (fe_data.update_each & update_values)
        dof_sign = fe_data.dof_sign_change[dof_index];

      if (is_quad_dof)
        {
          /*
           * Find the face belonging to this dof_index. This is integer
           * division.
           */
          unsigned int face_index_from_dof_index =
            (dof_index - first_quad_index) / (n_dofs_per_quad);

          unsigned int local_quad_dof_index = dof_index % n_dofs_per_quad;

          // Correct the dof_sign if necessary
          if (adjust_quad_dof_sign_for_face_orientation(
                local_quad_dof_index,
                face_index_from_dof_index,
                cell->face_orientation(face_index_from_dof_index),
                cell->face_flip(face_index_from_dof_index),
                cell->face_rotation(face_index_from_dof_index)))
            dof_sign = -1.0;
        }

      const MappingKind mapping_kind = get_mapping_kind(dof_index);

      const unsigned int first =
        output_data.shape_function_to_row_table
          [dof_index * this->n_components() +
           this->get_nonzero_components(dof_index).first_selected_component()];

      if (fe_data.update_each & update_values)
        {
          switch (mapping_kind)
            {
              case mapping_none:
                {
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_values(first + d, k) =
                        fe_data.shape_values[dof_index][k + offset][d];
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
                                                    dof_index,
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
                                                    dof_index,
                                                    offset,
                                                    n_q_points),
                                    mapping_piola,
                                    mapping_internal,
                                    transformed_shape_values);
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_values(first + d, k) =
                        dof_sign * transformed_shape_values[k][d];
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
                                                    dof_index,
                                                    offset,
                                                    n_q_points),
                                    mapping_covariant,
                                    mapping_internal,
                                    transformed_shape_values);

                  for (unsigned int k = 0; k < n_q_points; ++k)
                    for (unsigned int d = 0; d < dim; ++d)
                      output_data.shape_values(first + d, k) =
                        dof_sign * transformed_shape_values[k][d];

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
                  mapping.transform(make_array_view(fe_data.shape_grads,
                                                    dof_index,
                                                    offset,
                                                    n_q_points),
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
                  mapping.transform(make_array_view(fe_data.shape_grads,
                                                    dof_index,
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
                        transformed_shape_grads[k][d];

                  break;
                }

              case mapping_contravariant:
                {
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    fe_data.untransformed_shape_grads[k + offset] =
                      fe_data.shape_grads[dof_index][k + offset];

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
                      fe_data.shape_grads[dof_index][k + offset];

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
                        dof_sign * transformed_shape_grads[k][d];

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
                      fe_data.shape_grads[dof_index][k + offset];

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
                        dof_sign * transformed_shape_grads[k][d];

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
                                                    dof_index,
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
                      fe_data.shape_grad_grads[dof_index][k + offset];

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
                      fe_data.shape_grad_grads[dof_index][k + offset];

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
                      fe_data.shape_grad_grads[dof_index][k + offset];

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
                        dof_sign * transformed_shape_hessians[k][d];

                  break;
                }

              case mapping_nedelec:
                {
                  for (unsigned int k = 0; k < n_q_points; ++k)
                    fe_data.untransformed_shape_hessian_tensors[k + offset] =
                      fe_data.shape_grad_grads[dof_index][k + offset];

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
                        dof_sign * transformed_shape_hessians[k][d];

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

  for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
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
