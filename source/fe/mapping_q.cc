// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/array_view.h>
#include <deal.II/base/derivative_form.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q_internal.h>

#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <boost/container/small_vector.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>


DEAL_II_NAMESPACE_OPEN


template <int dim, int spacedim>
MappingQ<dim, spacedim>::InternalData::InternalData(
  const unsigned int polynomial_degree)
  : polynomial_degree(polynomial_degree)
  , n_shape_functions(Utilities::fixed_power<dim>(polynomial_degree + 1))
  , line_support_points(QGaussLobatto<1>(polynomial_degree + 1))
  , tensor_product_quadrature(false)
  , output_data(nullptr)
{}



template <int dim, int spacedim>
std::size_t
MappingQ<dim, spacedim>::InternalData::memory_consumption() const
{
  return (
    Mapping<dim, spacedim>::InternalDataBase::memory_consumption() +
    MemoryConsumption::memory_consumption(quadrature_points) +
    MemoryConsumption::memory_consumption(unit_tangentials) +
    MemoryConsumption::memory_consumption(aux) +
    MemoryConsumption::memory_consumption(mapping_support_points) +
    MemoryConsumption::memory_consumption(cell_of_current_support_points) +
    MemoryConsumption::memory_consumption(volume_elements) +
    MemoryConsumption::memory_consumption(polynomial_degree) +
    MemoryConsumption::memory_consumption(n_shape_functions));
}



template <int dim, int spacedim>
void
MappingQ<dim, spacedim>::InternalData::reinit(const UpdateFlags update_flags,
                                              const Quadrature<dim> &quadrature)
{
  // store the flags in the internal data object so we can access them
  // in fill_fe_*_values()
  this->update_each = update_flags;

  const unsigned int n_q_points = quadrature.size();

  if (this->update_each & update_volume_elements)
    volume_elements.resize(n_q_points);

  tensor_product_quadrature = quadrature.is_tensor_product();

  // use of MatrixFree only for higher order elements and with more than one
  // point where tensor products do not make sense
  if (polynomial_degree < 2 || n_q_points == 1)
    tensor_product_quadrature = false;

  if (dim > 1)
    {
      // find out if the one-dimensional formula is the same
      // in all directions
      if (tensor_product_quadrature)
        {
          const std::array<Quadrature<1>, dim> &quad_array =
            quadrature.get_tensor_basis();
          for (unsigned int i = 1; i < dim && tensor_product_quadrature; ++i)
            {
              if (quad_array[i - 1].size() != quad_array[i].size())
                {
                  tensor_product_quadrature = false;
                  break;
                }
              else
                {
                  const std::vector<Point<1>> &points_1 =
                    quad_array[i - 1].get_points();
                  const std::vector<Point<1>> &points_2 =
                    quad_array[i].get_points();
                  const std::vector<double> &weights_1 =
                    quad_array[i - 1].get_weights();
                  const std::vector<double> &weights_2 =
                    quad_array[i].get_weights();
                  for (unsigned int j = 0; j < quad_array[i].size(); ++j)
                    {
                      if (std::abs(points_1[j][0] - points_2[j][0]) > 1.e-10 ||
                          std::abs(weights_1[j] - weights_2[j]) > 1.e-10)
                        {
                          tensor_product_quadrature = false;
                          break;
                        }
                    }
                }
            }

          if (tensor_product_quadrature)
            {
              // use a 1d FE_DGQ and adjust the hierarchic -> lexicographic
              // numbering manually (building an FE_Q<dim> is relatively
              // expensive due to constraints)
              const FE_DGQ<1> fe(polynomial_degree);
              shape_info.reinit(quadrature.get_tensor_basis()[0], fe);
              shape_info.lexicographic_numbering =
                FETools::lexicographic_to_hierarchic_numbering<dim>(
                  polynomial_degree);
              shape_info.n_q_points = n_q_points;
              shape_info.dofs_per_component_on_cell =
                Utilities::pow(polynomial_degree + 1, dim);
            }
        }
    }
}



template <int dim, int spacedim>
void
MappingQ<dim, spacedim>::InternalData::initialize_face(
  const UpdateFlags      update_flags,
  const Quadrature<dim> &quadrature,
  const unsigned int     n_original_q_points)
{
  reinit(update_flags, quadrature);

  quadrature_points = quadrature.get_points();

  if (dim > 1 && tensor_product_quadrature)
    {
      constexpr unsigned int facedim = dim - 1;
      const FE_DGQ<1>        fe(polynomial_degree);
      shape_info.reinit(quadrature.get_tensor_basis()[0], fe);
      shape_info.lexicographic_numbering =
        FETools::lexicographic_to_hierarchic_numbering<facedim>(
          polynomial_degree);
      shape_info.n_q_points = n_original_q_points;
      shape_info.dofs_per_component_on_cell =
        Utilities::pow(polynomial_degree + 1, dim);
    }

  if (dim > 1)
    {
      if (this->update_each &
          (update_boundary_forms | update_normal_vectors | update_JxW_values))
        {
          aux.resize(dim - 1);
          aux[0].resize(n_original_q_points);
          if (dim > 2)
            aux[1].resize(n_original_q_points);

          // Compute tangentials to the unit cell.
          for (const unsigned int i : GeometryInfo<dim>::face_indices())
            {
              unit_tangentials[i].resize(n_original_q_points);
              std::fill(unit_tangentials[i].begin(),
                        unit_tangentials[i].end(),
                        GeometryInfo<dim>::unit_tangential_vectors[i][0]);
              if (dim > 2)
                {
                  unit_tangentials[GeometryInfo<dim>::faces_per_cell + i]
                    .resize(n_original_q_points);
                  std::fill(
                    unit_tangentials[GeometryInfo<dim>::faces_per_cell + i]
                      .begin(),
                    unit_tangentials[GeometryInfo<dim>::faces_per_cell + i]
                      .end(),
                    GeometryInfo<dim>::unit_tangential_vectors[i][1]);
                }
            }
        }
    }
}



template <int dim, int spacedim>
MappingQ<dim, spacedim>::MappingQ(const unsigned int p)
  : polynomial_degree(p)
  , line_support_points(
      QGaussLobatto<1>(this->polynomial_degree + 1).get_points())
  , polynomials_1d(
      Polynomials::generate_complete_Lagrange_basis(line_support_points))
  , renumber_lexicographic_to_hierarchic(
      FETools::lexicographic_to_hierarchic_numbering<dim>(p))
  , unit_cell_support_points(
      internal::MappingQImplementation::unit_support_points<dim>(
        line_support_points,
        renumber_lexicographic_to_hierarchic))
  , support_point_weights_perimeter_to_interior(
      internal::MappingQImplementation::
        compute_support_point_weights_perimeter_to_interior(
          this->polynomial_degree,
          dim))
  , support_point_weights_cell(
      internal::MappingQImplementation::compute_support_point_weights_cell<dim>(
        this->polynomial_degree))
{
  Assert(p >= 1,
         ExcMessage("It only makes sense to create polynomial mappings "
                    "with a polynomial degree greater or equal to one."));
}



template <int dim, int spacedim>
MappingQ<dim, spacedim>::MappingQ(const MappingQ<dim, spacedim> &mapping)
  : polynomial_degree(mapping.polynomial_degree)
  , line_support_points(mapping.line_support_points)
  , polynomials_1d(mapping.polynomials_1d)
  , renumber_lexicographic_to_hierarchic(
      mapping.renumber_lexicographic_to_hierarchic)
  , unit_cell_support_points(mapping.unit_cell_support_points)
  , support_point_weights_perimeter_to_interior(
      mapping.support_point_weights_perimeter_to_interior)
  , support_point_weights_cell(mapping.support_point_weights_cell)
{}



template <int dim, int spacedim>
std::unique_ptr<Mapping<dim, spacedim>>
MappingQ<dim, spacedim>::clone() const
{
  return std::make_unique<MappingQ<dim, spacedim>>(*this);
}



template <int dim, int spacedim>
unsigned int
MappingQ<dim, spacedim>::get_degree() const
{
  return polynomial_degree;
}



template <int dim, int spacedim>
Point<spacedim>
MappingQ<dim, spacedim>::transform_unit_to_real_cell(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const Point<dim>                                           &p) const
{
  if (polynomial_degree == 1)
    {
      const auto vertices = this->get_vertices(cell);
      return Point<spacedim>(
        internal::evaluate_tensor_product_value_linear(vertices.data(), p));
    }
  else
    return Point<spacedim>(internal::evaluate_tensor_product_value(
      polynomials_1d,
      make_const_array_view(this->compute_mapping_support_points(cell)),
      p,
      polynomials_1d.size() == 2,
      renumber_lexicographic_to_hierarchic));
}


// In the code below, GCC tries to instantiate MappingQ<3,4> when
// seeing which of the overloaded versions of
// do_transform_real_to_unit_cell_internal() to call. This leads to bad
// error messages and, generally, nothing very good. Avoid this by ensuring
// that this class exists, but does not have an inner InternalData
// type, thereby ruling out the codim-1 version of the function
// below when doing overload resolution.
template <>
class MappingQ<3, 4>
{};



// visual studio freaks out when trying to determine if
// do_transform_real_to_unit_cell_internal with dim=3 and spacedim=4 is a good
// candidate. So instead of letting the compiler pick the correct overload, we
// use template specialization to make sure we pick up the right function to
// call:

template <int dim, int spacedim>
Point<dim>
MappingQ<dim, spacedim>::transform_real_to_unit_cell_internal(
  const typename Triangulation<dim, spacedim>::cell_iterator &,
  const Point<spacedim> &,
  const Point<dim> &) const
{
  // default implementation (should never be called)
  DEAL_II_ASSERT_UNREACHABLE();
  return {};
}



template <>
Point<1>
MappingQ<1, 1>::transform_real_to_unit_cell_internal(
  const Triangulation<1, 1>::cell_iterator &cell,
  const Point<1>                           &p,
  const Point<1>                           &initial_p_unit) const
{
  // dispatch to the various specializations for spacedim=dim,
  // spacedim=dim+1, etc
  if (polynomial_degree == 1)
    {
      const auto vertices = this->get_vertices(cell);
      return internal::MappingQImplementation::
        do_transform_real_to_unit_cell_internal<1>(
          p,
          initial_p_unit,
          ArrayView<const Point<1>>(vertices.data(), vertices.size()),
          polynomials_1d,
          renumber_lexicographic_to_hierarchic);
    }
  else
    return internal::MappingQImplementation::
      do_transform_real_to_unit_cell_internal<1>(
        p,
        initial_p_unit,
        make_const_array_view(this->compute_mapping_support_points(cell)),
        polynomials_1d,
        renumber_lexicographic_to_hierarchic);
}



template <>
Point<2>
MappingQ<2, 2>::transform_real_to_unit_cell_internal(
  const Triangulation<2, 2>::cell_iterator &cell,
  const Point<2>                           &p,
  const Point<2>                           &initial_p_unit) const
{
  if (polynomial_degree == 1)
    {
      const auto vertices = this->get_vertices(cell);
      return internal::MappingQImplementation::
        do_transform_real_to_unit_cell_internal<2>(
          p,
          initial_p_unit,
          ArrayView<const Point<2>>(vertices.data(), vertices.size()),
          polynomials_1d,
          renumber_lexicographic_to_hierarchic);
    }
  else
    return internal::MappingQImplementation::
      do_transform_real_to_unit_cell_internal<2>(
        p,
        initial_p_unit,
        make_const_array_view(this->compute_mapping_support_points(cell)),
        polynomials_1d,
        renumber_lexicographic_to_hierarchic);
}



template <>
Point<3>
MappingQ<3, 3>::transform_real_to_unit_cell_internal(
  const Triangulation<3, 3>::cell_iterator &cell,
  const Point<3>                           &p,
  const Point<3>                           &initial_p_unit) const
{
  if (polynomial_degree == 1)
    {
      const auto vertices = this->get_vertices(cell);
      return internal::MappingQImplementation::
        do_transform_real_to_unit_cell_internal<3>(
          p,
          initial_p_unit,
          ArrayView<const Point<3>>(vertices.data(), vertices.size()),
          polynomials_1d,
          renumber_lexicographic_to_hierarchic);
    }
  else
    return internal::MappingQImplementation::
      do_transform_real_to_unit_cell_internal<3>(
        p,
        initial_p_unit,
        make_const_array_view(this->compute_mapping_support_points(cell)),
        polynomials_1d,
        renumber_lexicographic_to_hierarchic);
}



template <>
Point<1>
MappingQ<1, 2>::transform_real_to_unit_cell_internal(
  const Triangulation<1, 2>::cell_iterator &cell,
  const Point<2>                           &p,
  const Point<1>                           &initial_p_unit) const
{
  const int dim      = 1;
  const int spacedim = 2;

  const Quadrature<dim> point_quadrature(initial_p_unit);

  UpdateFlags update_flags = update_quadrature_points | update_jacobians;
  if (spacedim > dim)
    update_flags |= update_jacobian_grads;
  auto mdata = Utilities::dynamic_unique_cast<InternalData>(
    get_data(update_flags, point_quadrature));

  mdata->mapping_support_points = this->compute_mapping_support_points(cell);

  // dispatch to the various specializations for spacedim=dim,
  // spacedim=dim+1, etc
  return internal::MappingQImplementation::
    do_transform_real_to_unit_cell_internal_codim1<1>(
      p,
      initial_p_unit,
      make_const_array_view(mdata->mapping_support_points),
      polynomials_1d,
      renumber_lexicographic_to_hierarchic);
}



template <>
Point<2>
MappingQ<2, 3>::transform_real_to_unit_cell_internal(
  const Triangulation<2, 3>::cell_iterator &cell,
  const Point<3>                           &p,
  const Point<2>                           &initial_p_unit) const
{
  const int dim      = 2;
  const int spacedim = 3;

  const Quadrature<dim> point_quadrature(initial_p_unit);

  UpdateFlags update_flags = update_quadrature_points | update_jacobians;
  if (spacedim > dim)
    update_flags |= update_jacobian_grads;
  auto mdata = Utilities::dynamic_unique_cast<InternalData>(
    get_data(update_flags, point_quadrature));

  mdata->mapping_support_points = this->compute_mapping_support_points(cell);

  // dispatch to the various specializations for spacedim=dim,
  // spacedim=dim+1, etc
  return internal::MappingQImplementation::
    do_transform_real_to_unit_cell_internal_codim1<2>(
      p,
      initial_p_unit,
      make_const_array_view(mdata->mapping_support_points),
      polynomials_1d,
      renumber_lexicographic_to_hierarchic);
}



template <>
Point<1>
MappingQ<1, 3>::transform_real_to_unit_cell_internal(
  const Triangulation<1, 3>::cell_iterator &,
  const Point<3> &,
  const Point<1> &) const
{
  DEAL_II_NOT_IMPLEMENTED();
  return {};
}



template <int dim, int spacedim>
Point<dim>
MappingQ<dim, spacedim>::transform_real_to_unit_cell(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const Point<spacedim>                                      &p) const
{
  // Use an exact formula if one is available. this is only the case
  // for Q1 mappings in 1d, and in 2d if dim==spacedim
  if ((polynomial_degree == 1) &&
      ((dim == 1) || ((dim == 2) && (dim == spacedim))))
    {
      // The dimension-dependent algorithms are much faster (about 25-45x in
      // 2d) but fail most of the time when the given point (p) is not in the
      // cell. The dimension-independent Newton algorithm given below is
      // slower, but more robust (though it still sometimes fails). Therefore
      // this function implements the following strategy based on the
      // p's dimension:
      //
      // * In 1d this mapping is linear, so the mapping is always invertible
      //   (and the exact formula is known) as long as the cell has non-zero
      //   length.
      // * In 2d the exact (quadratic) formula is called first. If either the
      //   exact formula does not succeed (negative discriminant in the
      //   quadratic formula) or succeeds but finds a solution outside of the
      //   unit cell, then the Newton solver is called. The rationale for the
      //   second choice is that the exact formula may provide two different
      //   answers when mapping a point outside of the real cell, but the
      //   Newton solver (if it converges) will only return one answer.
      //   Otherwise the exact formula successfully found a point in the unit
      //   cell and that value is returned.
      // * In 3d there is no (known to the authors) exact formula, so the Newton
      //   algorithm is used.
      const auto vertices_ = this->get_vertices(cell);

      std::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell>
        vertices;
      for (unsigned int i = 0; i < vertices.size(); ++i)
        vertices[i] = vertices_[i];

      try
        {
          switch (dim)
            {
              case 1:
                {
                  // formula not subject to any issues in 1d
                  if (spacedim == 1)
                    return internal::MappingQ1::transform_real_to_unit_cell(
                      vertices, p);
                  else
                    break;
                }

              case 2:
                {
                  const Point<dim> point =
                    internal::MappingQ1::transform_real_to_unit_cell(vertices,
                                                                     p);

                  // formula not guaranteed to work for points outside of
                  // the cell. only take the computed point if it lies
                  // inside the reference cell
                  const double eps = 1e-15;
                  if (-eps <= point[1] && point[1] <= 1 + eps &&
                      -eps <= point[0] && point[0] <= 1 + eps)
                    {
                      return point;
                    }
                  else
                    break;
                }

              default:
                {
                  // we should not get here, based on the if-condition at the
                  // top
                  DEAL_II_ASSERT_UNREACHABLE();
                }
            }
        }
      catch (
        const typename Mapping<spacedim, spacedim>::ExcTransformationFailed &)
        {
          // simply fall through and continue on to the standard Newton code
        }
    }
  else
    {
      // we can't use an explicit formula,
    }


  // Find the initial value for the Newton iteration by a normal
  // projection to the least square plane determined by the vertices
  // of the cell
  Point<dim> initial_p_unit;
  if (this->preserves_vertex_locations())
    {
      initial_p_unit = cell->real_to_unit_cell_affine_approximation(p);
      // in 1d with spacedim > 1 the affine approximation is exact
      if (dim == 1 && polynomial_degree == 1)
        return initial_p_unit;
    }
  else
    {
      // else, we simply use the mid point
      for (unsigned int d = 0; d < dim; ++d)
        initial_p_unit[d] = 0.5;
    }

  // perform the Newton iteration and return the result. note that this
  // statement may throw an exception, which we simply pass up to the caller
  const Point<dim> p_unit =
    this->transform_real_to_unit_cell_internal(cell, p, initial_p_unit);
  AssertThrow(p_unit[0] != std::numeric_limits<double>::lowest(),
              (typename Mapping<dim, spacedim>::ExcTransformationFailed()));
  return p_unit;
}



template <int dim, int spacedim>
void
MappingQ<dim, spacedim>::transform_points_real_to_unit_cell(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const ArrayView<const Point<spacedim>>                     &real_points,
  const ArrayView<Point<dim>>                                &unit_points) const
{
  // Go to base class functions for dim < spacedim because it is not yet
  // implemented with optimized code.
  if (dim < spacedim)
    {
      Mapping<dim, spacedim>::transform_points_real_to_unit_cell(cell,
                                                                 real_points,
                                                                 unit_points);
      return;
    }

  AssertDimension(real_points.size(), unit_points.size());
  std::vector<Point<spacedim>> support_points_higher_order;
  boost::container::small_vector<Point<spacedim>,
#ifndef _MSC_VER
                                 ReferenceCells::max_n_vertices<dim>()
#else
                                 GeometryInfo<dim>::vertices_per_cell
#endif
                                 >
    vertices;
  if (polynomial_degree == 1)
    vertices = this->get_vertices(cell);
  else
    support_points_higher_order = this->compute_mapping_support_points(cell);
  const ArrayView<const Point<spacedim>> support_points(
    polynomial_degree == 1 ? vertices.data() :
                             support_points_higher_order.data(),
    Utilities::pow(polynomial_degree + 1, dim));

  // From the given (high-order) support points, now only pick the first
  // 2^dim points and construct an affine approximation from those.
  internal::MappingQImplementation::InverseQuadraticApproximation<dim, spacedim>
    inverse_approximation(support_points, unit_cell_support_points);

  const unsigned int n_points = real_points.size();
  const unsigned int n_lanes  = VectorizedArray<double>::size();

  // Use the more heavy VectorizedArray code path if there is more than
  // one point left to compute
  for (unsigned int i = 0; i < n_points; i += n_lanes)
    if (n_points - i > 1)
      {
        Point<spacedim, VectorizedArray<double>> p_vec;
        for (unsigned int j = 0; j < n_lanes; ++j)
          if (i + j < n_points)
            for (unsigned int d = 0; d < spacedim; ++d)
              p_vec[d][j] = real_points[i + j][d];
          else
            for (unsigned int d = 0; d < spacedim; ++d)
              p_vec[d][j] = real_points[i][d];

        Point<dim, VectorizedArray<double>> unit_point =
          internal::MappingQImplementation::
            do_transform_real_to_unit_cell_internal<dim, spacedim>(
              p_vec,
              inverse_approximation.compute(p_vec),
              support_points,
              polynomials_1d,
              renumber_lexicographic_to_hierarchic);

        // If the vectorized computation failed, it could be that only some of
        // the lanes failed but others would have succeeded if we had let them
        // compute alone without interference (like negative Jacobian
        // determinants) from other SIMD lanes. Repeat the computation in this
        // unlikely case with scalar arguments.
        for (unsigned int j = 0; j < n_lanes && i + j < n_points; ++j)
          if (unit_point[0][j] != std::numeric_limits<double>::lowest())
            for (unsigned int d = 0; d < dim; ++d)
              unit_points[i + j][d] = unit_point[d][j];
          else
            unit_points[i + j] = internal::MappingQImplementation::
              do_transform_real_to_unit_cell_internal<dim, spacedim>(
                real_points[i + j],
                inverse_approximation.compute(real_points[i + j]),
                support_points,
                polynomials_1d,
                renumber_lexicographic_to_hierarchic);
      }
    else
      unit_points[i] = internal::MappingQImplementation::
        do_transform_real_to_unit_cell_internal<dim, spacedim>(
          real_points[i],
          inverse_approximation.compute(real_points[i]),
          support_points,
          polynomials_1d,
          renumber_lexicographic_to_hierarchic);
}



template <int dim, int spacedim>
UpdateFlags
MappingQ<dim, spacedim>::requires_update_flags(const UpdateFlags in) const
{
  // add flags if the respective quantities are necessary to compute
  // what we need. note that some flags appear in both the conditions
  // and in subsequent set operations. this leads to some circular
  // logic. the only way to treat this is to iterate. since there are
  // 5 if-clauses in the loop, it will take at most 5 iterations to
  // converge. do them:
  UpdateFlags out = in;
  for (unsigned int i = 0; i < 5; ++i)
    {
      // The following is a little incorrect:
      // If not applied on a face,
      // update_boundary_forms does not
      // make sense. On the other hand,
      // it is necessary on a
      // face. Currently,
      // update_boundary_forms is simply
      // ignored for the interior of a
      // cell.
      if (out & (update_JxW_values | update_normal_vectors))
        out |= update_boundary_forms;

      if (out &
          (update_covariant_transformation | update_jacobian_grads |
           update_jacobians | update_boundary_forms | update_normal_vectors))
        out |= update_contravariant_transformation;

      if (out &
          (update_inverse_jacobians | update_jacobian_pushed_forward_grads |
           update_jacobian_pushed_forward_2nd_derivatives |
           update_jacobian_pushed_forward_3rd_derivatives))
        out |= update_covariant_transformation;

      // The contravariant transformation is used in the Piola
      // transformation, which requires the determinant of the Jacobi
      // matrix of the transformation.  Because we have no way of
      // knowing here whether the finite element wants to use the
      // contravariant or the Piola transforms, we add the JxW values
      // to the list of flags to be updated for each cell.
      if (out & update_contravariant_transformation)
        out |= update_volume_elements;

      // the same is true when computing normal vectors: they require
      // the determinant of the Jacobian
      if (out & update_normal_vectors)
        out |= update_volume_elements;
    }

  return out;
}



template <int dim, int spacedim>
std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
MappingQ<dim, spacedim>::get_data(const UpdateFlags      update_flags,
                                  const Quadrature<dim> &q) const
{
  std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase> data_ptr =
    std::make_unique<InternalData>(polynomial_degree);
  data_ptr->reinit(update_flags, q);
  return data_ptr;
}



template <int dim, int spacedim>
std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
MappingQ<dim, spacedim>::get_face_data(
  const UpdateFlags               update_flags,
  const hp::QCollection<dim - 1> &quadrature) const
{
  AssertDimension(quadrature.size(), 1);

  std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase> data_ptr =
    std::make_unique<InternalData>(polynomial_degree);
  auto &data = dynamic_cast<InternalData &>(*data_ptr);
  data.initialize_face(this->requires_update_flags(update_flags),
                       QProjector<dim>::project_to_all_faces(
                         ReferenceCells::get_hypercube<dim>(), quadrature[0]),
                       quadrature[0].size());

  return data_ptr;
}



template <int dim, int spacedim>
std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
MappingQ<dim, spacedim>::get_subface_data(
  const UpdateFlags          update_flags,
  const Quadrature<dim - 1> &quadrature) const
{
  std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase> data_ptr =
    std::make_unique<InternalData>(polynomial_degree);
  auto &data = dynamic_cast<InternalData &>(*data_ptr);
  data.initialize_face(this->requires_update_flags(update_flags),
                       QProjector<dim>::project_to_all_subfaces(
                         ReferenceCells::get_hypercube<dim>(), quadrature),
                       quadrature.size());

  return data_ptr;
}



template <int dim, int spacedim>
CellSimilarity::Similarity
MappingQ<dim, spacedim>::fill_fe_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const CellSimilarity::Similarity                            cell_similarity,
  const Quadrature<dim>                                      &quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  // ensure that the following static_cast is really correct:
  Assert(dynamic_cast<const InternalData *>(&internal_data) != nullptr,
         ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(internal_data);
  data.output_data         = &output_data;

  const unsigned int n_q_points = quadrature.size();

  // recompute the support points of the transformation of this
  // cell. we tried to be clever here in an earlier version of the
  // library by checking whether the cell is the same as the one we
  // had visited last, but it turns out to be difficult to determine
  // that because a cell for the purposes of a mapping is
  // characterized not just by its (triangulation, level, index)
  // triple, but also by the locations of its vertices, the manifold
  // object attached to the cell and all of its bounding faces/edges,
  // etc. to reliably test that the "cell" we are on is, therefore,
  // not easily done
  if (polynomial_degree == 1)
    {
      data.mapping_support_points.resize(GeometryInfo<dim>::vertices_per_cell);
      const auto vertices = this->get_vertices(cell);
      for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
        data.mapping_support_points[i] = vertices[i];
    }
  else
    data.mapping_support_points = this->compute_mapping_support_points(cell);

  data.cell_of_current_support_points = cell;

  // if the order of the mapping is greater than 1, then do not reuse any cell
  // similarity information. This is necessary because the cell similarity
  // value is computed with just cell vertices and does not take into account
  // cell curvature.
  const CellSimilarity::Similarity computed_cell_similarity =
    (polynomial_degree == 1 && this->preserves_vertex_locations() ?
       cell_similarity :
       CellSimilarity::none);

  if (dim > 1 && data.tensor_product_quadrature)
    {
      internal::MappingQImplementation::
        maybe_update_q_points_Jacobians_and_grads_tensor<dim, spacedim>(
          computed_cell_similarity,
          data,
          output_data.quadrature_points,
          output_data.jacobians,
          output_data.inverse_jacobians,
          output_data.jacobian_grads);
    }
  else
    {
      internal::MappingQImplementation::maybe_update_q_points_Jacobians_generic(
        computed_cell_similarity,
        data,
        make_array_view(quadrature.get_points()),
        polynomials_1d,
        renumber_lexicographic_to_hierarchic,
        output_data.quadrature_points,
        output_data.jacobians,
        output_data.inverse_jacobians);

      internal::MappingQImplementation::maybe_update_jacobian_grads<dim,
                                                                    spacedim>(
        computed_cell_similarity,
        data,
        make_array_view(quadrature.get_points()),
        polynomials_1d,
        renumber_lexicographic_to_hierarchic,
        output_data.jacobian_grads);
    }

  internal::MappingQImplementation::maybe_update_jacobian_pushed_forward_grads<
    dim,
    spacedim>(computed_cell_similarity,
              data,
              make_array_view(quadrature.get_points()),
              polynomials_1d,
              renumber_lexicographic_to_hierarchic,
              output_data.jacobian_pushed_forward_grads);

  internal::MappingQImplementation::maybe_update_jacobian_2nd_derivatives<
    dim,
    spacedim>(computed_cell_similarity,
              data,
              make_array_view(quadrature.get_points()),
              polynomials_1d,
              renumber_lexicographic_to_hierarchic,
              output_data.jacobian_2nd_derivatives);

  internal::MappingQImplementation::
    maybe_update_jacobian_pushed_forward_2nd_derivatives<dim, spacedim>(
      computed_cell_similarity,
      data,
      make_array_view(quadrature.get_points()),
      polynomials_1d,
      renumber_lexicographic_to_hierarchic,
      output_data.jacobian_pushed_forward_2nd_derivatives);

  internal::MappingQImplementation::maybe_update_jacobian_3rd_derivatives<
    dim,
    spacedim>(computed_cell_similarity,
              data,
              make_array_view(quadrature.get_points()),
              polynomials_1d,
              renumber_lexicographic_to_hierarchic,
              output_data.jacobian_3rd_derivatives);

  internal::MappingQImplementation::
    maybe_update_jacobian_pushed_forward_3rd_derivatives<dim, spacedim>(
      computed_cell_similarity,
      data,
      make_array_view(quadrature.get_points()),
      polynomials_1d,
      renumber_lexicographic_to_hierarchic,
      output_data.jacobian_pushed_forward_3rd_derivatives);

  const UpdateFlags          update_flags = data.update_each;
  const std::vector<double> &weights      = quadrature.get_weights();

  // Multiply quadrature weights by absolute value of Jacobian determinants or
  // the area element g=sqrt(DX^t DX) in case of codim > 0

  if (update_flags & (update_normal_vectors | update_JxW_values))
    {
      Assert(!(update_flags & update_JxW_values) ||
               (output_data.JxW_values.size() == n_q_points),
             ExcDimensionMismatch(output_data.JxW_values.size(), n_q_points));

      Assert(!(update_flags & update_normal_vectors) ||
               (output_data.normal_vectors.size() == n_q_points),
             ExcDimensionMismatch(output_data.normal_vectors.size(),
                                  n_q_points));


      if (computed_cell_similarity != CellSimilarity::translation)
        for (unsigned int point = 0; point < n_q_points; ++point)
          {
            if (dim == spacedim)
              {
                const double det = data.volume_elements[point];

                // check for distorted cells.

                // TODO: this allows for anisotropies of up to 1e6 in 3d and
                // 1e12 in 2d. might want to find a finer
                // (dimension-independent) criterion
                Assert(det >
                         1e-12 * Utilities::fixed_power<dim>(
                                   cell->diameter() / std::sqrt(double(dim))),
                       (typename Mapping<dim, spacedim>::ExcDistortedMappedCell(
                         cell->center(), det, point)));

                output_data.JxW_values[point] = weights[point] * det;
              }
            // if dim==spacedim, then there is no cell normal to
            // compute. since this is for FEValues (and not FEFaceValues),
            // there are also no face normals to compute
            else // codim>0 case
              {
                Tensor<1, spacedim> DX_t[dim];
                for (unsigned int i = 0; i < spacedim; ++i)
                  for (unsigned int j = 0; j < dim; ++j)
                    DX_t[j][i] = output_data.jacobians[point][i][j];

                Tensor<2, dim> G; // First fundamental form
                for (unsigned int i = 0; i < dim; ++i)
                  for (unsigned int j = 0; j < dim; ++j)
                    G[i][j] = DX_t[i] * DX_t[j];

                if (update_flags & update_JxW_values)
                  output_data.JxW_values[point] =
                    std::sqrt(determinant(G)) * weights[point];

                if (computed_cell_similarity ==
                    CellSimilarity::inverted_translation)
                  {
                    // we only need to flip the normal
                    if (update_flags & update_normal_vectors)
                      output_data.normal_vectors[point] *= -1.;
                  }
                else
                  {
                    if (update_flags & update_normal_vectors)
                      {
                        Assert(spacedim == dim + 1,
                               ExcMessage(
                                 "There is no (unique) cell normal for " +
                                 Utilities::int_to_string(dim) +
                                 "-dimensional cells in " +
                                 Utilities::int_to_string(spacedim) +
                                 "-dimensional space. This only works if the "
                                 "space dimension is one greater than the "
                                 "dimensionality of the mesh cells."));

                        if (dim == 1)
                          output_data.normal_vectors[point] =
                            cross_product_2d(-DX_t[0]);
                        else // dim == 2
                          output_data.normal_vectors[point] =
                            cross_product_3d(DX_t[0], DX_t[1]);

                        output_data.normal_vectors[point] /=
                          output_data.normal_vectors[point].norm();

                        if (cell->direction_flag() == false)
                          output_data.normal_vectors[point] *= -1.;
                      }
                  }
              } // codim>0 case
          }
    }

  return computed_cell_similarity;
}



template <int dim, int spacedim>
void
MappingQ<dim, spacedim>::fill_fe_face_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const hp::QCollection<dim - 1>                             &quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  AssertDimension(quadrature.size(), 1);

  // ensure that the following cast is really correct:
  Assert((dynamic_cast<const InternalData *>(&internal_data) != nullptr),
         ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(internal_data);
  data.output_data         = &output_data;

  // if necessary, recompute the support points of the transformation of this
  // cell (note that we need to first check the triangulation pointer, since
  // otherwise the second test might trigger an exception if the
  // triangulations are not the same)
  if ((data.mapping_support_points.empty()) ||
      (&cell->get_triangulation() !=
       &data.cell_of_current_support_points->get_triangulation()) ||
      (cell != data.cell_of_current_support_points))
    {
      if (polynomial_degree == 1)
        {
          data.mapping_support_points.resize(
            GeometryInfo<dim>::vertices_per_cell);
          const auto vertices = this->get_vertices(cell);
          for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell;
               ++i)
            data.mapping_support_points[i] = vertices[i];
        }
      else
        data.mapping_support_points =
          this->compute_mapping_support_points(cell);
      data.cell_of_current_support_points = cell;
    }

  internal::MappingQImplementation::do_fill_fe_face_values(
    *this,
    cell,
    face_no,
    numbers::invalid_unsigned_int,
    QProjector<dim>::DataSetDescriptor::face(
      ReferenceCells::get_hypercube<dim>(),
      face_no,
      cell->combined_face_orientation(face_no),
      quadrature[0].size()),
    quadrature[0],
    data,
    polynomials_1d,
    renumber_lexicographic_to_hierarchic,
    output_data);
}



template <int dim, int spacedim>
void
MappingQ<dim, spacedim>::fill_fe_subface_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const unsigned int                                          subface_no,
  const Quadrature<dim - 1>                                  &quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  // ensure that the following cast is really correct:
  Assert((dynamic_cast<const InternalData *>(&internal_data) != nullptr),
         ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(internal_data);
  data.output_data         = &output_data;

  // if necessary, recompute the support points of the transformation of this
  // cell (note that we need to first check the triangulation pointer, since
  // otherwise the second test might trigger an exception if the
  // triangulations are not the same)
  if ((data.mapping_support_points.empty()) ||
      (&cell->get_triangulation() !=
       &data.cell_of_current_support_points->get_triangulation()) ||
      (cell != data.cell_of_current_support_points))
    {
      if (polynomial_degree == 1)
        {
          data.mapping_support_points.resize(
            GeometryInfo<dim>::vertices_per_cell);
          const auto vertices = this->get_vertices(cell);
          for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell;
               ++i)
            data.mapping_support_points[i] = vertices[i];
        }
      else
        data.mapping_support_points =
          this->compute_mapping_support_points(cell);
      data.cell_of_current_support_points = cell;
    }

  internal::MappingQImplementation::do_fill_fe_face_values(
    *this,
    cell,
    face_no,
    subface_no,
    QProjector<dim>::DataSetDescriptor::subface(
      ReferenceCells::get_hypercube<dim>(),
      face_no,
      subface_no,
      cell->combined_face_orientation(face_no),
      quadrature.size(),
      cell->subface_case(face_no)),
    quadrature,
    data,
    polynomials_1d,
    renumber_lexicographic_to_hierarchic,
    output_data);
}



template <int dim, int spacedim>
void
MappingQ<dim, spacedim>::fill_fe_immersed_surface_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const NonMatching::ImmersedSurfaceQuadrature<dim>          &quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  Assert(dim == spacedim, ExcNotImplemented());

  // ensure that the following static_cast is really correct:
  Assert(dynamic_cast<const InternalData *>(&internal_data) != nullptr,
         ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(internal_data);
  data.output_data         = &output_data;

  const unsigned int n_q_points = quadrature.size();

  if (polynomial_degree == 1)
    {
      data.mapping_support_points.resize(GeometryInfo<dim>::vertices_per_cell);
      const auto vertices = this->get_vertices(cell);
      for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
        data.mapping_support_points[i] = vertices[i];
    }
  else
    data.mapping_support_points = this->compute_mapping_support_points(cell);
  data.cell_of_current_support_points = cell;

  internal::MappingQImplementation::maybe_update_q_points_Jacobians_generic(
    CellSimilarity::none,
    data,
    make_array_view(quadrature.get_points()),
    polynomials_1d,
    renumber_lexicographic_to_hierarchic,
    output_data.quadrature_points,
    output_data.jacobians,
    output_data.inverse_jacobians);

  internal::MappingQImplementation::maybe_update_jacobian_grads<dim, spacedim>(
    CellSimilarity::none,
    data,
    make_array_view(quadrature.get_points()),
    polynomials_1d,
    renumber_lexicographic_to_hierarchic,
    output_data.jacobian_grads);

  internal::MappingQImplementation::maybe_update_jacobian_pushed_forward_grads<
    dim,
    spacedim>(CellSimilarity::none,
              data,
              make_array_view(quadrature.get_points()),
              polynomials_1d,
              renumber_lexicographic_to_hierarchic,
              output_data.jacobian_pushed_forward_grads);

  internal::MappingQImplementation::maybe_update_jacobian_2nd_derivatives<
    dim,
    spacedim>(CellSimilarity::none,
              data,
              make_array_view(quadrature.get_points()),
              polynomials_1d,
              renumber_lexicographic_to_hierarchic,
              output_data.jacobian_2nd_derivatives);

  internal::MappingQImplementation::
    maybe_update_jacobian_pushed_forward_2nd_derivatives<dim, spacedim>(
      CellSimilarity::none,
      data,
      make_array_view(quadrature.get_points()),
      polynomials_1d,
      renumber_lexicographic_to_hierarchic,
      output_data.jacobian_pushed_forward_2nd_derivatives);

  internal::MappingQImplementation::maybe_update_jacobian_3rd_derivatives<
    dim,
    spacedim>(CellSimilarity::none,
              data,
              make_array_view(quadrature.get_points()),
              polynomials_1d,
              renumber_lexicographic_to_hierarchic,
              output_data.jacobian_3rd_derivatives);

  internal::MappingQImplementation::
    maybe_update_jacobian_pushed_forward_3rd_derivatives<dim, spacedim>(
      CellSimilarity::none,
      data,
      make_array_view(quadrature.get_points()),
      polynomials_1d,
      renumber_lexicographic_to_hierarchic,
      output_data.jacobian_pushed_forward_3rd_derivatives);

  const UpdateFlags          update_flags = data.update_each;
  const std::vector<double> &weights      = quadrature.get_weights();

  if ((update_flags & (update_normal_vectors | update_JxW_values)) != 0u)
    {
      AssertDimension(output_data.JxW_values.size(), n_q_points);

      Assert(!(update_flags & update_normal_vectors) ||
               (output_data.normal_vectors.size() == n_q_points),
             ExcDimensionMismatch(output_data.normal_vectors.size(),
                                  n_q_points));


      for (unsigned int point = 0; point < n_q_points; ++point)
        {
          const double det = data.volume_elements[point];

          // check for distorted cells.

          // TODO: this allows for anisotropies of up to 1e6 in 3d and
          // 1e12 in 2d. might want to find a finer
          // (dimension-independent) criterion
          Assert(det > 1e-12 * Utilities::fixed_power<dim>(
                                 cell->diameter() / std::sqrt(double(dim))),
                 (typename Mapping<dim, spacedim>::ExcDistortedMappedCell(
                   cell->center(), det, point)));

          // The normals are n = J^{-T} * \hat{n} before normalizing.
          Tensor<1, spacedim> normal;
          for (unsigned int d = 0; d < spacedim; d++)
            normal[d] = output_data.inverse_jacobians[point].transpose()[d] *
                        quadrature.normal_vector(point);

          output_data.JxW_values[point] = weights[point] * det * normal.norm();

          if ((update_flags & update_normal_vectors) != 0u)
            {
              normal /= normal.norm();
              output_data.normal_vectors[point] = normal;
            }
        }
    }
}



template <int dim, int spacedim>
void
MappingQ<dim, spacedim>::fill_mapping_data_for_generic_points(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const ArrayView<const Point<dim>>                          &unit_points,
  const UpdateFlags                                           update_flags,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  if (update_flags == update_default)
    return;

  Assert(update_flags & update_inverse_jacobians ||
           update_flags & update_jacobians ||
           update_flags & update_quadrature_points,
         ExcNotImplemented());

  output_data.initialize(unit_points.size(), update_flags);

  auto internal_data =
    this->get_data(update_flags,
                   Quadrature<dim>(std::vector<Point<dim>>(unit_points.begin(),
                                                           unit_points.end())));
  const InternalData &data = static_cast<const InternalData &>(*internal_data);
  data.output_data         = &output_data;
  if (polynomial_degree == 1)
    {
      data.mapping_support_points.resize(GeometryInfo<dim>::vertices_per_cell);
      const auto vertices = this->get_vertices(cell);
      for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
        data.mapping_support_points[i] = vertices[i];
    }
  else
    data.mapping_support_points = this->compute_mapping_support_points(cell);

  internal::MappingQImplementation::maybe_update_q_points_Jacobians_generic(
    CellSimilarity::none,
    data,
    unit_points,
    polynomials_1d,
    renumber_lexicographic_to_hierarchic,
    output_data.quadrature_points,
    output_data.jacobians,
    output_data.inverse_jacobians);
}



template <int dim, int spacedim>
void
MappingQ<dim, spacedim>::fill_mapping_data_for_face_quadrature(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const Quadrature<dim - 1>                                  &face_quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  if (face_quadrature.get_points().empty())
    return;

  // ensure that the following static_cast is really correct:
  Assert(dynamic_cast<const InternalData *>(&internal_data) != nullptr,
         ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(internal_data);

  if (polynomial_degree == 1)
    {
      data.mapping_support_points.resize(GeometryInfo<dim>::vertices_per_cell);
      const auto vertices = this->get_vertices(cell);
      for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
        data.mapping_support_points[i] = vertices[i];
    }
  else
    data.mapping_support_points = this->compute_mapping_support_points(cell);
  data.output_data = &output_data;

  internal::MappingQImplementation::do_fill_fe_face_values(
    *this,
    cell,
    face_no,
    numbers::invalid_unsigned_int,
    QProjector<dim>::DataSetDescriptor::cell(),
    face_quadrature,
    data,
    polynomials_1d,
    renumber_lexicographic_to_hierarchic,
    output_data);
}



template <int dim, int spacedim>
void
MappingQ<dim, spacedim>::transform(
  const ArrayView<const Tensor<1, dim>>                   &input,
  const MappingKind                                        mapping_kind,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  const ArrayView<Tensor<1, spacedim>>                    &output) const
{
  internal::MappingQImplementation::transform_fields(input,
                                                     mapping_kind,
                                                     mapping_data,
                                                     output);
}



template <int dim, int spacedim>
void
MappingQ<dim, spacedim>::transform(
  const ArrayView<const DerivativeForm<1, dim, spacedim>> &input,
  const MappingKind                                        mapping_kind,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  const ArrayView<Tensor<2, spacedim>>                    &output) const
{
  internal::MappingQImplementation::transform_differential_forms(input,
                                                                 mapping_kind,
                                                                 mapping_data,
                                                                 output);
}



template <int dim, int spacedim>
void
MappingQ<dim, spacedim>::transform(
  const ArrayView<const Tensor<2, dim>>                   &input,
  const MappingKind                                        mapping_kind,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  const ArrayView<Tensor<2, spacedim>>                    &output) const
{
  switch (mapping_kind)
    {
      case mapping_contravariant:
        internal::MappingQImplementation::transform_fields(input,
                                                           mapping_kind,
                                                           mapping_data,
                                                           output);
        return;

      case mapping_piola_gradient:
      case mapping_contravariant_gradient:
      case mapping_covariant_gradient:
        internal::MappingQImplementation::transform_gradients(input,
                                                              mapping_kind,
                                                              mapping_data,
                                                              output);
        return;
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
}



template <int dim, int spacedim>
void
MappingQ<dim, spacedim>::transform(
  const ArrayView<const DerivativeForm<2, dim, spacedim>> &input,
  const MappingKind                                        mapping_kind,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  const ArrayView<Tensor<3, spacedim>>                    &output) const
{
  AssertDimension(input.size(), output.size());
  Assert(dynamic_cast<const InternalData *>(&mapping_data) != nullptr,
         ExcInternalError());
  const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &data = *static_cast<const InternalData &>(mapping_data).output_data;

  switch (mapping_kind)
    {
      case mapping_covariant_gradient:
        {
          Assert(!data.inverse_jacobians.empty(),
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_covariant_transformation"));

          for (unsigned int q = 0; q < output.size(); ++q)
            {
              const DerivativeForm<1, dim, spacedim> covariant =
                data.inverse_jacobians[q].transpose();
              output[q] =
                internal::apply_covariant_gradient(covariant, input[q]);
            }
          return;
        }

      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
}



template <int dim, int spacedim>
void
MappingQ<dim, spacedim>::transform(
  const ArrayView<const Tensor<3, dim>>                   &input,
  const MappingKind                                        mapping_kind,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  const ArrayView<Tensor<3, spacedim>>                    &output) const
{
  switch (mapping_kind)
    {
      case mapping_piola_hessian:
      case mapping_contravariant_hessian:
      case mapping_covariant_hessian:
        internal::MappingQImplementation::transform_hessians(input,
                                                             mapping_kind,
                                                             mapping_data,
                                                             output);
        return;
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
}



template <int dim, int spacedim>
void
MappingQ<dim, spacedim>::add_line_support_points(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  std::vector<Point<spacedim>>                               &a) const
{
  // if we only need the midpoint, then ask for it.
  if (this->polynomial_degree == 2)
    {
      for (unsigned int line_no = 0;
           line_no < GeometryInfo<dim>::lines_per_cell;
           ++line_no)
        {
          const typename Triangulation<dim, spacedim>::line_iterator line =
            (dim == 1 ?
               static_cast<
                 typename Triangulation<dim, spacedim>::line_iterator>(cell) :
               cell->line(line_no));

          const Manifold<dim, spacedim> &manifold =
            ((line->manifold_id() == numbers::flat_manifold_id) &&
                 (dim < spacedim) ?
               cell->get_manifold() :
               line->get_manifold());
          a.push_back(manifold.get_new_point_on_line(line));
        }
    }
  else
    // otherwise call the more complicated functions and ask for inner points
    // from the manifold description
    {
      std::vector<Point<spacedim>> tmp_points;
      for (unsigned int line_no = 0;
           line_no < GeometryInfo<dim>::lines_per_cell;
           ++line_no)
        {
          const typename Triangulation<dim, spacedim>::line_iterator line =
            (dim == 1 ?
               static_cast<
                 typename Triangulation<dim, spacedim>::line_iterator>(cell) :
               cell->line(line_no));

          const Manifold<dim, spacedim> &manifold =
            ((line->manifold_id() == numbers::flat_manifold_id) &&
                 (dim < spacedim) ?
               cell->get_manifold() :
               line->get_manifold());

          const auto reference_cell = ReferenceCells::get_hypercube<dim>();
          const std::array<Point<spacedim>, 2> vertices{
            {cell->vertex(reference_cell.line_to_cell_vertices(line_no, 0)),
             cell->vertex(reference_cell.line_to_cell_vertices(line_no, 1))}};

          const std::size_t n_rows =
            support_point_weights_perimeter_to_interior[0].size(0);
          a.resize(a.size() + n_rows);
          auto a_view = make_array_view(a.end() - n_rows, a.end());
          manifold.get_new_points(
            make_array_view(vertices.begin(), vertices.end()),
            support_point_weights_perimeter_to_interior[0],
            a_view);
        }
    }
}



template <>
void
MappingQ<3, 3>::add_quad_support_points(
  const Triangulation<3, 3>::cell_iterator &cell,
  std::vector<Point<3>>                    &a) const
{
  const unsigned int faces_per_cell = GeometryInfo<3>::faces_per_cell;

  // used if face quad at boundary or entirely in the interior of the domain
  std::vector<Point<3>> tmp_points;

  // loop over all faces and collect points on them
  for (unsigned int face_no = 0; face_no < faces_per_cell; ++face_no)
    {
      const Triangulation<3>::face_iterator face = cell->face(face_no);

      if constexpr (running_in_debug_mode())
        {
          const bool face_orientation = cell->face_orientation(face_no),
                     face_flip        = cell->face_flip(face_no),
                     face_rotation    = cell->face_rotation(face_no);
          const unsigned int vertices_per_face =
                               GeometryInfo<3>::vertices_per_face,
                             lines_per_face = GeometryInfo<3>::lines_per_face;

          // some sanity checks up front
          for (unsigned int i = 0; i < vertices_per_face; ++i)
            Assert(face->vertex_index(i) ==
                     cell->vertex_index(GeometryInfo<3>::face_to_cell_vertices(
                       face_no, i, face_orientation, face_flip, face_rotation)),
                   ExcInternalError());

          // indices of the lines that bound a face are given by
          // GeometryInfo<3>:: face_to_cell_lines
          for (unsigned int i = 0; i < lines_per_face; ++i)
            Assert(face->line(i) ==
                     cell->line(GeometryInfo<3>::face_to_cell_lines(
                       face_no, i, face_orientation, face_flip, face_rotation)),
                   ExcInternalError());
        }
      // extract the points surrounding a quad from the points
      // already computed. First get the 4 vertices and then the points on
      // the four lines
      boost::container::small_vector<Point<3>, 200> tmp_points(
        GeometryInfo<2>::vertices_per_cell +
        GeometryInfo<2>::lines_per_cell * (polynomial_degree - 1));
      for (const unsigned int v : GeometryInfo<2>::vertex_indices())
        tmp_points[v] = a[GeometryInfo<3>::face_to_cell_vertices(face_no, v)];
      if (polynomial_degree > 1)
        for (unsigned int line = 0; line < GeometryInfo<2>::lines_per_cell;
             ++line)
          for (unsigned int i = 0; i < polynomial_degree - 1; ++i)
            tmp_points[4 + line * (polynomial_degree - 1) + i] =
              a[GeometryInfo<3>::vertices_per_cell +
                (polynomial_degree - 1) *
                  GeometryInfo<3>::face_to_cell_lines(face_no, line) +
                i];

      const std::size_t n_rows =
        support_point_weights_perimeter_to_interior[1].size(0);
      a.resize(a.size() + n_rows);
      auto a_view = make_array_view(a.end() - n_rows, a.end());
      face->get_manifold().get_new_points(
        make_array_view(tmp_points.begin(), tmp_points.end()),
        support_point_weights_perimeter_to_interior[1],
        a_view);
    }
}



template <>
void
MappingQ<2, 3>::add_quad_support_points(
  const Triangulation<2, 3>::cell_iterator &cell,
  std::vector<Point<3>>                    &a) const
{
  std::array<Point<3>, GeometryInfo<2>::vertices_per_cell> vertices;
  for (const unsigned int i : GeometryInfo<2>::vertex_indices())
    vertices[i] = cell->vertex(i);

  Table<2, double> weights(Utilities::fixed_power<2>(polynomial_degree - 1),
                           GeometryInfo<2>::vertices_per_cell);
  for (unsigned int q = 0, q2 = 0; q2 < polynomial_degree - 1; ++q2)
    for (unsigned int q1 = 0; q1 < polynomial_degree - 1; ++q1, ++q)
      {
        const Point<2> point(line_support_points[q1 + 1][0],
                             line_support_points[q2 + 1][0]);
        for (const unsigned int i : GeometryInfo<2>::vertex_indices())
          weights(q, i) = GeometryInfo<2>::d_linear_shape_function(point, i);
      }

  const std::size_t n_rows = weights.size(0);
  a.resize(a.size() + n_rows);
  auto a_view = make_array_view(a.end() - n_rows, a.end());
  cell->get_manifold().get_new_points(
    make_array_view(vertices.begin(), vertices.end()), weights, a_view);
}



template <int dim, int spacedim>
void
MappingQ<dim, spacedim>::add_quad_support_points(
  const typename Triangulation<dim, spacedim>::cell_iterator &,
  std::vector<Point<spacedim>> &) const
{
  DEAL_II_ASSERT_UNREACHABLE();
}



template <int dim, int spacedim>
std::vector<Point<spacedim>>
MappingQ<dim, spacedim>::compute_mapping_support_points(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell) const
{
  // get the vertices first
  std::vector<Point<spacedim>> a;
  a.reserve(Utilities::fixed_power<dim>(polynomial_degree + 1));
  for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
    a.push_back(cell->vertex(i));

  if (this->polynomial_degree > 1)
    {
      // check if all entities have the same manifold id which is when we can
      // simply ask the manifold for all points. the transfinite manifold can
      // do the interpolation better than this class, so if we detect that we
      // do not have to change anything here
      Assert(dim <= 3, ExcImpossibleInDim(dim));
      bool all_manifold_ids_are_equal = (dim == spacedim);
      if (all_manifold_ids_are_equal &&
          dynamic_cast<const TransfiniteInterpolationManifold<dim, spacedim> *>(
            &cell->get_manifold()) == nullptr)
        {
          for (auto f : GeometryInfo<dim>::face_indices())
            if (&cell->face(f)->get_manifold() != &cell->get_manifold())
              all_manifold_ids_are_equal = false;

          if (dim == 3)
            for (unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell; ++l)
              if (&cell->line(l)->get_manifold() != &cell->get_manifold())
                all_manifold_ids_are_equal = false;
        }

      if (all_manifold_ids_are_equal)
        {
          const std::size_t n_rows = support_point_weights_cell.size(0);
          a.resize(a.size() + n_rows);
          auto a_view = make_array_view(a.end() - n_rows, a.end());
          cell->get_manifold().get_new_points(make_array_view(a.begin(),
                                                              a.end() - n_rows),
                                              support_point_weights_cell,
                                              a_view);
        }
      else
        switch (dim)
          {
            case 1:
              add_line_support_points(cell, a);
              break;
            case 2:
              // in 2d, add the points on the four bounding lines to the
              // exterior (outer) points
              add_line_support_points(cell, a);

              // then get the interior support points
              if (dim != spacedim)
                add_quad_support_points(cell, a);
              else
                {
                  const std::size_t n_rows =
                    support_point_weights_perimeter_to_interior[1].size(0);
                  a.resize(a.size() + n_rows);
                  auto a_view = make_array_view(a.end() - n_rows, a.end());
                  cell->get_manifold().get_new_points(
                    make_array_view(a.begin(), a.end() - n_rows),
                    support_point_weights_perimeter_to_interior[1],
                    a_view);
                }
              break;

            case 3:
              // in 3d also add the points located on the boundary faces
              add_line_support_points(cell, a);
              add_quad_support_points(cell, a);

              // then compute the interior points
              {
                const std::size_t n_rows =
                  support_point_weights_perimeter_to_interior[2].size(0);
                a.resize(a.size() + n_rows);
                auto a_view = make_array_view(a.end() - n_rows, a.end());
                cell->get_manifold().get_new_points(
                  make_array_view(a.begin(), a.end() - n_rows),
                  support_point_weights_perimeter_to_interior[2],
                  a_view);
              }
              break;

            default:
              DEAL_II_NOT_IMPLEMENTED();
              break;
          }
    }

  return a;
}



template <int dim, int spacedim>
BoundingBox<spacedim>
MappingQ<dim, spacedim>::get_bounding_box(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell) const
{
  return BoundingBox<spacedim>(this->compute_mapping_support_points(cell));
}



template <int dim, int spacedim>
bool
MappingQ<dim, spacedim>::is_compatible_with(
  const ReferenceCell &reference_cell) const
{
  Assert(dim == reference_cell.get_dimension(),
         ExcMessage("The dimension of your mapping (" +
                    Utilities::to_string(dim) +
                    ") and the reference cell cell_type (" +
                    Utilities::to_string(reference_cell.get_dimension()) +
                    " ) do not agree."));

  return reference_cell.is_hyper_cube();
}



//--------------------------- Explicit instantiations -----------------------
#include "fe/mapping_q.inst"


DEAL_II_NAMESPACE_CLOSE
