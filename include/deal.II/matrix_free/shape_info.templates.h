// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2019 by the deal.II authors
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

#ifndef dealii_matrix_free_shape_info_templates_h
#define dealii_matrix_free_shape_info_templates_h


#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_piecewise.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/utilities.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_poly.h>
#include <deal.II/fe/fe_q_dg0.h>

#include <deal.II/matrix_free/shape_info.h>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace MatrixFreeFunctions
  {
    // ----------------- actual ShapeInfo functions --------------------

    template <typename Number>
    Number
    get_first_array_element(const Number a)
    {
      return a;
    }

    template <typename Number, int width>
    Number
    get_first_array_element(const VectorizedArray<Number, width> a)
    {
      return a[0];
    }

    template <typename Number>
    ShapeInfo<Number>::ShapeInfo()
      : element_type(tensor_general)
      , fe_degree(numbers::invalid_unsigned_int)
      , n_q_points_1d(0)
      , n_q_points(0)
      , dofs_per_component_on_cell(0)
      , n_q_points_face(0)
      , dofs_per_component_on_face(0)
      , nodal_at_cell_boundaries(false)
    {}



    template <typename Number>
    template <int dim>
    void
    ShapeInfo<Number>::reinit(const Quadrature<1> &     quad,
                              const FiniteElement<dim> &fe_in,
                              const unsigned int        base_element_number)
    {
      this->quadrature             = quad;
      const FiniteElement<dim> *fe = &fe_in.base_element(base_element_number);

      Assert(fe->n_components() == 1,
             ExcMessage("FEEvaluation only works for scalar finite elements."));

      fe_degree     = fe->degree;
      n_q_points_1d = quad.size();

      const unsigned int n_dofs_1d = std::min(fe->dofs_per_cell, fe_degree + 1);

      // renumber (this is necessary for FE_Q, for example, since there the
      // vertex DoFs come first, which is incompatible with the lexicographic
      // ordering necessary to apply tensor products efficiently)
      std::vector<unsigned int> scalar_lexicographic;
      Point<dim>                unit_point;
      {
        // find numbering to lexicographic
        Assert(fe->n_components() == 1,
               ExcMessage("Expected a scalar element"));

        const FE_Poly<TensorProductPolynomials<dim>, dim, dim> *fe_poly =
          dynamic_cast<
            const FE_Poly<TensorProductPolynomials<dim>, dim, dim> *>(fe);

        const FE_Poly<
          TensorProductPolynomials<dim,
                                   Polynomials::PiecewisePolynomial<double>>,
          dim,
          dim> *fe_poly_piece =
          dynamic_cast<const FE_Poly<
            TensorProductPolynomials<dim,
                                     Polynomials::PiecewisePolynomial<double>>,
            dim,
            dim> *>(fe);

        const FE_DGP<dim> *fe_dgp = dynamic_cast<const FE_DGP<dim> *>(fe);

        const FE_Q_DG0<dim> *fe_q_dg0 = dynamic_cast<const FE_Q_DG0<dim> *>(fe);

        element_type = tensor_general;
        if (fe_poly != nullptr)
          scalar_lexicographic = fe_poly->get_poly_space_numbering_inverse();
        else if (fe_poly_piece != nullptr)
          scalar_lexicographic =
            fe_poly_piece->get_poly_space_numbering_inverse();
        else if (fe_dgp != nullptr)
          {
            scalar_lexicographic.resize(fe_dgp->dofs_per_cell);
            for (unsigned int i = 0; i < fe_dgp->dofs_per_cell; ++i)
              scalar_lexicographic[i] = i;
            element_type = truncated_tensor;
          }
        else if (fe_q_dg0 != nullptr)
          {
            scalar_lexicographic = fe_q_dg0->get_poly_space_numbering_inverse();
            element_type         = tensor_symmetric_plus_dg0;
          }
        else if (fe->dofs_per_cell == 0)
          {
            // FE_Nothing case -> nothing to do here
          }
        else
          Assert(false, ExcNotImplemented());

        // Finally store the renumbering into the member variable of this
        // class
        if (fe_in.n_components() == 1)
          lexicographic_numbering = scalar_lexicographic;
        else
          {
            // have more than one component, get the inverse
            // permutation, invert it, sort the components one after one,
            // and invert back
            std::vector<unsigned int> scalar_inv =
              Utilities::invert_permutation(scalar_lexicographic);
            std::vector<unsigned int> lexicographic(
              fe_in.dofs_per_cell, numbers::invalid_unsigned_int);
            unsigned int components_before = 0;
            for (unsigned int e = 0; e < base_element_number; ++e)
              components_before += fe_in.element_multiplicity(e);
            for (unsigned int comp = 0;
                 comp < fe_in.element_multiplicity(base_element_number);
                 ++comp)
              for (unsigned int i = 0; i < scalar_inv.size(); ++i)
                lexicographic[fe_in.component_to_system_index(
                  comp + components_before, i)] =
                  scalar_inv.size() * comp + scalar_inv[i];

            // invert numbering again. Need to do it manually because we might
            // have undefined blocks
            lexicographic_numbering.resize(fe_in.element_multiplicity(
                                             base_element_number) *
                                             fe->dofs_per_cell,
                                           numbers::invalid_unsigned_int);
            for (unsigned int i = 0; i < lexicographic.size(); ++i)
              if (lexicographic[i] != numbers::invalid_unsigned_int)
                {
                  AssertIndexRange(lexicographic[i],
                                   lexicographic_numbering.size());
                  lexicographic_numbering[lexicographic[i]] = i;
                }
          }

        // to evaluate 1D polynomials, evaluate along the line with the first
        // unit support point, assuming that fe.shape_value(0,unit_point) ==
        // 1. otherwise, need other entry point (e.g. generating a 1D element
        // by reading the name, as done before r29356)
        if (fe->has_support_points())
          unit_point = fe->get_unit_support_points()[scalar_lexicographic[0]];
        Assert(fe->dofs_per_cell == 0 ||
                 std::abs(fe->shape_value(scalar_lexicographic[0], unit_point) -
                          1) < 1e-13,
               ExcInternalError("Could not decode 1D shape functions for the "
                                "element " +
                                fe->get_name()));
      }

      n_q_points = Utilities::fixed_power<dim>(n_q_points_1d);
      n_q_points_face =
        dim > 1 ? Utilities::fixed_power<dim - 1>(n_q_points_1d) : 1;
      dofs_per_component_on_cell = fe->dofs_per_cell;
      dofs_per_component_on_face =
        dim > 1 ? Utilities::fixed_power<dim - 1>(fe_degree + 1) : 1;

      const unsigned int array_size = n_dofs_1d * n_q_points_1d;
      this->shape_gradients.resize_fast(array_size);
      this->shape_values.resize_fast(array_size);
      this->shape_hessians.resize_fast(array_size);

      this->shape_data_on_face[0].resize(3 * n_dofs_1d);
      this->shape_data_on_face[1].resize(3 * n_dofs_1d);
      this->values_within_subface[0].resize(array_size);
      this->values_within_subface[1].resize(array_size);
      this->gradients_within_subface[0].resize(array_size);
      this->gradients_within_subface[1].resize(array_size);
      this->hessians_within_subface[0].resize(array_size);
      this->hessians_within_subface[1].resize(array_size);

      for (unsigned int i = 0; i < n_dofs_1d; ++i)
        {
          // need to reorder from hierarchical to lexicographic to get the
          // DoFs correct
          const unsigned int my_i = scalar_lexicographic[i];
          for (unsigned int q = 0; q < n_q_points_1d; ++q)
            {
              Point<dim> q_point = unit_point;
              q_point[0]         = quad.get_points()[q][0];

              shape_values[i * n_q_points_1d + q] =
                fe->shape_value(my_i, q_point);
              shape_gradients[i * n_q_points_1d + q] =
                fe->shape_grad(my_i, q_point)[0];
              shape_hessians[i * n_q_points_1d + q] =
                fe->shape_grad_grad(my_i, q_point)[0][0];

              // evaluate basis functions on the two 1D subfaces (i.e., at the
              // positions divided by one half and shifted by one half,
              // respectively)
              q_point[0] *= 0.5;
              values_within_subface[0][i * n_q_points_1d + q] =
                fe->shape_value(my_i, q_point);
              gradients_within_subface[0][i * n_q_points_1d + q] =
                fe->shape_grad(my_i, q_point)[0];
              hessians_within_subface[0][i * n_q_points_1d + q] =
                fe->shape_grad_grad(my_i, q_point)[0][0];
              q_point[0] += 0.5;
              values_within_subface[1][i * n_q_points_1d + q] =
                fe->shape_value(my_i, q_point);
              gradients_within_subface[1][i * n_q_points_1d + q] =
                fe->shape_grad(my_i, q_point)[0];
              hessians_within_subface[1][i * n_q_points_1d + q] =
                fe->shape_grad_grad(my_i, q_point)[0][0];
            }

          // evaluate basis functions on the 1D faces, i.e., in zero and one
          Point<dim> q_point             = unit_point;
          q_point[0]                     = 0;
          this->shape_data_on_face[0][i] = fe->shape_value(my_i, q_point);
          this->shape_data_on_face[0][i + n_dofs_1d] =
            fe->shape_grad(my_i, q_point)[0];
          this->shape_data_on_face[0][i + 2 * n_dofs_1d] =
            fe->shape_grad_grad(my_i, q_point)[0][0];
          q_point[0]                     = 1;
          this->shape_data_on_face[1][i] = fe->shape_value(my_i, q_point);
          this->shape_data_on_face[1][i + n_dofs_1d] =
            fe->shape_grad(my_i, q_point)[0];
          this->shape_data_on_face[1][i + 2 * n_dofs_1d] =
            fe->shape_grad_grad(my_i, q_point)[0][0];
        }

      // get gradient and Hessian transformation matrix for the polynomial
      // space associated with the quadrature rule (collocation space). We
      // need to avoid the case with more than a few hundreds of quadrature
      // points when the Lagrange polynomials constructed in
      // FE_DGQArbitraryNodes underflow.
      if (n_q_points_1d < 200)
        {
          const unsigned int stride = (n_q_points_1d + 1) / 2;
          shape_gradients_collocation_eo.resize(n_q_points_1d * stride);
          shape_hessians_collocation_eo.resize(n_q_points_1d * stride);
          FE_DGQArbitraryNodes<1> fe(quad.get_points());
          for (unsigned int i = 0; i < n_q_points_1d / 2; ++i)
            for (unsigned int q = 0; q < stride; ++q)
              {
                shape_gradients_collocation_eo[i * stride + q] =
                  0.5 *
                  (fe.shape_grad(i, quad.get_points()[q])[0] +
                   fe.shape_grad(i,
                                 quad.get_points()[n_q_points_1d - 1 - q])[0]);
                shape_gradients_collocation_eo[(n_q_points_1d - 1 - i) *
                                                 stride +
                                               q] =
                  0.5 *
                  (fe.shape_grad(i, quad.get_points()[q])[0] -
                   fe.shape_grad(i,
                                 quad.get_points()[n_q_points_1d - 1 - q])[0]);
                shape_hessians_collocation_eo[i * stride + q] =
                  0.5 * (fe.shape_grad_grad(i, quad.get_points()[q])[0][0] +
                         fe.shape_grad_grad(
                           i, quad.get_points()[n_q_points_1d - 1 - q])[0][0]);
                shape_hessians_collocation_eo[(n_q_points_1d - 1 - i) * stride +
                                              q] =
                  0.5 * (fe.shape_grad_grad(i, quad.get_points()[q])[0][0] -
                         fe.shape_grad_grad(
                           i, quad.get_points()[n_q_points_1d - 1 - q])[0][0]);
              }
          if (n_q_points_1d % 2 == 1)
            for (unsigned int q = 0; q < stride; ++q)
              {
                shape_gradients_collocation_eo[n_q_points_1d / 2 * stride + q] =
                  fe.shape_grad(n_q_points_1d / 2, quad.get_points()[q])[0];
                shape_hessians_collocation_eo[n_q_points_1d / 2 * stride + q] =
                  fe.shape_grad_grad(n_q_points_1d / 2,
                                     quad.get_points()[q])[0][0];
              }
        }

      if (element_type == tensor_general &&
          check_1d_shapes_symmetric(n_q_points_1d))
        {
          if (check_1d_shapes_collocation())
            element_type = tensor_symmetric_collocation;
          else
            element_type = tensor_symmetric;
          if (n_dofs_1d > 2 && element_type == tensor_symmetric)
            {
              // check if we are a Hermite type
              element_type = tensor_symmetric_hermite;
              for (unsigned int i = 1; i < n_dofs_1d; ++i)
                if (std::abs(get_first_array_element(
                      shape_data_on_face[0][i])) > 1e-12)
                  element_type = tensor_symmetric;
              for (unsigned int i = 2; i < n_dofs_1d; ++i)
                if (std::abs(get_first_array_element(
                      shape_data_on_face[0][n_dofs_1d + i])) > 1e-12)
                  element_type = tensor_symmetric;
            }
        }
      else if (element_type == tensor_symmetric_plus_dg0)
        check_1d_shapes_symmetric(n_q_points_1d);

      nodal_at_cell_boundaries = true;
      for (unsigned int i = 1; i < n_dofs_1d; ++i)
        if (std::abs(get_first_array_element(shape_data_on_face[0][i])) >
              1e-13 ||
            std::abs(get_first_array_element(shape_data_on_face[1][i - 1])) >
              1e-13)
          nodal_at_cell_boundaries = false;

      if (nodal_at_cell_boundaries == true)
        {
          face_to_cell_index_nodal.reinit(GeometryInfo<dim>::faces_per_cell,
                                          dofs_per_component_on_face);
          for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
            {
              const unsigned int direction = f / 2;
              const unsigned int stride =
                direction < dim - 1 ? (fe_degree + 1) : 1;
              int shift = 1;
              for (unsigned int d = 0; d < direction; ++d)
                shift *= fe_degree + 1;
              const unsigned int offset = (f % 2) * fe_degree * shift;

              if (direction == 0 || direction == dim - 1)
                for (unsigned int i = 0; i < dofs_per_component_on_face; ++i)
                  face_to_cell_index_nodal(f, i) = offset + i * stride;
              else
                // local coordinate system on faces 2 and 3 is zx in
                // deal.II, not xz as expected for tensor products -> adjust
                // that here
                for (unsigned int j = 0; j <= fe_degree; ++j)
                  for (unsigned int i = 0; i <= fe_degree; ++i)
                    {
                      const unsigned int ind =
                        offset + j * dofs_per_component_on_face + i;
                      AssertIndexRange(ind, dofs_per_component_on_cell);
                      const unsigned int l           = i * (fe_degree + 1) + j;
                      face_to_cell_index_nodal(f, l) = ind;
                    }
            }
        }

      if (element_type == tensor_symmetric_hermite)
        {
          face_to_cell_index_hermite.reinit(GeometryInfo<dim>::faces_per_cell,
                                            2 * dofs_per_component_on_face);
          for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
            {
              const unsigned int direction = f / 2;
              const unsigned int stride =
                direction < dim - 1 ? (fe_degree + 1) : 1;
              int shift = 1;
              for (unsigned int d = 0; d < direction; ++d)
                shift *= fe_degree + 1;
              const unsigned int offset = (f % 2) * fe_degree * shift;
              if (f % 2 == 1)
                shift = -shift;

              if (direction == 0 || direction == dim - 1)
                for (unsigned int i = 0; i < dofs_per_component_on_face; ++i)
                  {
                    face_to_cell_index_hermite(f, 2 * i) = offset + i * stride;
                    face_to_cell_index_hermite(f, 2 * i + 1) =
                      offset + i * stride + shift;
                  }
              else
                // local coordinate system on faces 2 and 3 is zx in
                // deal.II, not xz as expected for tensor products -> adjust
                // that here
                for (unsigned int j = 0; j <= fe_degree; ++j)
                  for (unsigned int i = 0; i <= fe_degree; ++i)
                    {
                      const unsigned int ind =
                        offset + j * dofs_per_component_on_face + i;
                      AssertIndexRange(ind, dofs_per_component_on_cell);
                      const unsigned int l = i * (fe_degree + 1) + j;
                      face_to_cell_index_hermite(f, 2 * l)     = ind;
                      face_to_cell_index_hermite(f, 2 * l + 1) = ind + shift;
                    }
            }
        }
    }



    template <typename Number>
    bool
    ShapeInfo<Number>::check_1d_shapes_symmetric(
      const unsigned int n_q_points_1d)
    {
      if (dofs_per_component_on_cell == 0)
        return false;

      const double zero_tol =
        std::is_same<Number, double>::value == true ? 1e-12 : 1e-7;
      // symmetry for values
      const unsigned int n_dofs_1d = fe_degree + 1;
      for (unsigned int i = 0; i < (n_dofs_1d + 1) / 2; ++i)
        for (unsigned int j = 0; j < n_q_points_1d; ++j)
          if (std::abs(get_first_array_element(
                shape_values[i * n_q_points_1d + j] -
                shape_values[(n_dofs_1d - i) * n_q_points_1d - j - 1])) >
              std::max(zero_tol,
                       zero_tol * std::abs(get_first_array_element(
                                    shape_values[i * n_q_points_1d + j]))))
            return false;

      // shape values should be zero at x=0.5 for all basis functions except
      // for the middle one
      if (n_q_points_1d % 2 == 1 && n_dofs_1d % 2 == 1)
        {
          for (unsigned int i = 0; i < n_dofs_1d / 2; ++i)
            if (std::abs(get_first_array_element(
                  shape_values[i * n_q_points_1d + n_q_points_1d / 2])) >
                zero_tol)
              return false;
        }

      // skew-symmetry for gradient, zero of middle basis function in middle
      // quadrature point. Multiply tolerance by degree of the element to
      // the power of 1.5 to get a suitable gradient scaling
      const double zero_tol_gradient =
        zero_tol * std::sqrt(fe_degree + 1.) * (fe_degree + 1);
      for (unsigned int i = 0; i < (n_dofs_1d + 1) / 2; ++i)
        for (unsigned int j = 0; j < n_q_points_1d; ++j)
          if (std::abs(get_first_array_element(
                shape_gradients[i * n_q_points_1d + j] +
                shape_gradients[(n_dofs_1d - i) * n_q_points_1d - j - 1])) >
              zero_tol_gradient)
            return false;
      if (n_dofs_1d % 2 == 1 && n_q_points_1d % 2 == 1)
        if (std::abs(get_first_array_element(
              shape_gradients[(n_dofs_1d / 2) * n_q_points_1d +
                              (n_q_points_1d / 2)])) > zero_tol_gradient)
          return false;

      // symmetry for Hessian. Multiply tolerance by degree^3 of the element
      // to get a suitable Hessian scaling
      const double zero_tol_hessian =
        zero_tol * (fe_degree + 1) * (fe_degree + 1) * (fe_degree + 1);
      for (unsigned int i = 0; i < (n_dofs_1d + 1) / 2; ++i)
        for (unsigned int j = 0; j < n_q_points_1d; ++j)
          if (std::abs(get_first_array_element(
                shape_hessians[i * n_q_points_1d + j] -
                shape_hessians[(n_dofs_1d - i) * n_q_points_1d - j - 1])) >
              zero_tol_hessian)
            return false;

      const unsigned int stride = (n_q_points_1d + 1) / 2;
      shape_values_eo.resize((fe_degree + 1) * stride);
      shape_gradients_eo.resize((fe_degree + 1) * stride);
      shape_hessians_eo.resize((fe_degree + 1) * stride);
      for (unsigned int i = 0; i < (fe_degree + 1) / 2; ++i)
        for (unsigned int q = 0; q < stride; ++q)
          {
            shape_values_eo[i * stride + q] =
              0.5 * (shape_values[i * n_q_points_1d + q] +
                     shape_values[i * n_q_points_1d + n_q_points_1d - 1 - q]);
            shape_values_eo[(fe_degree - i) * stride + q] =
              0.5 * (shape_values[i * n_q_points_1d + q] -
                     shape_values[i * n_q_points_1d + n_q_points_1d - 1 - q]);

            shape_gradients_eo[i * stride + q] =
              0.5 *
              (shape_gradients[i * n_q_points_1d + q] +
               shape_gradients[i * n_q_points_1d + n_q_points_1d - 1 - q]);
            shape_gradients_eo[(fe_degree - i) * stride + q] =
              0.5 *
              (shape_gradients[i * n_q_points_1d + q] -
               shape_gradients[i * n_q_points_1d + n_q_points_1d - 1 - q]);

            shape_hessians_eo[i * stride + q] =
              0.5 * (shape_hessians[i * n_q_points_1d + q] +
                     shape_hessians[i * n_q_points_1d + n_q_points_1d - 1 - q]);
            shape_hessians_eo[(fe_degree - i) * stride + q] =
              0.5 * (shape_hessians[i * n_q_points_1d + q] -
                     shape_hessians[i * n_q_points_1d + n_q_points_1d - 1 - q]);
          }
      if (fe_degree % 2 == 0)
        for (unsigned int q = 0; q < stride; ++q)
          {
            shape_values_eo[fe_degree / 2 * stride + q] =
              shape_values[(fe_degree / 2) * n_q_points_1d + q];
            shape_gradients_eo[fe_degree / 2 * stride + q] =
              shape_gradients[(fe_degree / 2) * n_q_points_1d + q];
            shape_hessians_eo[fe_degree / 2 * stride + q] =
              shape_hessians[(fe_degree / 2) * n_q_points_1d + q];
          }

      return true;
    }



    template <typename Number>
    bool
    ShapeInfo<Number>::check_1d_shapes_collocation() const
    {
      if (dofs_per_component_on_cell != n_q_points)
        return false;

      const double zero_tol =
        std::is_same<Number, double>::value == true ? 1e-12 : 1e-7;
      // check: identity operation for shape values
      const unsigned int n_points_1d = fe_degree + 1;
      for (unsigned int i = 0; i < n_points_1d; ++i)
        for (unsigned int j = 0; j < n_points_1d; ++j)
          if (i != j)
            {
              if (std::abs(get_first_array_element(
                    shape_values[i * n_points_1d + j])) > zero_tol)
                return false;
            }
          else
            {
              if (std::abs(
                    get_first_array_element(shape_values[i * n_points_1d + j]) -
                    1.) > zero_tol)
                return false;
            }
      return true;
    }



    template <typename Number>
    std::size_t
    ShapeInfo<Number>::memory_consumption() const
    {
      std::size_t memory = sizeof(*this);
      memory += MemoryConsumption::memory_consumption(shape_values);
      memory += MemoryConsumption::memory_consumption(shape_gradients);
      memory += MemoryConsumption::memory_consumption(shape_hessians);
      memory += MemoryConsumption::memory_consumption(shape_values_eo);
      memory += MemoryConsumption::memory_consumption(shape_gradients_eo);
      memory += MemoryConsumption::memory_consumption(shape_hessians_eo);
      memory +=
        MemoryConsumption::memory_consumption(shape_gradients_collocation_eo);
      memory +=
        MemoryConsumption::memory_consumption(shape_hessians_collocation_eo);
      for (unsigned int i = 0; i < 2; ++i)
        {
          memory +=
            MemoryConsumption::memory_consumption(shape_data_on_face[i]);
          memory +=
            MemoryConsumption::memory_consumption(values_within_subface[i]);
          memory +=
            MemoryConsumption::memory_consumption(gradients_within_subface[i]);
        }
      return memory;
    }

    // end of functions for ShapeInfo

  } // end of namespace MatrixFreeFunctions
} // end of namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif
