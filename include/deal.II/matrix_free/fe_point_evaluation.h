// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_point_evaluation_h
#define dealii_fe_point_evaluation_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/matrix_free/evaluation_flags.h>
#include <deal.II/matrix_free/evaluation_kernels_face.h>
#include <deal.II/matrix_free/mapping_info_storage.h>
#include <deal.II/matrix_free/shape_info.h>
#include <deal.II/matrix_free/tensor_product_point_kernels.h>

#include <deal.II/non_matching/mapping_info.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace FEPointEvaluation
  {
    DeclException1(
      ExcFEPointEvaluationAccessToUninitializedMappingField,
      std::string,
      << "You are requesting information from an FEPointEvaluationBase "
      << "object for which this kind of information has not been computed. "
      << "What information these objects compute is determined by the update_* "
      << "flags you pass to MappingInfo() in the Constructor. Here, "
      << "the operation you are attempting requires the <" << arg1
      << "> flag to be set, but it was apparently not specified "
      << "upon initialization.");

    /**
     * Struct to distinguish between the value and gradient types of different
     * numbers of components used by the FlexibleEvaluator class.
     */
    template <int dim,
              int spacedim,
              int n_components,
              typename Number,
              typename Enable = void>
    struct EvaluatorTypeTraits
    {
      using ScalarNumber =
        typename internal::VectorizedArrayTrait<Number>::value_type;
      using VectorizedArrayType =
        typename dealii::internal::VectorizedArrayTrait<
          Number>::vectorized_value_type;
      using value_type        = Tensor<1, n_components, Number>;
      using scalar_value_type = Tensor<1, n_components, ScalarNumber>;
      using vectorized_value_type =
        Tensor<1, n_components, VectorizedArrayType>;
      using unit_gradient_type =
        Tensor<1, n_components, Tensor<1, dim, Number>>;
      using real_gradient_type = std::conditional_t<
        n_components == spacedim,
        Tensor<2, spacedim, Number>,
        Tensor<1, n_components, Tensor<1, spacedim, Number>>>;
      using scalar_unit_gradient_type =
        Tensor<1, n_components, Tensor<1, dim, ScalarNumber>>;
      using vectorized_unit_gradient_type =
        Tensor<1, n_components, Tensor<1, dim, VectorizedArrayType>>;
      using interface_vectorized_unit_gradient_type =
        Tensor<1, dim, Tensor<1, n_components, VectorizedArrayType>>;

      static void
      read_value(const ScalarNumber vector_entry,
                 const unsigned int component,
                 scalar_value_type &result)
      {
        AssertIndexRange(component, n_components);
        result[component] = vector_entry;
      }

      static scalar_value_type
      sum_value(const scalar_value_type &result)
      {
        return result;
      }

      static scalar_value_type
      sum_value(const vectorized_value_type &result)
      {
        scalar_value_type result_scalar = {};

        for (unsigned int c = 0; c < n_components; ++c)
          result_scalar[c] = result[c].sum();

        return result_scalar;
      }

      static ScalarNumber
      sum_value(const unsigned int           component,
                const vectorized_value_type &result)
      {
        AssertIndexRange(component, n_components);
        return result[component].sum();
      }

      static void
      set_gradient(const interface_vectorized_unit_gradient_type &value,
                   const unsigned int                             vector_lane,
                   unit_gradient_type                            &result)
      {
        for (unsigned int i = 0; i < n_components; ++i)
          for (unsigned int d = 0; d < dim; ++d)
            result[i][d] =
              internal::VectorizedArrayTrait<Number>::get_from_vectorized(
                value[d][i], vector_lane);
      }

      static void
      get_gradient(interface_vectorized_unit_gradient_type &value,
                   const unsigned int                       vector_lane,
                   const unit_gradient_type                &result)
      {
        for (unsigned int i = 0; i < n_components; ++i)
          for (unsigned int d = 0; d < dim; ++d)
            internal::VectorizedArrayTrait<Number>::get_from_vectorized(
              value[d][i], vector_lane) = result[i][d];
      }

      static void
      get_gradient(interface_vectorized_unit_gradient_type &value,
                   const unsigned int                       vector_lane,
                   const DerivativeForm<1, dim, n_components, Number> &result)
      {
        for (unsigned int i = 0; i < n_components; ++i)
          for (unsigned int d = 0; d < dim; ++d)
            internal::VectorizedArrayTrait<Number>::get_from_vectorized(
              value[d][i], vector_lane) = result[i][d];
      }

      static void
      set_zero_gradient(real_gradient_type &value,
                        const unsigned int  vector_lane)
      {
        for (unsigned int i = 0; i < n_components; ++i)
          for (unsigned int d = 0; d < spacedim; ++d)
            internal::VectorizedArrayTrait<Number>::get(value[i][d],
                                                        vector_lane) = 0.;
      }

      static void
      set_value(const vectorized_value_type &value,
                const unsigned int           vector_lane,
                scalar_value_type           &result)
      {
        for (unsigned int i = 0; i < n_components; ++i)
          result[i] = value[i][vector_lane];
      }

      static void
      set_value(const vectorized_value_type &value,
                const unsigned int,
                vectorized_value_type &result)
      {
        result = value;
      }

      static void
      get_value(vectorized_value_type   &value,
                const unsigned int       vector_lane,
                const scalar_value_type &result)
      {
        for (unsigned int i = 0; i < n_components; ++i)
          value[i][vector_lane] = result[i];
      }

      static void
      get_value(vectorized_value_type &value,
                const unsigned int,
                const vectorized_value_type &result)
      {
        value = result;
      }

      static void
      set_zero_value(value_type &value, const unsigned int vector_lane)
      {
        for (unsigned int i = 0; i < n_components; ++i)
          internal::VectorizedArrayTrait<Number>::get(value[i], vector_lane) =
            0.;
      }

      static void
      access(value_type         &value,
             const unsigned int  vector_lane,
             const unsigned int  component,
             const ScalarNumber &shape_value)
      {
        internal::VectorizedArrayTrait<Number>::get(value[component],
                                                    vector_lane) += shape_value;
      }

      static ScalarNumber
      access(const value_type  &value,
             const unsigned int vector_lane,
             const unsigned int component)
      {
        return internal::VectorizedArrayTrait<Number>::get(value[component],
                                                           vector_lane);
      }

      static void
      access(real_gradient_type                      &value,
             const unsigned int                       vector_lane,
             const unsigned int                       component,
             const Tensor<1, spacedim, ScalarNumber> &shape_gradient)
      {
        for (unsigned int d = 0; d < spacedim; ++d)
          internal::VectorizedArrayTrait<Number>::get(value[component][d],
                                                      vector_lane) +=
            shape_gradient[d];
      }

      static Tensor<1, spacedim, ScalarNumber>
      access(const real_gradient_type &value,
             const unsigned int        vector_lane,
             const unsigned int        component)
      {
        Tensor<1, spacedim, ScalarNumber> result;
        for (unsigned int d = 0; d < spacedim; ++d)
          result[d] =
            internal::VectorizedArrayTrait<Number>::get(value[component][d],
                                                        vector_lane);
        return result;
      }
    };

    template <int dim, int spacedim, typename Number>
    struct EvaluatorTypeTraits<dim, spacedim, 1, Number>
    {
      using ScalarNumber =
        typename internal::VectorizedArrayTrait<Number>::value_type;
      using VectorizedArrayType =
        typename dealii::internal::VectorizedArrayTrait<
          Number>::vectorized_value_type;
      using value_type                    = Number;
      using scalar_value_type             = ScalarNumber;
      using vectorized_value_type         = VectorizedArrayType;
      using unit_gradient_type            = Tensor<1, dim, Number>;
      using real_gradient_type            = Tensor<1, spacedim, Number>;
      using scalar_unit_gradient_type     = Tensor<1, dim, ScalarNumber>;
      using vectorized_unit_gradient_type = Tensor<1, dim, VectorizedArrayType>;
      using interface_vectorized_unit_gradient_type =
        vectorized_unit_gradient_type;

      static void
      read_value(const ScalarNumber vector_entry,
                 const unsigned int,
                 scalar_value_type &result)
      {
        result = vector_entry;
      }

      static scalar_value_type
      sum_value(const scalar_value_type &result)
      {
        return result;
      }

      static scalar_value_type
      sum_value(const vectorized_value_type &result)
      {
        return result.sum();
      }

      static ScalarNumber
      sum_value(const unsigned int, const vectorized_value_type &result)
      {
        return result.sum();
      }

      static void
      set_gradient(const vectorized_unit_gradient_type &value,
                   const unsigned int                   vector_lane,
                   scalar_unit_gradient_type           &result)
      {
        for (unsigned int d = 0; d < dim; ++d)
          result[d] = value[d][vector_lane];
      }

      static void
      set_gradient(const vectorized_unit_gradient_type &value,
                   const unsigned int,
                   vectorized_unit_gradient_type &result)
      {
        result = value;
      }

      static void
      get_gradient(vectorized_unit_gradient_type   &value,
                   const unsigned int               vector_lane,
                   const scalar_unit_gradient_type &result)
      {
        for (unsigned int d = 0; d < dim; ++d)
          value[d][vector_lane] = result[d];
      }

      static void
      get_gradient(vectorized_unit_gradient_type &value,
                   const unsigned int,
                   const vectorized_unit_gradient_type &result)
      {
        value = result;
      }

      static void
      set_zero_gradient(real_gradient_type &value,
                        const unsigned int  vector_lane)
      {
        for (unsigned int d = 0; d < spacedim; ++d)
          internal::VectorizedArrayTrait<Number>::get(value[d], vector_lane) =
            0.;
      }

      static void
      set_value(const vectorized_value_type &value,
                const unsigned int           vector_lane,
                scalar_value_type           &result)
      {
        result = value[vector_lane];
      }

      static void
      set_value(const vectorized_value_type &value,
                const unsigned int,
                vectorized_value_type &result)
      {
        result = value;
      }

      static void
      get_value(vectorized_value_type   &value,
                const unsigned int       vector_lane,
                const scalar_value_type &result)
      {
        value[vector_lane] = result;
      }

      static void
      get_value(vectorized_value_type &value,
                const unsigned int,
                const vectorized_value_type &result)
      {
        value = result;
      }

      static void
      set_zero_value(value_type &value, const unsigned int vector_lane)
      {
        internal::VectorizedArrayTrait<Number>::get(value, vector_lane) = 0.;
      }

      static void
      access(value_type        &value,
             const unsigned int vector_lane,
             const unsigned int,
             const ScalarNumber &shape_value)
      {
        internal::VectorizedArrayTrait<Number>::get(value, vector_lane) +=
          shape_value;
      }

      static ScalarNumber
      access(const value_type  &value,
             const unsigned int vector_lane,
             const unsigned int)
      {
        return internal::VectorizedArrayTrait<Number>::get(value, vector_lane);
      }

      static void
      access(real_gradient_type &value,
             const unsigned int  vector_lane,
             const unsigned int,
             const Tensor<1, spacedim, ScalarNumber> &shape_gradient)
      {
        for (unsigned int d = 0; d < spacedim; ++d)
          internal::VectorizedArrayTrait<Number>::get(value[d], vector_lane) +=
            shape_gradient[d];
      }

      static Tensor<1, spacedim, ScalarNumber>
      access(const real_gradient_type &value,
             const unsigned int        vector_lane,
             const unsigned int)
      {
        Tensor<1, spacedim, ScalarNumber> result;
        for (unsigned int d = 0; d < spacedim; ++d)
          result[d] =
            internal::VectorizedArrayTrait<Number>::get(value[d], vector_lane);
        return result;
      }
    };

    template <int dim, typename Number>
    struct EvaluatorTypeTraits<dim,
                               dim,
                               dim,
                               Number,
                               std::enable_if_t<dim != 1>>
    {
      using ScalarNumber =
        typename internal::VectorizedArrayTrait<Number>::value_type;
      using VectorizedArrayType =
        typename dealii::internal::VectorizedArrayTrait<
          Number>::vectorized_value_type;
      using value_type                    = Tensor<1, dim, Number>;
      using scalar_value_type             = Tensor<1, dim, ScalarNumber>;
      using vectorized_value_type         = Tensor<1, dim, VectorizedArrayType>;
      using unit_gradient_type            = Tensor<2, dim, Number>;
      using real_gradient_type            = unit_gradient_type;
      using scalar_unit_gradient_type     = Tensor<2, dim, ScalarNumber>;
      using vectorized_unit_gradient_type = Tensor<2, dim, VectorizedArrayType>;
      using interface_vectorized_unit_gradient_type =
        Tensor<1, dim, Tensor<1, dim, VectorizedArrayType>>;

      static void
      read_value(const ScalarNumber vector_entry,
                 const unsigned int component,
                 scalar_value_type &result)
      {
        AssertIndexRange(component, dim);
        result[component] = vector_entry;
      }

      static scalar_value_type
      sum_value(const scalar_value_type &result)
      {
        return result;
      }

      static scalar_value_type
      sum_value(const vectorized_value_type &result)
      {
        scalar_value_type result_scalar = {};

        for (unsigned int c = 0; c < dim; ++c)
          result_scalar[c] = result[c].sum();

        return result_scalar;
      }

      static ScalarNumber
      sum_value(const unsigned int           component,
                const vectorized_value_type &result)
      {
        AssertIndexRange(component, dim);
        return result[component].sum();
      }

      static void
      set_gradient(const interface_vectorized_unit_gradient_type &value,
                   const unsigned int                             vector_lane,
                   unit_gradient_type                            &result)
      {
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int d = 0; d < dim; ++d)
            result[i][d] =
              internal::VectorizedArrayTrait<Number>::get_from_vectorized(
                value[d][i], vector_lane);
      }

      static void
      get_gradient(interface_vectorized_unit_gradient_type &value,
                   const unsigned int                       vector_lane,
                   const unit_gradient_type                &result)
      {
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int d = 0; d < dim; ++d)
            internal::VectorizedArrayTrait<Number>::get_from_vectorized(
              value[d][i], vector_lane) = result[i][d];
      }

      static void
      set_zero_gradient(unit_gradient_type &value,
                        const unsigned int  vector_lane)
      {
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int d = 0; d < dim; ++d)
            internal::VectorizedArrayTrait<Number>::get(value[i][d],
                                                        vector_lane) = 0.;
      }

      static void
      set_value(const vectorized_value_type &value,
                const unsigned int           vector_lane,
                scalar_value_type           &result)
      {
        for (unsigned int i = 0; i < dim; ++i)
          result[i] = value[i][vector_lane];
      }

      static void
      set_value(const vectorized_value_type &value,
                const unsigned int,
                vectorized_value_type &result)
      {
        result = value;
      }

      static void
      get_value(vectorized_value_type   &value,
                const unsigned int       vector_lane,
                const scalar_value_type &result)
      {
        for (unsigned int i = 0; i < dim; ++i)
          value[i][vector_lane] = result[i];
      }

      static void
      get_value(vectorized_value_type &value,
                const unsigned int,
                const vectorized_value_type &result)
      {
        value = result;
      }

      static void
      set_zero_value(value_type &value, const unsigned int vector_lane)
      {
        for (unsigned int i = 0; i < dim; ++i)
          internal::VectorizedArrayTrait<Number>::get(value[i], vector_lane) =
            0.;
      }

      static void
      access(value_type         &value,
             const unsigned int  vector_lane,
             const unsigned int  component,
             const ScalarNumber &shape_value)
      {
        internal::VectorizedArrayTrait<Number>::get(value[component],
                                                    vector_lane) += shape_value;
      }

      static ScalarNumber
      access(const value_type  &value,
             const unsigned int vector_lane,
             const unsigned int component)
      {
        return internal::VectorizedArrayTrait<Number>::get(value[component],
                                                           vector_lane);
      }

      static void
      access(real_gradient_type                 &value,
             const unsigned int                  vector_lane,
             const unsigned int                  component,
             const Tensor<1, dim, ScalarNumber> &shape_gradient)
      {
        for (unsigned int d = 0; d < dim; ++d)
          internal::VectorizedArrayTrait<Number>::get(value[component][d],
                                                      vector_lane) +=
            shape_gradient[d];
      }

      static Tensor<1, dim, ScalarNumber>
      access(const real_gradient_type &value,
             const unsigned int        vector_lane,
             const unsigned int        component)
      {
        Tensor<1, dim, ScalarNumber> result;
        for (unsigned int d = 0; d < dim; ++d)
          result[d] =
            internal::VectorizedArrayTrait<Number>::get(value[component][d],
                                                        vector_lane);
        return result;
      }
    };

    template <int dim, int spacedim>
    bool
    is_fast_path_supported(const FiniteElement<dim, spacedim> &fe,
                           const unsigned int base_element_number);

    template <int dim, int spacedim>
    bool
    is_fast_path_supported(const Mapping<dim, spacedim> &mapping);

    template <int dim, int spacedim>
    std::vector<Polynomials::Polynomial<double>>
    get_polynomial_space(const FiniteElement<dim, spacedim> &fe);
  } // namespace FEPointEvaluation
} // namespace internal



/**
 * Base class of FEPointEvaluation and FEFacePointEvaluation. This class needs
 * usually not be called in user code and does not have any public
 * constructor. The usage is through the class
 * FEPointEvaluation/FEFacePointEvaluation instead.
 */
template <int n_components_,
          int dim,
          int spacedim    = dim,
          typename Number = double>
class FEPointEvaluationBase
{
public:
  static constexpr unsigned int dimension    = dim;
  static constexpr unsigned int n_components = n_components_;

  using number_type = Number;

  using ScalarNumber =
    typename internal::VectorizedArrayTrait<Number>::value_type;
  using VectorizedArrayType = typename dealii::internal::VectorizedArrayTrait<
    Number>::vectorized_value_type;
  using ETT = typename internal::FEPointEvaluation::
    EvaluatorTypeTraits<dim, spacedim, n_components, Number>;
  using value_type            = typename ETT::value_type;
  using scalar_value_type     = typename ETT::scalar_value_type;
  using vectorized_value_type = typename ETT::vectorized_value_type;
  using gradient_type         = typename ETT::real_gradient_type;
  using interface_vectorized_unit_gradient_type =
    typename ETT::interface_vectorized_unit_gradient_type;

protected:
  /**
   * Constructor.
   *
   * @param mapping The Mapping class describing the actual geometry of a cell
   * passed to the evaluate() function.
   *
   * @param fe The FiniteElement object that is used for the evaluation, which
   * is typically the same on all cells to be evaluated.
   *
   * @param update_flags Specify the quantities to be computed by the mapping
   * during the call of reinit(). During evaluate() or integrate(), this data
   * is queried to produce the desired result (e.g., the gradient of a finite
   * element solution).
   *
   * @param first_selected_component For multi-component FiniteElement
   * objects, this parameter allows to select a range of `n_components`
   * components starting from this parameter.
   */
  FEPointEvaluationBase(const Mapping<dim, spacedim>       &mapping,
                        const FiniteElement<dim, spacedim> &fe,
                        const UpdateFlags                   update_flags,
                        const unsigned int first_selected_component = 0);

  /**
   * Constructor to make the present class able to re-use the geometry
   * data also used by other `FEPointEvaluationBase` objects.
   *
   * @param mapping_info The MappingInfo class describes the geometry-related
   * data for evaluating finite-element solutions. This object enables to
   * construct such an object on the outside, possibly re-using it between
   * several objects or between several calls to the same cell and unit points.
   *
   * @param fe The FiniteElement object that is used for the evaluation, which
   * is typically the same on all cells to be evaluated.
   *
   * @param first_selected_component For multi-component FiniteElement
   * objects, this parameter allows to select a range of `n_components`
   * components starting from this parameter.
   *
   * @param is_interior Defines if interior or exterior. Only makes sense for
   * faces.
   */
  FEPointEvaluationBase(
    const NonMatching::MappingInfo<dim, spacedim, Number> &mapping_info,
    const FiniteElement<dim, spacedim>                    &fe,
    const unsigned int first_selected_component = 0,
    const bool         is_interior              = true);

  /**
   * Copy constructor.
   */
  FEPointEvaluationBase(FEPointEvaluationBase &other) noexcept;

  /**
   * Move constructor.
   */
  FEPointEvaluationBase(FEPointEvaluationBase &&other) noexcept;

public:
  /**
   * Return the value at quadrature point number @p point_index after a call to
   * FEPointEvaluation::evaluate() with EvaluationFlags::values set, or
   * the value that has been stored there with a call to
   * FEPointEvaluationBase::submit_value(). If the object is vector-valued, a
   * vector-valued return argument is given.
   */
  const value_type &
  get_value(const unsigned int point_index) const;

  /**
   * Write a value to the field containing the values on points
   * with component point_index. Access to the same field as through
   * get_value(). If applied before the function
   * FEPointEvaluation::integrate() with EvaluationFlags::values set is
   * called, this specifies the value which is tested by all basis function on
   * the current cell and integrated over.
   */
  void
  submit_value(const value_type &value, const unsigned int point_index);

  /**
   * Return the gradient in real coordinates at the point with index
   * `point_index` after a call to FEPointEvaluation::evaluate() with
   * EvaluationFlags::gradients set, or the gradient that has been stored there
   * with a call to FEPointEvaluationBase::submit_gradient(). The gradient in
   * real coordinates is obtained by taking the unit gradient (also accessible
   * via get_unit_gradient()) and applying the inverse Jacobian of the mapping.
   * If the object is vector-valued, a vector-valued return argument is given.
   */
  const gradient_type &
  get_gradient(const unsigned int point_index) const;

  /**
   * Write a contribution that is tested by the gradient to the field
   * containing the values on points with the given `point_index`. Access to
   * the same field as through get_gradient(). If applied before the function
   * FEPointEvaluation::integrate(EvaluationFlags::gradients) is called,
   * this specifies what is tested by all basis function gradients on the
   * current cell and integrated over.
   */
  void
  submit_gradient(const gradient_type &, const unsigned int point_index);

  /**
   * Return the divergence in real coordinates at the point with index
   * `point_index` after a call to FEPointEvaluation::evaluate() with
   * EvaluationFlags::gradients set, or the divergence that has been stored
   * there with a call to FEPointEvaluationBase::submit_divergence(). This
   * functions only makes sense for a vector field with dim components.
   */
  Number
  get_divergence(const unsigned int point_index) const;

  /**
   * Write a contribution that is tested by the divergence to the field
   * containing the values on points with the given `point_index`. Access to
   * the same field as through get_divergence(). If applied before the function
   * FEPointEvaluation::integrate(EvaluationFlags::gradients) is called,
   * this specifies what is tested by all basis function gradients on the
   * current cell and integrated over.
   *
   * @note This operation writes the data to the same field as
   * submit_gradient(). As a consequence, only one of these functions can be
   * used. In case several terms of this kind appear in a weak form, the
   * contribution of a potential call to this function must be added into the
   * diagonal of the rank-2 tensor contribution passed to submit_gradient(),
   * in order not to overwrite information.
   */
  void
  submit_divergence(const Number &value, const unsigned int point_index);

  /**
   * Return the curl in real coordinates at the point with index
   * `point_index` after a call to FEPointEvaluation::evaluate() with
   * EvaluationFlags::gradients set. This functions only makes sense for a
   * vector field with dim components and dim > 1.
   */
  Tensor<1, (dim == 2 ? 1 : dim), Number>
  get_curl(const unsigned int point_index) const;

  /**
   * Return the Jacobian of the transformation on the current cell with the
   * given point index. Prerequisite: This class needs to be constructed with
   * UpdateFlags containing `update_jacobian`.
   */
  DerivativeForm<1, dim, spacedim, Number>
  jacobian(const unsigned int point_index) const;

  /**
   * Return the inverse of the Jacobian of the transformation on the current
   * cell with the given point index. Prerequisite: This class needs to be
   * constructed with UpdateFlags containing `update_inverse_jacobian` or
   * `update_gradients`.
   */
  DerivativeForm<1, spacedim, dim, Number>
  inverse_jacobian(const unsigned int point_index) const;

  /**
   * Return the Jacobian determinant multiplied by the quadrature weight. This
   * class or the MappingInfo object passed to this function needs to be
   * constructed with UpdateFlags containing `update_JxW_values`.
   */
  Number
  JxW(const unsigned int point_index) const;

  /**
   * Return the position in real coordinates of the given point index among
   * the points passed to reinit().
   *
   * @deprecated Use the function quadrature_point() instead.
   */
  DEAL_II_DEPRECATED Point<spacedim, Number>
                     real_point(const unsigned int point_index) const;

  /**
   * Return the position in real coordinates of the given point index among
   * the points passed to reinit().
   */
  Point<spacedim, Number>
  quadrature_point(const unsigned int point_index) const;

  /**
   * Return the position in unit/reference coordinates of the given point
   * index, i.e., the respective point passed to the reinit() function.
   */
  Point<dim, Number>
  unit_point(const unsigned int point_index) const;

  /**
￼   * Take values collected at quadrature points via the submit_value()
    * function, multiply by the Jacobian determinant
    * and quadrature weights (JxW) and sum the values for all quadrature
￼   * points on the cell. The result is a scalar, representing the integral
￼   * of the function over the cell.
￼   */
  scalar_value_type
  integrate_value() const;

  /**
   * Return an object that can be thought of as an array containing all indices
   * from zero to n_quadrature_points. This allows to write code using
   * range-based for loops.
   */
  inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  quadrature_point_indices() const;

  /**
   * Returns how many lanes of a quadrature batch are active.
   */
  unsigned int
  n_active_entries_per_quadrature_batch(unsigned int q);

protected:
  static constexpr std::size_t n_lanes_user_interface =
    internal::VectorizedArrayTrait<Number>::width();
  static constexpr std::size_t n_lanes_internal =
    internal::VectorizedArrayTrait<VectorizedArrayType>::width();
  static constexpr std::size_t stride =
    internal::VectorizedArrayTrait<Number>::stride();

  /**
   * Common setup function for both constructors. Does the setup for both fast
   * and slow path.
   *
   * @param first_selected_component For multi-component FiniteElement
   * objects, this parameter allows to select a range of `n_components`
   * components starting from this parameter.
   */
  void
  setup(const unsigned int first_selected_component);

  /**
   * Shared functionality of all @p reinit() functions. Resizes data fields and
   * precomputes the @p shapes vector, holding the evaluation of 1D basis
   * functions of tensor product polynomials, if necessary.
   */
  template <bool is_face, bool is_linear>
  void
  do_reinit();

  /**
   * Number of quadrature batches of the current cell/face.
   */
  const unsigned int n_q_batches;

  /**
   * Number of quadrature points/batches of the current cell/face.
   */
  const unsigned int n_q_points;

  /**
   * Number of quadrature points of the current cell/face.
   */
  const unsigned int n_q_points_scalar;

  /**
   * Pointer to the Mapping object passed to the constructor.
   */
  ObserverPointer<const Mapping<dim, spacedim>> mapping;

  /**
   * Pointer to the FiniteElement object passed to the constructor.
   */
  ObserverPointer<const FiniteElement<dim, spacedim>> fe;

  /**
   * Description of the 1d polynomial basis for tensor product elements used
   * for the fast path of this class using tensor product evaluators.
   */
  std::vector<Polynomials::Polynomial<double>> poly;

  /**
   * Store whether the linear path should be used.
   */
  bool use_linear_path;

  /**
   * Renumbering between the unknowns of unknowns implied by the FiniteElement
   * class and a lexicographic numbering used for the tensorized code path.
   */
  std::vector<unsigned int> renumber;

  /**
   * Temporary array to store the `solution_values` passed to the evaluate()
   * function in a format compatible with the tensor product evaluators. For
   * vector-valued setups, this array uses a `Tensor<1, n_components,
   * ScalarNumber>` type to collect the unknowns for a particular basis
   * function.
   */
  std::vector<scalar_value_type> solution_renumbered;

  /**
   * Temporary array to store a vectorized version of the `solution_values`
   * computed during `integrate()` in a format compatible with the tensor
   * product evaluators. For vector-valued setups, this array uses a
   * `Tensor<1, n_components, VectorizedArrayType>` format.
   */
  AlignedVector<vectorized_value_type> solution_renumbered_vectorized;

  /**
   * Temporary array for the face path (scalar).
   */
  AlignedVector<ScalarNumber> scratch_data_scalar;

  /**
   * Temporary array to store the values at the points.
   */
  std::vector<value_type> values;

  /**
   * Temporary array to store the gradients in real coordinates at the points.
   */
  std::vector<gradient_type> gradients;

  /**
   * Pointer to first unit point batch of current cell/face from MappingInfo,
   * set internally during do_reinit().
   */
  const Point<dim, VectorizedArrayType> *unit_point_ptr;

  /**
   * Pointer to first unit point batch of current face from MappingInfo,
   * set internally during do_reinit(). Needed for face path.
   */
  const Point<dim - 1, VectorizedArrayType> *unit_point_faces_ptr;

  /**
   * Pointer to real point of first quadrature point of current cell/face from
   * MappingInfo, set internally during do_reinit().
   */
  const Point<spacedim, Number> *real_point_ptr;

  /**
   * Pointer to Jacobian of first quadrature point of current cell/face from
   * MappingInfo, set internally during do_reinit().
   */
  const DerivativeForm<1, dim, spacedim, Number> *jacobian_ptr;

  /**
   * Pointer to inverse Jacobian of first quadrature point of current cell/face
   * from MappingInfo, set internally during do_reinit().
   */
  const DerivativeForm<1, spacedim, dim, Number> *inverse_jacobian_ptr;

  /**
   * Pointer to normal vector of first quadrature point of current cell/face
   * from MappingInfo, set internally during do_reinit().
   */
  const Tensor<1, spacedim, Number> *normal_ptr;

  /**
   * Pointer to Jacobian determinant times quadrature weight of first quadrature
   * point of current cell/face from MappingInfo, set internally during
   * do_reinit().
   */
  const Number *JxW_ptr;

  /**
   * Cell type describing the geometry of the cell and compression of jacobians.
   */
  internal::MatrixFreeFunctions::GeometryType cell_type;

  /**
   * Number of unknowns per component, i.e., number of unique basis functions,
   * for the chosen FiniteElement (or base element).
   */
  unsigned int dofs_per_component;

  /**
   * Number of unknowns per component, i.e., number of unique basis functions,
   * for a restriction to the face of the chosen FiniteElement (or base
   * element). This means a (dim-1)-dimensional basis.
   */
  unsigned int dofs_per_component_face;

  /**
   * Scalar ShapeInfo object needed for face path.
   */
  internal::MatrixFreeFunctions::ShapeInfo<ScalarNumber> shape_info;

  /**
   * The first selected component in the active base element.
   */
  unsigned int component_in_base_element;

  /**
   * For complicated FiniteElement objects this variable informs us about
   * which unknowns actually carry degrees of freedom in the selected
   * components.
   */
  std::vector<std::array<bool, n_components>> nonzero_shape_function_component;

  /**
   * The desired update flags for the evaluation.
   */
  const UpdateFlags update_flags;

  /**
   * The FEValues object underlying the slow evaluation path.
   */
  std::shared_ptr<FEValues<dim, spacedim>> fe_values;

  /**
   * Pointer to mapping info on the fly computed during reinit.
   */
  std::unique_ptr<NonMatching::MappingInfo<dim, spacedim, Number>>
    mapping_info_on_the_fly;

  /**
   * Pointer to currently used mapping info (either on the fly or external
   * precomputed).
   */
  ObserverPointer<const NonMatching::MappingInfo<dim, spacedim, Number>>
    mapping_info;

  /**
   * The current cell index to access mapping data from mapping info.
   */
  unsigned int current_cell_index;

  /**
   * The current face number to access mapping data from mapping info.
   */
  unsigned int current_face_number;

  /**
   * Bool indicating if fast path is chosen.
   */
  bool fast_path;

  /**
   * Bool indicating if class needs to call reinit() inside
   * evaluate()/integrate() functions, which is the case when the present
   * class does not own the MappingInfo object but shares evaluation points
   * with another object.
   */
  bool must_reinitialize_pointers;

  /**
   * Vector containing tensor product shape functions evaluated (during
   * reinit()) at the vectorized unit points.
   */
  AlignedVector<dealii::ndarray<VectorizedArrayType, 2, dim>> shapes;

  /**
   * Vector containing tensor product shape functions evaluated (during
   * reinit()) at the vectorized unit points on faces.
   */
  AlignedVector<dealii::ndarray<VectorizedArrayType, 2, dim - 1>> shapes_faces;

  const bool is_interior;
};



/**
 * This class provides an interface to the evaluation of interpolated solution
 * values and gradients on cells on arbitrary reference point positions. These
 * points can change from cell to cell, both with respect to their quantity as
 * well to the location. The two typical use cases are evaluations on
 * non-matching grids and particle simulations.
 *
 * The use of this class is similar to FEValues or FEEvaluation: The class is
 * first initialized to a cell by calling `FEPointEvaluation::reinit(cell,
 * unit_points)`, with the main difference to the other concepts that the
 * underlying points in reference coordinates need to be passed along. Then,
 * upon call to evaluate() or integrate(), the user can compute information at
 * the give points. Eventually, the access functions get_value() or
 * get_gradient() allow to query this information at a specific point index.
 *
 * The functionality is similar to creating an FEValues object with a
 * Quadrature object on the `unit_points` on every cell separately and then
 * calling FEValues::get_function_values or FEValues::get_function_gradients,
 * and for some elements and mappings this is what actually happens
 * internally. For specific combinations of Mapping and FiniteElement
 * realizations, however, there is a much more efficient implementation that
 * avoids the memory allocation and other expensive start-up cost of
 * FEValues. Currently, the functionality is specialized for mappings derived
 * from MappingQ and MappingCartesian and for finite elements with tensor
 * product structure that work with the
 * @ref matrixfree
 * topic. In those cases, the cost implied
 * by this class is similar (or sometimes even somewhat lower) than using
 * `FEValues::reinit(cell)` followed by `FEValues::get_function_gradients`.
 */
template <int n_components_,
          int dim,
          int spacedim    = dim,
          typename Number = double>
class FEPointEvaluation
  : public FEPointEvaluationBase<n_components_, dim, spacedim, Number>
{
public:
  static constexpr unsigned int dimension    = dim;
  static constexpr unsigned int n_components = n_components_;

  using number_type = Number;

  using ScalarNumber =
    typename internal::VectorizedArrayTrait<Number>::value_type;
  using VectorizedArrayType = typename dealii::internal::VectorizedArrayTrait<
    Number>::vectorized_value_type;
  using ETT = typename internal::FEPointEvaluation::
    EvaluatorTypeTraits<dim, spacedim, n_components, Number>;
  using value_type            = typename ETT::value_type;
  using scalar_value_type     = typename ETT::scalar_value_type;
  using vectorized_value_type = typename ETT::vectorized_value_type;
  using unit_gradient_type    = typename ETT::unit_gradient_type;
  using gradient_type         = typename ETT::real_gradient_type;
  using interface_vectorized_unit_gradient_type =
    typename ETT::interface_vectorized_unit_gradient_type;

  /**
   * Constructor.
   *
   * @param mapping The Mapping class describing the actual geometry of a cell
   * passed to the evaluate() function.
   *
   * @param fe The FiniteElement object that is used for the evaluation, which
   * is typically the same on all cells to be evaluated.
   *
   * @param update_flags Specify the quantities to be computed by the mapping
   * during the call of reinit(). During evaluate() or integrate(), this data
   * is queried to produce the desired result (e.g., the gradient of a finite
   * element solution).
   *
   * @param first_selected_component For multi-component FiniteElement
   * objects, this parameter allows to select a range of `n_components`
   * components starting from this parameter.
   *
   * @param force_lexicographic_numbering By default, this class uses the DoF
   * numbering detected from the FiniteElement passed to this class (which is
   * lexicographic for e.g. FE_DGQ but hierarchic for FE_Q). However, it is
   * possible to force a lexicographic numbering of DoFs instead by setting this
   * parameter to true, which is useful if this class is used together with
   * FEEvaluation, for instance. This is because DoF values within FEEvaluation
   * are stored in lexicographic numbering.
   */
  FEPointEvaluation(const Mapping<dim, spacedim>       &mapping,
                    const FiniteElement<dim, spacedim> &fe,
                    const UpdateFlags                   update_flags,
                    const unsigned int first_selected_component      = 0,
                    const bool         force_lexicographic_numbering = false);

  /**
   * Constructor to make the present class able to re-use the geometry
   * data also used by other `FEPointEvaluation` objects.
   *
   * @param mapping_info The MappingInfo class describes the geometry-related
   * data for evaluating finite-element solutions. This object enables to
   * construct such an object on the outside, possibly re-using it between
   * several objects or between several calls to the same cell and unit points.
   *
   * @param fe The FiniteElement object that is used for the evaluation, which
   * is typically the same on all cells to be evaluated.
   *
   * @param first_selected_component For multi-component FiniteElement
   * objects, this parameter allows to select a range of `n_components`
   * components starting from this parameter.
   *
   * @param force_lexicographic_numbering By default, this class uses the DoF
   * numbering detected from the FiniteElement passed to this class (which is
   * lexicographic for e.g. FE_DGQ but hierarchic for FE_Q). However, it is
   * possible to force a lexicographic numbering of DoFs instead by setting this
   * parameter to true, which is useful if this class is used together with
   * FEEvaluation, for instance. This is because DoF values within FEEvaluation
   * are stored in lexicographic numbering.
   */
  FEPointEvaluation(
    const NonMatching::MappingInfo<dim, spacedim, Number> &mapping_info,
    const FiniteElement<dim, spacedim>                    &fe,
    const unsigned int first_selected_component      = 0,
    const bool         force_lexicographic_numbering = false);

  /**
   * Set up the mapping information for the given cell, e.g., by computing the
   * Jacobian of the mapping for the given points if gradients of the functions
   * are requested.
   *
   * @param[in] cell An iterator to the current cell
   *
   * @param[in] unit_points List of points in the reference locations of the
   * current cell where the FiniteElement object should be
   * evaluated/integrated in the evaluate() and integrate() functions.
   */
  void
  reinit(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
         const ArrayView<const Point<dim>> &unit_points);

  /**
   * Reinitialize the evaluator to point to the correct precomputed mapping of
   * the single cell in the MappingInfo object.
   */
  void
  reinit();

  /**
   * Reinitialize the evaluator to point to the correct precomputed mapping of
   * the cell in the MappingInfo object.
   */
  void
  reinit(const unsigned int cell_index);


  /**
   * This function interpolates the finite element solution, represented by
   * `solution_values`, on the cell and `unit_points` passed to reinit().
   *
   * @param[in] solution_values This array is supposed to contain the unknown
   * values on the element read out by
   * `FEEvaluation::read_dof_values(global_vector)`.
   *
   * @param[in] evaluation_flags Flags specifying which quantities should be
   * evaluated at the points.
   */
  template <std::size_t stride_view>
  void
  evaluate(
    const StridedArrayView<const ScalarNumber, stride_view> &solution_values,
    const EvaluationFlags::EvaluationFlags                  &evaluation_flags);

  /**
   * This function interpolates the finite element solution, represented by
   * `solution_values`, on the cell and `unit_points` passed to reinit().
   *
   * @param[in] solution_values This array is supposed to contain the unknown
   * values on the element as returned by `cell->get_dof_values(global_vector,
   * solution_values)`.
   *
   * @param[in] evaluation_flags Flags specifying which quantities should be
   * evaluated at the points.
   */
  void
  evaluate(const ArrayView<const ScalarNumber>    &solution_values,
           const EvaluationFlags::EvaluationFlags &evaluation_flags);

  /**
   * This function multiplies the quantities passed in by previous
   * submit_value() or submit_gradient() calls by the value or gradient of the
   * test functions, and performs summation over all given points multiplied be
   * the Jacobian determinant times the quadrature weight (JxW).
   *
   * @param[out] solution_values This array will contain the result of the
   * integral, which can be used during
   * `FEEvaluation::set_dof_values(global_vector)` or
   * `FEEvaluation::distribute_local_to_global(global_vector)`. Note
   * that for multi-component systems where only some of the components are
   * selected by the present class, the entries in `solution_values` not touched
   * by this class will be set to zero.
   *
   * @param[in] integration_flags Flags specifying which quantities should be
   * integrated at the points.
   *
   * @param[in] sum_into_values Flag specifying if the integrated values
   * should be summed into the solution values. For the default value
   * `sum_into_values=false` every value of @p solution_values is zeroed out.
   *
   */
  template <std::size_t stride_view>
  void
  integrate(const StridedArrayView<ScalarNumber, stride_view> &solution_values,
            const EvaluationFlags::EvaluationFlags &integration_flags,
            const bool                              sum_into_values = false);

  /**
   * This function multiplies the quantities passed in by previous
   * submit_value() or submit_gradient() calls by the value or gradient of the
   * test functions, and performs summation over all given points multiplied be
   * the Jacobian determinant times the quadrature weight (JxW).
   *
   * @param[out] solution_values This array will contain the result of the
   * integral, which can be used to during
   * `cell->set_dof_values(solution_values, global_vector)` or
   * `cell->distribute_local_to_global(solution_values, global_vector)`. Note
   * that for multi-component systems where only some of the components are
   * selected by the present class, the entries in `solution_values` not touched
   * by this class will be set to zero.
   *
   * @param[in] integration_flags Flags specifying which quantities should be
   * integrated at the points.
   *
   * @param[in] sum_into_values Flag specifying if the integrated values
   * should be summed into the solution values. For the default value
   * `sum_into_values=false` every value of @p solution_values is zeroed out.
   *
   */
  void
  integrate(const ArrayView<ScalarNumber>          &solution_values,
            const EvaluationFlags::EvaluationFlags &integration_flags,
            const bool                              sum_into_values = false);

  /**
   * This function multiplies the quantities passed in by previous
   * submit_value() or submit_gradient() calls by the value or gradient of the
   * test functions, and performs summation over all given points. This is
   * similar to the integration of a bilinear form in terms of the test
   * function, with the difference that this formula does not include a `JxW`
   * factor (in contrast to the integrate function of this class). This allows
   * the class to naturally embed point information (e.g. particles) into a
   * finite element formulation.
   *
   * @param[out] solution_values This array will contain the result of the
   * integral, which can be used during
   * `FEEvaluation::set_dof_values(global_vector)` or
   * `FEEvaluation::distribute_local_to_global(global_vector)`. Note
   * that for multi-component systems where only some of the components are
   * selected by the present class, the entries in `solution_values` not touched
   * by this class will be set to zero.
   *
   * @param[in] integration_flags Flags specifying which quantities should be
   * integrated at the points.
   *
   * @param[in] sum_into_values Flag specifying if the integrated values
   * should be summed into the solution values. For the default value
   * `sum_into_values=false` every value of @p solution_values is zeroed out.
   *
   */
  template <std::size_t stride_view>
  void
  test_and_sum(
    const StridedArrayView<ScalarNumber, stride_view> &solution_values,
    const EvaluationFlags::EvaluationFlags            &integration_flags,
    const bool                                         sum_into_values = false);

  /**
   * This function multiplies the quantities passed in by previous
   * submit_value() or submit_gradient() calls by the value or gradient of the
   * test functions, and performs summation over all given points. This is
   * similar to the integration of a bilinear form in terms of the test
   * function, with the difference that this formula does not include a `JxW`
   * factor (in contrast to the integrate function of this class). This allows
   * the class to naturally embed point information (e.g. particles) into a
   * finite element formulation.
   *
   * @param[out] solution_values This array will contain the result of the
   * integral, which can be used during
   * `cell->set_dof_values(solution_values, global_vector)` or
   * `cell->distribute_local_to_global(solution_values, global_vector)`. Note
   * that for multi-component systems where only some of the components are
   * selected by the present class, the entries in `solution_values` not touched
   * by this class will be set to zero.
   *
   * @param[in] integration_flags Flags specifying which quantities should be
   * integrated at the points.
   *
   * @param[in] sum_into_values Flag specifying if the integrated values
   * should be summed into the solution values. For the default value
   * `sum_into_values=false` every value of @p solution_values is zeroed out.
   *
   */
  void
  test_and_sum(const ArrayView<ScalarNumber>          &solution_values,
               const EvaluationFlags::EvaluationFlags &integration_flags,
               const bool                              sum_into_values = false);

  /**
   * Return the normal vector. This class or the MappingInfo object passed to
   * this function needs to be constructed with UpdateFlags containing
   * `update_normal_vectors`.
   */
  Tensor<1, spacedim, Number>
  normal_vector(const unsigned int point_index) const;

  /**
   * Return the normal derivative in real coordinates at the point with index
   * `point_index` after a call to FEPointEvaluation::evaluate() with
   * EvaluationFlags::gradients set.
   */
  const value_type
  get_normal_derivative(const unsigned int point_index) const;

  /**
   * Write a contribution that is tested by the normal derivative to the field
   * containing the values on points with the given `point_index`. Access to
   * the same field as through set_gradient()/get_gradient.
   */
  void
  submit_normal_derivative(const value_type &, const unsigned int point_index);

private:
  static constexpr std::size_t n_lanes_user_interface =
    internal::VectorizedArrayTrait<Number>::width();
  static constexpr std::size_t n_lanes_internal =
    internal::VectorizedArrayTrait<VectorizedArrayType>::width();
  static constexpr std::size_t stride =
    internal::VectorizedArrayTrait<Number>::stride();

  const bool lexicographic_numbering;

  /**
   * Resizes necessary data fields, reads in and renumbers solution values.
   * Interpolates onto face if face path is selected.
   */
  template <bool is_linear, std::size_t stride_view>
  void
  prepare_evaluate_fast(
    const StridedArrayView<const ScalarNumber, stride_view> &solution_values);

  /**
   * Evaluates the actual interpolation on the cell or face for a quadrature
   * batch.
   */
  template <bool is_linear, std::size_t stride_view>
  void
  compute_evaluate_fast(
    const StridedArrayView<const ScalarNumber, stride_view> &solution_values,
    const EvaluationFlags::EvaluationFlags                  &evaluation_flags,
    const unsigned int                                       n_shapes,
    const unsigned int                                       qb,
    vectorized_value_type                                   &value,
    interface_vectorized_unit_gradient_type                 &gradient);

  /**
   * Fast path of the evaluate function.
   */
  template <bool is_linear, std::size_t stride_view>
  void
  evaluate_fast(
    const StridedArrayView<const ScalarNumber, stride_view> &solution_values,
    const EvaluationFlags::EvaluationFlags                  &evaluation_flags);

  /**
   * Slow path of the evaluate function using FEValues.
   */
  template <std::size_t stride_view>
  void
  evaluate_slow(
    const StridedArrayView<const ScalarNumber, stride_view> &solution_values,
    const EvaluationFlags::EvaluationFlags                  &evaluation_flags);

  /**
   * Integrates the product of the data passed in by submit_value() and
   * submit_gradient() with the values or gradients of test functions on the
   * cell or face for a given quadrature batch.
   */
  template <bool is_linear>
  void
  compute_integrate_fast(
    const EvaluationFlags::EvaluationFlags       &integration_flags,
    const unsigned int                            n_shapes,
    const unsigned int                            qb,
    const vectorized_value_type                   value,
    const interface_vectorized_unit_gradient_type gradient,
    vectorized_value_type *solution_values_vectorized_linear);

  /**
   * Addition across the lanes of VectorizedArray as accumulated by the
   * compute_integrate_fast_function(), writing the sum into the result vector.
   * Applies face contributions to cell contributions for face path.
   */
  template <bool is_linear, std::size_t stride_view>
  void
  finish_integrate_fast(
    const StridedArrayView<ScalarNumber, stride_view> &solution_values,
    vectorized_value_type *solution_values_vectorized_linear,
    const bool             sum_into_values);

  /**
   * Fast path of the integrate function.
   */
  template <bool do_JxW, bool is_linear, std::size_t stride_view>
  void
  integrate_fast(
    const StridedArrayView<ScalarNumber, stride_view> &solution_values,
    const EvaluationFlags::EvaluationFlags            &integration_flags,
    const bool                                         sum_into_values);

  /**
   * Slow path of the integrate function using FEValues.
   */
  template <bool do_JxW, std::size_t stride_view>
  void
  integrate_slow(
    const StridedArrayView<ScalarNumber, stride_view> &solution_values,
    const EvaluationFlags::EvaluationFlags            &integration_flags,
    const bool                                         sum_into_values);

  /**
   * Implementation of the integrate/test_and_sum function.
   */
  template <bool do_JxW, std::size_t stride_view>
  void
  do_integrate(
    const StridedArrayView<ScalarNumber, stride_view> &solution_values,
    const EvaluationFlags::EvaluationFlags            &integration_flags,
    const bool                                         sum_into_values);

  /**
   * Internal function to initialize the pointers of this class when an
   * external MappingInfo has already queried the mapping for the relevant
   * information.
   */
  void
  internal_reinit_single_cell_state_mapping_info();
};



/**
 * This class provides an interface to the evaluation of interpolated solution
 * values and gradients on faces on arbitrary reference point positions. These
 * points can change from face to face, both with respect to their quantity as
 * well to the location. A typical use case is evaluations on non-matching
 * grids.
 *
 * The use of this class is similar to FEEvaluation: In the constructor, a
 * reference to a NonMatching::MappingInfo object is passed, where the
 * quadrature points in reference position is stored together with the mapping
 * information. The class is then reinitialized to a cell by calling
 * `FEFacePointEvaluation::reinit(face_index)` or
 * `FEFacePointEvaluation::reinit(cell_index, face_number)`. Then, upon call to
 * evaluate() or integrate(), the user can compute information at the given
 * points. Eventually, the access functions get_value() or get_gradient() allow
 * to query this information at a specific point index.
 */
template <int n_components_,
          int dim,
          int spacedim    = dim,
          typename Number = double>
class FEFacePointEvaluation
  : public FEPointEvaluationBase<n_components_, dim, spacedim, Number>
{
public:
  static constexpr unsigned int dimension    = dim;
  static constexpr unsigned int n_components = n_components_;

  using number_type = Number;

  using ScalarNumber =
    typename internal::VectorizedArrayTrait<Number>::value_type;
  using VectorizedArrayType = typename dealii::internal::VectorizedArrayTrait<
    Number>::vectorized_value_type;
  using ETT = typename internal::FEPointEvaluation::
    EvaluatorTypeTraits<dim, spacedim, n_components, Number>;
  using value_type            = typename ETT::value_type;
  using scalar_value_type     = typename ETT::scalar_value_type;
  using vectorized_value_type = typename ETT::vectorized_value_type;
  using unit_gradient_type    = typename ETT::unit_gradient_type;
  using gradient_type         = typename ETT::real_gradient_type;
  using interface_vectorized_unit_gradient_type =
    typename ETT::interface_vectorized_unit_gradient_type;

  /**
   * Constructor. Allows to select if interior or exterior face is selected.
   */
  FEFacePointEvaluation(
    const NonMatching::MappingInfo<dim, spacedim, Number> &mapping_info,
    const FiniteElement<dim, spacedim>                    &fe,
    const bool                                             is_interior = true,
    const unsigned int first_selected_component                        = 0);

  /**
   * Reinitialize the evaluator to point to the correct precomputed mapping of
   * the face in the MappingInfo object. Used in element-centric loops (ECL).
   */
  void
  reinit(const unsigned int cell_index, const unsigned int face_number);

  /**
   * Reinitialize the evaluator to point to the correct precomputed mapping of
   * the face in the MappingInfo object. Used in face-centric loops (FCL).
   */
  void
  reinit(const unsigned int face_index);

  /**
   * This function interpolates the finite element solution, represented by
   * `solution_values`, on the cell and `unit_points` passed to reinit().
   *
   * @param[in] solution_values This array is supposed to contain the unknown
   * values on the element read out by
   * `FEEvaluation::read_dof_values(global_vector)`.
   *
   * @param[in] evaluation_flags Flags specifying which quantities should be
   * evaluated at the points.
   */
  template <std::size_t stride_view>
  void
  evaluate(
    const StridedArrayView<const ScalarNumber, stride_view> &solution_values,
    const EvaluationFlags::EvaluationFlags                  &evaluation_flags);

  /**
   * This function interpolates the finite element solution, represented by
   * `solution_values`, on the cell and `unit_points` passed to reinit().
   *
   * @param[in] solution_values This array is supposed to contain the unknown
   * values on the element as returned by `cell->get_dof_values(global_vector,
   * solution_values)`.
   *
   * @param[in] evaluation_flags Flags specifying which quantities should be
   * evaluated at the points.
   */
  void
  evaluate(const ArrayView<const ScalarNumber>    &solution_values,
           const EvaluationFlags::EvaluationFlags &evaluation_flags);

  /**
   * This function multiplies the quantities passed in by previous
   * submit_value() or submit_gradient() calls by the value or gradient of the
   * test functions, and performs summation over all given points multiplied be
   * the Jacobian determinant times the quadrature weight (JxW).
   *
   * @param[out] solution_values This array will contain the result of the
   * integral, which can be used during
   * `FEEvaluation::set_dof_values(global_vector)` or
   * `FEEvaluation::distribute_local_to_global(global_vector)`. Note
   * that for multi-component systems where only some of the components are
   * selected by the present class, the entries in `solution_values` not touched
   * by this class will be set to zero.
   *
   * @param[in] integration_flags Flags specifying which quantities should be
   * integrated at the points.
   *
   * @param[in] sum_into_values Flag specifying if the integrated values
   * should be summed into the solution values. For the default value
   * `sum_into_values=false` every value of @p solution_values is zeroed out.
   *
   */
  template <std::size_t stride_view>
  void
  integrate(const StridedArrayView<ScalarNumber, stride_view> &solution_values,
            const EvaluationFlags::EvaluationFlags &integration_flags,
            const bool                              sum_into_values = false);

  /**
   * This function multiplies the quantities passed in by previous
   * submit_value() or submit_gradient() calls by the value or gradient of the
   * test functions, and performs summation over all given points multiplied be
   * the Jacobian determinant times the quadrature weight (JxW).
   *
   * @param[out] solution_values This array will contain the result of the
   * integral, which can be used to during
   * `cell->set_dof_values(solution_values, global_vector)` or
   * `cell->distribute_local_to_global(solution_values, global_vector)`. Note
   * that for multi-component systems where only some of the components are
   * selected by the present class, the entries in `solution_values` not touched
   * by this class will be set to zero.
   *
   * @param[in] integration_flags Flags specifying which quantities should be
   * integrated at the points.
   *
   * @param[in] sum_into_values Flag specifying if the integrated values
   * should be summed into the solution values. For the default value
   * `sum_into_values=false` every value of @p solution_values is zeroed out.
   *
   */
  void
  integrate(const ArrayView<ScalarNumber>          &solution_values,
            const EvaluationFlags::EvaluationFlags &integration_flags,
            const bool                              sum_into_values = false);

  /**
   * This function multiplies the quantities passed in by previous
   * submit_value() or submit_gradient() calls by the value or gradient of the
   * test functions, and performs summation over all given points multiplied be
   * the Jacobian determinant times the quadrature weight (JxW).
   *
   * @param[out] solution_values This array will contain the result of the
   * integral, which can be used during
   * `FEEvaluation::set_dof_values(global_vector)` or
   * `FEEvaluation::distribute_local_to_global(global_vector)`. Note
   * that for multi-component systems where only some of the components are
   * selected by the present class, the entries in `solution_values` not touched
   * by this class will be set to zero.
   *
   * @param[in] integration_flags Flags specifying which quantities should be
   * integrated at the points.
   *
   * @param[in] sum_into_values Flag specifying if the integrated values
   * should be summed into the solution values. For the default value
   * `sum_into_values=false` every value of @p solution_values is zeroed out.
   *
   */
  template <std::size_t stride_view>
  void
  test_and_sum(
    const StridedArrayView<ScalarNumber, stride_view> &solution_values,
    const EvaluationFlags::EvaluationFlags            &integration_flags,
    const bool                                         sum_into_values = false);

  /**
   * This function multiplies the quantities passed in by previous
   * submit_value() or submit_gradient() calls by the value or gradient of the
   * test functions, and performs summation over all given points multiplied be
   * the Jacobian determinant times the quadrature weight (JxW).
   *
   * @param[out] solution_values This array will contain the result of the
   * integral, which can be used to during
   * `cell->set_dof_values(solution_values, global_vector)` or
   * `cell->distribute_local_to_global(solution_values, global_vector)`. Note
   * that for multi-component systems where only some of the components are
   * selected by the present class, the entries in `solution_values` not touched
   * by this class will be set to zero.
   *
   * @param[in] integration_flags Flags specifying which quantities should be
   * integrated at the points.
   *
   * @param[in] sum_into_values Flag specifying if the integrated values
   * should be summed into the solution values. For the default value
   * `sum_into_values=false` every value of @p solution_values is zeroed out.
   *
   */
  void
  test_and_sum(const ArrayView<ScalarNumber>          &solution_values,
               const EvaluationFlags::EvaluationFlags &integration_flags,
               const bool                              sum_into_values = false);

  /**
   * Evaluate values and gradients in face for the selected face (lane) of the
   * batch. Default stride into the face dofs is width of
   * VectorizedArray<selected_floating_point_type> which is the default
   * vectorization over faces for FEFaceEvaluation.
   */
  template <int stride_face_dof = VectorizedArrayType::size()>
  void
  evaluate_in_face(const ScalarNumber                     *face_dof_values,
                   const EvaluationFlags::EvaluationFlags &evaluation_flags);

  /**
   * Integrate values and gradients in face for the selected face (lane) of the
   * batch. Default stride into the face dofs is width of
   * VectorizedArray<selected_floating_point_type> which is the default
   * vectorization over faces for FEFaceEvaluation.
   */
  template <int stride_face_dof = VectorizedArrayType::size()>
  void
  integrate_in_face(ScalarNumber                           *face_dof_values,
                    const EvaluationFlags::EvaluationFlags &integration_flags,
                    const bool sum_into_values = false);

  /**
   * Return the normal vector. This class or the MappingInfo object passed to
   * this function needs to be constructed with UpdateFlags containing
   * `update_normal_vectors`.
   */
  Tensor<1, spacedim, Number>
  normal_vector(const unsigned int point_index) const;

  /**
   * Return the normal derivative in real coordinates at the point with index
   * `point_index` after a call to FEFacePointEvaluation::evaluate() with
   * EvaluationFlags::gradients set.
   */
  const value_type
  get_normal_derivative(const unsigned int point_index) const;

  /**
   * Write a contribution that is tested by the normal derivative to the field
   * containing the values on points with the given `point_index`. Access to
   * the same field as through set_gradient()/get_gradient.
   */
  void
  submit_normal_derivative(const value_type &, const unsigned int point_index);

private:
  static constexpr std::size_t n_lanes_user_interface =
    internal::VectorizedArrayTrait<Number>::width();
  static constexpr std::size_t n_lanes_internal =
    internal::VectorizedArrayTrait<VectorizedArrayType>::width();
  static constexpr std::size_t stride =
    internal::VectorizedArrayTrait<Number>::stride();

  template <bool is_linear, std::size_t stride_view>
  void
  do_evaluate(
    const StridedArrayView<const ScalarNumber, stride_view> &solution_values,
    const EvaluationFlags::EvaluationFlags                  &evaluation_flags);

  template <bool do_JxW, bool is_linear, std::size_t stride_view>
  void
  do_integrate(
    const StridedArrayView<ScalarNumber, stride_view> &solution_values,
    const EvaluationFlags::EvaluationFlags            &integration_flags,
    const bool                                         sum_into_values);

  /**
   * Actually does the evaluation templated on the chosen code path (linear or
   * higher order).
   */
  template <bool is_linear, int stride_face_dof>
  void
  do_evaluate_in_face(const ScalarNumber                     *face_dof_values,
                      const EvaluationFlags::EvaluationFlags &evaluation_flags);

  /**
   * Actually does the integration templated on the chosen code path (linear or
   * higher order).
   */
  template <bool do_JxW, bool is_linear, int stride_face_dof>
  void
  do_integrate_in_face(
    ScalarNumber                           *face_dof_values,
    const EvaluationFlags::EvaluationFlags &integration_flags,
    const bool                              sum_into_values);
};



// ----------------------- template and inline function ----------------------


template <int n_components_, int dim, int spacedim, typename Number>
FEPointEvaluationBase<n_components_, dim, spacedim, Number>::
  FEPointEvaluationBase(const Mapping<dim, spacedim>       &mapping,
                        const FiniteElement<dim, spacedim> &fe,
                        const UpdateFlags                   update_flags,
                        const unsigned int first_selected_component)
  : n_q_batches(numbers::invalid_unsigned_int)
  , n_q_points(numbers::invalid_unsigned_int)
  , n_q_points_scalar(numbers::invalid_unsigned_int)
  , mapping(&mapping)
  , fe(&fe)
  , JxW_ptr(nullptr)
  , update_flags(update_flags)
  , mapping_info_on_the_fly(
      std::make_unique<NonMatching::MappingInfo<dim, spacedim, Number>>(
        mapping,
        update_flags))
  , mapping_info(mapping_info_on_the_fly.get())
  , current_cell_index(numbers::invalid_unsigned_int)
  , current_face_number(numbers::invalid_unsigned_int)
  , must_reinitialize_pointers(false)
  , is_interior(true)
{
  setup(first_selected_component);
}



template <int n_components_, int dim, int spacedim, typename Number>
FEPointEvaluationBase<n_components_, dim, spacedim, Number>::
  FEPointEvaluationBase(
    const NonMatching::MappingInfo<dim, spacedim, Number> &mapping_info,
    const FiniteElement<dim, spacedim>                    &fe,
    const unsigned int first_selected_component,
    const bool         is_interior)
  : n_q_batches(numbers::invalid_unsigned_int)
  , n_q_points(numbers::invalid_unsigned_int)
  , n_q_points_scalar(numbers::invalid_unsigned_int)
  , mapping(&mapping_info.get_mapping())
  , fe(&fe)
  , JxW_ptr(nullptr)
  , update_flags(mapping_info.get_update_flags())
  , mapping_info(&mapping_info)
  , current_cell_index(numbers::invalid_unsigned_int)
  , current_face_number(numbers::invalid_unsigned_int)
  , must_reinitialize_pointers(true)
  , is_interior(is_interior)
{
  setup(first_selected_component);
}



template <int n_components_, int dim, int spacedim, typename Number>
FEPointEvaluationBase<n_components_, dim, spacedim, Number>::
  FEPointEvaluationBase(
    FEPointEvaluationBase<n_components_, dim, spacedim, Number> &other) noexcept
  : n_q_batches(other.n_q_batches)
  , n_q_points(other.n_q_points)
  , n_q_points_scalar(other.n_q_points_scalar)
  , mapping(other.mapping)
  , fe(other.fe)
  , poly(other.poly)
  , use_linear_path(other.use_linear_path)
  , renumber(other.renumber)
  , solution_renumbered(other.solution_renumbered)
  , solution_renumbered_vectorized(other.solution_renumbered_vectorized)
  , values(other.values)
  , gradients(other.gradients)
  , dofs_per_component(other.dofs_per_component)
  , dofs_per_component_face(other.dofs_per_component_face)
  , component_in_base_element(other.component_in_base_element)
  , nonzero_shape_function_component(other.nonzero_shape_function_component)
  , update_flags(other.update_flags)
  , fe_values(other.fe_values)
  , mapping_info_on_the_fly(
      other.mapping_info_on_the_fly ?
        std::make_unique<NonMatching::MappingInfo<dim, spacedim, Number>>(
          *mapping,
          update_flags) :
        nullptr)
  , mapping_info(other.mapping_info)
  , current_cell_index(other.current_cell_index)
  , current_face_number(other.current_face_number)
  , fast_path(other.fast_path)
  , must_reinitialize_pointers(true)
  , is_interior(other.is_interior)
{}



template <int n_components_, int dim, int spacedim, typename Number>
FEPointEvaluationBase<n_components_, dim, spacedim, Number>::
  FEPointEvaluationBase(FEPointEvaluationBase &&other) noexcept = default;



template <int n_components_, int dim, int spacedim, typename Number>
void
FEPointEvaluationBase<n_components_, dim, spacedim, Number>::setup(
  const unsigned int first_selected_component)
{
  AssertIndexRange(first_selected_component + n_components,
                   fe->n_components() + 1);

  shapes.reserve(100);

  bool         same_base_element   = true;
  unsigned int base_element_number = 0;
  component_in_base_element        = 0;
  unsigned int component           = 0;
  for (; base_element_number < fe->n_base_elements(); ++base_element_number)
    if (component + fe->element_multiplicity(base_element_number) >
        first_selected_component)
      {
        if (first_selected_component + n_components >
            component + fe->element_multiplicity(base_element_number))
          same_base_element = false;
        component_in_base_element = first_selected_component - component;
        break;
      }
    else
      component += fe->element_multiplicity(base_element_number);

  if (internal::FEPointEvaluation::is_fast_path_supported(*mapping) &&
      internal::FEPointEvaluation::is_fast_path_supported(
        *fe, base_element_number) &&
      same_base_element)
    {
      shape_info.reinit(QMidpoint<1>(), *fe, base_element_number);
      renumber                = shape_info.lexicographic_numbering;
      dofs_per_component      = shape_info.dofs_per_component_on_cell;
      dofs_per_component_face = shape_info.dofs_per_component_on_face;
      poly = internal::FEPointEvaluation::get_polynomial_space(
        fe->base_element(base_element_number));

      bool is_lexicographic = true;
      for (unsigned int i = 0; i < renumber.size(); ++i)
        if (i != renumber[i])
          is_lexicographic = false;

      if (is_lexicographic)
        renumber.clear();

      use_linear_path = (poly.size() == 2 && poly[0].value(0.) == 1. &&
                         poly[0].value(1.) == 0. && poly[1].value(0.) == 0. &&
                         poly[1].value(1.) == 1.) &&
                        (fe->n_components() == n_components);

      const unsigned int size_face = 3 * dofs_per_component_face * n_components;
      const unsigned int size_cell = dofs_per_component * n_components;
      scratch_data_scalar.resize(size_face + size_cell);

      solution_renumbered.resize(dofs_per_component);
      solution_renumbered_vectorized.resize(dofs_per_component);

      fast_path = true;
    }
  else
    {
      nonzero_shape_function_component.resize(fe->n_dofs_per_cell());
      for (unsigned int d = 0; d < n_components; ++d)
        {
          const unsigned int component = first_selected_component + d;
          for (unsigned int i = 0; i < fe->n_dofs_per_cell(); ++i)
            {
              const bool is_primitive =
                fe->is_primitive() || fe->is_primitive(i);
              if (is_primitive)
                nonzero_shape_function_component[i][d] =
                  (component == fe->system_to_component_index(i).first);
              else
                nonzero_shape_function_component[i][d] =
                  (fe->get_nonzero_components(i)[component] == true);
            }
        }

      fast_path = false;
    }
}



template <int n_components_, int dim, int spacedim, typename Number>
template <bool is_face, bool is_linear>
inline void
FEPointEvaluationBase<n_components_, dim, spacedim, Number>::do_reinit()
{
  const unsigned int geometry_index =
    mapping_info->template compute_geometry_index_offset<is_face>(
      current_cell_index, current_face_number);

  cell_type = mapping_info->get_cell_type(geometry_index);

  const_cast<unsigned int &>(n_q_points_scalar) =
    mapping_info->get_n_q_points_unvectorized(geometry_index);

  // round up n_q_points_scalar / n_lanes_internal
  const_cast<unsigned int &>(n_q_batches) =
    (n_q_points_scalar + n_lanes_internal - 1) / n_lanes_internal;

  const unsigned int n_q_points_before = n_q_points;

  const_cast<unsigned int &>(n_q_points) =
    (stride == 1) ? n_q_batches : n_q_points_scalar;

  if (n_q_points != n_q_points_before)
    {
      if (update_flags & update_values)
        values.resize(n_q_points);
      if (update_flags & update_gradients)
        gradients.resize(n_q_points);
    }

  if (n_q_points == 0)
    return;

  // set unit point pointer
  const unsigned int unit_point_offset =
    mapping_info->compute_unit_point_index_offset(geometry_index);

  if (is_face)
    unit_point_faces_ptr =
      mapping_info->get_unit_point_faces(unit_point_offset);
  else
    unit_point_ptr = mapping_info->get_unit_point(unit_point_offset);

  // set data pointers
  const unsigned int data_offset =
    mapping_info->compute_data_index_offset(geometry_index);
  const unsigned int compressed_data_offset =
    mapping_info->compute_compressed_data_index_offset(geometry_index);
  if constexpr (running_in_debug_mode())
    {
      const UpdateFlags update_flags_mapping =
        mapping_info->get_update_flags_mapping();
      if (update_flags_mapping & UpdateFlags::update_quadrature_points)
        real_point_ptr = mapping_info->get_real_point(data_offset);
      if (update_flags_mapping & UpdateFlags::update_jacobians)
        jacobian_ptr =
          mapping_info->get_jacobian(compressed_data_offset, is_interior);
      if (update_flags_mapping & UpdateFlags::update_inverse_jacobians)
        inverse_jacobian_ptr =
          mapping_info->get_inverse_jacobian(compressed_data_offset,
                                             is_interior);
      if (update_flags_mapping & UpdateFlags::update_normal_vectors)
        normal_ptr = mapping_info->get_normal_vector(data_offset);
      if (update_flags_mapping & UpdateFlags::update_JxW_values)
        JxW_ptr = mapping_info->get_JxW(data_offset);
    }
  else
    {
      real_point_ptr = mapping_info->get_real_point(data_offset);
      jacobian_ptr =
        mapping_info->get_jacobian(compressed_data_offset, is_interior);
      inverse_jacobian_ptr =
        mapping_info->get_inverse_jacobian(compressed_data_offset, is_interior);
      normal_ptr = mapping_info->get_normal_vector(data_offset);
      JxW_ptr    = mapping_info->get_JxW(data_offset);
    }

  if (!is_linear && fast_path)
    {
      const std::size_t n_shapes = poly.size();
      if (is_face)
        shapes_faces.resize_fast(n_q_batches * n_shapes);
      else
        shapes.resize_fast(n_q_batches * n_shapes);

      for (unsigned int qb = 0; qb < n_q_batches; ++qb)
        if (is_face)
          {
            if (dim > 1)
              {
                internal::compute_values_of_array(
                  shapes_faces.data() + qb * n_shapes,
                  poly,
                  unit_point_faces_ptr[qb],
                  update_flags & UpdateFlags::update_gradients ? 1 : 0);
              }
          }
        else
          {
            if (update_flags & UpdateFlags::update_gradients)
              {
                internal::compute_values_of_array(shapes.data() + qb * n_shapes,
                                                  poly,
                                                  unit_point_ptr[qb],
                                                  1);
              }
            else if (qb + 1 < n_q_batches)
              {
                // Use function with reduced overhead to compute for two
                // points at once
                internal::compute_values_of_array_in_pairs(
                  shapes.data() + qb * n_shapes,
                  poly,
                  unit_point_ptr[qb],
                  unit_point_ptr[qb + 1]);
                ++qb;
              }
            else
              {
                internal::compute_values_of_array(shapes.data() + qb * n_shapes,
                                                  poly,
                                                  unit_point_ptr[qb],
                                                  0);
              }
          }
    }
}



template <int n_components_, int dim, int spacedim, typename Number>
inline const typename FEPointEvaluationBase<n_components_,
                                            dim,
                                            spacedim,
                                            Number>::value_type &
FEPointEvaluationBase<n_components_, dim, spacedim, Number>::get_value(
  const unsigned int point_index) const
{
  AssertIndexRange(point_index, values.size());
  return values[point_index];
}



template <int n_components_, int dim, int spacedim, typename Number>
inline const typename FEPointEvaluationBase<n_components_,
                                            dim,
                                            spacedim,
                                            Number>::gradient_type &
FEPointEvaluationBase<n_components_, dim, spacedim, Number>::get_gradient(
  const unsigned int point_index) const
{
  AssertIndexRange(point_index, gradients.size());
  return gradients[point_index];
}



template <int n_components_, int dim, int spacedim, typename Number>
inline Number
FEPointEvaluationBase<n_components_, dim, spacedim, Number>::get_divergence(
  const unsigned int point_index) const
{
  static_assert(n_components == dim,
                "Only makes sense for a vector field with dim components");

  AssertIndexRange(point_index, values.size());
  return trace(gradients[point_index]);
}



template <int n_components_, int dim, int spacedim, typename Number>
inline void
FEPointEvaluationBase<n_components_, dim, spacedim, Number>::submit_value(
  const value_type  &value,
  const unsigned int point_index)
{
  AssertIndexRange(point_index, n_q_points);
  values[point_index] = value;
}



template <int n_components_, int dim, int spacedim, typename Number>
inline void
FEPointEvaluationBase<n_components_, dim, spacedim, Number>::submit_gradient(
  const gradient_type &gradient,
  const unsigned int   point_index)
{
  AssertIndexRange(point_index, n_q_points);
  gradients[point_index] = gradient;
}



template <int n_components_, int dim, int spacedim, typename Number>
inline void
FEPointEvaluationBase<n_components_, dim, spacedim, Number>::submit_divergence(
  const Number      &value,
  const unsigned int point_index)
{
  static_assert(n_components == dim,
                "Only makes sense for a vector field with dim components");

  AssertIndexRange(point_index, n_q_points);
  gradients[point_index] = gradient_type();
  for (unsigned int d = 0; d < dim; ++d)
    gradients[point_index][d][d] = value;
}



template <int n_components_, int dim, int spacedim, typename Number>
Tensor<1, (dim == 2 ? 1 : dim), Number>
FEPointEvaluationBase<n_components_, dim, spacedim, Number>::get_curl(
  const unsigned int point_index) const
{
  static_assert(
    dim > 1 && n_components == dim,
    "Only makes sense for a vector field with dim components and dim > 1");

  const Tensor<2, dim, Number>            grad = get_gradient(point_index);
  Tensor<1, (dim == 2 ? 1 : dim), Number> curl;
  switch (dim)
    {
      case 2:
        curl[0] = grad[1][0] - grad[0][1];
        break;
      case 3:
        curl[0] = grad[2][1] - grad[1][2];
        curl[1] = grad[0][2] - grad[2][0];
        curl[2] = grad[1][0] - grad[0][1];
        break;
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
  return curl;
}



template <int n_components_, int dim, int spacedim, typename Number>
inline DerivativeForm<1, dim, spacedim, Number>
FEPointEvaluationBase<n_components_, dim, spacedim, Number>::jacobian(
  const unsigned int point_index) const
{
  AssertIndexRange(point_index, n_q_points);
  Assert(jacobian_ptr != nullptr,
         internal::FEPointEvaluation::
           ExcFEPointEvaluationAccessToUninitializedMappingField(
             "update_jacobians"));
  return jacobian_ptr[cell_type <= ::dealii::internal::MatrixFreeFunctions::
                                     GeometryType::affine ?
                        0 :
                        point_index];
}



template <int n_components_, int dim, int spacedim, typename Number>
inline DerivativeForm<1, spacedim, dim, Number>
FEPointEvaluationBase<n_components_, dim, spacedim, Number>::inverse_jacobian(
  const unsigned int point_index) const
{
  AssertIndexRange(point_index, n_q_points);
  Assert(inverse_jacobian_ptr != nullptr,
         internal::FEPointEvaluation::
           ExcFEPointEvaluationAccessToUninitializedMappingField(
             "update_inverse_jacobians"));
  return inverse_jacobian_ptr
    [cell_type <=
         ::dealii::internal::MatrixFreeFunctions::GeometryType::affine ?
       0 :
       point_index];
}



template <int n_components_, int dim, int spacedim, typename Number>
inline Number
FEPointEvaluationBase<n_components_, dim, spacedim, Number>::JxW(
  const unsigned int point_index) const
{
  AssertIndexRange(point_index, n_q_points);
  Assert(JxW_ptr != nullptr,
         internal::FEPointEvaluation::
           ExcFEPointEvaluationAccessToUninitializedMappingField(
             "update_JxW_values"));
  return JxW_ptr[point_index];
}



template <int n_components_, int dim, int spacedim, typename Number>
inline Point<spacedim, Number>
FEPointEvaluationBase<n_components_, dim, spacedim, Number>::real_point(
  const unsigned int point_index) const
{
  return quadrature_point(point_index);
}



template <int n_components_, int dim, int spacedim, typename Number>
inline Point<spacedim, Number>
FEPointEvaluationBase<n_components_, dim, spacedim, Number>::quadrature_point(
  const unsigned int point_index) const
{
  AssertIndexRange(point_index, n_q_points);
  Assert(real_point_ptr != nullptr,
         internal::FEPointEvaluation::
           ExcFEPointEvaluationAccessToUninitializedMappingField(
             "update_quadrature_points"));
  return real_point_ptr[point_index];
}



template <int n_components_, int dim, int spacedim, typename Number>
inline Point<dim, Number>
FEPointEvaluationBase<n_components_, dim, spacedim, Number>::unit_point(
  const unsigned int point_index) const
{
  AssertIndexRange(point_index, n_q_points);
  Assert(unit_point_ptr != nullptr, ExcMessage("unit_point_ptr is not set!"));
  Point<dim, Number> unit_point;
  for (unsigned int d = 0; d < dim; ++d)
    unit_point[d] = internal::VectorizedArrayTrait<Number>::get_from_vectorized(
      unit_point_ptr[point_index / stride][d], point_index % stride);
  return unit_point;
}



template <int n_components_, int dim, int spacedim, typename Number>
inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
FEPointEvaluationBase<n_components_, dim, spacedim, Number>::
  quadrature_point_indices() const
{
  return std_cxx20::ranges::iota_view<unsigned int, unsigned int>(0U,
                                                                  n_q_points);
}



template <int n_components_, int dim, int spacedim, typename Number>
FEPointEvaluation<n_components_, dim, spacedim, Number>::FEPointEvaluation(
  const NonMatching::MappingInfo<dim, spacedim, Number> &mapping_info,
  const FiniteElement<dim, spacedim>                    &fe,
  const unsigned int first_selected_component,
  const bool         force_lexicographic_numbering)
  : FEPointEvaluationBase<n_components_, dim, spacedim, Number>(
      mapping_info,
      fe,
      first_selected_component)
  , lexicographic_numbering(force_lexicographic_numbering ||
                            this->renumber.empty())
{}



template <int n_components_, int dim, int spacedim, typename Number>
FEPointEvaluation<n_components_, dim, spacedim, Number>::FEPointEvaluation(
  const Mapping<dim, spacedim>       &mapping,
  const FiniteElement<dim, spacedim> &fe,
  const UpdateFlags                   update_flags,
  const unsigned int                  first_selected_component,
  const bool                          force_lexicographic_numbering)
  : FEPointEvaluationBase<n_components_, dim, spacedim, Number>(
      mapping,
      fe,
      update_flags,
      first_selected_component)
  , lexicographic_numbering(force_lexicographic_numbering ||
                            this->renumber.empty())
{}



template <int n_components_, int dim, int spacedim, typename Number>
inline void
FEPointEvaluation<n_components_, dim, spacedim, Number>::
  internal_reinit_single_cell_state_mapping_info()
{
  this->current_cell_index  = numbers::invalid_unsigned_int;
  this->current_face_number = numbers::invalid_unsigned_int;

  if (this->use_linear_path)
    this->template do_reinit<false, true>();
  else
    this->template do_reinit<false, false>();
}



template <int n_components_, int dim, int spacedim, typename Number>
inline void
FEPointEvaluation<n_components_, dim, spacedim, Number>::reinit()
{
  internal_reinit_single_cell_state_mapping_info();
  this->must_reinitialize_pointers = false;
}



template <int n_components_, int dim, int spacedim, typename Number>
inline void
FEPointEvaluation<n_components_, dim, spacedim, Number>::reinit(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const ArrayView<const Point<dim>>                          &unit_points)
{
  // reinit is only allowed for mapping computation on the fly
  AssertThrow(this->mapping_info_on_the_fly.get() != nullptr,
              ExcNotImplemented());

  this->mapping_info_on_the_fly->reinit(cell, unit_points);
  this->must_reinitialize_pointers = false;

  if (!this->fast_path)
    {
      this->fe_values = std::make_shared<FEValues<dim, spacedim>>(
        *this->mapping,
        *this->fe,
        Quadrature<dim>(
          std::vector<Point<dim>>(unit_points.begin(), unit_points.end())),
        this->update_flags);
      this->fe_values->reinit(cell);
    }

  if (this->use_linear_path)
    this->template do_reinit<false, true>();
  else
    this->template do_reinit<false, false>();
}



template <int n_components_, int dim, int spacedim, typename Number>
inline void
FEPointEvaluation<n_components_, dim, spacedim, Number>::reinit(
  const unsigned int cell_index)
{
  this->current_cell_index         = cell_index;
  this->current_face_number        = numbers::invalid_unsigned_int;
  this->must_reinitialize_pointers = false;

  if (this->use_linear_path)
    this->template do_reinit<false, true>();
  else
    this->template do_reinit<false, false>();

  if (!this->fast_path)
    {
      std::vector<Point<dim>> unit_points(this->n_q_points_scalar);

      for (unsigned int v = 0; v < this->n_q_points_scalar; ++v)
        for (unsigned int d = 0; d < dim; ++d)
          unit_points[v][d] =
            this->unit_point_ptr[v / n_lanes_internal][d][v % n_lanes_internal];

      this->fe_values = std::make_shared<FEValues<dim, spacedim>>(
        *this->mapping,
        *this->fe,
        Quadrature<dim>(
          std::vector<Point<dim>>(unit_points.begin(), unit_points.end())),
        this->update_flags);

      this->fe_values->reinit(
        this->mapping_info->get_cell_iterator(this->current_cell_index));
    }
}



template <int n_components_, int dim, int spacedim, typename Number>
template <std::size_t stride_view>
void
FEPointEvaluation<n_components_, dim, spacedim, Number>::evaluate(
  const StridedArrayView<const ScalarNumber, stride_view> &solution_values,
  const EvaluationFlags::EvaluationFlags                  &evaluation_flags)
{
  Assert(!(evaluation_flags & EvaluationFlags::hessians), ExcNotImplemented());

  if (!((evaluation_flags & EvaluationFlags::values) ||
        (evaluation_flags & EvaluationFlags::gradients))) // no evaluation flags
    return;

  if (this->must_reinitialize_pointers)
    internal_reinit_single_cell_state_mapping_info();

  if (this->n_q_points == 0)
    return;

  AssertDimension(solution_values.size(), this->fe->dofs_per_cell);
  if (this->fast_path)
    {
      if (this->use_linear_path)
        evaluate_fast<true>(solution_values, evaluation_flags);
      else
        evaluate_fast<false>(solution_values, evaluation_flags);
    }
  else
    evaluate_slow(solution_values, evaluation_flags);
}



template <int n_components_, int dim, int spacedim, typename Number>
void
FEPointEvaluation<n_components_, dim, spacedim, Number>::evaluate(
  const ArrayView<const ScalarNumber>    &solution_values,
  const EvaluationFlags::EvaluationFlags &evaluation_flags)
{
  evaluate(StridedArrayView<const ScalarNumber, 1>(solution_values.data(),
                                                   solution_values.size()),
           evaluation_flags);
}



template <int n_components_, int dim, int spacedim, typename Number>
template <std::size_t stride_view>
void
FEPointEvaluation<n_components_, dim, spacedim, Number>::integrate(
  const StridedArrayView<ScalarNumber, stride_view> &solution_values,
  const EvaluationFlags::EvaluationFlags            &integration_flags,
  const bool                                         sum_into_values)
{
  do_integrate<true>(solution_values, integration_flags, sum_into_values);
}



template <int n_components_, int dim, int spacedim, typename Number>
void
FEPointEvaluation<n_components_, dim, spacedim, Number>::integrate(
  const ArrayView<ScalarNumber>          &solution_values,
  const EvaluationFlags::EvaluationFlags &integration_flags,
  const bool                              sum_into_values)
{
  integrate(StridedArrayView<ScalarNumber, 1>(solution_values.data(),
                                              solution_values.size()),
            integration_flags,
            sum_into_values);
}



template <int n_components_, int dim, int spacedim, typename Number>
inline typename FEPointEvaluationBase<n_components_, dim, spacedim, Number>::
  scalar_value_type
  FEPointEvaluationBase<n_components_, dim, spacedim, Number>::integrate_value()
    const
{
  value_type return_value = {};

  for (const auto point_index : this->quadrature_point_indices())
    return_value += values[point_index] * this->JxW(point_index);

  return ETT::sum_value(return_value);
}



template <int n_components_, int dim, int spacedim, typename Number>
unsigned int
FEPointEvaluationBase<n_components_, dim, spacedim, Number>::
  n_active_entries_per_quadrature_batch(unsigned int q)
{
  Assert(stride == 1,
         ExcMessage(
           "Calling this function only makes sense in fully vectorized mode."));
  if (q == n_q_batches - 1)
    {
      const unsigned int n_filled_lanes =
        n_q_points_scalar & (n_lanes_user_interface - 1);
      if (n_filled_lanes == 0)
        return n_lanes_user_interface;
      else
        return n_filled_lanes;
    }
  else
    return n_lanes_user_interface;
}



template <int n_components_, int dim, int spacedim, typename Number>
template <std::size_t stride_view>
void
FEPointEvaluation<n_components_, dim, spacedim, Number>::test_and_sum(
  const StridedArrayView<ScalarNumber, stride_view> &solution_values,
  const EvaluationFlags::EvaluationFlags            &integration_flags,
  const bool                                         sum_into_values)
{
  do_integrate<false>(solution_values, integration_flags, sum_into_values);
}



template <int n_components_, int dim, int spacedim, typename Number>
void
FEPointEvaluation<n_components_, dim, spacedim, Number>::test_and_sum(
  const ArrayView<ScalarNumber>          &solution_values,
  const EvaluationFlags::EvaluationFlags &integration_flags,
  const bool                              sum_into_values)
{
  test_and_sum(StridedArrayView<ScalarNumber, 1>(solution_values.data(),
                                                 solution_values.size()),
               integration_flags,
               sum_into_values);
}



template <int n_components_, int dim, int spacedim, typename Number>
template <bool is_linear, std::size_t stride_view>
inline void
FEPointEvaluation<n_components_, dim, spacedim, Number>::prepare_evaluate_fast(
  const StridedArrayView<const ScalarNumber, stride_view> &solution_values)
{
  const unsigned int dofs_per_comp =
    is_linear ? Utilities::pow(2, dim) : this->dofs_per_component;

  for (unsigned int comp = 0; comp < n_components; ++comp)
    {
      const std::size_t offset =
        (this->component_in_base_element + comp) * dofs_per_comp;

      if ((is_linear && n_components == 1) || lexicographic_numbering)
        {
          for (unsigned int i = 0; i < dofs_per_comp; ++i)
            ETT::read_value(solution_values[i + offset],
                            comp,
                            this->solution_renumbered[i]);
        }
      else
        {
          const unsigned int *renumber_ptr = this->renumber.data() + offset;
          for (unsigned int i = 0; i < dofs_per_comp; ++i)
            ETT::read_value(solution_values[renumber_ptr[i]],
                            comp,
                            this->solution_renumbered[i]);
        }
    }
}



template <int n_components_, int dim, int spacedim, typename Number>
template <bool is_linear, std::size_t stride_view>
inline void
FEPointEvaluation<n_components_, dim, spacedim, Number>::compute_evaluate_fast(
  const StridedArrayView<const ScalarNumber, stride_view> &solution_values,
  const EvaluationFlags::EvaluationFlags                  &evaluation_flags,
  const unsigned int                                       n_shapes,
  const unsigned int                                       qb,
  vectorized_value_type                                   &value,
  interface_vectorized_unit_gradient_type                 &gradient)
{
  if (evaluation_flags & EvaluationFlags::gradients)
    {
      std::array<vectorized_value_type, dim + 1> result;
      if constexpr (is_linear)
        {
          if constexpr (n_components == 1)
            result =
              internal::evaluate_tensor_product_value_and_gradient_linear<
                dim,
                scalar_value_type,
                VectorizedArrayType,
                1,
                stride_view>(solution_values.data(), this->unit_point_ptr[qb]);
          else
            result =
              internal::evaluate_tensor_product_value_and_gradient_linear(
                this->solution_renumbered.data(), this->unit_point_ptr[qb]);
        }
      else
        result = internal::evaluate_tensor_product_value_and_gradient_shapes<
          dim,
          scalar_value_type,
          VectorizedArrayType,
          1,
          false>(this->shapes.data() + qb * n_shapes,
                 n_shapes,
                 this->solution_renumbered.data());
      gradient[0] = result[0];
      if (dim > 1)
        gradient[1] = result[1];
      if (dim > 2)
        gradient[2] = result[2];
      value = result[dim];
    }
  else
    {
      if constexpr (is_linear)
        {
          if constexpr (n_components == 1)
            value = internal::evaluate_tensor_product_value_linear<
              dim,
              scalar_value_type,
              VectorizedArrayType,
              stride_view>(solution_values.data(), this->unit_point_ptr[qb]);
          else
            value = internal::evaluate_tensor_product_value_linear(
              this->solution_renumbered.data(), this->unit_point_ptr[qb]);
        }
      else
        value =
          internal::evaluate_tensor_product_value_shapes<dim,
                                                         scalar_value_type,
                                                         VectorizedArrayType,
                                                         false>(
            this->shapes.data() + qb * n_shapes,
            n_shapes,
            this->solution_renumbered.data());
    }
}



template <int n_components_, int dim, int spacedim, typename Number>
template <bool is_linear, std::size_t stride_view>
inline void
FEPointEvaluation<n_components_, dim, spacedim, Number>::evaluate_fast(
  const StridedArrayView<const ScalarNumber, stride_view> &solution_values,
  const EvaluationFlags::EvaluationFlags                  &evaluation_flags)
{
  if (!(is_linear && n_components == 1))
    prepare_evaluate_fast<is_linear>(solution_values);

  // loop over quadrature batches qb
  const unsigned int n_shapes = is_linear ? 2 : this->poly.size();

  for (unsigned int qb = 0; qb < this->n_q_batches; ++qb)
    {
      vectorized_value_type                   value;
      interface_vectorized_unit_gradient_type gradient;

      compute_evaluate_fast<is_linear>(
        solution_values, evaluation_flags, n_shapes, qb, value, gradient);

      if (evaluation_flags & EvaluationFlags::values)
        {
          for (unsigned int v = 0, offset = qb * stride;
               v < stride && (stride == 1 || offset < this->n_q_points_scalar);
               ++v, ++offset)
            ETT::set_value(value, v, this->values[offset]);
        }
      if (evaluation_flags & EvaluationFlags::gradients)
        {
          Assert(this->update_flags & update_gradients ||
                   this->update_flags & update_inverse_jacobians,
                 ExcNotInitialized());

          for (unsigned int v = 0, offset = qb * stride;
               v < stride && (stride == 1 || offset < this->n_q_points_scalar);
               ++v, ++offset)
            {
              unit_gradient_type unit_gradient;
              ETT::set_gradient(gradient, v, unit_gradient);
              this->gradients[offset] =
                this->cell_type <=
                    internal::MatrixFreeFunctions::GeometryType::cartesian ?
                  apply_diagonal_transformation(
                    this->inverse_jacobian_ptr[0].transpose(), unit_gradient) :
                  apply_transformation(
                    this
                      ->inverse_jacobian_ptr[this->cell_type <=
                                                 internal::MatrixFreeFunctions::
                                                   GeometryType::affine ?
                                               0 :
                                               offset]
                      .transpose(),
                    unit_gradient);
            }
        }
    }
}



template <int n_components_, int dim, int spacedim, typename Number>
template <std::size_t stride_view>
inline void
FEPointEvaluation<n_components_, dim, spacedim, Number>::evaluate_slow(
  const StridedArrayView<const ScalarNumber, stride_view> &solution_values,
  const EvaluationFlags::EvaluationFlags                  &evaluation_flags)
{
  // slow path with FEValues
  Assert(this->fe_values.get() != nullptr,
         ExcMessage(
           "Not initialized. Please call FEPointEvaluation::reinit()!"));

  const std::size_t n_points = this->fe_values->get_quadrature().size();

  if (evaluation_flags & EvaluationFlags::values)
    {
      this->values.resize(this->n_q_points);
      std::fill(this->values.begin(), this->values.end(), value_type());
      for (unsigned int i = 0; i < this->fe->n_dofs_per_cell(); ++i)
        {
          const ScalarNumber value = solution_values[i];
          for (unsigned int d = 0; d < n_components; ++d)
            if (this->nonzero_shape_function_component[i][d] &&
                (this->fe->is_primitive(i) || this->fe->is_primitive()))
              for (unsigned int qb = 0, q = 0; q < n_points;
                   ++qb, q += n_lanes_user_interface)
                for (unsigned int v = 0;
                     v < n_lanes_user_interface && q + v < n_points;
                     ++v)
                  ETT::access(this->values[qb],
                              v,
                              d,
                              this->fe_values->shape_value(i, q + v) * value);
            else if (this->nonzero_shape_function_component[i][d])
              for (unsigned int qb = 0, q = 0; q < n_points;
                   ++qb, q += n_lanes_user_interface)
                for (unsigned int v = 0;
                     v < n_lanes_user_interface && q + v < n_points;
                     ++v)
                  ETT::access(this->values[qb],
                              v,
                              d,
                              this->fe_values->shape_value_component(i,
                                                                     q + v,
                                                                     d) *
                                value);
        }
    }

  if (evaluation_flags & EvaluationFlags::gradients)
    {
      this->gradients.resize(this->n_q_points);
      std::fill(this->gradients.begin(),
                this->gradients.end(),
                gradient_type());
      for (unsigned int i = 0; i < this->fe->n_dofs_per_cell(); ++i)
        {
          const ScalarNumber value = solution_values[i];
          for (unsigned int d = 0; d < n_components; ++d)
            if (this->nonzero_shape_function_component[i][d] &&
                (this->fe->is_primitive(i) || this->fe->is_primitive()))
              for (unsigned int qb = 0, q = 0; q < n_points;
                   ++qb, q += n_lanes_user_interface)
                for (unsigned int v = 0;
                     v < n_lanes_user_interface && q + v < n_points;
                     ++v)
                  ETT::access(this->gradients[qb],
                              v,
                              d,
                              this->fe_values->shape_grad(i, q + v) * value);
            else if (this->nonzero_shape_function_component[i][d])
              for (unsigned int qb = 0, q = 0; q < n_points;
                   ++qb, q += n_lanes_user_interface)
                for (unsigned int v = 0;
                     v < n_lanes_user_interface && q + v < n_points;
                     ++v)
                  ETT::access(
                    this->gradients[qb],
                    v,
                    d,
                    this->fe_values->shape_grad_component(i, q + v, d) * value);
        }
    }
}



template <int n_components_, int dim, int spacedim, typename Number>
template <bool is_linear>
inline void
FEPointEvaluation<n_components_, dim, spacedim, Number>::compute_integrate_fast(
  const EvaluationFlags::EvaluationFlags       &integration_flags,
  const unsigned int                            n_shapes,
  const unsigned int                            qb,
  const vectorized_value_type                   value,
  const interface_vectorized_unit_gradient_type gradient,
  vectorized_value_type *solution_values_vectorized_linear)
{
  if (integration_flags & EvaluationFlags::gradients)
    internal::integrate_tensor_product_value_and_gradient<
      is_linear,
      dim,
      VectorizedArrayType,
      vectorized_value_type>(this->shapes.data() + qb * n_shapes,
                             n_shapes,
                             &value,
                             gradient,
                             is_linear ?
                               solution_values_vectorized_linear :
                               this->solution_renumbered_vectorized.data(),
                             this->unit_point_ptr[qb],
                             qb != 0);
  else
    internal::integrate_tensor_product_value<is_linear,
                                             dim,
                                             VectorizedArrayType,
                                             vectorized_value_type>(
      this->shapes.data() + qb * n_shapes,
      n_shapes,
      value,
      is_linear ? solution_values_vectorized_linear :
                  this->solution_renumbered_vectorized.data(),
      this->unit_point_ptr[qb],
      qb != 0);
}



template <int n_components_, int dim, int spacedim, typename Number>
template <bool is_linear, std::size_t stride_view>
inline void
FEPointEvaluation<n_components_, dim, spacedim, Number>::finish_integrate_fast(
  const StridedArrayView<ScalarNumber, stride_view> &solution_values,
  vectorized_value_type *solution_values_vectorized_linear,
  const bool             sum_into_values)
{
  if (!sum_into_values && this->fe->n_components() > n_components)
    for (unsigned int i = 0; i < solution_values.size(); ++i)
      solution_values[i] = 0;

  const unsigned int dofs_per_comp =
    is_linear ? Utilities::pow(2, dim) : this->dofs_per_component;

  for (unsigned int comp = 0; comp < n_components; ++comp)
    {
      const std::size_t offset =
        (this->component_in_base_element + comp) * dofs_per_comp;

      if (is_linear || lexicographic_numbering)
        {
          for (unsigned int i = 0; i < dofs_per_comp; ++i)
            if (sum_into_values)
              solution_values[i + offset] +=
                ETT::sum_value(comp,
                               is_linear ?
                                 *(solution_values_vectorized_linear + i) :
                                 this->solution_renumbered_vectorized[i]);
            else
              solution_values[i + offset] =
                ETT::sum_value(comp,
                               is_linear ?
                                 *(solution_values_vectorized_linear + i) :
                                 this->solution_renumbered_vectorized[i]);
        }
      else
        {
          const unsigned int *renumber_ptr = this->renumber.data() + offset;
          for (unsigned int i = 0; i < dofs_per_comp; ++i)
            if (sum_into_values)
              solution_values[renumber_ptr[i]] +=
                ETT::sum_value(comp, this->solution_renumbered_vectorized[i]);
            else
              solution_values[renumber_ptr[i]] =
                ETT::sum_value(comp, this->solution_renumbered_vectorized[i]);
        }
    }
}



template <int n_components_, int dim, int spacedim, typename Number>
template <bool do_JxW, bool is_linear, std::size_t stride_view>
inline void
FEPointEvaluation<n_components_, dim, spacedim, Number>::integrate_fast(
  const StridedArrayView<ScalarNumber, stride_view> &solution_values,
  const EvaluationFlags::EvaluationFlags            &integration_flags,
  const bool                                         sum_into_values)
{
  // zero out lanes of incomplete last quadrature point batch
  if constexpr (stride == 1)
    if (const unsigned int n_filled_lanes =
          this->n_q_points_scalar & (n_lanes_internal - 1);
        n_filled_lanes > 0)
      {
        if (integration_flags & EvaluationFlags::values)
          for (unsigned int v = n_filled_lanes; v < n_lanes_internal; ++v)
            ETT::set_zero_value(this->values.back(), v);
        if (integration_flags & EvaluationFlags::gradients)
          for (unsigned int v = n_filled_lanes; v < n_lanes_internal; ++v)
            ETT::set_zero_gradient(this->gradients.back(), v);
      }

  std::array<vectorized_value_type, is_linear ? Utilities::pow(2, dim) : 0>
    solution_values_vectorized_linear = {};

  // loop over quadrature batches qb
  const unsigned int n_shapes = is_linear ? 2 : this->poly.size();

  const bool cartesian_cell =
    this->cell_type <= internal::MatrixFreeFunctions::GeometryType::cartesian;
  const bool affine_cell =
    this->cell_type <= internal::MatrixFreeFunctions::GeometryType::affine;
  for (unsigned int qb = 0; qb < this->n_q_batches; ++qb)
    {
      vectorized_value_type                 value = {};
      Tensor<1, dim, vectorized_value_type> gradient;

      if (integration_flags & EvaluationFlags::values)
        for (unsigned int v = 0, offset = qb * stride;
             v < stride && (stride == 1 || offset < this->n_q_points_scalar);
             ++v, ++offset)
          ETT::get_value(value,
                         v,
                         do_JxW ? this->values[offset] * this->JxW_ptr[offset] :
                                  this->values[offset]);

      if (integration_flags & EvaluationFlags::gradients)
        for (unsigned int v = 0, offset = qb * stride;
             v < stride && (stride == 1 || offset < this->n_q_points_scalar);
             ++v, ++offset)
          {
            const gradient_type grad_w =
              do_JxW ? this->gradients[offset] * this->JxW_ptr[offset] :
                       this->gradients[offset];
            ETT::get_gradient(
              gradient,
              v,
              cartesian_cell ?
                apply_diagonal_transformation(this->inverse_jacobian_ptr[0],
                                              grad_w) :
                apply_transformation(
                  this->inverse_jacobian_ptr[affine_cell ? 0 : offset],
                  grad_w));
          }

      compute_integrate_fast<is_linear>(
        integration_flags,
        n_shapes,
        qb,
        value,
        gradient,
        solution_values_vectorized_linear.data());
    }

  // add between the lanes and write into the result
  finish_integrate_fast<is_linear>(solution_values,
                                   solution_values_vectorized_linear.data(),
                                   sum_into_values);
}



template <int n_components_, int dim, int spacedim, typename Number>
template <bool do_JxW, std::size_t stride_view>
inline void
FEPointEvaluation<n_components_, dim, spacedim, Number>::integrate_slow(
  const StridedArrayView<ScalarNumber, stride_view> &solution_values,
  const EvaluationFlags::EvaluationFlags            &integration_flags,
  const bool                                         sum_into_values)
{
  // slow path with FEValues
  Assert(this->fe_values.get() != nullptr,
         ExcMessage(
           "Not initialized. Please call FEPointEvaluation::reinit()!"));
  if (!sum_into_values)
    for (unsigned int i = 0; i < solution_values.size(); ++i)
      solution_values[i] = 0;

  const std::size_t n_points = this->fe_values->get_quadrature().size();

  if (integration_flags & EvaluationFlags::values)
    {
      AssertIndexRange(this->n_q_points, this->values.size() + 1);
      for (unsigned int i = 0; i < this->fe->n_dofs_per_cell(); ++i)
        {
          for (unsigned int d = 0; d < n_components; ++d)
            if (this->nonzero_shape_function_component[i][d] &&
                (this->fe->is_primitive(i) || this->fe->is_primitive()))
              for (unsigned int qb = 0, q = 0; q < n_points;
                   ++qb, q += n_lanes_user_interface)
                for (unsigned int v = 0;
                     v < n_lanes_user_interface && q + v < n_points;
                     ++v)
                  solution_values[i] +=
                    this->fe_values->shape_value(i, q + v) *
                    ETT::access(this->values[qb], v, d) *
                    (do_JxW ? this->fe_values->JxW(q + v) : 1.);
            else if (this->nonzero_shape_function_component[i][d])
              for (unsigned int qb = 0, q = 0; q < n_points;
                   ++qb, q += n_lanes_user_interface)
                for (unsigned int v = 0;
                     v < n_lanes_user_interface && q + v < n_points;
                     ++v)
                  solution_values[i] +=
                    this->fe_values->shape_value_component(i, q + v, d) *
                    ETT::access(this->values[qb], v, d) *
                    (do_JxW ? this->fe_values->JxW(q + v) : 1.);
        }
    }

  if (integration_flags & EvaluationFlags::gradients)
    {
      AssertIndexRange(this->n_q_points, this->gradients.size() + 1);
      for (unsigned int i = 0; i < this->fe->n_dofs_per_cell(); ++i)
        {
          for (unsigned int d = 0; d < n_components; ++d)
            if (this->nonzero_shape_function_component[i][d] &&
                (this->fe->is_primitive(i) || this->fe->is_primitive()))
              for (unsigned int qb = 0, q = 0; q < n_points;
                   ++qb, q += n_lanes_user_interface)
                for (unsigned int v = 0;
                     v < n_lanes_user_interface && q + v < n_points;
                     ++v)
                  solution_values[i] +=
                    this->fe_values->shape_grad(i, q + v) *
                    ETT::access(this->gradients[qb], v, d) *
                    (do_JxW ? this->fe_values->JxW(q + v) : 1.);
            else if (this->nonzero_shape_function_component[i][d])
              for (unsigned int qb = 0, q = 0; q < n_points;
                   ++qb, q += n_lanes_user_interface)
                for (unsigned int v = 0;
                     v < n_lanes_user_interface && q + v < n_points;
                     ++v)
                  solution_values[i] +=
                    this->fe_values->shape_grad_component(i, q + v, d) *
                    ETT::access(this->gradients[qb], v, d) *
                    (do_JxW ? this->fe_values->JxW(q + v) : 1.);
        }
    }
}



template <int n_components_, int dim, int spacedim, typename Number>
template <bool do_JxW, std::size_t stride_view>
void
FEPointEvaluation<n_components_, dim, spacedim, Number>::do_integrate(
  const StridedArrayView<ScalarNumber, stride_view> &solution_values,
  const EvaluationFlags::EvaluationFlags            &integration_flags,
  const bool                                         sum_into_values)
{
  if (this->must_reinitialize_pointers)
    internal_reinit_single_cell_state_mapping_info();

  Assert(!(integration_flags & EvaluationFlags::hessians), ExcNotImplemented());

  if (this->n_q_points == 0 || // no evaluation points provided
      !((integration_flags & EvaluationFlags::values) ||
        (integration_flags &
         EvaluationFlags::gradients))) // no integration flags
    {
      if (!sum_into_values)
        for (unsigned int i = 0; i < solution_values.size(); ++i)
          solution_values[i] = 0;
      return;
    }

  Assert(
    !do_JxW || this->JxW_ptr != nullptr,
    ExcMessage(
      "JxW pointer is not set! If you do not want to integrate() use test_and_sum()"));

  AssertDimension(solution_values.size(), this->fe->dofs_per_cell);
  if (this->fast_path)
    {
      if (this->use_linear_path)
        integrate_fast<do_JxW, true>(solution_values,
                                     integration_flags,
                                     sum_into_values);
      else
        integrate_fast<do_JxW, false>(solution_values,
                                      integration_flags,
                                      sum_into_values);
    }
  else
    integrate_slow<do_JxW>(solution_values, integration_flags, sum_into_values);
}



template <int n_components_, int dim, int spacedim, typename Number>
inline Tensor<1, spacedim, Number>
FEPointEvaluation<n_components_, dim, spacedim, Number>::normal_vector(
  const unsigned int point_index) const
{
  AssertIndexRange(point_index, this->n_q_points);
  Assert(this->normal_ptr != nullptr,
         internal::FEPointEvaluation::
           ExcFEPointEvaluationAccessToUninitializedMappingField(
             "update_normal_vectors"));
  return this->normal_ptr[point_index];
}



template <int n_components_, int dim, int spacedim, typename Number>
inline const typename FEPointEvaluation<n_components_,
                                        dim,
                                        spacedim,
                                        Number>::value_type
FEPointEvaluation<n_components_, dim, spacedim, Number>::get_normal_derivative(
  const unsigned int point_index) const
{
  AssertIndexRange(point_index, this->gradients.size());

  value_type normal_derivative;
  if constexpr (n_components == 1)
    normal_derivative =
      this->gradients[point_index] * normal_vector(point_index);
  else
    for (unsigned int comp = 0; comp < n_components; ++comp)
      normal_derivative[comp] =
        this->gradients[point_index][comp] * normal_vector(point_index);

  return normal_derivative;
}



template <int n_components_, int dim, int spacedim, typename Number>
inline void
FEPointEvaluation<n_components_, dim, spacedim, Number>::
  submit_normal_derivative(const value_type  &value,
                           const unsigned int point_index)
{
  AssertIndexRange(point_index, this->gradients.size());
  if constexpr (n_components == 1)
    this->gradients[point_index] = value * normal_vector(point_index);
  else
    for (unsigned int comp = 0; comp < n_components; ++comp)
      this->gradients[point_index][comp] =
        value[comp] * normal_vector(point_index);
}



template <int n_components_, int dim, int spacedim, typename Number>
FEFacePointEvaluation<n_components_, dim, spacedim, Number>::
  FEFacePointEvaluation(
    const NonMatching::MappingInfo<dim, spacedim, Number> &mapping_info,
    const FiniteElement<dim, spacedim>                    &fe,
    const bool                                             is_interior,
    const unsigned int first_selected_component)
  : FEPointEvaluationBase<n_components_, dim, spacedim, Number>(
      mapping_info,
      fe,
      first_selected_component,
      is_interior)
{}



template <int n_components_, int dim, int spacedim, typename Number>
inline void
FEFacePointEvaluation<n_components_, dim, spacedim, Number>::reinit(
  const unsigned int cell_index,
  const unsigned int face_number)
{
  this->current_cell_index         = cell_index;
  this->current_face_number        = face_number;
  this->must_reinitialize_pointers = false;

  if (this->use_linear_path)
    this->template do_reinit<true, true>();
  else
    this->template do_reinit<true, false>();
}



template <int n_components_, int dim, int spacedim, typename Number>
inline void
FEFacePointEvaluation<n_components_, dim, spacedim, Number>::reinit(
  const unsigned int face_index)
{
  this->current_cell_index = face_index;
  this->current_face_number =
    this->mapping_info->get_face_number(face_index, this->is_interior);
  this->must_reinitialize_pointers = false;

  if (this->use_linear_path)
    this->template do_reinit<true, true>();
  else
    this->template do_reinit<true, false>();
}



template <int n_components_, int dim, int spacedim, typename Number>
template <std::size_t stride_view>
void
FEFacePointEvaluation<n_components_, dim, spacedim, Number>::evaluate(
  const StridedArrayView<const ScalarNumber, stride_view> &solution_values,
  const EvaluationFlags::EvaluationFlags                  &evaluation_flags)
{
  Assert(!this->must_reinitialize_pointers,
         ExcMessage("Object has not been reinitialized!"));

  if (this->n_q_points == 0)
    return;

  Assert(!(evaluation_flags & EvaluationFlags::hessians), ExcNotImplemented());

  if (!((evaluation_flags & EvaluationFlags::values) ||
        (evaluation_flags & EvaluationFlags::gradients))) // no evaluation flags
    return;

  AssertDimension(solution_values.size(), this->fe->dofs_per_cell);

  if (this->use_linear_path)
    do_evaluate<true>(solution_values, evaluation_flags);
  else
    do_evaluate<false>(solution_values, evaluation_flags);
}



template <int n_components_, int dim, int spacedim, typename Number>
void
FEFacePointEvaluation<n_components_, dim, spacedim, Number>::evaluate(
  const ArrayView<const ScalarNumber>    &solution_values,
  const EvaluationFlags::EvaluationFlags &evaluation_flags)
{
  evaluate(StridedArrayView<const ScalarNumber, 1>(solution_values.data(),
                                                   solution_values.size()),
           evaluation_flags);
}



template <int n_components_, int dim, int spacedim, typename Number>
template <bool is_linear, std::size_t stride_view>
void
FEFacePointEvaluation<n_components_, dim, spacedim, Number>::do_evaluate(
  const StridedArrayView<const ScalarNumber, stride_view> &solution_values,
  const EvaluationFlags::EvaluationFlags                  &evaluation_flags)
{
  const unsigned int dofs_per_comp =
    is_linear ? Utilities::pow(2, dim) : this->dofs_per_component;

  const ScalarNumber *input;
  if (stride_view == 1 && this->component_in_base_element == 0 &&
      (is_linear || this->renumber.empty()))
    input = solution_values.data();
  else
    {
      for (unsigned int comp = 0; comp < n_components; ++comp)
        {
          const std::size_t offset =
            (this->component_in_base_element + comp) * dofs_per_comp;

          if (is_linear || this->renumber.empty())
            {
              for (unsigned int i = 0; i < dofs_per_comp; ++i)
                this->scratch_data_scalar[i + comp * dofs_per_comp] =
                  solution_values[i + offset];
            }
          else
            {
              const unsigned int *renumber_ptr = this->renumber.data() + offset;
              for (unsigned int i = 0; i < dofs_per_comp; ++i)
                this->scratch_data_scalar[i + comp * dofs_per_comp] =
                  solution_values[renumber_ptr[i]];
            }
        }
      input = this->scratch_data_scalar.data();
    }

  ScalarNumber *output =
    this->scratch_data_scalar.begin() + dofs_per_comp * n_components;

  internal::FEFaceNormalEvaluationImpl<dim, is_linear ? 1 : -1, ScalarNumber>::
    template interpolate<true, false>(n_components,
                                      evaluation_flags,
                                      this->shape_info,
                                      input,
                                      output,
                                      this->current_face_number);

  do_evaluate_in_face<is_linear, 1>(output, evaluation_flags);
}



template <int n_components_, int dim, int spacedim, typename Number>
template <std::size_t stride_view>
void
FEFacePointEvaluation<n_components_, dim, spacedim, Number>::integrate(
  const StridedArrayView<ScalarNumber, stride_view> &solution_values,
  const EvaluationFlags::EvaluationFlags            &integration_flags,
  const bool                                         sum_into_values)
{
  Assert(!this->must_reinitialize_pointers,
         ExcMessage("Object has not been reinitialized!"));

  Assert(!(integration_flags & EvaluationFlags::hessians), ExcNotImplemented());

  if (this->n_q_points == 0 || // no evaluation points provided
      !((integration_flags & EvaluationFlags::values) ||
        (integration_flags &
         EvaluationFlags::gradients))) // no integration flags
    {
      if (!sum_into_values)
        for (unsigned int i = 0; i < solution_values.size(); ++i)
          solution_values[i] = 0;
      return;
    }

  AssertDimension(solution_values.size(), this->fe->dofs_per_cell);

  if (this->use_linear_path)
    do_integrate<true, true>(solution_values,
                             integration_flags,
                             sum_into_values);
  else
    do_integrate<true, false>(solution_values,
                              integration_flags,
                              sum_into_values);
}



template <int n_components_, int dim, int spacedim, typename Number>
void
FEFacePointEvaluation<n_components_, dim, spacedim, Number>::integrate(
  const ArrayView<ScalarNumber>          &solution_values,
  const EvaluationFlags::EvaluationFlags &integration_flags,
  const bool                              sum_into_values)
{
  integrate(StridedArrayView<ScalarNumber, 1>(solution_values.data(),
                                              solution_values.size()),
            integration_flags,
            sum_into_values);
}



template <int n_components_, int dim, int spacedim, typename Number>
template <std::size_t stride_view>
void
FEFacePointEvaluation<n_components_, dim, spacedim, Number>::test_and_sum(
  const StridedArrayView<ScalarNumber, stride_view> &solution_values,
  const EvaluationFlags::EvaluationFlags            &integration_flags,
  const bool                                         sum_into_values)
{
  Assert(!this->must_reinitialize_pointers,
         ExcMessage("Object has not been reinitialized!"));

  Assert(!(integration_flags & EvaluationFlags::hessians), ExcNotImplemented());

  if (this->n_q_points == 0 || // no evaluation points provided
      !((integration_flags & EvaluationFlags::values) ||
        (integration_flags &
         EvaluationFlags::gradients))) // no integration flags
    {
      if (!sum_into_values)
        for (unsigned int i = 0; i < solution_values.size(); ++i)
          solution_values[i] = 0;
      return;
    }

  AssertDimension(solution_values.size(), this->fe->dofs_per_cell);

  if (this->use_linear_path)
    do_integrate<false, true>(solution_values,
                              integration_flags,
                              sum_into_values);
  else
    do_integrate<false, false>(solution_values,
                               integration_flags,
                               sum_into_values);
}



template <int n_components_, int dim, int spacedim, typename Number>
void
FEFacePointEvaluation<n_components_, dim, spacedim, Number>::test_and_sum(
  const ArrayView<ScalarNumber>          &solution_values,
  const EvaluationFlags::EvaluationFlags &integration_flags,
  const bool                              sum_into_values)
{
  test_and_sum(StridedArrayView<ScalarNumber, 1>(solution_values.data(),
                                                 solution_values.size()),
               integration_flags,
               sum_into_values);
}



template <int n_components_, int dim, int spacedim, typename Number>
template <bool do_JxW, bool is_linear, std::size_t stride_view>
void
FEFacePointEvaluation<n_components_, dim, spacedim, Number>::do_integrate(
  const StridedArrayView<ScalarNumber, stride_view> &solution_values,
  const EvaluationFlags::EvaluationFlags            &integration_flags,
  const bool                                         sum_into_values)
{
  if (!sum_into_values && this->fe->n_components() > n_components)
    for (unsigned int i = 0; i < solution_values.size(); ++i)
      solution_values[i] = 0;

  do_integrate_in_face<do_JxW, is_linear, 1>(this->scratch_data_scalar.begin(),
                                             integration_flags,
                                             false);

  ScalarNumber *input = this->scratch_data_scalar.begin();

  if (stride_view == 1 && this->component_in_base_element == 0 &&
      (is_linear || this->renumber.empty()))
    {
      if (sum_into_values)
        internal::
          FEFaceNormalEvaluationImpl<dim, is_linear ? 1 : -1, ScalarNumber>::
            template interpolate<false, true>(n_components,
                                              integration_flags,
                                              this->shape_info,
                                              input,
                                              solution_values.data(),
                                              this->current_face_number);
      else
        internal::
          FEFaceNormalEvaluationImpl<dim, is_linear ? 1 : -1, ScalarNumber>::
            template interpolate<false, false>(n_components,
                                               integration_flags,
                                               this->shape_info,
                                               input,
                                               solution_values.data(),
                                               this->current_face_number);
    }
  else
    {
      const unsigned int dofs_per_comp_face =
        is_linear ? Utilities::pow(2, dim - 1) : this->dofs_per_component_face;

      const unsigned int size_input = 3 * dofs_per_comp_face * n_components;
      ScalarNumber      *output     = input + size_input;

      internal::
        FEFaceNormalEvaluationImpl<dim, is_linear ? 1 : -1, ScalarNumber>::
          template interpolate<false, false>(n_components,
                                             integration_flags,
                                             this->shape_info,
                                             input,
                                             output,
                                             this->current_face_number);

      const unsigned int dofs_per_comp =
        is_linear ? Utilities::pow(2, dim) : this->dofs_per_component;

      for (unsigned int comp = 0; comp < n_components; ++comp)
        {
          const std::size_t offset =
            (this->component_in_base_element + comp) * dofs_per_comp;

          if (is_linear || this->renumber.empty())
            {
              for (unsigned int i = 0; i < dofs_per_comp; ++i)
                if (sum_into_values)
                  solution_values[i + offset] +=
                    output[i + comp * dofs_per_comp];
                else
                  solution_values[i + offset] =
                    output[i + comp * dofs_per_comp];
            }
          else
            {
              const unsigned int *renumber_ptr = this->renumber.data() + offset;
              for (unsigned int i = 0; i < dofs_per_comp; ++i)
                if (sum_into_values)
                  solution_values[renumber_ptr[i]] +=
                    output[i + comp * dofs_per_comp];
                else
                  solution_values[renumber_ptr[i]] =
                    output[i + comp * dofs_per_comp];
            }
        }
    }
}



template <int n_components_, int dim, int spacedim, typename Number>
template <int stride_face_dof>
void
FEFacePointEvaluation<n_components_, dim, spacedim, Number>::evaluate_in_face(
  const ScalarNumber                     *face_dof_values,
  const EvaluationFlags::EvaluationFlags &evaluation_flags)
{
  if (this->use_linear_path)
    do_evaluate_in_face<true, stride_face_dof>(face_dof_values,
                                               evaluation_flags);
  else
    do_evaluate_in_face<false, stride_face_dof>(face_dof_values,
                                                evaluation_flags);
}



template <int n_components_, int dim, int spacedim, typename Number>
template <bool is_linear, int stride_face_dof>
inline void
FEFacePointEvaluation<n_components_, dim, spacedim, Number>::
  do_evaluate_in_face(const ScalarNumber                     *face_dof_values,
                      const EvaluationFlags::EvaluationFlags &evaluation_flags)
{
  const scalar_value_type *face_dof_values_ptr;
  if constexpr (n_components == 1)
    face_dof_values_ptr = face_dof_values;
  else
    {
      const unsigned int dofs_per_comp_face =
        is_linear ? Utilities::pow(2, dim - 1) : this->dofs_per_component_face;
      for (unsigned int comp = 0; comp < n_components; ++comp)
        for (unsigned int i = 0; i < 2 * dofs_per_comp_face; ++i)
          ETT::read_value(face_dof_values[(i + comp * 3 * dofs_per_comp_face) *
                                          stride_face_dof],
                          comp,
                          this->solution_renumbered[i]);

      face_dof_values_ptr = this->solution_renumbered.data();
    }

  constexpr int stride_face_dof_actual =
    n_components == 1 ? stride_face_dof : 1;

  // loop over quadrature batches qb
  const unsigned int n_shapes = is_linear ? 2 : this->poly.size();

  for (unsigned int qb = 0; qb < this->n_q_batches; ++qb)
    {
      vectorized_value_type                   value;
      interface_vectorized_unit_gradient_type gradient;

      if (evaluation_flags & EvaluationFlags::gradients)
        {
          const std::array<vectorized_value_type, dim + 1> interpolated_value =
            is_linear ?
              internal::evaluate_tensor_product_value_and_gradient_linear<
                dim - 1,
                scalar_value_type,
                VectorizedArrayType,
                2,
                stride_face_dof_actual>(face_dof_values_ptr,
                                        this->unit_point_faces_ptr[qb]) :
              internal::evaluate_tensor_product_value_and_gradient_shapes<
                dim - 1,
                scalar_value_type,
                VectorizedArrayType,
                2,
                false,
                stride_face_dof_actual>(this->shapes_faces.data() +
                                          qb * n_shapes,
                                        n_shapes,
                                        face_dof_values_ptr);

          value = interpolated_value[dim - 1];
          // reorder derivative from tangential/normal derivatives into tensor
          // in physical coordinates
          if (this->current_face_number / 2 == 0)
            {
              gradient[0] = interpolated_value[dim];
              if (dim > 1)
                gradient[1] = interpolated_value[0];
              if (dim > 2)
                gradient[2] = interpolated_value[1];
            }
          else if (this->current_face_number / 2 == 1)
            {
              if (dim > 1)
                gradient[1] = interpolated_value[dim];
              if (dim == 3)
                {
                  gradient[0] = interpolated_value[1];
                  gradient[2] = interpolated_value[0];
                }
              else if (dim == 2)
                gradient[0] = interpolated_value[0];
              else
                DEAL_II_ASSERT_UNREACHABLE();
            }
          else if (this->current_face_number / 2 == 2)
            {
              if (dim > 2)
                {
                  gradient[0] = interpolated_value[0];
                  gradient[1] = interpolated_value[1];
                  gradient[2] = interpolated_value[dim];
                }
              else
                DEAL_II_ASSERT_UNREACHABLE();
            }
          else
            DEAL_II_ASSERT_UNREACHABLE();
        }
      else
        {
          value = is_linear ?
                    internal::evaluate_tensor_product_value_linear<
                      dim - 1,
                      scalar_value_type,
                      VectorizedArrayType,
                      stride_face_dof_actual>(face_dof_values_ptr,
                                              this->unit_point_faces_ptr[qb]) :
                    internal::evaluate_tensor_product_value_shapes<
                      dim - 1,
                      scalar_value_type,
                      VectorizedArrayType,
                      false,
                      stride_face_dof_actual>(this->shapes_faces.data() +
                                                qb * n_shapes,
                                              n_shapes,
                                              face_dof_values_ptr);
        }

      if (evaluation_flags & EvaluationFlags::values)
        {
          for (unsigned int v = 0, offset = qb * stride;
               v < stride && (stride == 1 || offset < this->n_q_points_scalar);
               ++v, ++offset)
            ETT::set_value(value, v, this->values[offset]);
        }
      if (evaluation_flags & EvaluationFlags::gradients)
        {
          Assert(this->update_flags & update_gradients ||
                   this->update_flags & update_inverse_jacobians,
                 ExcNotInitialized());

          for (unsigned int v = 0, offset = qb * stride;
               v < stride && (stride == 1 || offset < this->n_q_points_scalar);
               ++v, ++offset)
            {
              unit_gradient_type unit_gradient;
              ETT::set_gradient(gradient, v, unit_gradient);
              this->gradients[offset] =
                this->cell_type <=
                    internal::MatrixFreeFunctions::GeometryType::cartesian ?
                  apply_diagonal_transformation(
                    this->inverse_jacobian_ptr[0].transpose(), unit_gradient) :
                  apply_transformation(
                    this
                      ->inverse_jacobian_ptr[this->cell_type <=
                                                 internal::MatrixFreeFunctions::
                                                   GeometryType::affine ?
                                               0 :
                                               offset]
                      .transpose(),
                    unit_gradient);
            }
        }
    }
}



template <int n_components_, int dim, int spacedim, typename Number>
template <int stride_face_dof>
void
FEFacePointEvaluation<n_components_, dim, spacedim, Number>::integrate_in_face(
  ScalarNumber                           *face_dof_values,
  const EvaluationFlags::EvaluationFlags &integration_flags,
  const bool                              sum_into_values)
{
  if (this->use_linear_path)
    do_integrate_in_face<true, true, stride_face_dof>(face_dof_values,
                                                      integration_flags,
                                                      sum_into_values);
  else
    do_integrate_in_face<true, false, stride_face_dof>(face_dof_values,
                                                       integration_flags,
                                                       sum_into_values);
}



template <int n_components_, int dim, int spacedim, typename Number>
template <bool do_JxW, bool is_linear, int stride_face_dof>
inline void
FEFacePointEvaluation<n_components_, dim, spacedim, Number>::
  do_integrate_in_face(
    ScalarNumber                           *face_dof_values,
    const EvaluationFlags::EvaluationFlags &integration_flags,
    const bool                              sum_into_values)
{
  // zero out lanes of incomplete last quadrature point batch
  if constexpr (stride == 1)
    if (const unsigned int n_filled_lanes =
          this->n_q_points_scalar & (n_lanes_internal - 1);
        n_filled_lanes > 0)
      {
        if (integration_flags & EvaluationFlags::values)
          for (unsigned int v = n_filled_lanes; v < n_lanes_internal; ++v)
            ETT::set_zero_value(this->values.back(), v);
        if (integration_flags & EvaluationFlags::gradients)
          for (unsigned int v = n_filled_lanes; v < n_lanes_internal; ++v)
            ETT::set_zero_gradient(this->gradients.back(), v);
      }

  std::array<vectorized_value_type,
             is_linear ? 2 * Utilities::pow(2, dim - 1) : 0>
    solution_values_vectorized_linear = {};

  // loop over quadrature batches qb
  const unsigned int n_shapes = is_linear ? 2 : this->poly.size();

  const bool cartesian_cell =
    this->cell_type <= internal::MatrixFreeFunctions::GeometryType::cartesian;
  const bool affine_cell =
    this->cell_type <= internal::MatrixFreeFunctions::GeometryType::affine;
  for (unsigned int qb = 0; qb < this->n_q_batches; ++qb)
    {
      vectorized_value_type                 value = {};
      Tensor<1, dim, vectorized_value_type> gradient;

      if (integration_flags & EvaluationFlags::values)
        for (unsigned int v = 0, offset = qb * stride;
             v < stride && (stride == 1 || offset < this->n_q_points_scalar);
             ++v, ++offset)
          ETT::get_value(value,
                         v,
                         do_JxW ? this->values[offset] * this->JxW_ptr[offset] :
                                  this->values[offset]);

      if (integration_flags & EvaluationFlags::gradients)
        for (unsigned int v = 0, offset = qb * stride;
             v < stride && (stride == 1 || offset < this->n_q_points_scalar);
             ++v, ++offset)
          {
            const auto grad_w =
              do_JxW ? this->gradients[offset] * this->JxW_ptr[offset] :
                       this->gradients[offset];
            ETT::get_gradient(
              gradient,
              v,
              cartesian_cell ?
                apply_diagonal_transformation(this->inverse_jacobian_ptr[0],
                                              grad_w) :
                apply_transformation(
                  this->inverse_jacobian_ptr[affine_cell ? 0 : offset],
                  grad_w));
          }

      if (integration_flags & EvaluationFlags::gradients)
        {
          std::array<vectorized_value_type, 2>      value_face = {};
          Tensor<1, dim - 1, vectorized_value_type> gradient_in_face;

          value_face[0] = value;
          // fill derivative in physical coordinates into tangential/normal
          // derivatives
          if (this->current_face_number / 2 == 0)
            {
              value_face[1] = gradient[0];
              if (dim > 1)
                gradient_in_face[0] = gradient[1];
              if (dim > 2)
                gradient_in_face[1] = gradient[2];
            }
          else if (this->current_face_number / 2 == 1)
            {
              if (dim > 1)
                value_face[1] = gradient[1];
              if (dim == 3)
                {
                  gradient_in_face[0] = gradient[2];
                  gradient_in_face[1] = gradient[0];
                }
              else if (dim == 2)
                gradient_in_face[0] = gradient[0];
              else
                DEAL_II_ASSERT_UNREACHABLE();
            }
          else if (this->current_face_number / 2 == 2)
            {
              if (dim > 2)
                {
                  value_face[1]       = gradient[2];
                  gradient_in_face[0] = gradient[0];
                  gradient_in_face[1] = gradient[1];
                }
              else
                DEAL_II_ASSERT_UNREACHABLE();
            }
          else
            DEAL_II_ASSERT_UNREACHABLE();

          internal::integrate_tensor_product_value_and_gradient<
            is_linear,
            dim - 1,
            VectorizedArrayType,
            vectorized_value_type,
            2>(this->shapes_faces.data() + qb * n_shapes,
               n_shapes,
               value_face.data(),
               gradient_in_face,
               is_linear ? solution_values_vectorized_linear.data() :
                           this->solution_renumbered_vectorized.data(),
               this->unit_point_faces_ptr[qb],
               qb != 0);
        }
      else
        internal::integrate_tensor_product_value<is_linear,
                                                 dim - 1,
                                                 VectorizedArrayType,
                                                 vectorized_value_type>(
          this->shapes_faces.data() + qb * n_shapes,
          n_shapes,
          value,
          is_linear ? solution_values_vectorized_linear.data() :
                      this->solution_renumbered_vectorized.data(),
          this->unit_point_faces_ptr[qb],
          qb != 0);
    }

  const unsigned int dofs_per_comp_face =
    is_linear ? Utilities::pow(2, dim - 1) : this->dofs_per_component_face;

  for (unsigned int comp = 0; comp < n_components; ++comp)
    for (unsigned int i = 0; i < 2 * dofs_per_comp_face; ++i)
      if (sum_into_values)
        face_dof_values[(i + comp * 3 * dofs_per_comp_face) *
                        stride_face_dof] +=
          ETT::sum_value(comp,
                         is_linear ?
                           *(solution_values_vectorized_linear.data() + i) :
                           this->solution_renumbered_vectorized[i]);
      else
        face_dof_values[(i + comp * 3 * dofs_per_comp_face) * stride_face_dof] =
          ETT::sum_value(comp,
                         is_linear ?
                           *(solution_values_vectorized_linear.data() + i) :
                           this->solution_renumbered_vectorized[i]);
}



template <int n_components_, int dim, int spacedim, typename Number>
inline Tensor<1, spacedim, Number>
FEFacePointEvaluation<n_components_, dim, spacedim, Number>::normal_vector(
  const unsigned int point_index) const
{
  AssertIndexRange(point_index, this->n_q_points);
  Assert(this->normal_ptr != nullptr,
         internal::FEPointEvaluation::
           ExcFEPointEvaluationAccessToUninitializedMappingField(
             "update_normal_vectors"));
  if (this->cell_type <= dealii::internal::MatrixFreeFunctions::affine)
    {
      Tensor<1, spacedim, Number> normal;
      for (unsigned int d = 0; d < dim; ++d)
        normal[d] =
          internal::VectorizedArrayTrait<Number>::get(this->normal_ptr[0][d],
                                                      0);

      return normal;
    }
  else
    {
      return this->normal_ptr[point_index];
    }
}



template <int n_components_, int dim, int spacedim, typename Number>
inline const typename FEFacePointEvaluation<n_components_,
                                            dim,
                                            spacedim,
                                            Number>::value_type
FEFacePointEvaluation<n_components_, dim, spacedim, Number>::
  get_normal_derivative(const unsigned int point_index) const
{
  AssertIndexRange(point_index, this->gradients.size());

  value_type normal_derivative;
  if constexpr (n_components == 1)
    normal_derivative =
      this->gradients[point_index] * normal_vector(point_index);
  else
    for (unsigned int comp = 0; comp < n_components; ++comp)
      normal_derivative[comp] =
        this->gradients[point_index][comp] * normal_vector(point_index);

  return normal_derivative;
}



template <int n_components_, int dim, int spacedim, typename Number>
inline void
FEFacePointEvaluation<n_components_, dim, spacedim, Number>::
  submit_normal_derivative(const value_type  &value,
                           const unsigned int point_index)
{
  AssertIndexRange(point_index, this->gradients.size());
  if constexpr (n_components == 1)
    this->gradients[point_index] = value * normal_vector(point_index);
  else
    for (unsigned int comp = 0; comp < n_components; ++comp)
      this->gradients[point_index][comp] =
        value[comp] * normal_vector(point_index);
}

DEAL_II_NAMESPACE_CLOSE

#endif
