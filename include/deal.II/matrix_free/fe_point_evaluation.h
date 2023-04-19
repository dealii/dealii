// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2022 by the deal.II authors
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
#include <deal.II/matrix_free/shape_info.h>
#include <deal.II/matrix_free/tensor_product_kernels.h>

#include <deal.II/non_matching/mapping_info.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace FEPointEvaluation
  {
    template <typename Number>
    struct VectorizedPointsTypeTraits
    {
      using value_type                    = VectorizedArray<Number>;
      static constexpr std::size_t stride = value_type::size();
    };

    template <typename Number, std::size_t width_>
    struct VectorizedPointsTypeTraits<VectorizedArray<Number, width_>>
    {
      using value_type                    = VectorizedArray<Number, width_>;
      static constexpr std::size_t stride = 1;
    };


    /**
     * Struct to distinguish between the value and gradient types of different
     * numbers of components used by the FlexibleEvaluator class.
     */
    template <int dim, int n_components, typename Number>
    struct EvaluatorTypeTraits
    {
      using ScalarNumber =
        typename internal::VectorizedArrayTrait<Number>::value_type;
      using VectorizedArrayType =
        typename internal::FEPointEvaluation::VectorizedPointsTypeTraits<
          Number>::value_type;
      using value_type        = Tensor<1, n_components, Number>;
      using scalar_value_type = Tensor<1, n_components, ScalarNumber>;
      using vectorized_value_type =
        Tensor<1, n_components, VectorizedArrayType>;
      using gradient_type = Tensor<1, n_components, Tensor<1, dim, Number>>;
      using scalar_gradient_type =
        Tensor<1, n_components, Tensor<1, dim, ScalarNumber>>;
      using vectorized_gradient_type =
        Tensor<1, n_components, Tensor<1, dim, VectorizedArrayType>>;
      using interface_vectorized_gradient_type =
        Tensor<1, dim, Tensor<1, n_components, VectorizedArrayType>>;

      static void
      read_value(const Number       vector_entry,
                 const unsigned int component,
                 value_type &       result)
      {
        AssertIndexRange(component, n_components);
        result[component] = vector_entry;
      }

      static void
      write_value(Number &           vector_entry,
                  const unsigned int component,
                  const value_type & result)
      {
        AssertIndexRange(component, n_components);
        vector_entry = result[component];
      }

      static void
      set_gradient(const interface_vectorized_gradient_type &value,
                   const unsigned int                        vector_lane,
                   scalar_gradient_type &                    result)
      {
        for (unsigned int i = 0; i < n_components; ++i)
          for (unsigned int d = 0; d < dim; ++d)
            result[i][d] = value[d][i][vector_lane];
      }

      static void
      set_gradient(const interface_vectorized_gradient_type &value,
                   const unsigned int,
                   vectorized_gradient_type &result)
      {
        for (unsigned int i = 0; i < n_components; ++i)
          for (unsigned int d = 0; d < dim; ++d)
            result[i][d] = value[d][i];
      }

      static void
      get_gradient(interface_vectorized_gradient_type &value,
                   const unsigned int                  vector_lane,
                   const scalar_gradient_type &        result)
      {
        for (unsigned int i = 0; i < n_components; ++i)
          for (unsigned int d = 0; d < dim; ++d)
            value[d][i][vector_lane] = result[i][d];
      }

      static void
      get_gradient(interface_vectorized_gradient_type &value,
                   const unsigned int,
                   const vectorized_gradient_type &result)
      {
        for (unsigned int i = 0; i < n_components; ++i)
          for (unsigned int d = 0; d < dim; ++d)
            value[d][i] = result[i][d];
      }

      static void
      set_zero_gradient(gradient_type &value, const unsigned int vector_lane)
      {
        for (unsigned int i = 0; i < n_components; ++i)
          for (unsigned int d = 0; d < dim; ++d)
            internal::VectorizedArrayTrait<Number>::get(value[i][d],
                                                        vector_lane) = 0.;
      }

      static void
      set_value(const vectorized_value_type &value,
                const unsigned int           vector_lane,
                scalar_value_type &          result)
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
      get_value(vectorized_value_type &  value,
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
      access(value_type &        value,
             const unsigned int  vector_lane,
             const unsigned int  component,
             const ScalarNumber &shape_value)
      {
        internal::VectorizedArrayTrait<Number>::get(value[component],
                                                    vector_lane) += shape_value;
      }

      static ScalarNumber
      access(const value_type & value,
             const unsigned int vector_lane,
             const unsigned int component)
      {
        return internal::VectorizedArrayTrait<Number>::get(value[component],
                                                           vector_lane);
      }

      static void
      access(gradient_type &                     value,
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
      access(const gradient_type &value,
             const unsigned int   vector_lane,
             const unsigned int   component)
      {
        Tensor<1, dim, ScalarNumber> result;
        for (unsigned int d = 0; d < dim; ++d)
          result[d] =
            internal::VectorizedArrayTrait<Number>::get(value[component][d],
                                                        vector_lane);
        return result;
      }
    };

    template <int dim, typename Number>
    struct EvaluatorTypeTraits<dim, 1, Number>
    {
      using ScalarNumber =
        typename internal::VectorizedArrayTrait<Number>::value_type;
      using VectorizedArrayType =
        typename internal::FEPointEvaluation::VectorizedPointsTypeTraits<
          Number>::value_type;
      using value_type               = Number;
      using scalar_value_type        = ScalarNumber;
      using vectorized_value_type    = VectorizedArrayType;
      using gradient_type            = Tensor<1, dim, Number>;
      using scalar_gradient_type     = Tensor<1, dim, ScalarNumber>;
      using vectorized_gradient_type = Tensor<1, dim, VectorizedArrayType>;

      static void
      read_value(const Number vector_entry,
                 const unsigned int,
                 value_type &result)
      {
        result = vector_entry;
      }

      static void
      write_value(Number &vector_entry,
                  const unsigned int,
                  const value_type &result)
      {
        vector_entry = result;
      }

      static void
      set_gradient(const vectorized_gradient_type &value,
                   const unsigned int              vector_lane,
                   scalar_gradient_type &          result)
      {
        for (unsigned int d = 0; d < dim; ++d)
          result[d] = value[d][vector_lane];
      }

      static void
      set_gradient(const vectorized_gradient_type &value,
                   const unsigned int,
                   vectorized_gradient_type &result)
      {
        result = value;
      }

      static void
      get_gradient(vectorized_gradient_type &  value,
                   const unsigned int          vector_lane,
                   const scalar_gradient_type &result)
      {
        for (unsigned int d = 0; d < dim; ++d)
          value[d][vector_lane] = result[d];
      }

      static void
      get_gradient(vectorized_gradient_type &value,
                   const unsigned int,
                   const vectorized_gradient_type &result)
      {
        value = result;
      }

      static void
      set_zero_gradient(gradient_type &value, const unsigned int vector_lane)
      {
        for (unsigned int d = 0; d < dim; ++d)
          internal::VectorizedArrayTrait<Number>::get(value[d], vector_lane) =
            0.;
      }

      static void
      set_value(const vectorized_value_type &value,
                const unsigned int           vector_lane,
                scalar_value_type &          result)
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
      get_value(vectorized_value_type &  value,
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
      access(value_type &       value,
             const unsigned int vector_lane,
             const unsigned int,
             const ScalarNumber &shape_value)
      {
        internal::VectorizedArrayTrait<Number>::get(value, vector_lane) +=
          shape_value;
      }

      static ScalarNumber
      access(const value_type & value,
             const unsigned int vector_lane,
             const unsigned int)
      {
        return internal::VectorizedArrayTrait<Number>::get(value, vector_lane);
      }

      static void
      access(gradient_type &    value,
             const unsigned int vector_lane,
             const unsigned int,
             const scalar_gradient_type &shape_gradient)
      {
        for (unsigned int d = 0; d < dim; ++d)
          internal::VectorizedArrayTrait<Number>::get(value[d], vector_lane) +=
            shape_gradient[d];
      }

      static scalar_gradient_type
      access(const gradient_type &value,
             const unsigned int   vector_lane,
             const unsigned int)
      {
        scalar_gradient_type result;
        for (unsigned int d = 0; d < dim; ++d)
          result[d] =
            internal::VectorizedArrayTrait<Number>::get(value[d], vector_lane);
        return result;
      }
    };

    template <int dim, typename Number>
    struct EvaluatorTypeTraits<dim, dim, Number>
    {
      using ScalarNumber =
        typename internal::VectorizedArrayTrait<Number>::value_type;
      using VectorizedArrayType =
        typename internal::FEPointEvaluation::VectorizedPointsTypeTraits<
          Number>::value_type;
      using value_type               = Tensor<1, dim, Number>;
      using scalar_value_type        = Tensor<1, dim, ScalarNumber>;
      using vectorized_value_type    = Tensor<1, dim, VectorizedArrayType>;
      using gradient_type            = Tensor<2, dim, Number>;
      using scalar_gradient_type     = Tensor<2, dim, ScalarNumber>;
      using vectorized_gradient_type = Tensor<2, dim, VectorizedArrayType>;
      using interface_vectorized_gradient_type =
        Tensor<1, dim, Tensor<1, dim, VectorizedArrayType>>;

      static void
      read_value(const Number       vector_entry,
                 const unsigned int component,
                 value_type &       result)
      {
        result[component] = vector_entry;
      }

      static void
      write_value(Number &           vector_entry,
                  const unsigned int component,
                  const value_type & result)
      {
        vector_entry = result[component];
      }

      static void
      set_gradient(const interface_vectorized_gradient_type &value,
                   const unsigned int                        vector_lane,
                   scalar_gradient_type &                    result)
      {
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int d = 0; d < dim; ++d)
            result[i][d] = value[d][i][vector_lane];
      }

      static void
      set_gradient(const interface_vectorized_gradient_type &value,
                   const unsigned int,
                   vectorized_gradient_type &result)
      {
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int d = 0; d < dim; ++d)
            result[i][d] = value[d][i];
      }

      static void
      get_gradient(interface_vectorized_gradient_type &value,
                   const unsigned int                  vector_lane,
                   const scalar_gradient_type &        result)
      {
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int d = 0; d < dim; ++d)
            value[d][i][vector_lane] = result[i][d];
      }

      static void
      get_gradient(interface_vectorized_gradient_type &value,
                   const unsigned int,
                   const vectorized_gradient_type &result)
      {
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int d = 0; d < dim; ++d)
            value[d][i] = result[i][d];
      }

      static void
      set_zero_gradient(gradient_type &value, const unsigned int vector_lane)
      {
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int d = 0; d < dim; ++d)
            internal::VectorizedArrayTrait<Number>::get(value[i][d],
                                                        vector_lane) = 0.;
      }

      static void
      set_value(const vectorized_value_type &value,
                const unsigned int           vector_lane,
                scalar_value_type &          result)
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
      get_value(vectorized_value_type &  value,
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
      access(value_type &        value,
             const unsigned int  vector_lane,
             const unsigned int  component,
             const ScalarNumber &shape_value)
      {
        internal::VectorizedArrayTrait<Number>::get(value[component],
                                                    vector_lane) += shape_value;
      }

      static ScalarNumber
      access(const value_type & value,
             const unsigned int vector_lane,
             const unsigned int component)
      {
        return internal::VectorizedArrayTrait<Number>::get(value[component],
                                                           vector_lane);
      }

      static void
      access(gradient_type &                     value,
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
      access(const gradient_type &value,
             const unsigned int   vector_lane,
             const unsigned int   component)
      {
        Tensor<1, dim, ScalarNumber> result;
        for (unsigned int d = 0; d < dim; ++d)
          result[d] =
            internal::VectorizedArrayTrait<Number>::get(value[component][d],
                                                        vector_lane);
        return result;
      }
    };

    template <typename Number>
    struct EvaluatorTypeTraits<1, 1, Number>
    {
      using ScalarNumber =
        typename internal::VectorizedArrayTrait<Number>::value_type;
      using VectorizedArrayType =
        typename internal::FEPointEvaluation::VectorizedPointsTypeTraits<
          Number>::value_type;
      using value_type               = Number;
      using scalar_value_type        = ScalarNumber;
      using vectorized_value_type    = VectorizedArrayType;
      using gradient_type            = Tensor<1, 1, Number>;
      using scalar_gradient_type     = Tensor<1, 1, ScalarNumber>;
      using vectorized_gradient_type = Tensor<1, 1, VectorizedArrayType>;

      static void
      read_value(const Number vector_entry,
                 const unsigned int,
                 value_type &result)
      {
        result = vector_entry;
      }

      static void
      write_value(Number &vector_entry,
                  const unsigned int,
                  const value_type &result)
      {
        vector_entry = result;
      }

      static void
      set_gradient(const vectorized_gradient_type &value,
                   const unsigned int              vector_lane,
                   scalar_gradient_type &          result)
      {
        result[0] = value[0][vector_lane];
      }

      static void
      set_gradient(const vectorized_gradient_type &value,
                   const unsigned int,
                   vectorized_gradient_type &result)
      {
        result = value;
      }

      static void
      get_gradient(vectorized_gradient_type &  value,
                   const unsigned int          vector_lane,
                   const scalar_gradient_type &result)
      {
        value[0][vector_lane] = result[0];
      }

      static void
      get_gradient(vectorized_gradient_type &value,
                   const unsigned int,
                   const vectorized_gradient_type &result)
      {
        value = result;
      }

      static void
      set_zero_gradient(gradient_type &value, const unsigned int vector_lane)
      {
        internal::VectorizedArrayTrait<Number>::get(value[0], vector_lane) = 0.;
      }

      static void
      set_value(const vectorized_value_type &value,
                const unsigned int           vector_lane,
                scalar_value_type &          result)
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
      get_value(vectorized_value_type &  value,
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
      access(value_type &       value,
             const unsigned int vector_lane,
             const unsigned int,
             const ScalarNumber &shape_value)
      {
        internal::VectorizedArrayTrait<Number>::get(value, vector_lane) +=
          shape_value;
      }

      static ScalarNumber
      access(const value_type & value,
             const unsigned int vector_lane,
             const unsigned int)
      {
        return internal::VectorizedArrayTrait<Number>::get(value, vector_lane);
      }

      static void
      access(gradient_type &    value,
             const unsigned int vector_lane,
             const unsigned int,
             const scalar_gradient_type &shape_gradient)
      {
        internal::VectorizedArrayTrait<Number>::get(value[0], vector_lane) +=
          shape_gradient[0];
      }

      static scalar_gradient_type
      access(const gradient_type &value,
             const unsigned int   vector_lane,
             const unsigned int)
      {
        scalar_gradient_type result;
        result[0] =
          internal::VectorizedArrayTrait<Number>::get(value[0], vector_lane);
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
 * module. In those cases, the cost implied
 * by this class is similar (or sometimes even somewhat lower) than using
 * `FEValues::reinit(cell)` followed by `FEValues::get_function_gradients`.
 */
template <int n_components_,
          int dim,
          int spacedim    = dim,
          typename Number = double>
class FEPointEvaluation
{
public:
  static constexpr unsigned int dimension    = dim;
  static constexpr unsigned int n_components = n_components_;

  using number_type = Number;

  using ScalarNumber =
    typename internal::VectorizedArrayTrait<Number>::value_type;
  using VectorizedArrayType =
    typename internal::FEPointEvaluation::VectorizedPointsTypeTraits<
      Number>::value_type;
  using ETT = typename internal::FEPointEvaluation::
    EvaluatorTypeTraits<dim, n_components, Number>;
  using value_type            = typename ETT::value_type;
  using scalar_value_type     = typename ETT::scalar_value_type;
  using vectorized_value_type = typename ETT::vectorized_value_type;
  using gradient_type         = typename ETT::gradient_type;

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
  FEPointEvaluation(const Mapping<dim> &      mapping,
                    const FiniteElement<dim> &fe,
                    const UpdateFlags         update_flags,
                    const unsigned int        first_selected_component = 0);

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
   */
  FEPointEvaluation(NonMatching::MappingInfo<dim, spacedim> &mapping_info,
                    const FiniteElement<dim> &               fe,
                    const unsigned int first_selected_component = 0);

  /**
   * Copy constructor.
   */
  FEPointEvaluation(FEPointEvaluation &other) noexcept;

  /**
   * Move constructor.
   */
  FEPointEvaluation(FEPointEvaluation &&other) noexcept;

  /**
   * Destructor.
   */
  ~FEPointEvaluation();

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
   * Reinitialize the evaluator to point to the correct precomputed mapping of
   * the face in the MappingInfo object.
   */
  void
  reinit(const unsigned int cell_index, const unsigned int face_number);

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
  evaluate(const ArrayView<const ScalarNumber> &   solution_values,
           const EvaluationFlags::EvaluationFlags &evaluation_flags);

  /**
   * This function multiplies the quantities passed in by previous
   * submit_value() or submit_gradient() calls by the value or gradient of the
   * test functions, and performs summation over all given points. This is
   * similar to the integration of a bilinear form in terms of the test
   * function, with the difference that this formula does not include a `JxW`
   * factor. This allows the class to naturally embed point information
   * (e.g. particles) into a finite element formulation. Of course, by
   * multiplication of a `JxW` information of the data given to
   * submit_value(), the integration can also be represented by this class.
   *
   * @param[out] solution_values This array will contain the result of the
   * integral, which can be used to during
   * `cell->set_dof_values(solution_values, global_vector)` or
   * `cell->distribute_local_to_global(solution_values, global_vector)`. Note
   * that for multi-component systems where only some of the components are
   * selected by the present class, the entries not touched by this class will
   * be zeroed out.
   *
   * @param[in] integration_flags Flags specifying which quantities should be
   * integrated at the points.
   *
   */
  void
  integrate(const ArrayView<ScalarNumber> &         solution_values,
            const EvaluationFlags::EvaluationFlags &integration_flags);

  /**
   * Return the value at quadrature point number @p point_index after a call to
   * FEPointEvaluation::evaluate() with EvaluationFlags::values set, or
   * the value that has been stored there with a call to
   * FEPointEvaluation::submit_value(). If the object is vector-valued, a
   * vector-valued return argument is given.
   */
  const value_type &
  get_value(const unsigned int point_index) const;

  /**
   * Write a value to the field containing the values on points
   * with component point_index. Access to the same field as through
   * get_value(). If applied before the function FEPointEvaluation::integrate()
   * with EvaluationFlags::values set is called, this specifies the value
   * which is tested by all basis function on the current cell and
   * integrated over.
   */
  void
  submit_value(const value_type &value, const unsigned int point_index);

  /**
   * Return the gradient in real coordinates at the point with index
   * `point_index` after a call to FEPointEvaluation::evaluate() with
   * EvaluationFlags::gradients set, or the gradient that has been stored there
   * with a call to FEPointEvaluation::submit_gradient(). The gradient in real
   * coordinates is obtained by taking the unit gradient (also accessible via
   * get_unit_gradient()) and applying the inverse Jacobian of the mapping. If
   * the object is vector-valued, a vector-valued return argument is given.
   */
  const gradient_type &
  get_gradient(const unsigned int point_index) const;

  /**
   * Return the gradient in unit coordinates at the point with index
   * `point_index` after a call to FEPointEvaluation::evaluate() with
   * EvaluationFlags::gradients set, or the gradient that has been stored there
   * with a call to FEPointEvaluation::submit_gradient(). If the object is
   * vector-valued, a vector-valued return argument is given. Note that when
   * vectorization is enabled, values from several points are grouped
   * together.
   */
  const gradient_type &
  get_unit_gradient(const unsigned int point_index) const;

  /**
   * Write a contribution that is tested by the gradient to the field
   * containing the values on points with the given `point_index`. Access to
   * the same field as through get_gradient(). If applied before the function
   * FEPointEvaluation::integrate(EvaluationFlags::gradients) is called, this
   * specifies what is tested by all basis function gradients on the current
   * cell and integrated over.
   */
  void
  submit_gradient(const gradient_type &, const unsigned int point_index);

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
   * Return the normal vector. This class or the MappingInfo object passed to
   * this function needs to be constructed with UpdateFlags containing
   * `update_normal_vectors`.
   */
  Tensor<1, spacedim, Number>
  normal_vector(const unsigned int point_index) const;

  /**
   * Return the position in real coordinates of the given point index among
   * the points passed to reinit().
   */
  Point<spacedim, Number>
  real_point(const unsigned int point_index) const;

  /**
   * Return the position in unit/reference coordinates of the given point
   * index, i.e., the respective point passed to the reinit() function.
   */
  Point<dim, Number>
  unit_point(const unsigned int point_index) const;

  /**
   * Return an object that can be thought of as an array containing all indices
   * from zero to n_quadrature_points. This allows to write code using
   * range-based for loops.
   */
  inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  quadrature_point_indices() const;

private:
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
  void
  do_reinit();

  /**
   * Fast path of the evaluate function.
   */
  void
  evaluate_fast(const ArrayView<const ScalarNumber> &   solution_values,
                const EvaluationFlags::EvaluationFlags &evaluation_flags);

  /**
   * Slow path of the evaluate function using FEValues.
   */
  void
  evaluate_slow(const ArrayView<const ScalarNumber> &   solution_values,
                const EvaluationFlags::EvaluationFlags &evaluation_flags);

  /**
   * Fast path of the integrate function.
   */
  void
  integrate_fast(const ArrayView<ScalarNumber> &         solution_values,
                 const EvaluationFlags::EvaluationFlags &integration_flags);

  /**
   * Slow path of the integrate function using FEValues.
   */
  void
  integrate_slow(const ArrayView<ScalarNumber> &         solution_values,
                 const EvaluationFlags::EvaluationFlags &integration_flags);

  /**
   * Number of quadrature points of the current cell/face.
   */
  const unsigned int n_q_points;

  /**
   * Number of active quadrature points of the last quadrature point batch.
   */
  const unsigned int n_filled_lanes_last_batch;

  /**
   * Pointer to the Mapping object passed to the constructor.
   */
  SmartPointer<const Mapping<dim, spacedim>> mapping;

  /**
   * Pointer to the FiniteElement object passed to the constructor.
   */
  SmartPointer<const FiniteElement<dim>> fe;

  /**
   * Description of the 1d polynomial basis for tensor product elements used
   * for the fast path of this class using tensor product evaluators.
   */
  std::vector<Polynomials::Polynomial<double>> poly;

  /**
   * Store whether the polynomials are linear with nodes at 0 and 1.
   */
  bool polynomials_are_hat_functions;

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
   * Temporary array to store the values at the points.
   */
  std::vector<value_type> values;

  /**
   * Temporary array to store the gradients in unit coordinates at the points.
   */
  std::vector<gradient_type> unit_gradients;

  /**
   * Temporary array to store the gradients in real coordinates at the points.
   */
  std::vector<gradient_type> gradients;

  /**
   * Number of unknowns per component, i.e., number of unique basis functions,
   * for the chosen FiniteElement (or base element).
   */
  unsigned int dofs_per_component;

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
  std::unique_ptr<NonMatching::MappingInfo<dim, spacedim>>
    mapping_info_on_the_fly;

  /**
   * Pointer to currently used mapping info (either on the fly or external
   * precomputed).
   */
  SmartPointer<NonMatching::MappingInfo<dim, spacedim>> mapping_info;

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
   * Connection to NonMatching::MappingInfo to check wheter mapping data
   * has been invalidated.
   */
  boost::signals2::connection connection_is_reinitialized;

  /**
   * Bool indicating if class is reinitialized and data vectors a resized.
   */
  bool is_reinitialized;

  /**
   * Vector containing tensor product shape functions evaluated (during
   * reinit()) at the vectorized unit points.
   */
  AlignedVector<dealii::ndarray<VectorizedArrayType, 2, dim>> shapes;
};

// ----------------------- template and inline function ----------------------


template <int n_components, int dim, int spacedim, typename Number>
FEPointEvaluation<n_components, dim, spacedim, Number>::FEPointEvaluation(
  const Mapping<dim> &      mapping,
  const FiniteElement<dim> &fe,
  const UpdateFlags         update_flags,
  const unsigned int        first_selected_component)
  : n_q_points(numbers::invalid_unsigned_int)
  , n_filled_lanes_last_batch(numbers::invalid_unsigned_int)
  , mapping(&mapping)
  , fe(&fe)
  , update_flags(update_flags)
  , mapping_info_on_the_fly(
      std::make_unique<NonMatching::MappingInfo<dim, spacedim>>(mapping,
                                                                update_flags))
  , mapping_info(mapping_info_on_the_fly.get())
  , current_cell_index(numbers::invalid_unsigned_int)
  , current_face_number(numbers::invalid_unsigned_int)
  , is_reinitialized(false)
{
  setup(first_selected_component);
}



template <int n_components, int dim, int spacedim, typename Number>
FEPointEvaluation<n_components, dim, spacedim, Number>::FEPointEvaluation(
  NonMatching::MappingInfo<dim, spacedim> &mapping_info,
  const FiniteElement<dim> &               fe,
  const unsigned int                       first_selected_component)
  : n_q_points(numbers::invalid_unsigned_int)
  , n_filled_lanes_last_batch(numbers::invalid_unsigned_int)
  , mapping(&mapping_info.get_mapping())
  , fe(&fe)
  , update_flags(mapping_info.get_update_flags())
  , mapping_info(&mapping_info)
  , current_cell_index(numbers::invalid_unsigned_int)
  , current_face_number(numbers::invalid_unsigned_int)
  , is_reinitialized(false)
{
  setup(first_selected_component);
  connection_is_reinitialized = mapping_info.connect_is_reinitialized(
    [this]() { this->is_reinitialized = false; });
}



template <int n_components_, int dim, int spacedim, typename Number>
FEPointEvaluation<n_components_, dim, spacedim, Number>::FEPointEvaluation(
  FEPointEvaluation<n_components_, dim, spacedim, Number> &other) noexcept
  : n_q_points(other.n_q_points)
  , n_filled_lanes_last_batch(numbers::invalid_unsigned_int)
  , mapping(other.mapping)
  , fe(other.fe)
  , poly(other.poly)
  , polynomials_are_hat_functions(other.polynomials_are_hat_functions)
  , renumber(other.renumber)
  , solution_renumbered(other.solution_renumbered)
  , solution_renumbered_vectorized(other.solution_renumbered_vectorized)
  , values(other.values)
  , unit_gradients(other.unit_gradients)
  , gradients(other.gradients)
  , dofs_per_component(other.dofs_per_component)
  , component_in_base_element(other.component_in_base_element)
  , nonzero_shape_function_component(other.nonzero_shape_function_component)
  , update_flags(other.update_flags)
  , fe_values(other.fe_values)
  , mapping_info_on_the_fly(
      other.mapping_info_on_the_fly ?
        std::make_unique<NonMatching::MappingInfo<dim, spacedim>>(
          *mapping,
          update_flags) :
        nullptr)
  , mapping_info(other.mapping_info)
  , current_cell_index(other.current_cell_index)
  , current_face_number(other.current_face_number)
  , fast_path(other.fast_path)
  , is_reinitialized(false)
  , shapes(other.shapes)
{
  connection_is_reinitialized = mapping_info->connect_is_reinitialized(
    [this]() { this->is_reinitialized = false; });
}



template <int n_components_, int dim, int spacedim, typename Number>
FEPointEvaluation<n_components_, dim, spacedim, Number>::FEPointEvaluation(
  FEPointEvaluation<n_components_, dim, spacedim, Number> &&other) noexcept
  : n_q_points(other.n_q_points)
  , n_filled_lanes_last_batch(numbers::invalid_unsigned_int)
  , mapping(other.mapping)
  , fe(other.fe)
  , poly(other.poly)
  , polynomials_are_hat_functions(other.polynomials_are_hat_functions)
  , renumber(other.renumber)
  , solution_renumbered(other.solution_renumbered)
  , solution_renumbered_vectorized(other.solution_renumbered_vectorized)
  , values(other.values)
  , unit_gradients(other.unit_gradients)
  , gradients(other.gradients)
  , dofs_per_component(other.dofs_per_component)
  , component_in_base_element(other.component_in_base_element)
  , nonzero_shape_function_component(other.nonzero_shape_function_component)
  , update_flags(other.update_flags)
  , fe_values(other.fe_values)
  , mapping_info_on_the_fly(std::move(other.mapping_info_on_the_fly))
  , mapping_info(other.mapping_info)
  , current_cell_index(other.current_cell_index)
  , current_face_number(other.current_face_number)
  , fast_path(other.fast_path)
  , is_reinitialized(false)
  , shapes(other.shapes)
{
  connection_is_reinitialized = mapping_info->connect_is_reinitialized(
    [this]() { this->is_reinitialized = false; });
}



template <int n_components, int dim, int spacedim, typename Number>
FEPointEvaluation<n_components, dim, spacedim, Number>::~FEPointEvaluation()
{
  connection_is_reinitialized.disconnect();
}



template <int n_components, int dim, int spacedim, typename Number>
void
FEPointEvaluation<n_components, dim, spacedim, Number>::setup(
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
      internal::MatrixFreeFunctions::ShapeInfo<double> shape_info;

      shape_info.reinit(QMidpoint<1>(), *fe, base_element_number);
      renumber           = shape_info.lexicographic_numbering;
      dofs_per_component = shape_info.dofs_per_component_on_cell;
      poly               = internal::FEPointEvaluation::get_polynomial_space(
        fe->base_element(base_element_number));

      polynomials_are_hat_functions =
        (poly.size() == 2 && poly[0].value(0.) == 1. &&
         poly[0].value(1.) == 0. && poly[1].value(0.) == 0. &&
         poly[1].value(1.) == 1.);

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



template <int n_components, int dim, int spacedim, typename Number>
void
FEPointEvaluation<n_components, dim, spacedim, Number>::reinit(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const ArrayView<const Point<dim>> &                         unit_points)
{
  // reinit is only allowed for mapping computation on the fly
  AssertThrow(mapping_info_on_the_fly.get() != nullptr, ExcNotImplemented());

  mapping_info->reinit(cell, unit_points);

  if (!fast_path)
    {
      fe_values = std::make_shared<FEValues<dim, spacedim>>(
        *mapping,
        *fe,
        Quadrature<dim>(
          std::vector<Point<dim>>(unit_points.begin(), unit_points.end())),
        update_flags);
      fe_values->reinit(cell);
    }

  do_reinit();
}



template <int n_components, int dim, int spacedim, typename Number>
void
FEPointEvaluation<n_components, dim, spacedim, Number>::reinit()
{
  current_cell_index  = numbers::invalid_unsigned_int;
  current_face_number = numbers::invalid_unsigned_int;

  do_reinit();
}



template <int n_components, int dim, int spacedim, typename Number>
void
FEPointEvaluation<n_components, dim, spacedim, Number>::reinit(
  const unsigned int cell_index)
{
  current_cell_index  = cell_index;
  current_face_number = numbers::invalid_unsigned_int;

  do_reinit();
}



template <int n_components, int dim, int spacedim, typename Number>
void
FEPointEvaluation<n_components, dim, spacedim, Number>::reinit(
  const unsigned int cell_index,
  const unsigned int face_number)
{
  current_cell_index  = cell_index;
  current_face_number = face_number;

  do_reinit();
}



template <int n_components, int dim, int spacedim, typename Number>
void
FEPointEvaluation<n_components, dim, spacedim, Number>::do_reinit()
{
  const auto unit_points =
    mapping_info->get_unit_points(current_cell_index, current_face_number);

  const unsigned int n_q_points_unvectorized = unit_points.size();

  if (std::is_same<ScalarNumber, Number>::value)
    {
      const_cast<unsigned int &>(n_q_points) = unit_points.size();
    }
  else
    {
      const unsigned int n_lanes =
        internal::VectorizedArrayTrait<VectorizedArrayType>::width;
      const_cast<unsigned int &>(n_filled_lanes_last_batch) =
        n_q_points_unvectorized % n_lanes;
      const_cast<unsigned int &>(n_q_points) =
        n_q_points_unvectorized / n_lanes;
      if (n_filled_lanes_last_batch > 0)
        ++const_cast<unsigned int &>(n_q_points);
    }

  if (update_flags & update_values)
    values.resize(n_q_points, numbers::signaling_nan<value_type>());
  if (update_flags & update_gradients)
    gradients.resize(n_q_points, numbers::signaling_nan<gradient_type>());

  if (!polynomials_are_hat_functions)
    {
      const std::size_t n_points = unit_points.size();
      const std::size_t n_lanes  = VectorizedArrayType::size();
      const std::size_t n_batches =
        n_points / n_lanes + (n_points % n_lanes > 0 ? 1 : 0);
      const std::size_t n_shapes = poly.size();
      shapes.resize_fast(n_batches * n_shapes);
      for (unsigned int i = 0, qb = 0; i < n_points; i += n_lanes, ++qb)
        {
          // convert to vectorized format
          Point<dim, VectorizedArrayType> vectorized_points;
          for (unsigned int j = 0; j < n_lanes && i + j < n_points; ++j)
            for (unsigned int d = 0; d < dim; ++d)
              vectorized_points[d][j] = unit_points[i + j][d];

          auto view =
            make_array_view(shapes.begin() + qb * n_shapes,
                            shapes.begin() + (qb * n_shapes + n_shapes));

          internal::compute_values_of_array(view, poly, vectorized_points);
        }
    }

  is_reinitialized = true;
}



template <int n_components, int dim, int spacedim, typename Number>
inline void
FEPointEvaluation<n_components, dim, spacedim, Number>::evaluate_fast(
  const ArrayView<const ScalarNumber> &   solution_values,
  const EvaluationFlags::EvaluationFlags &evaluation_flags)
{
  // fast path with tensor product evaluation
  if (solution_renumbered.size() != dofs_per_component)
    solution_renumbered.resize(dofs_per_component);
  for (unsigned int comp = 0; comp < n_components; ++comp)
    for (unsigned int i = 0; i < dofs_per_component; ++i)
      internal::FEPointEvaluation::
        EvaluatorTypeTraits<dim, n_components, ScalarNumber>::read_value(
          solution_values[renumber[(component_in_base_element + comp) *
                                     dofs_per_component +
                                   i]],
          comp,
          solution_renumbered[i]);

  // unit gradients are currently only implemented with the fast tensor
  // path
  unit_gradients.resize(n_q_points, numbers::signaling_nan<gradient_type>());

  const auto unit_points =
    mapping_info->get_unit_points(current_cell_index, current_face_number);
  const auto &mapping_data =
    mapping_info->get_mapping_data(current_cell_index, current_face_number);

  const std::size_t n_points = unit_points.size();
  const std::size_t n_lanes =
    internal::VectorizedArrayTrait<VectorizedArrayType>::width;
  const std::size_t stride =
    internal::FEPointEvaluation::VectorizedPointsTypeTraits<Number>::stride;

  // loop over quadrature batches qb / points q
  for (unsigned int qb = 0, q = 0; q < n_points; ++qb, q += n_lanes)
    {
      // compute
      const unsigned int n_shapes     = poly.size();
      const auto         val_and_grad = [&]() {
        if (polynomials_are_hat_functions)
          {
            // convert quadrature points to vectorized format
            Point<dim, VectorizedArrayType> vectorized_points;
            for (unsigned int v = 0; v < n_lanes && q + v < n_points; ++v)
              for (unsigned int d = 0; d < dim; ++d)
                vectorized_points[d][v] = unit_points[q + v][d];

            return internal::evaluate_tensor_product_value_and_gradient_linear(
              poly, solution_renumbered, vectorized_points);
          }
        else
          return internal::evaluate_tensor_product_value_and_gradient_shapes<
            dim,
            scalar_value_type,
            VectorizedArrayType>(make_array_view(shapes.begin() + qb * n_shapes,
                                                 shapes.begin() +
                                                   (qb * n_shapes + n_shapes)),
                                 poly.size(),
                                 solution_renumbered);
      }();

      if (evaluation_flags & EvaluationFlags::values)
        {
          for (unsigned int v = 0; v < stride && q + v < n_points; ++v)
            internal::FEPointEvaluation::
              EvaluatorTypeTraits<dim, n_components, Number>::set_value(
                val_and_grad.first, v, values[qb * stride + v]);
        }
      if (evaluation_flags & EvaluationFlags::gradients)
        {
          Assert(update_flags & update_gradients ||
                   update_flags & update_inverse_jacobians,
                 ExcNotInitialized());
          if (std::is_same<Number, ScalarNumber>::value)
            {
              // convert back to standard format
              for (unsigned int v = 0; v < n_lanes && q + v < n_points; ++v)
                {
                  internal::FEPointEvaluation::EvaluatorTypeTraits<
                    dim,
                    n_components,
                    Number>::set_gradient(val_and_grad.second,
                                          v,
                                          unit_gradients[q + v]);
                  gradients[q + v] = apply_transformation(
                    mapping_data.inverse_jacobians[q + v].transpose(),
                    unit_gradients[q + v]);
                }
            }
          else
            {
              // convert to vectorized format
              DerivativeForm<1, spacedim, dim, Number>
                vectorized_inverse_transposed_jacobians;
              for (unsigned int v = 0; v < n_lanes && q + v < n_points; ++v)
                for (unsigned int d = 0; d < dim; ++d)
                  for (unsigned int s = 0; s < spacedim; ++s)
                    internal::VectorizedArrayTrait<Number>::get(
                      vectorized_inverse_transposed_jacobians[s][d], v) =
                      mapping_data.inverse_jacobians[q + v].transpose()[s][d];

              internal::FEPointEvaluation::
                EvaluatorTypeTraits<dim, n_components, Number>::set_gradient(
                  val_and_grad.second, 0, unit_gradients[qb]);
              gradients[qb] =
                apply_transformation(vectorized_inverse_transposed_jacobians,
                                     unit_gradients[qb]);
            }
        }
    }
}



template <int n_components, int dim, int spacedim, typename Number>
inline void
FEPointEvaluation<n_components, dim, spacedim, Number>::evaluate_slow(
  const ArrayView<const ScalarNumber> &   solution_values,
  const EvaluationFlags::EvaluationFlags &evaluation_flags)
{
  // slow path with FEValues
  Assert(fe_values.get() != nullptr,
         ExcMessage(
           "Not initialized. Please call FEPointEvaluation::reinit()!"));

  const std::size_t n_points = fe_values->get_quadrature().size();
  const std::size_t n_lanes  = internal::VectorizedArrayTrait<Number>::width;

  if (evaluation_flags & EvaluationFlags::values)
    {
      values.resize(n_q_points);
      std::fill(values.begin(), values.end(), value_type());
      for (unsigned int i = 0; i < fe->n_dofs_per_cell(); ++i)
        {
          const ScalarNumber value = solution_values[i];
          for (unsigned int d = 0; d < n_components; ++d)
            if (nonzero_shape_function_component[i][d] &&
                (fe->is_primitive(i) || fe->is_primitive()))
              for (unsigned int qb = 0, q = 0; q < n_points; ++qb, q += n_lanes)
                for (unsigned int v = 0;
                     v < (q + n_lanes > n_points ? n_filled_lanes_last_batch :
                                                   n_lanes);
                     ++v)
                  internal::FEPointEvaluation::
                    EvaluatorTypeTraits<dim, n_components, Number>::access(
                      values[qb],
                      v,
                      d,
                      fe_values->shape_value(i, q + v) * value);
            else if (nonzero_shape_function_component[i][d])
              for (unsigned int qb = 0, q = 0; q < n_points; ++qb, q += n_lanes)
                for (unsigned int v = 0;
                     v < (q + n_lanes > n_points ? n_filled_lanes_last_batch :
                                                   n_lanes);
                     ++v)
                  internal::FEPointEvaluation::
                    EvaluatorTypeTraits<dim, n_components, Number>::access(
                      values[qb],
                      v,
                      d,
                      fe_values->shape_value_component(i, q + v, d) * value);
        }
    }

  if (evaluation_flags & EvaluationFlags::gradients)
    {
      gradients.resize(n_q_points);
      std::fill(gradients.begin(), gradients.end(), gradient_type());
      for (unsigned int i = 0; i < fe->n_dofs_per_cell(); ++i)
        {
          const ScalarNumber value = solution_values[i];
          for (unsigned int d = 0; d < n_components; ++d)
            if (nonzero_shape_function_component[i][d] &&
                (fe->is_primitive(i) || fe->is_primitive()))
              for (unsigned int qb = 0, q = 0; q < n_points; ++qb, q += n_lanes)
                for (unsigned int v = 0;
                     v < (q + n_lanes > n_points ? n_filled_lanes_last_batch :
                                                   n_lanes);
                     ++v)
                  internal::FEPointEvaluation::
                    EvaluatorTypeTraits<dim, n_components, Number>::access(
                      gradients[qb],
                      v,
                      d,
                      fe_values->shape_grad(i, q + v) * value);
            else if (nonzero_shape_function_component[i][d])
              for (unsigned int qb = 0, q = 0; q < n_points; ++qb, q += n_lanes)
                for (unsigned int v = 0;
                     v < (q + n_lanes > n_points ? n_filled_lanes_last_batch :
                                                   n_lanes);
                     ++v)
                  internal::FEPointEvaluation::
                    EvaluatorTypeTraits<dim, n_components, Number>::access(
                      gradients[qb],
                      v,
                      d,
                      fe_values->shape_grad_component(i, q + v, d) * value);
        }
    }
}



template <int n_components, int dim, int spacedim, typename Number>
void
FEPointEvaluation<n_components, dim, spacedim, Number>::evaluate(
  const ArrayView<const ScalarNumber> &   solution_values,
  const EvaluationFlags::EvaluationFlags &evaluation_flags)
{
  if (!is_reinitialized)
    reinit();

  if (n_q_points == 0)
    return;

  Assert(!(evaluation_flags & EvaluationFlags::hessians), ExcNotImplemented());

  if (!((evaluation_flags & EvaluationFlags::values) ||
        (evaluation_flags & EvaluationFlags::gradients))) // no evaluation flags
    return;

  AssertDimension(solution_values.size(), fe->dofs_per_cell);
  if (fast_path)
    evaluate_fast(solution_values, evaluation_flags);
  else
    evaluate_slow(solution_values, evaluation_flags);
}



template <int n_components, int dim, int spacedim, typename Number>
inline void
FEPointEvaluation<n_components, dim, spacedim, Number>::integrate_fast(
  const ArrayView<ScalarNumber> &         solution_values,
  const EvaluationFlags::EvaluationFlags &integration_flags)

{
  // fast path with tensor product integration
  if (solution_renumbered_vectorized.size() != dofs_per_component)
    solution_renumbered_vectorized.resize(dofs_per_component);
  // zero content
  solution_renumbered_vectorized.fill(
    typename internal::FEPointEvaluation::EvaluatorTypeTraits<
      dim,
      n_components,
      VectorizedArrayType>::value_type());

  const auto unit_points =
    mapping_info->get_unit_points(current_cell_index, current_face_number);
  const auto &mapping_data =
    mapping_info->get_mapping_data(current_cell_index, current_face_number);

  const std::size_t n_points = unit_points.size();
  const std::size_t n_lanes =
    internal::VectorizedArrayTrait<VectorizedArrayType>::width;
  const std::size_t stride =
    internal::FEPointEvaluation::VectorizedPointsTypeTraits<Number>::stride;

  // loop over quadrature batches qb / points q
  for (unsigned int qb = 0, q = 0; q < n_points; ++qb, q += n_lanes)
    {
      const bool incomplete_last_batch =
        (qb == (n_q_points - 1)) && (n_filled_lanes_last_batch > 0);

      typename internal::ProductTypeNoPoint<value_type,
                                            VectorizedArrayType>::type value =
        {};
      Tensor<1,
             dim,
             typename internal::ProductTypeNoPoint<value_type,
                                                   VectorizedArrayType>::type>
        gradient;

      if (integration_flags & EvaluationFlags::values)
        {
          // zero out lanes of incomplete last quadrature point batch
          if (incomplete_last_batch)
            {
              for (unsigned int v = n_filled_lanes_last_batch; v < n_lanes; ++v)
                internal::FEPointEvaluation::EvaluatorTypeTraits<
                  dim,
                  n_components,
                  Number>::set_zero_value(values[qb], v);
            }
          for (unsigned int v = 0; v < stride && q + v < n_points; ++v)
            internal::FEPointEvaluation::
              EvaluatorTypeTraits<dim, n_components, Number>::get_value(
                value, v, values[qb * stride + v]);
        }
      if (integration_flags & EvaluationFlags::gradients)
        {
          if (std::is_same<Number, ScalarNumber>::value)
            {
              for (unsigned int v = 0; v < n_lanes && q + v < n_points; ++v)
                {
                  gradients[q + v] =
                    apply_transformation(mapping_data.inverse_jacobians[q + v],
                                         gradients[q + v]);
                  internal::FEPointEvaluation::EvaluatorTypeTraits<
                    dim,
                    n_components,
                    Number>::get_gradient(gradient, v, gradients[q + v]);
                }
            }
          else
            {
              // convert to vectorized format
              DerivativeForm<1, spacedim, dim, Number>
                vectorized_inverse_jacobians;
              for (unsigned int v = 0; v < n_lanes && q + v < n_points; ++v)
                for (unsigned int d = 0; d < dim; ++d)
                  for (unsigned int s = 0; s < spacedim; ++s)
                    internal::VectorizedArrayTrait<Number>::get(
                      vectorized_inverse_jacobians[s][d], v) =
                      mapping_data.inverse_jacobians[q + v][s][d];

              // zero out lanes of incomplete last quadrature point batch
              if (incomplete_last_batch)
                {
                  for (unsigned int v = n_filled_lanes_last_batch; v < n_lanes;
                       ++v)
                    internal::FEPointEvaluation::EvaluatorTypeTraits<
                      dim,
                      n_components,
                      Number>::set_zero_gradient(gradients[qb], v);
                }

              gradients[qb] = apply_transformation(vectorized_inverse_jacobians,
                                                   gradients[qb]);
              internal::FEPointEvaluation::
                EvaluatorTypeTraits<dim, n_components, Number>::get_gradient(
                  gradient, 0, gradients[qb]);
            }
        }

      // compute
      const unsigned int n_shapes = poly.size();
      if (polynomials_are_hat_functions)
        {
          // convert to vectorized format
          Point<dim, VectorizedArrayType> vectorized_points;
          for (unsigned int j = 0; j < n_lanes && q + j < n_points; ++j)
            for (unsigned int d = 0; d < dim; ++d)
              vectorized_points[d][j] = unit_points[q + j][d];

          internal::integrate_add_tensor_product_value_and_gradient_linear(
            poly,
            value,
            gradient,
            solution_renumbered_vectorized,
            vectorized_points);
        }
      else
        internal::integrate_add_tensor_product_value_and_gradient_shapes<
          dim,
          VectorizedArrayType,
          typename internal::ProductTypeNoPoint<value_type,
                                                VectorizedArrayType>::type>(
          make_array_view(shapes.begin() + qb * n_shapes,
                          shapes.begin() + (qb * n_shapes + n_shapes)),
          n_shapes,
          value,
          gradient,
          solution_renumbered_vectorized);
    }

  // add between the lanes and write into the result
  std::fill(solution_values.begin(), solution_values.end(), ScalarNumber());
  for (unsigned int comp = 0; comp < n_components; ++comp)
    for (unsigned int i = 0; i < dofs_per_component; ++i)
      {
        VectorizedArrayType result;
        internal::FEPointEvaluation::EvaluatorTypeTraits<
          dim,
          n_components,
          VectorizedArrayType>::write_value(result,
                                            comp,
                                            solution_renumbered_vectorized[i]);
        for (unsigned int lane = n_lanes / 2; lane > 0; lane /= 2)
          for (unsigned int j = 0; j < lane; ++j)
            result[j] += result[lane + j];
        solution_values[renumber[comp * dofs_per_component + i]] = result[0];
      }
}



template <int n_components, int dim, int spacedim, typename Number>
inline void
FEPointEvaluation<n_components, dim, spacedim, Number>::integrate_slow(
  const ArrayView<ScalarNumber> &         solution_values,
  const EvaluationFlags::EvaluationFlags &integration_flags)
{
  // slow path with FEValues
  Assert(fe_values.get() != nullptr,
         ExcMessage(
           "Not initialized. Please call FEPointEvaluation::reinit()!"));
  std::fill(solution_values.begin(), solution_values.end(), 0.0);

  const std::size_t n_points = fe_values->get_quadrature().size();
  const std::size_t n_lanes  = internal::VectorizedArrayTrait<Number>::width;

  if (integration_flags & EvaluationFlags::values)
    {
      AssertIndexRange(n_q_points, values.size() + 1);
      for (unsigned int i = 0; i < fe->n_dofs_per_cell(); ++i)
        {
          for (unsigned int d = 0; d < n_components; ++d)
            if (nonzero_shape_function_component[i][d] &&
                (fe->is_primitive(i) || fe->is_primitive()))
              for (unsigned int qb = 0, q = 0; q < n_points; ++qb, q += n_lanes)
                for (unsigned int v = 0;
                     v < (q + n_lanes > n_points ? n_filled_lanes_last_batch :
                                                   n_lanes);
                     ++v)
                  solution_values[i] +=
                    fe_values->shape_value(i, q + v) *
                    internal::FEPointEvaluation::
                      EvaluatorTypeTraits<dim, n_components, Number>::access(
                        values[qb], v, d);
            else if (nonzero_shape_function_component[i][d])
              for (unsigned int qb = 0, q = 0; q < n_points; ++qb, q += n_lanes)
                for (unsigned int v = 0;
                     v < (q + n_lanes > n_points ? n_filled_lanes_last_batch :
                                                   n_lanes);
                     ++v)
                  solution_values[i] +=
                    fe_values->shape_value_component(i, q + v, d) *
                    internal::FEPointEvaluation::
                      EvaluatorTypeTraits<dim, n_components, Number>::access(
                        values[qb], v, d);
        }
    }

  if (integration_flags & EvaluationFlags::gradients)
    {
      AssertIndexRange(n_q_points, gradients.size() + 1);
      for (unsigned int i = 0; i < fe->n_dofs_per_cell(); ++i)
        {
          for (unsigned int d = 0; d < n_components; ++d)
            if (nonzero_shape_function_component[i][d] &&
                (fe->is_primitive(i) || fe->is_primitive()))
              for (unsigned int qb = 0, q = 0; q < n_points; ++qb, q += n_lanes)
                for (unsigned int v = 0;
                     v < (q + n_lanes > n_points ? n_filled_lanes_last_batch :
                                                   n_lanes);
                     ++v)
                  solution_values[i] +=
                    fe_values->shape_grad(i, q + v) *
                    internal::FEPointEvaluation::
                      EvaluatorTypeTraits<dim, n_components, Number>::access(
                        gradients[qb], v, d);
            else if (nonzero_shape_function_component[i][d])
              for (unsigned int qb = 0, q = 0; q < n_points; ++qb, q += n_lanes)
                for (unsigned int v = 0;
                     v < (q + n_lanes > n_points ? n_filled_lanes_last_batch :
                                                   n_lanes);
                     ++v)
                  solution_values[i] +=
                    fe_values->shape_grad_component(i, q + v, d) *
                    internal::FEPointEvaluation::
                      EvaluatorTypeTraits<dim, n_components, Number>::access(
                        gradients[qb], v, d);
        }
    }
}



template <int n_components, int dim, int spacedim, typename Number>
void
FEPointEvaluation<n_components, dim, spacedim, Number>::integrate(
  const ArrayView<ScalarNumber> &         solution_values,
  const EvaluationFlags::EvaluationFlags &integration_flags)
{
  if (!is_reinitialized)
    reinit();

  if (n_q_points == 0) // no evaluation points provided
    {
      std::fill(solution_values.begin(), solution_values.end(), 0.0);
      return;
    }

  Assert(!(integration_flags & EvaluationFlags::hessians), ExcNotImplemented());

  if (!((integration_flags & EvaluationFlags::values) ||
        (integration_flags &
         EvaluationFlags::gradients))) // no integration flags
    {
      std::fill(solution_values.begin(), solution_values.end(), 0.0);
      return;
    }

  AssertDimension(solution_values.size(), fe->dofs_per_cell);
  if (fast_path)
    integrate_fast(solution_values, integration_flags);
  else
    integrate_slow(solution_values, integration_flags);
}



template <int n_components, int dim, int spacedim, typename Number>
inline const typename FEPointEvaluation<n_components, dim, spacedim, Number>::
  value_type &
  FEPointEvaluation<n_components, dim, spacedim, Number>::get_value(
    const unsigned int point_index) const
{
  AssertIndexRange(point_index, values.size());
  return values[point_index];
}



template <int n_components, int dim, int spacedim, typename Number>
inline const typename FEPointEvaluation<n_components, dim, spacedim, Number>::
  gradient_type &
  FEPointEvaluation<n_components, dim, spacedim, Number>::get_gradient(
    const unsigned int point_index) const
{
  AssertIndexRange(point_index, gradients.size());
  return gradients[point_index];
}



template <int n_components, int dim, int spacedim, typename Number>
inline const typename FEPointEvaluation<n_components, dim, spacedim, Number>::
  gradient_type &
  FEPointEvaluation<n_components, dim, spacedim, Number>::get_unit_gradient(
    const unsigned int point_index) const
{
  Assert(fast_path,
         ExcMessage("Unit gradients are currently only implemented for tensor "
                    "product finite elements combined with MappingQ "
                    "mappings"));
  AssertIndexRange(point_index, unit_gradients.size());
  return unit_gradients[point_index];
}



template <int n_components, int dim, int spacedim, typename Number>
inline void
FEPointEvaluation<n_components, dim, spacedim, Number>::submit_value(
  const value_type & value,
  const unsigned int point_index)
{
  AssertIndexRange(point_index, n_q_points);
  values[point_index] = value;
}



template <int n_components, int dim, int spacedim, typename Number>
inline void
FEPointEvaluation<n_components, dim, spacedim, Number>::submit_gradient(
  const gradient_type &gradient,
  const unsigned int   point_index)
{
  AssertIndexRange(point_index, n_q_points);
  gradients[point_index] = gradient;
}



template <int n_components, int dim, int spacedim, typename Number>
inline DerivativeForm<1, dim, spacedim, Number>
FEPointEvaluation<n_components, dim, spacedim, Number>::jacobian(
  const unsigned int point_index) const
{
  const auto &mapping_data =
    mapping_info->get_mapping_data(current_cell_index, current_face_number);

  const unsigned int n_lanes = internal::VectorizedArrayTrait<Number>::width;
  DerivativeForm<1, dim, spacedim, Number> vectorized_jacobian;
  for (unsigned int v = 0;
       v < n_lanes && point_index * n_lanes + v < mapping_data.jacobians.size();
       ++v)
    for (unsigned int d = 0; d < dim; ++d)
      for (unsigned int s = 0; s < spacedim; ++s)
        internal::VectorizedArrayTrait<Number>::get(vectorized_jacobian[d][s],
                                                    v) =
          mapping_data.jacobians[point_index * n_lanes + v][d][s];
  return vectorized_jacobian;
}



template <int n_components, int dim, int spacedim, typename Number>
inline DerivativeForm<1, spacedim, dim, Number>
FEPointEvaluation<n_components, dim, spacedim, Number>::inverse_jacobian(
  const unsigned int point_index) const
{
  const auto &mapping_data =
    mapping_info->get_mapping_data(current_cell_index, current_face_number);

  const unsigned int n_lanes = internal::VectorizedArrayTrait<Number>::width;
  DerivativeForm<1, dim, spacedim, Number> vectorized_inverse_jacobian;
  for (unsigned int v = 0;
       v < n_lanes &&
       point_index * n_lanes + v < mapping_data.inverse_jacobians.size();
       ++v)
    for (unsigned int d = 0; d < dim; ++d)
      for (unsigned int s = 0; s < spacedim; ++s)
        internal::VectorizedArrayTrait<Number>::get(
          vectorized_inverse_jacobian[d][s], v) =
          mapping_data.inverse_jacobians[point_index * n_lanes + v][d][s];
  return vectorized_inverse_jacobian;
}



template <int n_components, int dim, int spacedim, typename Number>
inline Number
FEPointEvaluation<n_components, dim, spacedim, Number>::JxW(
  const unsigned int point_index) const
{
  const auto &mapping_data =
    mapping_info->get_mapping_data(current_cell_index, current_face_number);

  const unsigned int n_lanes = internal::VectorizedArrayTrait<Number>::width;
  Number             vectorized_JxW = 0;
  for (unsigned int v = 0; v < n_lanes && point_index * n_lanes + v <
                                            mapping_data.JxW_values.size();
       ++v)
    internal::VectorizedArrayTrait<Number>::get(vectorized_JxW, v) =
      mapping_data.JxW_values[point_index * n_lanes + v];
  return vectorized_JxW;
}


template <int n_components, int dim, int spacedim, typename Number>
inline Tensor<1, spacedim, Number>
FEPointEvaluation<n_components, dim, spacedim, Number>::normal_vector(
  const unsigned int point_index) const
{
  const auto &mapping_data =
    mapping_info->get_mapping_data(current_cell_index, current_face_number);

  const unsigned int n_lanes = internal::VectorizedArrayTrait<Number>::width;
  Tensor<1, spacedim, Number> vectorized_normal_vectors;
  for (unsigned int v = 0; v < n_lanes && point_index * n_lanes + v <
                                            mapping_data.normal_vectors.size();
       ++v)
    for (unsigned int s = 0; s < spacedim; ++s)
      internal::VectorizedArrayTrait<Number>::get(vectorized_normal_vectors[s],
                                                  v) =
        mapping_data.normal_vectors[point_index * n_lanes + v][s];
  return vectorized_normal_vectors;
}



template <int n_components, int dim, int spacedim, typename Number>
inline Point<spacedim, Number>
FEPointEvaluation<n_components, dim, spacedim, Number>::real_point(
  const unsigned int point_index) const
{
  const auto &mapping_data =
    mapping_info->get_mapping_data(current_cell_index, current_face_number);

  const unsigned int n_lanes = internal::VectorizedArrayTrait<Number>::width;
  Point<spacedim, Number> vectorized_real_point;
  for (unsigned int v = 0;
       v < n_lanes &&
       point_index * n_lanes + v < mapping_data.quadrature_points.size();
       ++v)
    for (unsigned int s = 0; s < spacedim; ++s)
      internal::VectorizedArrayTrait<Number>::get(vectorized_real_point[s], v) =
        mapping_data.quadrature_points[point_index * n_lanes + v][s];
  return vectorized_real_point;
}



template <int n_components, int dim, int spacedim, typename Number>
inline Point<dim, Number>
FEPointEvaluation<n_components, dim, spacedim, Number>::unit_point(
  const unsigned int point_index) const
{
  const auto unit_points =
    mapping_info->get_unit_points(current_cell_index, current_face_number);
  const unsigned int n_lanes = internal::VectorizedArrayTrait<Number>::width;
  Point<spacedim, Number> vectorized_unit_point;
  for (unsigned int v = 0;
       v < n_lanes && point_index * n_lanes + v < unit_points.size();
       ++v)
    for (unsigned int d = 0; d < dim; ++d)
      internal::VectorizedArrayTrait<Number>::get(vectorized_unit_point[d], v) =
        unit_points[point_index * n_lanes + v][d];
  return vectorized_unit_point;
}



template <int n_components, int dim, int spacedim, typename Number>
inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
FEPointEvaluation<n_components, dim, spacedim, Number>::
  quadrature_point_indices() const
{
  return {0U, n_q_points};
}

DEAL_II_NAMESPACE_CLOSE

#endif
