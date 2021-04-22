/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2020 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 */

#ifndef dealii_fe_point_evaluation_h
#define dealii_fe_point_evaluation_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/matrix_free/evaluation_flags.h>
#include <deal.II/matrix_free/shape_info.h>
#include <deal.II/matrix_free/tensor_product_kernels.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace FEPointEvaluation
  {
    /**
     * Struct to distinguish between the value and gradient types of different
     * numbers of components used by the FlexibleEvaluator class.
     */
    template <int dim, int n_components>
    struct EvaluatorTypeTraits
    {
      using value_type    = Tensor<1, n_components>;
      using gradient_type = Tensor<1, n_components, Tensor<1, dim>>;

      static void
      read_value(const double       vector_entry,
                 const unsigned int component,
                 value_type &       result)
      {
        AssertIndexRange(component, n_components);
        result[component] = vector_entry;
      }

      static void
      write_value(double &           vector_entry,
                  const unsigned int component,
                  const value_type & result)
      {
        AssertIndexRange(component, n_components);
        vector_entry = result[component];
      }

      static void
      set_gradient(
        const Tensor<1, dim, Tensor<1, n_components, VectorizedArray<double>>>
          &                value,
        const unsigned int vector_lane,
        gradient_type &    result)
      {
        for (unsigned int i = 0; i < n_components; ++i)
          for (unsigned int d = 0; d < dim; ++d)
            result[i][d] = value[d][i][vector_lane];
      }

      static void get_gradient(
        Tensor<1, dim, Tensor<1, n_components, VectorizedArray<double>>> &value,
        const unsigned int   vector_lane,
        const gradient_type &result)
      {
        for (unsigned int i = 0; i < n_components; ++i)
          for (unsigned int d = 0; d < dim; ++d)
            value[d][i][vector_lane] = result[i][d];
      }

      static void
      set_value(const Tensor<1, n_components, VectorizedArray<double>> &value,
                const unsigned int vector_lane,
                value_type &       result)
      {
        for (unsigned int i = 0; i < n_components; ++i)
          result[i] = value[i][vector_lane];
      }

      static void
        get_value(Tensor<1, n_components, VectorizedArray<double>> &value,
                  const unsigned int                                vector_lane,
                  const value_type &                                result)
      {
        for (unsigned int i = 0; i < n_components; ++i)
          value[i][vector_lane] = result[i];
      }

      template <typename Number>
      static Number &access(Tensor<1, n_components, Number> &value,
                            const unsigned int               component)
      {
        return value[component];
      }

      template <typename Number>
      static const Number &
      access(const Tensor<1, n_components, Number> &value,
             const unsigned int                     component)
      {
        return value[component];
      }
    };

    template <int dim>
    struct EvaluatorTypeTraits<dim, 1>
    {
      using value_type    = double;
      using gradient_type = Tensor<1, dim>;

      static void
      read_value(const double vector_entry,
                 const unsigned int,
                 value_type &result)
      {
        result = vector_entry;
      }

      static void
      write_value(double &vector_entry,
                  const unsigned int,
                  const value_type &result)
      {
        vector_entry = result;
      }

      static void
      set_gradient(const Tensor<1, dim, VectorizedArray<double>> &value,
                   const unsigned int                             vector_lane,
                   gradient_type &                                result)
      {
        for (unsigned int d = 0; d < dim; ++d)
          result[d] = value[d][vector_lane];
      }

      static void get_gradient(Tensor<1, dim, VectorizedArray<double>> &value,
                               const unsigned int   vector_lane,
                               const gradient_type &result)
      {
        for (unsigned int d = 0; d < dim; ++d)
          value[d][vector_lane] = result[d];
      }

      static void
      set_value(const VectorizedArray<double> &value,
                const unsigned int             vector_lane,
                value_type &                   result)
      {
        result = value[vector_lane];
      }

      static void
      get_value(VectorizedArray<double> &value,
                const unsigned int       vector_lane,
                const value_type &       result)
      {
        value[vector_lane] = result;
      }

      template <typename Number>
      static Number &
      access(Number &value, const unsigned int)
      {
        return value;
      }

      template <typename Number>
      static const Number &
      access(const Number &value, const unsigned int)
      {
        return value;
      }
    };

    template <int dim>
    struct EvaluatorTypeTraits<dim, dim>
    {
      using value_type    = Tensor<1, dim>;
      using gradient_type = Tensor<2, dim>;

      static void
      read_value(const double       vector_entry,
                 const unsigned int component,
                 value_type &       result)
      {
        result[component] = vector_entry;
      }

      static void
      write_value(double &           vector_entry,
                  const unsigned int component,
                  const value_type & result)
      {
        vector_entry = result[component];
      }

      static void
      set_gradient(
        const Tensor<1, dim, Tensor<1, dim, VectorizedArray<double>>> &value,
        const unsigned int vector_lane,
        gradient_type &    result)
      {
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int d = 0; d < dim; ++d)
            result[i][d] = value[d][i][vector_lane];
      }

      static void get_gradient(
        Tensor<1, dim, Tensor<1, dim, VectorizedArray<double>>> &value,
        const unsigned int                                       vector_lane,
        const gradient_type &                                    result)
      {
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int d = 0; d < dim; ++d)
            value[d][i][vector_lane] = result[i][d];
      }

      static void
      set_value(const Tensor<1, dim, VectorizedArray<double>> &value,
                const unsigned int                             vector_lane,
                value_type &                                   result)
      {
        for (unsigned int i = 0; i < dim; ++i)
          result[i] = value[i][vector_lane];
      }

      static void get_value(Tensor<1, dim, VectorizedArray<double>> &value,
                            const unsigned int vector_lane,
                            const value_type & result)
      {
        for (unsigned int i = 0; i < dim; ++i)
          value[i][vector_lane] = result[i];
      }

      static double &
      access(value_type &value, const unsigned int component)
      {
        return value[component];
      }

      static const double &
      access(const value_type &value, const unsigned int component)
      {
        return value[component];
      }

      static Tensor<1, dim> &
      access(gradient_type &value, const unsigned int component)
      {
        return value[component];
      }

      static const Tensor<1, dim> &
      access(const gradient_type &value, const unsigned int component)
      {
        return value[component];
      }
    };

    template <>
    struct EvaluatorTypeTraits<1, 1>
    {
      using value_type    = double;
      using gradient_type = Tensor<1, 1>;

      static void
      read_value(const double vector_entry,
                 const unsigned int,
                 value_type &result)
      {
        result = vector_entry;
      }

      static void
      write_value(double &vector_entry,
                  const unsigned int,
                  const value_type &result)
      {
        vector_entry = result;
      }

      static void
      set_gradient(const Tensor<1, 1, VectorizedArray<double>> &value,
                   const unsigned int                           vector_lane,
                   gradient_type &                              result)
      {
        result[0] = value[0][vector_lane];
      }

      static void get_gradient(Tensor<1, 1, VectorizedArray<double>> &value,
                               const unsigned int   vector_lane,
                               const gradient_type &result)
      {
        value[0][vector_lane] = result[0];
      }

      static void
      set_value(const VectorizedArray<double> &value,
                const unsigned int             vector_lane,
                value_type &                   result)
      {
        result = value[vector_lane];
      }

      static void
      get_value(VectorizedArray<double> &value,
                const unsigned int       vector_lane,
                const value_type &       result)
      {
        value[vector_lane] = result;
      }

      template <typename Number>
      static Number &
      access(Number &value, const unsigned int)
      {
        return value;
      }

      template <typename Number>
      static const Number &
      access(const Number &value, const unsigned int)
      {
        return value;
      }
    };
  } // namespace FEPointEvaluation
} // namespace internal



/**
 * This class provides an interface to the evaluation of interpolated solution
 * values and gradients on cells on arbitrary reference point positions. These
 * points can change from cell to cell, both with respect to their quantity as
 * well to the location. The two typical use cases are evaluations on
 * non-matching grids and particle simulations.
 *
 * The functionality is similar to creating an FEValues object with a
 * Quadrature object on the `unit_points` on every cell separately and then
 * calling FEValues::get_function_values or FEValues::get_function_gradients,
 * and for some elements and mappings this is what actually happens
 * internally. For specific combinations of Mapping and FiniteElement
 * realizations, however, there is a much more efficient implementation that
 * avoids the memory allocation and other expensive start-up cost of
 * FEValues. Currently, the functionality is specialized for mappings derived
 * from MappingQGeneric and for finite elements with tensor product structure
 * that work with the @ref matrixfree module. In those cases, the cost implied
 * by this class is similar (or sometimes even somewhat lower) than using
 * `FEValues::reinit(cell)` followed by `FEValues::get_function_gradients`.
 */
template <int n_components, int dim, int spacedim = dim>
class FEPointEvaluation
{
public:
  using value_type = typename internal::FEPointEvaluation::
    EvaluatorTypeTraits<dim, n_components>::value_type;
  using gradient_type = typename internal::FEPointEvaluation::
    EvaluatorTypeTraits<dim, n_components>::gradient_type;

  /**
   * Constructor.
   *
   * @param mapping The Mapping class describing the actual geometry of a cell
   * passed to the evaluate() function.
   *
   * @param fe The FiniteElement object that is used for the evaluation, which
   * is typically the same on all cells to be evaluated.
   *
   * @param first_selected_component For multi-component FiniteElement
   * objects, this parameter allows to select a range of `n_components`
   * components starting from this parameter.
   */
  FEPointEvaluation(const Mapping<dim> &      mapping,
                    const FiniteElement<dim> &fe,
                    const unsigned int        first_selected_component = 0);

  /**
   * This is one of the main functions in this class and the one that does the
   * heavy lifting.
   *
   * @param[in] cell An iterator to the current cell in question
   *
   * @param[in] unit_points List of points in the reference locations of the
   * current cell where the FiniteElement object should be evaluated
   *
   * @param[in] solution_values This array is supposed to contain the unknown
   * values on the element as returned by `cell->get_dof_values(global_vector,
   * solution_values)`.
   *
   * @param[in] evaluation_flags Flags specifying which quantities should be
   * evaluated at the points.
   */
  void
  evaluate(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
           const ArrayView<const Point<dim>> &     unit_points,
           const ArrayView<const double> &         solution_values,
           const EvaluationFlags::EvaluationFlags &evaluation_flags);

  /**
   * This is one of the main functions in this class and the one that does the
   * heavy lifting.
   *
   * @param[in] cell An iterator to the current cell in question
   *
   * @param[in] unit_points List of points in the reference locations of the
   * current cell where the FiniteElement object should be evaluated
   *
   * @param[out] solution_values This array is filled and can be used to
   * during `cell->set_dof_values(global_vector, solution_values)`.
   *
   * @param[in] integration_flags Flags specifying which quantities should be
   * integrated at the points.
   *
   */
  void
  integrate(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
            const ArrayView<const Point<dim>> &     unit_points,
            const ArrayView<double> &               solution_values,
            const EvaluationFlags::EvaluationFlags &integration_flags);

  /**
   * Return the value at quadrature point number @p q_point after a call to
   * FEPointEvaluation::evaluate() with EvaluationFlags::value set, or
   * the value that has been stored there with a call to
   * FEPointEvaluation::submit_value(). If the object is vector-valued, a
   * vector-valued return argument is given. Note that when
   * vectorization is enabled, values from several points are grouped together.
   */
  const value_type &
  get_value(const unsigned int q_point) const;

  /**
   * Write a value to the field containing the values on points
   * with component q_point. Access to the same field as through get_value().
   * If applied before the function FEPointEvaluation::integrate()
   * with EvaluationFlags::values set is called, this specifies the value
   * which is tested by all basis function on the current cell and
   * integrated over.
   */
  void
  submit_value(const value_type &value, const unsigned int q_point);

  /**
   * Return the gradient at quadrature point number @p q_point after a call to
   * FEPointEvaluation::evaluate() with EvaluationFlags::gradient set, or
   * the gradient that has been stored there with a call to
   * FEPointEvaluation::submit_gradient(). If the object is vector-valued, a
   * vector-valued return argument is given. Note that when
   * vectorization is enabled, values from several points are grouped together.
   */
  const gradient_type &
  get_gradient(const unsigned int q_point) const;

  /**
   * Write a contribution that is tested by the gradient to the field
   * containing the values on points with component q_point. Access to the same
   * field as through get_gradient(). If applied before the function
   * FEPointEvaluation::integrate(EvaluationFlags::gradients) is called,
   * this specifies what is tested by all basis function gradients on the
   * current cell and integrated over.
   */
  void
  submit_gradient(const gradient_type &, const unsigned int q_point);

private:
  /**
   * Pointer to the Mapping object passed to the constructor.
   */
  SmartPointer<const Mapping<dim>> mapping;

  /**
   * Pointer to MappingQGeneric class that enables the fast path of this
   * class.
   */
  const MappingQGeneric<dim, spacedim> *mapping_q_generic;

  /**
   * Pointer to the FiniteElement object passed to the constructor.
   */
  SmartPointer<const FiniteElement<dim>> fe;

  /**
   * Description of the 1D polynomial basis for tensor product elements used
   * for the fast path of this class using tensor product evaluators.
   */
  std::vector<Polynomials::Polynomial<double>> poly;

  /**
   * Renumbering between the unknowns of unknowns implied by the FiniteElement
   * class and a lexicographic numbering used for the tensorized code path.
   */
  std::vector<unsigned int> renumber;

  /**
   * Temporary array to store the `solution_values` passed to the evaluate()
   * function in a format compatible with the tensor product evaluators. For
   * vector-valued setups, this array uses a `Tensor<1, n_components>` type to
   * collect the unknowns for a particular basis function.
   */
  std::vector<value_type> solution_renumbered;

  /**
   * Temporary array to store the values at the points.
   */
  std::vector<value_type> values;

  /**
   * Temporary array to store the gradients at the points.
   */
  std::vector<gradient_type> gradients;

  /**
   * Number of unknowns per component, i.e., number of unique basis functions,
   * for the chosen FiniteElement (or base element).
   */
  unsigned int dofs_per_component;

  /**
   * For complicated FiniteElement objects this variable informs us about
   * which unknowns actually carry degrees of freedom in the selected
   * components.
   */
  std::vector<std::array<bool, n_components>> nonzero_shape_function_component;
};

// ----------------------- template and inline function ----------------------


template <int n_components, int dim, int spacedim>
FEPointEvaluation<n_components, dim, spacedim>::FEPointEvaluation(
  const Mapping<dim> &      mapping,
  const FiniteElement<dim> &fe,
  const unsigned int        first_selected_component)
  : mapping(&mapping)
  , mapping_q_generic(
      dynamic_cast<const MappingQGeneric<dim, spacedim> *>(&mapping))
  , fe(&fe)
{
  if (mapping_q_generic != nullptr &&
      internal::MatrixFreeFunctions::ShapeInfo<double>::is_supported(fe))
    {
      internal::MatrixFreeFunctions::ShapeInfo<double> shape_info;
      unsigned int                                     base_element_number = 0;
      unsigned int                                     component           = 0;
      for (; base_element_number < fe.n_base_elements(); ++base_element_number)
        if (component + fe.element_multiplicity(base_element_number) >
            first_selected_component)
          {
            Assert(first_selected_component + n_components <=
                     component + fe.element_multiplicity(base_element_number),
                   ExcMessage("You selected an evaluation which crosses "
                              "different base elements, a case not supported"));
            break;
          }
        else
          component += fe.element_multiplicity(base_element_number);

      shape_info.reinit(QMidpoint<1>(), fe, base_element_number);
      renumber           = shape_info.lexicographic_numbering;
      dofs_per_component = shape_info.dofs_per_component_on_cell;
      poly               = Polynomials::generate_complete_Lagrange_basis(
        QGaussLobatto<1>(shape_info.data[0].fe_degree + 1).get_points());
    }
  if (true /*TODO: as long as the fast path of integrate() is not working*/)
    {
      nonzero_shape_function_component.resize(fe.n_dofs_per_cell());
      for (unsigned int d = 0; d < n_components; ++d)
        {
          const unsigned int component = first_selected_component + d;
          for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
            {
              const bool is_primitive = fe.is_primitive() || fe.is_primitive(i);
              if (is_primitive)
                nonzero_shape_function_component[i][d] =
                  (component == fe.system_to_component_index(i).first);
              else
                nonzero_shape_function_component[i][d] =
                  (fe.get_nonzero_components(i)[component] == true);
            }
        }
    }
}



template <int n_components, int dim, int spacedim>
void
FEPointEvaluation<n_components, dim, spacedim>::evaluate(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const ArrayView<const Point<dim>> &                         unit_points,
  const ArrayView<const double> &                             solution_values,
  const EvaluationFlags::EvaluationFlags &                    evaluation_flag)
{
  AssertDimension(solution_values.size(), fe->dofs_per_cell);
  if (((evaluation_flag & EvaluationFlags::values) ||
       (evaluation_flag & EvaluationFlags::gradients)) &&
      !poly.empty())
    {
      // fast path with tensor product evaluation
      if (solution_renumbered.size() != dofs_per_component)
        solution_renumbered.resize(dofs_per_component);
      for (unsigned int comp = 0; comp < n_components; ++comp)
        for (unsigned int i = 0; i < dofs_per_component; ++i)
          internal::FEPointEvaluation::EvaluatorTypeTraits<dim, n_components>::
            read_value(solution_values[renumber[comp * dofs_per_component + i]],
                       comp,
                       solution_renumbered[i]);

      if (evaluation_flag & EvaluationFlags::values)
        values.resize(unit_points.size());
      if (evaluation_flag & EvaluationFlags::gradients)
        gradients.resize(unit_points.size());

      const std::size_t n_points = unit_points.size();
      const std::size_t n_lanes  = VectorizedArray<double>::size();
      for (unsigned int i = 0; i < n_points; i += n_lanes)
        {
          // convert to vectorized format
          Point<dim, VectorizedArray<double>> vectorized_points;
          for (unsigned int j = 0; j < n_lanes && i + j < n_points; ++j)
            for (unsigned int d = 0; d < dim; ++d)
              vectorized_points[d][j] = unit_points[i + j][d];

          // compute
          const auto val_and_grad =
            internal::evaluate_tensor_product_value_and_gradient(
              poly, solution_renumbered, vectorized_points, poly.size() == 2);

          // convert back to standard format
          if (evaluation_flag & EvaluationFlags::values)
            for (unsigned int j = 0; j < n_lanes && i + j < n_points; ++j)
              internal::FEPointEvaluation::EvaluatorTypeTraits<
                dim,
                n_components>::set_value(val_and_grad.first, j, values[i + j]);
          if (evaluation_flag & EvaluationFlags::gradients)
            for (unsigned int j = 0; j < n_lanes && i + j < n_points; ++j)
              internal::FEPointEvaluation::
                EvaluatorTypeTraits<dim, n_components>::set_gradient(
                  val_and_grad.second, j, gradients[i + j]);
        }

      // let mapping compute the transformation
      if (evaluation_flag & EvaluationFlags::gradients)
        {
          Assert(mapping_q_generic != nullptr, ExcInternalError());
          mapping_q_generic->transform_variable(
            cell,
            mapping_covariant,
            unit_points,
            ArrayView<const gradient_type>(gradients.data(), gradients.size()),
            ArrayView<gradient_type>(gradients.data(), gradients.size()));
        }
    }
  else if ((evaluation_flag & EvaluationFlags::values) ||
           (evaluation_flag & EvaluationFlags::gradients))
    {
      // slow path with FEValues
      const UpdateFlags flags =
        ((evaluation_flag & EvaluationFlags::values) ? update_values :
                                                       update_default) |
        ((evaluation_flag & EvaluationFlags::gradients) ? update_gradients :
                                                          update_default);
      FEValues<dim, spacedim> fe_values(
        *mapping,
        *fe,
        Quadrature<dim>(
          std::vector<Point<dim>>(unit_points.begin(), unit_points.end())),
        flags);
      fe_values.reinit(cell);

      if (evaluation_flag & EvaluationFlags::values)
        {
          values.resize(unit_points.size());
          std::fill(values.begin(), values.end(), value_type());
          for (unsigned int i = 0; i < fe->n_dofs_per_cell(); ++i)
            {
              const double value = solution_values[i];
              for (unsigned int d = 0; d < n_components; ++d)
                if (nonzero_shape_function_component[i][d] &&
                    (fe->is_primitive(i) || fe->is_primitive()))
                  for (unsigned int q = 0; q < unit_points.size(); ++q)
                    internal::FEPointEvaluation::
                      EvaluatorTypeTraits<dim, n_components>::access(values[q],
                                                                     d) +=
                      fe_values.shape_value(i, q) * value;
                else if (nonzero_shape_function_component[i][d])
                  for (unsigned int q = 0; q < unit_points.size(); ++q)
                    internal::FEPointEvaluation::
                      EvaluatorTypeTraits<dim, n_components>::access(values[q],
                                                                     d) +=
                      fe_values.shape_value_component(i, q, d) * value;
            }
        }

      if (evaluation_flag & EvaluationFlags::gradients)
        {
          gradients.resize(unit_points.size());
          std::fill(gradients.begin(), gradients.end(), gradient_type());
          for (unsigned int i = 0; i < fe->n_dofs_per_cell(); ++i)
            {
              const double value = solution_values[i];
              for (unsigned int d = 0; d < n_components; ++d)
                if (nonzero_shape_function_component[i][d] &&
                    (fe->is_primitive(i) || fe->is_primitive()))
                  for (unsigned int q = 0; q < unit_points.size(); ++q)
                    internal::FEPointEvaluation::
                      EvaluatorTypeTraits<dim, n_components>::access(
                        gradients[q], d) += fe_values.shape_grad(i, q) * value;
                else if (nonzero_shape_function_component[i][d])
                  for (unsigned int q = 0; q < unit_points.size(); ++q)
                    internal::FEPointEvaluation::EvaluatorTypeTraits<
                      dim,
                      n_components>::access(gradients[q], d) +=
                      fe_values.shape_grad_component(i, q, d) * value;
            }
        }
    }
}



template <int n_components, int dim, int spacedim>
void
FEPointEvaluation<n_components, dim, spacedim>::integrate(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const ArrayView<const Point<dim>> &                         unit_points,
  const ArrayView<double> &                                   solution_values,
  const EvaluationFlags::EvaluationFlags &                    integration_flags)
{
  AssertDimension(solution_values.size(), fe->dofs_per_cell);
  if (false /*TODO*/ && (((integration_flags & EvaluationFlags::values) ||
                          (integration_flags & EvaluationFlags::gradients)) &&
                         !poly.empty()))
    {
      Assert(false, ExcNotImplemented());
      // fast path with tensor product integration

      if (solution_renumbered.size() != dofs_per_component)
        solution_renumbered.resize(dofs_per_component);

      // let mapping compute the transformation
      if (integration_flags & EvaluationFlags::gradients)
        {
          Assert(mapping_q_generic != nullptr, ExcInternalError());

          mapping_q_generic->transform_variable(
            cell,
            mapping_covariant,
            unit_points,
            ArrayView<const gradient_type>(gradients.data(), gradients.size()),
            ArrayView<gradient_type>(gradients.data(), gradients.size()));
        }

      if (integration_flags & EvaluationFlags::values)
        AssertIndexRange(unit_points.size(), values.size() + 1);
      if (integration_flags & EvaluationFlags::gradients)
        AssertIndexRange(unit_points.size(), gradients.size() + 1);

      const std::size_t n_points = unit_points.size();
      const std::size_t n_lanes  = VectorizedArray<double>::size();
      for (unsigned int i = 0; i < n_points; i += n_lanes)
        {
          // convert to vectorized format
          Point<dim, VectorizedArray<double>> vectorized_points;
          for (unsigned int j = 0; j < n_lanes && i + j < n_points; ++j)
            for (unsigned int d = 0; d < dim; ++d)
              vectorized_points[d][j] = unit_points[i + j][d];

          typename internal::ProductTypeNoPoint<value_type,
                                                VectorizedArray<double>>::type
            value;
          Tensor<1,
                 dim,
                 typename internal::ProductTypeNoPoint<
                   value_type,
                   VectorizedArray<double>>::type>
            gradient;

          // convert back to standard format
          if (integration_flags & EvaluationFlags::values)
            for (unsigned int j = 0; j < n_lanes && i + j < n_points; ++j)
              internal::FEPointEvaluation::EvaluatorTypeTraits<
                dim,
                n_components>::get_value(value, j, values[i + j]);
          if (integration_flags & EvaluationFlags::gradients)
            for (unsigned int j = 0; j < n_lanes && i + j < n_points; ++j)
              internal::FEPointEvaluation::EvaluatorTypeTraits<
                dim,
                n_components>::get_gradient(gradient, j, gradients[i + j]);

          // compute
          internal::integrate_tensor_product_value_and_gradient(
            poly,
            solution_renumbered,
            value,
            gradient,
            vectorized_points,
            poly.size() == 2);
        }

      for (unsigned int comp = 0; comp < n_components; ++comp)
        for (unsigned int i = 0; i < dofs_per_component; ++i)
          internal::FEPointEvaluation::EvaluatorTypeTraits<dim, n_components>::
            write_value(
              solution_values[renumber[comp * dofs_per_component + i]],
              comp,
              solution_renumbered[i]);
    }
  else if ((integration_flags & EvaluationFlags::values) ||
           (integration_flags & EvaluationFlags::gradients))
    {
      // slow path with FEValues
      const UpdateFlags flags =
        ((integration_flags & EvaluationFlags::values) ? update_values :
                                                         update_default) |
        ((integration_flags & EvaluationFlags::gradients) ? update_gradients :
                                                            update_default);
      FEValues<dim, spacedim> fe_values(
        *mapping,
        *fe,
        Quadrature<dim>(
          std::vector<Point<dim>>(unit_points.begin(), unit_points.end())),
        flags);
      fe_values.reinit(cell);

      std::fill(solution_values.begin(), solution_values.end(), 0.0);

      if (integration_flags & EvaluationFlags::values)
        {
          AssertIndexRange(unit_points.size(), values.size() + 1);
          for (unsigned int i = 0; i < fe->n_dofs_per_cell(); ++i)
            {
              for (unsigned int d = 0; d < n_components; ++d)
                if (nonzero_shape_function_component[i][d] &&
                    (fe->is_primitive(i) || fe->is_primitive()))
                  for (unsigned int q = 0; q < unit_points.size(); ++q)
                    solution_values[i] +=
                      fe_values.shape_value(i, q) *
                      internal::FEPointEvaluation::EvaluatorTypeTraits<
                        dim,
                        n_components>::access(values[q], d);
                else if (nonzero_shape_function_component[i][d])
                  for (unsigned int q = 0; q < unit_points.size(); ++q)
                    solution_values[i] +=
                      fe_values.shape_value_component(i, q, d) *
                      internal::FEPointEvaluation::EvaluatorTypeTraits<
                        dim,
                        n_components>::access(values[q], d);
            }
        }

      if (integration_flags & EvaluationFlags::gradients)
        {
          AssertIndexRange(unit_points.size(), gradients.size() + 1);
          for (unsigned int i = 0; i < fe->n_dofs_per_cell(); ++i)
            {
              for (unsigned int d = 0; d < n_components; ++d)
                if (nonzero_shape_function_component[i][d] &&
                    (fe->is_primitive(i) || fe->is_primitive()))
                  for (unsigned int q = 0; q < unit_points.size(); ++q)
                    solution_values[i] +=
                      fe_values.shape_grad(i, q) *
                      internal::FEPointEvaluation::EvaluatorTypeTraits<
                        dim,
                        n_components>::access(gradients[q], d);
                else if (nonzero_shape_function_component[i][d])
                  for (unsigned int q = 0; q < unit_points.size(); ++q)
                    solution_values[i] +=
                      fe_values.shape_grad_component(i, q, d) *
                      internal::FEPointEvaluation::EvaluatorTypeTraits<
                        dim,
                        n_components>::access(gradients[q], d);
            }
        }
    }
}



template <int n_components, int dim, int spacedim>
inline const typename FEPointEvaluation<n_components, dim, spacedim>::value_type
  &
  FEPointEvaluation<n_components, dim, spacedim>::get_value(
    const unsigned int q_point) const
{
  AssertIndexRange(q_point, values.size());
  return values[q_point];
}



template <int n_components, int dim, int spacedim>
inline const typename FEPointEvaluation<n_components, dim, spacedim>::
  gradient_type &
  FEPointEvaluation<n_components, dim, spacedim>::get_gradient(
    const unsigned int q_point) const
{
  AssertIndexRange(q_point, gradients.size());
  return gradients[q_point];
}



template <int n_components, int dim, int spacedim>
inline void
FEPointEvaluation<n_components, dim, spacedim>::submit_value(
  const value_type & value,
  const unsigned int q_point)
{
  if (values.size() <= q_point)
    values.resize(values.size() + 1);

  values[q_point] = value;
}



template <int n_components, int dim, int spacedim>
inline void
FEPointEvaluation<n_components, dim, spacedim>::submit_gradient(
  const gradient_type &gradient,
  const unsigned int   q_point)
{
  if (gradients.size() <= q_point)
    gradients.resize(gradients.size() + 1);

  gradients[q_point] = gradient;
}

DEAL_II_NAMESPACE_CLOSE

#endif
