// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
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

#ifndef dealii_fe_coupling_values_h
#define dealii_fe_coupling_values_h

#include <deal.II/base/config.h>

#include <deal.II/algorithms/general_data_storage.h>

#include <deal.II/base/subscriptor.h>
#include <deal.II/base/thread_local_storage.h>
#include <deal.II/base/utilities.h>

#include <deal.II/fe/fe_values_views.h>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
template <int dim, int spacedim>
class FEValuesBase;
#endif

namespace FEValuesViews
{
  /**
   * A class that stores the renumbering of degrees of freedom and quadrature
   * points for the RenumberedView class. This class is common to all possible
   * renumbered views, and therefore it can be shared between all views that use
   * the same set of renumbering vectors.
   *
   * The renumbering is stored in two vectors, one for the degrees of freedom
   * and one for the quadrature points, and it is used by the
   * FEValuesViews::RenumberedView class.
   */
  struct RenumberingData
  {
    /**
     * Construct a new renumbering data object.
     *
     * The @p dof_renumbering vector is used to renumber the degrees of freedom,
     * while the @p quadrature_renumbering vector is used to renumber the
     * quadrature points. An empty renumbering vector simply means that no
     * renumbering is performed.
     *
     * @note The renumbering vectors are *not* required to match the dimensions
     * of the underlying view, i.e., you could use this class to only run over
     * half of the degrees of freedom, and integrating only over half of the
     * current cell by selecting a subset of quadrature points, if you wish to
     * do so, or to run over some of the degrees of freedom more than once with
     * regard to the underlying view. The important part is that every entry of
     * each renumbering vector is a legal index within the underlying view
     * object. We allow the dof renumbering vector to contain
     * `numbers::invalid_unsigned_int` values, which are ignored, and produce
     * zero values, gradients, hessians, etc. for the corresponding shape
     * function.
     *
     * An example of the renumbering ordering is the following. Assume that you
     * have an FE_Q(1) finite element in dimension two (i.e., four degrees of
     * freedom), and that you are manually taking care of an additional degree
     * of freedom (say, a custom shape function in the middle of the cell),
     * which, for your own convenience, should be numbered locally in the middle
     * of the others. You would like your loops on degrees of freedom to run
     * over all degrees of freedom, and to have the same structure as if the
     * additional degree of freedom was not there. Then the renumbering data
     * object should look like this:
     * @code
     * ...
     * RenumberingData data({{0, 1, numbers::invalid_unsigned_int, 2, 3}});
     * @endcode
     *
     * Using this RenumberingData object, the RenumberedView will return zero
     * when we ask for values, gradients, etc., of the shape function with index
     * two.
     */
    RenumberingData(
      const unsigned int n_inner_dofs = numbers::invalid_unsigned_int,
      const unsigned int n_inner_quadrature_points =
        numbers::invalid_unsigned_int,
      const std::vector<unsigned int> &dof_renumbering        = {},
      const std::vector<unsigned int> &quadrature_renumbering = {});

    /**
     * Make sure we never copy this object, for performance reasons.
     */
    RenumberingData(const RenumberingData &other) = delete;

    /**
     * The number of dofs in the underlying view (before any renumbering).
     */
    const unsigned int n_inner_dofs;

    /**
     * The number of dofs in the renumbered view.
     */
    const unsigned int n_dofs;

    /**
     * The number of quadrature points in the underlying view (before any
     * renumbering).
     */
    const unsigned int n_inner_quadrature_points;

    /**
     * The number of quadrature points in the renumbered view.
     */
    const unsigned int n_quadrature_points;

    /**
     * The renumbering of degrees of freedom.
     */
    const std::vector<unsigned int> dof_renumbering;

    /**
     * The renumbering of quadrature points.
     */
    const std::vector<unsigned int> quadrature_renumbering;

    /**
     * General data storage to store temporary vectors.
     *
     * When the renumbering vectors are non empty, the RenumberedView class may
     * need to construct temporary vectors to store the values of the solution
     * and/or of its gradients with the sizes given by the underlying view
     * object. Unfortunately, we don't know beforehand the Number types with
     * which the vectors will be instantiated, so we cannot use a simple cache
     * internally, and we use a GeneralDataStorage object to avoid allocating
     * memory at each call.
     */
    mutable Threads::ThreadLocalStorage<GeneralDataStorage> data_storage;
  };

  /**
   * A class that provides a renumbered view to a given FEValuesViews object.
   *
   * In general, the order of the degrees of freedom and quadrature points
   * follows the one of the FEValues object, which itself uses the numbering
   * provided by the FiniteElement and Quadrature objects it uses. However, in
   * some cases, it is convenient to group together degrees of freedom and
   * quadrature points in a different order, to select only a subset of the
   * degrees of freedom, or to combine two different sets of degrees of freedom
   * together. This class provides a view to a given FEValuesViews object, where
   * the degrees of freedom and quadrature points are renumbered according to
   * the given RenumberingData object (see there for a documentation of how the
   * renumbering is interpreted).
   *
   * Users will typically not use this class directly, but rather pass an
   * extractor from the FEValuesExtractors namespace to the FECouplingValues
   * class, which returns an object of this type. This is the same mechanism
   * used in the FEValues classes, where passing an extractor returns an
   * FEValuesViews object, and the user rarely instantiates an object of this
   * type.
   */
  template <typename ViewType>
  class RenumberedView
  {
  public:
    /**
     * An alias for the data type of values of the underlying view.
     */
    using value_type = typename ViewType::value_type;

    /**
     * An alias for the data type of gradients of the underlying view.
     */
    using gradient_type = typename ViewType::gradient_type;

    /**
     * An alias for the data type of the product of a @p Number and the values
     * of the underlying view type.
     */
    template <typename Number>
    using solution_value_type =
      typename ViewType::template solution_value_type<Number>;

    /**
     * An alias for the data type of the product of a @p Number and the gradients
     * of the underlying view type.
     */
    template <typename Number>
    using solution_gradient_type =
      typename ViewType::template solution_gradient_type<Number>;

    /**
     * Construct a new RenumberedView object.
     *
     * The renumbering information is taken from the class RenumberingData (see
     * there for a documentation of how the renumbering is interpreted).
     *
     * @note The renumbering information is stored as a reference in this
     * object, so you have to make sure that the RenumberingData object lives
     * longer than this object.
     *
     * @param view The underlying FEValuesViews object.
     * @param data A RenumberingData object, containing renumbering information.
     */
    RenumberedView(const ViewType &view, const RenumberingData &data);

    /**
     * Return the value of the underlying view, for the shape function and
     * quadrature point selected by the arguments.
     *
     * @param shape_function Number of the shape function to be evaluated. Note
     * that this number runs from zero to size of the renumbering vector
     * provided at construction time, or to `dofs_per_cell`, if the renumbering
     * vector is empty.
     *
     * @param q_point Number of the quadrature point at which the function is to
     * be evaluated. Note that this number runs from zero to the size of the
     * renumbering vector provided at construction time, or to
     * `n_quadrature_points`, if the renumbering vector is empty.
     *
     * @dealiiRequiresUpdateFlags{update_values}
     */
    value_type
    value(const unsigned int shape_function, const unsigned int q_point) const;

    /**
     * Return the gradient of the underlying view, for the shape function and
     * quadrature point selected by the arguments.
     *
     * @param shape_function Number of the shape function to be evaluated. This
     * number runs from zero to size of the renumbering vector provided at
     * construction time, or to `dofs_per_cell`, if the renumbering vector is
     * empty.
     *
     * @param q_point Number of the quadrature point at which the function is to
     * be evaluated. This number runs from zero to the size of the renumbering
     * vector provided at construction time, or to `n_quadrature_points`, if the
     * renumbering vector is empty.
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    gradient_type
    gradient(const unsigned int shape_function,
             const unsigned int q_point) const;

    /**
     * Return the values of the underlying view characterized by
     * <tt>fe_function</tt> at the renumbered quadrature points.
     *
     * This function is the equivalent of the FEValuesBase::get_function_values
     * function but it only works on the selected view.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the values of shape functions (i.e., @p value_type) times the
     * type used to store the values of the unknowns $U_j$ of your finite
     * element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_values}
     */
    template <typename Number>
    void
    get_function_values(const ReadVector<Number>                 &fe_function,
                        std::vector<solution_value_type<Number>> &values) const;

    /**
     * Same as above, but using a vector of renumbered local degree-of-freedom
     * values. In other words, instead of extracting the nodal values of the
     * degrees of freedom located on the current cell from a global vector
     * associated with a DoFHandler object (as the function above does), this
     * function instead takes these local nodal values through its first
     * argument.
     *
     * @param[in] dof_values A vector of local nodal values. This vector must
     *   have a length equal to the size of the dof renumbering vector.
     *
     * @param[out] values A vector of values of the given finite element field,
     *   at the renumbered quadrature points on the current object.
     *
     * @tparam InputVector The @p InputVector type must allow creation of an
     *   ArrayView object from it; this is satisfied by the `std::vector` class,
     *   among others.
     *
     * @dealiiRequiresUpdateFlags{update_values}
     */
    template <class InputVector>
    void
    get_function_values_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_value_type<typename InputVector::value_type>>
        &values) const;

    /**
     * Return the gradients of the underlying view characterized by
     * <tt>fe_function</tt> at the renumbered quadrature points.
     *
     * This function is the equivalent of the
     * FEValuesBase::get_function_gradients function but it only works on the
     * selected view.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the gradients of shape functions (i.e., @p value_type) times the
     * type used to store the gradients of the unknowns $U_j$ of your finite
     * element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    template <typename Number>
    void
    get_function_gradients(
      const ReadVector<Number>                    &fe_function,
      std::vector<solution_gradient_type<Number>> &gradients) const;

    /**
     * Same as above, but using a vector of renumbered local degree-of-freedom
     * gradients. In other words, instead of extracting the nodal values of the
     * degrees of freedom located on the current cell from a global vector
     * associated with a DoFHandler object (as the function above does), this
     * function instead takes these local nodal values through its first
     * argument.
     *
     * @param[in] dof_values A vector of local nodal values. This vector must
     *   have a length equal to the size of the dof renumbering vector.
     *
     * @param[out] gradients A vector of gradients of the given finite element
     *   field, at the renumbered quadrature points on the current object.
     *
     * @tparam InputVector The @p InputVector type must allow creation of an
     *   ArrayView object from it; this is satisfied by the `std::vector` class,
     *   among others.
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    template <class InputVector>
    void
    get_function_gradients_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_gradient_type<typename InputVector::value_type>>
        &gradients) const;

  private:
    /**
     * The data structure that stores the renumbering of the degrees of freedom
     * and of the quadrature points.
     */
    const RenumberingData &data;


    /**
     * Helper function that constructs a unique name for a container, based on a
     * string prefix, on its size, and on the type stored in the container.
     *
     * When the renumbering vectors are non empty, this class may need to
     * construct temporary vectors to store the values of the solution and/or of
     * its gradients with the sizes given by the underlying view object.
     * Unfortunately, we don't know before hand the Number types with which the
     * vectors will be instantiated, so we cannot use a simple cache internally,
     * and we use a GeneralDataStorage object to avoid allocating memory at each
     * call.
     *
     * This function constructs a unique name for temporary containers that will
     * be stored upon the first request in the internal GeneralDataStorage
     * object, to access the
     */
    template <typename Number>
    std::string
    get_unique_container_name(const std::string &prefix,
                              const unsigned int size,
                              const Number      &exemplar_number) const;

    /**
     * Produce an inner vector compatible with the inner view, after copying the
     * values with the correct numbering from the outer vector.
     */
    template <typename InputVector>
    const InputVector &
    outer_to_inner_dofs(const InputVector &outer_vector) const;

    /**
     * Produce an inner vector compatible with the inner view, and zero out its
     * entries if necessary.
     */
    template <typename ValueType>
    std::vector<ValueType> &
    outer_to_inner_values(std::vector<ValueType> &outer_values) const;

    /**
     * Return the outer argument renumbered according to the quadrature
     * renumbering. The values in the inner values are copied to the right
     * position in the outer vector.
     */
    template <typename ValueType>
    void
    inner_to_outer_values(const std::vector<ValueType> &inner_values,
                          std::vector<ValueType>       &outer_values) const;

    /**
     * Store a reference to the underlying view.
     */
    const ViewType &view;
  };
} // namespace FEValuesViews


#ifndef DOXYGEN


/*------------------------ Inline functions: namespace FEValuesViews --------*/

namespace FEValuesViews
{
  RenumberingData::RenumberingData(
    const unsigned int               n_inner_dofs,
    const unsigned int               n_inner_quadrature_points,
    const std::vector<unsigned int> &dof_renumbering,
    const std::vector<unsigned int> &quadrature_renumbering)
    : n_inner_dofs(n_inner_dofs)
    , n_dofs(dof_renumbering.empty() ? n_inner_dofs : dof_renumbering.size())
    , n_inner_quadrature_points(n_inner_quadrature_points)
    , n_quadrature_points(quadrature_renumbering.empty() ?
                            n_inner_quadrature_points :
                            quadrature_renumbering.size())
    , dof_renumbering(dof_renumbering)
    , quadrature_renumbering(quadrature_renumbering)
  {
// Check that the renumbering vectors are valid.
#  ifdef DEBUG
    // While for dofs we admit invalid values, this is not the case for
    // quadrature points.
    for (const auto i : dof_renumbering)
      Assert(i < n_inner_dofs || i == numbers::invalid_unsigned_int,
             ExcIndexRange(i, 0, n_inner_dofs));

    for (const auto q : quadrature_renumbering)
      AssertIndexRange(q, n_inner_quadrature_points);
#  endif
  }



  template <typename ViewType>
  RenumberedView<ViewType>::RenumberedView(const ViewType        &view,
                                           const RenumberingData &data)
    : view(view)
    , data(data)
  {}



  template <typename ViewType>
  template <typename Number>
  inline std::string
  RenumberedView<ViewType>::get_unique_container_name(
    const std::string &prefix,
    const unsigned int size,
    const Number      &exemplar_number) const
  {
    return prefix + "_" + Utilities::int_to_string(size) + "_" +
           Utilities::type_to_string(exemplar_number);
  }



  template <typename ViewType>
  typename RenumberedView<ViewType>::value_type
  RenumberedView<ViewType>::value(const unsigned int shape_function,
                                  const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, data.n_dofs);
    AssertIndexRange(q_point, data.n_quadrature_points);

    const auto inner_shape_function = data.dof_renumbering.empty() ?
                                        shape_function :
                                        data.dof_renumbering[shape_function];
    const auto inner_q_point        = data.quadrature_renumbering.empty() ?
                                        q_point :
                                        data.quadrature_renumbering[q_point];
    if (inner_shape_function == numbers::invalid_unsigned_int)
      return value_type(0);
    else
      {
        AssertIndexRange(inner_shape_function, data.n_inner_dofs);
        AssertIndexRange(inner_q_point, data.n_inner_quadrature_points);
        return view.value(inner_shape_function, inner_q_point);
      }
  }



  template <typename ViewType>
  typename RenumberedView<ViewType>::gradient_type
  RenumberedView<ViewType>::gradient(const unsigned int shape_function,
                                     const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, data.n_dofs);
    AssertIndexRange(q_point, data.n_quadrature_points);

    const auto inner_shape_function = data.dof_renumbering.empty() ?
                                        shape_function :
                                        data.dof_renumbering[shape_function];
    const auto inner_q_point        = data.quadrature_renumbering.empty() ?
                                        q_point :
                                        data.quadrature_renumbering[q_point];
    if (inner_shape_function == numbers::invalid_unsigned_int)
      return gradient_type();
    else
      return view.gradient(inner_shape_function, inner_q_point);
  }



  template <typename ViewType>
  template <typename ValueType>
  std::vector<ValueType> &
  RenumberedView<ViewType>::outer_to_inner_values(
    std::vector<ValueType> &outer_values) const
  {
    AssertDimension(outer_values.size(), data.n_quadrature_points);
    if (data.quadrature_renumbering.empty())
      {
        return outer_values;
      }
    else
      {
        const auto name =
          get_unique_container_name("RenumberedView::outer_to_inner_values",
                                    data.n_inner_quadrature_points,
                                    outer_values[0]);
        auto &inner_values =
          data.data_storage.get()
            .template get_or_add_object_with_name<std::vector<ValueType>>(
              name, data.n_inner_quadrature_points);
        return inner_values;
      }
  }



  template <typename ViewType>
  template <typename VectorType>
  const VectorType &
  RenumberedView<ViewType>::outer_to_inner_dofs(
    const VectorType &outer_dofs) const
  {
    AssertDimension(outer_dofs.size(), data.n_dofs);
    if (data.dof_renumbering.empty())
      {
        return outer_dofs;
      }
    else
      {
        const auto name =
          get_unique_container_name("RenumberedView::outer_to_inner_dofs",
                                    data.n_inner_dofs,
                                    outer_dofs[0]);

        auto &inner_dofs = data.data_storage.get()
                             .template get_or_add_object_with_name<VectorType>(
                               name, data.n_inner_dofs);
        for (unsigned int i = 0; i < data.n_dofs; ++i)
          {
            const auto inner_i = data.dof_renumbering[i];
            if (inner_i != numbers::invalid_unsigned_int)
              {
                AssertIndexRange(inner_i, data.n_inner_dofs);
                inner_dofs[inner_i] = outer_dofs[i];
              }
          }
        return inner_dofs;
      }
  }



  template <typename ViewType>
  template <typename ValueType>
  void
  RenumberedView<ViewType>::inner_to_outer_values(
    const std::vector<ValueType> &inner_values,
    std::vector<ValueType>       &outer_values) const
  {
    AssertDimension(outer_values.size(), data.n_quadrature_points);
    AssertDimension(inner_values.size(), data.n_inner_quadrature_points);
    if (data.quadrature_renumbering.empty())
      {
        Assert(&inner_values == &outer_values, ExcInternalError());
        return;
      }
    for (unsigned int i = 0; i < data.quadrature_renumbering.size(); ++i)
      {
        outer_values[i] = inner_values[data.quadrature_renumbering[i]];
      }
  }



  template <typename ViewType>
  template <typename Number>
  void
  RenumberedView<ViewType>::get_function_values(
    const ReadVector<Number>                 &fe_function,
    std::vector<solution_value_type<Number>> &values) const
  {
    auto &inner_values = outer_to_inner_values(values);
    view.get_function_values(fe_function, inner_values);
    inner_to_outer_values(inner_values, values);
  }



  template <typename ViewType>
  template <typename InputVector>
  void
  RenumberedView<ViewType>::get_function_values_from_local_dof_values(
    const InputVector &dof_values,
    std::vector<solution_value_type<typename InputVector::value_type>> &values)
    const
  {
    const auto &inner_dof_values = outer_to_inner_dofs(dof_values);
    auto       &inner_values     = outer_to_inner_values(values);

    view.get_function_values_from_local_dof_values(inner_dof_values,
                                                   inner_values);
    inner_to_outer_values(inner_values, values);
  }


  template <typename ViewType>
  template <typename Number>
  void
  RenumberedView<ViewType>::get_function_gradients(
    const ReadVector<Number>                    &fe_function,
    std::vector<solution_gradient_type<Number>> &gradients) const
  {
    auto &inner_gradients = outer_to_inner_values(gradients);
    view.get_function_gradients(fe_function, inner_gradients);
    inner_to_outer_values(inner_gradients, gradients);
  }



  template <typename ViewType>
  template <typename InputVector>
  void
  RenumberedView<ViewType>::get_function_gradients_from_local_dof_values(
    const InputVector &dof_values,
    std::vector<solution_gradient_type<typename InputVector::value_type>>
      &gradients) const
  {
    const auto &inner_dof_values = outer_to_inner_dofs(dof_values);
    auto       &inner_gradients  = outer_to_inner_values(gradients);

    view.get_function_gradients_from_local_dof_values(inner_dof_values,
                                                      inner_gradients);
    inner_to_outer_values(inner_gradients, gradients);
  }
} // namespace FEValuesViews

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
