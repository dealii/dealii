// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_interface_values_h
#define dealii_fe_interface_values_h

#include <deal.II/base/config.h>

#include <deal.II/base/std_cxx20/iota_view.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/q_collection.h>

DEAL_II_NAMESPACE_OPEN

#ifndef DOXYGEN
template <int dim, int spacedim>
class FEInterfaceValues;
#endif

/**
 * Namespace for views you get from accessing FEInterfaceValues using an
 * extractor.
 */
namespace FEInterfaceViews
{
  /**
   * The base class for the views.
   */
  template <int dim, int spacedim = dim>
  class Base
  {
  public:
    /**
     * The constructor.
     */
    Base(const FEInterfaceValues<dim, spacedim> &fe_interface);

  protected:
    /**
     * Store a pointer to the FEInterfaceValues instance.
     */
    const FEInterfaceValues<dim, spacedim> *fe_interface;

    /**
     * Get the local degree-of-freedom values associated with the internally
     * initialized cell interface.
     */
    template <class InputVector, class OutputVector>
    void
    get_local_dof_values(const InputVector &dof_values,
                         OutputVector      &local_dof_values) const;
  };



  /**
   * The view of a scalar variable for FEInterfaceValues.
   */
  template <int dim, int spacedim = dim>
  class Scalar : public Base<dim, spacedim>
  {
  public:
    /**
     * This is the type returned for values.
     */
    using value_type = double;

    /**
     * This is the type returned for gradients, for example from
     * average_of_gradients().
     */
    using gradient_type =
      typename FEValuesViews::Scalar<dim, spacedim>::gradient_type;

    /**
     * This is the type returned for hessians, for example from
     * jump_in_hessians().
     */
    using hessian_type =
      typename FEValuesViews::Scalar<dim, spacedim>::hessian_type;

    /**
     * This is the type returned for third derivatives, for example from
     * jump_in_hessians().
     */
    using third_derivative_type =
      typename FEValuesViews::Scalar<dim, spacedim>::third_derivative_type;

    /**
     * An alias for the data type of the product of a @p Number and the
     * values of the view this class provides. This is the data type of
     * scalar components of a finite element field whose degrees of
     * freedom are described by a vector with elements of type @p Number.
     */
    template <typename Number>
    using solution_value_type = typename ProductType<Number, value_type>::type;

    /**
     * An alias for the data type of the product of a @p Number and the
     * gradients of the view this class provides. This is the data type of
     * scalar components of a finite element field whose degrees of
     * freedom are described by a vector with elements of type @p Number.
     */
    template <typename Number>
    using solution_gradient_type =
      typename ProductType<Number, gradient_type>::type;

    /**
     * An alias for the data type of the product of a @p Number and the
     * Hessians of the view this class provides. This is the data type of
     * scalar components of a finite element field whose degrees of
     * freedom are described by a vector with elements of type @p Number.
     */
    template <typename Number>
    using solution_hessian_type =
      typename ProductType<Number, hessian_type>::type;

    /**
     * An alias for the data type of the product of a @p Number and the
     * third derivatives of the view this class provides. This is the data type
     * of scalar components of a finite element field whose degrees of
     * freedom are described by a vector with elements of type @p Number.
     */
    template <typename Number>
    using solution_third_derivative_type =
      typename ProductType<Number, third_derivative_type>::type;

    /**
     * Constructor for an object that represents a single scalar component
     */
    Scalar(const FEInterfaceValues<dim, spacedim> &fe_interface,
           const unsigned int                      component);

    /**
     * @name Access to shape functions
     */
    /** @{ */

    /**
     * Return the value of the shape function
     * with interface dof index @p interface_dof_index in
     * quadrature point @p q_point of the component selected by this view.
     *
     * The argument @p here_or_there selects between the upstream value and
     * the downstream value as defined by the direction of the normal vector
     * in this quadrature point. If @p here_or_there is true, the shape
     * functions from the first cell of the interface is used.
     *
     * In other words, this function returns the limit of the value of the shape
     * function in the given quadrature point when approaching it from one of
     * the two cells of the interface.
     *
     * @note This function is typically used to pick the upstream or downstream
     * value based on a direction. This can be achieved by using
     * <code>(direction * normal)>0</code> as the first argument of this
     * function.
     */
    value_type
    value(const bool         here_or_there,
          const unsigned int interface_dof_index,
          const unsigned int q_point) const;

    /** @} */

    /**
     * @name Access to jumps in shape functions and their derivatives
     */
    /** @{ */

    /**
     * Return the jump $\jump{u}=u_1 - u_2$ on the interface for the shape
     * function
     * @p interface_dof_index in the quadrature point @p q_point
     * of the component selected by this view.
     *
     * @note The name of the function is supposed to be read as "the jump
     *   (singular) in the values (plural: one or two possible values) of
     *   the shape function (singular)".
     */
    value_type
    jump_in_values(const unsigned int interface_dof_index,
                   const unsigned int q_point) const;

    /**
     * Return the jump of the gradient $\jump{nabla u}$ on the interface for
     * the shape
     * function @p interface_dof_index in the quadrature point @p q_point
     * of the component selected by this view.
     *
     * @note The name of the function is supposed to be read as "the jump
     *   (singular) in the gradients (plural: one or two possible gradients)
     *   of the shape function (singular)".
     */
    gradient_type
    jump_in_gradients(const unsigned int interface_dof_index,
                      const unsigned int q_point) const;

    /**
     * Return the jump in the gradient $\jump{\nabla u}=\nabla u_{\text{cell0}}
     * - \nabla u_{\text{cell1}}$ on the interface for the shape function @p
     * interface_dof_index at the quadrature point @p q_point of
     * the component selected by this view.
     *
     * @note The name of the function is supposed to be read as "the jump
     *   (singular) in the Hessians (plural: one or two possible values
     *   for the second derivative) of the shape function (singular)".
     */
    hessian_type
    jump_in_hessians(const unsigned int interface_dof_index,
                     const unsigned int q_point) const;

    /**
     * Return the jump in the third derivative $\jump{\nabla^3 u} = \nabla^3
     * u_{\text{cell0}} - \nabla^3 u_{\text{cell1}}$ on the interface for the
     * shape function @p interface_dof_index at the quadrature point @p q_point of
     * the component selected by this view.
     *
     * @note The name of the function is supposed to be read as "the jump
     *   (singular) in the third derivatives (plural: one or two possible values
     *   for the third derivative) of the shape function (singular)".
     */
    third_derivative_type
    jump_in_third_derivatives(const unsigned int interface_dof_index,
                              const unsigned int q_point) const;

    /** @} */

    /**
     * @name Access to the average of shape functions and their derivatives
     */
    /** @{ */

    /**
     * Return the average value $\average{u}=\frac{1}{2}(u_1 + u_2)$ on the
     * interface for the shape
     * function @p interface_dof_index in the quadrature point @p q_point
     * of the component selected by this view.
     *
     * @note The name of the function is supposed to be read as "the average
     *   (singular) of the values (plural: one or two possible values) of
     *   the shape function (singular)".
     */
    value_type
    average_of_values(const unsigned int interface_dof_index,
                      const unsigned int q_point) const;

    /**
     * Return the average of the gradient $\average{\nabla u}$ on the interface
     * for the shape
     * function @p interface_dof_index in the quadrature point @p q_point
     * of the component selected by this view.
     *
     * @note The name of the function is supposed to be read as "the average
     *   (singular) of the gradients (plural: one or two possible values of
     *   the derivative) of the shape function (singular)".
     */
    gradient_type
    average_of_gradients(const unsigned int interface_dof_index,
                         const unsigned int q_point) const;

    /**
     * Return the average of the Hessian $\average{\nabla^2 u} =
     * \frac{1}{2}\nabla^2 u_{\text{cell0}} + \frac{1}{2} \nabla^2
     * u_{\text{cell1}}$ on the interface
     * for the shape function @p interface_dof_index at the quadrature point @p
     * q_point of the component selected by this view.
     *
     * @note The name of the function is supposed to be read as "the average
     *   (singular) in the Hessians (plural: one or two possible values of
     *   the second derivative) of the shape function (singular)".
     */
    hessian_type
    average_of_hessians(const unsigned int interface_dof_index,
                        const unsigned int q_point) const;

    /** @} */

    /**
     * @name Access to values of global finite element fields
     */
    /** @{ */

    /**
     * Return the values of the selected scalar component of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell interface selected the last time
     * the <tt>reinit</tt> function of the FEInterfaceValues object was called.
     *
     * The argument @p here_or_there selects between the value on cell 0 (here, @p true)
     * and cell 1 (there, @p false). You can also interpret it as "upstream" (@p true)
     * and "downstream" (@p false) as defined by the direction of the normal
     * vector in this quadrature point. If @p here_or_there is true, the values
     * from the first cell of the interface is used.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the values of shape functions (i.e., @p value_type) times the
     * type used to store the values of the unknowns $U_j$ of your finite
     * element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_values}
     */
    template <class InputVector>
    void
    get_function_values(
      const bool         here_or_there,
      const InputVector &fe_function,
      std::vector<solution_value_type<typename InputVector::value_type>>
        &values) const;

    /**
     * Same as above, but using a vector of local degree-of-freedom values. In
     * other words, instead of extracting the nodal values of the degrees of
     * freedom located on the current cell interface from a global vector
     * associated with a DoFHandler object (as the function above does), this
     * function instead takes these local nodal values through its first
     * argument.
     *
     * @param[in] here_or_there Same as the one in the above function.
     *
     * @param[in] local_dof_values A vector of local nodal values. This vector
     *   must have a length equal to number of DoFs on the current cell, and
     * must be ordered in the same order as degrees of freedom are numbered on
     *   the reference cell.
     *
     * @param[out] values A vector of values of the given finite element field,
     *   at the quadrature points on the current object.
     *
     * @tparam InputVector The @p InputVector type must allow creation
     *   of an ArrayView object from it; this is satisfied by the
     *   `std::vector` class, among others.
     */
    template <class InputVector>
    void
    get_function_values_from_local_dof_values(
      const bool         here_or_there,
      const InputVector &local_dof_values,
      std::vector<solution_value_type<typename InputVector::value_type>>
        &values) const;

    /** @} */

    /**
     * @name Access to jumps in global finite element fields
     */
    /** @{ */

    /**
     * Return the jump in the values of the selected scalar component of the
     * finite element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell interface selected the last time
     * the <tt>reinit</tt> function of the FEInterfaceValues object was called.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the values of shape functions (i.e., @p value_type) times the
     * type used to store the values of the unknowns $U_j$ of your finite
     * element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_values}
     */
    template <class InputVector>
    void
    get_jump_in_function_values(
      const InputVector &fe_function,
      std::vector<solution_value_type<typename InputVector::value_type>>
        &values) const;

    /**
     * This function relates to get_jump_in_function_values() in the same way
     * as get_function_values_from_local_dof_values() relates to
     * get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more
     * information.
     */
    template <class InputVector>
    void
    get_jump_in_function_values_from_local_dof_values(
      const InputVector &local_dof_values,
      std::vector<solution_value_type<typename InputVector::value_type>>
        &values) const;

    /**
     * Return the jump in the gradients of the selected scalar components of the
     * finite element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell interface selected the last time
     * the <tt>reinit</tt> function of the FEInterfaceValues object was called.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the gradients of shape functions (i.e., @p gradient_type)
     * times the type used to store the values of the unknowns $U_j$ of your
     * finite element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    template <class InputVector>
    void
    get_jump_in_function_gradients(
      const InputVector &fe_function,
      std::vector<solution_gradient_type<typename InputVector::value_type>>
        &gradients) const;

    /**
     * This function relates to get_jump_in_function_gradients() in the same way
     * as get_function_values_from_local_dof_values() relates to
     * get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more
     * information.
     */
    template <class InputVector>
    void
    get_jump_in_function_gradients_from_local_dof_values(
      const InputVector &local_dof_values,
      std::vector<solution_gradient_type<typename InputVector::value_type>>
        &gradients) const;

    /**
     * Return the jump in the Hessians of the selected scalar component of the
     * finite element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEInterfaceValues object was called.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the Hessians of shape functions (i.e., @p hessian_type) times
     * the type used to store the values of the unknowns $U_j$ of your finite
     * element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_hessians}
     */
    template <class InputVector>
    void
    get_jump_in_function_hessians(
      const InputVector &fe_function,
      std::vector<solution_hessian_type<typename InputVector::value_type>>
        &hessians) const;

    /**
     * This function relates to get_jump_in_function_hessians() in the same way
     * as get_function_values_from_local_dof_values() relates to
     * get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more
     * information.
     */
    template <class InputVector>
    void
    get_jump_in_function_hessians_from_local_dof_values(
      const InputVector &local_dof_values,
      std::vector<solution_hessian_type<typename InputVector::value_type>>
        &hessians) const;

    /**
     * Return the jump in the third derivatives of the selected scalar component
     * of the finite element function characterized by <tt>fe_function</tt> at
     * the quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEInterfaceValues object was called.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the third derivatives of shape functions (i.e., @p
     * third_derivative_type) times the type used to store the values of the
     * unknowns $U_j$ of your finite element vector $U$ (represented by the @p
     * fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_third_derivatives}
     */
    template <class InputVector>
    void
    get_jump_in_function_third_derivatives(
      const InputVector &fe_function,
      std::vector<
        solution_third_derivative_type<typename InputVector::value_type>>
        &third_derivatives) const;

    /**
     * This function relates to get_jump_in_function_third_derivatives() in the
     * same way as get_function_values_from_local_dof_values() relates
     * to get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more
     * information.
     */
    template <class InputVector>
    void
    get_jump_in_function_third_derivatives_from_local_dof_values(
      const InputVector &local_dof_values,
      std::vector<
        solution_third_derivative_type<typename InputVector::value_type>>
        &third_derivatives) const;

    /** @} */

    /**
     * @name Access to the average of global finite element fields
     */
    /** @{ */

    /**
     * Return the average of the values of the selected scalar component of the
     * finite element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell interface selected the last time
     * the <tt>reinit</tt> function of the FEInterfaceValues object was called.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the values of shape functions (i.e., @p value_type) times the
     * type used to store the values of the unknowns $U_j$ of your finite
     * element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_values}
     */
    template <class InputVector>
    void
    get_average_of_function_values(
      const InputVector &fe_function,
      std::vector<solution_value_type<typename InputVector::value_type>>
        &values) const;

    /**
     * This function relates to get_average_of_function_values() in the same way
     * as get_function_values_from_local_dof_values() relates to
     * get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more
     * information.
     */
    template <class InputVector>
    void
    get_average_of_function_values_from_local_dof_values(
      const InputVector &local_dof_values,
      std::vector<solution_value_type<typename InputVector::value_type>>
        &values) const;

    /**
     * Return the average of the gradients of the selected scalar components of
     * the finite element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell interface selected the last time
     * the <tt>reinit</tt> function of the FEInterfaceValues object was called.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the gradients of shape functions (i.e., @p gradient_type)
     * times the type used to store the values of the unknowns $U_j$ of your
     * finite element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    template <class InputVector>
    void
    get_average_of_function_gradients(
      const InputVector &fe_function,
      std::vector<solution_gradient_type<typename InputVector::value_type>>
        &gradients) const;

    /**
     * This function relates to get_average_of_function_gradients() in the same
     * way as get_function_values_from_local_dof_values() relates to
     * get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more
     * information.
     */
    template <class InputVector>
    void
    get_average_of_function_gradients_from_local_dof_values(
      const InputVector &local_dof_values,
      std::vector<solution_gradient_type<typename InputVector::value_type>>
        &gradients) const;

    /**
     * Return the average of the Hessians of the selected scalar component of
     * the finite element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEInterfaceValues object was called.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the Hessians of shape functions (i.e., @p hessian_type) times
     * the type used to store the values of the unknowns $U_j$ of your finite
     * element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_hessians}
     */
    template <class InputVector>
    void
    get_average_of_function_hessians(
      const InputVector &fe_function,
      std::vector<solution_hessian_type<typename InputVector::value_type>>
        &hessians) const;

    /**
     * This function relates to get_average_of_function_hessians() in the same
     * way as get_function_values_from_local_dof_values() relates to
     * get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more
     * information.
     */
    template <class InputVector>
    void
    get_average_of_function_hessians_from_local_dof_values(
      const InputVector &local_dof_values,
      std::vector<solution_hessian_type<typename InputVector::value_type>>
        &hessians) const;

    /** @} */

  private:
    /**
     * The extractor for this view.
     */
    const FEValuesExtractors::Scalar extractor;
  };



  /**
   * The view of a vector-valued variable for FEInterfaceValues.
   */
  template <int dim, int spacedim = dim>
  class Vector : public Base<dim, spacedim>
  {
  public:
    /**
     * This is the type returned for values.
     */
    using value_type =
      typename FEValuesViews::Vector<dim, spacedim>::value_type;

    /**
     * This is the type returned for gradients, for example from
     * average_of_gradients().
     */
    using gradient_type =
      typename FEValuesViews::Vector<dim, spacedim>::gradient_type;

    /**
     * An alias for the type of second derivatives of the view this class
     * represents. Here, for a set of <code>dim</code> components of the
     * finite element, the Hessian is a <code>Tensor@<3,dim@></code>.
     */
    using hessian_type =
      typename FEValuesViews::Vector<dim, spacedim>::hessian_type;

    /**
     * An alias for the type of third derivatives of the view this class
     * represents. Here, for a set of <code>dim</code> components of the
     * finite element, the third derivative is a <code>Tensor@<4,dim@></code>.
     */
    using third_derivative_type =
      typename FEValuesViews::Vector<dim, spacedim>::third_derivative_type;

    /**
     * An alias for the data type of the product of a @p Number and the
     * values of the view this class provides. This is the data type of
     * vector components of a finite element field whose degrees of
     * freedom are described by a vector with elements of type @p Number.
     */
    template <typename Number>
    using solution_value_type = typename ProductType<Number, value_type>::type;

    /**
     * An alias for the data type of the product of a @p Number and the
     * gradients of the view this class provides. This is the data type of
     * vector components of a finite element field whose degrees of
     * freedom are described by a vector with elements of type @p Number.
     */
    template <typename Number>
    using solution_gradient_type =
      typename ProductType<Number, gradient_type>::type;

    /**
     * An alias for the data type of the product of a @p Number and the
     * Hessians of the view this class provides. This is the data type of
     * vector components of a finite element field whose degrees of
     * freedom are described by a vector with elements of type @p Number.
     */
    template <typename Number>
    using solution_hessian_type =
      typename ProductType<Number, hessian_type>::type;

    /**
     * An alias for the data type of the product of a @p Number and the
     * third derivatives of the view this class provides. This is the data type
     * of vector components of a finite element field whose degrees of
     * freedom are described by a vector with elements of type @p Number.
     */
    template <typename Number>
    using solution_third_derivative_type =
      typename ProductType<Number, third_derivative_type>::type;

    /**
     * Constructor for an object that represents a vector component
     */
    Vector(const FEInterfaceValues<dim, spacedim> &fe_interface,
           const unsigned int                      first_vector_component);

    /**
     * @name Access to shape functions
     */
    /** @{ */

    /**
     * Return the value of the vector components selected by this view
     * with interface dof index @p interface_dof_index in
     * quadrature point @p q_point.
     *
     * The argument @p here_or_there selects between the upstream value and
     * the downstream value as defined by the direction of the normal vector
     * in this quadrature point. If @p here_or_there is true, the shape
     * functions from the first cell of the interface is used.
     *
     * In other words, this function returns the limit of the value of the shape
     * function in the given quadrature point when approaching it from one of
     * the two cells of the interface.
     *
     * @note This function is typically used to pick the upstream or downstream
     * value based on a direction. This can be achieved by using
     * <code>(direction * normal)>0</code> as the first argument of this
     * function.
     */
    value_type
    value(const bool         here_or_there,
          const unsigned int interface_dof_index,
          const unsigned int q_point) const;

    /** @} */

    /**
     * @name Access to jumps in shape functions and their derivatives
     */
    /** @{ */

    /**
     * Return the jump vector $[\mathbf{u}]=\mathbf{u_1} - \mathbf{u_2}$ on the
     * interface for the shape function
     * @p interface_dof_index in the quadrature point @p q_point.
     *
     * @note The name of the function is supposed to be read as "the jump
     *   (singular) in the values (plural: one or two possible values) of
     *   the shape function (singular)".
     */
    value_type
    jump_in_values(const unsigned int interface_dof_index,
                   const unsigned int q_point) const;

    /**
     * Return the jump of the gradient (a tensor of rank 2) $\jump{\nabla
     * \mathbf{u}}$ on the interface for the shape
     * function @p interface_dof_index in the quadrature point @p q_point.
     *
     * @note The name of the function is supposed to be read as "the jump
     *   (singular) in the gradients (plural: one or two possible gradients)
     *   of the shape function (singular)".
     */
    gradient_type
    jump_in_gradients(const unsigned int interface_dof_index,
                      const unsigned int q_point) const;

    /**
     * The same as above.
     *
     * @deprecated Use the average_of_gradients() function instead.
     */
    gradient_type
    jump_gradient(const unsigned int interface_dof_index,
                  const unsigned int q_point) const;

    /**
     * Return the jump in the gradient $\jump{\nabla u}=\nabla u_{\text{cell0}}
     * - \nabla u_{\text{cell1}}$ on the interface for the shape function @p
     * interface_dof_index at the quadrature point @p q_point of
     * the component selected by this view.
     *
     * @note The name of the function is supposed to be read as "the jump
     *   (singular) in the Hessians (plural: one or two possible values
     *   for the second derivative) of the shape function (singular)".
     */
    hessian_type
    jump_in_hessians(const unsigned int interface_dof_index,
                     const unsigned int q_point) const;

    /**
     * The same as above.
     *
     * @deprecated Use the average_of_hessians() function instead.
     */
    hessian_type
    jump_hessian(const unsigned int interface_dof_index,
                 const unsigned int q_point) const;

    /**
     * Return the jump in the third derivative $\jump{\nabla^3 u} = \nabla^3
     * u_{\text{cell0}} - \nabla^3 u_{\text{cell1}}$ on the interface for the
     * shape function @p interface_dof_index at the quadrature point @p q_point of
     * the component selected by this view.
     *
     * @note The name of the function is supposed to be read as "the jump
     *   (singular) in the third derivatives (plural: one or two possible values
     *   for the third derivative) of the shape function (singular)".
     */
    third_derivative_type
    jump_in_third_derivatives(const unsigned int interface_dof_index,
                              const unsigned int q_point) const;

    /** @} */

    /**
     * @name Access to the average of shape functions and their derivatives
     */
    /** @{ */

    /**
     * Return the average vector $\average{\mathbf{u}}=\frac{1}{2}(\mathbf{u_1}
     * + \mathbf{u_2})$ on the interface for the shape
     * function @p interface_dof_index in the quadrature point @p q_point.
     *
     * @note The name of the function is supposed to be read as "the average
     *   (singular) of the values (plural: one or two possible values) of
     *   the shape function (singular)".
     */
    value_type
    average_of_values(const unsigned int interface_dof_index,
                      const unsigned int q_point) const;

    /**
     * Return the average of the gradient (a tensor of rank 2) $\average{\nabla
     * \mathbf{u}}$ on the interface for the shape
     * function @p interface_dof_index in the quadrature point @p q_point.
     *
     * @note The name of the function is supposed to be read as "the average
     *   (singular) of the gradients (plural: one or two possible values
     *   of the derivative) of the shape function (singular)".
     */
    gradient_type
    average_of_gradients(const unsigned int interface_dof_index,
                         const unsigned int q_point) const;

    /**
     * Return the average of the Hessian $\average{\nabla^2 u} =
     * \frac{1}{2}\nabla^2 u_{\text{cell0}} + \frac{1}{2} \nabla^2
     * u_{\text{cell1}}$ on the interface
     * for the shape function @p interface_dof_index at the quadrature point @p
     * q_point of the component selected by this view.
     *
     * @note The name of the function is supposed to be read as "the average
     *   (singular) in the Hessians (plural: one or two possible values of
     *   the second derivative) of the shape function (singular)".
     */
    hessian_type
    average_of_hessians(const unsigned int interface_dof_index,
                        const unsigned int q_point) const;

    /**
     * The same as above.
     *
     * @deprecated Use the average_of_hessians() function instead.
     */
    hessian_type
    average_hessian(const unsigned int interface_dof_index,
                    const unsigned int q_point) const;

    /** @} */

    /**
     * @name Access to values of global finite element fields
     */
    /** @{ */

    /**
     * Return the values of the selected vector component of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell interface selected the last time
     * the <tt>reinit</tt> function of the FEInterfaceValues object was called.
     *
     * The argument @p here_or_there selects between the value on cell 0 (here, @p true)
     * and cell 1 (there, @p false). You can also interpret it as "upstream" (@p true)
     * and "downstream" (@p false) as defined by the direction of the normal
     * vector in this quadrature point. If @p here_or_there is true, the values
     * from the first cell of the interface is used.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the values of shape functions (i.e., @p value_type) times the
     * type used to store the values of the unknowns $U_j$ of your finite
     * element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_values}
     */
    template <class InputVector>
    void
    get_function_values(
      const bool         here_or_there,
      const InputVector &fe_function,
      std::vector<solution_value_type<typename InputVector::value_type>>
        &values) const;

    /**
     * Same as above, but using a vector of local degree-of-freedom values. In
     * other words, instead of extracting the nodal values of the degrees of
     * freedom located on the current cell interface from a global vector
     * associated with a DoFHandler object (as the function above does), this
     * function instead takes these local nodal values through its first
     * argument.
     *
     * @param[in] here_or_there Same as the one in the above function.
     *
     * @param[in] local_dof_values A vector of local nodal values. This vector
     *   must have a length equal to number of DoFs on the current cell, and
     * must be ordered in the same order as degrees of freedom are numbered on
     *   the reference cell.
     *
     * @param[out] values A vector of values of the given finite element field,
     *   at the quadrature points on the current object.
     *
     * @tparam InputVector The @p InputVector type must allow creation
     *   of an ArrayView object from it; this is satisfied by the
     *   `std::vector` class, among others.
     */
    template <class InputVector>
    void
    get_function_values_from_local_dof_values(
      const bool         here_or_there,
      const InputVector &local_dof_values,
      std::vector<solution_value_type<typename InputVector::value_type>>
        &values) const;

    /** @} */

    /**
     * @name Access to jumps in global finite element fields
     */
    /** @{ */

    /**
     * Return the jump in the values of the selected vector component of the
     * finite element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell interface selected the last time
     * the <tt>reinit</tt> function of the FEInterfaceValues object was called.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the values of shape functions (i.e., @p value_type) times the
     * type used to store the values of the unknowns $U_j$ of your finite
     * element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_values}
     */
    template <class InputVector>
    void
    get_jump_in_function_values(
      const InputVector &fe_function,
      std::vector<solution_value_type<typename InputVector::value_type>>
        &values) const;

    /**
     * This function relates to get_jump_in_function_values() in the same way
     * as get_function_values_from_local_dof_values() relates to
     * get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more
     * information.
     */
    template <class InputVector>
    void
    get_jump_in_function_values_from_local_dof_values(
      const InputVector &local_dof_values,
      std::vector<solution_value_type<typename InputVector::value_type>>
        &values) const;

    /**
     * Return the jump in the gradients of the selected vector components of the
     * finite element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell interface selected the last time
     * the <tt>reinit</tt> function of the FEInterfaceValues object was called.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the gradients of shape functions (i.e., @p gradient_type)
     * times the type used to store the values of the unknowns $U_j$ of your
     * finite element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    template <class InputVector>
    void
    get_jump_in_function_gradients(
      const InputVector &fe_function,
      std::vector<solution_gradient_type<typename InputVector::value_type>>
        &gradients) const;

    /**
     * This function relates to get_jump_in_function_gradients() in the same way
     * as get_function_values_from_local_dof_values() relates to
     * get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more
     * information.
     */
    template <class InputVector>
    void
    get_jump_in_function_gradients_from_local_dof_values(
      const InputVector &local_dof_values,
      std::vector<solution_gradient_type<typename InputVector::value_type>>
        &gradients) const;

    /**
     * Return the jump in the Hessians of the selected vector component of the
     * finite element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEInterfaceValues object was called.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the Hessians of shape functions (i.e., @p hessian_type) times
     * the type used to store the values of the unknowns $U_j$ of your finite
     * element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_hessians}
     */
    template <class InputVector>
    void
    get_jump_in_function_hessians(
      const InputVector &fe_function,
      std::vector<solution_hessian_type<typename InputVector::value_type>>
        &hessians) const;

    /**
     * This function relates to get_jump_in_function_hessians() in the same way
     * as get_function_values_from_local_dof_values() relates to
     * get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more
     * information.
     */
    template <class InputVector>
    void
    get_jump_in_function_hessians_from_local_dof_values(
      const InputVector &local_dof_values,
      std::vector<solution_hessian_type<typename InputVector::value_type>>
        &hessians) const;

    /**
     * Return the jump in the third derivatives of the selected vector component
     * of the finite element function characterized by <tt>fe_function</tt> at
     * the quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEInterfaceValues object was called.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the third derivatives of shape functions (i.e., @p
     * third_derivative_type) times the type used to store the values of the
     * unknowns $U_j$ of your finite element vector $U$ (represented by the @p
     * fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_third_derivatives}
     */
    template <class InputVector>
    void
    get_jump_in_function_third_derivatives(
      const InputVector &fe_function,
      std::vector<
        solution_third_derivative_type<typename InputVector::value_type>>
        &third_derivatives) const;

    /**
     * This function relates to get_jump_in_function_third_derivatives() in the
     * same way as get_function_values_from_local_dof_values() relates
     * to get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more
     * information.
     */
    template <class InputVector>
    void
    get_jump_in_function_third_derivatives_from_local_dof_values(
      const InputVector &local_dof_values,
      std::vector<
        solution_third_derivative_type<typename InputVector::value_type>>
        &third_derivatives) const;

    /** @} */

    /**
     * @name Access to the average of global finite element fields
     */
    /** @{ */

    /**
     * Return the average of the values of the selected vector component of the
     * finite element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell interface selected the last time
     * the <tt>reinit</tt> function of the FEInterfaceValues object was called.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the values of shape functions (i.e., @p value_type) times the
     * type used to store the values of the unknowns $U_j$ of your finite
     * element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_values}
     */
    template <class InputVector>
    void
    get_average_of_function_values(
      const InputVector &fe_function,
      std::vector<solution_value_type<typename InputVector::value_type>>
        &values) const;

    /**
     * This function relates to get_average_of_function_values() in the same way
     * as get_function_values_from_local_dof_values() relates to
     * get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more
     * information.
     */
    template <class InputVector>
    void
    get_average_of_function_values_from_local_dof_values(
      const InputVector &local_dof_values,
      std::vector<solution_value_type<typename InputVector::value_type>>
        &values) const;

    /**
     * Return the average of the gradients of the selected vector components of
     * the finite element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell interface selected the last time
     * the <tt>reinit</tt> function of the FEInterfaceValues object was called.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the gradients of shape functions (i.e., @p gradient_type)
     * times the type used to store the values of the unknowns $U_j$ of your
     * finite element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    template <class InputVector>
    void
    get_average_of_function_gradients(
      const InputVector &fe_function,
      std::vector<solution_gradient_type<typename InputVector::value_type>>
        &gradients) const;

    /**
     * This function relates to get_average_of_function_gradients() in the same
     * way as get_function_values_from_local_dof_values() relates to
     * get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more
     * information.
     */
    template <class InputVector>
    void
    get_average_of_function_gradients_from_local_dof_values(
      const InputVector &local_dof_values,
      std::vector<solution_gradient_type<typename InputVector::value_type>>
        &gradients) const;

    /**
     * Return the average of the Hessians of the selected vector component of
     * the finite element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEInterfaceValues object was called.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the Hessians of shape functions (i.e., @p hessian_type) times
     * the type used to store the values of the unknowns $U_j$ of your finite
     * element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_hessians}
     */
    template <class InputVector>
    void
    get_average_of_function_hessians(
      const InputVector &fe_function,
      std::vector<solution_hessian_type<typename InputVector::value_type>>
        &hessians) const;

    /**
     * This function relates to get_average_of_function_hessians() in the same
     * way as get_function_values_from_local_dof_values() relates to
     * get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more
     * information.
     */
    template <class InputVector>
    void
    get_average_of_function_hessians_from_local_dof_values(
      const InputVector &local_dof_values,
      std::vector<solution_hessian_type<typename InputVector::value_type>>
        &hessians) const;

    /** @} */

  private:
    /**
     * The extractor for this view.
     */
    const FEValuesExtractors::Vector extractor;
  };
} // namespace FEInterfaceViews


namespace internal
{
  namespace FEInterfaceViews
  {
    /**
     * A class whose specialization is used to define what FEInterfaceViews
     * object corresponds to the given FEValuesExtractors object.
     */
    template <int dim, int spacedim, typename Extractor>
    struct ViewType
    {};

    /**
     * A class whose specialization is used to define what FEInterfaceViews
     * object corresponds to the given FEValuesExtractors object.
     *
     * When using FEValuesExtractors::Scalar, the corresponding view is an
     * FEInterfaceViews::Scalar<dim, spacedim>.
     */
    template <int dim, int spacedim>
    struct ViewType<dim, spacedim, FEValuesExtractors::Scalar>
    {
      using type = typename dealii::FEInterfaceViews::Scalar<dim, spacedim>;
    };

    /**
     * A class whose specialization is used to define what FEInterfaceViews
     * object corresponds to the given FEValuesExtractors object.
     *
     * When using FEValuesExtractors::Vector, the corresponding view is an
     * FEInterfaceViews::Vector<dim, spacedim>.
     */
    template <int dim, int spacedim>
    struct ViewType<dim, spacedim, FEValuesExtractors::Vector>
    {
      using type = typename dealii::FEInterfaceViews::Vector<dim, spacedim>;
    };
  } // namespace FEInterfaceViews
} // namespace internal

namespace FEInterfaceViews
{
  /**
   * A templated alias that associates to a given Extractor class
   * the corresponding view in FEInterfaceViews.
   */
  template <int dim, int spacedim, typename Extractor>
  using View = typename dealii::internal::FEInterfaceViews::
    ViewType<dim, spacedim, Extractor>::type;
} // namespace FEInterfaceViews



/**
 * FEInterfaceValues is a data structure to access and assemble finite element
 * data on interfaces between two cells of a mesh.
 *
 * It provides a way to access averages, jump terms, and similar operations used
 * in Discontinuous Galerkin methods on a face between two neighboring cells.
 * This allows the computation of typical mesh-dependent linear or bilinear
 * forms in a similar way as FEValues does for cells and FEFaceValues does for
 * faces. In
 * the literature, the faces between neighboring cells are called "inner
 * interfaces" or "facets".
 *
 * Internally, this class provides an abstraction for two FEFaceValues
 * objects (or FESubfaceValues when using adaptive refinement). The class
 * introduces a new "interface dof index" that walks over
 * the union of the dof indices of the two FEFaceValues objects. Helper
 * functions allow translating between the new "interface dof index" and the
 * corresponding "cell index" (0 for the first cell, 1 for the second cell)
 * and "dof index" within that cell.
 *
 * The class is made to be used inside MeshWorker::mesh_loop(). It is intended
 * to be a low level replacement for MeshWorker and LocalIntegrators and a
 * higher level abstraction compared to assembling face terms manually.
 */
template <int dim, int spacedim = dim>
class FEInterfaceValues
{
public:
  /**
   * Number of quadrature points.
   */
  const unsigned int n_quadrature_points;

  /**
   * Construct the FEInterfaceValues with a single FiniteElement (same on both
   * sides of the facet). The FEFaceValues objects will be initialized with
   * the given @p mapping, @p quadrature, and @p update_flags.
   */
  FEInterfaceValues(const Mapping<dim, spacedim>       &mapping,
                    const FiniteElement<dim, spacedim> &fe,
                    const Quadrature<dim - 1>          &quadrature,
                    const UpdateFlags                   update_flags);

  /**
   * The same as above but taking a collection of quadrature rules
   * so that different quadrature rules can be assigned to different
   * faces.
   */
  FEInterfaceValues(const Mapping<dim, spacedim>       &mapping,
                    const FiniteElement<dim, spacedim> &fe,
                    const hp::QCollection<dim - 1>     &quadrature,
                    const UpdateFlags                   update_flags);

  /**
   * Construct the FEInterfaceValues with a single FiniteElement and
   * a Q1 Mapping.
   *
   * See the constructor above.
   */
  FEInterfaceValues(const FiniteElement<dim, spacedim> &fe,
                    const Quadrature<dim - 1>          &quadrature,
                    const UpdateFlags                   update_flags);

  /**
   * Construct the FEInterfaceValues object with different FiniteElements
   * assigned to different faces.
   */
  FEInterfaceValues(
    const hp::MappingCollection<dim, spacedim> &mapping_collection,
    const hp::FECollection<dim, spacedim>      &fe_collection,
    const hp::QCollection<dim - 1>             &quadrature_collection,
    const UpdateFlags                           update_flags);

  /**
   * Same as above, but using the default linear mapping.
   */
  FEInterfaceValues(const hp::FECollection<dim, spacedim> &fe_collection,
                    const hp::QCollection<dim - 1> &quadrature_collection,
                    const UpdateFlags               update_flags);

  /**
   * Re-initialize this object to be used on a new interface given by two faces
   * of two neighboring cells. The `cell` and `cell_neighbor` cells will be
   * referred to through `cell_index` zero and one after this call in all places
   * where one needs to identify the two cells adjacent to the interface.
   *
   * Use numbers::invalid_unsigned_int for @p sub_face_no or @p
   * sub_face_no_neighbor to indicate that you want to work on the entire face,
   * not a sub-face.
   *
   * The arguments (including their order) are identical to the @p face_worker
   * arguments in MeshWorker::mesh_loop().
   *
   * In order to do integration on a face or sub-face, this object will have to
   * choose what quadrature formula to use. This is simple if you initialized
   * the current FEInterfaceValues object with just a single FiniteElement and
   * a single (face) Quadrature object, because in that case there is only one
   * element and quadrature that will be used for all interfaces. But it is not
   * so simple in the hp case where there may be different elements used on
   * different cells, and different quadrature formulas should be used on
   * different faces; one may also want to use a different mapping for different
   * faces. As a consequence, you would have initialized the current object
   * with a hp::FECollection, hp::QCollection, and possible an
   * hp::MappingCollection object. The question then is: For a given face or
   * subface, which quadrature and mapping should be used? The following
   * decision tree will then be used:
   * -# If the `q_index` and `mapping_index` arguments to this function are
   *    explicitly specified (rather than leaving them at their default
   *    values), then these indices will be used to select which element
   *    of the hp::QCollection and hp::MappingCollection passed to the
   *    constructor should serve as the quadrature and mapping to be used.
   * -# If one of these arguments is left at its default value, then the
   *    function will need to choose a quadrature and/or mapping that is
   *    appropriate for the two finite element spaces used on the two cells
   *    adjacent to the current interface. As the first choice, if the
   *    quadrature or mapping collection we are considering has only one
   *    element, then that is clearly the one that should be used.
   * -# If the quadrature or mapping collection have multiple elements,
   *    then we need to dig further. For quadrature objects, we can
   *    compare whether the two quadrature objects that correspond to
   *    the `active_fe_index` values of the two adjacent cells are
   *    identical (i.e., have quadrature points at the same locations,
   *    and have the same weights). If this is so, then it does not
   *    matter which one of the two we take, and we choose one or the
   *    other.
   * -# If this has still not helped, we try to find out
   *    which of the two finite element spaces on the two adjacent cells
   *    is "larger" (say, if you had used $Q_2$
   *    and $Q_4$ elements on the two adjacent cells, then the $Q_4$
   *    element is the larger one); the determination of which space
   *    is "larger" is made using the hp::FECollection::find_dominated_fe()
   *    function, which is not necessarily intended for this kind of query,
   *    but yields a result that serves just fine for our purposes here.
   *    We then operate on the assumption that the quadrature object
   *    associated with the "larger" of the two spaces is the appropriate
   *    one to use for the face that separates these two spaces.
   *    - If this function returns that one of the two elements in question
   *      is dominated by the other, then presumably it is "larger" one and
   *      we take the quadrature formula and mapping that corresponds
   *      to this "larger" element is. For example, for the $Q_2$ element
   *      mentioned above, one would generally use a QGauss(3) quadrature
   *      formula, whereas for the $Q_4$ element, one would use QGauss(5).
   *      To integrate jump and average terms on the interface between
   *      cells using these two elements, QGauss(5) is appropriate. Because,
   *      typically, people will order elements
   *      in the hp::FECollection in the same order as the quadrature and
   *      mapping objects in hp::QCollection and hp::MappingCollection, this
   *      function will use the index of the "larger" element in the
   *      hp::FECollection to also index into the hp::QCollection and
   *      hp::MappingCollection to retrieve quadrature and mapping objects
   *      appropriate for the current face.
   *    - There are cases where neither element dominates the other. For
   *      example, if one uses $Q_2\times Q_1$ and $Q_1\times Q_2$ elements on
   *      neighboring cells, neither of the two spaces dominates the other --
   *      or, in the context of the current function, neither space is "larger"
   *      than the other. In that case, there is no way for the current function
   *      to determine quadrature and mapping objects associated with the two
   *      elements are the appropriate ones. If that happens, you will get an
   *      error -- and the only way to avoid the error is to explicitly specify
   *      for these interfaces which quadrature and mapping objects you want
   *      to use, by providing non-default values for the @p q_index and
   *      @p mapping_index arguments to this function.
   *
   * @param[in] cell An iterator to the first cell adjacent to the interface.
   * @param[in] face_no An integer identifying which face of the first cell the
   *   interface is on.
   * @param[in] sub_face_no An integer identifying the subface (child) of the
   *   face (identified by the previous two arguments) that the interface
   *   corresponds to. If equal to numbers::invalid_unsigned_int, then the
   *   interface is considered to be the entire face.
   * @param[in] cell_neighbor An iterator to the second cell adjacent to
   *   the interface. The type of this iterator does not have to equal that
   *   of `cell`, but must be convertible to it. This allows using an
   *   active cell iterator for `cell`, and `cell->neighbor(f)` for
   *   `cell_neighbor`, since the return type of `cell->neighbor(f)` is
   *   simply a cell iterator (not necessarily an active cell iterator).
   * @param[in] face_no_neighbor Like `face_no`, just for the neighboring
   *   cell.
   * @param[in] sub_face_no_neighbor Like `sub_face_no`, just for the
   *   neighboring cell.
   * @param[in] q_index The index of the quadrature object within the
   *   hp::QCollection passed to the constructor to use on the current
   *   interface. See the documentation above what happens if this argument
   *   is not explicitly provided but left at its default value.
   * @param[in] mapping_index The index of the mapping object within the
   *   hp::MappingCollection passed to the constructor to use on the current
   *   interface. See the documentation above what happens if this argument
   *   is not explicitly provided but left at its default value.
   * @param[in] fe_index If left at its default, use the finite element within
   *   the hp::FECollection passed to the constructor as given by the dominating
   *   finite element across the interface (only used if the FEInterface object
   *   is initialized with an hp::FECollection, an hp::QCollection, and possibly
   *   an hp::MappingCollection).
   * @param[in] fe_index_neighbor Active fe index of neighboring cell. Useful
   *   if hp capabilities are used and for non-DoFHandler iterators.
   */
  template <typename CellIteratorType, typename CellNeighborIteratorType>
  void
  reinit(const CellIteratorType         &cell,
         const unsigned int              face_no,
         const unsigned int              sub_face_no,
         const CellNeighborIteratorType &cell_neighbor,
         const unsigned int              face_no_neighbor,
         const unsigned int              sub_face_no_neighbor,
         const unsigned int q_index           = numbers::invalid_unsigned_int,
         const unsigned int mapping_index     = numbers::invalid_unsigned_int,
         const unsigned int fe_index          = numbers::invalid_unsigned_int,
         const unsigned int fe_index_neighbor = numbers::invalid_unsigned_int);

  /**
   * Re-initialize this object to be used on an interface given by a single face
   * @p face_no of the cell @p cell. This is useful to use FEInterfaceValues
   * on boundaries of the domain.
   *
   * As a consequence, members like jump() will assume a value of zero for the
   * values on the "other" side. Note that no sub_face_number is needed as a
   * boundary face can not neighbor a finer cell.
   *
   * After calling this function at_boundary() will return true.
   *
   * @param[in] cell An iterator to the first cell adjacent to the interface.
   * @param[in] face_no An integer identifying which face of the first cell the
   *   interface is on.
   * @param[in] q_index This argument selects which quadrature formula to use
   *   See the discussion in the documentation of the other reinit() function
   *   for what happens when this argument is left at its default value.
   * @param[in] mapping_index This argument selects which mapping to use
   *   See the discussion in the documentation of the other reinit() function
   *   for what happens when this argument is left at its default value.
   * @param[in] fe_index If left at its default, use the finite element within
   * the hp::FECollection passed to the constructor as given by the dominating
   * finite element across the interface (only used if the FEInterface object
   * is initialized with an hp::FECollection, an hp::QCollection, and possibly
   * an hp::MappingCollection).
   */
  template <typename CellIteratorType>
  void
  reinit(const CellIteratorType &cell,
         const unsigned int      face_no,
         const unsigned int      q_index       = numbers::invalid_unsigned_int,
         const unsigned int      mapping_index = numbers::invalid_unsigned_int,
         const unsigned int      fe_index      = numbers::invalid_unsigned_int);

  /**
   * Return a reference to the FEFaceValues or FESubfaceValues object
   * of the specified cell of the interface.
   *
   * The @p cell_index is either 0 or 1 and corresponds to the cell index
   * returned by interface_dof_to_cell_and_dof_index().
   */
  const FEFaceValuesBase<dim, spacedim> &
  get_fe_face_values(const unsigned int cell_index) const;

  /**
   * Constant reference to the selected mapping object.
   */
  const Mapping<dim, spacedim> &
  get_mapping() const;

  /**
   * Constant reference to the selected finite element object.
   */
  const FiniteElement<dim, spacedim> &
  get_fe() const;

  /**
   * Return a reference to the quadrature object in use.
   */
  const Quadrature<dim - 1> &
  get_quadrature() const;

  /**
   * Return a reference to the used mapping.
   */
  const hp::MappingCollection<dim, spacedim> &
  get_mapping_collection() const;

  /**
   * Return a reference to the selected finite element object.
   */
  const hp::FECollection<dim, spacedim> &
  get_fe_collection() const;

  /**
   * Return a reference to the face quadrature object in use.
   */
  const hp::QCollection<dim - 1> &
  get_quadrature_collection() const;

  /**
   * Returns a boolean indicating whether or not this FEInterfaceValues object
   * has hp-capabilities enabled.
   */
  bool
  has_hp_capabilities() const;

  /**
   * Return the update flags set.
   */
  UpdateFlags
  get_update_flags() const;

  /**
   * Return a triangulation iterator to the current cell of the interface.
   *
   * The @p cell_index is either 0 or 1 and corresponds to the cell index
   * returned by interface_dof_to_cell_and_dof_index().
   */
  typename Triangulation<dim, spacedim>::cell_iterator
  get_cell(const unsigned int cell_index) const;

  /**
   * Return the number of the face on the interface selected the last time
   * the reinit() function was called.
   *
   * The @p cell_index is either 0 or 1 and corresponds to the cell index
   * returned by interface_dof_to_cell_and_dof_index().
   */
  unsigned int
  get_face_number(const unsigned int cell_index) const;

  /**
   * @name Functions to query information on a given interface
   * @{
   */

  /**
   * Return if the current interface is a boundary face or an internal
   * face with two adjacent cells.
   *
   * See the corresponding reinit() functions for details.
   */
  bool
  at_boundary() const;

  /**
   * Mapped quadrature weight. This value equals the
   * mapped surface element times the weight of the quadrature
   * point.
   *
   * You can think of the quantity returned by this function as the
   * surface element $ds$ in the integral that we implement here by
   * quadrature.
   *
   * @dealiiRequiresUpdateFlags{update_JxW_values}
   */
  double
  JxW(const unsigned int quadrature_point) const;

  /**
   * Return the vector of JxW values for each quadrature point.
   *
   * @dealiiRequiresUpdateFlags{update_JxW_values}
   */
  const std::vector<double> &
  get_JxW_values() const;

  /**
   * Return the normal in a given quadrature point.
   *
   * The normal points in outwards direction as seen from the first cell of
   * this interface.
   *
   * @deprecated Use the function normal_vector().
   *
   * @dealiiRequiresUpdateFlags{update_normal_vectors}
   */
  DEAL_II_DEPRECATED_EARLY_WITH_COMMENT("Use the function normal_vector().")
  Tensor<1, spacedim>
  normal(const unsigned int q_point_index) const;

  /**
   * Return the normal in a given quadrature point.
   *
   * The normal points in outwards direction as seen from the first cell of
   * this interface.
   *
   * @dealiiRequiresUpdateFlags{update_normal_vectors}
   */
  Tensor<1, spacedim>
  normal_vector(const unsigned int q_point_index) const;

  /**
   * Return the normal vector of the interface in each quadrature point.
   *
   * The return value is identical to get_fe_face_values(0).get_normal_vectors()
   * and therefore, are outside normal vectors from the perspective of the
   * first cell of this interface.
   *
   * @dealiiRequiresUpdateFlags{update_normal_vectors}
   */
  const std::vector<Tensor<1, spacedim>> &
  get_normal_vectors() const;

  /**
   * Return an object that can be thought of as an array containing all
   * indices from zero to `n_quadrature_points`. This allows to write code
   * using range-based `for` loops.
   *
   * @see CPP11
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  quadrature_point_indices() const;

  /**
   * Return the location of the <tt>q_point</tt>th quadrature point in
   * real space.
   *
   * @dealiiRequiresUpdateFlags{update_quadrature_points}
   */
  const Point<spacedim> &
  quadrature_point(const unsigned int q_point) const;

  /**
   * Return a reference to the quadrature points in real space.
   *
   * @dealiiRequiresUpdateFlags{update_quadrature_points}
   */
  const std::vector<Point<spacedim>> &
  get_quadrature_points() const;

  /**
   * Return the number of DoFs (or shape functions) on the current interface.
   *
   * @note This number is only available after a call to reinit() and can change
   * from one call to reinit() to the next. For example, on a boundary interface
   * it is equal to the number of dofs of the single FEFaceValues object, while
   * it is twice that for an interior interface for a DG element. For a
   * continuous element, it is slightly smaller because the two cells on the
   * interface share some of the dofs.
   */
  unsigned
  n_current_interface_dofs() const;

  /**
   * Return an object that can be thought of as an array containing all
   * indices from zero (inclusive) to `n_current_interface_dofs()` (exclusive).
   * This allows one to write code using range-based `for` loops of the
   * following kind:
   * @code
   *   FEInterfaceValues<dim> fe_iv (...);
   *   FullMatrix<double>     cell_matrix (...);
   *
   *   for (auto &cell : dof_handler.active_cell_iterators())
   *     {
   *       cell_matrix = 0;
   *       for (const auto &face : cell->face_iterators())
   *         {
   *           fe_iv.values.reinit(cell, face, ...);
   *           for (const auto q : fe_iv.quadrature_point_indices())
   *             for (const auto i : fe_iv.dof_indices())
   *               for (const auto j : fe_iv.dof_indices())
   *                 cell_matrix(i,j) += ...; // Do something for DoF indices
   *                                          // (i,j) at quadrature point q
   *         }
   *     }
   * @endcode
   * Here, we are looping over all degrees of freedom on all cell interfaces,
   * with `i` and `j` taking on all valid indices for interface degrees of
   * freedom, as defined by the finite element passed to `fe_iv`.
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  dof_indices() const;

  /**
   * Return the set of joint DoF indices. This includes indices from both cells.
   * If reinit was called with an active cell iterator, the indices are based
   * on the active indices (returned by `DoFCellAccessor::get_dof_indices()` ),
   * in case of level cell (that is, if is_level_cell() return true )
   * the mg dof indices are returned.
   *
   * @note This function is only available after a call to reinit() and can
   * change from one call to reinit() to the next.
   */
  std::vector<types::global_dof_index>
  get_interface_dof_indices() const;

  /**
   * Convert an interface dof index into the corresponding local DoF indices of
   * the two cells. If an interface DoF is only active on one of the
   * cells, the other index will be numbers::invalid_unsigned_int.
   *
   * For discontinuous finite elements, each interface dof is located on exactly
   * one side of the interface and, consequently, only one of the two values
   * returned is valid (i.e., different from numbers::invalid_unsigned_int).
   *
   * @note This function is only available after a call to reinit() and the
   * returned values may change from one call to reinit() to the next.
   */
  std::array<unsigned int, 2>
  interface_dof_to_dof_indices(const unsigned int interface_dof_index) const;

  /**
   * @}
   */

  /**
   * @name Access to shape functions and their derivatives
   * @{
   */

  /**
   * Return component @p component of the value of the shape function
   * with interface dof index @p interface_dof_index in
   * quadrature point @p q_point.
   *
   * The argument @p here_or_there selects between the value on cell 0 (here, @p true)
   * and cell 1 (there, @p false). You can also interpret it as "upstream" (@p true)
   * and "downstream" (@p false) as defined by the direction of the normal
   * vector
   * in this quadrature point. If @p here_or_there is true, the shape
   * functions from the first cell of the interface is used.
   *
   * In other words, this function returns the limit of the value of the shape
   * function in the given quadrature point when approaching it from one of the
   * two cells of the interface.
   *
   * @note This function is typically used to pick the upstream or downstream
   * value based on a direction. This can be achieved by using
   * <code>(direction * normal)>0</code> as the first argument of this
   * function.
   */
  double
  shape_value(const bool         here_or_there,
              const unsigned int interface_dof_index,
              const unsigned int q_point,
              const unsigned int component = 0) const;

  /**
   * Return component @p component of the gradient of the shape function
   * with interface dof index @p interface_dof_index in
   * quadrature point @p q_point.
   *
   * Similar to @p shape_value() , the argument @p here_or_there selects between
   * the gradient on cell 0 (here, @p true) and cell 1 (there, @p false)
   * at this quadrature point.
   */
  Tensor<1, spacedim>
  shape_grad(const bool         here_or_there,
             const unsigned int interface_dof_index,
             const unsigned int q_point,
             const unsigned int component = 0) const;

  /**
   * @}
   */

  /**
   * @name Access to jumps in shape function values and their derivatives
   * @{
   */

  /**
   * Return the jump $\jump{u}=u_{\text{cell0}} - u_{\text{cell1}}$ on the
   * interface
   * for the shape function @p interface_dof_index at the quadrature point
   * @p q_point of component @p component.
   *
   * Note that one can define the jump in
   * different ways (the value "there" minus the value "here", or the other way
   * around; both are used in the finite element literature). The definition
   * here uses "value here minus value there", as seen from the first cell.
   *
   * If this is a boundary face (at_boundary() returns true), then
   * $\jump{u}=u_{\text{cell0}}$, that is "the value here (minus zero)".
   *
   * @note The name of the function is supposed to be read as "the jump
   *   (singular) in the values (plural: one or two possible values)
   *   of the shape function (singular)".
   */
  double
  jump_in_shape_values(const unsigned int interface_dof_index,
                       const unsigned int q_point,
                       const unsigned int component = 0) const;

  /**
   * Return the jump in the gradient $\jump{\nabla u}=\nabla u_{\text{cell0}} -
   * \nabla u_{\text{cell1}}$ on the interface for the shape function @p
   * interface_dof_index at the quadrature point @p q_point of component @p
   * component.
   *
   * If this is a boundary face (at_boundary() returns true), then
   * $\jump{\nabla u}=\nabla u_{\text{cell0}}$.
   *
   * @note The name of the function is supposed to be read as "the jump
   *   (singular) in the gradients (plural: one or two possible gradients)
   *   of the shape function (singular)".
   */
  Tensor<1, spacedim>
  jump_in_shape_gradients(const unsigned int interface_dof_index,
                          const unsigned int q_point,
                          const unsigned int component = 0) const;

  /**
   * Return the jump in the Hessian $\jump{\nabla^2 u} = \nabla^2
   * u_{\text{cell0}} - \nabla^2 u_{\text{cell1}}$ on the interface for the
   * shape function
   * @p interface_dof_index at the quadrature point @p q_point of component
   * @p component.
   *
   * If this is a boundary face (at_boundary() returns true), then
   * $\jump{\nabla^2 u} = \nabla^2 u_{\text{cell0}}$.
   *
   * @note The name of the function is supposed to be read as "the jump
   *   (singular) in the Hessians (plural: one or two possible values
   *   for the derivative) of the shape function (singular)".
   */
  Tensor<2, spacedim>
  jump_in_shape_hessians(const unsigned int interface_dof_index,
                         const unsigned int q_point,
                         const unsigned int component = 0) const;

  /**
   * Return the jump in the third derivative $\jump{\nabla^3 u} = \nabla^3
   * u_{\text{cell0}} - \nabla^3 u_{\text{cell1}}$ on the interface for the
   * shape function @p interface_dof_index at the quadrature point @p q_point of
   * component @p component.
   *
   * If this is a boundary face (at_boundary() returns true), then
   * $\jump{\nabla^3 u} = \nabla^3 u_{\text{cell0}}$.
   *
   * @note The name of the function is supposed to be read as "the jump
   *   (singular) in the third derivatives (plural: one or two possible values
   *   for the derivative) of the shape function (singular)".
   */
  Tensor<3, spacedim>
  jump_in_shape_3rd_derivatives(const unsigned int interface_dof_index,
                                const unsigned int q_point,
                                const unsigned int component = 0) const;

  /**
   * @}
   */

  /**
   * @name Access to the average of shape function values and their derivatives
   * @{
   */

  /**
   * Return the average $\average{u}=\frac{1}{2}u_{\text{cell0}} +
   * \frac{1}{2}u_{\text{cell1}}$ on the interface
   * for the shape function @p interface_dof_index at the quadrature point
   * @p q_point of component @p component.
   *
   * If this is a boundary face (at_boundary() returns true), then
   * $\average{u}=u_{\text{cell0}}$.
   *
   * @note The name of the function is supposed to be read as "the average
   *   (singular) of the values (plural: one or two possible values)
   *   of the shape function (singular)".
   */
  double
  average_of_shape_values(const unsigned int interface_dof_index,
                          const unsigned int q_point,
                          const unsigned int component = 0) const;

  /**
   * Return the average of the gradient $\average{\nabla u} = \frac{1}{2}\nabla
   * u_{\text{cell0}} + \frac{1}{2} \nabla u_{\text{cell1}}$ on the interface
   * for the shape function @p interface_dof_index at the quadrature point @p
   * q_point of component @p component.
   *
   * If this is a boundary face (at_boundary() returns true), then
   * $\average{\nabla u}=\nabla u_{\text{cell0}}$.
   *
   * @note The name of the function is supposed to be read as "the average
   *   (singular) of the gradients (plural: one or two possible values
   *   for the gradient) of the shape function (singular)".
   */
  Tensor<1, spacedim>
  average_of_shape_gradients(const unsigned int interface_dof_index,
                             const unsigned int q_point,
                             const unsigned int component = 0) const;

  /**
   * Return the average of the Hessian $\average{\nabla^2 u} =
   * \frac{1}{2}\nabla^2 u_{\text{cell0}} + \frac{1}{2} \nabla^2
   * u_{\text{cell1}}$ on the interface
   * for the shape function @p interface_dof_index at the quadrature point @p
   * q_point of component @p component.
   *
   * If this is a boundary face (at_boundary() returns true), then
   * $\average{\nabla^2 u}=\nabla^2 u_{\text{cell0}}$.
   *
   * @note The name of the function is supposed to be read as "the average
   *   (singular) of the Hessians (plural: one or two possible values
   *   for the second derivatives) of the shape function (singular)".
   */
  Tensor<2, spacedim>
  average_of_shape_hessians(const unsigned int interface_dof_index,
                            const unsigned int q_point,
                            const unsigned int component = 0) const;

  /**
   * @}
   */



  /**
   * @name Access to jumps in the function values and derivatives
   * @{
   */

  /**
   * Return the jump in the values of the
   * finite element function characterized by <tt>fe_function</tt> at the
   * quadrature points of the cell interface selected the last time
   * the <tt>reinit</tt> function of the FEInterfaceValues object was called.
   *
   * @dealiiRequiresUpdateFlags{update_values}
   */
  template <class InputVector>
  void
  get_jump_in_function_values(
    const InputVector                             &fe_function,
    std::vector<typename InputVector::value_type> &values) const;

  /**
   * Return the jump in the gradients of the
   * finite element function characterized by <tt>fe_function</tt> at the
   * quadrature points of the cell interface selected the last time
   * the <tt>reinit</tt> function of the FEInterfaceValues object was called.
   *
   * @dealiiRequiresUpdateFlags{update_gradients}
   */
  template <class InputVector>
  void
  get_jump_in_function_gradients(
    const InputVector &fe_function,
    std::vector<Tensor<1, spacedim, typename InputVector::value_type>>
      &gradients) const;

  /**
   * Return the jump in the Hessians of the
   * finite element function characterized by <tt>fe_function</tt> at the
   * quadrature points of the cell interface selected the last time
   * the <tt>reinit</tt> function of the FEInterfaceValues object was called.
   * @dealiiRequiresUpdateFlags{update_hessians}
   */
  template <class InputVector>
  void
  get_jump_in_function_hessians(
    const InputVector &fe_function,
    std::vector<Tensor<2, spacedim, typename InputVector::value_type>>
      &hessians) const;

  /**
   * Return the jump in the third derivatives of the
   * finite element function characterized by <tt>fe_function</tt> at
   * the quadrature points of the cell interface selected the last time
   * the <tt>reinit</tt> function of the FEInterfaceValues object was called.
   *
   * @dealiiRequiresUpdateFlags{update_third_derivatives}
   */
  template <class InputVector>
  void
  get_jump_in_function_third_derivatives(
    const InputVector &fe_function,
    std::vector<Tensor<3, spacedim, typename InputVector::value_type>>
      &third_derivatives) const;

  /** @} */

  /**
   * @name Access to the average of the function values and derivatives
   */
  /** @{ */

  /**
   * Return the average of the values of the
   * finite element function characterized by <tt>fe_function</tt> at the
   * quadrature points of the cell interface selected the last time
   * the <tt>reinit</tt> function of the FEInterfaceValues object was called.
   *
   * @dealiiRequiresUpdateFlags{update_values}
   */
  template <class InputVector>
  void
  get_average_of_function_values(
    const InputVector                             &fe_function,
    std::vector<typename InputVector::value_type> &values) const;

  /**
   * Return the average of the gradients of the
   * finite element function characterized by <tt>fe_function</tt> at the
   * quadrature points of the cell interface selected the last time
   * the <tt>reinit</tt> function of the FEInterfaceValues object was called.
   * @dealiiRequiresUpdateFlags{update_gradients}
   */
  template <class InputVector>
  void
  get_average_of_function_gradients(
    const InputVector &fe_function,
    std::vector<Tensor<1, spacedim, typename InputVector::value_type>>
      &gradients) const;

  /**
   * Return the average of the Hessians of the
   * finite element function characterized by <tt>fe_function</tt> at the
   * quadrature points of the cell interface selected the last time
   * the <tt>reinit</tt> function of the FEInterfaceValues object was called.
   * @dealiiRequiresUpdateFlags{update_hessians}
   */
  template <class InputVector>
  void
  get_average_of_function_hessians(
    const InputVector &fe_function,
    std::vector<Tensor<2, spacedim, typename InputVector::value_type>>
      &hessians) const;

  /**
   * @}
   */



  /**
   * @name Extractors Methods to extract individual components
   * @{
   */

  /**
   * Create a view of the current FEInterfaceValues object that represents a
   * particular scalar component of the possibly vector-valued finite element.
   * The concept of views is explained in the documentation of the namespace
   * FEValuesViews.
   */
  FEInterfaceViews::Scalar<dim, spacedim>
  operator[](const FEValuesExtractors::Scalar &scalar) const;

  /**
   * Create a view of the current FEInterfaceValues object that represents a set
   * of <code>dim</code> scalar components (i.e. a vector) of the vector-valued
   * finite element. The concept of views is explained in the documentation of
   * the namespace FEValuesViews.
   */
  FEInterfaceViews::Vector<dim, spacedim>
  operator[](const FEValuesExtractors::Vector &vector) const;

  /**
   * @}
   */

private:
  /**
   * The list of DoF indices for the current interface, filled in reinit().
   */
  std::vector<types::global_dof_index> interface_dof_indices;

  /**
   * The mapping from interface dof to the two local dof indices of the
   * FEFaceValues objects. If an interface DoF is only active on one of the
   * cells, the other one will have numbers::invalid_unsigned_int.
   */
  std::vector<std::array<unsigned int, 2>> dofmap;

  /**
   * Pointer to internal_fe_face_values or internal_fe_subface_values,
   * respectively as determined in reinit().
   */
  FEFaceValuesBase<dim, spacedim> *fe_face_values;

  /**
   * Pointer to internal_fe_face_values_neighbor,
   * internal_fe_subface_values_neighbor, or nullptr, respectively
   * as determined in reinit().
   */
  FEFaceValuesBase<dim, spacedim> *fe_face_values_neighbor;

  /**
   * @name Data that supports the standard FE implementation
   */
  /** @{ */ // non-hp data

  /**
   * The FEFaceValues object for the current cell.
   */
  std::unique_ptr<FEFaceValues<dim, spacedim>> internal_fe_face_values;

  /**
   * The FEFaceValues object for the current cell if the cell is refined.
   */
  std::unique_ptr<FESubfaceValues<dim, spacedim>> internal_fe_subface_values;

  /**
   * The FEFaceValues object for the neighboring cell.
   */
  std::unique_ptr<FEFaceValues<dim, spacedim>> internal_fe_face_values_neighbor;

  /**
   * The FEFaceValues object for the neighboring cell if the cell is refined.
   */
  std::unique_ptr<FESubfaceValues<dim, spacedim>>
    internal_fe_subface_values_neighbor;

  /** @} */ // non-hp data

  /**
   * @name Data that supports the hp-FE implementation
   */
  /** @{ */ // hp data

  /**
   * An hp::FEValues object for the FEFaceValues on the
   * present cell.
   */
  std::unique_ptr<hp::FEFaceValues<dim, spacedim>> internal_hp_fe_face_values;

  /**
   * An hp::FEValues object for the FESubfaceValues on the
   * present cell.
   */
  std::unique_ptr<hp::FESubfaceValues<dim, spacedim>>
    internal_hp_fe_subface_values;

  /**
   * An hp::FEValues object for the FEFaceValues on the
   * neighbor of the present cell.
   */
  std::unique_ptr<hp::FEFaceValues<dim, spacedim>>
    internal_hp_fe_face_values_neighbor;

  /**
   * An hp::FEValues object for the FESubfaceValues on the
   * neighboring cell.
   */
  std::unique_ptr<hp::FESubfaceValues<dim, spacedim>>
    internal_hp_fe_subface_values_neighbor;

  /**
   * Exception used when a certain feature doesn't make sense when
   * FEInterfaceValues does has hp-capabilities enabled.
   *
   * @ingroup Exceptions
   */
  DeclExceptionMsg(ExcOnlyAvailableWithoutHP,
                   "The current function doesn't make sense when used with a "
                   "FEInterfaceValues object with hp-capabilities.");

  /**
   * Exception used when a certain feature doesn't make sense when
   * FEInterfaceValues does not have hp-capabilities enabled.
   *
   * @ingroup Exceptions
   */
  DeclExceptionMsg(ExcOnlyAvailableWithHP,
                   "The current function doesn't make sense when used with a "
                   "FEInterfaceValues object without hp-capabilities.");

  /** @} */ // hp data

  /*
   * Make the view classes friends of this class, since they
   * access internal data.
   */
  template <int, int>
  friend class FEInterfaceViews::Scalar;
  template <int, int>
  friend class FEInterfaceViews::Vector;
};



#ifndef DOXYGEN

/*---------------------- Inline functions ---------------------*/

template <int dim, int spacedim>
FEInterfaceValues<dim, spacedim>::FEInterfaceValues(
  const Mapping<dim, spacedim>       &mapping,
  const FiniteElement<dim, spacedim> &fe,
  const Quadrature<dim - 1>          &quadrature,
  const UpdateFlags                   update_flags)
  : n_quadrature_points(quadrature.size())
  , fe_face_values(nullptr)
  , fe_face_values_neighbor(nullptr)
  , internal_fe_face_values(
      std::make_unique<FEFaceValues<dim, spacedim>>(mapping,
                                                    fe,
                                                    quadrature,
                                                    update_flags))
  , internal_fe_subface_values(
      std::make_unique<FESubfaceValues<dim, spacedim>>(mapping,
                                                       fe,
                                                       quadrature,
                                                       update_flags))
  , internal_fe_face_values_neighbor(
      std::make_unique<FEFaceValues<dim, spacedim>>(mapping,
                                                    fe,
                                                    quadrature,
                                                    update_flags))
  , internal_fe_subface_values_neighbor(
      std::make_unique<FESubfaceValues<dim, spacedim>>(mapping,
                                                       fe,
                                                       quadrature,
                                                       update_flags))
{}



template <int dim, int spacedim>
FEInterfaceValues<dim, spacedim>::FEInterfaceValues(
  const FiniteElement<dim, spacedim> &fe,
  const Quadrature<dim - 1>          &quadrature,
  const UpdateFlags                   update_flags)
  : FEInterfaceValues(
      fe.reference_cell().template get_default_linear_mapping<dim, spacedim>(),
      fe,
      quadrature,
      update_flags)
{}



template <int dim, int spacedim>
FEInterfaceValues<dim, spacedim>::FEInterfaceValues(
  const Mapping<dim, spacedim>       &mapping,
  const FiniteElement<dim, spacedim> &fe,
  const hp::QCollection<dim - 1>     &quadrature,
  const UpdateFlags                   update_flags)
  : n_quadrature_points(quadrature.max_n_quadrature_points())
  , fe_face_values(nullptr)
  , fe_face_values_neighbor(nullptr)
  , internal_fe_face_values(
      std::make_unique<FEFaceValues<dim, spacedim>>(mapping,
                                                    fe,
                                                    quadrature,
                                                    update_flags))
  , internal_fe_subface_values(
      std::make_unique<FESubfaceValues<dim, spacedim>>(mapping,
                                                       fe,
                                                       quadrature,
                                                       update_flags))
  , internal_fe_face_values_neighbor(
      std::make_unique<FEFaceValues<dim, spacedim>>(mapping,
                                                    fe,
                                                    quadrature[0],
                                                    update_flags))
  , internal_fe_subface_values_neighbor(
      std::make_unique<FESubfaceValues<dim, spacedim>>(mapping,
                                                       fe,
                                                       quadrature[0],
                                                       update_flags))
{}



template <int dim, int spacedim>
FEInterfaceValues<dim, spacedim>::FEInterfaceValues(
  const hp::MappingCollection<dim, spacedim> &mapping_collection,
  const hp::FECollection<dim, spacedim>      &fe_collection,
  const hp::QCollection<dim - 1>             &quadrature_collection,
  const UpdateFlags                           update_flags)
  : n_quadrature_points(quadrature_collection.max_n_quadrature_points())
  , fe_face_values(nullptr)
  , fe_face_values_neighbor(nullptr)
  , internal_hp_fe_face_values(
      std::make_unique<hp::FEFaceValues<dim, spacedim>>(mapping_collection,
                                                        fe_collection,
                                                        quadrature_collection,
                                                        update_flags))
  , internal_hp_fe_subface_values(
      std::make_unique<hp::FESubfaceValues<dim, spacedim>>(
        mapping_collection,
        fe_collection,
        quadrature_collection,
        update_flags))
  , internal_hp_fe_face_values_neighbor(
      std::make_unique<hp::FEFaceValues<dim, spacedim>>(mapping_collection,
                                                        fe_collection,
                                                        quadrature_collection,
                                                        update_flags))
  , internal_hp_fe_subface_values_neighbor(
      std::make_unique<hp::FESubfaceValues<dim, spacedim>>(
        mapping_collection,
        fe_collection,
        quadrature_collection,
        update_flags))
{
  AssertDimension(dim, spacedim);
}



template <int dim, int spacedim>
FEInterfaceValues<dim, spacedim>::FEInterfaceValues(
  const hp::FECollection<dim, spacedim> &fe_collection,
  const hp::QCollection<dim - 1>        &quadrature_collection,
  const UpdateFlags                      update_flags)
  : FEInterfaceValues(fe_collection.get_reference_cell_default_linear_mapping(),
                      fe_collection,
                      quadrature_collection,
                      update_flags)
{}



template <int dim, int spacedim>
template <typename CellIteratorType, typename CellNeighborIteratorType>
void
FEInterfaceValues<dim, spacedim>::reinit(
  const CellIteratorType         &cell,
  const unsigned int              face_no,
  const unsigned int              sub_face_no,
  const CellNeighborIteratorType &cell_neighbor,
  const unsigned int              face_no_neighbor,
  const unsigned int              sub_face_no_neighbor,
  const unsigned int              q_index,
  const unsigned int              mapping_index,
  const unsigned int              fe_index_in,
  const unsigned int              fe_index_neighbor_in)
{
  Assert(internal_fe_face_values || internal_hp_fe_face_values,
         ExcNotInitialized());

  constexpr bool is_dof_cell_accessor =
    std::is_same_v<DoFCellAccessor<dim, spacedim, true>,
                   typename CellIteratorType::AccessorType> ||
    std::is_same_v<DoFCellAccessor<dim, spacedim, false>,
                   typename CellIteratorType::AccessorType>;

  constexpr bool is_dof_cell_accessor_neighbor =
    std::is_same_v<DoFCellAccessor<dim, spacedim, true>,
                   typename CellNeighborIteratorType::AccessorType> ||
    std::is_same_v<DoFCellAccessor<dim, spacedim, false>,
                   typename CellNeighborIteratorType::AccessorType>;

  unsigned int active_fe_index          = 0;
  unsigned int active_fe_index_neighbor = 0;

  if (internal_fe_face_values)
    {
      if (sub_face_no == numbers::invalid_unsigned_int)
        {
          internal_fe_face_values->reinit(cell, face_no);
          fe_face_values = internal_fe_face_values.get();
        }
      else
        {
          internal_fe_subface_values->reinit(cell, face_no, sub_face_no);
          fe_face_values = internal_fe_subface_values.get();
        }
      if (sub_face_no_neighbor == numbers::invalid_unsigned_int)
        {
          internal_fe_face_values_neighbor->reinit(cell_neighbor,
                                                   face_no_neighbor);
          fe_face_values_neighbor = internal_fe_face_values_neighbor.get();
        }
      else
        {
          internal_fe_subface_values_neighbor->reinit(cell_neighbor,
                                                      face_no_neighbor,
                                                      sub_face_no_neighbor);
          fe_face_values_neighbor = internal_fe_subface_values_neighbor.get();
        }

      AssertDimension(fe_face_values->n_quadrature_points,
                      fe_face_values_neighbor->n_quadrature_points);

      const_cast<unsigned int &>(this->n_quadrature_points) =
        fe_face_values->n_quadrature_points;
    }
  else if (internal_hp_fe_face_values)
    {
      active_fe_index = fe_index_in;
      active_fe_index_neighbor =
        (fe_index_neighbor_in != numbers::invalid_unsigned_int) ?
          fe_index_neighbor_in :
          fe_index_in;

      if (active_fe_index == numbers::invalid_unsigned_int)
        {
          if constexpr (is_dof_cell_accessor)
            active_fe_index = cell->active_fe_index();
          else
            active_fe_index = 0;
        }

      if (active_fe_index_neighbor == numbers::invalid_unsigned_int)
        {
          if constexpr (is_dof_cell_accessor_neighbor)
            active_fe_index_neighbor = cell_neighbor->active_fe_index();
          else
            active_fe_index_neighbor = 0;
        }

      unsigned int used_q_index       = q_index;
      unsigned int used_mapping_index = mapping_index;

      // First check. If there is only one element in a collection, and if none
      // had been specified explicitly, then that's clearly the one to take:
      if (used_q_index == numbers::invalid_unsigned_int)
        if (internal_hp_fe_face_values->get_quadrature_collection().size() == 1)
          used_q_index = 0;

      if (used_mapping_index == numbers::invalid_unsigned_int)
        if (internal_hp_fe_face_values->get_mapping_collection().size() == 1)
          used_mapping_index = 0;

      // Second check: See if the two quadrature objects are the same, because
      // in that case it does not matter which one we use. Unfortunately, we
      // currently have no way of testing that two mapping objects are the
      // same :-(
      if (used_q_index == numbers::invalid_unsigned_int)
        if (internal_hp_fe_face_values
              ->get_quadrature_collection()[active_fe_index] ==
            internal_hp_fe_face_values
              ->get_quadrature_collection()[active_fe_index_neighbor])
          used_q_index = active_fe_index;

      // Third check, if the above did not already suffice. We see if we
      // can get somewhere via the dominated's finite element index.
      if ((used_q_index == numbers::invalid_unsigned_int) ||
          (used_mapping_index == numbers::invalid_unsigned_int))
        {
          const unsigned int dominated_fe_index =
            ((used_q_index == numbers::invalid_unsigned_int) ||
                 (used_mapping_index == numbers::invalid_unsigned_int) ?
               internal_hp_fe_face_values->get_fe_collection()
                 .find_dominated_fe(
                   {active_fe_index, active_fe_index_neighbor}) :
               numbers::invalid_unsigned_int);

          if (used_q_index == numbers::invalid_unsigned_int)
            {
              Assert(dominated_fe_index != numbers::invalid_fe_index,
                     ExcMessage(
                       "You called this function with 'q_index' left at its "
                       "default value, but this can only work if one of "
                       "the two finite elements adjacent to this face "
                       "dominates the other. See the documentation "
                       "of this function for more information of how "
                       "to deal with this situation."));
              used_q_index = dominated_fe_index;
            }

          if (used_mapping_index == numbers::invalid_unsigned_int)
            {
              Assert(dominated_fe_index != numbers::invalid_fe_index,
                     ExcMessage(
                       "You called this function with 'mapping_index' left "
                       "at its default value, but this can only work if one "
                       "of the two finite elements adjacent to this face "
                       "dominates the other. See the documentation "
                       "of this function for more information of how "
                       "to deal with this situation."));
              used_mapping_index = dominated_fe_index;
            }
        }

      // Same as if above, but when hp is enabled.
      if (sub_face_no == numbers::invalid_unsigned_int)
        {
          internal_hp_fe_face_values->reinit(
            cell, face_no, used_q_index, used_mapping_index, active_fe_index);
          fe_face_values = &const_cast<FEFaceValues<dim, spacedim> &>(
            internal_hp_fe_face_values->get_present_fe_values());
        }
      else
        {
          internal_hp_fe_subface_values->reinit(cell,
                                                face_no,
                                                sub_face_no,
                                                used_q_index,
                                                used_mapping_index,
                                                active_fe_index);

          fe_face_values = &const_cast<FESubfaceValues<dim, spacedim> &>(
            internal_hp_fe_subface_values->get_present_fe_values());
        }
      if (sub_face_no_neighbor == numbers::invalid_unsigned_int)
        {
          internal_hp_fe_face_values_neighbor->reinit(cell_neighbor,
                                                      face_no_neighbor,
                                                      used_q_index,
                                                      used_mapping_index,
                                                      active_fe_index_neighbor);

          fe_face_values_neighbor = &const_cast<FEFaceValues<dim, spacedim> &>(
            internal_hp_fe_face_values_neighbor->get_present_fe_values());
        }
      else
        {
          internal_hp_fe_subface_values_neighbor->reinit(
            cell_neighbor,
            face_no_neighbor,
            sub_face_no_neighbor,
            used_q_index,
            used_mapping_index,
            active_fe_index_neighbor);

          fe_face_values_neighbor =
            &const_cast<FESubfaceValues<dim, spacedim> &>(
              internal_hp_fe_subface_values_neighbor->get_present_fe_values());
        }

      AssertDimension(fe_face_values->n_quadrature_points,
                      fe_face_values_neighbor->n_quadrature_points);

      const_cast<unsigned int &>(this->n_quadrature_points) =
        fe_face_values->n_quadrature_points;
    }

  // Set up dof mapping and remove duplicates (for continuous elements).
  if constexpr (is_dof_cell_accessor_neighbor && is_dof_cell_accessor)
    {
      // Get dof indices first:
      std::vector<types::global_dof_index> v(
        fe_face_values->get_fe().n_dofs_per_cell());
      cell->get_active_or_mg_dof_indices(v);
      std::vector<types::global_dof_index> v2(
        fe_face_values_neighbor->get_fe().n_dofs_per_cell());
      cell_neighbor->get_active_or_mg_dof_indices(v2);

      // Fill a map from the global dof index to the left and right
      // local index.
      std::map<types::global_dof_index, std::pair<unsigned int, unsigned int>>
                                            tempmap;
      std::pair<unsigned int, unsigned int> invalid_entry(
        numbers::invalid_unsigned_int, numbers::invalid_unsigned_int);

      for (unsigned int i = 0; i < v.size(); ++i)
        {
          // If not already existing, add an invalid entry:
          auto result = tempmap.insert(std::make_pair(v[i], invalid_entry));
          result.first->second.first = i;
        }

      for (unsigned int i = 0; i < v2.size(); ++i)
        {
          // If not already existing, add an invalid entry:
          auto result = tempmap.insert(std::make_pair(v2[i], invalid_entry));
          result.first->second.second = i;
        }

      // Transfer from the map to the sorted std::vectors.
      dofmap.resize(tempmap.size());
      interface_dof_indices.resize(tempmap.size());
      unsigned int idx = 0;
      for (auto &x : tempmap)
        {
          interface_dof_indices[idx] = x.first;
          dofmap[idx]                = {{x.second.first, x.second.second}};
          ++idx;
        }
    }
  else
    {
      const unsigned int n_dofs_per_cell_1 = fe_face_values->dofs_per_cell;
      const unsigned int n_dofs_per_cell_2 =
        fe_face_values_neighbor->dofs_per_cell;

      interface_dof_indices.resize(n_dofs_per_cell_1 + n_dofs_per_cell_2);
      dofmap.resize(n_dofs_per_cell_1 + n_dofs_per_cell_2);

      for (unsigned int i = 0; i < n_dofs_per_cell_1; ++i)
        {
          interface_dof_indices[i] = numbers::invalid_dof_index;
          dofmap[i]                = {{i, numbers::invalid_unsigned_int}};
        }

      for (unsigned int i = 0; i < n_dofs_per_cell_2; ++i)
        {
          interface_dof_indices[i + n_dofs_per_cell_1] =
            numbers::invalid_dof_index;
          dofmap[i + n_dofs_per_cell_1] = {{numbers::invalid_unsigned_int, i}};
        }
    }
}



template <int dim, int spacedim>
template <typename CellIteratorType>
void
FEInterfaceValues<dim, spacedim>::reinit(const CellIteratorType &cell,
                                         const unsigned int      face_no,
                                         const unsigned int      q_index,
                                         const unsigned int      mapping_index,
                                         const unsigned int      fe_index)
{
  Assert(internal_fe_face_values || internal_hp_fe_face_values,
         ExcNotInitialized());

  if (internal_fe_face_values)
    {
      Assert((q_index == 0 || q_index == numbers::invalid_unsigned_int),
             ExcNotImplemented());
      Assert((mapping_index == 0 ||
              mapping_index == numbers::invalid_unsigned_int),
             ExcNotImplemented());
      Assert((fe_index == 0 || fe_index == numbers::invalid_unsigned_int),
             ExcNotImplemented());

      internal_fe_face_values->reinit(cell, face_no);
      fe_face_values          = internal_fe_face_values.get();
      fe_face_values_neighbor = nullptr;
    }
  else if (internal_hp_fe_face_values)
    {
      internal_hp_fe_face_values->reinit(
        cell, face_no, q_index, mapping_index, fe_index);
      fe_face_values = &const_cast<FEFaceValues<dim> &>(
        internal_hp_fe_face_values->get_present_fe_values());
      fe_face_values_neighbor = nullptr;
    }

  interface_dof_indices.resize(fe_face_values->get_fe().n_dofs_per_cell());

  if constexpr (std::is_same_v<typename CellIteratorType::AccessorType,
                               DoFCellAccessor<dim, spacedim, true>> ||
                std::is_same_v<typename CellIteratorType::AccessorType,
                               DoFCellAccessor<dim, spacedim, false>>)
    {
      cell->get_active_or_mg_dof_indices(interface_dof_indices);
    }
  else
    {
      for (auto &i : interface_dof_indices)
        i = numbers::invalid_dof_index;
    }

  dofmap.resize(interface_dof_indices.size());
  for (unsigned int i = 0; i < interface_dof_indices.size(); ++i)
    dofmap[i] = {{i, numbers::invalid_unsigned_int}};
}



template <int dim, int spacedim>
inline double
FEInterfaceValues<dim, spacedim>::JxW(const unsigned int q) const
{
  Assert(fe_face_values != nullptr,
         ExcMessage("This call requires a call to reinit() first."));
  return fe_face_values->JxW(q);
}



template <int dim, int spacedim>
const std::vector<double> &
FEInterfaceValues<dim, spacedim>::get_JxW_values() const
{
  Assert(fe_face_values != nullptr,
         ExcMessage("This call requires a call to reinit() first."));
  return fe_face_values->get_JxW_values();
}



template <int dim, int spacedim>
const std::vector<Tensor<1, spacedim>> &
FEInterfaceValues<dim, spacedim>::get_normal_vectors() const
{
  Assert(fe_face_values != nullptr,
         ExcMessage("This call requires a call to reinit() first."));
  return fe_face_values->get_normal_vectors();
}



template <int dim, int spacedim>
const Mapping<dim, spacedim> &
FEInterfaceValues<dim, spacedim>::get_mapping() const
{
  Assert(!has_hp_capabilities(), ExcOnlyAvailableWithoutHP());
  return internal_fe_face_values->get_mapping();
}



template <int dim, int spacedim>
const FiniteElement<dim, spacedim> &
FEInterfaceValues<dim, spacedim>::get_fe() const
{
  Assert(!has_hp_capabilities(), ExcOnlyAvailableWithoutHP());
  return internal_fe_face_values->get_fe();
}



template <int dim, int spacedim>
const Quadrature<dim - 1> &
FEInterfaceValues<dim, spacedim>::get_quadrature() const
{
  Assert(!has_hp_capabilities(), ExcOnlyAvailableWithoutHP());
  return internal_fe_face_values->get_quadrature();
}



template <int dim, int spacedim>
const hp::MappingCollection<dim, spacedim> &
FEInterfaceValues<dim, spacedim>::get_mapping_collection() const
{
  Assert(has_hp_capabilities(), ExcOnlyAvailableWithHP());
  return internal_hp_fe_face_values->get_mapping_collection();
}



template <int dim, int spacedim>
const hp::FECollection<dim, spacedim> &
FEInterfaceValues<dim, spacedim>::get_fe_collection() const
{
  Assert(has_hp_capabilities(), ExcOnlyAvailableWithHP());
  return internal_hp_fe_face_values->get_fe_collection();
}



template <int dim, int spacedim>
const hp::QCollection<dim - 1> &
FEInterfaceValues<dim, spacedim>::get_quadrature_collection() const
{
  Assert(has_hp_capabilities(), ExcOnlyAvailableWithHP());
  return internal_hp_fe_face_values->get_quadrature_collection();
}



template <int dim, int spacedim>
bool
FEInterfaceValues<dim, spacedim>::has_hp_capabilities() const
{
  if (internal_hp_fe_face_values || internal_hp_fe_subface_values ||
      internal_hp_fe_face_values_neighbor ||
      internal_hp_fe_subface_values_neighbor)
    {
      Assert(!internal_fe_face_values, ExcInternalError());
      Assert(!internal_fe_subface_values, ExcInternalError());
      Assert(!internal_fe_face_values_neighbor, ExcInternalError());
      Assert(!internal_fe_subface_values_neighbor, ExcInternalError());

      return true;
    }

  Assert(internal_fe_face_values || internal_fe_subface_values ||
           internal_fe_face_values_neighbor ||
           internal_fe_subface_values_neighbor,
         ExcInternalError());
  Assert(!internal_hp_fe_face_values, ExcInternalError());
  Assert(!internal_hp_fe_subface_values, ExcInternalError());
  Assert(!internal_hp_fe_face_values_neighbor, ExcInternalError());
  Assert(!internal_hp_fe_subface_values_neighbor, ExcInternalError());

  return false;
}



template <int dim, int spacedim>
inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
FEInterfaceValues<dim, spacedim>::quadrature_point_indices() const
{
  return std_cxx20::ranges::iota_view<unsigned int, unsigned int>(
    0U, n_quadrature_points);
}



template <int dim, int spacedim>
const Point<spacedim> &
FEInterfaceValues<dim, spacedim>::quadrature_point(
  const unsigned int q_point) const
{
  Assert(fe_face_values != nullptr,
         ExcMessage("This call requires a call to reinit() first."));
  return fe_face_values->quadrature_point(q_point);
}



template <int dim, int spacedim>
const std::vector<Point<spacedim>> &
FEInterfaceValues<dim, spacedim>::get_quadrature_points() const
{
  Assert(fe_face_values != nullptr,
         ExcMessage("This call requires a call to reinit() first."));
  return fe_face_values->get_quadrature_points();
}



template <int dim, int spacedim>
UpdateFlags
FEInterfaceValues<dim, spacedim>::get_update_flags() const
{
  if (has_hp_capabilities())
    return internal_hp_fe_face_values->get_update_flags();
  else
    return internal_fe_face_values->get_update_flags();
}



template <int dim, int spacedim>
typename Triangulation<dim, spacedim>::cell_iterator
FEInterfaceValues<dim, spacedim>::get_cell(const unsigned int cell_index) const
{
  return get_fe_face_values(cell_index).get_cell();
}



template <int dim, int spacedim>
inline unsigned int
FEInterfaceValues<dim, spacedim>::get_face_number(
  const unsigned int cell_index) const
{
  return get_fe_face_values(cell_index).get_face_number();
}



template <int dim, int spacedim>
unsigned
FEInterfaceValues<dim, spacedim>::n_current_interface_dofs() const
{
  Assert(
    interface_dof_indices.size() > 0,
    ExcMessage(
      "n_current_interface_dofs() is only available after a call to reinit()."));
  return interface_dof_indices.size();
}



template <int dim, int spacedim>
inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
FEInterfaceValues<dim, spacedim>::dof_indices() const
{
  return {0U, n_current_interface_dofs()};
}



template <int dim, int spacedim>
bool
FEInterfaceValues<dim, spacedim>::at_boundary() const
{
  return fe_face_values_neighbor == nullptr;
}



template <int dim, int spacedim>
std::vector<types::global_dof_index>
FEInterfaceValues<dim, spacedim>::get_interface_dof_indices() const
{
  return interface_dof_indices;
}



template <int dim, int spacedim>
std::array<unsigned int, 2>
FEInterfaceValues<dim, spacedim>::interface_dof_to_dof_indices(
  const unsigned int interface_dof_index) const
{
  AssertIndexRange(interface_dof_index, n_current_interface_dofs());
  return dofmap[interface_dof_index];
}



template <int dim, int spacedim>
const FEFaceValuesBase<dim, spacedim> &
FEInterfaceValues<dim, spacedim>::get_fe_face_values(
  const unsigned int cell_index) const
{
  AssertIndexRange(cell_index, 2);
  Assert(
    cell_index == 0 || !at_boundary(),
    ExcMessage(
      "You are on a boundary, so you can only ask for the first FEFaceValues object."));

  return (cell_index == 0) ? *fe_face_values : *fe_face_values_neighbor;
}



template <int dim, int spacedim>
Tensor<1, spacedim>
FEInterfaceValues<dim, spacedim>::normal(const unsigned int q_point_index) const
{
  return normal_vector(q_point_index);
}



template <int dim, int spacedim>
Tensor<1, spacedim>
FEInterfaceValues<dim, spacedim>::normal_vector(
  const unsigned int q_point_index) const
{
  return fe_face_values->normal_vector(q_point_index);
}



template <int dim, int spacedim>
double
FEInterfaceValues<dim, spacedim>::shape_value(
  const bool         here_or_there,
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const auto dof_pair = dofmap[interface_dof_index];

  if (here_or_there && dof_pair[0] != numbers::invalid_unsigned_int)
    return get_fe_face_values(0).shape_value_component(dof_pair[0],
                                                       q_point,
                                                       component);
  if (!here_or_there && dof_pair[1] != numbers::invalid_unsigned_int)
    return get_fe_face_values(1).shape_value_component(dof_pair[1],
                                                       q_point,
                                                       component);

  return 0.0;
}


template <int dim, int spacedim>
Tensor<1, spacedim>
FEInterfaceValues<dim, spacedim>::shape_grad(
  const bool         here_or_there,
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const auto dof_pair = dofmap[interface_dof_index];

  Tensor<1, dim> value;

  if (here_or_there && dof_pair[0] != numbers::invalid_unsigned_int)
    value = get_fe_face_values(0).shape_grad_component(dof_pair[0],
                                                       q_point,
                                                       component);
  if (!here_or_there && dof_pair[1] != numbers::invalid_unsigned_int)
    value = get_fe_face_values(1).shape_grad_component(dof_pair[1],
                                                       q_point,
                                                       component);

  return value;
}

template <int dim, int spacedim>
double
FEInterfaceValues<dim, spacedim>::jump_in_shape_values(
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const auto dof_pair = dofmap[interface_dof_index];

  double value = 0.0;

  if (dof_pair[0] != numbers::invalid_unsigned_int)
    value += get_fe_face_values(0).shape_value_component(dof_pair[0],
                                                         q_point,
                                                         component);
  if (dof_pair[1] != numbers::invalid_unsigned_int)
    value -= get_fe_face_values(1).shape_value_component(dof_pair[1],
                                                         q_point,
                                                         component);
  return value;
}



template <int dim, int spacedim>
double
FEInterfaceValues<dim, spacedim>::average_of_shape_values(
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const auto dof_pair = dofmap[interface_dof_index];

  if (at_boundary())
    return get_fe_face_values(0).shape_value_component(dof_pair[0],
                                                       q_point,
                                                       component);

  double value = 0.0;

  if (dof_pair[0] != numbers::invalid_unsigned_int)
    value += 0.5 * get_fe_face_values(0).shape_value_component(dof_pair[0],
                                                               q_point,
                                                               component);
  if (dof_pair[1] != numbers::invalid_unsigned_int)
    value += 0.5 * get_fe_face_values(1).shape_value_component(dof_pair[1],
                                                               q_point,
                                                               component);

  return value;
}



template <int dim, int spacedim>
Tensor<1, spacedim>
FEInterfaceValues<dim, spacedim>::average_of_shape_gradients(
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const auto dof_pair = dofmap[interface_dof_index];

  if (at_boundary())
    return get_fe_face_values(0).shape_grad_component(dof_pair[0],
                                                      q_point,
                                                      component);

  Tensor<1, spacedim> value;

  if (dof_pair[0] != numbers::invalid_unsigned_int)
    value += 0.5 * get_fe_face_values(0).shape_grad_component(dof_pair[0],
                                                              q_point,
                                                              component);
  if (dof_pair[1] != numbers::invalid_unsigned_int)
    value += 0.5 * get_fe_face_values(1).shape_grad_component(dof_pair[1],
                                                              q_point,
                                                              component);

  return value;
}



template <int dim, int spacedim>
Tensor<2, spacedim>
FEInterfaceValues<dim, spacedim>::average_of_shape_hessians(
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const auto dof_pair = dofmap[interface_dof_index];

  if (at_boundary())
    return get_fe_face_values(0).shape_hessian_component(dof_pair[0],
                                                         q_point,
                                                         component);

  Tensor<2, spacedim> value;

  if (dof_pair[0] != numbers::invalid_unsigned_int)
    value += 0.5 * get_fe_face_values(0).shape_hessian_component(dof_pair[0],
                                                                 q_point,
                                                                 component);
  if (dof_pair[1] != numbers::invalid_unsigned_int)
    value += 0.5 * get_fe_face_values(1).shape_hessian_component(dof_pair[1],
                                                                 q_point,
                                                                 component);

  return value;
}



template <int dim, int spacedim>
Tensor<1, spacedim>
FEInterfaceValues<dim, spacedim>::jump_in_shape_gradients(
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const auto dof_pair = dofmap[interface_dof_index];

  if (at_boundary())
    return get_fe_face_values(0).shape_grad_component(dof_pair[0],
                                                      q_point,
                                                      component);

  Tensor<1, spacedim> value;

  if (dof_pair[0] != numbers::invalid_unsigned_int)
    value += get_fe_face_values(0).shape_grad_component(dof_pair[0],
                                                        q_point,
                                                        component);
  if (dof_pair[1] != numbers::invalid_unsigned_int)
    value -= get_fe_face_values(1).shape_grad_component(dof_pair[1],
                                                        q_point,
                                                        component);

  return value;
}



template <int dim, int spacedim>
Tensor<2, spacedim>
FEInterfaceValues<dim, spacedim>::jump_in_shape_hessians(
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const auto dof_pair = dofmap[interface_dof_index];

  if (at_boundary())
    return get_fe_face_values(0).shape_hessian_component(dof_pair[0],
                                                         q_point,
                                                         component);

  Tensor<2, spacedim> value;

  if (dof_pair[0] != numbers::invalid_unsigned_int)
    value += get_fe_face_values(0).shape_hessian_component(dof_pair[0],
                                                           q_point,
                                                           component);
  if (dof_pair[1] != numbers::invalid_unsigned_int)
    value -= get_fe_face_values(1).shape_hessian_component(dof_pair[1],
                                                           q_point,
                                                           component);

  return value;
}



template <int dim, int spacedim>
Tensor<3, spacedim>
FEInterfaceValues<dim, spacedim>::jump_in_shape_3rd_derivatives(
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const auto dof_pair = dofmap[interface_dof_index];

  if (at_boundary())
    return get_fe_face_values(0).shape_3rd_derivative_component(dof_pair[0],
                                                                q_point,
                                                                component);

  Tensor<3, spacedim> value;

  if (dof_pair[0] != numbers::invalid_unsigned_int)
    value += get_fe_face_values(0).shape_3rd_derivative_component(dof_pair[0],
                                                                  q_point,
                                                                  component);
  if (dof_pair[1] != numbers::invalid_unsigned_int)
    value -= get_fe_face_values(1).shape_3rd_derivative_component(dof_pair[1],
                                                                  q_point,
                                                                  component);

  return value;
}



template <int dim, int spacedim>
template <class InputVector>
void
FEInterfaceValues<dim, spacedim>::get_jump_in_function_values(
  const InputVector                             &fe_function,
  std::vector<typename InputVector::value_type> &values) const
{
  AssertDimension(values.size(), n_quadrature_points);

  const FEValuesExtractors::Scalar scalar(0);
  this->operator[](scalar).get_jump_in_function_values(fe_function, values);
}



template <int dim, int spacedim>
template <class InputVector>
void
FEInterfaceValues<dim, spacedim>::get_jump_in_function_gradients(
  const InputVector &fe_function,
  std::vector<Tensor<1, spacedim, typename InputVector::value_type>> &gradients)
  const
{
  AssertDimension(gradients.size(), n_quadrature_points);

  const FEValuesExtractors::Scalar scalar(0);
  this->operator[](scalar).get_jump_in_function_gradients(fe_function,
                                                          gradients);
}



template <int dim, int spacedim>
template <class InputVector>
void
FEInterfaceValues<dim, spacedim>::get_jump_in_function_hessians(
  const InputVector &fe_function,
  std::vector<Tensor<2, spacedim, typename InputVector::value_type>> &hessians)
  const
{
  AssertDimension(hessians.size(), n_quadrature_points);

  const FEValuesExtractors::Scalar scalar(0);
  this->operator[](scalar).get_jump_in_function_hessians(fe_function, hessians);
}



template <int dim, int spacedim>
template <class InputVector>
void
FEInterfaceValues<dim, spacedim>::get_jump_in_function_third_derivatives(
  const InputVector &fe_function,
  std::vector<Tensor<3, spacedim, typename InputVector::value_type>>
    &third_derivatives) const
{
  AssertDimension(third_derivatives.size(), n_quadrature_points);

  const FEValuesExtractors::Scalar scalar(0);
  this->operator[](scalar).get_jump_in_function_third_derivatives(
    fe_function, third_derivatives);
}



template <int dim, int spacedim>
template <class InputVector>
void
FEInterfaceValues<dim, spacedim>::get_average_of_function_values(
  const InputVector                             &fe_function,
  std::vector<typename InputVector::value_type> &values) const
{
  AssertDimension(values.size(), n_quadrature_points);

  const FEValuesExtractors::Scalar scalar(0);
  this->operator[](scalar).get_average_of_function_values(fe_function, values);
}



template <int dim, int spacedim>
template <class InputVector>
void
FEInterfaceValues<dim, spacedim>::get_average_of_function_gradients(
  const InputVector &fe_function,
  std::vector<Tensor<1, spacedim, typename InputVector::value_type>> &gradients)
  const
{
  AssertDimension(gradients.size(), n_quadrature_points);

  const FEValuesExtractors::Scalar scalar(0);
  this->operator[](scalar).get_average_of_function_gradients(fe_function,
                                                             gradients);
}



template <int dim, int spacedim>
template <class InputVector>
void
FEInterfaceValues<dim, spacedim>::get_average_of_function_hessians(
  const InputVector &fe_function,
  std::vector<Tensor<2, spacedim, typename InputVector::value_type>> &hessians)
  const
{
  AssertDimension(hessians.size(), n_quadrature_points);

  const FEValuesExtractors::Scalar scalar(0);
  this->operator[](scalar).get_average_of_function_hessians(fe_function,
                                                            hessians);
}



/*------------ Inline functions: FEInterfaceValues------------*/
template <int dim, int spacedim>
inline FEInterfaceViews::Scalar<dim, spacedim>
FEInterfaceValues<dim, spacedim>::operator[](
  const FEValuesExtractors::Scalar &scalar) const
{
  const unsigned int n_components =
    (this->has_hp_capabilities() ? this->get_fe_collection().n_components() :
                                   this->get_fe().n_components());
  (void)n_components;
  AssertIndexRange(scalar.component, n_components);
  return FEInterfaceViews::Scalar<dim, spacedim>(*this, scalar.component);
}



template <int dim, int spacedim>
inline FEInterfaceViews::Vector<dim, spacedim>
FEInterfaceValues<dim, spacedim>::operator[](
  const FEValuesExtractors::Vector &vector) const
{
  const unsigned int n_components =
    (this->has_hp_capabilities() ? this->get_fe_collection().n_components() :
                                   this->get_fe().n_components());
  const unsigned int n_vectors =
    (n_components >= Tensor<1, spacedim>::n_independent_components ?
       n_components - Tensor<1, spacedim>::n_independent_components + 1 :
       0);
  (void)n_components;
  (void)n_vectors;
  AssertIndexRange(vector.first_vector_component, n_vectors);
  return FEInterfaceViews::Vector<dim, spacedim>(*this,
                                                 vector.first_vector_component);
}



namespace FEInterfaceViews
{
  template <int dim, int spacedim>
  Base<dim, spacedim>::Base(
    const FEInterfaceValues<dim, spacedim> &fe_interface)
    : fe_interface(&fe_interface)
  {}



  template <int dim, int spacedim>
  Scalar<dim, spacedim>::Scalar(
    const FEInterfaceValues<dim, spacedim> &fe_interface,
    const unsigned int                      component)
    : Base<dim, spacedim>(fe_interface)
    , extractor(component)
  {}



  template <int dim, int spacedim>
  template <class InputVector, class OutputVector>
  void
  Base<dim, spacedim>::get_local_dof_values(
    const InputVector &dof_values,
    OutputVector      &local_dof_values) const
  {
    const auto &interface_dof_indices =
      this->fe_interface->get_interface_dof_indices();

    AssertDimension(interface_dof_indices.size(), local_dof_values.size());

    for (const unsigned int i : this->fe_interface->dof_indices())
      local_dof_values[i] = dof_values(interface_dof_indices[i]);
  }



  template <int dim, int spacedim>
  typename Scalar<dim, spacedim>::value_type
  Scalar<dim, spacedim>::value(const bool         here_or_there,
                               const unsigned int interface_dof_index,
                               const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (here_or_there && dof_pair[0] != numbers::invalid_unsigned_int)
      return (*(this->fe_interface->fe_face_values))[extractor].value(
        dof_pair[0], q_point);

    if (!here_or_there && dof_pair[1] != numbers::invalid_unsigned_int)
      return (*(this->fe_interface->fe_face_values_neighbor))[extractor].value(
        dof_pair[1], q_point);

    return 0.0;
  }



  template <int dim, int spacedim>
  typename Scalar<dim, spacedim>::value_type
  Scalar<dim, spacedim>::jump_in_values(const unsigned int interface_dof_index,
                                        const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    value_type value = 0.0;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        (*(this->fe_interface->fe_face_values))[extractor].value(dof_pair[0],
                                                                 q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value -=
        (*(this->fe_interface->fe_face_values_neighbor))[extractor].value(
          dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  typename Scalar<dim, spacedim>::value_type
  Scalar<dim, spacedim>::average_of_values(
    const unsigned int interface_dof_index,
    const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (this->fe_interface->at_boundary())
      return (*(this->fe_interface->fe_face_values))[extractor].value(
        dof_pair[0], q_point);

    value_type value = 0.0;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        0.5 *
        (*(this->fe_interface->fe_face_values))[extractor].value(dof_pair[0],
                                                                 q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value +=
        0.5 * (*(this->fe_interface->fe_face_values_neighbor))[extractor].value(
                dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  typename Scalar<dim, spacedim>::gradient_type
  Scalar<dim, spacedim>::average_of_gradients(
    const unsigned int interface_dof_index,
    const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (this->fe_interface->at_boundary())
      return (*(this->fe_interface->fe_face_values))[extractor].gradient(
        dof_pair[0], q_point);

    gradient_type value;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        0.5 *
        (*(this->fe_interface->fe_face_values))[extractor].gradient(dof_pair[0],
                                                                    q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value += 0.5 * (*(this->fe_interface->fe_face_values_neighbor))[extractor]
                       .gradient(dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  typename Scalar<dim, spacedim>::gradient_type
  Scalar<dim, spacedim>::jump_in_gradients(
    const unsigned int interface_dof_index,
    const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (this->fe_interface->at_boundary())
      return (*(this->fe_interface->fe_face_values))[extractor].gradient(
        dof_pair[0], q_point);

    gradient_type value;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        (*(this->fe_interface->fe_face_values))[extractor].gradient(dof_pair[0],
                                                                    q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value -=
        (*(this->fe_interface->fe_face_values_neighbor))[extractor].gradient(
          dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  typename Scalar<dim, spacedim>::hessian_type
  Scalar<dim, spacedim>::average_of_hessians(
    const unsigned int interface_dof_index,
    const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (this->fe_interface->at_boundary())
      return (*(this->fe_interface->fe_face_values))[extractor].hessian(
        dof_pair[0], q_point);

    hessian_type value;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        0.5 *
        (*(this->fe_interface->fe_face_values))[extractor].hessian(dof_pair[0],
                                                                   q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value += 0.5 * (*(this->fe_interface->fe_face_values_neighbor))[extractor]
                       .hessian(dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  typename Scalar<dim, spacedim>::third_derivative_type
  Scalar<dim, spacedim>::jump_in_third_derivatives(
    const unsigned int interface_dof_index,
    const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (this->fe_interface->at_boundary())
      return (*(this->fe_interface->fe_face_values))[extractor]
        .third_derivative(dof_pair[0], q_point);

    third_derivative_type value;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        (*(this->fe_interface->fe_face_values))[extractor].third_derivative(
          dof_pair[0], q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value -= (*(this->fe_interface->fe_face_values_neighbor))[extractor]
                 .third_derivative(dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  typename Scalar<dim, spacedim>::hessian_type
  Scalar<dim, spacedim>::jump_in_hessians(
    const unsigned int interface_dof_index,
    const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (this->fe_interface->at_boundary())
      return (*(this->fe_interface->fe_face_values))[extractor].hessian(
        dof_pair[0], q_point);

    hessian_type value;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        (*(this->fe_interface->fe_face_values))[extractor].hessian(dof_pair[0],
                                                                   q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value -=
        (*(this->fe_interface->fe_face_values_neighbor))[extractor].hessian(
          dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Scalar<dim, spacedim>::get_function_values_from_local_dof_values(
    const bool         here_or_there,
    const InputVector &local_dof_values,
    std::vector<solution_value_type<typename InputVector::value_type>> &values)
    const
  {
    AssertDimension(values.size(), this->fe_interface->n_quadrature_points);

    for (const auto dof_index : this->fe_interface->dof_indices())
      for (const auto q_index : this->fe_interface->quadrature_point_indices())
        {
          if (dof_index == 0)
            values[q_index] = 0.;

          values[q_index] += local_dof_values[dof_index] *
                             value(here_or_there, dof_index, q_index);
        }
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Scalar<dim, spacedim>::get_function_values(
    const bool         here_or_there,
    const InputVector &fe_function,
    std::vector<solution_value_type<typename InputVector::value_type>> &values)
    const
  {
    std::vector<typename InputVector::value_type> local_dof_values(
      this->fe_interface->n_current_interface_dofs());
    this->get_local_dof_values(fe_function, local_dof_values);

    get_function_values_from_local_dof_values(here_or_there,
                                              local_dof_values,
                                              values);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Scalar<dim, spacedim>::get_jump_in_function_values_from_local_dof_values(
    const InputVector &local_dof_values,
    std::vector<solution_value_type<typename InputVector::value_type>> &values)
    const
  {
    AssertDimension(values.size(), this->fe_interface->n_quadrature_points);

    for (const auto dof_index : this->fe_interface->dof_indices())
      for (const auto q_index : this->fe_interface->quadrature_point_indices())
        {
          if (dof_index == 0)
            values[q_index] = 0.;

          values[q_index] +=
            local_dof_values[dof_index] * jump_in_values(dof_index, q_index);
        }
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Scalar<dim, spacedim>::get_jump_in_function_values(
    const InputVector &fe_function,
    std::vector<solution_value_type<typename InputVector::value_type>> &values)
    const
  {
    std::vector<typename InputVector::value_type> local_dof_values(
      this->fe_interface->n_current_interface_dofs());
    this->get_local_dof_values(fe_function, local_dof_values);

    get_jump_in_function_values_from_local_dof_values(local_dof_values, values);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Scalar<dim, spacedim>::get_jump_in_function_gradients_from_local_dof_values(
    const InputVector &local_dof_values,
    std::vector<solution_gradient_type<typename InputVector::value_type>>
      &gradients) const
  {
    AssertDimension(gradients.size(), this->fe_interface->n_quadrature_points);

    for (const auto dof_index : this->fe_interface->dof_indices())
      for (const auto q_index : this->fe_interface->quadrature_point_indices())
        {
          if (dof_index == 0)
            gradients[q_index] = 0.;

          gradients[q_index] +=
            local_dof_values[dof_index] * jump_in_gradients(dof_index, q_index);
        }
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Scalar<dim, spacedim>::get_jump_in_function_gradients(
    const InputVector &fe_function,
    std::vector<solution_gradient_type<typename InputVector::value_type>>
      &gradients) const
  {
    std::vector<typename InputVector::value_type> local_dof_values(
      this->fe_interface->n_current_interface_dofs());
    this->get_local_dof_values(fe_function, local_dof_values);

    get_jump_in_function_gradients_from_local_dof_values(local_dof_values,
                                                         gradients);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Scalar<dim, spacedim>::get_average_of_function_values_from_local_dof_values(
    const InputVector &local_dof_values,
    std::vector<solution_value_type<typename InputVector::value_type>> &values)
    const
  {
    AssertDimension(values.size(), this->fe_interface->n_quadrature_points);

    for (const auto dof_index : this->fe_interface->dof_indices())
      for (const auto q_index : this->fe_interface->quadrature_point_indices())
        {
          if (dof_index == 0)
            values[q_index] = 0.;

          values[q_index] +=
            local_dof_values[dof_index] * average_of_values(dof_index, q_index);
        }
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Scalar<dim, spacedim>::get_average_of_function_values(
    const InputVector &fe_function,
    std::vector<solution_value_type<typename InputVector::value_type>> &values)
    const
  {
    std::vector<typename InputVector::value_type> local_dof_values(
      this->fe_interface->n_current_interface_dofs());
    this->get_local_dof_values(fe_function, local_dof_values);

    get_average_of_function_values_from_local_dof_values(local_dof_values,
                                                         values);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Scalar<dim, spacedim>::
    get_average_of_function_gradients_from_local_dof_values(
      const InputVector &local_dof_values,
      std::vector<solution_gradient_type<typename InputVector::value_type>>
        &gradients) const
  {
    AssertDimension(gradients.size(), this->fe_interface->n_quadrature_points);

    for (const auto dof_index : this->fe_interface->dof_indices())
      for (const auto q_index : this->fe_interface->quadrature_point_indices())
        {
          if (dof_index == 0)
            gradients[q_index] = 0.;

          gradients[q_index] += local_dof_values[dof_index] *
                                average_of_gradients(dof_index, q_index);
        }
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Scalar<dim, spacedim>::get_average_of_function_gradients(
    const InputVector &fe_function,
    std::vector<solution_gradient_type<typename InputVector::value_type>>
      &gradients) const
  {
    std::vector<typename InputVector::value_type> local_dof_values(
      this->fe_interface->n_current_interface_dofs());
    this->get_local_dof_values(fe_function, local_dof_values);

    get_average_of_function_gradients_from_local_dof_values(local_dof_values,
                                                            gradients);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Scalar<dim, spacedim>::get_jump_in_function_hessians_from_local_dof_values(
    const InputVector &local_dof_values,
    std::vector<solution_hessian_type<typename InputVector::value_type>>
      &hessians) const
  {
    AssertDimension(hessians.size(), this->fe_interface->n_quadrature_points);

    for (const auto dof_index : this->fe_interface->dof_indices())
      for (const auto q_index : this->fe_interface->quadrature_point_indices())
        {
          if (dof_index == 0)
            hessians[q_index] = 0.;

          hessians[q_index] +=
            local_dof_values[dof_index] * jump_in_hessians(dof_index, q_index);
        }
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Scalar<dim, spacedim>::get_jump_in_function_hessians(
    const InputVector &fe_function,
    std::vector<solution_hessian_type<typename InputVector::value_type>>
      &hessians) const
  {
    std::vector<typename InputVector::value_type> local_dof_values(
      this->fe_interface->n_current_interface_dofs());
    this->get_local_dof_values(fe_function, local_dof_values);

    get_jump_in_function_hessians_from_local_dof_values(local_dof_values,
                                                        hessians);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Scalar<dim, spacedim>::get_average_of_function_hessians_from_local_dof_values(
    const InputVector &local_dof_values,
    std::vector<solution_hessian_type<typename InputVector::value_type>>
      &hessians) const
  {
    AssertDimension(hessians.size(), this->fe_interface->n_quadrature_points);

    for (const unsigned int dof_index : this->fe_interface->dof_indices())
      for (const auto q_index : this->fe_interface->quadrature_point_indices())
        {
          if (dof_index == 0)
            hessians[q_index] = 0.;

          hessians[q_index] += local_dof_values[dof_index] *
                               average_of_hessians(dof_index, q_index);
        }
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Scalar<dim, spacedim>::get_average_of_function_hessians(
    const InputVector &fe_function,
    std::vector<solution_hessian_type<typename InputVector::value_type>>
      &hessians) const
  {
    std::vector<typename InputVector::value_type> local_dof_values(
      this->fe_interface->n_current_interface_dofs());
    this->get_local_dof_values(fe_function, local_dof_values);

    get_average_of_function_hessians_from_local_dof_values(local_dof_values,
                                                           hessians);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Scalar<dim, spacedim>::
    get_jump_in_function_third_derivatives_from_local_dof_values(
      const InputVector &local_dof_values,
      std::vector<
        solution_third_derivative_type<typename InputVector::value_type>>
        &third_derivatives) const
  {
    AssertDimension(third_derivatives.size(),
                    this->fe_interface->n_quadrature_points);

    for (const unsigned int dof_index : this->fe_interface->dof_indices())
      for (const auto q_index : this->fe_interface->quadrature_point_indices())
        {
          if (dof_index == 0)
            third_derivatives[q_index] = 0.;

          third_derivatives[q_index] +=
            local_dof_values[dof_index] *
            jump_in_third_derivatives(dof_index, q_index);
        }
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Scalar<dim, spacedim>::get_jump_in_function_third_derivatives(
    const InputVector &fe_function,
    std::vector<
      solution_third_derivative_type<typename InputVector::value_type>>
      &third_derivatives) const
  {
    std::vector<typename InputVector::value_type> local_dof_values(
      this->fe_interface->n_current_interface_dofs());
    this->get_local_dof_values(fe_function, local_dof_values);

    get_jump_in_function_third_derivatives_from_local_dof_values(
      local_dof_values, third_derivatives);
  }



  template <int dim, int spacedim>
  Vector<dim, spacedim>::Vector(
    const FEInterfaceValues<dim, spacedim> &fe_interface,
    const unsigned int                      first_vector_component)
    : Base<dim, spacedim>(fe_interface)
    , extractor(first_vector_component)
  {}



  template <int dim, int spacedim>
  typename Vector<dim, spacedim>::value_type
  Vector<dim, spacedim>::value(const bool         here_or_there,
                               const unsigned int interface_dof_index,
                               const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (here_or_there && dof_pair[0] != numbers::invalid_unsigned_int)
      return (*(this->fe_interface->fe_face_values))[extractor].value(
        dof_pair[0], q_point);

    if (!here_or_there && dof_pair[1] != numbers::invalid_unsigned_int)
      return (*(this->fe_interface->fe_face_values_neighbor))[extractor].value(
        dof_pair[1], q_point);

    return value_type();
  }



  template <int dim, int spacedim>
  typename Vector<dim, spacedim>::value_type
  Vector<dim, spacedim>::jump_in_values(const unsigned int interface_dof_index,
                                        const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    value_type value;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        (*(this->fe_interface->fe_face_values))[extractor].value(dof_pair[0],
                                                                 q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value -=
        (*(this->fe_interface->fe_face_values_neighbor))[extractor].value(
          dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  typename Vector<dim, spacedim>::value_type
  Vector<dim, spacedim>::average_of_values(
    const unsigned int interface_dof_index,
    const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (this->fe_interface->at_boundary())
      return (*(this->fe_interface->fe_face_values))[extractor].value(
        dof_pair[0], q_point);

    value_type value;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        0.5 *
        (*(this->fe_interface->fe_face_values))[extractor].value(dof_pair[0],
                                                                 q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value +=
        0.5 * (*(this->fe_interface->fe_face_values_neighbor))[extractor].value(
                dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  typename Vector<dim, spacedim>::gradient_type
  Vector<dim, spacedim>::average_of_gradients(
    const unsigned int interface_dof_index,
    const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (this->fe_interface->at_boundary())
      return (*(this->fe_interface->fe_face_values))[extractor].gradient(
        dof_pair[0], q_point);

    gradient_type value;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        0.5 *
        (*(this->fe_interface->fe_face_values))[extractor].gradient(dof_pair[0],
                                                                    q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value += 0.5 * (*(this->fe_interface->fe_face_values_neighbor))[extractor]
                       .gradient(dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  typename Vector<dim, spacedim>::gradient_type
  Vector<dim, spacedim>::jump_in_gradients(
    const unsigned int interface_dof_index,
    const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (this->fe_interface->at_boundary())
      return (*(this->fe_interface->fe_face_values))[extractor].gradient(
        dof_pair[0], q_point);

    gradient_type value;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        (*(this->fe_interface->fe_face_values))[extractor].gradient(dof_pair[0],
                                                                    q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value -=
        (*(this->fe_interface->fe_face_values_neighbor))[extractor].gradient(
          dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  typename Vector<dim, spacedim>::gradient_type
  Vector<dim, spacedim>::jump_gradient(const unsigned int interface_dof_index,
                                       const unsigned int q_point) const
  {
    return jump_in_gradients(interface_dof_index, q_point);
  }



  template <int dim, int spacedim>
  typename Vector<dim, spacedim>::hessian_type
  Vector<dim, spacedim>::average_of_hessians(
    const unsigned int interface_dof_index,
    const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (this->fe_interface->at_boundary())
      return (*(this->fe_interface->fe_face_values))[extractor].hessian(
        dof_pair[0], q_point);

    hessian_type value;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        0.5 *
        (*(this->fe_interface->fe_face_values))[extractor].hessian(dof_pair[0],
                                                                   q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value += 0.5 * (*(this->fe_interface->fe_face_values_neighbor))[extractor]
                       .hessian(dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  typename Vector<dim, spacedim>::hessian_type
  Vector<dim, spacedim>::average_hessian(const unsigned int interface_dof_index,
                                         const unsigned int q_point) const
  {
    return average_of_hessians(interface_dof_index, q_point);
  }



  template <int dim, int spacedim>
  typename Vector<dim, spacedim>::hessian_type
  Vector<dim, spacedim>::jump_in_hessians(
    const unsigned int interface_dof_index,
    const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (this->fe_interface->at_boundary())
      return (*(this->fe_interface->fe_face_values))[extractor].hessian(
        dof_pair[0], q_point);

    hessian_type value;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        (*(this->fe_interface->fe_face_values))[extractor].hessian(dof_pair[0],
                                                                   q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value -=
        (*(this->fe_interface->fe_face_values_neighbor))[extractor].hessian(
          dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  typename Vector<dim, spacedim>::hessian_type
  Vector<dim, spacedim>::jump_hessian(const unsigned int interface_dof_index,
                                      const unsigned int q_point) const
  {
    return jump_in_hessians(interface_dof_index, q_point);
  }



  template <int dim, int spacedim>
  typename Vector<dim, spacedim>::third_derivative_type
  Vector<dim, spacedim>::jump_in_third_derivatives(
    const unsigned int interface_dof_index,
    const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (this->fe_interface->at_boundary())
      return (*(this->fe_interface->fe_face_values))[extractor]
        .third_derivative(dof_pair[0], q_point);

    third_derivative_type value;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        (*(this->fe_interface->fe_face_values))[extractor].third_derivative(
          dof_pair[0], q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value -= (*(this->fe_interface->fe_face_values_neighbor))[extractor]
                 .third_derivative(dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Vector<dim, spacedim>::get_function_values_from_local_dof_values(
    const bool         here_or_there,
    const InputVector &local_dof_values,
    std::vector<solution_value_type<typename InputVector::value_type>> &values)
    const
  {
    AssertDimension(values.size(), this->fe_interface->n_quadrature_points);

    for (const auto dof_index : this->fe_interface->dof_indices())
      for (const auto q_index : this->fe_interface->quadrature_point_indices())
        {
          if (dof_index == 0)
            values[q_index] = 0.;

          values[q_index] += local_dof_values[dof_index] *
                             value(here_or_there, dof_index, q_index);
        }
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Vector<dim, spacedim>::get_function_values(
    const bool         here_or_there,
    const InputVector &fe_function,
    std::vector<solution_value_type<typename InputVector::value_type>> &values)
    const
  {
    std::vector<typename InputVector::value_type> local_dof_values(
      this->fe_interface->n_current_interface_dofs());
    this->get_local_dof_values(fe_function, local_dof_values);

    get_function_values_from_local_dof_values(here_or_there,
                                              local_dof_values,
                                              values);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Vector<dim, spacedim>::get_jump_in_function_values_from_local_dof_values(
    const InputVector &local_dof_values,
    std::vector<solution_value_type<typename InputVector::value_type>> &values)
    const
  {
    AssertDimension(values.size(), this->fe_interface->n_quadrature_points);

    for (const auto dof_index : this->fe_interface->dof_indices())
      for (const auto q_index : this->fe_interface->quadrature_point_indices())
        {
          if (dof_index == 0)
            values[q_index] = 0.;

          values[q_index] +=
            local_dof_values[dof_index] * jump_in_values(dof_index, q_index);
        }
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Vector<dim, spacedim>::get_jump_in_function_values(
    const InputVector &fe_function,
    std::vector<solution_value_type<typename InputVector::value_type>> &values)
    const
  {
    std::vector<typename InputVector::value_type> local_dof_values(
      this->fe_interface->n_current_interface_dofs());
    this->get_local_dof_values(fe_function, local_dof_values);

    get_jump_in_function_values_from_local_dof_values(local_dof_values, values);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Vector<dim, spacedim>::get_jump_in_function_gradients_from_local_dof_values(
    const InputVector &local_dof_values,
    std::vector<solution_gradient_type<typename InputVector::value_type>>
      &gradients) const
  {
    AssertDimension(gradients.size(), this->fe_interface->n_quadrature_points);

    for (const auto dof_index : this->fe_interface->dof_indices())
      for (const auto q_index : this->fe_interface->quadrature_point_indices())
        {
          if (dof_index == 0)
            gradients[q_index] = 0.;

          gradients[q_index] +=
            local_dof_values[dof_index] * jump_in_gradients(dof_index, q_index);
        }
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Vector<dim, spacedim>::get_jump_in_function_gradients(
    const InputVector &fe_function,
    std::vector<solution_gradient_type<typename InputVector::value_type>>
      &gradients) const
  {
    std::vector<typename InputVector::value_type> local_dof_values(
      this->fe_interface->n_current_interface_dofs());
    this->get_local_dof_values(fe_function, local_dof_values);

    get_jump_in_function_gradients_from_local_dof_values(local_dof_values,
                                                         gradients);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Vector<dim, spacedim>::get_average_of_function_values_from_local_dof_values(
    const InputVector &local_dof_values,
    std::vector<solution_value_type<typename InputVector::value_type>> &values)
    const
  {
    AssertDimension(values.size(), this->fe_interface->n_quadrature_points);

    for (const auto dof_index : this->fe_interface->dof_indices())
      for (const auto q_index : this->fe_interface->quadrature_point_indices())
        {
          if (dof_index == 0)
            values[q_index] = 0.;

          values[q_index] +=
            local_dof_values[dof_index] * average_of_values(dof_index, q_index);
        }
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Vector<dim, spacedim>::get_average_of_function_values(
    const InputVector &fe_function,
    std::vector<solution_value_type<typename InputVector::value_type>> &values)
    const
  {
    std::vector<typename InputVector::value_type> local_dof_values(
      this->fe_interface->n_current_interface_dofs());
    this->get_local_dof_values(fe_function, local_dof_values);

    get_average_of_function_values_from_local_dof_values(local_dof_values,
                                                         values);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Vector<dim, spacedim>::
    get_average_of_function_gradients_from_local_dof_values(
      const InputVector &local_dof_values,
      std::vector<solution_gradient_type<typename InputVector::value_type>>
        &gradients) const
  {
    AssertDimension(gradients.size(), this->fe_interface->n_quadrature_points);

    for (const auto dof_index : this->fe_interface->dof_indices())
      for (const auto q_index : this->fe_interface->quadrature_point_indices())
        {
          if (dof_index == 0)
            gradients[q_index] = 0.;

          gradients[q_index] += local_dof_values[dof_index] *
                                average_of_gradients(dof_index, q_index);
        }
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Vector<dim, spacedim>::get_average_of_function_gradients(
    const InputVector &fe_function,
    std::vector<solution_gradient_type<typename InputVector::value_type>>
      &gradients) const
  {
    std::vector<typename InputVector::value_type> local_dof_values(
      this->fe_interface->n_current_interface_dofs());
    this->get_local_dof_values(fe_function, local_dof_values);

    get_average_of_function_gradients_from_local_dof_values(local_dof_values,
                                                            gradients);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Vector<dim, spacedim>::get_jump_in_function_hessians_from_local_dof_values(
    const InputVector &local_dof_values,
    std::vector<solution_hessian_type<typename InputVector::value_type>>
      &hessians) const
  {
    AssertDimension(hessians.size(), this->fe_interface->n_quadrature_points);

    for (const auto dof_index : this->fe_interface->dof_indices())
      for (const auto q_index : this->fe_interface->quadrature_point_indices())
        {
          if (dof_index == 0)
            hessians[q_index] = 0.;

          hessians[q_index] +=
            local_dof_values[dof_index] * jump_in_hessians(dof_index, q_index);
        }
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Vector<dim, spacedim>::get_jump_in_function_hessians(
    const InputVector &fe_function,
    std::vector<solution_hessian_type<typename InputVector::value_type>>
      &hessians) const
  {
    std::vector<typename InputVector::value_type> local_dof_values(
      this->fe_interface->n_current_interface_dofs());
    this->get_local_dof_values(fe_function, local_dof_values);

    get_jump_in_function_hessians_from_local_dof_values(local_dof_values,
                                                        hessians);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Vector<dim, spacedim>::get_average_of_function_hessians_from_local_dof_values(
    const InputVector &local_dof_values,
    std::vector<solution_hessian_type<typename InputVector::value_type>>
      &hessians) const
  {
    AssertDimension(hessians.size(), this->fe_interface->n_quadrature_points);

    for (const unsigned int dof_index : this->fe_interface->dof_indices())
      for (const auto q_index : this->fe_interface->quadrature_point_indices())
        {
          if (dof_index == 0)
            hessians[q_index] = 0.;

          hessians[q_index] += local_dof_values[dof_index] *
                               average_of_hessians(dof_index, q_index);
        }
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Vector<dim, spacedim>::get_average_of_function_hessians(
    const InputVector &fe_function,
    std::vector<solution_hessian_type<typename InputVector::value_type>>
      &hessians) const
  {
    std::vector<typename InputVector::value_type> local_dof_values(
      this->fe_interface->n_current_interface_dofs());
    this->get_local_dof_values(fe_function, local_dof_values);

    get_average_of_function_hessians_from_local_dof_values(local_dof_values,
                                                           hessians);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Vector<dim, spacedim>::
    get_jump_in_function_third_derivatives_from_local_dof_values(
      const InputVector &local_dof_values,
      std::vector<
        solution_third_derivative_type<typename InputVector::value_type>>
        &third_derivatives) const
  {
    AssertDimension(third_derivatives.size(),
                    this->fe_interface->n_quadrature_points);

    for (const unsigned int dof_index : this->fe_interface->dof_indices())
      for (const auto q_index : this->fe_interface->quadrature_point_indices())
        {
          if (dof_index == 0)
            third_derivatives[q_index] = 0.;

          third_derivatives[q_index] +=
            local_dof_values[dof_index] *
            jump_in_third_derivatives(dof_index, q_index);
        }
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Vector<dim, spacedim>::get_jump_in_function_third_derivatives(
    const InputVector &fe_function,
    std::vector<
      solution_third_derivative_type<typename InputVector::value_type>>
      &third_derivatives) const
  {
    std::vector<typename InputVector::value_type> local_dof_values(
      this->fe_interface->n_current_interface_dofs());
    this->get_local_dof_values(fe_function, local_dof_values);

    get_jump_in_function_third_derivatives_from_local_dof_values(
      local_dof_values, third_derivatives);
  }
} // namespace FEInterfaceViews

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
