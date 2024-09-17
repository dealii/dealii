// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_values_views_h
#define dealii_fe_values_views_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/lazy.h>
#include <deal.II/base/observer_pointer.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/fe_values_extractors.h>

#include <deal.II/lac/read_vector.h>

#include <type_traits>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
template <int dim, int spacedim = dim>
class FEValuesBase;
#endif

namespace internal
{
  /**
   * A class whose specialization is used to define what type the curl of a
   * vector valued function corresponds to.
   */
  template <int dim, typename NumberType = double>
  struct CurlType;

  /**
   * A class whose specialization is used to define what type the curl of a
   * vector valued function corresponds to.
   *
   * In 1d, the curl is a scalar.
   */
  template <typename NumberType>
  struct CurlType<1, NumberType>
  {
    using type = Tensor<1, 1, NumberType>;
  };

  /**
   * A class whose specialization is used to define what type the curl of a
   * vector valued function corresponds to.
   *
   * In 2d, the curl is a scalar.
   */
  template <typename NumberType>
  struct CurlType<2, NumberType>
  {
    using type = Tensor<1, 1, NumberType>;
  };

  /**
   * A class whose specialization is used to define what type the curl of a
   * vector valued function corresponds to.
   *
   * In 3d, the curl is a vector.
   */
  template <typename NumberType>
  struct CurlType<3, NumberType>
  {
    using type = Tensor<1, 3, NumberType>;
  };
} // namespace internal



/**
 * A namespace for "views" on a FEValues, FEFaceValues, or FESubfaceValues
 * object. A view represents only a certain part of the whole: whereas the
 * FEValues object represents <i>all</i> values, gradients, or second
 * derivatives of all components of a vector-valued element, views restrict
 * the attention to only a single component or a subset of components. You
 * typically get objects of classes defined in this namespace by applying
 * FEValuesExtractors objects to a FEValues, FEFaceValues or FESubfaceValues
 * objects using the square bracket operator.
 *
 * There are classes that present views for single scalar components, vector
 * components consisting of <code>dim</code> elements, and symmetric second
 * order tensor components consisting of <code>(dim*dim + dim)/2</code>
 * elements
 *
 * See the description of the
 * @ref vector_valued
 * topic for examples how to use the features of this namespace.
 *
 * @ingroup feaccess vector_valued
 */
namespace FEValuesViews
{
  /**
   * A class representing a view to a single scalar component of a possibly
   * vector-valued finite element. Views are discussed in the
   * @ref vector_valued
   * topic.
   *
   * You get an object of this type if you apply a FEValuesExtractors::Scalar
   * to an FEValues, FEFaceValues or FESubfaceValues object.
   *
   * @ingroup feaccess vector_valued
   */
  template <int dim, int spacedim = dim>
  class Scalar
  {
  public:
    /**
     * An alias for the data type of values of the view this class
     * represents. Since we deal with a single components, the value type is a
     * scalar double.
     */
    using value_type = double;

    /**
     * An alias for the type of gradients of the view this class represents.
     * Here, for a scalar component of the finite element, the gradient is a
     * <code>Tensor@<1,dim@></code>.
     */
    using gradient_type = dealii::Tensor<1, spacedim>;

    /**
     * An alias for the type of second derivatives of the view this class
     * represents. Here, for a scalar component of the finite element, the
     * Hessian is a <code>Tensor@<2,dim@></code>.
     */
    using hessian_type = dealii::Tensor<2, spacedim>;

    /**
     * An alias for the type of third derivatives of the view this class
     * represents. Here, for a scalar component of the finite element, the
     * Third derivative is a <code>Tensor@<3,dim@></code>.
     */
    using third_derivative_type = dealii::Tensor<3, spacedim>;

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
     * laplacians of the view this class provides. This is the data type of
     * scalar components of a finite element field whose degrees of
     * freedom are described by a vector with elements of type @p Number.
     */
    template <typename Number>
    using solution_laplacian_type =
      typename ProductType<Number, value_type>::type;

    /**
     * An alias for the data type of the product of a @p Number and the
     * hessians of the view this class provides. This is the data type of
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
     * A structure where for each shape function we pre-compute a bunch of
     * data that will make later accesses much cheaper.
     */
    struct ShapeFunctionData
    {
      /**
       * For each shape function, store whether the selected vector component
       * may be nonzero. For primitive shape functions we know for sure
       * whether a certain scalar component of a given shape function is
       * nonzero, whereas for non-primitive shape functions this may not be
       * entirely clear (e.g. for RT elements it depends on the shape of a
       * cell).
       */
      bool is_nonzero_shape_function_component;

      /**
       * For each shape function, store the row index within the shape_values,
       * shape_gradients, and shape_hessians tables (the column index is the
       * quadrature point index). If the shape function is primitive, then we
       * can get this information from the shape_function_to_row_table of the
       * FEValues object; otherwise, we have to work a bit harder to compute
       * this information.
       */
      unsigned int row_index;
    };

    /**
     * Default constructor. Creates an invalid object.
     */
    Scalar();

    /**
     * Constructor for an object that represents a single scalar component of
     * a FEValuesBase object (or of one of the classes derived from
     * FEValuesBase).
     */
    Scalar(const FEValuesBase<dim, spacedim> &fe_values_base,
           const unsigned int                 component);

    /**
     * Copy constructor. This is not a lightweight object so we don't allow
     * copying and generate a compile-time error if this function is called.
     */
    Scalar(const Scalar<dim, spacedim> &) = delete;

    /**
     * Move constructor.
     */
    // NOLINTNEXTLINE OSX does not compile with noexcept
    Scalar(Scalar<dim, spacedim> &&) = default;

    /**
     * Destructor.
     */
    ~Scalar() = default;

    /**
     * Copy operator. This is not a lightweight object so we don't allow
     * copying and generate a compile-time error if this function is called.
     */
    Scalar &
    operator=(const Scalar<dim, spacedim> &) = delete;

    /**
     * Move assignment operator.
     */
    Scalar &
    operator=(Scalar<dim, spacedim> &&) noexcept = default;

    /**
     * Return the value of the vector component selected by this view, for the
     * shape function and quadrature point selected by the arguments.
     *
     * @param shape_function Number of the shape function to be evaluated.
     * Note that this number runs from zero to dofs_per_cell, even in the case
     * of an FEFaceValues or FESubfaceValues object.
     *
     * @param q_point Number of the quadrature point at which function is to
     * be evaluated.
     *
     * @dealiiRequiresUpdateFlags{update_values}
     */
    value_type
    value(const unsigned int shape_function, const unsigned int q_point) const;

    /**
     * Return the gradient (a tensor of rank 1) of the vector component
     * selected by this view, for the shape function and quadrature point
     * selected by the arguments.
     *
     * @note The meaning of the arguments is as documented for the value()
     * function.
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    gradient_type
    gradient(const unsigned int shape_function,
             const unsigned int q_point) const;

    /**
     * Return the Hessian (the tensor of rank 2 of all second derivatives) of
     * the vector component selected by this view, for the shape function and
     * quadrature point selected by the arguments.
     *
     * @note The meaning of the arguments is as documented for the value()
     * function.
     *
     * @dealiiRequiresUpdateFlags{update_hessians}
     */
    hessian_type
    hessian(const unsigned int shape_function,
            const unsigned int q_point) const;

    /**
     * Return the tensor of rank 3 of all third derivatives of the vector
     * component selected by this view, for the shape function and quadrature
     * point selected by the arguments.
     *
     * @note The meaning of the arguments is as documented for the value()
     * function.
     *
     * @dealiiRequiresUpdateFlags{update_third_derivatives}
     */
    third_derivative_type
    third_derivative(const unsigned int shape_function,
                     const unsigned int q_point) const;

    /**
     * Return the values of the selected scalar component of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called.
     *
     * This function is the equivalent of the
     * FEValuesBase::get_function_values function but it only works on the
     * selected scalar component.
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
     * Same as above, but using a vector of local degree-of-freedom values. In
     * other words, instead of extracting the nodal values of the degrees of
     * freedom located on the current cell from a global vector associated with
     * a DoFHandler object (as the function above does), this function instead
     * takes these local nodal values through its first argument. A typical
     * way to obtain such a vector is by calling code such as
     * @code
     *   cell->get_dof_values (dof_values, local_dof_values);
     * @endcode
     * (See DoFCellAccessor::get_dof_values() for more information on this
     * function.) The point of the current function is then that one could
     * modify these local values first, for example by applying a limiter
     * or by ensuring that all nodal values are positive, before evaluating
     * the finite element field that corresponds to these local values on the
     * current cell. Another application is where one wants to postprocess
     * the solution on a cell into a different finite element space on every
     * cell, without actually creating a corresponding DoFHandler -- in that
     * case, all one would compute is a local representation of that
     * postprocessed function, characterized by its nodal values; this function
     * then allows the evaluation of that representation at quadrature points.
     *
     * @param[in] dof_values A vector of local nodal values. This vector must
     *   have a length equal to number of DoFs on the current cell, and must
     *   be ordered in the same order as degrees of freedom are numbered on
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
      const InputVector &dof_values,
      std::vector<solution_value_type<typename InputVector::value_type>>
        &values) const;

    /**
     * Return the gradients of the selected scalar component of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called.
     *
     * This function is the equivalent of the
     * FEValuesBase::get_function_gradients function but it only works on the
     * selected scalar component.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the gradients of shape functions (i.e., @p gradient_type)
     * times the type used to store the values of the unknowns $U_j$ of your
     * finite element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    template <typename Number>
    void
    get_function_gradients(
      const ReadVector<Number>                    &fe_function,
      std::vector<solution_gradient_type<Number>> &gradients) const;

    /**
     * This function relates to get_function_gradients() in the same way
     * as get_function_values_from_local_dof_values() relates to
     * get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more information.
     */
    template <class InputVector>
    void
    get_function_gradients_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_gradient_type<typename InputVector::value_type>>
        &gradients) const;

    /**
     * Return the Hessians of the selected scalar component of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called.
     *
     * This function is the equivalent of the
     * FEValuesBase::get_function_hessians function but it only works on the
     * selected scalar component.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the Hessians of shape functions (i.e., @p hessian_type) times
     * the type used to store the values of the unknowns $U_j$ of your finite
     * element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_hessians}
     */
    template <typename Number>
    void
    get_function_hessians(
      const ReadVector<Number>                   &fe_function,
      std::vector<solution_hessian_type<Number>> &hessians) const;

    /**
     * This function relates to get_function_hessians() in the same way
     * as get_function_values_from_local_dof_values() relates to
     * get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more information.
     */
    template <class InputVector>
    void
    get_function_hessians_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_hessian_type<typename InputVector::value_type>>
        &hessians) const;


    /**
     * Return the Laplacians of the selected scalar component of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called. The
     * Laplacians are the trace of the Hessians.
     *
     * This function is the equivalent of the
     * FEValuesBase::get_function_laplacians function but it only works on the
     * selected scalar component.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the Laplacians of shape functions (i.e., @p value_type) times
     * the type used to store the values of the unknowns $U_j$ of your finite
     * element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_hessians}
     */
    template <typename Number>
    void
    get_function_laplacians(
      const ReadVector<Number>                     &fe_function,
      std::vector<solution_laplacian_type<Number>> &laplacians) const;

    /**
     * This function relates to get_function_laplacians() in the same way
     * as get_function_values_from_local_dof_values() relates to
     * get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more information.
     */
    template <class InputVector>
    void
    get_function_laplacians_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_laplacian_type<typename InputVector::value_type>>
        &laplacians) const;


    /**
     * Return the third derivatives of the selected scalar component of the
     * finite element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called.
     *
     * This function is the equivalent of the
     * FEValuesBase::get_function_third_derivatives function but it only works
     * on the selected scalar component.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the third derivatives of shape functions (i.e., @p
     * third_derivative_type) times the type used to store the values of the
     * unknowns $U_j$ of your finite element vector $U$ (represented by the @p
     * fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_third_derivatives}
     */
    template <typename Number>
    void
    get_function_third_derivatives(
      const ReadVector<Number>                            &fe_function,
      std::vector<solution_third_derivative_type<Number>> &third_derivatives)
      const;

    /**
     * This function relates to get_function_third_derivatives() in the same way
     * as get_function_values_from_local_dof_values() relates to
     * get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more information.
     */
    template <class InputVector>
    void
    get_function_third_derivatives_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<
        solution_third_derivative_type<typename InputVector::value_type>>
        &third_derivatives) const;


  private:
    /**
     * A pointer to the FEValuesBase object we operate on.
     */
    ObserverPointer<const FEValuesBase<dim, spacedim>> fe_values;

    /**
     * The single scalar component this view represents of the FEValuesBase
     * object.
     */
    unsigned int component;

    /**
     * Store the data about shape functions.
     */
    std::vector<ShapeFunctionData> shape_function_data;
  };



  /**
   * A class representing a view to a set of <code>spacedim</code> components
   * forming a vector part of a vector-valued finite element. Views are
   * discussed in the
   * @ref vector_valued
   * topic.
   *
   * Note that in the current context, a vector is meant in the sense physics
   * uses it: it has <code>spacedim</code> components that behave in specific
   * ways under coordinate system transformations. Examples include velocity
   * or displacement fields. This is opposed to how mathematics uses the word
   * "vector" (and how we use this word in other contexts in the library, for
   * example in the Vector class), where it really stands for a collection of
   * numbers. An example of this latter use of the word could be the set of
   * concentrations of chemical species in a flame; however, these are really
   * just a collection of scalar variables, since they do not change if the
   * coordinate system is rotated, unlike the components of a velocity vector,
   * and consequently, this class should not be used for this context.
   *
   * This class allows to query the value, gradient and divergence of
   * (components of) shape functions and solutions representing vectors. The
   * gradient of a vector $d_{k}, 0\le k<\text{dim}$ is defined as $S_{ij} =
   * \frac{\partial d_{i}}{\partial x_j}, 0\le i,j<\text{dim}$.
   *
   * You get an object of this type if you apply a FEValuesExtractors::Vector
   * to an FEValues, FEFaceValues or FESubfaceValues object.
   *
   * @ingroup feaccess vector_valued
   */
  template <int dim, int spacedim = dim>
  class Vector
  {
  public:
    /**
     * An alias for the data type of values of the view this class
     * represents. Since we deal with a set of <code>dim</code> components,
     * the value type is a Tensor<1,spacedim>.
     */
    using value_type = dealii::Tensor<1, spacedim>;

    /**
     * An alias for the type of gradients of the view this class represents.
     * Here, for a set of <code>dim</code> components of the finite element,
     * the gradient is a <code>Tensor@<2,spacedim@></code>.
     *
     * See the general documentation of this class for how exactly the
     * gradient of a vector is defined.
     */
    using gradient_type = dealii::Tensor<2, spacedim>;

    /**
     * An alias for the type of symmetrized gradients of the view this class
     * represents. Here, for a set of <code>dim</code> components of the
     * finite element, the symmetrized gradient is a
     * <code>SymmetricTensor@<2,spacedim@></code>.
     *
     * The symmetric gradient of a vector field $\mathbf v$ is defined as
     * $\varepsilon(\mathbf v)=\frac 12 (\nabla \mathbf v + \nabla \mathbf
     * v^T)$.
     */
    using symmetric_gradient_type = dealii::SymmetricTensor<2, spacedim>;

    /**
     * An alias for the type of the divergence of the view this class
     * represents. Here, for a set of <code>dim</code> components of the
     * finite element, the divergence of course is a scalar.
     */
    using divergence_type = double;

    /**
     * An alias for the type of the curl of the view this class represents.
     * Here, for a set of <code>spacedim=2</code> components of the finite
     * element, the curl is a <code>Tensor@<1, 1@></code>. For
     * <code>spacedim=3</code> it is a <code>Tensor@<1, dim@></code>.
     */
    using curl_type = typename dealii::internal::CurlType<spacedim>::type;

    /**
     * An alias for the type of second derivatives of the view this class
     * represents. Here, for a set of <code>dim</code> components of the
     * finite element, the Hessian is a <code>Tensor@<3,dim@></code>.
     */
    using hessian_type = dealii::Tensor<3, spacedim>;

    /**
     * An alias for the type of third derivatives of the view this class
     * represents. Here, for a set of <code>dim</code> components of the
     * finite element, the third derivative is a <code>Tensor@<4,dim@></code>.
     */
    using third_derivative_type = dealii::Tensor<4, spacedim>;

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
     * symmetric gradients of the view this class provides. This is the data
     * type of vector components of a finite element field whose degrees of
     * freedom are described by a vector with elements of type @p Number.
     */
    template <typename Number>
    using solution_symmetric_gradient_type =
      typename ProductType<Number, symmetric_gradient_type>::type;

    /**
     * An alias for the data type of the product of a @p Number and the
     * divergences of the view this class provides. This is the data type of
     * vector components of a finite element field whose degrees of
     * freedom are described by a vector with elements of type @p Number.
     */
    template <typename Number>
    using solution_divergence_type =
      typename ProductType<Number, divergence_type>::type;

    /**
     * An alias for the data type of the product of a @p Number and the
     * laplacians of the view this class provides. This is the data type of
     * vector components of a finite element field whose degrees of
     * freedom are described by a vector with elements of type @p Number.
     */
    template <typename Number>
    using solution_laplacian_type =
      typename ProductType<Number, value_type>::type;

    /**
     * An alias for the data type of the product of a @p Number and the
     * curls of the view this class provides. This is the data type of
     * vector components of a finite element field whose degrees of
     * freedom are described by a vector with elements of type @p Number.
     */
    template <typename Number>
    using solution_curl_type = typename ProductType<Number, curl_type>::type;

    /**
     * An alias for the data type of the product of a @p Number and the
     * hessians of the view this class provides. This is the data type of
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
     * A structure where for each shape function we pre-compute a bunch of
     * data that will make later accesses much cheaper.
     */
    struct ShapeFunctionData
    {
      /**
       * For each pair (shape function,component within vector), store whether
       * the selected vector component may be nonzero. For primitive shape
       * functions we know for sure whether a certain scalar component of a
       * given shape function is nonzero, whereas for non-primitive shape
       * functions this may not be entirely clear (e.g. for RT elements it
       * depends on the shape of a cell).
       */
      bool is_nonzero_shape_function_component[spacedim];

      /**
       * For each pair (shape function, component within vector), store the
       * row index within the shape_values, shape_gradients, and
       * shape_hessians tables (the column index is the quadrature point
       * index). If the shape function is primitive, then we can get this
       * information from the shape_function_to_row_table of the FEValues
       * object; otherwise, we have to work a bit harder to compute this
       * information.
       */
      unsigned int row_index[spacedim];

      /**
       * For each shape function say the following: if only a single entry in
       * is_nonzero_shape_function_component for this shape function is
       * nonzero, then store the corresponding value of row_index and
       * single_nonzero_component_index represents the index between 0 and dim
       * for which it is attained. If multiple components are nonzero, then
       * store -1. If no components are nonzero then store -2.
       */
      int          single_nonzero_component;
      unsigned int single_nonzero_component_index;
    };

    /**
     * Default constructor. Creates an invalid object.
     */
    Vector();

    /**
     * Constructor for an object that represents dim components of a
     * FEValuesBase object (or of one of the classes derived from
     * FEValuesBase), representing a vector-valued variable.
     *
     * The second argument denotes the index of the first component of the
     * selected vector.
     */
    Vector(const FEValuesBase<dim, spacedim> &fe_values_base,
           const unsigned int                 first_vector_component);

    /**
     * Copy constructor. This is not a lightweight object so we don't allow
     * copying and generate a compile-time error if this function is called.
     */
    Vector(const Vector<dim, spacedim> &) = delete;

    /**
     * Move constructor.
     */
    // NOLINTNEXTLINE OSX does not compile with noexcept
    Vector(Vector<dim, spacedim> &&) = default;

    /**
     * Destructor.
     */
    ~Vector() = default;

    /**
     * Copy operator. This is not a lightweight object so we don't allow
     * copying and generate a compile-time error if this function is called.
     */
    Vector &
    operator=(const Vector<dim, spacedim> &) = delete;

    /**
     * Move assignment operator.
     */
    // NOLINTNEXTLINE OSX does not compile with noexcept
    Vector &
    operator=(Vector<dim, spacedim> &&) = default; // NOLINT

    /**
     * Return the value of the vector components selected by this view, for
     * the shape function and quadrature point selected by the arguments.
     * Here, since the view represents a vector-valued part of the FEValues
     * object with <code>dim</code> components, the return type is a tensor of
     * rank 1 with <code>dim</code> components.
     *
     * @param shape_function Number of the shape function to be evaluated.
     * Note that this number runs from zero to dofs_per_cell, even in the case
     * of an FEFaceValues or FESubfaceValues object.
     *
     * @param q_point Number of the quadrature point at which function is to
     * be evaluated.
     *
     * @dealiiRequiresUpdateFlags{update_values}
     */
    value_type
    value(const unsigned int shape_function, const unsigned int q_point) const;

    /**
     * Return the gradient (a tensor of rank 2) of the vector component
     * selected by this view, for the shape function and quadrature point
     * selected by the arguments.
     *
     * See the general documentation of this class for how exactly the
     * gradient of a vector is defined.
     *
     * @note The meaning of the arguments is as documented for the value()
     * function.
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    gradient_type
    gradient(const unsigned int shape_function,
             const unsigned int q_point) const;

    /**
     * Return the symmetric gradient (a symmetric tensor of rank 2) of the
     * vector component selected by this view, for the shape function and
     * quadrature point selected by the arguments.
     *
     * The symmetric gradient is defined as $\frac 12 [(\nabla \phi_i(x_q)) +
     * (\nabla \phi_i(x_q))^T]$, where $\phi_i$ represents the
     * <code>dim</code> components selected from the FEValuesBase object, and
     * $x_q$ is the location of the $q$-th quadrature point.
     *
     * @note The meaning of the arguments is as documented for the value()
     * function.
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    symmetric_gradient_type
    symmetric_gradient(const unsigned int shape_function,
                       const unsigned int q_point) const;

    /**
     * Return the scalar divergence of the vector components selected by this
     * view, for the shape function and quadrature point selected by the
     * arguments.
     *
     * @note The meaning of the arguments is as documented for the value()
     * function.
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    divergence_type
    divergence(const unsigned int shape_function,
               const unsigned int q_point) const;

    /**
     * Return the vector curl of the vector components selected by this view,
     * for the shape function and quadrature point selected by the arguments.
     * For 1d this function does not make any sense. Thus it is not
     * implemented for <code>spacedim=1</code>.  In 2d the curl is defined as
     * @f{equation*}{
     * \operatorname{curl}(u) \dealcoloneq \frac{du_2}{dx} -\frac{du_1}{dy},
     * @f}
     * whereas in 3d it is given by
     * @f{equation*}{
     * \operatorname{curl}(u) \dealcoloneq \left( \begin{array}{c}
     * \frac{du_3}{dy}-\frac{du_2}{dz}\\ \frac{du_1}{dz}-\frac{du_3}{dx}\\
     * \frac{du_2}{dx}-\frac{du_1}{dy} \end{array} \right).
     * @f}
     *
     * @note The meaning of the arguments is as documented for the value()
     * function.
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    curl_type
    curl(const unsigned int shape_function, const unsigned int q_point) const;

    /**
     * Return the Hessian (the tensor of rank 2 of all second derivatives) of
     * the vector components selected by this view, for the shape function and
     * quadrature point selected by the arguments.
     *
     * @note The meaning of the arguments is as documented for the value()
     * function.
     *
     * @dealiiRequiresUpdateFlags{update_hessians}
     */
    hessian_type
    hessian(const unsigned int shape_function,
            const unsigned int q_point) const;

    /**
     * Return the tensor of rank 3 of all third derivatives of the vector
     * components selected by this view, for the shape function and quadrature
     * point selected by the arguments.
     *
     * @note The meaning of the arguments is as documented for the value()
     * function.
     *
     * @dealiiRequiresUpdateFlags{update_3rd_derivatives}
     */
    third_derivative_type
    third_derivative(const unsigned int shape_function,
                     const unsigned int q_point) const;

    /**
     * Return the values of the selected vector components of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called.
     *
     * This function is the equivalent of the
     * FEValuesBase::get_function_values function but it only works on the
     * selected vector components.
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
     * Same as above, but using a vector of local degree-of-freedom values. In
     * other words, instead of extracting the nodal values of the degrees of
     * freedom located on the current cell from a global vector associated with
     * a DoFHandler object (as the function above does), this function instead
     * takes these local nodal values through its first argument. A typical
     * way to obtain such a vector is by calling code such as
     * @code
     *   cell->get_dof_values (dof_values, local_dof_values);
     * @endcode
     * (See DoFCellAccessor::get_dof_values() for more information on this
     * function.) The point of the current function is then that one could
     * modify these local values first, for example by applying a limiter
     * or by ensuring that all nodal values are positive, before evaluating
     * the finite element field that corresponds to these local values on the
     * current cell. Another application is where one wants to postprocess
     * the solution on a cell into a different finite element space on every
     * cell, without actually creating a corresponding DoFHandler -- in that
     * case, all one would compute is a local representation of that
     * postprocessed function, characterized by its nodal values; this function
     * then allows the evaluation of that representation at quadrature points.
     *
     * @param[in] dof_values A vector of local nodal values. This vector must
     *   have a length equal to number of DoFs on the current cell, and must
     *   be ordered in the same order as degrees of freedom are numbered on
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
      const InputVector &dof_values,
      std::vector<solution_value_type<typename InputVector::value_type>>
        &values) const;

    /**
     * Return the gradients of the selected vector components of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called.
     *
     * This function is the equivalent of the
     * FEValuesBase::get_function_gradients function but it only works on the
     * selected vector components.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the gradients of shape functions (i.e., @p gradient_type)
     * times the type used to store the values of the unknowns $U_j$ of your
     * finite element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    template <typename Number>
    void
    get_function_gradients(
      const ReadVector<Number>                    &fe_function,
      std::vector<solution_gradient_type<Number>> &gradients) const;

    /**
     * This function relates to get_function_gradients() in the same way
     * as get_function_values_from_local_dof_values() relates to
     * get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more information.
     */
    template <class InputVector>
    void
    get_function_gradients_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_gradient_type<typename InputVector::value_type>>
        &gradients) const;

    /**
     * Return the symmetrized gradients of the selected vector components of
     * the finite element function characterized by <tt>fe_function</tt> at
     * the quadrature points of the cell, face or subface selected the last
     * time the <tt>reinit</tt> function of the FEValues object was called.
     *
     * The symmetric gradient of a vector field $\mathbf v$ is defined as
     * $\varepsilon(\mathbf v)=\frac 12 (\nabla \mathbf v + \nabla \mathbf
     * v^T)$.
     *
     * @note There is no equivalent function such as
     * FEValuesBase::get_function_symmetric_gradients in the FEValues classes
     * but the information can be obtained from
     * FEValuesBase::get_function_gradients, of course.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the symmetric gradients of shape functions (i.e., @p
     * symmetric_gradient_type) times the type used to store the values of the
     * unknowns $U_j$ of your finite element vector $U$ (represented by the @p
     * fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    template <typename Number>
    void
    get_function_symmetric_gradients(
      const ReadVector<Number> &fe_function,
      std::vector<solution_symmetric_gradient_type<Number>>
        &symmetric_gradients) const;

    /**
     * This function relates to get_function_symmetric_gradients() in the same
     * way as get_function_values_from_local_dof_values() relates to
     * get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more information.
     */
    template <class InputVector>
    void
    get_function_symmetric_gradients_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<
        solution_symmetric_gradient_type<typename InputVector::value_type>>
        &symmetric_gradients) const;

    /**
     * Return the divergence of the selected vector components of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called.
     *
     * There is no equivalent function such as
     * FEValuesBase::get_function_divergences in the FEValues classes but the
     * information can be obtained from FEValuesBase::get_function_gradients,
     * of course.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the divergences of shape functions (i.e., @p divergence_type)
     * times the type used to store the values of the unknowns $U_j$ of your
     * finite element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    template <typename Number>
    void
    get_function_divergences(
      const ReadVector<Number>                      &fe_function,
      std::vector<solution_divergence_type<Number>> &divergences) const;

    /**
     * This function relates to get_function_divergences() in the same way
     * as get_function_values_from_local_dof_values() relates to
     * get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more information.
     */
    template <class InputVector>
    void
    get_function_divergences_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_divergence_type<typename InputVector::value_type>>
        &divergences) const;

    /**
     * Return the curl of the selected vector components of the finite element
     * function characterized by <tt>fe_function</tt> at the quadrature points
     * of the cell, face or subface selected the last time the <tt>reinit</tt>
     * function of the FEValues object was called.
     *
     * There is no equivalent function such as
     * FEValuesBase::get_function_curls in the FEValues classes but the
     * information can be obtained from FEValuesBase::get_function_gradients,
     * of course.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the curls of shape functions (i.e., @p curl_type) times the
     * type used to store the values of the unknowns $U_j$ of your finite
     * element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    template <typename Number>
    void
    get_function_curls(const ReadVector<Number>                &fe_function,
                       std::vector<solution_curl_type<Number>> &curls) const;

    /**
     * This function relates to get_function_curls() in the same way
     * as get_function_values_from_local_dof_values() relates to
     * get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more information.
     */
    template <class InputVector>
    void
    get_function_curls_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_curl_type<typename InputVector::value_type>> &curls)
      const;

    /**
     * Return the Hessians of the selected vector components of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called.
     *
     * This function is the equivalent of the
     * FEValuesBase::get_function_hessians function but it only works on the
     * selected vector components.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the Hessians of shape functions (i.e., @p hessian_type) times
     * the type used to store the values of the unknowns $U_j$ of your finite
     * element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_hessians}
     */
    template <typename Number>
    void
    get_function_hessians(
      const ReadVector<Number>                   &fe_function,
      std::vector<solution_hessian_type<Number>> &hessians) const;

    /**
     * This function relates to get_function_hessians() in the same way
     * as get_function_values_from_local_dof_values() relates to
     * get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more information.
     */
    template <class InputVector>
    void
    get_function_hessians_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_hessian_type<typename InputVector::value_type>>
        &hessians) const;

    /**
     * Return the Laplacians of the selected vector components of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called. The
     * Laplacians are the trace of the Hessians.
     *
     * This function is the equivalent of the
     * FEValuesBase::get_function_laplacians function but it only works on the
     * selected vector components.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the Laplacians of shape functions (i.e., @p laplacian_type)
     * times the type used to store the values of the unknowns $U_j$ of your
     * finite element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_hessians}
     */
    template <typename Number>
    void
    get_function_laplacians(
      const ReadVector<Number>                     &fe_function,
      std::vector<solution_laplacian_type<Number>> &laplacians) const;

    /**
     * This function relates to get_function_laplacians() in the same way
     * as get_function_values_from_local_dof_values() relates to
     * get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more information.
     */
    template <class InputVector>
    void
    get_function_laplacians_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_laplacian_type<typename InputVector::value_type>>
        &laplacians) const;

    /**
     * Return the third derivatives of the selected scalar component of the
     * finite element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called.
     *
     * This function is the equivalent of the
     * FEValuesBase::get_function_third_derivatives function but it only works
     * on the selected scalar component.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the third derivatives of shape functions (i.e., @p
     * third_derivative_type) times the type used to store the values of the
     * unknowns $U_j$ of your finite element vector $U$ (represented by the @p
     * fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_third_derivatives}
     */
    template <typename Number>
    void
    get_function_third_derivatives(
      const ReadVector<Number>                            &fe_function,
      std::vector<solution_third_derivative_type<Number>> &third_derivatives)
      const;

    /**
     * This function relates to get_function_third_derivatives() in the same way
     * as get_function_values_from_local_dof_values() relates to
     * get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more information.
     */
    template <class InputVector>
    void
    get_function_third_derivatives_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<
        solution_third_derivative_type<typename InputVector::value_type>>
        &third_derivatives) const;

  private:
    /**
     * A pointer to the FEValuesBase object we operate on.
     */
    ObserverPointer<const FEValuesBase<dim, spacedim>> fe_values;

    /**
     * The first component of the vector this view represents of the
     * FEValuesBase object.
     */
    unsigned int first_vector_component;

    /**
     * Store the data about shape functions.
     */
    std::vector<ShapeFunctionData> shape_function_data;
  };


  template <int rank, int dim, int spacedim = dim>
  class SymmetricTensor;

  /**
   * A class representing a view to a set of <code>(dim*dim + dim)/2</code>
   * components forming a symmetric second-order tensor from a vector-valued
   * finite element. Views are discussed in the
   * @ref vector_valued
   * topic.
   *
   * This class allows to query the value and divergence of (components of)
   * shape functions and solutions representing symmetric tensors. The
   * divergence of a symmetric tensor $S_{ij}, 0\le i,j<\text{dim}$ is defined
   * as $d_i = \sum_j \frac{\partial S_{ij}}{\partial x_j}, 0\le
   * i<\text{dim}$, which due to the symmetry of the tensor is also $d_i =
   * \sum_j \frac{\partial S_{ji}}{\partial x_j}$.  In other words, it due to
   * the symmetry of $S$ it does not matter whether we apply the nabla
   * operator by row or by column to get the divergence.
   *
   * You get an object of this type if you apply a
   * FEValuesExtractors::SymmetricTensor to an FEValues, FEFaceValues or
   * FESubfaceValues object.
   *
   * @ingroup feaccess vector_valued
   */
  template <int dim, int spacedim>
  class SymmetricTensor<2, dim, spacedim>
  {
  public:
    /**
     * An alias for the data type of values of the view this class
     * represents. Since we deal with a set of <code>(dim*dim + dim)/2</code>
     * components (i.e. the unique components of a symmetric second-order
     * tensor), the value type is a SymmetricTensor<2,spacedim>.
     */
    using value_type = dealii::SymmetricTensor<2, spacedim>;

    /**
     * An alias for the type of the divergence of the view this class
     * represents. Here, for a set of <code>(dim*dim + dim)/2</code> unique
     * components of the finite element representing a symmetric second-order
     * tensor, the divergence of course is a * <code>Tensor@<1,dim@></code>.
     *
     * See the general discussion of this class for a definition of the
     * divergence.
     */
    using divergence_type = dealii::Tensor<1, spacedim>;

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
     * divergences of the view this class provides. This is the data type of
     * vector components of a finite element field whose degrees of
     * freedom are described by a vector with elements of type @p Number.
     */
    template <typename Number>
    using solution_divergence_type =
      typename ProductType<Number, divergence_type>::type;


    /**
     * A structure where for each shape function we pre-compute a bunch of
     * data that will make later accesses much cheaper.
     */
    struct ShapeFunctionData
    {
      /**
       * For each pair (shape function,component within vector), store whether
       * the selected vector component may be nonzero. For primitive shape
       * functions we know for sure whether a certain scalar component of a
       * given shape function is nonzero, whereas for non-primitive shape
       * functions this may not be entirely clear (e.g. for RT elements it
       * depends on the shape of a cell).
       */
      bool is_nonzero_shape_function_component
        [value_type::n_independent_components];

      /**
       * For each pair (shape function, component within vector), store the
       * row index within the shape_values, shape_gradients, and
       * shape_hessians tables (the column index is the quadrature point
       * index). If the shape function is primitive, then we can get this
       * information from the shape_function_to_row_table of the FEValues
       * object; otherwise, we have to work a bit harder to compute this
       * information.
       */
      unsigned int row_index[value_type::n_independent_components];

      /**
       * For each shape function say the following: if only a single entry in
       * is_nonzero_shape_function_component for this shape function is
       * nonzero, then store the corresponding value of row_index and
       * single_nonzero_component_index represents the index between 0 and
       * (dim^2 + dim)/2 for which it is attained. If multiple components are
       * nonzero, then store -1. If no components are nonzero then store -2.
       */
      int single_nonzero_component;

      /**
       * Index of the @p single_nonzero_component .
       */
      unsigned int single_nonzero_component_index;
    };

    /**
     * Default constructor. Creates an invalid object.
     */
    SymmetricTensor();

    /**
     * Constructor for an object that represents <code>(dim*dim +
     * dim)/2</code> components of a FEValuesBase object (or of one of the
     * classes derived from FEValuesBase), representing the unique components
     * comprising a symmetric second- order tensor valued variable.
     *
     * The second argument denotes the index of the first component of the
     * selected symmetric second order tensor.
     */
    SymmetricTensor(const FEValuesBase<dim, spacedim> &fe_values_base,
                    const unsigned int                 first_tensor_component);

    /**
     * Copy constructor. This is not a lightweight object so we don't allow
     * copying and generate a compile-time error if this function is called.
     */
    SymmetricTensor(const SymmetricTensor<2, dim, spacedim> &) = delete;

    /**
     * Move constructor.
     */
    // NOLINTNEXTLINE OSX does not compile with noexcept
    SymmetricTensor(SymmetricTensor<2, dim, spacedim> &&) = default;

    /**
     * Copy operator. This is not a lightweight object so we don't allow
     * copying and generate a compile-time error if this function is called.
     */
    SymmetricTensor &
    operator=(const SymmetricTensor<2, dim, spacedim> &) = delete;

    /**
     * Move assignment operator.
     */
    SymmetricTensor &
    operator=(SymmetricTensor<2, dim, spacedim> &&) noexcept = default;

    /**
     * Return the value of the vector components selected by this view, for
     * the shape function and quadrature point selected by the arguments.
     * Here, since the view represents a vector-valued part of the FEValues
     * object with <code>(dim*dim + dim)/2</code> components (the unique
     * components of a symmetric second-order tensor), the return type is a
     * symmetric tensor of rank 2.
     *
     * @param shape_function Number of the shape function to be evaluated.
     * Note that this number runs from zero to dofs_per_cell, even in the case
     * of an FEFaceValues or FESubfaceValues object.
     *
     * @param q_point Number of the quadrature point at which function is to
     * be evaluated.
     *
     * @dealiiRequiresUpdateFlags{update_values}
     */
    value_type
    value(const unsigned int shape_function, const unsigned int q_point) const;

    /**
     * Return the vector divergence of the vector components selected by this
     * view, for the shape function and quadrature point selected by the
     * arguments.
     *
     * See the general discussion of this class for a definition of the
     * divergence.
     *
     * @note The meaning of the arguments is as documented for the value()
     * function.
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    divergence_type
    divergence(const unsigned int shape_function,
               const unsigned int q_point) const;

    /**
     * Return the values of the selected vector components of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called.
     *
     * This function is the equivalent of the
     * FEValuesBase::get_function_values function but it only works on the
     * selected vector components.
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
     * Same as above, but using a vector of local degree-of-freedom values. In
     * other words, instead of extracting the nodal values of the degrees of
     * freedom located on the current cell from a global vector associated with
     * a DoFHandler object (as the function above does), this function instead
     * takes these local nodal values through its first argument. A typical
     * way to obtain such a vector is by calling code such as
     * @code
     *   cell->get_dof_values (dof_values, local_dof_values);
     * @endcode
     * (See DoFCellAccessor::get_dof_values() for more information on this
     * function.) The point of the current function is then that one could
     * modify these local values first, for example by applying a limiter
     * or by ensuring that all nodal values are positive, before evaluating
     * the finite element field that corresponds to these local values on the
     * current cell. Another application is where one wants to postprocess
     * the solution on a cell into a different finite element space on every
     * cell, without actually creating a corresponding DoFHandler -- in that
     * case, all one would compute is a local representation of that
     * postprocessed function, characterized by its nodal values; this function
     * then allows the evaluation of that representation at quadrature points.
     *
     * @param[in] dof_values A vector of local nodal values. This vector must
     *   have a length equal to number of DoFs on the current cell, and must
     *   be ordered in the same order as degrees of freedom are numbered on
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
      const InputVector &dof_values,
      std::vector<solution_value_type<typename InputVector::value_type>>
        &values) const;

    /**
     * Return the divergence of the selected vector components of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called.
     *
     * There is no equivalent function such as
     * FEValuesBase::get_function_divergences in the FEValues classes but the
     * information can be obtained from FEValuesBase::get_function_gradients,
     * of course.
     *
     * See the general discussion of this class for a definition of the
     * divergence.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the divergences of shape functions (i.e., @p divergence_type)
     * times the type used to store the values of the unknowns $U_j$ of your
     * finite element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    template <typename Number>
    void
    get_function_divergences(
      const ReadVector<Number>                      &fe_function,
      std::vector<solution_divergence_type<Number>> &divergences) const;

    /**
     * This function relates to get_function_divergences() in the same way
     * as get_function_values_from_local_dof_values() relates to
     * get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more information.
     */
    template <class InputVector>
    void
    get_function_divergences_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_divergence_type<typename InputVector::value_type>>
        &divergences) const;

  private:
    /**
     * A pointer to the FEValuesBase object we operate on.
     */
    ObserverPointer<const FEValuesBase<dim, spacedim>> fe_values;

    /**
     * The first component of the vector this view represents of the
     * FEValuesBase object.
     */
    unsigned int first_tensor_component;

    /**
     * Store the data about shape functions.
     */
    std::vector<ShapeFunctionData> shape_function_data;
  };


  template <int rank, int dim, int spacedim = dim>
  class Tensor;

  /**
   * A class representing a view to a set of <code>dim*dim</code> components
   * forming a second-order tensor from a vector-valued finite element. Views
   * are discussed in the
   * @ref vector_valued
   * topic.
   *
   * This class allows to query the value, gradient and divergence of
   * (components of) shape functions and solutions representing tensors. The
   * divergence of a tensor $T_{ij},\, 0\le i,j<\text{dim}$ is defined as $d_i =
   * \sum_j \frac{\partial T_{ij}}{\partial x_j}, \, 0\le i<\text{dim}$, whereas
   * its gradient is $G_{ijk} = \frac{\partial T_{ij}}{\partial x_k}$.
   *
   * You get an object of this type if you apply a FEValuesExtractors::Tensor
   * to an FEValues, FEFaceValues or FESubfaceValues object.
   *
   * @ingroup feaccess vector_valued
   */
  template <int dim, int spacedim>
  class Tensor<2, dim, spacedim>
  {
  public:
    /**
     * Data type for what you get when you apply an extractor of this kind to
     * a vector-valued finite element.
     */
    using value_type = dealii::Tensor<2, spacedim>;

    /**
     * Data type for taking the divergence of a tensor: a vector.
     */
    using divergence_type = dealii::Tensor<1, spacedim>;

    /**
     * Data type for taking the gradient of a second order tensor: a third order
     * tensor.
     */
    using gradient_type = dealii::Tensor<3, spacedim>;

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
     * divergences of the view this class provides. This is the data type of
     * vector components of a finite element field whose degrees of
     * freedom are described by a vector with elements of type @p Number.
     */
    template <typename Number>
    using solution_divergence_type =
      typename ProductType<Number, divergence_type>::type;

    /**
     * An alias for the data type of the product of a @p Number and the
     * gradient of the view this class provides. This is the data type of
     * vector components of a finite element field whose degrees of
     * freedom are described by a vector with elements of type @p Number.
     */
    template <typename Number>
    using solution_gradient_type =
      typename ProductType<Number, gradient_type>::type;


    /**
     * A structure where for each shape function we pre-compute a bunch of
     * data that will make later accesses much cheaper.
     */
    struct ShapeFunctionData
    {
      /**
       * For each pair (shape function,component within vector), store whether
       * the selected vector component may be nonzero. For primitive shape
       * functions we know for sure whether a certain scalar component of a
       * given shape function is nonzero, whereas for non-primitive shape
       * functions this may not be entirely clear (e.g. for RT elements it
       * depends on the shape of a cell).
       */
      bool is_nonzero_shape_function_component
        [value_type::n_independent_components];

      /**
       * For each pair (shape function, component within vector), store the
       * row index within the shape_values, shape_gradients, and
       * shape_hessians tables (the column index is the quadrature point
       * index). If the shape function is primitive, then we can get this
       * information from the shape_function_to_row_table of the FEValues
       * object; otherwise, we have to work a bit harder to compute this
       * information.
       */
      unsigned int row_index[value_type::n_independent_components];

      /**
       * For each shape function say the following: if only a single entry in
       * is_nonzero_shape_function_component for this shape function is
       * nonzero, then store the corresponding value of row_index and
       * single_nonzero_component_index represents the index between 0 and
       * (dim^2) for which it is attained. If multiple components are nonzero,
       * then store -1. If no components are nonzero then store -2.
       */
      int single_nonzero_component;

      /**
       * Index of the @p single_nonzero_component .
       */
      unsigned int single_nonzero_component_index;
    };

    /**
     * Default constructor. Creates an invalid object.
     */
    Tensor();

    /**
     * Copy constructor. This is not a lightweight object so we don't allow
     * copying and generate a compile-time error if this function is called.
     */
    Tensor(const Tensor<2, dim, spacedim> &) = delete;

    /**
     * Move constructor.
     */
    // NOLINTNEXTLINE OSX does not compile with noexcept
    Tensor(Tensor<2, dim, spacedim> &&) = default;

    /**
     * Destructor.
     */
    ~Tensor() = default;

    /**
     * Constructor for an object that represents <code>(dim*dim)</code>
     * components of a FEValuesBase object (or of one of the classes derived
     * from FEValuesBase), representing the unique components comprising a
     * second-order tensor valued variable.
     *
     * The second argument denotes the index of the first component of the
     * selected symmetric second order tensor.
     */
    Tensor(const FEValuesBase<dim, spacedim> &fe_values_base,
           const unsigned int                 first_tensor_component);


    /**
     * Copy operator. This is not a lightweight object so we don't allow
     * copying and generate a compile-time error if this function is called.
     */
    Tensor &
    operator=(const Tensor<2, dim, spacedim> &) = delete;

    /**
     * Move assignment operator.
     */
    Tensor &
    operator=(Tensor<2, dim, spacedim> &&) = default; // NOLINT

    /**
     * Return the value of the vector components selected by this view, for
     * the shape function and quadrature point selected by the arguments.
     * Here, since the view represents a vector-valued part of the FEValues
     * object with <code>(dim*dim)</code> components (the unique components of
     * a second-order tensor), the return type is a tensor of rank 2.
     *
     * @param shape_function Number of the shape function to be evaluated.
     * Note that this number runs from zero to dofs_per_cell, even in the case
     * of an FEFaceValues or FESubfaceValues object.
     *
     * @param q_point Number of the quadrature point at which function is to
     * be evaluated.
     *
     * @dealiiRequiresUpdateFlags{update_values}
     */
    value_type
    value(const unsigned int shape_function, const unsigned int q_point) const;

    /**
     * Return the vector divergence of the vector components selected by this
     * view, for the shape function and quadrature point selected by the
     * arguments.
     *
     * See the general discussion of this class for a definition of the
     * divergence.
     *
     * @note The meaning of the arguments is as documented for the value()
     * function.
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    divergence_type
    divergence(const unsigned int shape_function,
               const unsigned int q_point) const;

    /**
     * Return the gradient (3-rd order tensor) of the vector components selected
     * by this view, for the shape function and quadrature point selected by the
     * arguments.
     *
     * See the general discussion of this class for a definition of the
     * gradient.
     *
     * @note The meaning of the arguments is as documented for the value()
     * function.
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    gradient_type
    gradient(const unsigned int shape_function,
             const unsigned int q_point) const;

    /**
     * Return the values of the selected vector components of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called.
     *
     * This function is the equivalent of the
     * FEValuesBase::get_function_values function but it only works on the
     * selected vector components.
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
     * Same as above, but using a vector of local degree-of-freedom values. In
     * other words, instead of extracting the nodal values of the degrees of
     * freedom located on the current cell from a global vector associated with
     * a DoFHandler object (as the function above does), this function instead
     * takes these local nodal values through its first argument. A typical
     * way to obtain such a vector is by calling code such as
     * @code
     *   cell->get_dof_values (dof_values, local_dof_values);
     * @endcode
     * (See DoFCellAccessor::get_dof_values() for more information on this
     * function.) The point of the current function is then that one could
     * modify these local values first, for example by applying a limiter
     * or by ensuring that all nodal values are positive, before evaluating
     * the finite element field that corresponds to these local values on the
     * current cell. Another application is where one wants to postprocess
     * the solution on a cell into a different finite element space on every
     * cell, without actually creating a corresponding DoFHandler -- in that
     * case, all one would compute is a local representation of that
     * postprocessed function, characterized by its nodal values; this function
     * then allows the evaluation of that representation at quadrature points.
     *
     * @param[in] dof_values A vector of local nodal values. This vector must
     *   have a length equal to number of DoFs on the current cell, and must
     *   be ordered in the same order as degrees of freedom are numbered on
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
      const InputVector &dof_values,
      std::vector<solution_value_type<typename InputVector::value_type>>
        &values) const;

    /**
     * Return the divergence of the selected vector components of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called.
     *
     * There is no equivalent function such as
     * FEValuesBase::get_function_divergences in the FEValues classes but the
     * information can be obtained from FEValuesBase::get_function_gradients,
     * of course.
     *
     * See the general discussion of this class for a definition of the
     * divergence.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the divergences of shape functions (i.e., @p divergence_type)
     * times the type used to store the values of the unknowns $U_j$ of your
     * finite element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    template <typename Number>
    void
    get_function_divergences(
      const ReadVector<Number>                      &fe_function,
      std::vector<solution_divergence_type<Number>> &divergences) const;

    /**
     * This function relates to get_function_divergences() in the same way
     * as get_function_values_from_local_dof_values() relates to
     * get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more information.
     */
    template <class InputVector>
    void
    get_function_divergences_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_divergence_type<typename InputVector::value_type>>
        &divergences) const;

    /**
     * Return the gradient of the selected vector components of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called.
     *
     * See the general discussion of this class for a definition of the
     * gradient.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the gradients of shape functions (i.e., @p gradient_type)
     * times the type used to store the values of the unknowns $U_j$ of your
     * finite element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    template <typename Number>
    void
    get_function_gradients(
      const ReadVector<Number>                    &fe_function,
      std::vector<solution_gradient_type<Number>> &gradients) const;

    /**
     * This function relates to get_function_gradients() in the same way
     * as get_function_values_from_local_dof_values() relates to
     * get_function_values(). See the documentation of
     * get_function_values_from_local_dof_values() for more information.
     */
    template <class InputVector>
    void
    get_function_gradients_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_gradient_type<typename InputVector::value_type>>
        &gradients) const;

  private:
    /**
     * A pointer to the FEValuesBase object we operate on.
     */
    ObserverPointer<const FEValuesBase<dim, spacedim>> fe_values;

    /**
     * The first component of the vector this view represents of the
     * FEValuesBase object.
     */
    unsigned int first_tensor_component;

    /**
     * Store the data about shape functions.
     */
    std::vector<ShapeFunctionData> shape_function_data;
  };

} // namespace FEValuesViews


namespace internal
{
  namespace FEValuesViews
  {
    /**
     * A class whose specialization is used to define what FEValuesViews
     * object corresponds to the given FEValuesExtractors object.
     */
    template <int dim, int spacedim, typename Extractor>
    struct ViewType
    {};

    /**
     * A class whose specialization is used to define what FEValuesViews
     * object corresponds to the given FEValuesExtractors object.
     *
     * When using FEValuesExtractors::Scalar, the corresponding view is an
     * FEValuesViews::Scalar<dim, spacedim>.
     */
    template <int dim, int spacedim>
    struct ViewType<dim, spacedim, FEValuesExtractors::Scalar>
    {
      using type = typename dealii::FEValuesViews::Scalar<dim, spacedim>;
    };

    /**
     * A class whose specialization is used to define what FEValuesViews
     * object corresponds to the given FEValuesExtractors object.
     *
     * When using FEValuesExtractors::Vector, the corresponding view is an
     * FEValuesViews::Vector<dim, spacedim>.
     */
    template <int dim, int spacedim>
    struct ViewType<dim, spacedim, FEValuesExtractors::Vector>
    {
      using type = typename dealii::FEValuesViews::Vector<dim, spacedim>;
    };

    /**
     * A class whose specialization is used to define what FEValuesViews
     * object corresponds to the given FEValuesExtractors object.
     *
     * When using FEValuesExtractors::Tensor<rank>, the corresponding view is an
     * FEValuesViews::Tensor<rank, dim, spacedim>.
     */
    template <int dim, int spacedim, int rank>
    struct ViewType<dim, spacedim, FEValuesExtractors::Tensor<rank>>
    {
      using type = typename dealii::FEValuesViews::Tensor<rank, dim, spacedim>;
    };

    /**
     * A class whose specialization is used to define what FEValuesViews
     * object corresponds to the given FEValuesExtractors object.
     *
     * When using FEValuesExtractors::SymmetricTensor<rank>, the corresponding
     * view is an FEValuesViews::SymmetricTensor<rank, dim, spacedim>.
     */
    template <int dim, int spacedim, int rank>
    struct ViewType<dim, spacedim, FEValuesExtractors::SymmetricTensor<rank>>
    {
      using type =
        typename dealii::FEValuesViews::SymmetricTensor<rank, dim, spacedim>;
    };

    /**
     * A class objects of which store a collection of FEValuesViews::Scalar,
     * FEValuesViews::Vector, etc object. The FEValuesBase class uses it to
     * generate all possible Views classes upon construction time; we do this
     * at construction time since the Views classes cache some information and
     * are therefore relatively expensive to create.
     */
    template <int dim, int spacedim>
    struct Cache
    {
      /**
       * Caches for scalar and vector, and symmetric second-order tensor
       * valued views.
       */
      std::vector<Lazy<dealii::FEValuesViews::Scalar<dim, spacedim>>> scalars;
      std::vector<Lazy<dealii::FEValuesViews::Vector<dim, spacedim>>> vectors;
      std::vector<
        Lazy<dealii::FEValuesViews::SymmetricTensor<2, dim, spacedim>>>
        symmetric_second_order_tensors;
      std::vector<Lazy<dealii::FEValuesViews::Tensor<2, dim, spacedim>>>
        second_order_tensors;

      /**
       * Constructor.
       */
      Cache(const FEValuesBase<dim, spacedim> &fe_values);
    };
  } // namespace FEValuesViews
} // namespace internal

namespace FEValuesViews
{
  /**
   * A templated alias that associates to a given Extractor class
   * the corresponding view in FEValuesViews.
   */
  template <int dim, int spacedim, typename Extractor>
  using View = typename dealii::internal::FEValuesViews::
    ViewType<dim, spacedim, Extractor>::type;
} // namespace FEValuesViews

#ifndef DOXYGEN

/*---------------- Inline functions: namespace FEValuesViews -----------------*/

namespace FEValuesViews
{
  template <int dim, int spacedim>
  inline typename Scalar<dim, spacedim>::value_type
  Scalar<dim, spacedim>::value(const unsigned int shape_function,
                               const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(
      fe_values->update_flags & update_values,
      ((typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
        "update_values"))));

    // an adaptation of the FEValuesBase::shape_value_component function
    // except that here we know the component as fixed and we have
    // pre-computed and cached a bunch of information. See the comments there.
    if (shape_function_data[shape_function].is_nonzero_shape_function_component)
      return fe_values->finite_element_output.shape_values(
        shape_function_data[shape_function].row_index, q_point);
    else
      return 0;
  }



  template <int dim, int spacedim>
  inline typename Scalar<dim, spacedim>::gradient_type
  Scalar<dim, spacedim>::gradient(const unsigned int shape_function,
                                  const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));

    // an adaptation of the FEValuesBase::shape_grad_component
    // function except that here we know the component as fixed and we have
    // pre-computed and cached a bunch of information. See the comments there.
    if (shape_function_data[shape_function].is_nonzero_shape_function_component)
      return fe_values->finite_element_output
        .shape_gradients[shape_function_data[shape_function].row_index]
                        [q_point];
    else
      return gradient_type();
  }



  template <int dim, int spacedim>
  inline typename Scalar<dim, spacedim>::hessian_type
  Scalar<dim, spacedim>::hessian(const unsigned int shape_function,
                                 const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_hessians,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_hessians")));

    // an adaptation of the FEValuesBase::shape_hessian_component
    // function except that here we know the component as fixed and we have
    // pre-computed and cached a bunch of information. See the comments there.
    if (shape_function_data[shape_function].is_nonzero_shape_function_component)
      return fe_values->finite_element_output
        .shape_hessians[shape_function_data[shape_function].row_index][q_point];
    else
      return hessian_type();
  }



  template <int dim, int spacedim>
  inline typename Scalar<dim, spacedim>::third_derivative_type
  Scalar<dim, spacedim>::third_derivative(const unsigned int shape_function,
                                          const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_3rd_derivatives,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_3rd_derivatives")));

    // an adaptation of the FEValuesBase::shape_3rdderivative_component
    // function except that here we know the component as fixed and we have
    // pre-computed and cached a bunch of information. See the comments there.
    if (shape_function_data[shape_function].is_nonzero_shape_function_component)
      return fe_values->finite_element_output
        .shape_3rd_derivatives[shape_function_data[shape_function].row_index]
                              [q_point];
    else
      return third_derivative_type();
  }



  template <int dim, int spacedim>
  inline typename Vector<dim, spacedim>::value_type
  Vector<dim, spacedim>::value(const unsigned int shape_function,
                               const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_values,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_values")));

    // same as for the scalar case except that we have one more index
    const int snc =
      shape_function_data[shape_function].single_nonzero_component;
    if (snc == -2)
      return value_type();
    else if (snc != -1)
      {
        value_type return_value;
        return_value[shape_function_data[shape_function]
                       .single_nonzero_component_index] =
          fe_values->finite_element_output.shape_values(snc, q_point);
        return return_value;
      }
    else
      {
        value_type return_value;
        for (unsigned int d = 0; d < dim; ++d)
          if (shape_function_data[shape_function]
                .is_nonzero_shape_function_component[d])
            return_value[d] = fe_values->finite_element_output.shape_values(
              shape_function_data[shape_function].row_index[d], q_point);

        return return_value;
      }
  }



  template <int dim, int spacedim>
  inline typename Vector<dim, spacedim>::gradient_type
  Vector<dim, spacedim>::gradient(const unsigned int shape_function,
                                  const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));

    // same as for the scalar case except that we have one more index
    const int snc =
      shape_function_data[shape_function].single_nonzero_component;
    if (snc == -2)
      return gradient_type();
    else if (snc != -1)
      {
        gradient_type return_value;
        return_value[shape_function_data[shape_function]
                       .single_nonzero_component_index] =
          fe_values->finite_element_output.shape_gradients[snc][q_point];
        return return_value;
      }
    else
      {
        gradient_type return_value;
        for (unsigned int d = 0; d < dim; ++d)
          if (shape_function_data[shape_function]
                .is_nonzero_shape_function_component[d])
            return_value[d] =
              fe_values->finite_element_output.shape_gradients
                [shape_function_data[shape_function].row_index[d]][q_point];

        return return_value;
      }
  }



  template <int dim, int spacedim>
  inline typename Vector<dim, spacedim>::divergence_type
  Vector<dim, spacedim>::divergence(const unsigned int shape_function,
                                    const unsigned int q_point) const
  {
    // this function works like in the case above
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));

    // same as for the scalar case except that we have one more index
    const int snc =
      shape_function_data[shape_function].single_nonzero_component;
    if (snc == -2)
      return divergence_type();
    else if (snc != -1)
      return fe_values->finite_element_output
        .shape_gradients[snc][q_point][shape_function_data[shape_function]
                                         .single_nonzero_component_index];
    else
      {
        divergence_type return_value = 0;
        for (unsigned int d = 0; d < dim; ++d)
          if (shape_function_data[shape_function]
                .is_nonzero_shape_function_component[d])
            return_value +=
              fe_values->finite_element_output.shape_gradients
                [shape_function_data[shape_function].row_index[d]][q_point][d];

        return return_value;
      }
  }



  template <int dim, int spacedim>
  inline typename Vector<dim, spacedim>::curl_type
  Vector<dim, spacedim>::curl(const unsigned int shape_function,
                              const unsigned int q_point) const
  {
    // this function works like in the case above

    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));
    // same as for the scalar case except that we have one more index
    const int snc =
      shape_function_data[shape_function].single_nonzero_component;

    if (snc == -2)
      return curl_type();

    else
      switch (dim)
        {
          case 1:
            {
              Assert(false,
                     ExcMessage(
                       "Computing the curl in 1d is not a useful operation"));
              return curl_type();
            }

          case 2:
            {
              if (snc != -1)
                {
                  curl_type return_value;

                  // the single nonzero component can only be zero or one in 2d
                  if (shape_function_data[shape_function]
                        .single_nonzero_component_index == 0)
                    return_value[0] =
                      -1.0 * fe_values->finite_element_output
                               .shape_gradients[snc][q_point][1];
                  else
                    return_value[0] = fe_values->finite_element_output
                                        .shape_gradients[snc][q_point][0];

                  return return_value;
                }

              else
                {
                  curl_type return_value;

                  return_value[0] = 0.0;

                  if (shape_function_data[shape_function]
                        .is_nonzero_shape_function_component[0])
                    return_value[0] -=
                      fe_values->finite_element_output
                        .shape_gradients[shape_function_data[shape_function]
                                           .row_index[0]][q_point][1];

                  if (shape_function_data[shape_function]
                        .is_nonzero_shape_function_component[1])
                    return_value[0] +=
                      fe_values->finite_element_output
                        .shape_gradients[shape_function_data[shape_function]
                                           .row_index[1]][q_point][0];

                  return return_value;
                }
            }

          case 3:
            {
              if (snc != -1)
                {
                  curl_type return_value;

                  switch (shape_function_data[shape_function]
                            .single_nonzero_component_index)
                    {
                      case 0:
                        {
                          return_value[0] = 0;
                          return_value[1] = fe_values->finite_element_output
                                              .shape_gradients[snc][q_point][2];
                          return_value[2] =
                            -1.0 * fe_values->finite_element_output
                                     .shape_gradients[snc][q_point][1];
                          return return_value;
                        }

                      case 1:
                        {
                          return_value[0] =
                            -1.0 * fe_values->finite_element_output
                                     .shape_gradients[snc][q_point][2];
                          return_value[1] = 0;
                          return_value[2] = fe_values->finite_element_output
                                              .shape_gradients[snc][q_point][0];
                          return return_value;
                        }

                      default:
                        {
                          return_value[0] = fe_values->finite_element_output
                                              .shape_gradients[snc][q_point][1];
                          return_value[1] =
                            -1.0 * fe_values->finite_element_output
                                     .shape_gradients[snc][q_point][0];
                          return_value[2] = 0;
                          return return_value;
                        }
                    }
                }

              else
                {
                  curl_type return_value;

                  for (unsigned int i = 0; i < dim; ++i)
                    return_value[i] = 0.0;

                  if (shape_function_data[shape_function]
                        .is_nonzero_shape_function_component[0])
                    {
                      return_value[1] +=
                        fe_values->finite_element_output
                          .shape_gradients[shape_function_data[shape_function]
                                             .row_index[0]][q_point][2];
                      return_value[2] -=
                        fe_values->finite_element_output
                          .shape_gradients[shape_function_data[shape_function]
                                             .row_index[0]][q_point][1];
                    }

                  if (shape_function_data[shape_function]
                        .is_nonzero_shape_function_component[1])
                    {
                      return_value[0] -=
                        fe_values->finite_element_output
                          .shape_gradients[shape_function_data[shape_function]
                                             .row_index[1]][q_point][2];
                      return_value[2] +=
                        fe_values->finite_element_output
                          .shape_gradients[shape_function_data[shape_function]
                                             .row_index[1]][q_point][0];
                    }

                  if (shape_function_data[shape_function]
                        .is_nonzero_shape_function_component[2])
                    {
                      return_value[0] +=
                        fe_values->finite_element_output
                          .shape_gradients[shape_function_data[shape_function]
                                             .row_index[2]][q_point][1];
                      return_value[1] -=
                        fe_values->finite_element_output
                          .shape_gradients[shape_function_data[shape_function]
                                             .row_index[2]][q_point][0];
                    }

                  return return_value;
                }
            }
        }
    // should not end up here
    DEAL_II_ASSERT_UNREACHABLE();
    return curl_type();
  }



  template <int dim, int spacedim>
  inline typename Vector<dim, spacedim>::hessian_type
  Vector<dim, spacedim>::hessian(const unsigned int shape_function,
                                 const unsigned int q_point) const
  {
    // this function works like in the case above
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_hessians,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_hessians")));

    // same as for the scalar case except that we have one more index
    const int snc =
      shape_function_data[shape_function].single_nonzero_component;
    if (snc == -2)
      return hessian_type();
    else if (snc != -1)
      {
        hessian_type return_value;
        return_value[shape_function_data[shape_function]
                       .single_nonzero_component_index] =
          fe_values->finite_element_output.shape_hessians[snc][q_point];
        return return_value;
      }
    else
      {
        hessian_type return_value;
        for (unsigned int d = 0; d < dim; ++d)
          if (shape_function_data[shape_function]
                .is_nonzero_shape_function_component[d])
            return_value[d] =
              fe_values->finite_element_output.shape_hessians
                [shape_function_data[shape_function].row_index[d]][q_point];

        return return_value;
      }
  }



  template <int dim, int spacedim>
  inline typename Vector<dim, spacedim>::third_derivative_type
  Vector<dim, spacedim>::third_derivative(const unsigned int shape_function,
                                          const unsigned int q_point) const
  {
    // this function works like in the case above
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_3rd_derivatives,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_3rd_derivatives")));

    // same as for the scalar case except that we have one more index
    const int snc =
      shape_function_data[shape_function].single_nonzero_component;
    if (snc == -2)
      return third_derivative_type();
    else if (snc != -1)
      {
        third_derivative_type return_value;
        return_value[shape_function_data[shape_function]
                       .single_nonzero_component_index] =
          fe_values->finite_element_output.shape_3rd_derivatives[snc][q_point];
        return return_value;
      }
    else
      {
        third_derivative_type return_value;
        for (unsigned int d = 0; d < dim; ++d)
          if (shape_function_data[shape_function]
                .is_nonzero_shape_function_component[d])
            return_value[d] =
              fe_values->finite_element_output.shape_3rd_derivatives
                [shape_function_data[shape_function].row_index[d]][q_point];

        return return_value;
      }
  }



  namespace internal
  {
    /**
     * Return the symmetrized version of a tensor whose n'th row equals the
     * second argument, with all other rows equal to zero.
     */
    inline dealii::SymmetricTensor<2, 1>
    symmetrize_single_row(const unsigned int n, const dealii::Tensor<1, 1> &t)
    {
      AssertIndexRange(n, 1);
      (void)n;

      return {{t[0]}};
    }



    inline dealii::SymmetricTensor<2, 2>
    symmetrize_single_row(const unsigned int n, const dealii::Tensor<1, 2> &t)
    {
      switch (n)
        {
          case 0:
            {
              return {{t[0], 0, t[1] / 2}};
            }
          case 1:
            {
              return {{0, t[1], t[0] / 2}};
            }
          default:
            {
              AssertIndexRange(n, 2);
              return {};
            }
        }
    }



    inline dealii::SymmetricTensor<2, 3>
    symmetrize_single_row(const unsigned int n, const dealii::Tensor<1, 3> &t)
    {
      switch (n)
        {
          case 0:
            {
              return {{t[0], 0, 0, t[1] / 2, t[2] / 2, 0}};
            }
          case 1:
            {
              return {{0, t[1], 0, t[0] / 2, 0, t[2] / 2}};
            }
          case 2:
            {
              return {{0, 0, t[2], 0, t[0] / 2, t[1] / 2}};
            }
          default:
            {
              AssertIndexRange(n, 3);
              return {};
            }
        }
    }
  } // namespace internal



  template <int dim, int spacedim>
  inline typename Vector<dim, spacedim>::symmetric_gradient_type
  Vector<dim, spacedim>::symmetric_gradient(const unsigned int shape_function,
                                            const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));

    // same as for the scalar case except that we have one more index
    const int snc =
      shape_function_data[shape_function].single_nonzero_component;
    if (snc == -2)
      return symmetric_gradient_type();
    else if (snc != -1)
      return internal::symmetrize_single_row(
        shape_function_data[shape_function].single_nonzero_component_index,
        fe_values->finite_element_output.shape_gradients[snc][q_point]);
    else
      {
        gradient_type return_value;
        for (unsigned int d = 0; d < dim; ++d)
          if (shape_function_data[shape_function]
                .is_nonzero_shape_function_component[d])
            return_value[d] =
              fe_values->finite_element_output.shape_gradients
                [shape_function_data[shape_function].row_index[d]][q_point];

        return symmetrize(return_value);
      }
  }



  template <int dim, int spacedim>
  inline typename SymmetricTensor<2, dim, spacedim>::value_type
  SymmetricTensor<2, dim, spacedim>::value(const unsigned int shape_function,
                                           const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_values,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_values")));

    // similar to the vector case where we have more then one index and we need
    // to convert between unrolled and component indexing for tensors
    const int snc =
      shape_function_data[shape_function].single_nonzero_component;

    if (snc == -2)
      {
        // shape function is zero for the selected components
        return value_type();
      }
    else if (snc != -1)
      {
        value_type         return_value;
        const unsigned int comp =
          shape_function_data[shape_function].single_nonzero_component_index;
        return_value[value_type::unrolled_to_component_indices(comp)] =
          fe_values->finite_element_output.shape_values(snc, q_point);
        return return_value;
      }
    else
      {
        value_type return_value;
        for (unsigned int d = 0; d < value_type::n_independent_components; ++d)
          if (shape_function_data[shape_function]
                .is_nonzero_shape_function_component[d])
            return_value[value_type::unrolled_to_component_indices(d)] =
              fe_values->finite_element_output.shape_values(
                shape_function_data[shape_function].row_index[d], q_point);
        return return_value;
      }
  }



  template <int dim, int spacedim>
  inline typename SymmetricTensor<2, dim, spacedim>::divergence_type
  SymmetricTensor<2, dim, spacedim>::divergence(
    const unsigned int shape_function,
    const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));

    const int snc =
      shape_function_data[shape_function].single_nonzero_component;

    if (snc == -2)
      {
        // shape function is zero for the selected components
        return divergence_type();
      }
    else if (snc != -1)
      {
        // we have a single non-zero component when the symmetric tensor is
        // represented in unrolled form. this implies we potentially have
        // two non-zero components when represented in component form!  we
        // will only have one non-zero entry if the non-zero component lies on
        // the diagonal of the tensor.
        //
        // the divergence of a second-order tensor is a first order tensor.
        //
        // assume the second-order tensor is A with components A_{ij}.  then
        // A_{ij} = A_{ji} and there is only one (if diagonal) or two non-zero
        // entries in the tensorial representation.  define the
        // divergence as:
        // b_i \dealcoloneq \dfrac{\partial phi_{ij}}{\partial x_j}.
        // (which is incidentally also
        // b_j \dealcoloneq \dfrac{\partial phi_{ij}}{\partial x_i}).
        // In both cases, a sum is implied.
        //
        // Now, we know the nonzero component in unrolled form: it is indicated
        // by 'snc'. we can figure out which tensor components belong to this:
        const unsigned int comp =
          shape_function_data[shape_function].single_nonzero_component_index;
        const unsigned int ii =
          value_type::unrolled_to_component_indices(comp)[0];
        const unsigned int jj =
          value_type::unrolled_to_component_indices(comp)[1];

        // given the form of the divergence above, if ii=jj there is only a
        // single nonzero component of the full tensor and the gradient
        // equals
        // b_ii \dealcoloneq \dfrac{\partial phi_{ii,ii}}{\partial x_ii}.
        // all other entries of 'b' are zero
        //
        // on the other hand, if ii!=jj, then there are two nonzero entries in
        // the full tensor and
        // b_ii \dealcoloneq \dfrac{\partial phi_{ii,jj}}{\partial x_ii}.
        // b_jj \dealcoloneq \dfrac{\partial phi_{ii,jj}}{\partial x_jj}.
        // again, all other entries of 'b' are zero
        const dealii::Tensor<1, spacedim> &phi_grad =
          fe_values->finite_element_output.shape_gradients[snc][q_point];

        divergence_type return_value;
        return_value[ii] = phi_grad[jj];

        if (ii != jj)
          return_value[jj] = phi_grad[ii];

        return return_value;
      }
    else
      {
        DEAL_II_NOT_IMPLEMENTED();
        divergence_type return_value;
        return return_value;
      }
  }



  template <int dim, int spacedim>
  inline typename Tensor<2, dim, spacedim>::value_type
  Tensor<2, dim, spacedim>::value(const unsigned int shape_function,
                                  const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_values,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_values")));

    // similar to the vector case where we have more then one index and we need
    // to convert between unrolled and component indexing for tensors
    const int snc =
      shape_function_data[shape_function].single_nonzero_component;

    if (snc == -2)
      {
        // shape function is zero for the selected components
        return value_type();
      }
    else if (snc != -1)
      {
        value_type         return_value;
        const unsigned int comp =
          shape_function_data[shape_function].single_nonzero_component_index;
        const TableIndices<2> indices =
          dealii::Tensor<2, spacedim>::unrolled_to_component_indices(comp);
        return_value[indices] =
          fe_values->finite_element_output.shape_values(snc, q_point);
        return return_value;
      }
    else
      {
        value_type return_value;
        for (unsigned int d = 0; d < dim * dim; ++d)
          if (shape_function_data[shape_function]
                .is_nonzero_shape_function_component[d])
            {
              const TableIndices<2> indices =
                dealii::Tensor<2, spacedim>::unrolled_to_component_indices(d);
              return_value[indices] =
                fe_values->finite_element_output.shape_values(
                  shape_function_data[shape_function].row_index[d], q_point);
            }
        return return_value;
      }
  }



  template <int dim, int spacedim>
  inline typename Tensor<2, dim, spacedim>::divergence_type
  Tensor<2, dim, spacedim>::divergence(const unsigned int shape_function,
                                       const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));

    const int snc =
      shape_function_data[shape_function].single_nonzero_component;

    if (snc == -2)
      {
        // shape function is zero for the selected components
        return divergence_type();
      }
    else if (snc != -1)
      {
        // we have a single non-zero component when the tensor is
        // represented in unrolled form.
        //
        // the divergence of a second-order tensor is a first order tensor.
        //
        // assume the second-order tensor is A with components A_{ij},
        // then divergence is d_i := \frac{\partial A_{ij}}{\partial x_j}
        //
        // Now, we know the nonzero component in unrolled form: it is indicated
        // by 'snc'. we can figure out which tensor components belong to this:
        const unsigned int comp =
          shape_function_data[shape_function].single_nonzero_component_index;
        const TableIndices<2> indices =
          dealii::Tensor<2, spacedim>::unrolled_to_component_indices(comp);
        const unsigned int ii = indices[0];
        const unsigned int jj = indices[1];

        const dealii::Tensor<1, spacedim> &phi_grad =
          fe_values->finite_element_output.shape_gradients[snc][q_point];

        divergence_type return_value;
        // note that we contract \nabla from the right
        return_value[ii] = phi_grad[jj];

        return return_value;
      }
    else
      {
        DEAL_II_NOT_IMPLEMENTED();
        divergence_type return_value;
        return return_value;
      }
  }



  template <int dim, int spacedim>
  inline typename Tensor<2, dim, spacedim>::gradient_type
  Tensor<2, dim, spacedim>::gradient(const unsigned int shape_function,
                                     const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, fe_values->fe->n_dofs_per_cell());
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));

    const int snc =
      shape_function_data[shape_function].single_nonzero_component;

    if (snc == -2)
      {
        // shape function is zero for the selected components
        return gradient_type();
      }
    else if (snc != -1)
      {
        // we have a single non-zero component when the tensor is
        // represented in unrolled form.
        //
        // the gradient of a second-order tensor is a third order tensor.
        //
        // assume the second-order tensor is A with components A_{ij},
        // then gradient is B_{ijk} := \frac{\partial A_{ij}}{\partial x_k}
        //
        // Now, we know the nonzero component in unrolled form: it is indicated
        // by 'snc'. we can figure out which tensor components belong to this:
        const unsigned int comp =
          shape_function_data[shape_function].single_nonzero_component_index;
        const TableIndices<2> indices =
          dealii::Tensor<2, spacedim>::unrolled_to_component_indices(comp);
        const unsigned int ii = indices[0];
        const unsigned int jj = indices[1];

        const dealii::Tensor<1, spacedim> &phi_grad =
          fe_values->finite_element_output.shape_gradients[snc][q_point];

        gradient_type return_value;
        return_value[ii][jj] = phi_grad;

        return return_value;
      }
    else
      {
        DEAL_II_NOT_IMPLEMENTED();
        gradient_type return_value;
        return return_value;
      }
  }
} // namespace FEValuesViews

#endif

DEAL_II_NAMESPACE_CLOSE

#endif
