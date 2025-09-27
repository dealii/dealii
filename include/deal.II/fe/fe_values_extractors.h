// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_values_extractors_h
#define dealii_fe_values_extractors_h


#include <deal.II/base/config.h>

#include <deal.II/base/types.h>

#include <string>


DEAL_II_NAMESPACE_OPEN


/**
 * A namespace in which we declare "extractors", i.e. classes that when used
 * as subscripts in operator[] expressions on FEValues, FEFaceValues, and
 * FESubfaceValues objects extract certain components of a vector-valued
 * element. The result of applying an extractor to these objects is an object
 * with corresponding type from the namespace FEValuesViews. There are
 * extractors for single scalar components, vector components consisting of
 * <code>dim</code> elements, and second order symmetric tensors consisting of
 * <code>(dim*dim + dim)/2</code> components, as well as second order
 * nonsymmetric tensors.
 *
 * One can think of extractors as the equivalent of an index, or an index range.
 * In the case of scalar extractors (i.e., the FEValuesExtractors::Scalar
 * class), creating an object like (see step-20 for this use)
 * @code
 *   const FEValuesExtractors::Scalar pressure(dim);
 * @endcode
 * can be thought of as creating a single index with value `dim`. By
 * itself, an index does not know what it is an index to, so it takes
 * the equivalent of an array to extract anything. Consequently,
 * assume that there is a finite element with at least `dim+1` vector
 * components (as indeed there is in step-20), and an FEValues object
 * that operates on it, then writing
 * @code
 *   fe_values[pressure]
 * @endcode
 * results in an object that represents the shape functions of only the
 * `dim`th component of the overall element. In the example, these
 * would be the values of the pressure shape functions, or more precisely:
 * the (scalar) pressure values of all shape functions (even for shape
 * functions that are not associated with the pressure, but for example
 * the velocity). In the example above, the result of using
 * `operator[]` on the `fe_values` object as shown is of type
 * FEValuesViews::Scalar.
 *
 * Likewise, when using
 * @code
 *   const FEValuesExtractors::Vector velocities(0);
 * @endcode
 * then the object so created can be thought of as an <i>index range</i>,
 * starting at zero and extending exactly `dim` components on. In Matlab
 * notation, one could write this as `0:dim-1`. Then, writing
 * @code
 *   fe_values[velocities]
 * @endcode
 * will result in an object that represents the values of a subset of
 * exactly `dim` vector components of the overall finite element, in
 * much the same way as writing `array(3:7)` in Matlab would return
 * an array of length 5 that has been extracted from the original
 * array by looking at indices 3 through 7 (inclusive).
 *
 * See the description of the
 * @ref vector_valued
 * topic for examples how to use the features of this namespace.
 *
 * @ingroup feaccess vector_valued
 */
namespace FEValuesExtractors
{
  /**
   * Extractor for a single scalar component of a vector-valued element. The
   * result of applying an object of this type to an FEValues, FEFaceValues or
   * FESubfaceValues object is of type FEValuesViews::Scalar. The concept of
   * extractors is defined in the documentation of the namespace
   * FEValuesExtractors and in the
   * @ref vector_valued
   * topic.
   *
   * @ingroup feaccess vector_valued
   */
  struct Scalar
  {
    /**
     * The selected scalar component of the vector.
     */
    unsigned int component;

    /**
     * Default constructor. Initialize the object with an invalid component.
     * This leads to an object that can not be used, but it allows objects of
     * this kind to be put into arrays that require a default constructor upon
     * resizing the array, and then later assigning a suitable object to each
     * element of the array.
     */
    constexpr Scalar();

    /**
     * Constructor. Take the selected vector component as argument.
     */
    constexpr Scalar(const unsigned int component);

    /**
     * Return a string that uniquely identifies this finite element extractor.
     */
    std::string
    get_name() const;
  };


  /**
   * Extractor for a vector of <code>spacedim</code> components of a
   * vector-valued element. The value of <code>spacedim</code> is defined by the
   * FEValues object the extractor is applied to. The result of applying an
   * object of this type to an FEValues, FEFaceValues or FESubfaceValues
   * object is of type FEValuesViews::Vector.
   *
   * The concept of extractors is defined in the documentation of the
   * namespace FEValuesExtractors and in the
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
   * @ingroup feaccess vector_valued
   */
  struct Vector
  {
    /**
     * The first component of the vector view.
     */
    unsigned int first_vector_component;

    /**
     * Default constructor. Initialize the object with an invalid component.
     * This leads to an object that can not be used, but it allows objects of
     * this kind to be put into arrays that require a default constructor upon
     * resizing the array, and then later assigning a suitable object to each
     * element of the array.
     */
    constexpr Vector();

    /**
     * Constructor. Take the first component of the selected vector inside the
     * FEValues object as argument.
     */
    constexpr Vector(const unsigned int first_vector_component);

    /**
     * Return a string that uniquely identifies this finite element extractor.
     */
    std::string
    get_name() const;
  };


  /**
   * Extractor for a symmetric tensor of a rank specified by the template
   * argument. For a second order symmetric tensor, this represents a
   * collection of <code>(dim*dim + dim)/2</code> components of a vector-valued
   * element. The value of <code>dim</code> is defined by the FEValues
   * object the extractor is applied to. The result of applying an object of
   * this type to an FEValues, FEFaceValues or FESubfaceValues object is of
   * type FEValuesViews::SymmetricTensor.
   *
   * The concept of extractors is defined in the documentation of the
   * namespace FEValuesExtractors and in the
   * @ref vector_valued
   * topic.
   *
   * @ingroup feaccess vector_valued
   */
  template <int rank>
  struct SymmetricTensor
  {
    /**
     * The first component of the tensor view.
     */
    unsigned int first_tensor_component;

    /**
     * Default constructor. Initialize the object with an invalid component.
     * This leads to an object that can not be used, but it allows objects of
     * this kind to be put into arrays that require a default constructor upon
     * resizing the array, and then later assigning a suitable object to each
     * element of the array.
     */
    constexpr SymmetricTensor();

    /**
     * Constructor. Take the first component of the selected tensor inside the
     * FEValues object as argument.
     */
    constexpr SymmetricTensor(const unsigned int first_tensor_component);

    /**
     * Return a string that uniquely identifies this finite element extractor.
     */
    std::string
    get_name() const;
  };


  /**
   * Extractor for a general tensor of a given rank specified by
   * the template argument. For a second order tensor, this represents a
   * collection of <code>(dim*dim)</code> components of a vector-valued
   * element. The value of <code>dim</code> is defined by the FEValues object
   * the extractor is applied to. The result of applying an object of this
   * type to an FEValues, FEFaceValues or FESubfaceValues object is of type
   * FEValuesViews::Tensor.
   *
   * The concept of extractors is defined in the documentation of the
   * namespace FEValuesExtractors and in the
   * @ref vector_valued
   * topic.
   *
   * @ingroup feaccess vector_valued
   */
  template <int rank>
  struct Tensor
  {
    /**
     * The first component of the tensor view.
     */
    unsigned int first_tensor_component;

    /**
     * Default constructor. Initialize the object with an invalid component.
     * This leads to an object that can not be used, but it allows objects of
     * this kind to be put into arrays that require a default constructor upon
     * resizing the array, and then later assigning a suitable object to each
     * element of the array.
     */
    constexpr Tensor();

    /**
     * Constructor. Take the first component of the selected tensor inside the
     * FEValues object as argument.
     */
    constexpr Tensor(const unsigned int first_tensor_component);

    /**
     * Return a string that uniquely identifies this finite element extractor.
     */
    std::string
    get_name() const;
  };

  /**
   * Helper struct to associate an extractor to the first FEValuesBase
   * sub-object of a FECouplingValues object.
   */
  template <typename Extractor>
  struct FirstCoupling
  {
    /**
     * Construct a new FirstCoupling object with the given extractor.
     */
    constexpr FirstCoupling(const Extractor &extractor);

    /**
     * The actual extractor object.
     */
    const Extractor extractor;
  };

  /**
   * Helper struct to associate an extractor to the second FEValuesBase
   * sub-object of a FECouplingValues object.
   */
  template <typename Extractor>
  struct SecondCoupling
  {
    /**
     * Construct a new SecondCoupling object with the given extractor.
     */
    constexpr SecondCoupling(const Extractor &extractor);

    /**
     * The actual extractor object.
     */
    const Extractor extractor;
  };
} // namespace FEValuesExtractors


/*-------------- Inline functions: namespace FEValuesExtractors -------------*/

namespace FEValuesExtractors
{
  constexpr inline Scalar::Scalar()
    : component(numbers::invalid_unsigned_int)
  {}



  constexpr inline Scalar::Scalar(const unsigned int component)
    : component(component)
  {}



  constexpr inline Vector::Vector()
    : first_vector_component(numbers::invalid_unsigned_int)
  {}


  constexpr inline Vector::Vector(const unsigned int first_vector_component)
    : first_vector_component(first_vector_component)
  {}


  template <int rank>
  constexpr inline SymmetricTensor<rank>::SymmetricTensor()
    : first_tensor_component(numbers::invalid_unsigned_int)
  {}


  template <int rank>
  constexpr inline SymmetricTensor<rank>::SymmetricTensor(
    const unsigned int first_tensor_component)
    : first_tensor_component(first_tensor_component)
  {}


  template <int rank>
  constexpr inline Tensor<rank>::Tensor()
    : first_tensor_component(numbers::invalid_unsigned_int)
  {}


  template <int rank>
  constexpr inline Tensor<rank>::Tensor(
    const unsigned int first_tensor_component)
    : first_tensor_component(first_tensor_component)
  {}


  template <typename Extractor>
  constexpr inline FirstCoupling<Extractor>::FirstCoupling(
    const Extractor &extractor)
    : extractor(extractor)
  {}


  template <typename Extractor>
  constexpr inline SecondCoupling<Extractor>::SecondCoupling(
    const Extractor &extractor)
    : extractor(extractor)
  {}
} // namespace FEValuesExtractors


DEAL_II_NAMESPACE_CLOSE

#endif
