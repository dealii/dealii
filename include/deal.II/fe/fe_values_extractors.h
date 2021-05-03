// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2019 by the deal.II authors
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

#ifndef dealii_fe_values_extractors_h
#define dealii_fe_values_extractors_h


#include <deal.II/base/config.h>

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
 * module for examples how to use the features of this namespace.
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
   * module.
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
    Scalar();

    /**
     * Constructor. Take the selected vector component as argument.
     */
    Scalar(const unsigned int component);

    /**
     * Return a string that uniquely identifies this finite element extractor.
     */
    std::string
    get_name() const;
  };


  /**
   * Extractor for a vector of <code>spacedim</code> components of a vector-
   * valued element. The value of <code>spacedim</code> is defined by the
   * FEValues object the extractor is applied to. The result of applying an
   * object of this type to an FEValues, FEFaceValues or FESubfaceValues
   * object is of type FEValuesViews::Vector.
   *
   * The concept of extractors is defined in the documentation of the
   * namespace FEValuesExtractors and in the
   * @ref vector_valued
   * module.
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
    Vector();

    /**
     * Constructor. Take the first component of the selected vector inside the
     * FEValues object as argument.
     */
    Vector(const unsigned int first_vector_component);

    /**
     * Return a string that uniquely identifies this finite element extractor.
     */
    std::string
    get_name() const;
  };


  /**
   * Extractor for a symmetric tensor of a rank specified by the template
   * argument. For a second order symmetric tensor, this represents a
   * collection of <code>(dim*dim + dim)/2</code> components of a vector-
   * valued element. The value of <code>dim</code> is defined by the FEValues
   * object the extractor is applied to. The result of applying an object of
   * this type to an FEValues, FEFaceValues or FESubfaceValues object is of
   * type FEValuesViews::SymmetricTensor.
   *
   * The concept of extractors is defined in the documentation of the
   * namespace FEValuesExtractors and in the
   * @ref vector_valued
   * module.
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
    SymmetricTensor();

    /**
     * Constructor. Take the first component of the selected tensor inside the
     * FEValues object as argument.
     */
    SymmetricTensor(const unsigned int first_tensor_component);

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
   * module.
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
    Tensor();

    /**
     * Constructor. Take the first component of the selected tensor inside the
     * FEValues object as argument.
     */
    Tensor(const unsigned int first_tensor_component);

    /**
     * Return a string that uniquely identifies this finite element extractor.
     */
    std::string
    get_name() const;
  };
} // namespace FEValuesExtractors


/*-------------- Inline functions: namespace FEValuesExtractors -------------*/

namespace FEValuesExtractors
{
  inline Scalar::Scalar()
    : component(numbers::invalid_unsigned_int)
  {}



  inline Scalar::Scalar(const unsigned int component)
    : component(component)
  {}



  inline Vector::Vector()
    : first_vector_component(numbers::invalid_unsigned_int)
  {}


  inline Vector::Vector(const unsigned int first_vector_component)
    : first_vector_component(first_vector_component)
  {}


  template <int rank>
  inline SymmetricTensor<rank>::SymmetricTensor()
    : first_tensor_component(numbers::invalid_unsigned_int)
  {}


  template <int rank>
  inline SymmetricTensor<rank>::SymmetricTensor(
    const unsigned int first_tensor_component)
    : first_tensor_component(first_tensor_component)
  {}


  template <int rank>
  inline Tensor<rank>::Tensor()
    : first_tensor_component(numbers::invalid_unsigned_int)
  {}


  template <int rank>
  inline Tensor<rank>::Tensor(const unsigned int first_tensor_component)
    : first_tensor_component(first_tensor_component)
  {}
} // namespace FEValuesExtractors



DEAL_II_NAMESPACE_CLOSE

#endif
