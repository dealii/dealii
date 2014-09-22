// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__fe_values_extractors_h
#define __deal2__fe_values_extractors_h


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
 * <code>(dim*dim + dim)/2</code> components
 *
 * See the description of the @ref vector_valued module for examples how to
 * use the features of this namespace.
 *
 * @ingroup feaccess vector_valued
 */
namespace FEValuesExtractors
{
  /**
   * Extractor for a single scalar component
   * of a vector-valued element. The result
   * of applying an object of this type to an
   * FEValues, FEFaceValues or
   * FESubfaceValues object is of type
   * FEValuesViews::Scalar. The concept of
   * extractors is defined in the
   * documentation of the namespace
   * FEValuesExtractors and in the @ref
   * vector_valued module.
   *
   * @ingroup feaccess vector_valued
   */
  struct Scalar
  {
    /**
     * The selected scalar component of the
     * vector.
     */
    unsigned int component;

    /**
     * Default constructor. Initialize the
     * object with an invalid component. This leads to
     * an object that can not be used, but it allows
     * objects of this kind to be put into arrays that
     * require a default constructor upon resizing the
     * array, and then later assigning a suitable
     * object to each element of the array.
     */
    Scalar ();

    /**
     * Constructor. Take the selected
     * vector component as argument.
     */
    Scalar (const unsigned int component);
  };


  /**
   * Extractor for a vector of
   * <code>spacedim</code> components of a
   * vector-valued element. The value of
   * <code>spacedim</code> is defined by the
   * FEValues object the extractor is applied
   * to. The result of applying an object of
   * this type to an FEValues, FEFaceValues
   * or FESubfaceValues object is of type
   * FEValuesViews::Vector.
   *
   * The concept of
   * extractors is defined in the
   * documentation of the namespace
   * FEValuesExtractors and in the @ref
   * vector_valued module.
   *
   * Note that in the current context, a
   * vector is meant in the sense physics
   * uses it: it has <code>spacedim</code>
   * components that behave in specific ways
   * under coordinate system
   * transformations. Examples include
   * velocity or displacement fields. This is
   * opposed to how mathematics uses the word
   * "vector" (and how we use this word in
   * other contexts in the library, for
   * example in the Vector class), where it
   * really stands for a collection of
   * numbers. An example of this latter use
   * of the word could be the set of
   * concentrations of chemical species in a
   * flame; however, these are really just a
   * collection of scalar variables, since
   * they do not change if the coordinate
   * system is rotated, unlike the components
   * of a velocity vector, and consequently,
   * this class should not be used for this
   * context.
   *
   * @ingroup feaccess vector_valued
   */
  struct Vector
  {
    /**
     * The first component of the vector
     * view.
     */
    unsigned int first_vector_component;

    /**
     * Default constructor. Initialize the
     * object with an invalid component. This leads to
     * an object that can not be used, but it allows
     * objects of this kind to be put into arrays that
     * require a default constructor upon resizing the
     * array, and then later assigning a suitable
     * object to each element of the array.
     */
    Vector ();

    /**
     * Constructor. Take the first
     * component of the selected vector
     * inside the FEValues object as
     * argument.
     */
    Vector (const unsigned int first_vector_component);
  };


  /**
   * Extractor for a symmetric tensor of a
   * rank specified by the template
   * argument. For a second order symmetric
   * tensor, this represents a collection of
   * <code>(dim*dim + dim)/2</code>
   * components of a vector-valued
   * element. The value of <code>dim</code>
   * is defined by the FEValues object the
   * extractor is applied to. The result of
   * applying an object of this type to an
   * FEValues, FEFaceValues or
   * FESubfaceValues object is of type
   * FEValuesViews::SymmetricTensor.
   *
   * The concept of
   * extractors is defined in the
   * documentation of the namespace
   * FEValuesExtractors and in the @ref
   * vector_valued module.
   *
   * @ingroup feaccess vector_valued
   *
   * @author Andrew McBride, 2009
   */
  template <int rank>
  struct SymmetricTensor
  {
    /**
     * The first component of the tensor
     * view.
     */
    unsigned int first_tensor_component;

    /**
     * Default constructor. Initialize the
     * object with an invalid component. This leads to
     * an object that can not be used, but it allows
     * objects of this kind to be put into arrays that
     * require a default constructor upon resizing the
     * array, and then later assigning a suitable
     * object to each element of the array.
     */
    SymmetricTensor ();

    /**
     * Constructor. Take the first
     * component of the selected tensor
     * inside the FEValues object as
     * argument.
     */
    SymmetricTensor (const unsigned int first_tensor_component);
  };


  /**
   * Extractor for a (possible non-)symmetric tensor of a
   * rank specified by the template
   * argument. For a second order
   * tensor, this represents a collection of
   * <code>(dim*dim)</code>
   * components of a vector-valued
   * element. The value of <code>dim</code>
   * is defined by the FEValues object the
   * extractor is applied to. The result of
   * applying an object of this type to an
   * FEValues, FEFaceValues or
   * FESubfaceValues object is of type
   * FEValuesViews::Tensor.
   *
   * The concept of
   * extractors is defined in the
   * documentation of the namespace
   * FEValuesExtractors and in the @ref
   * vector_valued module.
   *
   * @ingroup feaccess vector_valued
   *
   * @author Denis Davydov, 2013
   */
  template <int rank>
  struct Tensor
  {
    /**
     * The first component of the tensor
     * view.
     */
    unsigned int first_tensor_component;

    /**
     * Default constructor. Initialize the
     * object with an invalid component. This leads to
     * an object that can not be used, but it allows
     * objects of this kind to be put into arrays that
     * require a default constructor upon resizing the
     * array, and then later assigning a suitable
     * object to each element of the array.
     */
    Tensor ();

    /**
     * Constructor. Take the first
     * component of the selected tensor
     * inside the FEValues object as
     * argument.
     */
    Tensor (const unsigned int first_tensor_component);
  };
}


/*------------------------ Inline functions: namespace FEValuesExtractors --------*/

namespace FEValuesExtractors
{
  inline
  Scalar::Scalar ()
    :
    component (numbers::invalid_unsigned_int)
  {}



  inline
  Scalar::Scalar (const unsigned int component)
    :
    component (component)
  {}



  inline
  Vector::Vector ()
    :
    first_vector_component (numbers::invalid_unsigned_int)
  {}


  inline
  Vector::Vector (const unsigned int first_vector_component)
    :
    first_vector_component (first_vector_component)
  {}


  template <int rank>
  inline
  SymmetricTensor<rank>::SymmetricTensor ()
    :
    first_tensor_component(numbers::invalid_unsigned_int)
  {}


  template <int rank>
  inline
  SymmetricTensor<rank>::SymmetricTensor (const unsigned int first_tensor_component)
    :
    first_tensor_component (first_tensor_component)
  {}


  template <int rank>
  inline
  Tensor<rank>::Tensor ()
    :
    first_tensor_component(numbers::invalid_unsigned_int)
  {}


  template <int rank>
  inline
  Tensor<rank>::Tensor (const unsigned int first_tensor_component)
    :
    first_tensor_component (first_tensor_component)
  {}
}




DEAL_II_NAMESPACE_CLOSE

#endif
