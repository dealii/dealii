//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
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
                                        * Constructor. Take the first
                                        * component of the selected tensor
                                        * inside the FEValues object as
                                        * argument.
                                        */
      SymmetricTensor (const unsigned int first_tensor_component);
  };
}


DEAL_II_NAMESPACE_CLOSE

#endif
