// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2016 by the deal.II authors
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

#ifndef dealii__vector_space_vector_h
#define dealii__vector_space_vector_h

#include <deal.II/base/config.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/std_cxx11/shared_ptr.h>
#include <deal.II/lac/vector.h>


DEAL_II_NAMESPACE_OPEN

class CommunicationPatternBase;
class IndexSet;
template <typename Number> class ReadWriteVector;

namespace LinearAlgebra
{
  /*! @addtogroup Vectors
   *@{
   */

  /**
   * VectorSpaceVector is an abstract class that is used to define the
   * interface that vector classes need to implement when they want to
   * implement global operations. This class is complementary of
   * ReadWriteVector which allows the access of individual elements but does
   * not allow global operations.
   *
   * @author Bruno Turcksin, 2015.
   */
  template <typename Number>
  class VectorSpaceVector
  {
  public:
    typedef Number                                            value_type;
    typedef types::global_dof_index                           size_type;
    typedef typename numbers::NumberTraits<Number>::real_type real_type;


    /**
     * Multiply the entire vector by a fixed factor.
     */
    virtual VectorSpaceVector<Number> &operator*= (const Number factor) = 0;

    /**
     * Divide the entire vector by a fixed factor.
     */
    virtual VectorSpaceVector<Number> &operator/= (const Number factor) = 0;

    /**
     * Add the vector @p V to the present one.
     */
    virtual VectorSpaceVector<Number> &operator+= (const VectorSpaceVector<Number> &V) = 0;

    /**
     * Subtract the vector @p V from the present one.
     */
    virtual VectorSpaceVector<Number> &operator-= (const VectorSpaceVector<Number> &V) = 0;

    /**
     * Import all the elements present in the vector's IndexSet from the input
     * vector @p V. VectorOperation::values @p operation is used to decide if
     * the elements in @p V should be added to the current vector or replace the
     * current elements. The last parameter can be used if the same
     * communication pattern is used multiple times. This can be used to improve
     * performance.
     */
    virtual void import(const ReadWriteVector<Number> &V,
                        VectorOperation::values operation,
                        std_cxx11::shared_ptr<const CommunicationPatternBase> communication_pattern =
                          std_cxx11::shared_ptr<const CommunicationPatternBase> ()) = 0;

    /**
     * Return the scalar product of two vectors.
     */
    virtual Number operator* (const VectorSpaceVector<Number> &V) const = 0;

    /**
     * Add @p a to all components. Note that @p a is a scalar not a vector.
     */
    virtual void add(const Number a) = 0;

    /**
     * Simple addition of a multiple of a vector, i.e. <tt>*this += a*V</tt>.
     */
    virtual void add(const Number a, const VectorSpaceVector<Number> &V) = 0;

    /**
     * Multiple addition of scaled vectors, i.e. <tt>*this += a*V+b*W</tt>.
     */
    virtual void add(const Number a, const VectorSpaceVector<Number> &V,
                     const Number b, const VectorSpaceVector<Number> &W) = 0;

    /**
     * Scaling and simple addition of a multiple of a vector, i.e. <tt>*this =
     * s*(*this)+a*V</tt>.
     */
    virtual void sadd(const Number s, const Number a,
                      const VectorSpaceVector<Number> &V) = 0;

    /**
     * Scale each element of this vector by the corresponding element in the
     * argument. This function is mostly meant to simulate multiplication (and
     * immediate re-assignement) by a diagonal scaling matrix.
     */
    virtual void scale(const VectorSpaceVector<Number> &scaling_factors) = 0;

    /**
     * Assignement <tt>*this = a*V</tt>.
     */
    virtual void equ(const Number a, const VectorSpaceVector<Number> &V) = 0;

    /**
     * Return the l<sub>1</sub> norm of the vector (i.e., the sum of the
     * absolute values of all entries among all processors).
     */
    virtual real_type l1_norm() const = 0;

    /**
     * Return the l<sub>2</sub> norm of the vector (i.e., the square root of
     * the sum of the square of all entries among all processors).
     */
    virtual real_type l2_norm() const = 0;

    /**
     * Return the maximum norm of the vector (i.e., the maximum absolute value
     * among all entries and among all processors).
     */
    virtual real_type linfty_norm() const = 0;

    /**
     * Perform a combined operation of a vector addition and a subsequent
     * inner product, returning the value of the inner product. In other
     * words, the result of this function is the same as if the user called
     * @code
     * this->add(a, V);
     * return_value = *this * W;
     * @endcode
     *
     * The reason this function exists is that this operation involves less
     * memory transfer than calling the two functions separately. This method
     * only needs to load three vectors, @p this, @p V, @p W, whereas calling
     * separate methods means to load the calling vector @p this twice. Since
     * most vector operations are memory transfer limited, this reduces the
     * time by 25\% (or 50\% if @p W equals @p this).
     */
    virtual Number add_and_dot(const Number a,
                               const VectorSpaceVector<Number> &V,
                               const VectorSpaceVector<Number> &W) = 0;

    /**
     * Return the global size of the vector, equal to the sum of the number of
     * locally owned indices among all processors.
     */
    virtual size_type size() const = 0;

    /**
     * Return an index set that describes which elements of this vector are
     * owned by the current processor. As a consequence, the index sets
     * returned on different procesors if this is a distributed vector will
     * form disjoint sets that add up to the complete index set. Obviously, if
     * a vector is created on only one processor, then the result would
     * satisfy
     * @code
     *  vec.locally_owned_elements() == complete_index_set(vec.size())
     * @endcode
     */
    virtual dealii::IndexSet locally_owned_elements() const = 0;

    /**
     * Print the vector to the output stream @p out.
     */
    virtual void print(std::ostream &out,
                       const unsigned int precision=3,
                       const bool scientific=true,
                       const bool across=true) const = 0;

    /**
     * Return the memory consumption of this class in bytes.
     */
    virtual std::size_t memory_consumption() const = 0;
  };
  /*@}*/
}

DEAL_II_NAMESPACE_CLOSE

#endif
