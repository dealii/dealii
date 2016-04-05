// ---------------------------------------------------------------------
//
// Copyright (C)  2015 by the deal.II authors
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

#ifndef dealii__la_vector_h
#define dealii__la_vector_h


#include <deal.II/base/config.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/vector_space_vector.h>
#include <boost/serialization/array.hpp>
#include <boost/serialization/split_member.hpp>

#include <cstdio>
#include <iostream>
#include <cstring>
#include <vector>


DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  /*! @addtogroup Vectors
   *@{
   */

  /**
   * Numerical vector of data. This class derives from both
   * ::dealii::LinearAlgebra::ReadWriteVector and
   * ::dealii::LinearAlgebra::VectorSpaceVector. As opposed to the array of
   * the C++ standard library, this class implements an element of a vector
   * space suitable for numerical computations.
   *
   * @author Bruno Turcksin, 2015.
   */
  template <typename Number>
  class Vector : public ReadWriteVector<Number>, public VectorSpaceVector<Number>
  {
  public:
    typedef types::global_dof_index size_type;

    /**
     * Constructor. Create a vector of dimension zero.
     */
    Vector();

    /**
     * Copy constructor. Sets the dimension to that of the given vector and
     * copies all elements.
     */
    Vector(const Vector<Number> &V);

    /**
     * Constructor. Set dimension to @p n and initialize all elements with
     * zero.
     *
     * The constructor is made explicit to avoid accident like this:
     * <tt>v=0;</tt>. Presumably, the user wants to set every element of the
     * vector to zero, but instead, what happens is this call:
     * <tt>v=Vector@<Number@>(0);</tt>, i.e. the vector is replaced by one of
     * length zero.
     */
    explicit Vector(const size_type n);

    /**
     * Initialize the vector with a given range of values pointed to by the
     * iterators. This function exists in analogy to the @p std::vector class.
     */
    template <typename InputIterator>
    Vector(const InputIterator first, const InputIterator last);

    /**
     * Destructor, deallocates memory.
     */
    virtual ~Vector();

    /**
     * Copies the data of the input vector @p in_vector.
     */
    Vector<Number> &operator= (const Vector<Number> &in_vector);

    /**
     * Copies the data of the input vector @p in_vector.
     */
    template <typename Number2>
    Vector<Number> &operator= (const Vector<Number2> &in_vector);

    /**
     * Sets all elements of the vector to the scalar @p s. This operation is
     * only allowed if @p s is equal to zero.
     */
    Vector<Number> &operator= (const Number s);

    /**
     * Multiply the entire vector by a fixed factor.
     */
    virtual Vector<Number> &operator*= (const Number factor);

    /**
     * Divide the entire vector by a fixed factor.
     */
    virtual Vector<Number> &operator/= (const Number factor);

    /**
     * Add the vector @p V to the present one.
     */
    virtual Vector<Number> &operator+= (const VectorSpaceVector<Number> &V);

    /**
     * Substract the vector @p V from the present one.
     */
    virtual Vector<Number> &operator-= (const VectorSpaceVector<Number> &V);

    /**
     * Return the scalar product of two vectors.
     */
    virtual Number operator* (const VectorSpaceVector<Number> &V) const;

    /**
     * This function is not implemented and will throw an exception.
     */
    virtual void import(const ReadWriteVector<Number> &V,
                        VectorOperation::values operation,
                        std_cxx11::shared_ptr<const CommunicationPatternBase>
                        communication_pattern =
                          std_cxx11::shared_ptr<const CommunicationPatternBase>());

    /**
     * Add @p a to all components. Note that @p a is a scalar not a vector.
     */
    virtual void add(const Number a);

    /**
     * Simple addition of a multiple of a vector, i.e. <tt>*this += a*V</tt>.
     */
    virtual void add(const Number a, const VectorSpaceVector<Number> &V);

    /**
     * Multiple addition of a multiple of a vector, i.e. <tt>*this +=
     * a*V+b*W</tt>.
     */
    virtual void add(const Number a, const VectorSpaceVector<Number> &V,
                     const Number b, const VectorSpaceVector<Number> &W);

    /**
     * Scaling and simple addition of a multiple of a vector, i.e. <tt>*this =
     * s*(*this)+a*V</tt>.
     */
    virtual void sadd(const Number s, const Number a,
                      const VectorSpaceVector<Number> &V);

    /**
     * Scale each element of this vector by the corresponding element in the
     * argument. This function is mostly meant to simulate multiplication (and
     * immediate re-assignement) by a diagonal scaling matrix.
     */
    virtual void scale(const VectorSpaceVector<Number> &scaling_factors);

    /**
     * Assignement <tt>*this = a*V</tt>.
     */
    virtual void equ(const Number a, const VectorSpaceVector<Number> &V);

    /**
     * Return the l<sub>1</sub> norm of the vector (i.e., the sum of the
     * absolute values of all entries).
     */
    virtual typename VectorSpaceVector<Number>::real_type l1_norm() const;

    /**
     * Return the l<sub>2</sub> norm of the vector (i.e., the square root of
     * the sum of the square of all entries among all processors).
     */
    virtual typename VectorSpaceVector<Number>::real_type l2_norm() const;

    /**
     * Return the maximum norm of the vector (i.e., the maximum absolute value
     * among all entries and among all processors).
     */
    virtual typename VectorSpaceVector<Number>::real_type linfty_norm() const;

    /**
     * Perform a combined operation of a vector addition and a subsequent
     * inner product, returning the value of the inner product. In other
     * words, the result of this function is the same as if the user called
     * @code
     * this->add(a, V);
     * return_value = *this * W;
     * @endcode
     */
    virtual Number add_and_dot(const Number a,
                               const VectorSpaceVector<Number> &V,
                               const VectorSpaceVector<Number> &W);

    /**
     * Return the global size of the vector, equal to the sum of the number of
     * locally owned indices among all processors.
     */
    virtual size_type size() const;

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
    virtual dealii::IndexSet locally_owned_elements() const;

    /**
     * Prints the vector to the output stream @p out.
     */
    virtual void print(std::ostream &out,
                       const unsigned int precision=3,
                       const bool scientific=true,
                       const bool across=true) const;
    /**
     * Write the vector en bloc to a file. This is done in a binary mode, so
     * the output is neither readable by humans nor (probably) by other
     * computers using a different operating system or number format.
     */
    void block_write (std::ostream &out) const;

    /**
     * Read a vector en block from a file. This is done using the inverse
     * operations to the above function, so it is reasonably fast because the
     * bitstream is not interpreted.
     *
     * The vector is resized if necessary.
     *
     * A primitive form of error checking is performed which will recognize
     * the bluntest attempts to interpret some data as a vector stored bitwise
     * to a file, but not more.
     */
    void block_read (std::istream &in);

    /**
     * Returns the memory consumption of this class in bytes.
     */
    virtual std::size_t memory_consumption() const;

    /**
     * Attempt to perform an operation between two incompatible vector types.
     *
     * @ingroup Exceptions
     */
    DeclException0(ExcVectorTypeNotCompatible);

  private:
    /**
     * Compute the L1 norm in a recursive way by dividing the vector on
     * smaller and smaller intervals. This reduces the numerical error on
     * large vector.
     */
    typename VectorSpaceVector<Number>::real_type l1_norm_recursive(unsigned int i,
        unsigned int j) const;

    /**
     * Compute the squared L2 norm in a recursive way by dividing the vector
     * on smaller and smaller intervals. This reduces the numerical error on
     * large vector.
     */
    typename VectorSpaceVector<Number>::real_type l2_norm_squared_recursive(
      unsigned int i,
      unsigned int j) const;

    /**
     * Serialize the data of this object using boost. This function is
     * necessary to use boost::archive::text_iarchive and
     * boost::archive::text_oarchive.
     */
    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version);

    friend class boost::serialization::access;
  };

  /*@}*/
  /*----------------------- Inline functions ----------------------------------*/

  template <typename Number>
  inline
  Vector<Number>::Vector() {}



  template <typename Number>
  inline
  Vector<Number>::Vector(const Vector<Number> &V)
    :
    ReadWriteVector<Number>(V)
  {}



  template <typename Number>
  inline
  Vector<Number>::Vector(const size_type n)
    :
    ReadWriteVector<Number>(n)
  {}



  template <typename Number>
  template <typename InputIterator>
  inline
  Vector<Number>::Vector(const InputIterator first, const InputIterator last)
  {
    this->reinit(complete_index_set(std::distance (first, last)), true);
    std::copy(first, last, this->begin());
  }



  template <typename Number>
  inline
  Vector<Number>::~Vector() {}



  template <typename Number>
  inline
  Vector<Number> &Vector<Number>::operator= (const Vector<Number> &in_vector)
  {
    this->reinit(in_vector.size(),true);
    std::copy(in_vector.begin(), in_vector.end(), this->begin());

    return *this;
  }



  template <typename Number>
  template <typename Number2>
  inline
  Vector<Number> &Vector<Number>::operator= (const Vector<Number2> &in_vector)
  {
    this->reinit(in_vector.size(in_vector.size()),true);
    std::copy(in_vector.begin(), in_vector.end(), this->begin());

    return *this;
  }



  template <typename Number>
  inline
  Vector<Number> &Vector<Number>::operator= (const Number s)
  {
    Assert(s==static_cast<Number>(0), ExcMessage("Only 0 can be assigned to a vector."));
    (void) s;

    std::fill(this->begin(),this->end(),Number());

    return *this;
  }



  template <typename Number>
  inline
  void Vector<Number>::add(const Number a)
  {
    AssertIsFinite(a);
    for (unsigned int i=0; i<this->size(); ++i)
      this->val[i] += a;
  }



  template <typename Number>
  inline
  typename Vector<Number>::size_type Vector<Number>::size() const
  {
    return ReadWriteVector<Number>::size();
  }



  template <typename Number>
  inline
  dealii::IndexSet Vector<Number>::locally_owned_elements() const
  {
    return IndexSet(ReadWriteVector<Number>::get_stored_elements());
  }



  template <typename Number>
  inline
  void Vector<Number>::print(std::ostream &out,
                             const unsigned int precision,
                             const bool scientific,
                             const bool) const
  {
    ReadWriteVector<Number>::print(out, precision, scientific);
  }



  template <typename Number>
  template <typename Archive>
  inline
  void Vector<Number>::serialize(Archive &ar, const unsigned int)
  {
    unsigned int current_size = this->size();
    ar &static_cast<Subscriptor &>(*this);
    ar &this->stored_elements;
    // If necessary, resize the vector during a read operation
    if (this->size()!=current_size)
      this->reinit(this->size());
    ar &boost::serialization::make_array(this->val,this->size());
  }



  template <typename Number>
  inline
  std::size_t Vector<Number>::memory_consumption() const
  {
    return ReadWriteVector<Number>::memory_consumption();
  }
} // end of namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif
