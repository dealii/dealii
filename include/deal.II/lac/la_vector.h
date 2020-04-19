// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II authors
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

#ifndef dealii_la_vector_h
#define dealii_la_vector_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/logstream.h>

#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/vector_operation.h>
#include <deal.II/lac/vector_space_vector.h>
#include <deal.II/lac/vector_type_traits.h>

// boost::serialization::make_array used to be in array.hpp, but was
// moved to a different file in BOOST 1.64
#include <boost/version.hpp>
#if BOOST_VERSION >= 106400
#  include <boost/serialization/array_wrapper.hpp>
#else
#  include <boost/serialization/array.hpp>
#endif
#include <boost/serialization/split_member.hpp>

#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>


DEAL_II_NAMESPACE_OPEN

/**
 * A namespace for vector classes.
 *
 * This namespace contains various classes that provide wrappers to vector
 * classes from different external libraries like Trilinos (EPetra) or PETSc
 * and native implementations like LinearAlgebra::distributed::Vector.
 *
 * The different vector classes are derived from VectorSpaceVector to provide
 * a joint interface for vector space operations, are derived from
 * ReadWriteVector (or ReadWriteVector itself), or both. The separation of
 * vector space operations (like norms or vector additions) through
 * VectorSpaceVector and element access through ReadWriteVector are by design
 * and improve performance.
 */
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
  class Vector : public ReadWriteVector<Number>,
                 public VectorSpaceVector<Number>
  {
  public:
    using size_type  = types::global_dof_index;
    using value_type = typename ReadWriteVector<Number>::value_type;

    /**
     * Constructor. Create a vector of dimension zero.
     */
    Vector() = default;

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
     * Set the global size of the vector to @p size. The stored elements have
     * their index in [0,size).
     *
     * If the flag @p omit_zeroing_entries is set to false, the memory will be
     * initialized with zero, otherwise the memory will be untouched (and the
     * user must make sure to fill it with reasonable data before using it).
     */
    virtual void
    reinit(const size_type size,
           const bool      omit_zeroing_entries = false) override;

    /**
     * Uses the same IndexSet as the one of the input vector @p in_vector and
     * allocates memory for this vector.
     *
     * If the flag @p omit_zeroing_entries is set to false, the memory will be
     * initialized with zero, otherwise the memory will be untouched (and the
     * user must make sure to fill it with reasonable data before using it).
     */
    template <typename Number2>
    void
    reinit(const ReadWriteVector<Number2> &in_vector,
           const bool                      omit_zeroing_entries = false);

    /**
     * Initializes the vector. The indices are specified by @p
     * locally_stored_indices.
     *
     * If the flag @p omit_zeroing_entries is set to false, the memory will be
     * initialized with zero, otherwise the memory will be untouched (and the
     * user must make sure to fill it with reasonable data before using it).
     * locally_stored_indices.
     */
    virtual void
    reinit(const IndexSet &locally_stored_indices,
           const bool      omit_zeroing_entries = false) override;


    /**
     * Change the dimension to that of the vector V. The elements of V are not
     * copied.
     */
    virtual void
    reinit(const VectorSpaceVector<Number> &V,
           const bool omit_zeroing_entries = false) override;

    /**
     * Returns `false` as this is a serial vector.
     *
     * This functionality only needs to be called if using MPI based vectors and
     * exists in other objects for compatibility.
     */
    bool
    has_ghost_elements() const;

    /**
     * Copies the data of the input vector @p in_vector.
     */
    Vector<Number> &
    operator=(const Vector<Number> &in_vector);

    /**
     * Copies the data of the input vector @p in_vector.
     */
    template <typename Number2>
    Vector<Number> &
    operator=(const Vector<Number2> &in_vector);

    /**
     * Sets all elements of the vector to the scalar @p s. This operation is
     * only allowed if @p s is equal to zero.
     */
    virtual Vector<Number> &
    operator=(const Number s) override;

    /**
     * Multiply the entire vector by a fixed factor.
     */
    virtual Vector<Number> &
    operator*=(const Number factor) override;

    /**
     * Divide the entire vector by a fixed factor.
     */
    virtual Vector<Number> &
    operator/=(const Number factor) override;

    /**
     * Add the vector @p V to the present one.
     */
    virtual Vector<Number> &
    operator+=(const VectorSpaceVector<Number> &V) override;

    /**
     * Subtract the vector @p V from the present one.
     */
    virtual Vector<Number> &
    operator-=(const VectorSpaceVector<Number> &V) override;

    /**
     * Return the scalar product of two vectors.
     */
    virtual Number operator*(const VectorSpaceVector<Number> &V) const override;

    /**
     * This function is not implemented and will throw an exception.
     */
    virtual void
    import(
      const ReadWriteVector<Number> &                 V,
      VectorOperation::values                         operation,
      std::shared_ptr<const CommunicationPatternBase> communication_pattern =
        std::shared_ptr<const CommunicationPatternBase>()) override;

    /**
     * Add @p a to all components. Note that @p a is a scalar not a vector.
     */
    virtual void
    add(const Number a) override;

    /**
     * Simple addition of a multiple of a vector, i.e. <tt>*this += a*V</tt>.
     */
    virtual void
    add(const Number a, const VectorSpaceVector<Number> &V) override;

    /**
     * Multiple addition of a multiple of a vector, i.e. <tt>*this +=
     * a*V+b*W</tt>.
     */
    virtual void
    add(const Number                     a,
        const VectorSpaceVector<Number> &V,
        const Number                     b,
        const VectorSpaceVector<Number> &W) override;

    /**
     * Scaling and simple addition of a multiple of a vector, i.e. <tt>*this =
     * s*(*this)+a*V</tt>.
     */
    virtual void
    sadd(const Number                     s,
         const Number                     a,
         const VectorSpaceVector<Number> &V) override;

    /**
     * Scale each element of this vector by the corresponding element in the
     * argument. This function is mostly meant to simulate multiplication (and
     * immediate re-assignment) by a diagonal scaling matrix.
     */
    virtual void
    scale(const VectorSpaceVector<Number> &scaling_factors) override;

    /**
     * Assignment <tt>*this = a*V</tt>.
     */
    virtual void
    equ(const Number a, const VectorSpaceVector<Number> &V) override;

    /**
     * Return whether the vector contains only elements with value zero.
     */
    virtual bool
    all_zero() const override;

    /**
     * Return the mean value of all the entries of this vector.
     */
    virtual value_type
    mean_value() const override;

    /**
     * Return the l<sub>1</sub> norm of the vector (i.e., the sum of the
     * absolute values of all entries).
     */
    virtual typename VectorSpaceVector<Number>::real_type
    l1_norm() const override;

    /**
     * Return the l<sub>2</sub> norm of the vector (i.e., the square root of
     * the sum of the square of all entries among all processors).
     */
    virtual typename VectorSpaceVector<Number>::real_type
    l2_norm() const override;

    /**
     * Return the maximum norm of the vector (i.e., the maximum absolute value
     * among all entries and among all processors).
     */
    virtual typename VectorSpaceVector<Number>::real_type
    linfty_norm() const override;

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
     * most vector operations are memory transfer limited, this reduces the time
     * by 25\% (or 50\% if @p W equals @p this).
     *
     * For complex-valued vectors, the scalar product in the second step is
     * implemented as
     * $\left<v,w\right>=\sum_i v_i \bar{w_i}$.
     */
    virtual Number
    add_and_dot(const Number                     a,
                const VectorSpaceVector<Number> &V,
                const VectorSpaceVector<Number> &W) override;

    /**
     * Return the global size of the vector, equal to the sum of the number of
     * locally owned indices among all processors.
     */
    virtual size_type
    size() const override;

    /**
     * Return an index set that describes which elements of this vector are
     * owned by the current processor. As a consequence, the index sets
     * returned on different processors if this is a distributed vector will
     * form disjoint sets that add up to the complete index set. Obviously, if
     * a vector is created on only one processor, then the result would
     * satisfy
     * @code
     *  vec.locally_owned_elements() == complete_index_set(vec.size())
     * @endcode
     */
    virtual dealii::IndexSet
    locally_owned_elements() const override;

    /**
     * Print the vector to the output stream @p out.
     */
    virtual void
    print(std::ostream &     out,
          const unsigned int precision  = 3,
          const bool         scientific = true,
          const bool         across     = true) const override;

    /**
     * Print the vector to the output stream @p out in a format that can be
     * read by numpy::readtxt(). Note that the IndexSet is not printed but only
     * the values stored in the Vector. To load the vector in python just do
     * <code>
     * vector = numpy.loadtxt('my_vector.txt')
     * </code>
     */
    void
    print_as_numpy_array(std::ostream &     out,
                         const unsigned int precision = 9) const;

    /**
     * Write the vector en bloc to a file. This is done in a binary mode, so
     * the output is neither readable by humans nor (probably) by other
     * computers using a different operating system or number format.
     */
    void
    block_write(std::ostream &out) const;

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
    void
    block_read(std::istream &in);

    /**
     * Return the memory consumption of this class in bytes.
     */
    virtual std::size_t
    memory_consumption() const override;

    /**
     * Attempt to perform an operation between two incompatible vector types.
     *
     * @ingroup Exceptions
     */
    DeclException0(ExcVectorTypeNotCompatible);

  private:
    /**
     * Serialize the data of this object using boost. This function is
     * necessary to use boost::archive::text_iarchive and
     * boost::archive::text_oarchive.
     */
    template <typename Archive>
    void
    serialize(Archive &ar, const unsigned int version);

    friend class boost::serialization::access;

    // Make all other ReadWriteVector types friends.
    template <typename Number2>
    friend class Vector;
  };

  /*@}*/
  /*--------------------------- Inline functions ----------------------------*/

  template <typename Number>
  inline Vector<Number>::Vector(const Vector<Number> &V)
    : ReadWriteVector<Number>(V)
  {}



  template <typename Number>
  inline Vector<Number>::Vector(const size_type n)
    : ReadWriteVector<Number>(n)
  {}



  template <typename Number>
  template <typename InputIterator>
  inline Vector<Number>::Vector(const InputIterator first,
                                const InputIterator last)
  {
    this->reinit(complete_index_set(std::distance(first, last)), true);
    std::copy(first, last, this->begin());
  }



  template <typename Number>
  inline typename Vector<Number>::size_type
  Vector<Number>::size() const
  {
    return ReadWriteVector<Number>::size();
  }



  template <typename Number>
  inline dealii::IndexSet
  Vector<Number>::locally_owned_elements() const
  {
    return IndexSet(ReadWriteVector<Number>::get_stored_elements());
  }



  template <typename Number>
  inline void
  Vector<Number>::print(std::ostream &     out,
                        const unsigned int precision,
                        const bool         scientific,
                        const bool) const
  {
    ReadWriteVector<Number>::print(out, precision, scientific);
  }



  template <typename Number>
  template <typename Archive>
  inline void
  Vector<Number>::serialize(Archive &ar, const unsigned int)
  {
    size_type current_size = this->size();
    ar &static_cast<Subscriptor &>(*this);
    ar & this->stored_elements;
    // If necessary, resize the vector during a read operation
    if (this->size() != current_size)
      this->reinit(this->size());
    ar &boost::serialization::make_array(this->values.get(), this->size());
  }



  template <typename Number>
  inline std::size_t
  Vector<Number>::memory_consumption() const
  {
    return ReadWriteVector<Number>::memory_consumption();
  }
} // end of namespace LinearAlgebra


/**
 * Declare dealii::LinearAlgebra::Vector as serial vector.
 *
 * @author Uwe Koecher, 2017
 */
template <typename Number>
struct is_serial_vector<LinearAlgebra::Vector<Number>> : std::true_type
{};

#ifndef DOXYGEN
/*----------------------- Inline functions ----------------------------------*/

namespace LinearAlgebra
{
  template <typename Number>
  inline bool
  Vector<Number>::has_ghost_elements() const
  {
    return false;
  }
} // namespace LinearAlgebra

#endif


DEAL_II_NAMESPACE_CLOSE

#ifdef DEAL_II_MSVC
#  include <deal.II/lac/la_vector.templates.h>
#endif

#endif
