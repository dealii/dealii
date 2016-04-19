// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2016 by the deal.II authors
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

#ifndef dealii__vector_h
#define dealii__vector_h


#include <deal.II/base/config.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/index_set.h>
#include <boost/serialization/array.hpp>
#include <boost/serialization/split_member.hpp>

#include <cstdio>
#include <iostream>
#include <cstring>
#include <vector>

DEAL_II_NAMESPACE_OPEN


#ifdef DEAL_II_WITH_PETSC
namespace PETScWrappers
{
  class Vector;
  namespace MPI
  {
    class Vector;
  }
}
#endif

#ifdef DEAL_II_WITH_TRILINOS
namespace TrilinosWrappers
{
  namespace MPI
  {
    class Vector;
  }
  class Vector;
}
#endif

template<typename number> class LAPACKFullMatrix;

template <typename> class BlockVector;

template <typename> class VectorView;

namespace parallel
{
  namespace internal
  {
    class TBBPartitioner;
  }
}



/*! @addtogroup Vectors
 *@{
 */

/**
 * This enum keeps track of the current operation in parallel linear algebra
 * objects like Vectors and Matrices.
 *
 * It is used in the various compress() functions. They also exist in serial
 * codes for compatibility and are empty there.
 *
 * See
 * @ref GlossCompress "Compressing distributed objects"
 * for more information.
 */
struct VectorOperation
{
  enum values { unknown, insert, add };
};


/**
 * Numerical vector of data.  For this class there are different types of
 * functions available. The first type of function initializes the vector,
 * changes its size, or computes the norm of the vector in order to measure
 * its length in a suitable norm. The second type helps us to manipulate the
 * components of the vector. The third type defines the algebraic operations
 * for vectors, while the last type defines a few input and output functions.
 * As opposed to the array of the C++ standard library called @p vector (with
 * a lowercase "v"), this class implements an element of a vector space
 * suitable for numerical computations.
 *
 * @note Instantiations for this template are provided for <tt>@<float@>,
 * @<double@>, @<long double@>, @<std::complex@<float@>@>,
 * @<std::complex@<double@>@>, @<std::complex@<long double@>@></tt>; others
 * can be generated in application programs (see the section on
 * @ref Instantiations
 * in the manual).
 *
 * @author Guido Kanschat, Franz-Theo Suttmeier, Wolfgang Bangerth
 */
template <typename Number>
class Vector : public Subscriptor
{
public:
  /**
   * Declare standard types used in all containers. These types parallel those
   * in the <tt>C++</tt> standard libraries <tt>vector<...></tt> class.
   */
  typedef Number                                            value_type;
  typedef value_type                                       *pointer;
  typedef const value_type                                 *const_pointer;
  typedef value_type                                       *iterator;
  typedef const value_type                                 *const_iterator;
  typedef value_type                                       &reference;
  typedef const value_type                                 &const_reference;
  typedef types::global_dof_index                           size_type;

  /**
   * Declare a type that has holds real-valued numbers with the same precision
   * as the template argument to this class. If the template argument of this
   * class is a real data type, then real_type equals the template argument.
   * If the template argument is a std::complex type then real_type equals the
   * type underlying the complex numbers.
   *
   * This typedef is used to represent the return type of norms.
   */
  typedef typename numbers::NumberTraits<Number>::real_type real_type;

  /**
   * A variable that indicates whether this vector supports distributed data
   * storage. If true, then this vector also needs an appropriate compress()
   * function that allows communicating recent set or add operations to
   * individual elements to be communicated to other processors.
   *
   * For the current class, the variable equals false, since it does not
   * support parallel data storage.
   */
  static const bool supports_distributed_data = false;

public:

  /**
   * @name Basic object handling
   */
  //@{
  /**
   * Constructor. Create a vector of dimension zero.
   */
  Vector ();

  /**
   * Copy constructor. Sets the dimension to that of the given vector, and
   * copies all elements.
   *
   * We would like to make this constructor explicit, but standard containers
   * insist on using it implicitly.
   */
  Vector (const Vector<Number> &v);

#ifdef DEAL_II_WITH_CXX11
  /**
   * Move constructor. Creates a new vector by stealing the internal data of
   * the vector @p v.
   *
   * @note This constructor is only available if deal.II is configured with
   * C++11 support.
   */
  Vector (Vector<Number> &&v);
#endif


#ifndef DEAL_II_EXPLICIT_CONSTRUCTOR_BUG
  /**
   * Copy constructor taking a vector of another data type. This will fail if
   * there is no conversion path from @p OtherNumber to @p Number. Note that
   * you may lose accuracy when copying to a vector with data elements with
   * less accuracy.
   *
   * Older versions of gcc did not honor the @p explicit keyword on template
   * constructors. In such cases, it is easy to accidentally write code that
   * can be very inefficient, since the compiler starts performing hidden
   * conversions. To avoid this, this function is disabled if we have detected
   * a broken compiler during configuration.
   */
  template <typename OtherNumber>
  explicit
  Vector (const Vector<OtherNumber> &v);
#endif

#ifdef DEAL_II_WITH_PETSC
  /**
   * Another copy constructor: copy the values from a sequential PETSc wrapper
   * vector class. This copy constructor is only available if PETSc was
   * detected during configuration time.
   */
  Vector (const PETScWrappers::Vector &v);

  /**
   * Another copy constructor: copy the values from a parallel PETSc wrapper
   * vector class. This copy constructor is only available if PETSc was
   * detected during configuration time.
   *
   * Note that due to the communication model used in MPI, this operation can
   * only succeed if all processes do it at the same time. I.e., it is not
   * possible for only one process to obtain a copy of a parallel vector while
   * the other jobs do something else.
   */
  Vector (const PETScWrappers::MPI::Vector &v);
#endif

#ifdef DEAL_II_WITH_TRILINOS
  /**
   * Another copy constructor: copy the values from a Trilinos wrapper vector.
   * This copy constructor is only available if Trilinos was detected during
   * configuration time.
   *
   * Note that due to the communication model used in MPI, this operation can
   * only succeed if all processes do it at the same time. This means that it
   * is not possible for only one process to obtain a copy of a parallel
   * vector while the other jobs do something else. This call will rather
   * result in a copy of the vector on all processors.
   */
  Vector (const TrilinosWrappers::MPI::Vector &v);

  /**
   * Another copy constructor: copy the values from a localized Trilinos
   * wrapper vector. This copy constructor is only available if Trilinos was
   * detected during configuration time.
   */
  Vector (const TrilinosWrappers::Vector &v);
#endif

  /**
   * Constructor. Set dimension to @p n and initialize all elements with zero.
   *
   * The constructor is made explicit to avoid accidents like this:
   * <tt>v=0;</tt>. Presumably, the user wants to set every element of the
   * vector to zero, but instead, what happens is this call:
   * <tt>v=Vector@<number@>(0);</tt>, i.e. the vector is replaced by one of
   * length zero.
   */
  explicit Vector (const size_type n);

  /**
   * Initialize the vector with a given range of values pointed to by the
   * iterators. This function is there in analogy to the @p std::vector class.
   */
  template <typename InputIterator>
  Vector (const InputIterator first,
          const InputIterator last);

  /**
   * Destructor, deallocates memory. Made virtual to allow for derived classes
   * to behave properly.
   */
  virtual ~Vector ();

  /**
   * This function does nothing but is there for compatibility with the @p
   * PETScWrappers::Vector class.
   *
   * For the PETSc vector wrapper class, this function compresses the
   * underlying representation of the PETSc object, i.e. flushes the buffers
   * of the vector object if it has any. This function is necessary after
   * writing into a vector element-by-element and before anything else can be
   * done on it.
   *
   * However, for the implementation of this class, it is immaterial and thus
   * an empty function.
   */
  void compress (::dealii::VectorOperation::values operation
                 =::dealii::VectorOperation::unknown) const;

  /**
   * Change the dimension of the vector to @p N. The reserved memory for this
   * vector remains unchanged if possible, to make things faster; this may
   * waste some memory, so keep this in mind.  However, if <tt>N==0</tt> all
   * memory is freed, i.e. if you want to resize the vector and release the
   * memory not needed, you have to first call <tt>reinit(0)</tt> and then
   * <tt>reinit(N)</tt>. This cited behaviour is analogous to that of the
   * standard library containers.
   *
   * If @p omit_zeroing_entries is false, the vector is filled by zeros.
   * Otherwise, the elements are left an unspecified state.
   *
   * This function is virtual in order to allow for derived classes to handle
   * memory separately.
   */
  virtual void reinit (const size_type N,
                       const bool      omit_zeroing_entries=false);

  /**
   * Change the dimension to that of the vector @p V. The same applies as for
   * the other @p reinit function.
   *
   * The elements of @p V are not copied, i.e.  this function is the same as
   * calling <tt>reinit (V.size(), omit_zeroing_entries)</tt>.
   */
  template <typename Number2>
  void reinit (const Vector<Number2> &V,
               const bool            omit_zeroing_entries=false);

  /**
   * Swap the contents of this vector and the other vector @p v. One could do
   * this operation with a temporary variable and copying over the data
   * elements, but this function is significantly more efficient since it only
   * swaps the pointers to the data of the two vectors and therefore does not
   * need to allocate temporary storage and move data around.
   *
   * This function is analog to the the @p swap function of all C++ standard
   * containers. Also, there is a global function <tt>swap(u,v)</tt> that
   * simply calls <tt>u.swap(v)</tt>, again in analogy to standard functions.
   *
   * This function is virtual in order to allow for derived classes to handle
   * memory separately.
   */
  virtual void swap (Vector<Number> &v);

  /**
   * Set all components of the vector to the given number @p s. Simply pass
   * this down to the individual block objects, but we still need to declare
   * this function to make the example given in the discussion about making
   * the constructor explicit work.
   *
   * Since the semantics of assigning a scalar to a vector are not immediately
   * clear, this operator should really only be used if you want to set the
   * entire vector to zero. This allows the intuitive notation <tt>v=0</tt>.
   * Assigning other values is deprecated and may be disallowed in the future.
   *
   * @dealiiOperationIsMultithreaded
   */
  Vector<Number> &operator= (const Number s);

  /**
   * Copy the given vector. Resize the present vector if necessary.
   *
   * @dealiiOperationIsMultithreaded
   */
  Vector<Number> &operator= (const Vector<Number> &v);

#ifdef DEAL_II_WITH_CXX11
  /**
   * Move the given vector. This operator replaces the present vector with
   * the internal data of the vector @p v and resets @p v to the state it would
   * have after being newly default-constructed.
   *
   * @note This operator is only available if deal.II is configured with C++11
   * support.
   */
  Vector<Number> &operator= (Vector<Number> &&v);
#endif

  /**
   * Copy the given vector. Resize the present vector if necessary.
   *
   * @dealiiOperationIsMultithreaded
   */
  template <typename Number2>
  Vector<Number> &operator= (const Vector<Number2> &v);

  /**
   * Copy operator for assigning a block vector to a regular vector.
   */
  Vector<Number> &operator= (const BlockVector<Number> &v);

#ifdef DEAL_II_WITH_PETSC
  /**
   * Another copy operator: copy the values from a sequential PETSc wrapper
   * vector class. This operator is only available if PETSc was detected
   * during configuration time.
   */
  Vector<Number> &
  operator= (const PETScWrappers::Vector &v);

  /**
   * Another copy operator: copy the values from a parallel PETSc wrapper
   * vector class. This operator is only available if PETSc was detected
   * during configuration time.
   *
   * Note that due to the communication model used in MPI, this operation can
   * only succeed if all processes do it at the same time. I.e., it is not
   * possible for only one process to obtain a copy of a parallel vector while
   * the other jobs do something else.
   */
  Vector<Number> &
  operator= (const PETScWrappers::MPI::Vector &v);
#endif


#ifdef DEAL_II_WITH_TRILINOS
  /**
   * Another copy operator: copy the values from a (sequential or parallel,
   * depending on the underlying compiler) Trilinos wrapper vector class. This
   * operator is only available if Trilinos was detected during configuration
   * time.
   *
   * Note that due to the communication model used in MPI, this operation can
   * only succeed if all processes do it at the same time. I.e., it is not
   * possible for only one process to obtain a copy of a parallel vector while
   * the other jobs do something else.
   */
  Vector<Number> &
  operator= (const TrilinosWrappers::MPI::Vector &v);

  /**
   * Another copy operator: copy the values from a sequential Trilinos wrapper
   * vector class. This operator is only available if Trilinos was detected
   * during configuration time.
   */
  Vector<Number> &
  operator= (const TrilinosWrappers::Vector &v);
#endif

  /**
   * Test for equality. This function assumes that the present vector and the
   * one to compare with have the same size already, since comparing vectors
   * of different sizes makes not much sense anyway.
   */
  template <typename Number2>
  bool operator== (const Vector<Number2> &v) const;

  /**
   * Test for inequality. This function assumes that the present vector and
   * the one to compare with have the same size already, since comparing
   * vectors of different sizes makes not much sense anyway.
   */
  template <typename Number2>
  bool operator != (const Vector<Number2> &v) const;

  //@}


  /**
   * @name Scalar products, norms and related operations
   */
  //@{

  /**
   * Return the scalar product of two vectors.  The return type is the
   * underlying type of @p this vector, so the return type and the accuracy
   * with which it the result is computed depend on the order of the arguments
   * of this vector.
   *
   * For complex vectors, the scalar product is implemented as
   * $\left<v,w\right>=\sum_i v_i \bar{w_i}$.
   *
   * @dealiiOperationIsMultithreaded The algorithm uses pairwise summation
   * with the same order of summation in every run, which gives fully
   * repeatable results from one run to another.
   */
  template <typename Number2>
  Number operator * (const Vector<Number2> &V) const;

  /**
   * Return square of the $l_2$-norm.
   *
   * @dealiiOperationIsMultithreaded The algorithm uses pairwise summation
   * with the same order of summation in every run, which gives fully
   * repeatable results from one run to another.
   */
  real_type norm_sqr () const;

  /**
   * Mean value of the elements of this vector.
   *
   * @dealiiOperationIsMultithreaded The algorithm uses pairwise summation
   * with the same order of summation in every run, which gives fully
   * repeatable results from one run to another.
   */
  Number mean_value () const;

  /**
   * $l_1$-norm of the vector. The sum of the absolute values.
   *
   * @dealiiOperationIsMultithreaded The algorithm uses pairwise summation
   * with the same order of summation in every run, which gives fully
   * repeatable results from one run to another.
   */
  real_type l1_norm () const;

  /**
   * $l_2$-norm of the vector. The square root of the sum of the squares of
   * the elements.
   *
   * @dealiiOperationIsMultithreaded The algorithm uses pairwise summation
   * with the same order of summation in every run, which gives fully
   * repeatable results from one run to another.
   */
  real_type l2_norm () const;

  /**
   * $l_p$-norm of the vector. The pth root of the sum of the pth powers of
   * the absolute values of the elements.
   *
   * @dealiiOperationIsMultithreaded The algorithm uses pairwise summation
   * with the same order of summation in every run, which gives fully
   * repeatable results from one run to another.
   */
  real_type lp_norm (const real_type p) const;

  /**
   * Maximum absolute value of the elements.
   */
  real_type linfty_norm () const;

  /**
   * Performs a combined operation of a vector addition and a subsequent inner
   * product, returning the value of the inner product. In other words, the
   * result of this function is the same as if the user called
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
   * @dealiiOperationIsMultithreaded The algorithm uses pairwise summation
   * with the same order of summation in every run, which gives fully
   * repeatable results from one run to another.
   */
  Number add_and_dot (const Number          a,
                      const Vector<Number> &V,
                      const Vector<Number> &W);

  //@}


  /**
   * @name Data access
   */
  //@{

  /**
   * Make the @p Vector class a bit like the <tt>vector<></tt> class of the
   * C++ standard library by returning iterators to the start and end of the
   * elements of this vector.
   */
  iterator begin ();

  /**
   * Return constant iterator to the start of the vectors.
   */
  const_iterator begin () const;

  /**
   * Return an iterator pointing to the element past the end of the array.
   */
  iterator end ();

  /**
   * Return a constant iterator pointing to the element past the end of the
   * array.
   */
  const_iterator end () const;

  /**
   * Access the value of the @p ith component.
   */
  Number operator() (const size_type i) const;

  /**
   * Access the @p ith component as a writeable reference.
   */
  Number &operator() (const size_type i);

  /**
   * Access the value of the @p ith component.
   *
   * Exactly the same as operator().
   */
  Number operator[] (const size_type i) const;

  /**
   * Access the @p ith component as a writeable reference.
   *
   * Exactly the same as operator().
   */
  Number &operator[] (const size_type i);

  /**
   * A collective get operation: instead of getting individual elements of a
   * vector, this function allows to get a whole set of elements at once. The
   * indices of the elements to be read are stated in the first argument, the
   * corresponding values are returned in the second.
   */
  template <typename OtherNumber>
  void extract_subvector_to (const std::vector<size_type> &indices,
                             std::vector<OtherNumber> &values) const;

  /**
   * Just as the above, but with pointers. Useful in minimizing copying of
   * data around.
   */
  template <typename ForwardIterator, typename OutputIterator>
  void extract_subvector_to (ForwardIterator       indices_begin,
                             const ForwardIterator indices_end,
                             OutputIterator        values_begin) const;
  //@}


  /**
   * @name Modification of vectors
   */
  //@{

  /**
   * Add the given vector to the present one.
   *
   * @dealiiOperationIsMultithreaded
   */
  Vector<Number> &operator += (const Vector<Number> &V);

  /**
   * Subtract the given vector from the present one.
   *
   * @dealiiOperationIsMultithreaded
   */
  Vector<Number> &operator -= (const Vector<Number> &V);

  /**
   * A collective add operation: This function adds a whole set of values
   * stored in @p values to the vector components specified by @p indices.
   */
  template <typename OtherNumber>
  void add (const std::vector<size_type>   &indices,
            const std::vector<OtherNumber>  &values);

  /**
   * This is a second collective add operation. As a difference, this function
   * takes a deal.II vector of values.
   */
  template <typename OtherNumber>
  void add (const std::vector<size_type> &indices,
            const Vector<OtherNumber>    &values);

  /**
   * Take an address where <tt>n_elements</tt> are stored contiguously and add
   * them into the vector. Handles all cases which are not covered by the
   * other two <tt>add()</tt> functions above.
   */
  template <typename OtherNumber>
  void add (const size_type    n_elements,
            const size_type   *indices,
            const OtherNumber  *values);

  /**
   * Addition of @p s to all components. Note that @p s is a scalar and not a
   * vector.
   *
   * @dealiiOperationIsMultithreaded
   */
  void add (const Number s);

  /**
   * Simple vector addition, equal to the <tt>operator +=</tt>.
   *
   * @deprecated Use the <tt>operator +=</tt> instead.
   *
   * @dealiiOperationIsMultithreaded
   */
  void add (const Vector<Number> &V) DEAL_II_DEPRECATED;


  /**
   * Multiple addition of scaled vectors, i.e. <tt>*this += a*V+b*W</tt>.
   *
   * @dealiiOperationIsMultithreaded
   */
  void add (const Number a, const Vector<Number> &V,
            const Number b, const Vector<Number> &W);

  /**
   * Simple addition of a multiple of a vector, i.e. <tt>*this += a*V</tt>.
   *
   * @dealiiOperationIsMultithreaded
   */
  void add (const Number a, const Vector<Number> &V);

  /**
   * Scaling and simple vector addition, i.e.  <tt>*this = s*(*this)+V</tt>.
   *
   * @dealiiOperationIsMultithreaded
   */
  void sadd (const Number          s,
             const Vector<Number> &V);

  /**
   * Scaling and simple addition, i.e.  <tt>*this = s*(*this)+a*V</tt>.
   *
   * @dealiiOperationIsMultithreaded
   */
  void sadd (const Number          s,
             const Number          a,
             const Vector<Number> &V);

  /**
   * Scaling and multiple addition.
   *
   * This function is deprecated.
   *
   * @dealiiOperationIsMultithreaded
   */
  void sadd (const Number          s,
             const Number          a,
             const Vector<Number> &V,
             const Number          b,
             const Vector<Number> &W) DEAL_II_DEPRECATED;

  /**
   * Scaling and multiple addition.  <tt>*this = s*(*this)+a*V + b*W +
   * c*X</tt>.
   *
   * This function is deprecated.
   *
   * @dealiiOperationIsMultithreaded
   */
  void sadd (const Number          s,
             const Number          a,
             const Vector<Number> &V,
             const Number          b,
             const Vector<Number> &W,
             const Number          c,
             const Vector<Number> &X) DEAL_II_DEPRECATED;

  /**
   * Scale each element of the vector by a constant value.
   *
   * @dealiiOperationIsMultithreaded
   */
  Vector<Number> &operator *= (const Number factor);

  /**
   * Scale each element of the vector by the inverse of the given value.
   *
   * @dealiiOperationIsMultithreaded
   */
  Vector<Number> &operator /= (const Number factor);

  /**
   * Scale each element of this vector by the corresponding element in the
   * argument. This function is mostly meant to simulate multiplication (and
   * immediate re-assignment) by a diagonal scaling matrix.
   *
   * @dealiiOperationIsMultithreaded
   */
  void scale (const Vector<Number> &scaling_factors);

  /**
   * Scale each element of this vector by the corresponding element in the
   * argument. This function is mostly meant to simulate multiplication (and
   * immediate re-assignment) by a diagonal scaling matrix.
   */
  template <typename Number2>
  void scale (const Vector<Number2> &scaling_factors);

  /**
   * Assignment <tt>*this = a*u</tt>.
   *
   * @dealiiOperationIsMultithreaded
   */
  void equ (const Number a, const Vector<Number> &u);

  /**
   * Assignment <tt>*this = a*u</tt>.
   */
  template <typename Number2>
  void equ (const Number a, const Vector<Number2> &u);

  /**
   * Assignment <tt>*this = a*u + b*v</tt>.
   *
   * This function is deprecated.
   *
   * @dealiiOperationIsMultithreaded
   */
  void equ (const Number a, const Vector<Number> &u,
            const Number b, const Vector<Number> &v) DEAL_II_DEPRECATED;

  /**
   * Assignment <tt>*this = a*u + b*v + b*w</tt>.
   *
   * This function is deprecated.
   *
   * @dealiiOperationIsMultithreaded
   */
  void equ (const Number a, const Vector<Number> &u,
            const Number b, const Vector<Number> &v,
            const Number c, const Vector<Number> &w) DEAL_II_DEPRECATED;

  /**
   * Compute the elementwise ratio of the two given vectors, that is let
   * <tt>this[i] = a[i]/b[i]</tt>. This is useful for example if you want to
   * compute the cellwise ratio of true to estimated error.
   *
   * This vector is appropriately scaled to hold the result.
   *
   * If any of the <tt>b[i]</tt> is zero, the result is undefined. No attempt
   * is made to catch such situations.
   *
   * @dealiiOperationIsMultithreaded
   */
  void ratio (const Vector<Number> &a,
              const Vector<Number> &b) DEAL_II_DEPRECATED;

  /**
   * This function does nothing but is there for compatibility with the @p
   * PETScWrappers::Vector class.
   *
   * For the PETSc vector wrapper class, this function updates the ghost
   * values of the PETSc vector. This is necessary after any modification
   * before reading ghost values.
   *
   * However, for the implementation of this class, it is immaterial and thus
   * an empty function.
   */
  void update_ghost_values () const;
  //@}


  /**
   * @name Input and output
   */
  //@{
  /**
   * Output of vector in user-defined format. For complex-valued vectors, the
   * format should include specifiers for both the real and imaginary parts.
   *
   * This function is deprecated.
   */
  void print (const char *format = 0) const DEAL_II_DEPRECATED;

  /**
   * Print to a stream. @p precision denotes the desired precision with which
   * values shall be printed, @p scientific whether scientific notation shall
   * be used. If @p across is @p true then the vector is printed in a line,
   * while if @p false then the elements are printed on a separate line each.
   */
  void print (std::ostream &out,
              const unsigned int precision  = 3,
              const bool scientific = true,
              const bool across     = true) const;

  /**
   * Print to a LogStream. <tt>width</tt> is used as argument to the std::setw
   * manipulator, if printing across.  If @p across is @p true then the vector
   * is printed in a line, while if @p false then the elements are printed on
   * a separate line each.
   *
   * This function is deprecated.
   */
  void print (LogStream &out,
              const unsigned int width = 6,
              const bool across = true) const DEAL_II_DEPRECATED;

  /**
   * Write the vector en bloc to a file. This is done in a binary mode, so the
   * output is neither readable by humans nor (probably) by other computers
   * using a different operating system or number format.
   */
  void block_write (std::ostream &out) const;

  /**
   * Read a vector en block from a file. This is done using the inverse
   * operations to the above function, so it is reasonably fast because the
   * bitstream is not interpreted.
   *
   * The vector is resized if necessary.
   *
   * A primitive form of error checking is performed which will recognize the
   * bluntest attempts to interpret some data as a vector stored bitwise to a
   * file, but not more.
   */
  void block_read (std::istream &in);

  /**
   * Write the data of this object to a stream for the purpose of
   * serialization.
   */
  template <class Archive>
  void save (Archive &ar, const unsigned int version) const;

  /**
   * Read the data of this object from a stream for the purpose of
   * serialization.
   */
  template <class Archive>
  void load (Archive &ar, const unsigned int version);

  BOOST_SERIALIZATION_SPLIT_MEMBER()

  /**
   * @}
   */

  /**
   * @name Information about the object
   */
  //@{

  /**
   * Returns true if the given global index is in the local range of this
   * processor.  Since this is not a distributed vector the method always
   * returns true.
   */
  bool in_local_range (const size_type global_index) const;

  /**
   * Return an index set that describes which elements of this vector are
   * owned by the current processor. Note that this index set does not include
   * elements this vector may store locally as ghost elements but that are in
   * fact owned by another processor. As a consequence, the index sets
   * returned on different processors if this is a distributed vector will
   * form disjoint sets that add up to the complete index set. Obviously, if a
   * vector is created on only one processor, then the result would satisfy
   * @code
   *   vec.locally_owned_elements() == complete_index_set (vec.size())
   * @endcode
   *
   * Since the current data type does not support parallel data storage across
   * different processors, the returned index set is the complete index set.
   */
  IndexSet locally_owned_elements () const;

  /**
   * Return dimension of the vector.
   */
  std::size_t size () const;

  /**
   * Return whether the vector contains only elements with value zero. This
   * function is mainly for internal consistency checks and should seldom be
   * used when not in debug mode since it uses quite some time.
   */
  bool all_zero () const;

  /**
   * Return @p true if the vector has no negative entries, i.e. all entries
   * are zero or positive. This function is used, for example, to check
   * whether refinement indicators are really all positive (or zero).
   *
   * The function obviously only makes sense if the template argument of this
   * class is a real type. If it is a complex type, then an exception is
   * thrown.
   */
  bool is_non_negative () const;

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  std::size_t memory_consumption () const;
  //@}

protected:

  /**
   * Dimension. Actual number of components contained in the vector.  Get this
   * number by calling <tt>size()</tt>.
   */
  size_type vec_size;

  /**
   * Amount of memory actually reserved for this vector. This number may be
   * greater than @p vec_size if a @p reinit was called with less memory
   * requirements than the vector needed last time. At present @p reinit does
   * not free memory when the number of needed elements is reduced.
   */
  size_type max_vec_size;

  /**
   * Pointer to the array of elements of this vector.
   */
  Number *val;

  /**
   * For parallel loops with TBB, this member variable stores the affinity
   * information of loops.
   */
  mutable std_cxx11::shared_ptr<parallel::internal::TBBPartitioner> thread_loop_partitioner;

  /**
   * Make all other vector types friends.
   */
  template <typename Number2> friend class Vector;

  /**
   * LAPACK matrices need access to the data.
   */
  template <typename Number2> friend class LAPACKFullMatrix;

  /**
   * VectorView will access the pointer.
   */
  friend class VectorView<Number>;

private:

  /**
   * Allocate and align @p val along 64-byte boundaries. The size of the
   * allocated memory is determined by @p max_vec_size .
   */
  void allocate();

  /**
   * Deallocate @p val.
   */
  void deallocate();
};

/*@}*/
/*----------------------- Inline functions ----------------------------------*/


#ifndef DOXYGEN


template <typename Number>
inline
Vector<Number>::Vector ()
  :
  vec_size(0),
  max_vec_size(0),
  val(0)
{
  reinit(0);
}



template <typename Number>
template <typename InputIterator>
Vector<Number>::Vector (const InputIterator first, const InputIterator last)
  :
  vec_size (0),
  max_vec_size (0),
  val (0)
{
  // allocate memory. do not initialize it, as we will copy over to it in a
  // second
  reinit (std::distance (first, last), true);
  std::copy (first, last, begin());
}



template <typename Number>
inline
Vector<Number>::Vector (const size_type n)
  :
  vec_size(0),
  max_vec_size(0),
  val(0)
{
  reinit (n, false);
}



template <typename Number>
inline
Vector<Number>::~Vector ()
{
  if (val)
    {
      deallocate();
      val=0;
    }
}



template <typename Number>
inline
std::size_t Vector<Number>::size () const
{
  return vec_size;
}


template <typename Number>
inline
bool Vector<Number>::in_local_range
(const size_type) const
{
  return true;
}



template <typename Number>
inline
typename Vector<Number>::iterator
Vector<Number>::begin ()
{
  return &val[0];
}



template <typename Number>
inline
typename Vector<Number>::const_iterator
Vector<Number>::begin () const
{
  return &val[0];
}



template <typename Number>
inline
typename Vector<Number>::iterator
Vector<Number>::end ()
{
  return &val[vec_size];
}



template <typename Number>
inline
typename Vector<Number>::const_iterator
Vector<Number>::end () const
{
  return &val[vec_size];
}



template <typename Number>
inline
Number Vector<Number>::operator() (const size_type i) const
{
  Assert (i<vec_size, ExcIndexRange(i,0,vec_size));
  return val[i];
}



template <typename Number>
inline
Number &Vector<Number>::operator() (const size_type i)
{
  Assert (i<vec_size, ExcIndexRangeType<size_type>(i,0,vec_size));
  return val[i];
}



template <typename Number>
inline
Number Vector<Number>::operator[] (const size_type i) const
{
  return operator()(i);
}



template <typename Number>
inline
Number &Vector<Number>::operator[] (const size_type i)
{
  return operator()(i);
}



template <typename Number>
template <typename OtherNumber>
inline
void Vector<Number>::extract_subvector_to (const std::vector<size_type> &indices,
                                           std::vector<OtherNumber> &values) const
{
  for (size_type i = 0; i < indices.size(); ++i)
    values[i] = operator()(indices[i]);
}



template <typename Number>
template <typename ForwardIterator, typename OutputIterator>
inline
void Vector<Number>::extract_subvector_to (ForwardIterator          indices_begin,
                                           const ForwardIterator    indices_end,
                                           OutputIterator           values_begin) const
{
  while (indices_begin != indices_end)
    {
      *values_begin = operator()(*indices_begin);
      indices_begin++;
      values_begin++;
    }
}



template <typename Number>
inline
Vector<Number> &
Vector<Number>::operator /= (const Number factor)
{
  AssertIsFinite(factor);
  Assert (factor != Number(0.), ExcZero() );

  this->operator *= (Number(1.)/factor);
  return *this;
}



template <typename Number>
template <typename OtherNumber>
inline
void
Vector<Number>::add (const std::vector<size_type> &indices,
                     const std::vector<OtherNumber>  &values)
{
  Assert (indices.size() == values.size(),
          ExcDimensionMismatch(indices.size(), values.size()));
  add (indices.size(), &indices[0], &values[0]);
}



template <typename Number>
template <typename OtherNumber>
inline
void
Vector<Number>::add (const std::vector<size_type> &indices,
                     const Vector<OtherNumber>    &values)
{
  Assert (indices.size() == values.size(),
          ExcDimensionMismatch(indices.size(), values.size()));
  add (indices.size(), &indices[0], values.val);
}



template <typename Number>
template <typename OtherNumber>
inline
void
Vector<Number>::add (const size_type  n_indices,
                     const size_type *indices,
                     const OtherNumber  *values)
{
  for (size_type i=0; i<n_indices; ++i)
    {
      Assert (indices[i] < vec_size, ExcIndexRange(indices[i],0,vec_size));
      Assert (numbers::is_finite(values[i]),
              ExcMessage("The given value is not finite but either infinite or Not A Number (NaN)"));

      val[indices[i]] += values[i];
    }
}



template <typename Number>
template <typename Number2>
inline
bool
Vector<Number>::operator != (const Vector<Number2> &v) const
{
  return ! (*this == v);
}



template <typename Number>
inline
void
Vector<Number>::compress (::dealii::VectorOperation::values) const
{}



template <typename Number>
inline
void
Vector<Number>::update_ghost_values () const
{}



// Moved from vector.templates.h as an inline function by Luca Heltai
// on 2009/04/12 to prevent strange compiling errors, after making
// swap virtual.
template <typename Number>
inline
void
Vector<Number>::swap (Vector<Number> &v)
{
  std::swap (vec_size,     v.vec_size);
  std::swap (max_vec_size, v.max_vec_size);
  std::swap (val,          v.val);
}



template <typename Number>
template <class Archive>
inline
void
Vector<Number>::save (Archive &ar, const unsigned int) const
{
  // forward to serialization function in the base class.
  ar &static_cast<const Subscriptor &>(*this);

  ar &vec_size &max_vec_size ;
  ar &boost::serialization::make_array(val, max_vec_size);
}



template <typename Number>
template <class Archive>
inline
void
Vector<Number>::load (Archive &ar, const unsigned int)
{
  // get rid of previous content
  deallocate();

  // the load stuff again from the archive
  ar &static_cast<Subscriptor &>(*this);
  ar &vec_size &max_vec_size ;

  allocate();
  ar &boost::serialization::make_array(val, max_vec_size);
}

#endif


/*! @addtogroup Vectors
 *@{
 */


/**
 * Global function @p swap which overloads the default implementation of the
 * C++ standard library which uses a temporary object. The function simply
 * exchanges the data of the two vectors.
 *
 * @relates Vector
 * @author Wolfgang Bangerth, 2000
 */
template <typename Number>
inline
void swap (Vector<Number> &u, Vector<Number> &v)
{
  u.swap (v);
}


/**
 * Output operator writing a vector to a stream.
 */
template <typename number>
inline
std::ostream &
operator << (std::ostream &os, const Vector<number> &v)
{
  v.print(os);
  return os;
}

/**
 * Output operator writing a vector to a LogStream.
 */
template <typename number>
inline
LogStream &
operator << (LogStream &os, const Vector<number> &v)
{
  v.print(os);
  return os;
}


/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif
