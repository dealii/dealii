// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
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

#ifndef dealii_tensor_h
#define dealii_tensor_h

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/tensor_accessors.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/utilities.h>

#include <cmath>
#include <ostream>
#include <vector>


DEAL_II_NAMESPACE_OPEN

// Forward declarations:

template <int dim, typename Number> class Point;
template <int rank_, int dim, typename Number = double> class Tensor;
template <typename Number> class Vector;
template <typename Number> class VectorizedArray;

#ifndef DOXYGEN
// Overload invalid tensor types of negative rank that come up during
// overload resolution of operator* and related contraction variants.
template <int dim, typename Number>
class Tensor<-2, dim, Number>
{
};

template <int dim, typename Number>
class Tensor<-1, dim, Number>
{
};
#endif /* DOXYGEN */


/**
 * This class is a specialized version of the <tt>Tensor<rank,dim,Number></tt>
 * class. It handles tensors of rank zero, i.e. scalars. The second template
 * argument @p dim is ignored.
 *
 * This class exists because in some cases we want to construct objects of
 * type Tensor@<spacedim-dim,dim,Number@>, which should expand to scalars,
 * vectors, matrices, etc, depending on the values of the template arguments
 * @p dim and @p spacedim. We therefore need a class that acts as a scalar
 * (i.e. @p Number) for all purposes but is part of the Tensor template
 * family.
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * this tensor operates. This of course equals the number of coordinates that
 * identify a point and rank-1 tensor. Since the current object is a rank-0
 * tensor (a scalar), this template argument has no meaning for this class.
 *
 * @tparam Number The data type in which the tensor elements are to be stored.
 * This will, in almost all cases, simply be the default @p double, but there
 * are cases where one may want to store elements in a different (and always
 * scalar) type. It can be used to base tensors on @p float or @p complex
 * numbers or any other data type that implements basic arithmetic operations.
 * Another example would be a type that allows for Automatic Differentiation
 * (see, for example, the Sacado type used in step-33) and thereby can
 * generate analytic (spatial) derivatives of a function that takes a tensor
 * as argument.
 *
 * @ingroup geomprimitives
 * @author Wolfgang Bangerth, 2009, Matthias Maier, 2015
 */
template <int dim, typename Number>
class Tensor<0,dim,Number>
{
public:
  /**
   * Provide a way to get the dimension of an object without explicit
   * knowledge of it's data type. Implementation is this way instead of
   * providing a function <tt>dimension()</tt> because now it is possible to
   * get the dimension at compile time without the expansion and preevaluation
   * of an inlined function; the compiler may therefore produce more efficient
   * code and you may use this value to declare other data types.
   */
  static const unsigned int dimension = dim;

  /**
   * Publish the rank of this tensor to the outside world.
   */
  static const unsigned int rank = 0;

  /**
   * Number of independent components of a tensor of rank 0.
   */
  static const unsigned int n_independent_components = 1;

  /**
   * Declare a type that has holds real-valued numbers with the same precision
   * as the template argument to this class. For std::complex<number>, this
   * corresponds to type number, and it is equal to Number for all other
   * cases. See also the respective field in Vector<Number>.
   *
   * This typedef is used to represent the return type of norms.
   */
  typedef typename numbers::NumberTraits<Number>::real_type real_type;

  /**
   * Type of objects encapsulated by this container and returned by
   * operator[](). This is a scalar number type for a rank 0 tensor.
   */
  typedef Number value_type;

  /**
   * Declare an array type which can be used to initialize an object of this
   * type statically. In case of a tensor of rank 0 this is just the scalar
   * number type Number.
   */
  typedef Number array_type;

  /**
   * Constructor. Set to zero.
   *
   * @see CUDAWrappers
   */
  DEAL_II_CUDA_HOST_DEV Tensor ();

  /**
   * Constructor from tensors with different underlying scalar type. This
   * obviously requires that the @p OtherNumber type is convertible to @p
   * Number.
   */
  template <typename OtherNumber>
  Tensor (const Tensor<0,dim,OtherNumber> &initializer);

  /**
   * Constructor, where the data is copied from a C-style array.
   */
  template <typename OtherNumber>
  Tensor (const OtherNumber &initializer);

  /**
   * Return a pointer to the first element of the underlying storage.
   */
  Number *
  begin_raw();

  /**
   * Return a const pointer to the first element of the underlying storage.
   */
  const Number *
  begin_raw() const;

  /**
   * Return a pointer to the element past the end of the underlying storage.
   */
  Number *
  end_raw();

  /**
   * Return a const pointer to the element past the end of the underlying
   * storage.
   */
  const Number *
  end_raw() const;

  /**
   * Return a reference to the encapsulated Number object. Since rank-0
   * tensors are scalars, this is a natural operation.
   *
   * This is the non-const conversion operator that returns a writable
   * reference.
   *
   * @see CUDAWrappers
   */
  DEAL_II_CUDA_HOST_DEV operator Number &();

  /**
   * Return a reference to the encapsulated Number object. Since rank-0
   * tensors are scalars, this is a natural operation.
   *
   * This is the const conversion operator that returns a read-only reference.
   *
   * @see CUDAWrappers
   */
  DEAL_II_CUDA_HOST_DEV operator const Number &() const;

  /**
   * Assignment from tensors with different underlying scalar type. This
   * obviously requires that the @p OtherNumber type is convertible to @p
   * Number.
   */
  template <typename OtherNumber>
  Tensor &operator = (const Tensor<0,dim,OtherNumber> &rhs);

  /**
   * Assignment from tensors with same underlying scalar type.
   */
  Tensor &operator = (const Tensor<0,dim,Number> &rhs);

  /**
   * This operator assigns a scalar to a tensor. This obviously requires
   * that the @p OtherNumber type is convertible to @p Number.
   */
  template <typename OtherNumber>
  Tensor &operator = (const OtherNumber &d);

  /**
   * Test for equality of two tensors.
   */
  template <typename OtherNumber>
  bool operator == (const Tensor<0,dim,OtherNumber> &rhs) const;

  /**
   * Test for inequality of two tensors.
   */
  template <typename OtherNumber>
  bool operator != (const Tensor<0,dim,OtherNumber> &rhs) const;

  /**
   * Add another scalar
   */
  template <typename OtherNumber>
  Tensor &operator += (const Tensor<0,dim,OtherNumber> &rhs);

  /**
   * Subtract another scalar.
   */
  template <typename OtherNumber>
  Tensor &operator -= (const Tensor<0,dim,OtherNumber> &rhs);

  /**
   * Multiply the scalar with a <tt>factor</tt>.
   *
   * @see CUDAWrappers
   */
  template <typename OtherNumber>
  DEAL_II_CUDA_HOST_DEV Tensor &operator *= (const OtherNumber &factor);

  /**
   * Divide the scalar by <tt>factor</tt>.
   */
  template <typename OtherNumber>
  Tensor &operator /= (const OtherNumber &factor);

  /**
   * Tensor with inverted entries.
   */
  Tensor  operator - () const;

  /**
   * Reset all values to zero.
   *
   * Note that this is partly inconsistent with the semantics of the @p
   * clear() member functions of the standard library containers and of
   * several other classes within deal.II, which not only reset the values of
   * stored elements to zero, but release all memory and return the object
   * into a virginial state. However, since the size of objects of the present
   * type is determined by its template parameters, resizing is not an option,
   * and indeed the state where all elements have a zero value is the state
   * right after construction of such an object.
   */
  void clear ();

  /**
   * Return the Frobenius-norm of a tensor, i.e. the square root of the sum of
   * the absolute squares of all entries. For the present case of rank-1
   * tensors, this equals the usual <tt>l<sub>2</sub></tt> norm of the vector.
   */
  real_type norm () const;

  /**
   * Return the square of the Frobenius-norm of a tensor, i.e. the sum of the
   * absolute squares of all entries.
   *
   * @see CUDAWrappers
   */
  DEAL_II_CUDA_HOST_DEV real_type norm_square () const;

  /**
   * Read or write the data of this object to or from a stream for the purpose
   * of serialization
   */
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version);

  /**
   * Internal type declaration that is used to specialize the return type of
   * operator[]() for Tensor<1,dim,Number>
   */
  typedef Number tensor_type;

private:
  /**
   * The value of this scalar object.
   */
  Number value;

  /**
   * Internal helper function for unroll.
   */
  template <typename OtherNumber>
  void unroll_recursion(Vector<OtherNumber> &result,
                        unsigned int        &start_index) const;

  /**
   * Allow an arbitrary Tensor to access the underlying values.
   */
  template <int, int, typename> friend class Tensor;
};



/**
 * A general tensor class with an arbitrary rank, i.e. with an arbitrary
 * number of indices. The Tensor class provides an indexing operator and a bit
 * of infrastructure, but most functionality is recursively handed down to
 * tensors of rank 1 or put into external templated functions, e.g. the
 * <tt>contract</tt> family.
 *
 * Using this tensor class for objects of rank 2 has advantages over matrices
 * in many cases since the dimension is known to the compiler as well as the
 * location of the data. It is therefore possible to produce far more
 * efficient code than for matrices with runtime-dependent dimension. It also
 * makes the code easier to read because of the semantic difference between a
 * tensor (an object that relates to a coordinate system and has
 * transformation properties with regard to coordinate rotations and
 * transforms) and matrices (which we consider as operators on arbitrary
 * vector spaces related to linear algebra things).
 *
 * @tparam rank_ An integer that denotes the rank of this tensor. A rank-0
 * tensor is a scalar, a rank-1 tensor is a vector with @p dim components, a
 * rank-2 tensor is a matrix with dim-by-dim components, etc. There are
 * specializations of this class for rank-0 and rank-1 tensors. There is also
 * a related class SymmetricTensor for tensors of even rank whose elements are
 * symmetric.
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * this tensor operates. This of course equals the number of coordinates that
 * identify a point and rank-1 tensor.
 *
 * @tparam Number The data type in which the tensor elements are to be stored.
 * This will, in almost all cases, simply be the default @p double, but there
 * are cases where one may want to store elements in a different (and always
 * scalar) type. It can be used to base tensors on @p float or @p complex
 * numbers or any other data type that implements basic arithmetic operations.
 * Another example would be a type that allows for Automatic Differentiation
 * (see, for example, the Sacado type used in step-33) and thereby can
 * generate analytic (spatial) derivatives of a function that takes a tensor
 * as argument.
 *
 * @ingroup geomprimitives
 * @author Wolfgang Bangerth, 1998-2005, Matthias Maier, 2015
 */
template <int rank_, int dim, typename Number>
class Tensor
{
public:
  /**
   * Provide a way to get the dimension of an object without explicit
   * knowledge of it's data type. Implementation is this way instead of
   * providing a function <tt>dimension()</tt> because now it is possible to
   * get the dimension at compile time without the expansion and preevaluation
   * of an inlined function; the compiler may therefore produce more efficient
   * code and you may use this value to declare other data types.
   */
  static const unsigned int dimension = dim;

  /**
   * Publish the rank of this tensor to the outside world.
   */
  static const unsigned int rank = rank_;

  /**
   * Number of independent components of a tensor of current rank. This is dim
   * times the number of independent components of each sub-tensor.
   */
  static const unsigned int
  n_independent_components = Tensor<rank_-1,dim>::n_independent_components *dim;

  /**
   * Type of objects encapsulated by this container and returned by
   * operator[](). This is a tensor of lower rank for a general tensor, and a
   * scalar number type for Tensor<1,dim,Number>.
   */
  typedef typename Tensor<rank_-1,dim,Number>::tensor_type value_type;

  /**
   * Declare an array type which can be used to initialize an object of this
   * type statically.
   */
  typedef typename Tensor<rank_-1,dim,Number>::array_type
  array_type[(dim != 0) ? dim : 1];
  // ... avoid a compiler warning in case of dim == 0 and ensure that the
  // array always has positive size.

  /**
   * Constructor. Initialize all entries to zero.
   *
   * @see CUDAWrappers
   */
  DEAL_II_CUDA_HOST_DEV Tensor ();

  /**
   * Constructor, where the data is copied from a C-style array.
   */
  explicit Tensor (const array_type &initializer);

  /**
   * Constructor from tensors with different underlying scalar type. This
   * obviously requires that the @p OtherNumber type is convertible to @p
   * Number.
   */
  template <typename OtherNumber>
  Tensor (const Tensor<rank_,dim,OtherNumber> &initializer);

  /**
   * Constructor that converts from a "tensor of tensors".
   */
  template <typename OtherNumber>
  Tensor (const Tensor<1,dim,Tensor<rank_-1,dim,OtherNumber> > &initializer);

  /**
   * Conversion operator to tensor of tensors.
   */
  template <typename OtherNumber>
  operator Tensor<1,dim,Tensor<rank_-1,dim,OtherNumber> > () const;

  /**
   * Read-Write access operator.
   *
   * @see CUDAWrappers
   */
  DEAL_II_CUDA_HOST_DEV value_type &operator [] (const unsigned int i);

  /**
   * Read-only access operator.
   *
   * @see CUDAWrappers
   */
  DEAL_II_CUDA_HOST_DEV const value_type &operator[](const unsigned int i) const;

  /**
   * Read access using TableIndices <tt>indices</tt>
   */
  const Number &operator [] (const TableIndices<rank_> &indices) const;

  /**
   * Read and write access using TableIndices <tt>indices</tt>
   */
  Number &operator [] (const TableIndices<rank_> &indices);

  /**
   * Return a pointer to the first element of the underlying storage.
   */
  Number *
  begin_raw();

  /**
   * Return a const pointer to the first element of the underlying storage.
   */
  const Number *
  begin_raw() const;

  /**
   * Return a pointer to the element past the end of the underlying storage.
   */
  Number *
  end_raw();

  /**
   * Return a pointer to the element past the end of the underlying storage.
   */
  const Number *
  end_raw() const;

  /**
   * Assignment operator from tensors with different underlying scalar type.
   * This obviously requires that the @p OtherNumber type is convertible to @p
   * Number.
   */
  template <typename OtherNumber>
  Tensor &operator = (const Tensor<rank_,dim,OtherNumber> &rhs);

  /**
   * This operator assigns a scalar to a tensor. To avoid confusion with what
   * exactly it means to assign a scalar value to a tensor, zero is the only
   * value allowed for <tt>d</tt>, allowing the intuitive notation
   * <tt>t=0</tt> to reset all elements of the tensor to zero.
   */
  Tensor &operator = (const Number &d);

  /**
   * Test for equality of two tensors.
   */
  template <typename OtherNumber>
  bool operator == (const Tensor<rank_,dim,OtherNumber> &) const;

  /**
   * Test for inequality of two tensors.
   */
  template <typename OtherNumber>
  bool operator != (const Tensor<rank_,dim,OtherNumber> &) const;

  /**
   * Add another tensor.
   */
  template <typename OtherNumber>
  Tensor &operator += (const Tensor<rank_,dim,OtherNumber> &);

  /**
   * Subtract another tensor.
   */
  template <typename OtherNumber>
  Tensor &operator -= (const Tensor<rank_,dim,OtherNumber> &);

  /**
   * Scale the tensor by <tt>factor</tt>, i.e. multiply all components by
   * <tt>factor</tt>.
   *
   * @see CUDAWrappers
   */
  template <typename OtherNumber>
  DEAL_II_CUDA_HOST_DEV Tensor &operator *= (const OtherNumber &factor);

  /**
   * Scale the vector by <tt>1/factor</tt>.
   */
  template <typename OtherNumber>
  Tensor &operator /= (const OtherNumber &factor);

  /**
   * Unary minus operator. Negate all entries of a tensor.
   */
  Tensor  operator - () const;

  /**
   * Reset all values to zero.
   *
   * Note that this is partly inconsistent with the semantics of the @p
   * clear() member functions of the standard library containers and of
   * several other classes within deal.II, which not only reset the values of
   * stored elements to zero, but release all memory and return the object
   * into a virginial state. However, since the size of objects of the present
   * type is determined by its template parameters, resizing is not an option,
   * and indeed the state where all elements have a zero value is the state
   * right after construction of such an object.
   */
  void clear ();

  /**
   * Return the Frobenius-norm of a tensor, i.e. the square root of the sum of
   * the absolute squares of all entries. For the present case of rank-1
   * tensors, this equals the usual <tt>l<sub>2</sub></tt> norm of the vector.
   */

  typename numbers::NumberTraits<Number>::real_type norm() const;

  /**
   * Return the square of the Frobenius-norm of a tensor, i.e. the sum of the
   * absolute squares of all entries.
   *
   * @see CUDAWrappers
   */
  DEAL_II_CUDA_HOST_DEV typename numbers::NumberTraits<Number>::real_type norm_square() const;

  /**
   * Fill a vector with all tensor elements.
   *
   * This function unrolls all tensor entries into a single, linearly numbered
   * vector. As usual in C++, the rightmost index of the tensor marches
   * fastest.
   */
  template <typename OtherNumber>
  void unroll (Vector<OtherNumber> &result) const;

  /**
   * Return an unrolled index in the range [0,dim^rank-1] for the element of
   * the tensor indexed by the argument to the function.
   */
  static
  unsigned int
  component_to_unrolled_index(const TableIndices<rank_> &indices);

  /**
   * Opposite of  component_to_unrolled_index: For an index in the range
   * [0,dim^rank-1], return which set of indices it would correspond to.
   */
  static
  TableIndices<rank_> unrolled_to_component_indices(const unsigned int i);

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  static std::size_t memory_consumption ();

  /**
   * Read or write the data of this object to or from a stream for the purpose
   * of serialization
   */
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version);

  /**
   * Internal type declaration that is used to specialize the return type of
   * operator[]() for Tensor<1,dim,Number>
   */
  typedef Tensor<rank_, dim, Number> tensor_type;

private:
  /**
   * Array of tensors holding the subelements.
   */
  Tensor<rank_-1, dim, Number> values[(dim != 0) ? dim : 1];
  // ... avoid a compiler warning in case of dim == 0 and ensure that the
  // array always has positive size.

  /**
   * Internal helper function for unroll.
   */
  template <typename OtherNumber>
  void unroll_recursion(Vector<OtherNumber> &result,
                        unsigned int        &start_index) const;

  /**
   * Allow an arbitrary Tensor to access the underlying values.
   */
  template <int, int, typename> friend class Tensor;

  /**
   * Point is allowed access to the coordinates. This is supposed to improve
   * speed.
   */
  friend class Point<dim,Number>;
};


namespace internal
{
  /**
  * The structs below are needed since VectorizedArray<T1> is a POD-type without a constructor and
  * can be a template argument for Tensor<...,T2> where T2 would equal Tensor<1, dim, VectorizedArray >.
  * Internally, in previous versions of deal.II, Tensor<...,T2> would make use of the constructor
  * of T2 leading to a compile-time error. However simply adding a constructor for VectorizedArray<T1>
  * breaks the POD-idioms needed elsewhere. Calls to constructors of T2 subsequently got replaced by a
  * call to internal::NumberType<T2> which then determines the right function to use by template deduction.
  * A detailed discussion can be found at https://github.com/dealii/dealii/pull/3967 . Also see
  * numbers.h for another specialization.
   */
  template <int rank, int dim, typename T>
  struct NumberType<Tensor<rank,dim,T> >
  {
    static const Tensor<rank,dim,T> &value (const Tensor<rank,dim,T> &t)
    {
      return t;
    }

    static Tensor<rank,dim,T> value (const T &t)
    {
      Tensor<rank,dim,T> tmp;
      tmp=t;
      return tmp;
    }
  };

  template <int rank, int dim, typename T>
  struct NumberType<Tensor<rank,dim,VectorizedArray<T> > >
  {
    static const Tensor<rank,dim,VectorizedArray<T> > &value (const Tensor<rank,dim,VectorizedArray<T> > &t)
    {
      return t;
    }

    static Tensor<rank,dim,VectorizedArray<T> > value (const T &t)
    {
      Tensor<rank,dim,VectorizedArray<T> > tmp;
      tmp=internal::NumberType<VectorizedArray<T> >::value(t);
      return tmp;
    }

    static Tensor<rank,dim,VectorizedArray<T> > value (const VectorizedArray<T> &t)
    {
      Tensor<rank,dim,VectorizedArray<T> > tmp;
      tmp=t;
      return tmp;
    }
  };
}


/*---------------------- Inline functions: Tensor<0,dim> ---------------------*/


template <int dim,typename Number>
inline
DEAL_II_CUDA_HOST_DEV Tensor<0,dim,Number>::Tensor ()
// Some auto-differentiable numbers need explicit
// zero initialization.
  : value(internal::NumberType<Number>::value(0.0))
{
}



template <int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<0,dim,Number>::Tensor (const OtherNumber &initializer)
{
  value = internal::NumberType<Number>::value(initializer);
}



template <int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<0,dim,Number>::Tensor (const Tensor<0,dim,OtherNumber> &p)
{
  value = p.value;
}



template <int dim, typename Number>
inline
Number *
Tensor<0,dim,Number>::begin_raw()
{
  return std::addressof(value);
}



template <int dim, typename Number>
inline
const Number *
Tensor<0,dim,Number>::begin_raw() const
{
  return std::addressof(value);
}



template <int dim, typename Number>
inline
Number *
Tensor<0,dim,Number>::end_raw()
{
  return begin_raw()+n_independent_components;
}



template <int dim, typename Number>
inline
const Number *
Tensor<0,dim,Number>::end_raw() const
{
  return begin_raw()+n_independent_components;
}



template <int dim, typename Number>
inline
DEAL_II_CUDA_HOST_DEV Tensor<0,dim,Number>::operator Number &()
{
  // We cannot use Assert inside a CUDA kernel
#ifndef DEAL_II_WITH_CUDA
  Assert(dim != 0, ExcMessage("Cannot access an object of type Tensor<0,0,Number>"));
#endif
  return value;
}


template <int dim, typename Number>
inline
DEAL_II_CUDA_HOST_DEV Tensor<0,dim,Number>::operator const Number &() const
{
  // We cannot use Assert inside a CUDA kernel
#ifndef DEAL_II_WITH_CUDA
  Assert(dim != 0, ExcMessage("Cannot access an object of type Tensor<0,0,Number>"));
#endif
  return value;
}


template <int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<0,dim,Number> &Tensor<0,dim,Number>::operator = (const Tensor<0,dim,OtherNumber> &p)
{
  value = internal::NumberType<Number>::value(p);
  return *this;
}


template <int dim, typename Number>
inline
Tensor<0,dim,Number> &Tensor<0,dim,Number>::operator = (const Tensor<0,dim,Number> &p)
{
  value = p.value;
  return *this;
}


template <int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<0,dim,Number> &Tensor<0,dim,Number>::operator = (const OtherNumber &d)
{
  value = internal::NumberType<Number>::value(d);
  return *this;
}


template <int dim, typename Number>
template <typename OtherNumber>
inline
bool Tensor<0,dim,Number>::operator == (const Tensor<0,dim,OtherNumber> &p) const
{
  return (value == p.value);
}


template <int dim, typename Number>
template <typename OtherNumber>
inline
bool Tensor<0,dim,Number>::operator != (const Tensor<0,dim,OtherNumber> &p) const
{
  return !((*this) == p);
}


template <int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<0,dim,Number> &Tensor<0,dim,Number>::operator += (const Tensor<0,dim,OtherNumber> &p)
{
  value += p.value;
  return *this;
}


template <int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<0,dim,Number> &Tensor<0,dim,Number>::operator -= (const Tensor<0,dim,OtherNumber> &p)
{
  value -= p.value;
  return *this;
}


template <int dim, typename Number>
template <typename OtherNumber>
inline
DEAL_II_CUDA_HOST_DEV Tensor<0,dim,Number> &Tensor<0,dim,Number>::operator *= (const OtherNumber &s)
{
  value *= s;
  return *this;
}


template <int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<0,dim,Number> &Tensor<0,dim,Number>::operator /= (const OtherNumber &s)
{
  value /= s;
  return *this;
}


template <int dim, typename Number>
inline
Tensor<0,dim,Number> Tensor<0,dim,Number>::operator - () const
{
  return -value;
}


template <int dim, typename Number>
inline
typename Tensor<0,dim,Number>::real_type
Tensor<0,dim,Number>::norm () const
{
  Assert(dim != 0, ExcMessage("Cannot access an object of type Tensor<0,0,Number>"));
  return numbers::NumberTraits<Number>::abs (value);
}


template <int dim, typename Number>
inline
typename Tensor<0,dim,Number>::real_type
DEAL_II_CUDA_HOST_DEV Tensor<0,dim,Number>::norm_square () const
{
  Assert(dim != 0, ExcMessage("Cannot access an object of type Tensor<0,0,Number>"));
  return numbers::NumberTraits<Number>::abs_square (value);
}


template <int dim, typename Number>
template <typename OtherNumber>
inline
void
Tensor<0, dim, Number>::unroll_recursion (Vector<OtherNumber> &result,
                                          unsigned int        &index) const
{
  Assert(dim != 0, ExcMessage("Cannot unroll an object of type Tensor<0,0,Number>"));
  result[index] = value;
  ++index;
}


template <int dim, typename Number>
inline
void Tensor<0,dim,Number>::clear ()
{
  // Some auto-differentiable numbers need explicit
  // zero initialization.
  value = internal::NumberType<Number>::value(0.0);
}


template <int dim, typename Number>
template <class Archive>
inline
void Tensor<0,dim,Number>::serialize(Archive &ar, const unsigned int)
{
  ar &value;
}


/*-------------------- Inline functions: Tensor<rank,dim> --------------------*/


template <int rank_, int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE
DEAL_II_CUDA_HOST_DEV Tensor<rank_,dim,Number>::Tensor ()
{
  // All members of the c-style array values are already default initialized
  // and thus all values are already set to zero recursively.
}


template <int rank_, int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE
Tensor<rank_,dim,Number>::Tensor (const array_type &initializer)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] = Tensor<rank_-1, dim, Number>(initializer[i]);
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline DEAL_II_ALWAYS_INLINE
Tensor<rank_,dim,Number>::Tensor (const Tensor<rank_,dim,OtherNumber> &initializer)
{
  for (unsigned int i=0; i!=dim; ++i)
    values[i] = Tensor<rank_-1,dim,Number>(initializer[i]);
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline DEAL_II_ALWAYS_INLINE
Tensor<rank_,dim,Number>::Tensor
(const Tensor<1,dim,Tensor<rank_-1,dim,OtherNumber> > &initializer)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] = initializer[i];
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline DEAL_II_ALWAYS_INLINE
Tensor<rank_,dim,Number>::
operator Tensor<1,dim,Tensor<rank_-1,dim,OtherNumber> > () const
{
  return Tensor<1,dim,Tensor<rank_-1,dim,Number> > (values);
}



namespace internal
{
  namespace TensorSubscriptor
  {
    template <typename ArrayElementType, int dim>
    inline DEAL_II_ALWAYS_INLINE
    DEAL_II_CUDA_HOST_DEV
    ArrayElementType &
    subscript (ArrayElementType *values,
               const unsigned int i,
               std::integral_constant<int, dim>)
    {
      // We cannot use Assert in a CUDA kernel
#ifndef DEAL_II_WITH_CUDA
      Assert (i<dim, ExcIndexRange(i, 0, dim));
#endif
      return values[i];
    }


    template <typename ArrayElementType>
    DEAL_II_CUDA_HOST_DEV
    ArrayElementType &
    subscript (ArrayElementType *,
               const unsigned int,
               std::integral_constant<int, 0>)
    {
      // We cannot use Assert in a CUDA kernel
#ifndef DEAL_II_WITH_CUDA
      Assert(false, ExcMessage("Cannot access elements of an object of type Tensor<rank,0,Number>."));
#endif
      static ArrayElementType t;
      return t;
    }
  }
}


template <int rank_, int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE
DEAL_II_CUDA_HOST_DEV
typename Tensor<rank_,dim,Number>::value_type &
Tensor<rank_,dim,Number>::operator[] (const unsigned int i)
{
  return dealii::internal::TensorSubscriptor::subscript(values, i, std::integral_constant<int, dim>());
}


template <int rank_, int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE
DEAL_II_CUDA_HOST_DEV
const typename Tensor<rank_,dim,Number>::value_type &
Tensor<rank_,dim,Number>::operator[] (const unsigned int i) const
{
  return dealii::internal::TensorSubscriptor::subscript(values, i, std::integral_constant<int, dim>());
}


template <int rank_, int dim, typename Number>
inline
const Number &
Tensor<rank_,dim,Number>::operator[] (const TableIndices<rank_> &indices) const
{
  Assert(dim != 0, ExcMessage("Cannot access an object of type Tensor<rank_,0,Number>"));

  return TensorAccessors::extract<rank_>(*this, indices);
}



template <int rank_, int dim, typename Number>
inline
Number &
Tensor<rank_,dim,Number>::operator[] (const TableIndices<rank_> &indices)
{
  Assert(dim != 0, ExcMessage("Cannot access an object of type Tensor<rank_,0,Number>"));

  return TensorAccessors::extract<rank_>(*this, indices);
}



template <int rank_, int dim, typename Number>
inline
Number *
Tensor<rank_,dim,Number>::begin_raw()
{
  return std::addressof(this->operator[](this->unrolled_to_component_indices(0)));
}



template <int rank_, int dim, typename Number>
inline
const Number *
Tensor<rank_,dim,Number>::begin_raw() const
{
  return std::addressof(this->operator[](this->unrolled_to_component_indices(0)));
}



template <int rank_, int dim, typename Number>
inline
Number *
Tensor<rank_,dim,Number>::end_raw()
{
  return begin_raw()+n_independent_components;
}



template <int rank_, int dim, typename Number>
inline
const Number *
Tensor<rank_,dim,Number>::end_raw() const
{
  return begin_raw()+n_independent_components;
}



template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline DEAL_II_ALWAYS_INLINE
Tensor<rank_,dim,Number> &
Tensor<rank_,dim,Number>::operator = (const Tensor<rank_,dim,OtherNumber> &t)
{
  if (dim > 0)
    std::copy (&t.values[0], &t.values[0]+dim, &values[0]);
  return *this;
}


template <int rank_, int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE
Tensor<rank_,dim,Number> &
Tensor<rank_,dim,Number>::operator = (const Number &d)
{
  Assert (d == internal::NumberType<Number>::value(0.0),
          ExcMessage ("Only assignment with zero is allowed"));
  (void) d;

  for (unsigned int i=0; i<dim; ++i)
    values[i] = internal::NumberType<Number>::value(0.0);
  return *this;
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline
bool
Tensor<rank_,dim,Number>::operator == (const Tensor<rank_,dim,OtherNumber> &p) const
{
  for (unsigned int i=0; i<dim; ++i)
    if (values[i] != p.values[i])
      return false;
  return true;
}


// At some places in the library, we have Point<0> for formal reasons
// (e.g., we sometimes have Quadrature<dim-1> for faces, so we have
// Quadrature<0> for dim=1, and then we have Point<0>). To avoid warnings
// in the above function that the loop end check always fails, we
// implement this function here
template <>
template <>
inline
bool Tensor<1,0,double>::operator == (const Tensor<1,0,double> &) const
{
  return true;
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline
bool
Tensor<rank_,dim,Number>::operator != (const Tensor<rank_,dim,OtherNumber> &p) const
{
  return !((*this) == p);
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<rank_,dim,Number> &
Tensor<rank_,dim,Number>::operator += (const Tensor<rank_,dim,OtherNumber> &p)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] += p.values[i];
  return *this;
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<rank_,dim,Number> &
Tensor<rank_,dim,Number>::operator -= (const Tensor<rank_,dim,OtherNumber> &p)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] -= p.values[i];
  return *this;
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline
DEAL_II_CUDA_HOST_DEV
Tensor<rank_,dim,Number> &
Tensor<rank_,dim,Number>::operator *= (const OtherNumber &s)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] *= s;
  return *this;
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<rank_,dim,Number> &
Tensor<rank_,dim,Number>::operator /= (const OtherNumber &s)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] /= s;
  return *this;
}


template <int rank_, int dim, typename Number>
inline
Tensor<rank_,dim,Number>
Tensor<rank_,dim,Number>::operator - () const
{
  Tensor<rank_,dim,Number> tmp;

  for (unsigned int i=0; i<dim; ++i)
    tmp.values[i] = -values[i];

  return tmp;
}


template <int rank_, int dim, typename Number>
inline
typename numbers::NumberTraits<Number>::real_type
Tensor<rank_,dim,Number>::norm () const
{
  return std::sqrt (norm_square());
}


template <int rank_, int dim, typename Number>
inline
DEAL_II_CUDA_HOST_DEV
typename numbers::NumberTraits<Number>::real_type
Tensor<rank_,dim,Number>::norm_square () const
{
  typename numbers::NumberTraits<Number>::real_type s
    = internal::NumberType<typename numbers::NumberTraits<Number>::real_type>::value(0.0);
  for (unsigned int i=0; i<dim; ++i)
    s += values[i].norm_square();

  return s;
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline
void
Tensor<rank_, dim, Number>::unroll (Vector<OtherNumber> &result) const
{
  AssertDimension (result.size(),(Utilities::fixed_power<rank_, unsigned int>(dim)));

  unsigned int index = 0;
  unroll_recursion (result, index);
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline
void
Tensor<rank_, dim, Number>::unroll_recursion (Vector<OtherNumber> &result,
                                              unsigned int        &index) const
{
  for (unsigned int i=0; i<dim; ++i)
    values[i].unroll_recursion(result, index);
}


template <int rank_, int dim, typename Number>
inline
unsigned int
Tensor<rank_, dim, Number>::component_to_unrolled_index(const TableIndices<rank_> &indices)
{
  unsigned int index = 0;
  for (int r = 0; r < rank_; ++r)
    index = index * dim + indices[r];

  return index;
}


template <int rank_, int dim, typename Number>
inline
TableIndices<rank_>
Tensor<rank_, dim, Number>::unrolled_to_component_indices(const unsigned int i)
{
  Assert (i < n_independent_components,
          ExcIndexRange (i, 0, n_independent_components));

  TableIndices<rank_>   indices;

  unsigned int remainder = i;
  for (int r=rank_-1; r>=0; --r)
    {
      indices[r] = (remainder % dim);
      remainder /= dim;
    }
  Assert (remainder == 0, ExcInternalError());

  return indices;
}


template <int rank_, int dim, typename Number>
inline
void Tensor<rank_,dim,Number>::clear ()
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] = internal::NumberType<Number>::value(0.0);
}


template <int rank_, int dim, typename Number>
inline
std::size_t
Tensor<rank_,dim,Number>::memory_consumption ()
{
  return sizeof(Tensor<rank_,dim,Number>);
}


template <int rank_, int dim, typename Number>
template <class Archive>
inline
void
Tensor<rank_,dim,Number>::serialize(Archive &ar, const unsigned int)
{
  ar &values;
}


/* ----------------- Non-member functions operating on tensors. ------------ */

/**
 * @name Output functions for Tensor objects
 */
//@{

/**
 * Output operator for tensors. Print the elements consecutively, with a space
 * in between, two spaces between rank 1 subtensors, three between rank 2 and
 * so on.
 *
 * @relatesalso Tensor
 */
template <int rank_, int dim, typename Number>
inline
std::ostream &operator << (std::ostream &out, const Tensor<rank_,dim,Number> &p)
{
  for (unsigned int i = 0; i < dim; ++i)
    {
      out << p[i];
      if (i != dim - 1)
        out << ' ';
    }

  return out;
}


/**
 * Output operator for tensors of rank 0. Since such tensors are scalars, we
 * simply print this one value.
 *
 * @relatesalso Tensor<0,dim,Number>
 */
template <int dim, typename Number>
inline
std::ostream &operator << (std::ostream &out, const Tensor<0,dim,Number> &p)
{
  out << static_cast<const Number &>(p);
  return out;
}


//@}
/**
 * @name Vector space operations on Tensor objects:
 */
//@{


/**
 * Scalar multiplication of a tensor of rank 0 with an object from the left.
 *
 * This function unwraps the underlying @p Number stored in the Tensor and
 * multiplies @p object with it.
 *
 * @relatesalso Tensor<0,dim,Number>
 */
template <int dim, typename Number, typename Other>
inline DEAL_II_ALWAYS_INLINE
typename ProductType<Other, Number>::type
operator * (const Other                &object,
            const Tensor<0,dim,Number> &t)
{
  return object * static_cast<const Number &>(t);
}



/**
 * Scalar multiplication of a tensor of rank 0 with an object from the right.
 *
 * This function unwraps the underlying @p Number stored in the Tensor and
 * multiplies @p object with it.
 *
 * @relatesalso Tensor<0,dim,Number>
 */
template <int dim, typename Number, typename Other>
inline DEAL_II_ALWAYS_INLINE
typename ProductType<Number, Other>::type
operator * (const Tensor<0,dim,Number> &t,
            const Other                &object)
{
  return static_cast<const Number &>(t) * object;
}


/**
 * Scalar multiplication of two tensors of rank 0.
 *
 * This function unwraps the underlying objects of type @p Number and @p
 * OtherNumber that are stored within the Tensor and multiplies them. It
 * returns an unwrapped number of product type.
 *
 * @relatesalso Tensor<0,dim,Number>
 */
template <int dim, typename Number, typename OtherNumber>
inline DEAL_II_ALWAYS_INLINE
typename ProductType<Number, OtherNumber>::type
operator * (const Tensor<0, dim, Number>      &src1,
            const Tensor<0, dim, OtherNumber> &src2)
{
  return static_cast<const Number &>(src1) *
         static_cast<const OtherNumber &>(src2);
}


/**
 * Division of a tensor of rank 0 by a scalar number.
 *
 * @relatesalso Tensor<0,dim,Number>
 */
template <int dim, typename Number, typename OtherNumber>
inline DEAL_II_ALWAYS_INLINE
Tensor<0,dim,typename ProductType<Number, typename EnableIfScalar<OtherNumber>::type>::type>
operator / (const Tensor<0,dim,Number> &t,
            const OtherNumber          &factor)
{
  return static_cast<const Number &>(t) / factor;
}


/**
 * Add two tensors of rank 0.
 *
 * @relatesalso Tensor<0,dim,Number>
 */
template <int dim, typename Number, typename OtherNumber>
inline DEAL_II_ALWAYS_INLINE
Tensor<0, dim, typename ProductType<Number, OtherNumber>::type>
operator+ (const Tensor<0,dim,Number>      &p,
           const Tensor<0,dim,OtherNumber> &q)
{
  return static_cast<const Number &>(p) + static_cast<const OtherNumber &>(q);
}


/**
 * Subtract two tensors of rank 0.
 *
 * @relatesalso Tensor<0,dim,Number>
 */
template <int dim, typename Number, typename OtherNumber>
inline DEAL_II_ALWAYS_INLINE
Tensor<0, dim, typename ProductType<Number, OtherNumber>::type>
operator- (const Tensor<0,dim,Number>      &p,
           const Tensor<0,dim,OtherNumber> &q)
{
  return static_cast<const Number &>(p) - static_cast<const OtherNumber &>(q);
}


/**
 * Multiplication of a tensor of general rank with a scalar number from the
 * right.
 *
 * Only multiplication with a scalar number type (i.e., a floating point
 * number, a complex floating point number, etc.) is allowed, see the
 * documentation of EnableIfScalar for details.
 *
 * @relatesalso Tensor
 */
template <int rank, int dim,
          typename Number,
          typename OtherNumber>
inline DEAL_II_ALWAYS_INLINE
Tensor<rank,dim,typename ProductType<Number, typename EnableIfScalar<OtherNumber>::type>::type>
operator * (const Tensor<rank,dim,Number> &t,
            const OtherNumber             &factor)
{
  // recurse over the base objects
  Tensor<rank,dim,typename ProductType<Number,OtherNumber>::type> tt;
  for (unsigned int d=0; d<dim; ++d)
    tt[d] = t[d] * factor;
  return tt;
}


/**
 * Multiplication of a tensor of general rank with a scalar number from the
 * left.
 *
 * Only multiplication with a scalar number type (i.e., a floating point
 * number, a complex floating point number, etc.) is allowed, see the
 * documentation of EnableIfScalar for details.
 *
 * @relatesalso Tensor
 */
template <int rank, int dim,
          typename Number,
          typename OtherNumber>
inline DEAL_II_ALWAYS_INLINE
Tensor<rank,dim,typename ProductType<typename EnableIfScalar<Number>::type, OtherNumber>::type>
operator * (const Number                       &factor,
            const Tensor<rank,dim,OtherNumber> &t)
{
  // simply forward to the operator above
  return t * factor;
}


/**
 * Division of a tensor of general rank with a scalar number. See the
 * discussion on operator*() above for more information about template
 * arguments and the return type.
 *
 * @relatesalso Tensor
 */
template <int rank, int dim,
          typename Number,
          typename OtherNumber>
inline
Tensor<rank,dim,typename ProductType<Number, typename EnableIfScalar<OtherNumber>::type>::type>
operator / (const Tensor<rank,dim,Number> &t,
            const OtherNumber             &factor)
{
  // recurse over the base objects
  Tensor<rank,dim,typename ProductType<Number,OtherNumber>::type> tt;
  for (unsigned int d=0; d<dim; ++d)
    tt[d] = t[d] / factor;
  return tt;
}


/**
 * Addition of two tensors of general rank.
 *
 * @tparam rank The rank of both tensors.
 *
 * @relatesalso Tensor
 */
template <int rank, int dim, typename Number, typename OtherNumber>
inline DEAL_II_ALWAYS_INLINE
Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type>
operator+ (const Tensor<rank,dim,Number> &p,
           const Tensor<rank,dim,OtherNumber> &q)
{
  Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type> tmp (p);

  for (unsigned int i=0; i<dim; ++i)
    tmp[i] += q[i];

  return tmp;
}


/**
 * Subtraction of two tensors of general rank.
 *
 * @tparam rank The rank of both tensors.
 *
 * @relatesalso Tensor
 */
template <int rank, int dim, typename Number, typename OtherNumber>
inline DEAL_II_ALWAYS_INLINE
Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type>
operator- (const Tensor<rank,dim,Number> &p,
           const Tensor<rank,dim,OtherNumber> &q)
{
  Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type> tmp (p);

  for (unsigned int i=0; i<dim; ++i)
    tmp[i] -= q[i];

  return tmp;
}


//@}
/**
 * @name Contraction operations and the outer product for tensor objects
 */
//@{


/**
 * The dot product (single contraction) for tensors: Return a tensor of rank
 * $(\text{rank}_1 + \text{rank}_2 - 2)$ that is the contraction of the last
 * index of a tensor @p src1 of rank @p rank_1 with the first index of a
 * tensor @p src2 of rank @p rank_2:
 * @f[
 *   \text{result}_{i_1,\ldots,i_{r1},j_1,\ldots,j_{r2}}
 *   = \sum_{k}
 *     \text{left}_{i_1,\ldots,i_{r1}, k}
 *     \text{right}_{k, j_1,\ldots,j_{r2}}
 * @f]
 *
 * @note For the Tensor class, the multiplication operator only performs a
 * contraction over a single pair of indices. This is in contrast to the
 * multiplication operator for SymmetricTensor, which does the double
 * contraction.
 *
 * @note In case the contraction yields a tensor of rank 0 the scalar number
 * is returned as an unwrapped number type.
 *
 * @relatesalso Tensor
 * @author Matthias Maier, 2015
 */
template <int rank_1, int rank_2, int dim,
          typename Number, typename OtherNumber>
inline DEAL_II_ALWAYS_INLINE
typename Tensor<rank_1 + rank_2 - 2, dim, typename ProductType<Number, OtherNumber>::type>::tensor_type
operator * (const Tensor<rank_1, dim, Number> &src1,
            const Tensor<rank_2, dim, OtherNumber> &src2)
{
  typename Tensor<rank_1 + rank_2 - 2, dim, typename ProductType<Number, OtherNumber>::type>::tensor_type result;

  TensorAccessors::internal::ReorderedIndexView<0, rank_2, const Tensor<rank_2, dim, OtherNumber> >
  reordered = TensorAccessors::reordered_index_view<0, rank_2>(src2);
  TensorAccessors::contract<1, rank_1, rank_2, dim>(result, src1, reordered);

  return result;
}


/**
 * Generic contraction of a pair of indices of two tensors of arbitrary rank:
 * Return a tensor of rank $(\text{rank}_1 + \text{rank}_2 - 2)$ that is the
 * contraction of index @p index_1 of a tensor @p src1 of rank @p rank_1 with
 * the index @p index_2 of a tensor @p src2 of rank @p rank_2:
 * @f[
 *   \text{result}_{i_1,\ldots,i_{r1},j_1,\ldots,j_{r2}}
 *   = \sum_{k}
 *     \text{left}_{i_1,\ldots,k,\ldots,i_{r1}}
 *     \text{right}_{j_1,\ldots,k,\ldots,j_{r2}}
 * @f]
 *
 * If for example the first index (<code>index_1==0</code>) of a tensor
 * <code>t1</code> shall be contracted with the third index
 * (<code>index_2==2</code>) of a tensor <code>t2</code>, the invocation of
 * this function is
 * @code
 *   contract<0, 2>(t1, t2);
 * @endcode
 *
 * @note The position of the index is counted from 0, i.e.,
 * $0\le\text{index}_i<\text{range}_i$.
 *
 * @note In case the contraction yields a tensor of rank 0 the scalar number
 * is returned as an unwrapped number type.
 *
 * @relatesalso Tensor
 * @author Matthias Maier, 2015
 */
template <int index_1, int index_2,
          int rank_1, int rank_2, int dim,
          typename Number, typename OtherNumber>
inline DEAL_II_ALWAYS_INLINE
typename Tensor<rank_1 + rank_2 - 2, dim, typename ProductType<Number, OtherNumber>::type>::tensor_type
contract (const Tensor<rank_1, dim, Number> &src1,
          const Tensor<rank_2, dim, OtherNumber> &src2)
{
  Assert(0 <= index_1 && index_1 < rank_1,
         ExcMessage("The specified index_1 must lie within the range [0,rank_1)"));
  Assert(0 <= index_2 && index_2 < rank_2,
         ExcMessage("The specified index_2 must lie within the range [0,rank_2)"));

  using namespace TensorAccessors;
  using namespace TensorAccessors::internal;

  // Reorder index_1 to the end of src1:
  ReorderedIndexView<index_1, rank_1, const Tensor<rank_1, dim, Number> >
  reord_01 = reordered_index_view<index_1, rank_1>(src1);

  // Reorder index_2 to the end of src2:
  ReorderedIndexView<index_2, rank_2, const Tensor<rank_2, dim, OtherNumber> >
  reord_02 = reordered_index_view<index_2, rank_2>(src2);

  typename Tensor<rank_1 + rank_2 - 2, dim, typename ProductType<Number, OtherNumber>::type>::tensor_type
  result;
  TensorAccessors::contract<1, rank_1, rank_2, dim>(result, reord_01, reord_02);
  return result;
}


/**
 * Generic contraction of two pairs of indices of two tensors of arbitrary
 * rank: Return a tensor of rank $(\text{rank}_1 + \text{rank}_2 - 4)$ that is
 * the contraction of index @p index_1 with index @p index_2, and index @p
 * index_3 with index @p index_4 of a tensor @p src1 of rank @p rank_1 and a
 * tensor @p src2 of rank @p rank_2:
 * @f[
 *   \text{result}_{i_1,\ldots,i_{r1},j_1,\ldots,j_{r2}}
 *   = \sum_{k, l}
 *     \text{left}_{i_1,\ldots,k,\ldots,l,\ldots,i_{r1}}
 *     \text{right}_{j_1,\ldots,k,\ldots,l\ldots,j_{r2}}
 * @f]
 *
 * If for example the first index (<code>index_1==0</code>) shall be
 * contracted with the third index (<code>index_2==2</code>), and the second
 * index (<code>index_3==1</code>) with the first index
 * (<code>index_4==0</code>) the invocation of this function is this function
 * is
 * @code
 *   contract<0, 2, 1, 0>(t1, t2);
 * @endcode
 *
 * @note The position of the index is counted from 0, i.e.,
 * $0\le\text{index}_i<\text{range}_i$.
 *
 * @note In case the contraction yields a tensor of rank 0 the scalar number
 * is returned as an unwrapped number type.
 *
 * @relatesalso Tensor
 * @author Matthias Maier, 2015
 */
template <int index_1, int index_2, int index_3, int index_4,
          int rank_1, int rank_2, int dim,
          typename Number, typename OtherNumber>
inline
typename Tensor<rank_1 + rank_2 - 4, dim, typename ProductType<Number, OtherNumber>::type>::tensor_type
double_contract (const Tensor<rank_1, dim, Number> &src1,
                 const Tensor<rank_2, dim, OtherNumber> &src2)
{
  Assert(0 <= index_1 && index_1 < rank_1,
         ExcMessage("The specified index_1 must lie within the range [0,rank_1)"));
  Assert(0 <= index_3 && index_3 < rank_1,
         ExcMessage("The specified index_3 must lie within the range [0,rank_1)"));
  Assert(index_1 != index_3,
         ExcMessage("index_1 and index_3 must not be the same"));
  Assert(0 <= index_2 && index_2 < rank_2,
         ExcMessage("The specified index_2 must lie within the range [0,rank_2)"));
  Assert(0 <= index_4 && index_4 < rank_2,
         ExcMessage("The specified index_4 must lie within the range [0,rank_2)"));
  Assert(index_2 != index_4,
         ExcMessage("index_2 and index_4 must not be the same"));

  using namespace TensorAccessors;
  using namespace TensorAccessors::internal;

  // Reorder index_1 to the end of src1:
  ReorderedIndexView<index_1, rank_1, const Tensor<rank_1, dim, Number> >
  reord_1 = TensorAccessors::reordered_index_view<index_1, rank_1>(src1);

  // Reorder index_2 to the end of src2:
  ReorderedIndexView<index_2, rank_2, const Tensor<rank_2, dim, OtherNumber> >
  reord_2 = TensorAccessors::reordered_index_view<index_2, rank_2>(src2);

  // Now, reorder index_3 to the end of src1. We have to make sure to
  // preserve the orginial ordering: index_1 has been removed. If
  // index_3 > index_1, we have to use (index_3 - 1) instead:
  ReorderedIndexView<(index_3 < index_1 ? index_3 : index_3 - 1), rank_1, ReorderedIndexView<index_1, rank_1, const Tensor<rank_1, dim, Number> > >
  reord_3 = TensorAccessors::reordered_index_view<index_3 < index_1 ? index_3 : index_3 - 1, rank_1>(reord_1);

  // Now, reorder index_4 to the end of src2. We have to make sure to
  // preserve the orginial ordering: index_2 has been removed. If
  // index_4 > index_2, we have to use (index_4 - 1) instead:
  ReorderedIndexView<(index_4 < index_2 ? index_4 : index_4 - 1), rank_2, ReorderedIndexView<index_2, rank_2, const Tensor<rank_2, dim, OtherNumber> > >
  reord_4 = TensorAccessors::reordered_index_view<index_4 < index_2 ? index_4 : index_4 - 1, rank_2>(reord_2);

  typename Tensor<rank_1 + rank_2 - 4, dim, typename ProductType<Number, OtherNumber>::type>::tensor_type
  result;
  TensorAccessors::contract<2, rank_1, rank_2, dim>(result, reord_3, reord_4);
  return result;
}


/**
 * The scalar product, or (generalized) Frobenius inner product of two tensors
 * of equal rank: Return a scalar number that is the result of a full
 * contraction of a tensor @p left and @p right:
 * @f[
 *   \sum_{i_1,\ldots,i_r}
 *   \text{left}_{i_1,\ldots,i_r}
 *   \text{right}_{i_1,\ldots,i_r}
 * @f]
 *
 * @relatesalso Tensor
 * @author Matthias Maier, 2015
 */
template <int rank, int dim, typename Number, typename OtherNumber>
inline DEAL_II_ALWAYS_INLINE
typename ProductType<Number, OtherNumber>::type
scalar_product (const Tensor<rank, dim, Number> &left,
                const Tensor<rank, dim, OtherNumber> &right)
{
  typename ProductType<Number, OtherNumber>::type result;
  TensorAccessors::contract<rank, rank, rank, dim>(result, left, right);
  return result;
}


/**
 * Full contraction of three tensors: Return a scalar number that is the
 * result of a full contraction of a tensor @p left of rank @p rank_1, a
 * tensor @p middle of rank $(\text{rank}_1+\text{rank}_2)$ and a tensor @p
 * right of rank @p rank_2:
 * @f[
 *   \sum_{i_1,\ldots,i_{r1},j_1,\ldots,j_{r2}}
 *   \text{left}_{i_1,\ldots,i_{r1}}
 *   \text{middle}_{i_1,\ldots,i_{r1},j_1,\ldots,j_{r2}}
 *   \text{right}_{j_1,\ldots,j_{r2}}
 * @f]
 *
 * @note Each of the three input tensors can be either a Tensor or
 * SymmetricTensor.
 *
 * @relatesalso Tensor
 * @author Matthias Maier, 2015, Jean-Paul Pelteret 2017
 */
template <template <int, int, typename> class TensorT1,
          template <int, int, typename> class TensorT2,
          template <int, int, typename> class TensorT3,
          int rank_1, int rank_2, int dim,
          typename T1, typename T2, typename T3>
typename ProductType<T1, typename ProductType<T2, T3>::type>::type
contract3 (const TensorT1<rank_1, dim, T1>          &left,
           const TensorT2<rank_1 + rank_2, dim, T2> &middle,
           const TensorT3<rank_2, dim, T3>          &right)
{
  typedef typename ProductType<T1, typename ProductType<T2, T3>::type>::type
  return_type;
  return TensorAccessors::contract3<rank_1, rank_2, dim, return_type>(
           left, middle, right);
}


/**
 * The outer product of two tensors of @p rank_1 and @p rank_2: Returns a
 * tensor of rank $(\text{rank}_1 + \text{rank}_2)$:
 * @f[
 *   \text{result}_{i_1,\ldots,i_{r1},j_1,\ldots,j_{r2}}
 *   = \text{left}_{i_1,\ldots,i_{r1}}\,\text{right}_{j_1,\ldots,j_{r2}.}
 * @f]
 *
 * @relatesalso Tensor
 * @author Matthias Maier, 2015
 */
template <int rank_1, int rank_2, int dim,
          typename Number, typename OtherNumber>
inline
Tensor<rank_1 + rank_2, dim, typename ProductType<Number, OtherNumber>::type>
outer_product(const Tensor<rank_1, dim, Number> &src1,
              const Tensor<rank_2, dim, OtherNumber> &src2)
{
  typename Tensor<rank_1 + rank_2, dim, typename ProductType<Number, OtherNumber>::type>::tensor_type result;
  TensorAccessors::contract<0, rank_1, rank_2, dim>(result, src1, src2);
  return result;
}


//@}
/**
 * @name Special operations on tensors of rank 1
 */
//@{


/**
 * Return the cross product in 2d. This is just a rotation by 90 degrees
 * clockwise to compute the outer normal from a tangential vector. This
 * function is defined for all space dimensions to allow for dimension
 * independent programming (e.g. within switches over the space dimension),
 * but may only be called if the actual dimension of the arguments is two
 * (e.g. from the <tt>dim==2</tt> case in the switch).
 *
 * @relatesalso Tensor
 * @author Guido Kanschat, 2001
 */
template <int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE
Tensor<1,dim,Number>
cross_product_2d (const Tensor<1,dim,Number> &src)
{
  Assert (dim==2, ExcInternalError());

  Tensor<1, dim, Number> result;

  result[0] = src[1];
  result[1] = -src[0];

  return result;
}


/**
 * Return the cross product of 2 vectors in 3d. This function is defined for
 * all space dimensions to allow for dimension independent programming (e.g.
 * within switches over the space dimension), but may only be called if the
 * actual dimension of the arguments is three (e.g. from the <tt>dim==3</tt>
 * case in the switch).
 *
 * @relatesalso Tensor
 * @author Guido Kanschat, 2001
 */
template <int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE
Tensor<1,dim,Number>
cross_product_3d (const Tensor<1,dim,Number> &src1,
                  const Tensor<1,dim,Number> &src2)
{
  Assert (dim==3, ExcInternalError());

  Tensor<1, dim, Number> result;

  result[0] = src1[1]*src2[2] - src1[2]*src2[1];
  result[1] = src1[2]*src2[0] - src1[0]*src2[2];
  result[2] = src1[0]*src2[1] - src1[1]*src2[0];

  return result;
}


//@}
/**
 * @name Special operations on tensors of rank 2
 */
//@{


/**
 * Compute the determinant of a tensor or rank 2.
 *
 * @relatesalso Tensor
 * @author Wolfgang Bangerth, 2009
 */
template <int dim, typename Number>
inline
Number determinant (const Tensor<2,dim,Number> &t)
{
  // Compute the determinant using the Laplace expansion of the
  // determinant. We expand along the last row.
  Number det = internal::NumberType<Number>::value(0.0);

  for (unsigned int k=0; k<dim; ++k)
    {
      Tensor<2,dim-1,Number> minor;
      for (unsigned int i=0; i<dim-1; ++i)
        for (unsigned int j=0; j<dim-1; ++j)
          minor[i][j] = t[i][j<k ? j : j+1];

      const Number cofactor = ((k % 2 == 0) ? -1. : 1.) * determinant(minor);

      det += t[dim-1][k] * cofactor;
    }

  return ((dim % 2 == 0) ? 1. : -1.) * det;
}

/**
 * Specialization for dim==1.
 *
 * @relatesalso Tensor
 */
template <typename Number>
inline
Number determinant (const Tensor<2,1,Number> &t)
{
  return t[0][0];
}


/**
 * Compute and return the trace of a tensor of rank 2, i.e. the sum of its
 * diagonal entries.
 *
 * @relatesalso Tensor
 * @author Wolfgang Bangerth, 2001
 */
template <int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE
Number trace (const Tensor<2,dim,Number> &d)
{
  Number t=d[0][0];
  for (unsigned int i=1; i<dim; ++i)
    t += d[i][i];
  return t;
}


/**
 * Compute and return the inverse of the given tensor. Since the compiler can
 * perform the return value optimization, and since the size of the return
 * object is known, it is acceptable to return the result by value, rather
 * than by reference as a parameter.
 *
 * @relatesalso Tensor
 * @author Wolfgang Bangerth, 2000
 */
template <int dim, typename Number>
inline
Tensor<2,dim,Number>
invert (const Tensor<2,dim,Number> &)
{
  Number return_tensor [dim][dim];

  // if desired, take over the
  // inversion of a 4x4 tensor
  // from the FullMatrix
  AssertThrow (false, ExcNotImplemented());

  return Tensor<2,dim,Number>(return_tensor);
}


#ifndef DOXYGEN

template <typename Number>
inline
Tensor<2,1,Number>
invert (const Tensor<2,1,Number> &t)
{
  Number return_tensor [1][1];

  return_tensor[0][0] = 1.0/t[0][0];

  return Tensor<2,1,Number>(return_tensor);
}


template <typename Number>
inline
Tensor<2,2,Number>
invert (const Tensor<2,2,Number> &t)
{
  Tensor<2,2,Number> return_tensor;

  // this is Maple output,
  // thus a bit unstructured
  const Number inv_det_t = 1.0/(t[0][0]*t[1][1]-t[1][0]*t[0][1]);
  return_tensor[0][0] = t[1][1];
  return_tensor[0][1] = -t[0][1];
  return_tensor[1][0] = -t[1][0];
  return_tensor[1][1] = t[0][0];
  return_tensor *= inv_det_t;

  return return_tensor;
}


template <typename Number>
inline
Tensor<2,3,Number>
invert (const Tensor<2,3,Number> &t)
{
  Tensor<2,3,Number> return_tensor;

  const Number t4 = t[0][0]*t[1][1],
               t6 = t[0][0]*t[1][2],
               t8 = t[0][1]*t[1][0],
               t00 = t[0][2]*t[1][0],
               t01 = t[0][1]*t[2][0],
               t04 = t[0][2]*t[2][0],
               inv_det_t = 1.0/(t4*t[2][2]-t6*t[2][1]-t8*t[2][2]+
                                t00*t[2][1]+t01*t[1][2]-t04*t[1][1]);
  return_tensor[0][0] = t[1][1]*t[2][2]-t[1][2]*t[2][1];
  return_tensor[0][1] = t[0][2]*t[2][1]-t[0][1]*t[2][2];
  return_tensor[0][2] = t[0][1]*t[1][2]-t[0][2]*t[1][1];
  return_tensor[1][0] = t[1][2]*t[2][0]-t[1][0]*t[2][2];
  return_tensor[1][1] = t[0][0]*t[2][2]-t04;
  return_tensor[1][2] = t00-t6;
  return_tensor[2][0] = t[1][0]*t[2][1]-t[1][1]*t[2][0];
  return_tensor[2][1] = t01-t[0][0]*t[2][1];
  return_tensor[2][2] = t4-t8;
  return_tensor *= inv_det_t;

  return return_tensor;
}

#endif /* DOXYGEN */


/**
 * Return the transpose of the given tensor.
 *
 * @relatesalso Tensor
 * @author Wolfgang Bangerth, 2002
 */
template <int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE
Tensor<2,dim,Number>
transpose (const Tensor<2,dim,Number> &t)
{
  Tensor<2, dim, Number> tt;
  for (unsigned int i=0; i<dim; ++i)
    {
      tt[i][i] = t[i][i];
      for (unsigned int j=i+1; j<dim; ++j)
        {
          tt[i][j] = t[j][i];
          tt[j][i] = t[i][j];
        };
    }
  return tt;
}


/**
 * Return the adjugate of the given tensor of rank 2.
 * The adjugate of a tensor $\left(\bullet\right)$ is defined as
 * @f[
 *  \textrm{adj}\left(\bullet\right)
 *   := \textrm{det}\left(\bullet\right) \; \left(\bullet\right)^{-1} \; .
 * @f]
 *
 * @note This requires that the tensor is invertible.
 *
 * @relatesalso Tensor
 * @author Jean-Paul Pelteret, 2016
 */
template <int dim, typename Number>
inline
Tensor<2,dim,Number>
adjugate (const Tensor<2,dim,Number> &t)
{
  return determinant(t)*invert(t);
}


/**
 * Return the cofactor of the given tensor of rank 2.
 * The cofactor of a tensor $\left(\bullet\right)$ is defined as
 * @f[
 *  \textrm{cof}\left(\bullet\right)
 *   := \textrm{det}\left(\bullet\right) \; \left(\bullet\right)^{-T}
 *    = \left[ \textrm{adj}\left(\bullet\right) \right]^{T} \; .
 * @f]
 *
 * @note This requires that the tensor is invertible.
 *
 * @relatesalso Tensor
 * @author Jean-Paul Pelteret, 2016
 */
template <int dim, typename Number>
inline
Tensor<2,dim,Number>
cofactor (const Tensor<2,dim,Number> &t)
{
  return transpose(adjugate(t));
}


/**
 * Return the $l_1$ norm of the given rank-2 tensor, where $||t||_1 = \max_j
 * \sum_i |t_{ij}|$ (maximum of the sums over columns).
 *
 * @relatesalso Tensor
 * @author Wolfgang Bangerth, 2012
 */
template <int dim, typename Number>
inline
Number
l1_norm (const Tensor<2,dim,Number> &t)
{
  Number max = internal::NumberType<Number>::value(0.0);
  for (unsigned int j=0; j<dim; ++j)
    {
      Number sum = internal::NumberType<Number>::value(0.0);
      for (unsigned int i=0; i<dim; ++i)
        sum += std::fabs(t[i][j]);

      if (sum > max)
        max = sum;
    }

  return max;
}


/**
 * Return the $l_\infty$ norm of the given rank-2 tensor, where $||t||_\infty
 * = \max_i \sum_j |t_{ij}|$ (maximum of the sums over rows).
 *
 * @relatesalso Tensor
 * @author Wolfgang Bangerth, 2012
 */
template <int dim, typename Number>
inline
Number
linfty_norm (const Tensor<2,dim,Number> &t)
{
  Number max = internal::NumberType<Number>::value(0.0);
  for (unsigned int i=0; i<dim; ++i)
    {
      Number sum = internal::NumberType<Number>::value(0.0);
      for (unsigned int j=0; j<dim; ++j)
        sum += std::fabs(t[i][j]);

      if (sum > max)
        max = sum;
    }

  return max;
}

//@}

DEAL_II_NAMESPACE_CLOSE

// include deprecated non-member functions operating on Tensor
#include <deal.II/base/tensor_deprecated.h>

#endif
