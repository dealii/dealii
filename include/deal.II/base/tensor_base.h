// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2015 by the deal.II authors
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

#ifndef dealii__tensor_base_h
#define dealii__tensor_base_h


// single this file out from tensor.h, since we want to derive
// Point<dim,Number> from Tensor<1,dim,Number>. However, the point class will
// not need all the tensor stuff, so we don't want the whole tensor package to
// be included every time we use a point.


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/template_constraints.h>
#include <vector>

#include <cmath>
#include <ostream>

DEAL_II_NAMESPACE_OPEN

template <typename number> class Vector;

// forward declare Point and Tensor. This is the first definition of these
// classes and here we set the default Number type to double (this means that
// this file must be included when using something like Tensor<1,dim>, and
// Point and Tensor must not be forward declared without the number type
// specified)
template <int dim, typename Number> class Point;

// general template; specialized for rank==1; the general template is in
// tensor.h
template <int rank_, int dim, typename Number=double> class Tensor;
template <int dim, typename Number> class Tensor<0,dim,Number>;
template <int dim, typename Number> class Tensor<1,dim,Number>;



/**
 * This class is a specialized version of the
 * <tt>Tensor<rank,dim,Number></tt> class. It handles tensors of rank zero,
 * i.e. scalars. The second template argument @param dim is ignored.
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
  static const unsigned int rank      = 0;

  /**
   * Type of stored objects. This is a Number for a rank 1 tensor.
   */

  typedef Number value_type;

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
   * Declare an array type which can be used to initialize an object of this
   * type statically. In case of a a tensor of rank 0 this is just a scalar
   * number type
   */
  typedef Number array_type;

  /**
   * Constructor. Set to zero.
   */
  Tensor ();

  /**
   * Copy constructor.
   */
  Tensor (const Tensor<0,dim,Number> &initializer);

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
  Tensor (const OtherNumber initializer);

  /**
   * Conversion to Number. Since rank-0 tensors are scalars, this is a natural
   * operation.
   */
  operator Number () const;

  /**
   * Conversion to Number. Since rank-0 tensors are scalars, this is a natural
   * operation.
   *
   * This is the non-const conversion operator that returns a writable
   * reference.
   */
  operator Number &();

  /**
   * Copy assignment operator.
   */
  Tensor<0,dim,Number> &operator = (const Tensor<0,dim,Number> &rhs);

  /**
   * Assignment from tensors with different underlying scalar type.
   * This obviously requires that the @p OtherNumber type is convertible to @p
   * Number.
   */
  template <typename OtherNumber>
  Tensor<0,dim,Number> &operator = (const Tensor<0,dim,OtherNumber> &rhs);

  /**
   * Assignment operator.
   */
  template <typename OtherNumber>
  Tensor<0,dim,Number> &operator = (const OtherNumber d);

  /**
   * Test for equality of two tensors.
   */
  template<typename OtherNumber>
  bool operator == (const Tensor<0,dim,OtherNumber> &rhs) const;

  /**
   * Test for inequality of two tensors.
   */
  template<typename OtherNumber>
  bool operator != (const Tensor<0,dim,OtherNumber> &rhs) const;

  /**
   * Add another scalar
   */
  template<typename OtherNumber>
  Tensor<0,dim,Number> &operator += (const Tensor<0,dim,OtherNumber> &rhs);

  /**
   * Subtract another scalar.
   */
  template<typename OtherNumber>
  Tensor<0,dim,Number> &operator -= (const Tensor<0,dim,OtherNumber> &rhs);

  /**
   * Multiply the scalar with a <tt>factor</tt>.
   */
  template<typename OtherNumber>
  Tensor<0,dim,Number> &operator *= (const OtherNumber factor);

  /**
   * Divide the scalar by <tt>factor</tt>.
   */
  template<typename OtherNumber>
  Tensor<0,dim,Number> &operator /= (const OtherNumber factor);

  /**
   * Tensor with inverted entries.
   */
  Tensor<0,dim,Number>   operator - () const;

  /**
   * Return the Frobenius-norm of a tensor, i.e. the square root of the sum of
   * squares of all entries. For the present case of rank-1 tensors, this
   * equals the usual <tt>l<sub>2</sub></tt> norm of the vector.
   */
  real_type norm () const;

  /**
   * Return the square of the Frobenius-norm of a tensor, i.e. the square root
   * of the sum of squares of all entries.
   *
   * This function mainly exists because it makes computing the norm simpler
   * recursively, but may also be useful in other contexts.
   */
  real_type norm_square () const;

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
   * Only tensors with a positive dimension are implemented. This exception is
   * thrown by the constructor if the template argument <tt>dim</tt> is zero
   * or less.
   *
   * @ingroup Exceptions
   */
  DeclException1 (ExcDimTooSmall,
                  int,
                  << "dim must be positive, but was " << arg1);

  /**
   * Read or write the data of this object to or from a stream for the purpose
   * of serialization
   */
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version);

private:
  /**
   * The value of this scalar object.
   */
  Number value;
};



/**
 * This class is a specialized version of the <tt>Tensor<rank,dim,Number></tt>
 * class.  It handles tensors with one index, i.e. vectors, of fixed dimension
 * and provides the basis for the functionality needed for tensors of higher
 * rank.
 *
 * Within deal.II, the distinction between this class and its derived class
 * <tt>Point</tt> is that we use the <tt>Point</tt> class mainly to denote the
 * points that make up geometric objects. As such, they have a small number of
 * additional operations over general tensors of rank 1 for which we use the
 * <tt>Tensor<1,dim,Number></tt> class. In particular, there is a distance()
 * function to compute the Euclidean distance between two points in space.
 *
 * However, the <tt>Point</tt> class is really only used where the coordinates
 * of an object can be thought to possess the dimension of a length. For all
 * other uses, such as the gradient of a scalar function (which is a tensor of
 * rank 1, or vector, with as many elements as a point object, but with
 * different physical units), we use the <tt>Tensor<1,dim,Number></tt> class.
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * this tensor operates. This of course equals the number of coordinates that
 * identify a point and rank-1 tensor.
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
template <int dim,typename Number>
class Tensor<1,dim,Number>
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
  static const unsigned int rank = 1;

  /**
   * Number of independent components of a tensor of rank 1.
   */
  static const unsigned int
  n_independent_components = dim;

  /**
   * Type of stored objects. This is a Number for a rank 1 tensor.
   */

  typedef Number value_type;

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
   * Declare an array type which can be used to initialize statically an
   * object of this type.
   */
  // Avoid a bogus warning in case of dim==0, and always provide a type
  // with positive array size. The constructor will take care that no
  // Tensor with dim==0 will be constructed.
  typedef Number array_type[(dim!=0) ? dim : 1];

  /**
   * Constructor. Initialize all entries to zero if <tt>initialize==true</tt>;
   * this is the default behaviour.
   */
  explicit
  Tensor (const bool initialize = true);

  /**
   * Copy constructor.
   */
  Tensor (const Tensor<1,dim,Number> &initializer);

  /**
   * Copy constructor, where the data is copied from a C-style array.
   */
  Tensor (const array_type &initializer);

  /**
   * Copy constructor from tensors with different underlying scalar type. This
   * obviously requires that the @p OtherNumber type is convertible to @p
   * Number.
   */
  template <typename OtherNumber>
  Tensor (const Tensor<1,dim,OtherNumber> &initializer);

  /**
   * Read access to the <tt>index</tt>th coordinate.
   *
   * Note that the derived <tt>Point</tt> class also provides access through
   * the <tt>()</tt> operator for backcompatibility.
   */
  Number operator [] (const unsigned int index) const;

  /**
   * Read and write access to the <tt>index</tt>th coordinate.
   *
   * Note that the derived <tt>Point</tt> class also provides access through
   * the <tt>()</tt> operator for backcompatibility.
   */
  Number &operator [] (const unsigned int index);

  /**
   * Read access using TableIndices <tt>indices</tt>
   */
  Number operator [] (const TableIndices<1> &indices) const;

  /**
   * Read and write access using TableIndices <tt>indices</tt>
   */
  Number &operator [] (const TableIndices<1> &indices);

  /**
   * Copy assignment operator.
   */
  Tensor<1,dim,Number> &operator = (const Tensor<1,dim,Number> &rhs);

  /**
   * Assignment operator from tensors with different underlying scalar type.
   * This obviously requires that the @p OtherNumber type is convertible to @p
   * Number.
   */
  template <typename OtherNumber>
  Tensor<1,dim,Number> &operator = (const Tensor<1,dim,OtherNumber> &rhs);

  /**
   * This operator assigns a scalar to a tensor. To avoid confusion with what
   * exactly it means to assign a scalar value to a tensor, zero is the only
   * value allowed for <tt>d</tt>, allowing the intuitive notation
   * <tt>t=0</tt> to reset all elements of the tensor to zero.
   */
  template <typename OtherNumber>
  Tensor<1,dim,Number> &operator = (const OtherNumber d);

  /**
   * Test for equality of two tensors.
   */
  template <typename OtherNumber>
  bool operator == (const Tensor<1,dim,OtherNumber> &rhs) const;

  /**
   * Test for inequality of two tensors.
   */
  template <typename OtherNumber>
  bool operator != (const Tensor<1,dim,OtherNumber> &rhs) const;

  /**
   * Add another vector to this vector.
   */
  template <typename OtherNumber>
  Tensor<1,dim,Number> &operator += (const Tensor<1,dim,OtherNumber> &rhs);

  /**
   * Subtract another vector.
   */
  template <typename OtherNumber>
  Tensor<1,dim,Number> &operator -= (const Tensor<1,dim,OtherNumber> &rhs);

  /**
   * Scale the vector by <tt>factor</tt>, i.e., multiply all coordinates by
   * <tt>factor</tt>.
   */
  template <typename OtherNumber>
  Tensor<1,dim,Number> &operator *= (const OtherNumber factor);

  /**
   * Scale the vector by <tt>1/factor</tt>.
   */
  template <typename OtherNumber>
  Tensor<1,dim,Number> &operator /= (const OtherNumber factor);

  /**
   * Tensor with inverted entries.
   */
  Tensor<1,dim,Number> operator - () const;

  /**
   * Return the Frobenius-norm of a tensor, i.e. the square root of the sum of
   * squares of all entries. For the present case of rank-1 tensors, this
   * equals the usual <tt>l<sub>2</sub></tt> norm of the vector.
   */
  real_type norm () const;

  /**
   * Return the square of the Frobenius-norm of a tensor, i.e. the square root
   * of the sum of squares of all entries.
   *
   * This function mainly exists because it makes computing the norm simpler
   * recursively, but may also be useful in other contexts.
   */
  real_type norm_square () const;

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
   * Fill a vector with all tensor elements.
   *
   * This function unrolls all tensor entries into a single, linearly numbered
   * vector. As usual in C++, the rightmost index marches fastest.
   */
  template <typename Number2>
  void unroll (Vector<Number2> &result) const;

  /**
   * Returns an unrolled index in the range [0,dim-1] for the element of the
   * tensor indexed by the argument to the function.
   *
   * Given that this is a rank-1 object, the returned value is simply the
   * value of the only index stored by the argument.
   */
  static
  unsigned int
  component_to_unrolled_index(const TableIndices<1> &indices);

  /**
   * Opposite of  component_to_unrolled_index: For an index in the range
   * [0,dim-1], return which set of indices it would correspond to.
   *
   * Given that this is a rank-1 object, the returned set of indices consists
   * of only a single element with value equal to the argument to this
   * function.
   */
  static
  TableIndices<1> unrolled_to_component_indices(const unsigned int i);


  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  static std::size_t memory_consumption ();

  /**
   * Only tensors with a positive dimension are implemented. This exception is
   * thrown by the constructor if the template argument <tt>dim</tt> is zero
   * or less.
   *
   * @ingroup Exceptions
   */
  DeclException1 (ExcDimTooSmall,
                  int,
                  << "dim must be positive, but was " << arg1);

  /**
   * Read or write the data of this object to or from a stream for the purpose
   * of serialization
   */
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version);

private:
  /**
   * Store the values in a simple array.  For <tt>dim==0</tt> store one
   * element, because otherwise the compiler would choke.  We catch this case
   * in the constructor to disallow the creation of such an object.
   */
  array_type values;

  /**
   * Help function for unroll. If we have detected an access control bug in
   * the compiler, this function is declared public, otherwise private. Do not
   * attempt to use this function from outside in any case, even if it should
   * be public for your compiler.
   */
  template <typename Number2>
  void unroll_recursion (Vector<Number2> &result,
                         unsigned int    &start_index) const;

  /**
   * Make the following classes friends to this class. In principle, it would
   * suffice if otherrank==2, but that is not possible in C++ at present.
   *
   * Also, it would be sufficient to make the function unroll_loops a friend,
   * but that seems to be impossible as well.
   */
  template <int otherrank, int otherdim, typename OtherNumber>  friend class dealii::Tensor;

  /**
   * Point is allowed access to the coordinates. This is supposed to improve
   * speed.
   */
  friend class Point<dim,Number>;
};


/**
 * Prints the value of this scalar.
 */
template <int dim,typename Number>
std::ostream &operator << (std::ostream &out, const Tensor<0,dim,Number> &p);

/**
 * Prints the values of this tensor in the form <tt>x1 x2 x3 etc</tt>.
 */
template <int dim,typename Number>
std::ostream &operator << (std::ostream &out, const Tensor<1,dim,Number> &p);


#ifndef DOXYGEN



/*---------------------- Inline functions: Tensor<0,dim> ---------------------*/



template <int dim,typename Number>
inline
Tensor<0,dim,Number>::Tensor ()
{
  Assert (dim>0, ExcDimTooSmall(dim));

  value = value_type();
}



template <int dim, typename Number>
inline
Tensor<0,dim,Number>::Tensor (const Tensor<0,dim,Number> &p)
{
  Assert (dim>0, ExcDimTooSmall(dim));

  value = p.value;
}



template <int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<0,dim,Number>::Tensor (const OtherNumber initializer)
{
  Assert (dim>0, ExcDimTooSmall(dim));

  value = initializer;
}



template <int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<0,dim,Number>::Tensor (const Tensor<0,dim,OtherNumber> &p)
{
  Assert (dim>0, ExcDimTooSmall(dim));

  value = p.value;
}



template <int dim, typename Number>
inline
Tensor<0,dim,Number>::operator Number () const
{
  return value;
}



template <int dim, typename Number>
inline
Tensor<0,dim,Number>::operator Number &()
{
  return value;
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
Tensor<0,dim,Number> &Tensor<0,dim,Number>::operator = (const Tensor<0,dim,OtherNumber> &p)
{
  value = p.value;
  return *this;
}



template <int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<0,dim,Number> &Tensor<0,dim,Number>::operator = (const OtherNumber d)
{
  value = d;
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
Tensor<0,dim,Number> &Tensor<0,dim,Number>::operator *= (const OtherNumber s)
{
  value *= s;
  return *this;
}



template <int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<0,dim,Number> &Tensor<0,dim,Number>::operator /= (const OtherNumber s)
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
  return numbers::NumberTraits<Number>::abs (value);
}



template <int dim, typename Number>
inline
typename Tensor<0,dim,Number>::real_type
Tensor<0,dim,Number>::norm_square () const
{
  return numbers::NumberTraits<Number>::abs_square (value);
}



template <int dim, typename Number>
inline
void Tensor<0,dim,Number>::clear ()
{
  value = value_type();
}



template <int dim, typename Number>
template <class Archive>
inline
void Tensor<0,dim,Number>::serialize(Archive &ar, const unsigned int)
{
  ar &value;
}



#ifndef DEAL_II_WITH_CXX11

template <typename T, typename U, int rank, int dim>
struct ProductType<T,Tensor<rank,dim,U> >
{
  typedef Tensor<rank,dim,typename ProductType<T,U>::type> type;
};

template <typename T, typename U, int rank, int dim>
struct ProductType<Tensor<rank,dim,T>,U>
{
  typedef Tensor<rank,dim,typename ProductType<T,U>::type> type;
};

#endif



/**
 * TODO
 *
 * @relates Tensor
 * @relates EnableIfScalar
 */
template <int dim,
         typename Number,
         typename OtherNumber,
         typename = typename EnableIfScalar<OtherNumber>::type>
inline
Tensor<0,dim,typename ProductType<OtherNumber, Number>::type>
operator * (const OtherNumber           factor,
            const Tensor<0,dim,Number> &t)
{
  return factor * static_cast<Number>(t);
}



/**
 * TODO
 *
 * @relates Tensor
 * @relates EnableIfScalar
 */
template <int dim,
         typename Number,
         typename OtherNumber,
         typename = typename EnableIfScalar<OtherNumber>::type>
inline
Tensor<0,dim,typename ProductType<Number, OtherNumber>::type>
operator * (const Tensor<0,dim,Number> &t,
            const OtherNumber           factor)
{
  return static_cast<Number>(t) * factor;
}



/**
 * TODO
 *
 * @relates Tensor
 * @relates EnableIfScalar
 */
template <int dim,
         typename Number,
         typename OtherNumber,
         typename = typename EnableIfScalar<OtherNumber>::type>
inline
Tensor<0,dim,typename ProductType<Number, OtherNumber>::type>
operator / (const Tensor<0,dim,Number> &t,
            const OtherNumber           factor)
{
  return static_cast<Number>(t) / factor;
}



/**
 * Add two tensors of rank 0.
 *
 * @relates Tensor
 */
template <int dim, typename Number, typename OtherNumber>
inline
Tensor<0, dim, typename ProductType<Number, OtherNumber>::type>
operator+ (const Tensor<0,dim,Number> &p, const Tensor<0,dim,OtherNumber> &q)
{
  return static_cast<Number>(p) + static_cast<OtherNumber>(q);
}



/**
 * Subtract two tensors of rank 0.
 *
 * @relates Tensor
 */
template <int dim, typename Number, typename OtherNumber>
inline
Tensor<0, dim, typename ProductType<Number, OtherNumber>::type>
operator- (const Tensor<0,dim,Number> &p, const Tensor<0,dim,OtherNumber> &q)
{
  return static_cast<Number>(p) - static_cast<OtherNumber>(q);
}



/**
 * Returns the contraction of two Tensors of rank 0.
 *
 * @relates Tensor
 */
template <int dim, typename Number, typename OtherNumber>
inline
typename ProductType<Number, OtherNumber>::type
operator* (const Tensor<0,dim,Number> &p, const Tensor<0,dim,OtherNumber> &q)
{
  return static_cast<Number>(p) * static_cast<OtherNumber>(q);
}



/*---------------------- Inline functions: Tensor<1,dim> ---------------------*/



template <int dim, typename Number>
inline
Tensor<1,dim,Number>::Tensor (const bool initialize)
{
  if (initialize)
    // need to create an object Number() to initialize to zero to avoid
    // confusion with Tensor::operator=(scalar) when using something like
    // Tensor<1,dim,Tensor<1,dim,Number> >.
    for (unsigned int i=0; i!=dim; ++i)
      values[i] = Number();
}



template <int dim, typename Number>
inline
Tensor<1,dim,Number>::Tensor (const array_type &initializer)
{
  Assert (dim>0, ExcDimTooSmall(dim));

  for (unsigned int i=0; i<dim; ++i)
    values[i] = initializer[i];
}



template <int dim, typename Number>
inline
Tensor<1,dim,Number>::Tensor (const Tensor<1,dim,Number> &p)
{
  Assert (dim>0, ExcDimTooSmall(dim));

  for (unsigned int i=0; i<dim; ++i)
    values[i] = p.values[i];
}



template <int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<1,dim,Number>::Tensor (const Tensor<1,dim,OtherNumber> &p)
{
  Assert (dim>0, ExcDimTooSmall(dim));

  for (unsigned int i=0; i<dim; ++i)
    values[i] = Number(p.values[i]);
}



// At some places in the library, we have Point<0> for formal reasons
// (e.g., we sometimes have Quadrature<dim-1> for faces, so we have
// Quadrature<0> for dim=1, and then we have Point<0>). To avoid warnings
// in the above function that the loop end check always fails, we
// implement this function here
template <>
inline
Tensor<1,0,double>::Tensor (const Tensor<1,0,double> &)
{
}



template <int dim, typename Number>
inline
Number Tensor<1,dim,Number>::operator [] (const unsigned int index) const
{
  Assert (index<dim, ExcIndexRange (index, 0, dim));
  return values[index];
}



template <int dim, typename Number>
inline
Number &Tensor<1,dim,Number>::operator [] (const unsigned int index)
{
  Assert (index<dim, ExcIndexRange (index, 0, dim));
  return values[index];
}



template <int dim, typename Number>
inline
Number Tensor<1,dim,Number>::operator [] (const TableIndices<1> &indices) const
{
  Assert (indices[0]<dim, ExcIndexRange (indices[0], 0, dim));
  return values[indices[0]];
}



template <int dim, typename Number>
inline
Number &Tensor<1,dim,Number>::operator [] (const TableIndices<1> &indices)
{
  Assert (indices[0]<dim, ExcIndexRange (indices[0], 0, dim));
  return values[indices[0]];
}



template <int dim, typename Number>
inline
Tensor<1,dim,Number> &
Tensor<1,dim,Number>::operator = (const Tensor<1,dim,Number> &p)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] = p.values[i];

  return *this;
}



template <>
inline
Tensor<1,0,double> &Tensor<1,0,double>::operator = (const Tensor<1,0,double> &)
{
  // at some places in the library, we have Point<0> for formal reasons
  // (e.g., we sometimes have Quadrature<dim-1> for faces, so we have
  // Quadrature<0> for dim=1, and then we have Point<0>). To avoid warnings
  // in the above function that the loop end check always fails, we
  // implement this function here
  return *this;
}



template <int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<1,dim,Number> &
Tensor<1,dim,Number>::operator = (const Tensor<1,dim,OtherNumber> &p)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] = Number(p.values[i]);

  return *this;
}



template <int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<1,dim,Number> &Tensor<1,dim,Number>::operator = (const OtherNumber d)
{
  Assert (d == OtherNumber(), ExcMessage ("Only assignment with zero is allowed"));
  (void) d;

  for (unsigned int i=0; i<dim; ++i)
    values[i] = Number();

  return *this;
}



template <int dim, typename Number>
template <typename OtherNumber>
inline
bool Tensor<1,dim,Number>::operator == (const Tensor<1,dim,OtherNumber> &p) const
{
  for (unsigned int i=0; i<dim; ++i)
    if (values[i] != p.values[i])
      return false;
  return true;
}



template <>
template <>
inline
bool Tensor<1,0,double>::operator == (const Tensor<1,0,double> &) const
{
  return true;
}



template <int dim, typename Number>
template <typename OtherNumber>
inline
bool Tensor<1,dim,Number>::operator != (const Tensor<1,dim,OtherNumber> &p) const
{
  return !((*this) == p);
}



template <int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<1,dim,Number> &Tensor<1,dim,Number>::operator += (const Tensor<1,dim,OtherNumber> &p)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] += p.values[i];
  return *this;
}



template <int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<1,dim,Number> &Tensor<1,dim,Number>::operator -= (const Tensor<1,dim,OtherNumber> &p)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] -= p.values[i];
  return *this;
}



template <int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<1,dim,Number> &Tensor<1,dim,Number>::operator *= (const OtherNumber s)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] *= s;
  return *this;
}



template <int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<1,dim,Number> &Tensor<1,dim,Number>::operator /= (const OtherNumber s)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] /= s;
  return *this;
}



template <int dim, typename Number>
inline
Tensor<1,dim,Number> Tensor<1,dim,Number>::operator - () const
{
  Tensor<1,dim,Number> result (false);
  for (unsigned int i=0; i<dim; ++i)
    result.values[i] = -values[i];
  return result;
}



template <int dim, typename Number>
inline
typename Tensor<1,dim,Number>::real_type
Tensor<1,dim,Number>::norm () const
{
  return std::sqrt (norm_square());
}



template <int dim, typename Number>
inline
typename Tensor<1,dim,Number>::real_type
Tensor<1,dim,Number>::norm_square () const
{
  real_type s = numbers::NumberTraits<Number>::abs_square(values[0]);
  for (unsigned int i=1; i<dim; ++i)
    s += numbers::NumberTraits<Number>::abs_square(values[i]);

  return s;
}



template <int dim, typename Number>
template <typename Number2>
inline
void
Tensor<1,dim,Number>::unroll (Vector<Number2> &result) const
{
  Assert (result.size()==dim,
          ExcDimensionMismatch(dim, result.size()));

  unsigned int index = 0;
  unroll_recursion (result,index);
}



template<int dim, typename Number>
template <typename Number2>
inline
void
Tensor<1,dim,Number>::unroll_recursion (Vector<Number2> &result,
                                        unsigned int    &index) const
{
  for (unsigned int i=0; i<dim; ++i)
    result(index++) = operator[](i);
}


template <int dim, typename Number>
inline
unsigned int
Tensor<1, dim, Number>::component_to_unrolled_index (const TableIndices<1> &indices)
{
  return indices[0];
}

template <int dim, typename Number>
inline
TableIndices<1>
Tensor<1, dim, Number>::unrolled_to_component_indices (const unsigned int i)
{
  return TableIndices<1>(i);
}



template <int dim, typename Number>
inline
void Tensor<1,dim,Number>::clear ()
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] = 0;
}



template <int dim, typename Number>
inline
std::size_t
Tensor<1,dim,Number>::memory_consumption ()
{
  return sizeof(Tensor<1,dim,Number>);
}



template <int dim, typename Number>
template <class Archive>
inline
void Tensor<1,dim,Number>::serialize(Archive &ar, const unsigned int)
{
  ar &values;
}
#endif // DOXYGEN



/**
 * Output operator for tensors of rank 0. Since such tensors are scalars, we
 * simply print this one value.
 *
 * @relates Tensor<0,dim,Number>
 */
template <int dim, typename Number>
inline
std::ostream &operator << (std::ostream &out, const Tensor<0,dim,Number> &p)
{
  out << static_cast<Number>(p);
  return out;
}



/**
 * Output operator for tensors of rank 1. Print the elements consecutively,
 * with a space in between.
 *
 * @relates Tensor<1,dim,Number>
 */
template <int dim, typename Number>
inline
std::ostream &operator << (std::ostream &out, const Tensor<1,dim,Number> &p)
{
  for (unsigned int i=0; i<dim-1; ++i)
    out << p[i] << ' ';
  out << p[dim-1];

  return out;
}



/**
 * Output operator for tensors of rank 1 and dimension 1. This is implemented
 * specialized from the general template in order to avoid a compiler warning
 * that the loop is empty.
 *
 * @relates Tensor<1,dim,Number>
 */
inline
std::ostream &operator << (std::ostream &out, const Tensor<1,1,double> &p)
{
  out << p[0];

  return out;
}



/**
 * Multiplication of a tensor of rank with a scalar number from the right.
 *
 * The purpose of this operator is to enable only multiplication of a tensor
 * by a scalar number (i.e., a floating point number, a complex floating point
 * number, etc.). The function is written in a way that only allows the
 * compiler to consider the function if the second argument is indeed a scalar
 * number -- in other words, @p OtherNumber will not match, for example
 * <code>std::vector@<double@></code> as the product of a tensor and a vector
 * clearly would make no sense. The mechanism by which the compiler is
 * prohibited of considering this operator for multiplication with non-scalar
 * types are explained in the documentation of the EnableIfScalar class.
 *
 * The return type of the function is chosen so that it matches the types of
 * both the tensor and the scalar argument. For example, if you multiply a
 * <code>Tensor@<1,dim,double@></code> by <code>std::complex@<double@></code>,
 * then the result will be a
 * <code>Tensor@<1,dim,std::complex@<double@>@></code>. In other words, the
 * type with which the returned tensor stores its components equals the type
 * you would get if you multiplied an individual component of the input tensor
 * by the scalar factor.
 *
 * @relates Tensor
 * @relates EnableIfScalar
 */
template <int rank, int dim,
         typename Number,
         typename OtherNumber,
         typename = typename EnableIfScalar<OtherNumber>::type>
inline
Tensor<rank,dim,typename ProductType<Number, OtherNumber>::type>
operator * (const Tensor<rank,dim,Number> &t,
            const OtherNumber              factor)
{
  // recurse over the base objects
  Tensor<rank,dim,typename ProductType<Number,OtherNumber>::type> tt;
  for (unsigned int d=0; d<dim; ++d)
    tt[d] = t[d] * factor;
  return tt;
}



/**
 * Multiplication of a tensor of general rank with a scalar number from the
 * left. See the discussion with the operator with switched arguments for more
 * information about template arguments and the return type.
 *
 * @relates Tensor
 * @relates EnableIfScalar
 */
template <int rank, int dim,
         typename Number,
         typename OtherNumber,
         typename = typename EnableIfScalar<OtherNumber>::type>
inline
Tensor<rank,dim,typename ProductType<Number, OtherNumber>::type>
operator * (const Number                        factor,
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
 * @relates Tensor
 * @relates EnableIfScalar
 */
template <int rank, int dim,
         typename Number,
         typename OtherNumber,
         typename = typename EnableIfScalar<OtherNumber>::type>
inline
Tensor<rank,dim,typename ProductType<Number, OtherNumber>::type>
operator / (const Tensor<rank,dim,Number> &t,
            const OtherNumber              factor)
{
  // recurse over the base objects
  Tensor<rank,dim,typename ProductType<Number,OtherNumber>::type> tt;
  for (unsigned int d=0; d<dim; ++d)
    tt[d] = t[d] / factor;
  return tt;
}



/**
 * Addition of two tensors of general @tparam rank.
 *
 * @relates Tensor
 */
template <int rank, int dim, typename Number, typename OtherNumber>
inline
Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type>
operator+ (const Tensor<rank,dim,Number> &p, const Tensor<rank,dim,OtherNumber> &q)
{
  Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type> tmp (p);

  for (unsigned int i=0; i<dim; ++i)
    tmp[i] += q[i];

  return tmp;
}



/**
 * Subtraction of two tensors of general @tparam rank.
 *
 * @relates Tensor
 */
template <int rank, int dim, typename Number, typename OtherNumber>
inline
Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type>
operator- (const Tensor<rank,dim,Number> &p, const Tensor<rank,dim,OtherNumber> &q)
{
  Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type> tmp (p);

  for (unsigned int i=0; i<dim; ++i)
    tmp[i] -= q[i];

  return tmp;
}



/**
 * Contract a tensor of rank 1 with a tensor of rank 1. The result is
 * <tt>sum_j src1[j] src2[j]</tt>.
 *
 * @relates Tensor
 */
template <int dim, typename Number, typename OtherNumber>
inline
typename ProductType<Number,OtherNumber>::type
contract (const Tensor<1,dim,Number> &src1,
          const Tensor<1,dim,OtherNumber> &src2)
{
  typename ProductType<Number,OtherNumber>::type res
    = typename ProductType<Number,OtherNumber>::type();
  for (unsigned int i=0; i<dim; ++i)
    res += src1[i] * src2[i];

  return res;
}


/**
 * Multiplication operator performing a contraction of the last index of the
 * first argument and the first index of the second argument. This function
 * therefore does the same as the corresponding <tt>contract</tt> function,
 * but returns the result as a return value, rather than writing it into the
 * reference given as the first argument to the <tt>contract</tt> function.
 *
 * Note that for the <tt>Tensor</tt> class, the multiplication operator only
 * performs a contraction over a single pair of indices. This is in contrast
 * to the multiplication operator for symmetric tensors, which does the double
 * contraction.
 *
 * @relates Tensor
 */
template <int dim, typename Number, typename OtherNumber>
inline
typename ProductType<Number,OtherNumber>::type
operator * (const Tensor<1,dim,Number> &src1,
            const Tensor<1,dim,OtherNumber> &src2)
{
  return contract(src1, src2);
}



DEAL_II_NAMESPACE_CLOSE

#endif
