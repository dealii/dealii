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


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/utilities.h>
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

// general template; specialized for rank == 0
template <int rank_, int dim, typename Number = double> class Tensor;
template <int dim, typename Number> class Tensor<0,dim,Number>;



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
   * type statically. In case of a a tensor of rank 0 this is just the scalar
   * number type Number.
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
   * Return a reference to the encapsulated Number object. Since rank-0
   * tensors are scalars, this is a natural operation.
   *
   * This is the non-const conversion operator that returns a writable
   * reference.
   */
  operator Number &();

  /**
   * Return a reference to the encapsulated Number object. Since rank-0
   * tensors are scalars, this is a natural operation.
   *
   * This is the const conversion operator that returns a read-only
   * reference.
   */
  operator const Number &() const;

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
   * Return the Frobenius-norm of a tensor, i.e. the square root of the sum
   * of the absolute squares of all entries. For the present case of rank-1
   * tensors, this equals the usual <tt>l<sub>2</sub></tt> norm of the
   * vector.
   */
  real_type norm () const;

  /**
   * Return the square of the Frobenius-norm of a tensor, i.e. the sum of
   * the absolute squares of all entries.
   */
  real_type norm_square () const;

  /**
   * Read or write the data of this object to or from a stream for the purpose
   * of serialization
   */
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version);

private:
  /**
   * Internal type declaration that is used to specialize the return type
   * of operator[]() for Tensor<1,dim,Number>
   */
  typedef Number tensor_type;

  /**
   * The value of this scalar object.
   */
  Number value;

  /**
   * Help function for unroll.
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
   * Declare a type that holds real-valued numbers with the same precision
   * as the template argument to this class. For std::complex<number>, this
   * corresponds to type number, and it is equal to Number for all other
   * cases. See also the respective field in Vector<Number>.
   *
   * This typedef is used to represent the return type of norms.
   */
  typedef typename numbers::NumberTraits<Number>::real_type real_type;

  /**
   * Type of objects encapsulated by this container and returned by
   * operator[](). This is a tensor of lower rank for a general tensor, and
   * a scalar number type for Tensor<1,dim,Number>.
   */
  typedef typename Tensor<rank_-1,dim,Number>::tensor_type value_type;

  /**
   * Declare an array type which can be used to initialize an object of this
   * type statically.
   */
  typedef typename Tensor<rank_-1,dim,Number>::array_type
  array_type[(dim != 0) ? dim : 1];

  /**
   * Constructor. Initialize all entries to zero if
   * <tt>initialize==true</tt>; this is the default behaviour.
   */
  explicit
  Tensor (const bool initialize = true);

  /**
   * Copy constructor.
   */
  Tensor (const Tensor<rank_,dim,Number> &initializer);

  /**
   * Constructor, where the data is copied from a C-style array.
   */
  Tensor (const array_type &initializer);

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
   */
  value_type &operator [] (const unsigned int i);

  /**
   * Read-only access operator.
   */
  const value_type &operator[](const unsigned int i) const;

  /**
   * Read access using TableIndices <tt>indices</tt>
   */
  Number operator [] (const TableIndices<rank_> &indices) const;

  /**
   * Read and write access using TableIndices <tt>indices</tt>
   */
  Number &operator [] (const TableIndices<rank_> &indices);

  /**
   * Copy assignment operator.
   */
  Tensor &operator = (const Tensor<rank_,dim,Number> &rhs);

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
   *
   * @relates EnableIfScalar
   */
  template <typename OtherNumber,
            typename = typename EnableIfScalar<OtherNumber>::type>
  Tensor<rank_,dim,Number> &operator = (const OtherNumber d);

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
  Tensor<rank_,dim,Number> &operator += (const Tensor<rank_,dim,OtherNumber> &);

  /**
   * Subtract another tensor.
   */
  template <typename OtherNumber>
  Tensor<rank_,dim,Number> &operator -= (const Tensor<rank_,dim,OtherNumber> &);

  /**
   * Scale the tensor by <tt>factor</tt>, i.e. multiply all components by
   * <tt>factor</tt>.
   */
  template <typename OtherNumber>
  Tensor<rank_,dim,Number> &operator *= (const OtherNumber factor);

  /**
   * Scale the vector by <tt>1/factor</tt>.
   */
  template <typename OtherNumber>
  Tensor<rank_,dim,Number> &operator /= (const OtherNumber factor);

  /**
   * Unary minus operator. Negate all entries of a tensor.
   */
  Tensor<rank_,dim,Number>   operator - () const;

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
   * Return the Frobenius-norm of a tensor, i.e. the square root of the sum
   * of the absolute squares of all entries. For the present case of rank-1
   * tensors, this equals the usual <tt>l<sub>2</sub></tt> norm of the
   * vector.
   */
  real_type norm () const;

  /**
   * Return the square of the Frobenius-norm of a tensor, i.e. the sum of
   * the absolute squares of all entries.
   */
  real_type norm_square () const;

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
   * Returns an unrolled index in the range [0,dim^rank-1] for the element of
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
   * Exception.
   */
  DeclException1 (ExcInvalidTensorContractionIndex,
                  int,
                  << "You have requested contraction of tensors over index "
                  << arg1
                  << ", but this is not possible for tensors of the current type.");

private:
  /**
   * Internal type declaration that is used to specialize the return type
   * of operator[]() for Tensor<1,dim,Number>
   */
  typedef Tensor<rank_, dim, Number> tensor_type;

  /**
   * Array of tensors holding the subelements.
   */
  Tensor<rank_-1, dim, Number> values[(dim != 0) ? dim : 1];

  /**
   * Help function for unroll.
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


#ifndef DOXYGEN
/*---------------------- Inline functions: Tensor<0,dim> ---------------------*/


template <int dim,typename Number>
inline
Tensor<0,dim,Number>::Tensor ()
{
  value = value_type();
}


template <int dim, typename Number>
inline
Tensor<0,dim,Number>::Tensor (const Tensor<0,dim,Number> &p)
{
  Assert(dim != 0 || p.value == Number(),
         ExcMessage("Creation of a Tensor<0,0,Number> object with a non-zero scalar requested."));
  value = p.value;
}


template <int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<0,dim,Number>::Tensor (const OtherNumber initializer)
{
  Assert(dim != 0 || initializer == OtherNumber(),
         ExcMessage("Creation of a Tensor<0,0,Number> object with a non-zero scalar requested."));
  value = initializer;
}


template <int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<0,dim,Number>::Tensor (const Tensor<0,dim,OtherNumber> &p)
{
  Assert(dim != 0 || p.value == OtherNumber(),
         ExcMessage("Cannot return a non-zero scalar from a Tensor<0,0,Number> object."));
  value = p.value;
}


template <int dim, typename Number>
inline
Tensor<0,dim,Number>::operator Number &()
{
  Assert(dim != 0 || value == Number(),
         ExcMessage("Cannot return a non-zero scalar from a Tensor<0,0,Number> object."));
  return value;
}


template <int dim, typename Number>
inline
Tensor<0,dim,Number>::operator const Number &() const
{
  Assert(dim != 0 || value == Number(),
         ExcMessage("Cannot assign a non-zero scalar to a Tensor<0,0,Number> object."));
  return value;
}


template <int dim, typename Number>
inline
Tensor<0,dim,Number> &Tensor<0,dim,Number>::operator = (const Tensor<0,dim,Number> &p)
{
  Assert(dim != 0 || p.value == Number(),
         ExcMessage("Cannot assign a non-zero scalar to a Tensor<0,0,Number> object."));
  value = p.value;
  return *this;
}


template <int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<0,dim,Number> &Tensor<0,dim,Number>::operator = (const Tensor<0,dim,OtherNumber> &p)
{
  Assert(dim != 0 || p.value == OtherNumber(),
         ExcMessage("Cannot assign a non-zero scalar to a Tensor<0,0,Number> object."));
  value = p.value;
  return *this;
}


template <int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<0,dim,Number> &Tensor<0,dim,Number>::operator = (const OtherNumber d)
{
  Assert(dim != 0 || d == OtherNumber(),
         ExcMessage("Cannot assign a non-zero scalar to a Tensor<0,0,Number> object."));
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
template <typename OtherNumber>
inline
void
Tensor<0, dim, Number>::unroll_recursion (Vector<OtherNumber> &result,
                                          unsigned int        &index) const
{
  result[index] = value;
  ++index;
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


/*-------------------- Inline functions: Tensor<rank,dim> --------------------*/


namespace internal
{
  // TODO: Think about refactoring this into the TableIndices class as a
  // general, polymorphic for extracting an item out of an object with
  // nested identifiers.
  template<int rank_> struct TensorIndicesHelper
  {
    // used for implementing Tensor<rank,dim>::operator[] with TableIndices
    // tail recursive call to form up access to
    //   tensor[indices[0]][indices[1]]...[indices[rank_]]
    template<int rank, int dim, typename Number>
    static inline
    Number &extract(Tensor<rank_,dim,Number> &t, const TableIndices<rank> &indices)
    {
      Assert (indices[rank - rank_]<dim, ExcIndexRange (indices[rank - rank_], 0, dim));
      return TensorIndicesHelper<rank_ - 1>::template extract<rank, dim, Number>(
        t[indices[rank - rank_]], indices);
    }
  };

  template<> struct TensorIndicesHelper<1>
  {
    template<int rank, int dim, typename Number>
    static inline
    Number &extract(Tensor<1,dim,Number> &t, const TableIndices<rank> &indices)
    {
      Assert (indices[rank - 1]<dim, ExcIndexRange (indices[rank - 1], 0, dim));
      return t[indices[rank-1]];
    }
  };
} /* internal */


template <int rank_, int dim, typename Number>
inline
Tensor<rank_,dim,Number>::Tensor (const bool initialize)
{
  if (initialize)
    // need to create an object Number() to initialize to zero to avoid
    // confusion with Tensor::operator=(scalar) when using something like
    // Tensor<1,dim,Tensor<1,dim,Number> >.
    for (unsigned int i=0; i!=dim; ++i)
      values[i] = value_type();
}


template <int rank_, int dim, typename Number>
inline
Tensor<rank_,dim,Number>::Tensor (const Tensor<rank_,dim,Number> &initializer)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] = initializer[i];
}


template <int rank_, int dim, typename Number>
inline
Tensor<rank_,dim,Number>::Tensor (const array_type &initializer)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] = initializer[i];
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<rank_,dim,Number>::Tensor (const Tensor<rank_,dim,OtherNumber> &initializer)
{
  for (unsigned int i=0; i!=dim; ++i)
    values[i] = initializer[i];
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


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<rank_,dim,Number>::Tensor
(const Tensor<1,dim,Tensor<rank_-1,dim,OtherNumber> > &initializer)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] = initializer[i];
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<rank_,dim,Number>::operator
Tensor<1,dim,Tensor<rank_-1,dim,OtherNumber> > () const
{
  return Tensor<1,dim,Tensor<rank_-1,dim,Number> > (values);
}


template <int rank_, int dim, typename Number>
inline
typename Tensor<rank_,dim,Number>::value_type &
Tensor<rank_,dim,Number>::operator[] (const unsigned int i)
{
  Assert (i<dim, ExcIndexRange(i, 0, dim));
  return values[i];
}


template <int rank_, int dim, typename Number>
inline
const typename Tensor<rank_,dim,Number>::value_type &
Tensor<rank_,dim,Number>::operator[] (const unsigned int i) const
{
  Assert (i<dim, ExcIndexRange(i, 0, dim));
  return values[i];
}


template <int rank_, int dim, typename Number>
inline
Number
Tensor<rank_,dim,Number>::operator[] (const TableIndices<rank_> &indices) const
{
  Assert (indices[0]<dim, ExcIndexRange (indices[0], 0, dim));
  return internal::TensorIndicesHelper<rank_>::extract(*this, indices);
}


template <int rank_, int dim, typename Number>
inline
Number &
Tensor<rank_,dim,Number>::operator[] (const TableIndices<rank_> &indices)
{
  Assert (indices[0]<dim, ExcIndexRange (indices[0], 0, dim));
  return internal::TensorIndicesHelper<rank_>::extract(*this, indices);
}


template <int rank_, int dim, typename Number>
inline
Tensor<rank_,dim,Number> &
Tensor<rank_,dim,Number>::operator = (const Tensor<rank_,dim,Number> &t)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] = t.values[i];
  return *this;
}


// At some places in the library, we have Point<0> for formal reasons
// (e.g., we sometimes have Quadrature<dim-1> for faces, so we have
// Quadrature<0> for dim=1, and then we have Point<0>). To avoid warnings
// in the above function that the loop end check always fails, we
// implement this function here
template <>
inline
Tensor<1,0,double> &Tensor<1,0,double>::operator = (const Tensor<1,0,double> &)
{
  return *this;
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<rank_,dim,Number> &
Tensor<rank_,dim,Number>::operator = (const Tensor<rank_,dim,OtherNumber> &t)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] = t.values[i];
  return *this;
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber, typename>
inline
Tensor<rank_,dim,Number> &
Tensor<rank_,dim,Number>::operator = (const OtherNumber d)
{
  Assert (d == OtherNumber(), ExcMessage ("Only assignment with zero is allowed"));
  (void) d;

  for (unsigned int i=0; i<dim; ++i)
    values[i] = Number();
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
Tensor<rank_,dim,Number> &
Tensor<rank_,dim,Number>::operator *= (const OtherNumber s)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] *= s;
  return *this;
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<rank_,dim,Number> &
Tensor<rank_,dim,Number>::operator /= (const OtherNumber s)
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
typename Tensor<rank_,dim,Number>::real_type
Tensor<rank_,dim,Number>::norm () const
{
  return std::sqrt (norm_square());
}


template <int rank_, int dim, typename Number>
inline
typename Tensor<rank_,dim,Number>::real_type
Tensor<rank_,dim,Number>::norm_square () const
{
  real_type s = 0;
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
    values[i] = value_type();
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


#endif /* DOXYGEN */
/* ----------------- Non-member functions operating on tensors. ------------- */


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
 * @name Vector space operations on Tensor objects:
 */
//@{

/**
 * Multiplication of a tensor of general rank with a scalar number from the
 * right.
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

//@}


DEAL_II_NAMESPACE_CLOSE

#endif
