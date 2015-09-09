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

#ifndef dealii__tensor_h
#define dealii__tensor_h

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/utilities.h>

#include <cmath>
#include <ostream>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations:

template <int dim, typename Number> class Point;
template <int rank_, int dim, typename Number = double> class Tensor;



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

#ifdef DEAL_II_WITH_CXX11
  /**
   * Copy constructor.
   */
  Tensor (const Tensor<0,dim,Number> &initializer) = default;
  // always implicitly created, in case of C++11 we make this explicit (for
  // the sake of documentation).
#endif

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

#ifdef DEAL_II_WITH_CXX11
  /**
   * Copy assignment operator.
   */
  Tensor<0,dim,Number> &operator = (const Tensor<0,dim,Number> &rhs) = default;
  // always implicitly created, in case of C++11 we make this explicit (for
  // the sake of documentation).
#endif

  /**
   * Assignment from tensors with different underlying scalar type.
   * This obviously requires that the @p OtherNumber type is convertible to @p
   * Number.
   */
  template <typename OtherNumber>
  Tensor<0,dim,Number> &operator = (const Tensor<0,dim,OtherNumber> &rhs);

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
  // ... avoid a compiler warning in case of dim == 0 and ensure that the
  // array always has positive size.

  /**
   * Constructor. Initialize all entries to zero.
   */
  Tensor ();

#ifdef DEAL_II_WITH_CXX11
  /**
   * Copy constructor.
   */
  Tensor (const Tensor<rank_,dim,Number> &initializer) = default;
  // always implicitly created, in case of C++11 we make this explicit (for
  // the sake of documentation).
#endif

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

#ifdef DEAL_II_WITH_CXX11
  /**
   * Copy assignment operator.
   */
  Tensor &operator = (const Tensor<rank_,dim,Number> &rhs) = default;
  // always implicitly created, in case of C++11 we make this explicit (for
  // the sake of documentation).
#endif

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
  Tensor<rank_,dim,Number> &operator = (const Number d);

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
  // ... avoid a compiler warning in case of dim == 0 and ensure that the
  // array always has positive size.

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
  : value()
{
}


template <int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<0,dim,Number>::Tensor (const OtherNumber initializer)
{
  value = initializer;
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
Tensor<0,dim,Number>::operator Number &()
{
  Assert(dim != 0, ExcMessage("Cannot access an object of type Tensor<0,0,Number>"));
  return value;
}


template <int dim, typename Number>
inline
Tensor<0,dim,Number>::operator const Number &() const
{
  Assert(dim != 0, ExcMessage("Cannot access an object of type Tensor<0,0,Number>"));
  return value;
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
  Assert(dim != 0, ExcMessage("Cannot access an object of type Tensor<0,0,Number>"));
  return numbers::NumberTraits<Number>::abs (value);
}


template <int dim, typename Number>
inline
typename Tensor<0,dim,Number>::real_type
Tensor<0,dim,Number>::norm_square () const
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

    template<int rank, int dim, typename Number>
    static inline
    const Number &extract(const Tensor<rank_,dim,Number> &t, const TableIndices<rank> &indices)
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

    template<int rank, int dim, typename Number>
    static inline
    const Number &extract(const Tensor<1,dim,Number> &t, const TableIndices<rank> &indices)
    {
      Assert (indices[rank - 1]<dim, ExcIndexRange (indices[rank - 1], 0, dim));
      return t[indices[rank-1]];
    }
  };
} /* internal */


template <int rank_, int dim, typename Number>
inline
Tensor<rank_,dim,Number>::Tensor ()
{
  // All members of the c-style array values are already default initialized
  // and thus all values are already set to zero recursively.
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
Tensor<rank_,dim,Number>::
operator Tensor<1,dim,Tensor<rank_-1,dim,OtherNumber> > () const
{
  return Tensor<1,dim,Tensor<rank_-1,dim,Number> > (values);
}


template <int rank_, int dim, typename Number>
inline
typename Tensor<rank_,dim,Number>::value_type &
Tensor<rank_,dim,Number>::operator[] (const unsigned int i)
{
  Assert(dim != 0, ExcMessage("Cannot access an object of type Tensor<rank_,0,Number>"));
  Assert (i<dim, ExcIndexRange(i, 0, dim));
  return values[i];
}


template <int rank_, int dim, typename Number>
inline
const typename Tensor<rank_,dim,Number>::value_type &
Tensor<rank_,dim,Number>::operator[] (const unsigned int i) const
{
  Assert(dim != 0, ExcMessage("Cannot access an object of type Tensor<rank_,0,Number>"));
  Assert (i<dim, ExcIndexRange(i, 0, dim));
  return values[i];
}


template <int rank_, int dim, typename Number>
inline
Number
Tensor<rank_,dim,Number>::operator[] (const TableIndices<rank_> &indices) const
{
  Assert(dim != 0, ExcMessage("Cannot access an object of type Tensor<rank_,0,Number>"));
  Assert (indices[0]<dim, ExcIndexRange (indices[0], 0, dim));
  return internal::TensorIndicesHelper<rank_>::template extract<rank_, dim, Number>(*this, indices);
}


template <int rank_, int dim, typename Number>
inline
Number &
Tensor<rank_,dim,Number>::operator[] (const TableIndices<rank_> &indices)
{
  Assert(dim != 0, ExcMessage("Cannot access an object of type Tensor<rank_,0,Number>"));
  Assert (indices[0]<dim, ExcIndexRange (indices[0], 0, dim));
  return internal::TensorIndicesHelper<rank_>::template extract<rank_, dim, Number>(*this, indices);
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline
Tensor<rank_,dim,Number> &
Tensor<rank_,dim,Number>::operator = (const Tensor<rank_,dim,OtherNumber> &t)
{
  // avoid a warning with icc that complains about an empty loop. Even in
  // case of dim==0 we have an element available (as a workaround for
  // another set of warnings about an empty c-style array *sigh).
  // If dim==0, ust copy this bogus element:
  for (unsigned int i = 0; i < (dim != 0 ? dim : 1); ++i)
    values[i] = t.values[i];
  return *this;
}


template <int rank_, int dim, typename Number>
inline
Tensor<rank_,dim,Number> &
Tensor<rank_,dim,Number>::operator = (const Number d)
{
  Assert (d == Number(), ExcMessage ("Only assignment with zero is allowed"));
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
  real_type s = real_type();
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


/* ----------------- Non-member functions operating on tensors. ------------ */


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
 * @name Output functions for Tensor objects
 */
//@{

/**
 * Output operator for tensors. Print the elements consecutively, with a space
 * in between, two spaces between rank 1 subtensors, three between rank 2 and
 * so on.
 *
 * @relates Tensor
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
 * @relates Tensor<0,dim,Number>
 */
template <int dim, typename Number>
inline
std::ostream &operator << (std::ostream &out, const Tensor<0,dim,Number> &p)
{
  out << static_cast<Number>(p);
  return out;
}


//@}
/**
 * @name Vector space operations on Tensor objects:
 */
//@{

/**
 * Scalar multiplication of a tensor of rank 0 with a scalar from the left.
 *
 * @relates Tensor<0,dim,Number>
 * @relates EnableIfScalar
 */
template <int dim, typename Number, typename OtherNumber>
inline
Tensor<0,dim,typename ProductType<typename EnableIfScalar<OtherNumber>::type, Number>::type>
operator * (const OtherNumber           factor,
            const Tensor<0,dim,Number> &t)
{
  return factor * static_cast<Number>(t);
}


/**
 * Scalar multiplication of a tensor of rank 0 with a scalar from the right.
 *
 * @relates Tensor<0,dim,Number>
 * @relates EnableIfScalar
 */
template <int dim, typename Number, typename OtherNumber>
inline
Tensor<0,dim,typename ProductType<Number, typename EnableIfScalar<OtherNumber>::type>::type>
operator * (const Tensor<0,dim,Number> &t,
            const OtherNumber           factor)
{
  return static_cast<Number>(t) * factor;
}


/**
 * Division of a tensor of rank 0 by a scalar number.
 *
 * @relates Tensor<0,dim,Number>
 * @relates EnableIfScalar
 */
template <int dim, typename Number, typename OtherNumber>
inline
Tensor<0,dim,typename ProductType<Number, typename EnableIfScalar<OtherNumber>::type>::type>
operator / (const Tensor<0,dim,Number> &t,
            const OtherNumber           factor)
{
  return static_cast<Number>(t) / factor;
}


/**
 * Add two tensors of rank 0.
 *
 * @relates Tensor<0,dim,Number>
 */
template <int dim, typename Number, typename OtherNumber>
inline
Tensor<0, dim, typename ProductType<Number, OtherNumber>::type>
operator+ (const Tensor<0,dim,Number> &p, const Tensor<0,dim,OtherNumber> &q)
{
  return static_cast<const Number &>(p) + static_cast<const OtherNumber &>(q);
}


/**
 * Subtract two tensors of rank 0.
 *
 * @relates Tensor<0,dim,Number>
 */
template <int dim, typename Number, typename OtherNumber>
inline
Tensor<0, dim, typename ProductType<Number, OtherNumber>::type>
operator- (const Tensor<0,dim,Number> &p, const Tensor<0,dim,OtherNumber> &q)
{
  return static_cast<const Number &>(p) - static_cast<const OtherNumber &>(q);
}


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
          typename OtherNumber>
inline
Tensor<rank,dim,typename ProductType<Number, typename EnableIfScalar<OtherNumber>::type>::type>
operator * (const Tensor<rank,dim,Number> &t,
            const OtherNumber              factor)
{
  // recurse over the base objects
  Tensor<rank,dim,typename ProductType<Number,OtherNumber>::type> tt;
  for (unsigned int d=0; d<dim; ++d)
    tt[d] = t[d] * factor;
  return tt;
}


#ifdef DEAL_II_GCC_COMPLEX_CONV_BUG
template <int dim, typename Number, typename OtherNumber>
inline
Tensor<1,dim,typename ProductType<Number, typename EnableIfScalar<OtherNumber>::type>::type>
operator * (const Tensor<1,dim,Number> &t,
            const OtherNumber           factor)
{
  typedef typename ProductType<Number,OtherNumber>::type product_type;
  Tensor<1,dim,product_type> tt;
  for (unsigned int d=0; d<dim; ++d)
    tt[d] = product_type(t[d]) * product_type(factor);
  return tt;
}
#endif


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
          typename OtherNumber>
inline
Tensor<rank,dim,typename ProductType<typename EnableIfScalar<Number>::type, OtherNumber>::type>
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
          typename OtherNumber>
inline
Tensor<rank,dim,typename ProductType<Number, typename EnableIfScalar<OtherNumber>::type>::type>
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
/**
 * @name Contraction operations on Tensors
 */
//@{

/**
 * Returns the contraction of two Tensors of rank 0.
 *
 * @relates Tensor<0,dim,Number>
 */
template <int dim, typename Number, typename OtherNumber>
inline
typename ProductType<Number, OtherNumber>::type
operator* (const Tensor<0,dim,Number> &p, const Tensor<0,dim,OtherNumber> &q)
{
  return static_cast<const Number &>(p) * static_cast<const OtherNumber &>(q);
}

//@}


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


/**
 * Double contract two tensors of rank 2, thus computing the Frobenius inner
 * product <tt> sum<sub>i,j</sub> src1[i][j]*src2[i][j]</tt>.
 *
 * @relates Tensor
 * @author Guido Kanschat, 2000
 */
template <int dim, typename Number>
inline
Number double_contract (const Tensor<2, dim, Number> &src1,
                        const Tensor<2, dim, Number> &src2)
{
  Number res = 0.;
  for (unsigned int i=0; i<dim; ++i)
    res += contract(src1[i],src2[i]);

  return res;
}


/**
 * Contract a tensor of rank 2 with a tensor of rank 1. The result is
 * <tt>dest[i] = sum_j src1[i][j] src2[j]</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 1998
 */
template <int dim, typename Number>
inline
void contract (Tensor<1,dim,Number>       &dest,
               const Tensor<2,dim,Number> &src1,
               const Tensor<1,dim,Number> &src2)
{
  for (unsigned int i=0; i<dim; ++i)
    {
      dest[i] = src1[i][0] * src2[0];
      for (unsigned int j=1; j<dim; ++j)
        dest[i] += src1[i][j] * src2[j];
    }
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
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, typename Number>
Tensor<1,dim,Number>
operator * (const Tensor<2,dim,Number> &src1,
            const Tensor<1,dim,Number> &src2)
{
  Tensor<1,dim,Number> dest;
  for (unsigned int i=0; i<dim; ++i)
    {
      dest[i] = src1[i][0] * src2[0];
      for (unsigned int j=1; j<dim; ++j)
        dest[i] += src1[i][j] * src2[j];
    }
  return dest;
}


/**
 * Contract a tensor of rank 1 with a tensor of rank 2. The result is
 * <tt>dest[i] = sum_j src1[j] src2[j][i]</tt>.
 *
 * @relates Tensor
 * @author Guido Kanschat, 2001
 */
template <int dim, typename Number>
inline
void contract (Tensor<1,dim,Number>       &dest,
               const Tensor<1,dim,Number> &src1,
               const Tensor<2,dim,Number> &src2)
{
  for (unsigned int i=0; i<dim; ++i)
    {
      dest[i] = src1[0] * src2[0][i];
      for (unsigned int j=1; j<dim; ++j)
        dest[i] += src1[j] * src2[j][i];
    }
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
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, typename Number>
inline
Tensor<1,dim,Number>
operator * (const Tensor<1,dim,Number> &src1,
            const Tensor<2,dim,Number> &src2)
{
  Tensor<1,dim,Number> dest;
  for (unsigned int i=0; i<dim; ++i)
    {
      dest[i] = src1[0] * src2[0][i];
      for (unsigned int j=1; j<dim; ++j)
        dest[i] += src1[j] * src2[j][i];
    }
  return dest;
}


/**
 * Contract a tensor of rank 2 with a tensor of rank 2. The result is
 * <tt>dest[i][k] = sum_j src1[i][j] src2[j][k]</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 1998
 */
template <int dim, typename Number>
inline
void contract (Tensor<2,dim,Number>       &dest,
               const Tensor<2,dim,Number> &src1,
               const Tensor<2,dim,Number> &src2)
{
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      {
        dest[i][j] = src1[i][0] * src2[0][j];
        for (unsigned int k=1; k<dim; ++k)
          dest[i][j] += src1[i][k] * src2[k][j];
      }
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
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, typename Number>
inline
Tensor<2,dim,Number>
operator * (const Tensor<2,dim,Number> &src1,
            const Tensor<2,dim,Number> &src2)
{
  Tensor<2,dim,Number> dest;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        dest[i][j] += src1[i][k] * src2[k][j];
  return dest;
}


/**
 * Contract a tensor of rank 2 with a tensor of rank 2. The contraction is
 * performed over index <tt>index1</tt> of the first tensor, and
 * <tt>index2</tt> of the second tensor. Thus, if <tt>index1==2</tt>,
 * <tt>index2==1</tt>, the result is the usual contraction, but if for example
 * <tt>index1==1</tt>, <tt>index2==2</tt>, then the result is <tt>dest[i][k] =
 * sum_j src1[j][i] src2[k][j]</tt>.
 *
 * Note that the number of the index is counted from 1 on, not from zero as
 * usual.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 1998
 */
template <int dim, typename Number>
inline
void contract (Tensor<2,dim,Number>       &dest,
               const Tensor<2,dim,Number> &src1,   const unsigned int index1,
               const Tensor<2,dim,Number> &src2,   const unsigned int index2)
{
  dest.clear ();

  switch (index1)
    {
    case 1:
      switch (index2)
        {
        case 1:
          for (unsigned int i=0; i<dim; ++i)
            for (unsigned int j=0; j<dim; ++j)
              for (unsigned int k=0; k<dim; ++k)
                dest[i][j] += src1[k][i] * src2[k][j];
          break;
        case 2:
          for (unsigned int i=0; i<dim; ++i)
            for (unsigned int j=0; j<dim; ++j)
              for (unsigned int k=0; k<dim; ++k)
                dest[i][j] += src1[k][i] * src2[j][k];
          break;

        default:
          Assert (false,
                  (typename Tensor<2,dim,Number>::ExcInvalidTensorContractionIndex (index2)));
        };
      break;
    case 2:
      switch (index2)
        {
        case 1:
          for (unsigned int i=0; i<dim; ++i)
            for (unsigned int j=0; j<dim; ++j)
              for (unsigned int k=0; k<dim; ++k)
                dest[i][j] += src1[i][k] * src2[k][j];
          break;
        case 2:
          for (unsigned int i=0; i<dim; ++i)
            for (unsigned int j=0; j<dim; ++j)
              for (unsigned int k=0; k<dim; ++k)
                dest[i][j] += src1[i][k] * src2[j][k];
          break;

        default:
          Assert (false,
                  (typename Tensor<2,dim,Number>::ExcInvalidTensorContractionIndex (index2)));
        };
      break;

    default:
      Assert (false, (typename Tensor<2,dim,Number>::ExcInvalidTensorContractionIndex (index1)));
    };
}


/**
 * Contract a tensor of rank 3 with a tensor of rank 1. The contraction is
 * performed over index <tt>index1</tt> of the first tensor.
 *
 * Note that the number of the index is counted from 1 on, not from zero as
 * usual.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 1998
 */
template <int dim, typename Number>
inline
void contract (Tensor<2,dim,Number>       &dest,
               const Tensor<3,dim,Number> &src1,   const unsigned int index1,
               const Tensor<1,dim,Number> &src2)
{
  dest.clear ();

  switch (index1)
    {
    case 1:
      for (unsigned int i=0; i<dim; ++i)
        for (unsigned int j=0; j<dim; ++j)
          for (unsigned int k=0; k<dim; ++k)
            dest[i][j] += src1[k][i][j] * src2[k];
      break;

    case 2:
      for (unsigned int i=0; i<dim; ++i)
        for (unsigned int j=0; j<dim; ++j)
          for (unsigned int k=0; k<dim; ++k)
            dest[i][j] += src1[i][k][j] * src2[k];
      break;

    case 3:
      for (unsigned int i=0; i<dim; ++i)
        for (unsigned int j=0; j<dim; ++j)
          for (unsigned int k=0; k<dim; ++k)
            dest[i][j] += src1[i][j][k] * src2[k];
      break;

    default:
      Assert (false,
              (typename Tensor<2,dim,Number>::ExcInvalidTensorContractionIndex (index1)));
    };
}


/**
 * Contract a tensor of rank 3 with a tensor of rank 2. The result is
 * <tt>dest[i][j][l] = sum_k src1[i][j][k] src2[k][l]</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 1998
 */
template <int dim, typename Number>
inline
void contract (Tensor<3,dim,Number>       &dest,
               const Tensor<3,dim,Number> &src1,
               const Tensor<2,dim,Number> &src2)
{
  dest.clear ();
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        for (unsigned int l=0; l<dim; ++l)
          dest[i][j][k] += src1[i][j][l] * src2[l][k];
}


/**
 * Contract a tensor of rank 3 with a tensor of rank 2. The contraction is
 * performed over index <tt>index1</tt> of the first tensor, and
 * <tt>index2</tt> of the second tensor. Thus, if <tt>index1==3</tt>,
 * <tt>index2==1</tt>, the result is the usual contraction, but if for example
 * <tt>index1==1</tt>, <tt>index2==2</tt>, then the result is
 * <tt>dest[i][j][k] = sum_l src1[l][i][j] src2[k][l]</tt>.
 *
 * Note that the number of the index is counted from 1 on, not from zero as
 * usual.
 *
 * @relates Tensor
 */
template <int dim, typename Number>
inline
void contract (Tensor<3,dim,Number>       &dest,
               const Tensor<3,dim,Number> &src1, const unsigned int index1,
               const Tensor<2,dim,Number> &src2, const unsigned int index2)
{
  dest.clear ();

  switch (index1)
    {
    case 1:
      switch (index2)
        {
        case 1:
          for (unsigned int i=0; i<dim; ++i)
            for (unsigned int j=0; j<dim; ++j)
              for (unsigned int k=0; k<dim; ++k)
                for (unsigned int l=0; l<dim; ++l)
                  dest[i][j][k] += src1[l][i][j] * src2[l][k];
          break;
        case 2:
          for (unsigned int i=0; i<dim; ++i)
            for (unsigned int j=0; j<dim; ++j)
              for (unsigned int k=0; k<dim; ++k)
                for (unsigned int l=0; l<dim; ++l)
                  dest[i][j][k] += src1[l][i][j] * src2[k][l];
          break;
        default:
          Assert (false,
                  (typename Tensor<2,dim,Number>::ExcInvalidTensorContractionIndex (index2)));
        }

      break;
    case 2:
      switch (index2)
        {
        case 1:
          for (unsigned int i=0; i<dim; ++i)
            for (unsigned int j=0; j<dim; ++j)
              for (unsigned int k=0; k<dim; ++k)
                for (unsigned int l=0; l<dim; ++l)
                  dest[i][j][k] += src1[i][l][j] * src2[l][k];
          break;
        case 2:
          for (unsigned int i=0; i<dim; ++i)
            for (unsigned int j=0; j<dim; ++j)
              for (unsigned int k=0; k<dim; ++k)
                for (unsigned int l=0; l<dim; ++l)
                  dest[i][j][k] += src1[i][l][j] * src2[k][l];
          break;
        default:
          Assert (false,
                  (typename Tensor<2,dim,Number>::ExcInvalidTensorContractionIndex (index2)));
        }

      break;
    case 3:
      switch (index2)
        {
        case 1:
          for (unsigned int i=0; i<dim; ++i)
            for (unsigned int j=0; j<dim; ++j)
              for (unsigned int k=0; k<dim; ++k)
                for (unsigned int l=0; l<dim; ++l)
                  dest[i][j][k] += src1[i][j][l] * src2[l][k];
          break;
        case 2:
          for (unsigned int i=0; i<dim; ++i)
            for (unsigned int j=0; j<dim; ++j)
              for (unsigned int k=0; k<dim; ++k)
                for (unsigned int l=0; l<dim; ++l)
                  dest[i][j][k] += src1[i][j][l] * src2[k][l];
          break;
        default:
          Assert (false,
                  (typename Tensor<2,dim,Number>::ExcInvalidTensorContractionIndex (index2)));
        }

      break;
    default:
      Assert (false,
              (typename Tensor<3,dim,Number>::ExcInvalidTensorContractionIndex (index1)));
    }
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
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, typename Number>
inline
Tensor<3,dim,Number>
operator * (const Tensor<3,dim,Number> &src1,
            const Tensor<2,dim,Number> &src2)
{
  Tensor<3,dim,Number> dest;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        for (unsigned int l=0; l<dim; ++l)
          dest[i][j][k] += src1[i][j][l] * src2[l][k];
  return dest;
}


/**
 * Contract a tensor of rank 2 with a tensor of rank 3. The result is
 * <tt>dest[i][j][l] = sum_k src1[i][k] src2[k][j][l]</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 1998
 */
template <int dim, typename Number>
inline
void contract (Tensor<3,dim,Number>       &dest,
               const Tensor<2,dim,Number> &src1,
               const Tensor<3,dim,Number> &src2)
{
  dest.clear ();
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        for (unsigned int l=0; l<dim; ++l)
          dest[i][j][k] += src1[i][l] * src2[l][j][k];
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
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, typename Number>
inline
Tensor<3,dim,Number>
operator * (const Tensor<2,dim,Number> &src1,
            const Tensor<3,dim,Number> &src2)
{
  Tensor<3,dim,Number> dest;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        for (unsigned int l=0; l<dim; ++l)
          dest[i][j][k] += src1[i][l] * src2[l][j][k];
  return dest;
}


/**
 * Contract a tensor of rank 3 with a tensor of rank 3. The result is
 * <tt>dest[i][j][k][l] = sum_m src1[i][j][m] src2[m][k][l]</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 1998
 */
template <int dim, typename Number>
inline
Tensor<4,dim,Number>
operator * (const Tensor<3,dim,Number> &src1,
            const Tensor<3,dim,Number> &src2)
{
  Tensor<4,dim,Number> dest;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        for (unsigned int l=0; l<dim; ++l)
          for (unsigned int m=0; m<dim; ++m)
            dest[i][j][k][l] += src1[i][j][m] * src2[m][k][l];
  return dest;
}


/**
 * Contract the last two indices of <tt>src1</tt> with the two indices
 * <tt>src2</tt>, creating a rank-2 tensor. This is the matrix-vector product
 * analog operation between tensors of rank 4 and rank 2.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, typename Number>
inline
void double_contract (Tensor<2,dim,Number>       &dest,
                      const Tensor<4,dim,Number> &src1,
                      const Tensor<2,dim,Number> &src2)
{
  dest.clear ();
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        for (unsigned int l=0; l<dim; ++l)
          dest[i][j] += src1[i][j][k][l] * src2[k][l];
}


/**
 * Contract three tensors, corresponding to the matrix vector product
 * <i>u<sup>T</sup> A v</i>.
 *
 * @relates Tensor
 * @author Guido Kanschat, 2004
 */
template <int dim, typename Number>
inline
Number contract3 (const Tensor<1,dim,Number> &u,
                  const Tensor<2,dim,Number> &A,
                  const Tensor<1,dim,Number> &v)
{
  Number result = 0.;

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      result += u[i] * A[i][j] * v[j];
  return result;
}


/**
 * Compute the contraction of three tensors $s=\sum_{i,j,k}
 * a_{i}b_{ijk}c_{jk}$.
 *
 * @relates Tensor
 * @author Toby D. Young, 2011
 */
template <int dim, typename Number>
inline
Number
contract3 (const Tensor<1,dim,Number> &t1,
           const Tensor<3,dim,Number> &t2,
           const Tensor<2,dim,Number> &t3)
{
  Number s = 0;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        s += t1[i] * t2[i][j][k] * t3[j][k];
  return s;
}


/**
 * Compute the contraction of three tensors $s=\sum_{i,j,k}
 * a_{ij}b_{ijk}c_{k}$.
 *
 * @relates Tensor
 * @author Toby D. Young, 2011
 */
template <int dim, typename Number>
inline
Number
contract3 (const Tensor<2,dim,Number> &t1,
           const Tensor<3,dim,Number> &t2,
           const Tensor<1,dim,Number> &t3)
{
  Number s = 0;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        s += t1[i][j] * t2[i][j][k] * t3[k];
  return s;
}


/**
 * Compute the contraction of three tensors $s=\sum_{i,j,k,l}
 * a_{ij}b_{ijkl}c_{kl}$.
 *
 * @relates Tensor
 * @author Toby D. Young, 2011
 */
template <int dim, typename Number>
inline
Number
contract3 (const Tensor<2,dim,Number> &t1,
           const Tensor<4,dim,Number> &t2,
           const Tensor<2,dim,Number> &t3)
{
  Number s = 0;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        for (unsigned int l=0; l<dim; ++l)
          s += t1[i][j] * t2[i][j][k][l] * t3[k][l];
  return s;
}


/**
 * Form the outer product of two tensors of rank 1 and 1, i.e. <tt>dst[i][j] =
 * src1[i] * src2[j]</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2000
 */
template <int dim, typename Number>
void outer_product (Tensor<2,dim,Number>       &dst,
                    const Tensor<1,dim,Number> &src1,
                    const Tensor<1,dim,Number> &src2)
{
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      dst[i][j] = src1[i] * src2[j];
}


/**
 * Form the outer product of two tensors of rank 1 and 2, i.e.
 * <tt>dst[i][j][k] = src1[i] * src2[j][k]</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2000
 */
template <int dim, typename Number>
void outer_product (Tensor<3,dim,Number>       &dst,
                    const Tensor<1,dim,Number> &src1,
                    const Tensor<2,dim,Number> &src2)
{
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        dst[i][j][k] = src1[i] * src2[j][k];
}


/**
 * Form the outer product of two tensors of rank 2 and 1, i.e.
 * <tt>dst[i][j][k] = src1[i][j] * src2[k]</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2000
 */
template <int dim, typename Number>
void outer_product (Tensor<3,dim,Number>       &dst,
                    const Tensor<2,dim,Number> &src1,
                    const Tensor<1,dim,Number> &src2)
{
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        dst[i][j][k] = src1[i][j] * src2[k];
}


/**
 * Form the outer product of two tensors of rank 0 and 1, i.e. <tt>dst[i] =
 * src1 * src2[i]</tt>. Of course, this is only a scaling of <tt>src2</tt>,
 * but we consider this an outer product for completeness of these functions
 * and since this is sometimes needed when writing templates that depend on
 * the rank of a tensor, which may sometimes be zero (i.e. a scalar).
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2000
 */
template <int dim, typename Number>
void outer_product (Tensor<1,dim,Number>       &dst,
                    const Number                src1,
                    const Tensor<1,dim,Number> &src2)
{
  for (unsigned int i=0; i<dim; ++i)
    dst[i] = src1 * src2[i];
}



/**
 * Form the outer product of two tensors of rank 1 and 0, i.e. <tt>dst[i] =
 * src1[i] * src2</tt>. Of course, this is only a scaling of <tt>src1</tt>,
 * but we consider this an outer product for completeness of these functions
 * and since this is sometimes needed when writing templates that depend on
 * the rank of a tensor, which may sometimes be zero (i.e. a scalar).
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2000
 */
template <int dim, typename Number>
void outer_product (Tensor<1,dim,Number>       &dst,
                    const Tensor<1,dim,Number>  src1,
                    const Number         src2)
{
  for (unsigned int i=0; i<dim; ++i)
    dst[i] = src1[i] * src2;
}


/**
 * Cross-product in 2d. This is just a rotation by 90 degrees clockwise to
 * compute the outer normal from a tangential vector. This function is defined
 * for all space dimensions to allow for dimension independent programming
 * (e.g. within switches over the space dimension), but may only be called if
 * the actual dimension of the arguments is two (e.g. from the <tt>dim==2</tt>
 * case in the switch).
 *
 * @relates Tensor
 * @author Guido Kanschat, 2001
 */
template <int dim, typename Number>
inline
void
cross_product (Tensor<1,dim,Number>       &dst,
               const Tensor<1,dim,Number> &src)
{
  Assert (dim==2, ExcInternalError());

  dst[0] = src[1];
  dst[1] = -src[0];
}


/**
 * Cross-product of 2 vectors in 3d. This function is defined for all space
 * dimensions to allow for dimension independent programming (e.g. within
 * switches over the space dimension), but may only be called if the actual
 * dimension of the arguments is three (e.g. from the <tt>dim==3</tt> case in
 * the switch).
 *
 * @relates Tensor
 * @author Guido Kanschat, 2001
 */
template <int dim, typename Number>
inline
void
cross_product (Tensor<1,dim,Number>       &dst,
               const Tensor<1,dim,Number> &src1,
               const Tensor<1,dim,Number> &src2)
{
  Assert (dim==3, ExcInternalError());

  dst[0] = src1[1]*src2[2] - src1[2]*src2[1];
  dst[1] = src1[2]*src2[0] - src1[0]*src2[2];
  dst[2] = src1[0]*src2[1] - src1[1]*src2[0];
}


/**
 * Compute the scalar product $a:b=\sum_{i,j} a_{ij}b_{ij}$ between two
 * tensors $a,b$ of rank 2. We don't use <code>operator*</code> for this
 * operation since the product between two tensors is usually assumed to be
 * the contraction over the last index of the first tensor and the first index
 * of the second tensor, for example $(a\cdot b)_{ij}=\sum_k a_{ik}b_{kj}$.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2008
 */
template <int dim, typename Number>
inline
Number
scalar_product (const Tensor<2,dim,Number> &t1,
                const Tensor<2,dim,Number> &t2)
{
  Number s = 0;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      s += t1[i][j] * t2[i][j];
  return s;
}


/**
 * Compute the determinant of a tensor of arbitrary rank and dimension one.
 * Since this is a number, the return value is, of course, the number itself.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 1998
 */
template <int rank, typename Number>
inline
Number determinant (const Tensor<rank,1,Number> &t)
{
  return determinant(t[0]);
}



/**
 * Compute the determinant of a tensor of rank one and dimension one. Since
 * this is a number, the return value is, of course, the number itself.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 1998
 */
template <typename Number>
inline
Number determinant (const Tensor<1,1,Number> &t)
{
  return t[0];
}


/**
 * Compute the determinant of a tensor of rank two and dimension one. Since
 * this is a number, the return value is, of course, the number itself.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 1998
 */
template <typename Number>
inline
Number determinant (const Tensor<2,1,Number> &t)
{
  return t[0][0];
}



/**
 * Compute the determinant of a tensor or rank 2, here for <tt>dim==2</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 1998
 */
template <typename Number>
inline
Number determinant (const Tensor<2,2,Number> &t)
{
  return ((t[0][0] * t[1][1]) - (t[1][0] * t[0][1]));
}


/**
 * Compute the determinant of a tensor or rank 2, here for <tt>dim==3</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 1998
 */
template <typename Number>
inline
Number determinant (const Tensor<2,3,Number> &t)
{
  // use exactly the same expression with the
  // same order of operations as for the inverse
  // to let the compiler use common
  // subexpression elimination when using
  // determinant and inverse in nearby code
  const Number t4 = t[0][0]*t[1][1],
               t6 = t[0][0]*t[1][2],
               t8 = t[0][1]*t[1][0],
               t00 = t[0][2]*t[1][0],
               t01 = t[0][1]*t[2][0],
               t04 = t[0][2]*t[2][0],
               det = (t4*t[2][2]-t6*t[2][1]-t8*t[2][2]+
                      t00*t[2][1]+t01*t[1][2]-t04*t[1][1]);
  return det;
}


/**
 * Compute the determinant of a tensor or rank 2, here for all dimensions for
 * which no explicit specialization is available above.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2009
 */
template <int dim, typename Number>
inline
Number determinant (const Tensor<2,dim,Number> &t)
{
  // compute the determinant using the
  // Laplace expansion of the
  // determinant. this may not be the most
  // efficient algorithm, but it does for
  // small n.
  //
  // for some algorithmic simplicity, we use
  // the expansion along the last row
  Number det = 0;

  for (unsigned int k=0; k<dim; ++k)
    {
      Tensor<2,dim-1,Number> minor;
      for (unsigned int i=0; i<dim-1; ++i)
        for (unsigned int j=0; j<dim-1; ++j)
          minor[i][j] = t[i][j<k ? j : j+1];

      const Number cofactor = std::pow (-1., static_cast<Number>(k+1)) *
                              determinant (minor);

      det += t[dim-1][k] * cofactor;
    }

  return std::pow (-1., static_cast<Number>(dim)) * det;
}



/**
 * Compute and return the trace of a tensor of rank 2, i.e. the sum of its
 * diagonal entries.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2001
 */
template <int dim, typename Number>
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
 * @relates Tensor
 * @author Wolfgang Bangerth, 2000
 */
template <int dim, typename Number>
inline
Tensor<2,dim,Number>
invert (const Tensor<2,dim,Number> &t)
{
  Number return_tensor [dim][dim];
  switch (dim)
    {
    case 1:
      return_tensor[0][0] = 1.0/t[0][0];
      break;

    case 2:
      // this is Maple output,
      // thus a bit unstructured
    {
      const Number det = t[0][0]*t[1][1]-t[1][0]*t[0][1];
      const Number t4 = 1.0/det;
      return_tensor[0][0] = t[1][1]*t4;
      return_tensor[0][1] = -t[0][1]*t4;
      return_tensor[1][0] = -t[1][0]*t4;
      return_tensor[1][1] = t[0][0]*t4;
      break;
    }

    case 3:
    {
      const Number t4 = t[0][0]*t[1][1],
                   t6 = t[0][0]*t[1][2],
                   t8 = t[0][1]*t[1][0],
                   t00 = t[0][2]*t[1][0],
                   t01 = t[0][1]*t[2][0],
                   t04 = t[0][2]*t[2][0],
                   det = (t4*t[2][2]-t6*t[2][1]-t8*t[2][2]+
                          t00*t[2][1]+t01*t[1][2]-t04*t[1][1]),
                         t07 = 1.0/det;
      return_tensor[0][0] = (t[1][1]*t[2][2]-t[1][2]*t[2][1])*t07;
      return_tensor[0][1] = (t[0][2]*t[2][1]-t[0][1]*t[2][2])*t07;
      return_tensor[0][2] = (t[0][1]*t[1][2]-t[0][2]*t[1][1])*t07;
      return_tensor[1][0] = (t[1][2]*t[2][0]-t[1][0]*t[2][2])*t07;
      return_tensor[1][1] = (t[0][0]*t[2][2]-t04)*t07;
      return_tensor[1][2] = (t00-t6)*t07;
      return_tensor[2][0] = (t[1][0]*t[2][1]-t[1][1]*t[2][0])*t07;
      return_tensor[2][1] = (t01-t[0][0]*t[2][1])*t07;
      return_tensor[2][2] = (t4-t8)*t07;

      break;
    }

    // if desired, take over the
    // inversion of a 4x4 tensor
    // from the FullMatrix
    default:
      AssertThrow (false, ExcNotImplemented());
    }
  return Tensor<2,dim,Number>(return_tensor);
}



/**
 * Return the transpose of the given tensor. Since the compiler can perform
 * the return value optimization, and since the size of the return object is
 * known, it is acceptable to return the result by value, rather than by
 * reference as a parameter. Note that there are specializations of this
 * function for <tt>dim==1,2,3</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2002
 */
template <int dim, typename Number>
inline
Tensor<2,dim,Number>
transpose (const Tensor<2,dim,Number> &t)
{
  Number tt[dim][dim];
  for (unsigned int i=0; i<dim; ++i)
    {
      tt[i][i] = t[i][i];
      for (unsigned int j=i+1; j<dim; ++j)
        {
          tt[i][j] = t[j][i];
          tt[j][i] = t[i][j];
        };
    }
  return Tensor<2,dim,Number>(tt);
}

#ifndef DOXYGEN

/**
 * Return the transpose of the given tensor. This is the specialization of the
 * general template for <tt>dim==1</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2002
 */
template <typename Number>
inline
Tensor<2,1,Number>
transpose (const Tensor<2,1,Number> &t)
{
  return t;
}




/**
 * Return the transpose of the given tensor. This is the specialization of the
 * general template for <tt>dim==2</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2002
 */
template <typename Number>
inline
Tensor<2,2,Number>
transpose (const Tensor<2,2,Number> &t)
{
  const Number x[2][2] = {{t[0][0], t[1][0]}, {t[0][1], t[1][1]}};
  return Tensor<2,2,Number>(x);
}




/**
 * Return the transpose of the given tensor. This is the specialization of the
 * general template for <tt>dim==3</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2002
 */
template <typename Number>
inline
Tensor<2,3,Number>
transpose (const Tensor<2,3,Number> &t)
{
  const Number x[3][3] = {{t[0][0], t[1][0], t[2][0]},
    {t[0][1], t[1][1], t[2][1]},
    {t[0][2], t[1][2], t[2][2]}
  };
  return Tensor<2,3,Number>(x);
}

#endif // DOXYGEN


/**
 * Return the $l_1$ norm of the given rank-2 tensor, where $||t||_1 = \max_j
 * \sum_i |t_{ij}|$ (maximum of the sums over columns).
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2012
 */
template <int dim, typename Number>
inline
double
l1_norm (const Tensor<2,dim,Number> &t)
{
  double max = 0;
  for (unsigned int j=0; j<dim; ++j)
    {
      double sum = 0;
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
 * @relates Tensor
 * @author Wolfgang Bangerth, 2012
 */
template <int dim, typename Number>
inline
double
linfty_norm (const Tensor<2,dim,Number> &t)
{
  double max = 0;
  for (unsigned int i=0; i<dim; ++i)
    {
      double sum = 0;
      for (unsigned int j=0; j<dim; ++j)
        sum += std::fabs(t[i][j]);

      if (sum > max)
        max = sum;
    }

  return max;
}



DEAL_II_NAMESPACE_CLOSE

#endif
