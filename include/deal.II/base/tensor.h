// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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

#ifndef dealii_tensor_h
#define dealii_tensor_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/tensor_accessors.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/lapack_full_matrix.h>

#ifdef DEAL_II_WITH_ADOLC
#  include <adolc/adouble.h> // Taped double
#endif

#include <cmath>
#include <ostream>
#include <utility>
#include <vector>


DEAL_II_NAMESPACE_OPEN

// Forward declarations:
#ifndef DOXYGEN
template <int dim, typename Number>
class Point;
template <int rank_, int dim, typename Number = double>
class Tensor;
template <typename Number>
class Vector;
template <typename number>
class FullMatrix;
namespace Differentiation
{
  namespace SD
  {
    class Expression;
  }
} // namespace Differentiation
#endif


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
class Tensor<0, dim, Number>
{
public:
  static_assert(dim >= 0,
                "Tensors must have a dimension greater than or equal to one.");

  /**
   * Provide a way to get the dimension of an object without explicit
   * knowledge of it's data type. Implementation is this way instead of
   * providing a function <tt>dimension()</tt> because now it is possible to
   * get the dimension at compile time without the expansion and preevaluation
   * of an inlined function; the compiler may therefore produce more efficient
   * code and you may use this value to declare other data types.
   */
  static constexpr unsigned int dimension = dim;

  /**
   * Publish the rank of this tensor to the outside world.
   */
  static constexpr unsigned int rank = 0;

  /**
   * Number of independent components of a tensor of rank 0.
   */
  static constexpr unsigned int n_independent_components = 1;

  /**
   * Declare a type that has holds real-valued numbers with the same precision
   * as the template argument to this class. For std::complex<number>, this
   * corresponds to type number, and it is equal to Number for all other
   * cases. See also the respective field in Vector<Number>.
   *
   * This alias is used to represent the return type of norms.
   */
  using real_type = typename numbers::NumberTraits<Number>::real_type;

  /**
   * Type of objects encapsulated by this container and returned by
   * operator[](). This is a scalar number type for a rank 0 tensor.
   */
  using value_type = Number;

  /**
   * Declare an array type which can be used to initialize an object of this
   * type statically. In case of a tensor of rank 0 this is just the scalar
   * number type Number.
   */
  using array_type = Number;

  /**
   * Constructor. Set to zero.
   *
   * @note This function can also be used in CUDA device code.
   */
  constexpr DEAL_II_CUDA_HOST_DEV
  Tensor();

  /**
   * Constructor from tensors with different underlying scalar type. This
   * obviously requires that the @p OtherNumber type is convertible to @p
   * Number.
   *
   * @note This function can also be used in CUDA device code.
   */
  template <typename OtherNumber>
  constexpr DEAL_II_CUDA_HOST_DEV
  Tensor(const Tensor<0, dim, OtherNumber> &initializer);

  /**
   * Constructor, where the data is copied from a C-style array.
   *
   * @note This function can also be used in CUDA device code.
   */
  template <typename OtherNumber>
  constexpr DEAL_II_CUDA_HOST_DEV
  Tensor(const OtherNumber &initializer);

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
   * @note This function can also be used in CUDA device code.
   */
  DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV operator Number &();

  /**
   * Return a reference to the encapsulated Number object. Since rank-0
   * tensors are scalars, this is a natural operation.
   *
   * This is the const conversion operator that returns a read-only reference.
   *
   * @note This function can also be used in CUDA device code.
   */
  DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV operator const Number &() const;

  /**
   * Assignment from tensors with different underlying scalar type. This
   * obviously requires that the @p OtherNumber type is convertible to @p
   * Number.
   *
   * @note This function can also be used in CUDA device code.
   */
  template <typename OtherNumber>
  DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV Tensor &
                                          operator=(const Tensor<0, dim, OtherNumber> &rhs);

#ifdef __INTEL_COMPILER
  /**
   * Assignment from tensors with same underlying scalar type.
   * This is needed for ICC15 because it can't generate a suitable
   * copy constructor for Sacado::Rad::ADvar types automatically.
   * See https://github.com/dealii/dealii/pull/5865.
   *
   * @note This function can also be used in CUDA device code.
   */
  DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV Tensor &
                                          operator=(const Tensor<0, dim, Number> &rhs);
#endif

  /**
   * This operator assigns a scalar to a tensor. This obviously requires
   * that the @p OtherNumber type is convertible to @p Number.
   *
   * @note This function can also be used in CUDA device code.
   */
  template <typename OtherNumber>
  DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV Tensor &
                                          operator=(const OtherNumber &d);

  /**
   * Test for equality of two tensors.
   */
  template <typename OtherNumber>
  DEAL_II_CONSTEXPR bool
  operator==(const Tensor<0, dim, OtherNumber> &rhs) const;

  /**
   * Test for inequality of two tensors.
   */
  template <typename OtherNumber>
  constexpr bool
  operator!=(const Tensor<0, dim, OtherNumber> &rhs) const;

  /**
   * Add another scalar.
   *
   * @note This function can also be used in CUDA device code.
   */
  template <typename OtherNumber>
  DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV Tensor &
                                          operator+=(const Tensor<0, dim, OtherNumber> &rhs);

  /**
   * Subtract another scalar.
   *
   * @note This function can also be used in CUDA device code.
   */
  template <typename OtherNumber>
  DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV Tensor &
                                          operator-=(const Tensor<0, dim, OtherNumber> &rhs);

  /**
   * Multiply the scalar with a <tt>factor</tt>.
   *
   * @note This function can also be used in CUDA device code.
   */
  template <typename OtherNumber>
  DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV Tensor &
                                          operator*=(const OtherNumber &factor);

  /**
   * Divide the scalar by <tt>factor</tt>.
   *
   * @note This function can also be used in CUDA device code.
   */
  template <typename OtherNumber>
  DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV Tensor &
                                          operator/=(const OtherNumber &factor);

  /**
   * Tensor with inverted entries.
   *
   * @note This function can also be used in CUDA device code.
   */
  constexpr DEAL_II_CUDA_HOST_DEV Tensor
                                  operator-() const;

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
  DEAL_II_CONSTEXPR void
  clear();

  /**
   * Return the Frobenius-norm of a tensor, i.e. the square root of the sum of
   * the absolute squares of all entries. For the present case of rank-1
   * tensors, this equals the usual <tt>l<sub>2</sub></tt> norm of the vector.
   */
  real_type
  norm() const;

  /**
   * Return the square of the Frobenius-norm of a tensor, i.e. the sum of the
   * absolute squares of all entries.
   *
   * @note This function can also be used in CUDA device code.
   */
  DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV real_type
                                          norm_square() const;

  /**
   * Read or write the data of this object to or from a stream for the purpose
   * of serialization
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version);

  /**
   * Internal type declaration that is used to specialize the return type of
   * operator[]() for Tensor<1,dim,Number>
   */
  using tensor_type = Number;

private:
  /**
   * The value of this scalar object.
   */
  Number value;

  /**
   * Internal helper function for unroll.
   */
  template <typename OtherNumber>
  void
  unroll_recursion(Vector<OtherNumber> &result,
                   unsigned int &       start_index) const;

  // Allow an arbitrary Tensor to access the underlying values.
  template <int, int, typename>
  friend class Tensor;
};



/**
 * A general tensor class with an arbitrary rank, i.e. with an arbitrary
 * number of indices. The Tensor class provides an indexing operator and a bit
 * of infrastructure, but most functionality is recursively handed down to
 * tensors of rank 1 or put into external templated functions, e.g. the
 * <tt>contract</tt> family.
 *
 * The rank of a tensor specifies which types of physical quantities it can
 * represent:
 * <ul>
 *   <li> A rank-0 tensor is a scalar that can store quantities such as
 *     temperature or pressure. These scalar quantities are shown in this
 *     documentation as simple lower-case Latin letters e.g. $a, b, c, \dots$.
 *   </li>
 *   <li> A rank-1 tensor is a vector with @p dim components and it can
 *     represent vector quantities such as velocity, displacement, electric
 *     field, etc. They can also describe the gradient of a scalar field.
 *     The notation used for rank-1 tensors is bold-faced lower-case Latin
 *     letters e.g. $\mathbf a, \mathbf b, \mathbf c, \dots$.
 *     The components of a rank-1 tensor such as $\mathbf a$ are represented
 *     as $a_i$ where $i$ is an index between 0 and <tt>dim-1</tt>.
 *   </li>
 *   <li> A rank-2 tensor is a linear operator that can transform a vector
 *     into another vector. These tensors are similar to matrices with
 *     $\text{dim} \times \text{dim}$ components. There is a related class
 *     SymmetricTensor<2,dim> for tensors of rank 2 whose elements are
 *     symmetric. Rank-2 tensors are usually denoted by bold-faced upper-case
 *     Latin letters such as $\mathbf A, \mathbf B, \dots$ or bold-faced Greek
 *     letters for example $\boldsymbol{\varepsilon}, \boldsymbol{\sigma}$.
 *     The components of a rank 2 tensor such as $\mathbf A$ are shown with
 *     two indices $(i,j)$ as $A_{ij}$. These tensors usually describe the
 *     gradients of vector fields (deformation gradient, velocity gradient,
 *     etc.) or Hessians of scalar fields. Additionally, mechanical stress
 *     tensors are rank-2 tensors that map the unit normal vectors of internal
 *     surfaces into local traction (force per unit area) vectors.
 *   </li>
 *   <li> Tensors with ranks higher than 2 are similarly defined in a
 *     consistent manner. They have $\text{dim}^{\text{rank}}$ components and
 *     the number of indices required to identify a component equals
 *     <tt>rank</tt>. For rank-4 tensors, a symmetric variant called
 *     SymmetricTensor<4,dim> exists.
 *   </li>
 * </ul>
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
 * @tparam rank_ An integer that denotes the rank of this tensor. A
 * specialization of this class exists for rank-0 tensors.
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
  static_assert(rank_ >= 0,
                "Tensors must have a rank greater than or equal to one.");
  static_assert(dim >= 0,
                "Tensors must have a dimension greater than or equal to one.");
  /**
   * Provide a way to get the dimension of an object without explicit
   * knowledge of it's data type. Implementation is this way instead of
   * providing a function <tt>dimension()</tt> because now it is possible to
   * get the dimension at compile time without the expansion and preevaluation
   * of an inlined function; the compiler may therefore produce more efficient
   * code and you may use this value to declare other data types.
   */
  static constexpr unsigned int dimension = dim;

  /**
   * Publish the rank of this tensor to the outside world.
   */
  static constexpr unsigned int rank = rank_;

  /**
   * Number of independent components of a tensor of current rank. This is dim
   * times the number of independent components of each sub-tensor.
   */
  static constexpr unsigned int n_independent_components =
    Tensor<rank_ - 1, dim>::n_independent_components * dim;

  /**
   * Type of objects encapsulated by this container and returned by
   * operator[](). This is a tensor of lower rank for a general tensor, and a
   * scalar number type for Tensor<1,dim,Number>.
   */
  using value_type = typename Tensor<rank_ - 1, dim, Number>::tensor_type;

  /**
   * Declare an array type which can be used to initialize an object of this
   * type statically. For `dim == 0`, its size is 1. Otherwise, it is `dim`.
   */
  using array_type =
    typename Tensor<rank_ - 1, dim, Number>::array_type[(dim != 0) ? dim : 1];

  /**
   * Constructor. Initialize all entries to zero.
   *
   * @note This function can also be used in CUDA device code.
   */
  constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV
                                  Tensor()
#ifdef DEAL_II_MSVC
    : values{}
  {}
#else
    = default;
#endif

  /**
   * Constructor, where the data is copied from a C-style array.
   *
   * @note This function can also be used in CUDA device code.
   */
  constexpr DEAL_II_CUDA_HOST_DEV explicit Tensor(
    const array_type &initializer);

  /**
   * Constructor from tensors with different underlying scalar type. This
   * obviously requires that the @p OtherNumber type is convertible to @p
   * Number.
   *
   * @note This function can also be used in CUDA device code.
   */
  template <typename OtherNumber>
  constexpr DEAL_II_CUDA_HOST_DEV
  Tensor(const Tensor<rank_, dim, OtherNumber> &initializer);

  /**
   * Constructor that converts from a "tensor of tensors".
   */
  template <typename OtherNumber>
  constexpr Tensor(
    const Tensor<1, dim, Tensor<rank_ - 1, dim, OtherNumber>> &initializer);

  /**
   * Conversion operator to tensor of tensors.
   */
  template <typename OtherNumber>
  constexpr
  operator Tensor<1, dim, Tensor<rank_ - 1, dim, OtherNumber>>() const;

  /**
   * Read-Write access operator.
   *
   * @note This function can also be used in CUDA device code.
   */
  DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV value_type &
                                          operator[](const unsigned int i);

  /**
   * Read-only access operator.
   *
   * @note This function can also be used in CUDA device code.
   */
  constexpr DEAL_II_CUDA_HOST_DEV const value_type &
                                        operator[](const unsigned int i) const;

  /**
   * Read access using TableIndices <tt>indices</tt>
   */
  DEAL_II_CONSTEXPR const Number &
                          operator[](const TableIndices<rank_> &indices) const;

  /**
   * Read and write access using TableIndices <tt>indices</tt>
   */
  DEAL_II_CONSTEXPR Number &operator[](const TableIndices<rank_> &indices);

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
   *
   * @note This function can also be used in CUDA device code.
   */
  template <typename OtherNumber>
  DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV Tensor &
                                          operator=(const Tensor<rank_, dim, OtherNumber> &rhs);

  /**
   * This operator assigns a scalar to a tensor. To avoid confusion with what
   * exactly it means to assign a scalar value to a tensor, zero is the only
   * value allowed for <tt>d</tt>, allowing the intuitive notation
   * <tt>t=0</tt> to reset all elements of the tensor to zero.
   */
  DEAL_II_CONSTEXPR Tensor &
                    operator=(const Number &d);

  /**
   * Test for equality of two tensors.
   */
  template <typename OtherNumber>
  DEAL_II_CONSTEXPR bool
  operator==(const Tensor<rank_, dim, OtherNumber> &) const;

  /**
   * Test for inequality of two tensors.
   */
  template <typename OtherNumber>
  constexpr bool
  operator!=(const Tensor<rank_, dim, OtherNumber> &) const;

  /**
   * Add another tensor.
   *
   * @note This function can also be used in CUDA device code.
   */
  template <typename OtherNumber>
  DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV Tensor &
                                          operator+=(const Tensor<rank_, dim, OtherNumber> &);

  /**
   * Subtract another tensor.
   *
   * @note This function can also be used in CUDA device code.
   */
  template <typename OtherNumber>
  DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV Tensor &
                                          operator-=(const Tensor<rank_, dim, OtherNumber> &);

  /**
   * Scale the tensor by <tt>factor</tt>, i.e. multiply all components by
   * <tt>factor</tt>.
   *
   * @note This function can also be used in CUDA device code.
   */
  template <typename OtherNumber>
  DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV Tensor &
                                          operator*=(const OtherNumber &factor);

  /**
   * Scale the vector by <tt>1/factor</tt>.
   *
   * @note This function can also be used in CUDA device code.
   */
  template <typename OtherNumber>
  DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV Tensor &
                                          operator/=(const OtherNumber &factor);

  /**
   * Unary minus operator. Negate all entries of a tensor.
   *
   * @note This function can also be used in CUDA device code.
   */
  DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV Tensor
                                          operator-() const;

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
  DEAL_II_CONSTEXPR void
  clear();

  /**
   * Return the Frobenius-norm of a tensor, i.e. the square root of the sum of
   * the absolute squares of all entries. For the present case of rank-1
   * tensors, this equals the usual <tt>l<sub>2</sub></tt> norm of the vector.
   *
   * @note This function can also be used in CUDA device code.
   */
  DEAL_II_CUDA_HOST_DEV
  typename numbers::NumberTraits<Number>::real_type
  norm() const;

  /**
   * Return the square of the Frobenius-norm of a tensor, i.e. the sum of the
   * absolute squares of all entries.
   *
   * @note This function can also be used in CUDA device code.
   */
  DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV
    typename numbers::NumberTraits<Number>::real_type
    norm_square() const;

  /**
   * Fill a vector with all tensor elements.
   *
   * This function unrolls all tensor entries into a single, linearly numbered
   * vector. As usual in C++, the rightmost index of the tensor marches
   * fastest.
   */
  template <typename OtherNumber>
  void
  unroll(Vector<OtherNumber> &result) const;

  /**
   * Return an unrolled index in the range $[0,\text{dim}^{\text{rank}}-1]$
   * for the element of the tensor indexed by the argument to the function.
   */
  static DEAL_II_CONSTEXPR unsigned int
  component_to_unrolled_index(const TableIndices<rank_> &indices);

  /**
   * Opposite of  component_to_unrolled_index: For an index in the range
   * $[0, \text{dim}^{\text{rank}}-1]$, return which set of indices it would
   * correspond to.
   */
  static DEAL_II_CONSTEXPR TableIndices<rank_>
                           unrolled_to_component_indices(const unsigned int i);

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  static constexpr std::size_t
  memory_consumption();

  /**
   * Read or write the data of this object to or from a stream for the purpose
   * of serialization
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version);

  /**
   * Internal type declaration that is used to specialize the return type of
   * operator[]() for Tensor<1,dim,Number>
   */
  using tensor_type = Tensor<rank_, dim, Number>;

private:
  /**
   * Array of tensors holding the subelements.
   */
  Tensor<rank_ - 1, dim, Number> values[(dim != 0) ? dim : 1];
  // ... avoid a compiler warning in case of dim == 0 and ensure that the
  // array always has positive size.

  /**
   * Internal helper function for unroll.
   */
  template <typename OtherNumber>
  void
  unroll_recursion(Vector<OtherNumber> &result,
                   unsigned int &       start_index) const;

  /**
   * This constructor is for internal use. It provides a way
   * to create constexpr constructors for Tensor<rank, dim, Number>
   *
   * @note This function can also be used in CUDA device code.
   */
  template <typename ArrayLike, std::size_t... Indices>
  constexpr DEAL_II_CUDA_HOST_DEV
  Tensor(const ArrayLike &initializer, std::index_sequence<Indices...>);

  // Allow an arbitrary Tensor to access the underlying values.
  template <int, int, typename>
  friend class Tensor;

  // Point is allowed access to the coordinates. This is supposed to improve
  // speed.
  friend class Point<dim, Number>;
};


#ifndef DOXYGEN
namespace internal
{
  // Workaround: The following 4 overloads are necessary to be able to
  // compile the library with Apple Clang 8 and older. We should remove
  // these overloads again when we bump the minimal required version to
  // something later than clang-3.6 / Apple Clang 6.3.
  // - Jean-Paul Pelteret, Matthias Maier, Daniel Arndt 2020
  template <int rank, int dim, typename T, typename U>
  struct ProductTypeImpl<Tensor<rank, dim, T>, std::complex<U>>
  {
    using type =
      Tensor<rank, dim, std::complex<typename ProductType<T, U>::type>>;
  };

  template <int rank, int dim, typename T, typename U>
  struct ProductTypeImpl<Tensor<rank, dim, std::complex<T>>, std::complex<U>>
  {
    using type =
      Tensor<rank, dim, std::complex<typename ProductType<T, U>::type>>;
  };

  template <typename T, int rank, int dim, typename U>
  struct ProductTypeImpl<std::complex<T>, Tensor<rank, dim, U>>
  {
    using type =
      Tensor<rank, dim, std::complex<typename ProductType<T, U>::type>>;
  };

  template <int rank, int dim, typename T, typename U>
  struct ProductTypeImpl<std::complex<T>, Tensor<rank, dim, std::complex<U>>>
  {
    using type =
      Tensor<rank, dim, std::complex<typename ProductType<T, U>::type>>;
  };
  // end workaround

  /**
   * The structs below are needed to initialize nested Tensor objects.
   * Also see numbers.h for another specialization.
   */
  template <int rank, int dim, typename T>
  struct NumberType<Tensor<rank, dim, T>>
  {
    static constexpr DEAL_II_ALWAYS_INLINE const Tensor<rank, dim, T> &
                                                 value(const Tensor<rank, dim, T> &t)
    {
      return t;
    }

    static DEAL_II_CONSTEXPR DEAL_II_ALWAYS_INLINE Tensor<rank, dim, T>
                                                   value(const T &t)
    {
      Tensor<rank, dim, T> tmp;
      tmp = t;
      return tmp;
    }
  };
} // namespace internal


/*---------------------- Inline functions: Tensor<0,dim> ---------------------*/


template <int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV
                                Tensor<0, dim, Number>::Tensor()
  // Some auto-differentiable numbers need explicit
  // zero initialization such as adtl::adouble.
  : Tensor{0.0}
{}



template <int dim, typename Number>
template <typename OtherNumber>
constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV
                                Tensor<0, dim, Number>::Tensor(const OtherNumber &initializer)
  : value(internal::NumberType<Number>::value(initializer))
{}



template <int dim, typename Number>
template <typename OtherNumber>
constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV
                                Tensor<0, dim, Number>::Tensor(const Tensor<0, dim, OtherNumber> &p)
  : Tensor{p.value}
{}



template <int dim, typename Number>
inline Number *
Tensor<0, dim, Number>::begin_raw()
{
  return std::addressof(value);
}



template <int dim, typename Number>
inline const Number *
Tensor<0, dim, Number>::begin_raw() const
{
  return std::addressof(value);
}



template <int dim, typename Number>
inline Number *
Tensor<0, dim, Number>::end_raw()
{
  return begin_raw() + n_independent_components;
}



template <int dim, typename Number>
const Number *
Tensor<0, dim, Number>::end_raw() const
{
  return begin_raw() + n_independent_components;
}



template <int dim, typename Number>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
  DEAL_II_CUDA_HOST_DEV Tensor<0, dim, Number>::operator Number &()
{
  // We cannot use Assert inside a CUDA kernel
#  ifndef __CUDA_ARCH__
  Assert(dim != 0,
         ExcMessage("Cannot access an object of type Tensor<0,0,Number>"));
#  endif
  return value;
}


template <int dim, typename Number>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
  DEAL_II_CUDA_HOST_DEV Tensor<0, dim, Number>::operator const Number &() const
{
  // We cannot use Assert inside a CUDA kernel
#  ifndef __CUDA_ARCH__
  Assert(dim != 0,
         ExcMessage("Cannot access an object of type Tensor<0,0,Number>"));
#  endif
  return value;
}


template <int dim, typename Number>
template <typename OtherNumber>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
  DEAL_II_CUDA_HOST_DEV Tensor<0, dim, Number> &
  Tensor<0, dim, Number>::operator=(const Tensor<0, dim, OtherNumber> &p)
{
  value = internal::NumberType<Number>::value(p);
  return *this;
}


#  ifdef __INTEL_COMPILER
template <int dim, typename Number>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
  DEAL_II_CUDA_HOST_DEV Tensor<0, dim, Number> &
  Tensor<0, dim, Number>::operator=(const Tensor<0, dim, Number> &p)
{
  value = p.value;
  return *this;
}
#  endif


template <int dim, typename Number>
template <typename OtherNumber>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
  DEAL_II_CUDA_HOST_DEV Tensor<0, dim, Number> &
  Tensor<0, dim, Number>::operator=(const OtherNumber &d)
{
  value = internal::NumberType<Number>::value(d);
  return *this;
}


template <int dim, typename Number>
template <typename OtherNumber>
DEAL_II_CONSTEXPR inline bool
Tensor<0, dim, Number>::operator==(const Tensor<0, dim, OtherNumber> &p) const
{
#  if defined(DEAL_II_ADOLC_WITH_ADVANCED_BRANCHING)
  Assert(!(std::is_same<Number, adouble>::value ||
           std::is_same<OtherNumber, adouble>::value),
         ExcMessage(
           "The Tensor equality operator for ADOL-C taped numbers has not yet "
           "been extended to support advanced branching."));
#  endif

  return numbers::values_are_equal(value, p.value);
}


template <int dim, typename Number>
template <typename OtherNumber>
constexpr bool
Tensor<0, dim, Number>::operator!=(const Tensor<0, dim, OtherNumber> &p) const
{
  return !((*this) == p);
}


template <int dim, typename Number>
template <typename OtherNumber>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
  DEAL_II_CUDA_HOST_DEV Tensor<0, dim, Number> &
  Tensor<0, dim, Number>::operator+=(const Tensor<0, dim, OtherNumber> &p)
{
  value += p.value;
  return *this;
}


template <int dim, typename Number>
template <typename OtherNumber>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
  DEAL_II_CUDA_HOST_DEV Tensor<0, dim, Number> &
  Tensor<0, dim, Number>::operator-=(const Tensor<0, dim, OtherNumber> &p)
{
  value -= p.value;
  return *this;
}



namespace internal
{
  namespace ComplexWorkaround
  {
    template <typename Number, typename OtherNumber>
    DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV void
                                                   multiply_assign_scalar(Number &val, const OtherNumber &s)
    {
      val *= s;
    }

#  ifdef __CUDA_ARCH__
    template <typename Number, typename OtherNumber>
    DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV void
                                                   multiply_assign_scalar(std::complex<Number> &, const OtherNumber &)
    {
      printf("This function is not implemented for std::complex<Number>!\n");
      assert(false);
    }
#  endif
  } // namespace ComplexWorkaround
} // namespace internal


template <int dim, typename Number>
template <typename OtherNumber>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
  DEAL_II_CUDA_HOST_DEV Tensor<0, dim, Number> &
  Tensor<0, dim, Number>::operator*=(const OtherNumber &s)
{
  internal::ComplexWorkaround::multiply_assign_scalar(value, s);
  return *this;
}



template <int dim, typename Number>
template <typename OtherNumber>
DEAL_II_CONSTEXPR inline DEAL_II_CUDA_HOST_DEV Tensor<0, dim, Number> &
Tensor<0, dim, Number>::operator/=(const OtherNumber &s)
{
  value /= s;
  return *this;
}


template <int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV Tensor<0, dim, Number>
Tensor<0, dim, Number>::operator-() const
{
  return -value;
}


template <int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE typename Tensor<0, dim, Number>::real_type
Tensor<0, dim, Number>::norm() const
{
  Assert(dim != 0,
         ExcMessage("Cannot access an object of type Tensor<0,0,Number>"));
  return numbers::NumberTraits<Number>::abs(value);
}


template <int dim, typename Number>
DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV inline DEAL_II_ALWAYS_INLINE
  typename Tensor<0, dim, Number>::real_type
  Tensor<0, dim, Number>::norm_square() const
{
  // We cannot use Assert inside a CUDA kernel
#  ifndef __CUDA_ARCH__
  Assert(dim != 0,
         ExcMessage("Cannot access an object of type Tensor<0,0,Number>"));
#  endif
  return numbers::NumberTraits<Number>::abs_square(value);
}


template <int dim, typename Number>
template <typename OtherNumber>
inline void
Tensor<0, dim, Number>::unroll_recursion(Vector<OtherNumber> &result,
                                         unsigned int &       index) const
{
  Assert(dim != 0,
         ExcMessage("Cannot unroll an object of type Tensor<0,0,Number>"));
  result[index] = value;
  ++index;
}


template <int dim, typename Number>
DEAL_II_CONSTEXPR inline void
Tensor<0, dim, Number>::clear()
{
  // Some auto-differentiable numbers need explicit
  // zero initialization.
  value = internal::NumberType<Number>::value(0.0);
}


template <int dim, typename Number>
template <class Archive>
inline void
Tensor<0, dim, Number>::serialize(Archive &ar, const unsigned int)
{
  ar &value;
}


template <int dim, typename Number>
constexpr unsigned int Tensor<0, dim, Number>::n_independent_components;


/*-------------------- Inline functions: Tensor<rank,dim> --------------------*/

template <int rank_, int dim, typename Number>
template <typename ArrayLike, std::size_t... indices>
constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV
                                Tensor<rank_, dim, Number>::Tensor(const ArrayLike &initializer,
                                   std::index_sequence<indices...>)
  : values{Tensor<rank_ - 1, dim, Number>(initializer[indices])...}
{
  static_assert(sizeof...(indices) == dim,
                "dim should match the number of indices");
}


template <int rank_, int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV
                                Tensor<rank_, dim, Number>::Tensor(const array_type &initializer)
  : Tensor(initializer, std::make_index_sequence<dim>{})
{}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV
                                Tensor<rank_, dim, Number>::Tensor(
  const Tensor<rank_, dim, OtherNumber> &initializer)
  : Tensor(initializer, std::make_index_sequence<dim>{})
{}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
constexpr DEAL_II_ALWAYS_INLINE
Tensor<rank_, dim, Number>::Tensor(
  const Tensor<1, dim, Tensor<rank_ - 1, dim, OtherNumber>> &initializer)
  : Tensor(initializer, std::make_index_sequence<dim>{})
{}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
constexpr DEAL_II_ALWAYS_INLINE Tensor<rank_, dim, Number>::
                                operator Tensor<1, dim, Tensor<rank_ - 1, dim, OtherNumber>>() const
{
  return Tensor<1, dim, Tensor<rank_ - 1, dim, Number>>(values);
}



namespace internal
{
  namespace TensorSubscriptor
  {
    template <typename ArrayElementType, int dim>
    DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
      DEAL_II_CUDA_HOST_DEV ArrayElementType &
                            subscript(ArrayElementType * values,
                                      const unsigned int i,
                                      std::integral_constant<int, dim>)
    {
      // We cannot use Assert in a CUDA kernel
#  ifndef __CUDA_ARCH__
      AssertIndexRange(i, dim);
#  endif
      return values[i];
    }

    // The variables within this struct will be referenced in the next function.
    // It is a workaround that allows returning a reference to a static variable
    // while allowing constexpr evaluation of the function.
    // It has to be defined outside the function because constexpr functions
    // cannot define static variables
    template <typename ArrayElementType>
    struct Uninitialized
    {
      static ArrayElementType value;
    };

    template <typename Type>
    Type Uninitialized<Type>::value;

    template <typename ArrayElementType>
    DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
      DEAL_II_CUDA_HOST_DEV ArrayElementType &
                            subscript(ArrayElementType *,
                                      const unsigned int,
                                      std::integral_constant<int, 0>)
    {
      // We cannot use Assert in a CUDA kernel
#  ifndef __CUDA_ARCH__
      Assert(
        false,
        ExcMessage(
          "Cannot access elements of an object of type Tensor<rank,0,Number>."));
#  endif
      return Uninitialized<ArrayElementType>::value;
    }
  } // namespace TensorSubscriptor
} // namespace internal


template <int rank_, int dim, typename Number>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE     DEAL_II_CUDA_HOST_DEV
  typename Tensor<rank_, dim, Number>::value_type &Tensor<rank_, dim, Number>::
                                                   operator[](const unsigned int i)
{
  return dealii::internal::TensorSubscriptor::subscript(
    values, i, std::integral_constant<int, dim>());
}


template <int rank_, int dim, typename Number>
constexpr DEAL_II_ALWAYS_INLINE
    DEAL_II_CUDA_HOST_DEV const typename Tensor<rank_, dim, Number>::value_type &
    Tensor<rank_, dim, Number>::operator[](const unsigned int i) const
{
  return values[i];
}


template <int rank_, int dim, typename Number>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE const Number &
                                                     Tensor<rank_, dim, Number>::
                                                     operator[](const TableIndices<rank_> &indices) const
{
  Assert(dim != 0,
         ExcMessage("Cannot access an object of type Tensor<rank_,0,Number>"));

  return TensorAccessors::extract<rank_>(*this, indices);
}



template <int rank_, int dim, typename Number>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE Number &
  Tensor<rank_, dim, Number>::operator[](const TableIndices<rank_> &indices)
{
  Assert(dim != 0,
         ExcMessage("Cannot access an object of type Tensor<rank_,0,Number>"));

  return TensorAccessors::extract<rank_>(*this, indices);
}



template <int rank_, int dim, typename Number>
inline Number *
Tensor<rank_, dim, Number>::begin_raw()
{
  return std::addressof(
    this->operator[](this->unrolled_to_component_indices(0)));
}



template <int rank_, int dim, typename Number>
inline const Number *
Tensor<rank_, dim, Number>::begin_raw() const
{
  return std::addressof(
    this->operator[](this->unrolled_to_component_indices(0)));
}



template <int rank_, int dim, typename Number>
inline Number *
Tensor<rank_, dim, Number>::end_raw()
{
  return begin_raw() + n_independent_components;
}



template <int rank_, int dim, typename Number>
inline const Number *
Tensor<rank_, dim, Number>::end_raw() const
{
  return begin_raw() + n_independent_components;
}



template <int rank_, int dim, typename Number>
template <typename OtherNumber>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE Tensor<rank_, dim, Number> &
Tensor<rank_, dim, Number>::operator=(const Tensor<rank_, dim, OtherNumber> &t)
{
  // The following loop could be written more concisely using std::copy, but
  // that function is only constexpr from C++20 on.
  for (unsigned int i = 0; i < dim; ++i)
    values[i] = t.values[i];
  return *this;
}


template <int rank_, int dim, typename Number>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE Tensor<rank_, dim, Number> &
Tensor<rank_, dim, Number>::operator=(const Number &d)
{
  Assert(numbers::value_is_zero(d),
         ExcMessage("Only assignment with zero is allowed"));
  (void)d;

  for (unsigned int i = 0; i < dim; ++i)
    values[i] = internal::NumberType<Number>::value(0.0);
  return *this;
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
DEAL_II_CONSTEXPR inline bool
Tensor<rank_, dim, Number>::
operator==(const Tensor<rank_, dim, OtherNumber> &p) const
{
  for (unsigned int i = 0; i < dim; ++i)
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
DEAL_II_CONSTEXPR inline bool
Tensor<1, 0, double>::operator==(const Tensor<1, 0, double> &) const
{
  return true;
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
constexpr bool
Tensor<rank_, dim, Number>::
operator!=(const Tensor<rank_, dim, OtherNumber> &p) const
{
  return !((*this) == p);
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
  DEAL_II_CUDA_HOST_DEV Tensor<rank_, dim, Number> &
                        Tensor<rank_, dim, Number>::
                        operator+=(const Tensor<rank_, dim, OtherNumber> &p)
{
  for (unsigned int i = 0; i < dim; ++i)
    values[i] += p.values[i];
  return *this;
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
  DEAL_II_CUDA_HOST_DEV Tensor<rank_, dim, Number> &
                        Tensor<rank_, dim, Number>::
                        operator-=(const Tensor<rank_, dim, OtherNumber> &p)
{
  for (unsigned int i = 0; i < dim; ++i)
    values[i] -= p.values[i];
  return *this;
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
  DEAL_II_CUDA_HOST_DEV Tensor<rank_, dim, Number> &
  Tensor<rank_, dim, Number>::operator*=(const OtherNumber &s)
{
  for (unsigned int i = 0; i < dim; ++i)
    values[i] *= s;
  return *this;
}


namespace internal
{
  namespace TensorImplementation
  {
    template <int rank,
              int dim,
              typename Number,
              typename OtherNumber,
              typename std::enable_if<
                !std::is_integral<
                  typename ProductType<Number, OtherNumber>::type>::value &&
                  !std::is_same<Number, Differentiation::SD::Expression>::value,
                int>::type = 0>
    DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV inline DEAL_II_ALWAYS_INLINE void
                      division_operator(Tensor<rank, dim, Number> (&t)[dim],
                                        const OtherNumber &factor)
    {
      const Number inverse_factor = Number(1.) / factor;
      // recurse over the base objects
      for (unsigned int d = 0; d < dim; ++d)
        t[d] *= inverse_factor;
    }


    template <int rank,
              int dim,
              typename Number,
              typename OtherNumber,
              typename std::enable_if<
                std::is_integral<
                  typename ProductType<Number, OtherNumber>::type>::value ||
                  std::is_same<Number, Differentiation::SD::Expression>::value,
                int>::type = 0>
    DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV inline DEAL_II_ALWAYS_INLINE void
                      division_operator(Tensor<rank, dim, Number> (&t)[dim],
                                        const OtherNumber &factor)
    {
      // recurse over the base objects
      for (unsigned int d = 0; d < dim; ++d)
        t[d] /= factor;
    }
  } // namespace TensorImplementation
} // namespace internal


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
  DEAL_II_CUDA_HOST_DEV Tensor<rank_, dim, Number> &
  Tensor<rank_, dim, Number>::operator/=(const OtherNumber &s)
{
  internal::TensorImplementation::division_operator(values, s);
  return *this;
}


template <int rank_, int dim, typename Number>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
  DEAL_II_CUDA_HOST_DEV Tensor<rank_, dim, Number>
  Tensor<rank_, dim, Number>::operator-() const
{
  Tensor<rank_, dim, Number> tmp;

  for (unsigned int i = 0; i < dim; ++i)
    tmp.values[i] = -values[i];

  return tmp;
}


template <int rank_, int dim, typename Number>
inline typename numbers::NumberTraits<Number>::real_type
Tensor<rank_, dim, Number>::norm() const
{
  return std::sqrt(norm_square());
}


template <int rank_, int dim, typename Number>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV
  typename numbers::NumberTraits<Number>::real_type
  Tensor<rank_, dim, Number>::norm_square() const
{
  typename numbers::NumberTraits<Number>::real_type s = internal::NumberType<
    typename numbers::NumberTraits<Number>::real_type>::value(0.0);
  for (unsigned int i = 0; i < dim; ++i)
    s += values[i].norm_square();

  return s;
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline void
Tensor<rank_, dim, Number>::unroll(Vector<OtherNumber> &result) const
{
  AssertDimension(result.size(),
                  (Utilities::fixed_power<rank_, unsigned int>(dim)));

  unsigned int index = 0;
  unroll_recursion(result, index);
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline void
Tensor<rank_, dim, Number>::unroll_recursion(Vector<OtherNumber> &result,
                                             unsigned int &       index) const
{
  for (unsigned int i = 0; i < dim; ++i)
    values[i].unroll_recursion(result, index);
}


template <int rank_, int dim, typename Number>
DEAL_II_CONSTEXPR inline unsigned int
Tensor<rank_, dim, Number>::component_to_unrolled_index(
  const TableIndices<rank_> &indices)
{
  unsigned int index = 0;
  for (int r = 0; r < rank_; ++r)
    index = index * dim + indices[r];

  return index;
}



namespace internal
{
  // unrolled_to_component_indices is instantiated from DataOut for dim==0
  // and rank=2. Make sure we don't have compiler warnings.

  template <int dim>
  inline DEAL_II_CONSTEXPR unsigned int
  mod(const unsigned int x)
  {
    return x % dim;
  }

  template <>
  inline unsigned int
  mod<0>(const unsigned int x)
  {
    Assert(false, ExcInternalError());
    return x;
  }

  template <int dim>
  inline DEAL_II_CONSTEXPR unsigned int
  div(const unsigned int x)
  {
    return x / dim;
  }

  template <>
  inline unsigned int
  div<0>(const unsigned int x)
  {
    Assert(false, ExcInternalError());
    return x;
  }

} // namespace internal



template <int rank_, int dim, typename Number>
DEAL_II_CONSTEXPR inline TableIndices<rank_>
Tensor<rank_, dim, Number>::unrolled_to_component_indices(const unsigned int i)
{
  AssertIndexRange(i, n_independent_components);

  TableIndices<rank_> indices;

  unsigned int remainder = i;
  for (int r = rank_ - 1; r >= 0; --r)
    {
      indices[r] = internal::mod<dim>(remainder);
      remainder  = internal::div<dim>(remainder);
    }
  Assert(remainder == 0, ExcInternalError());

  return indices;
}


template <int rank_, int dim, typename Number>
DEAL_II_CONSTEXPR inline void
Tensor<rank_, dim, Number>::clear()
{
  for (unsigned int i = 0; i < dim; ++i)
    values[i] = internal::NumberType<Number>::value(0.0);
}


template <int rank_, int dim, typename Number>
constexpr std::size_t
Tensor<rank_, dim, Number>::memory_consumption()
{
  return sizeof(Tensor<rank_, dim, Number>);
}


template <int rank_, int dim, typename Number>
template <class Archive>
inline void
Tensor<rank_, dim, Number>::serialize(Archive &ar, const unsigned int)
{
  ar &values;
}


template <int rank_, int dim, typename Number>
constexpr unsigned int Tensor<rank_, dim, Number>::n_independent_components;

#endif // DOXYGEN

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
inline std::ostream &
operator<<(std::ostream &out, const Tensor<rank_, dim, Number> &p)
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
 * @relatesalso Tensor
 */
template <int dim, typename Number>
inline std::ostream &
operator<<(std::ostream &out, const Tensor<0, dim, Number> &p)
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
 * @note This function can also be used in CUDA device code.
 *
 * @relatesalso Tensor
 */
template <int dim, typename Number, typename Other>
DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV inline DEAL_II_ALWAYS_INLINE
  typename ProductType<Other, Number>::type
  operator*(const Other &object, const Tensor<0, dim, Number> &t)
{
  return object * static_cast<const Number &>(t);
}



/**
 * Scalar multiplication of a tensor of rank 0 with an object from the right.
 *
 * This function unwraps the underlying @p Number stored in the Tensor and
 * multiplies @p object with it.
 *
 * @note This function can also be used in CUDA device code.
 *
 * @relatesalso Tensor
 */
template <int dim, typename Number, typename Other>
DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV inline DEAL_II_ALWAYS_INLINE
  typename ProductType<Number, Other>::type
  operator*(const Tensor<0, dim, Number> &t, const Other &object)
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
 * @note This function can also be used in CUDA device code.
 *
 * @relatesalso Tensor
 */
template <int dim, typename Number, typename OtherNumber>
DEAL_II_CUDA_HOST_DEV constexpr DEAL_II_ALWAYS_INLINE
  typename ProductType<Number, OtherNumber>::type
  operator*(const Tensor<0, dim, Number> &     src1,
            const Tensor<0, dim, OtherNumber> &src2)
{
  return static_cast<const Number &>(src1) *
         static_cast<const OtherNumber &>(src2);
}


/**
 * Division of a tensor of rank 0 by a scalar number.
 *
 * @note This function can also be used in CUDA device code.
 *
 * @relatesalso Tensor
 */
template <int dim, typename Number, typename OtherNumber>
DEAL_II_CUDA_HOST_DEV constexpr DEAL_II_ALWAYS_INLINE
  Tensor<0,
         dim,
         typename ProductType<Number,
                              typename EnableIfScalar<OtherNumber>::type>::type>
  operator/(const Tensor<0, dim, Number> &t, const OtherNumber &factor)
{
  return static_cast<const Number &>(t) / factor;
}


/**
 * Add two tensors of rank 0.
 *
 * @note This function can also be used in CUDA device code.
 *
 * @relatesalso Tensor
 */
template <int dim, typename Number, typename OtherNumber>
constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV
                                Tensor<0, dim, typename ProductType<Number, OtherNumber>::type>
                                operator+(const Tensor<0, dim, Number> &     p,
            const Tensor<0, dim, OtherNumber> &q)
{
  return static_cast<const Number &>(p) + static_cast<const OtherNumber &>(q);
}


/**
 * Subtract two tensors of rank 0.
 *
 * @note This function can also be used in CUDA device code.
 *
 * @relatesalso Tensor
 */
template <int dim, typename Number, typename OtherNumber>
constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV
                                Tensor<0, dim, typename ProductType<Number, OtherNumber>::type>
                                operator-(const Tensor<0, dim, Number> &     p,
            const Tensor<0, dim, OtherNumber> &q)
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
 * @note This function can also be used in CUDA device code.
 *
 * @relatesalso Tensor
 */
template <int rank, int dim, typename Number, typename OtherNumber>
DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV inline DEAL_II_ALWAYS_INLINE
                  Tensor<rank,
         dim,
         typename ProductType<Number,
                              typename EnableIfScalar<OtherNumber>::type>::type>
                  operator*(const Tensor<rank, dim, Number> &t, const OtherNumber &factor)
{
  // recurse over the base objects
  Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type> tt;
  for (unsigned int d = 0; d < dim; ++d)
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
 * @note This function can also be used in CUDA device code.
 *
 * @relatesalso Tensor
 */
template <int rank, int dim, typename Number, typename OtherNumber>
DEAL_II_CUDA_HOST_DEV DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
                      Tensor<rank,
         dim,
         typename ProductType<typename EnableIfScalar<Number>::type,
                              OtherNumber>::type>
                      operator*(const Number &factor, const Tensor<rank, dim, OtherNumber> &t)
{
  // simply forward to the operator above
  return t * factor;
}


namespace internal
{
  namespace TensorImplementation
  {
    template <int rank,
              int dim,
              typename Number,
              typename OtherNumber,
              typename std::enable_if<
                !std::is_integral<
                  typename ProductType<Number, OtherNumber>::type>::value,
                int>::type = 0>
    DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV inline DEAL_II_ALWAYS_INLINE
                      Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type>
                      division_operator(const Tensor<rank, dim, Number> &t,
                                        const OtherNumber &              factor)
    {
      Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type> tt;
      const Number inverse_factor = Number(1.) / factor;
      // recurse over the base objects
      for (unsigned int d = 0; d < dim; ++d)
        tt[d] = t[d] * inverse_factor;
      return tt;
    }


    template <int rank,
              int dim,
              typename Number,
              typename OtherNumber,
              typename std::enable_if<
                std::is_integral<
                  typename ProductType<Number, OtherNumber>::type>::value,
                int>::type = 0>
    DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV inline DEAL_II_ALWAYS_INLINE
                      Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type>
                      division_operator(const Tensor<rank, dim, Number> &t,
                                        const OtherNumber &              factor)
    {
      Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type> tt;
      // recurse over the base objects
      for (unsigned int d = 0; d < dim; ++d)
        tt[d] = t[d] / factor;
      return tt;
    }
  } // namespace TensorImplementation
} // namespace internal


/**
 * Division of a tensor of general rank with a scalar number. See the
 * discussion on operator*() above for more information about template
 * arguments and the return type.
 *
 * @note This function can also be used in CUDA device code.
 *
 * @relatesalso Tensor
 */
template <int rank, int dim, typename Number, typename OtherNumber>
DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV inline DEAL_II_ALWAYS_INLINE
                  Tensor<rank,
         dim,
         typename ProductType<Number,
                              typename EnableIfScalar<OtherNumber>::type>::type>
                  operator/(const Tensor<rank, dim, Number> &t, const OtherNumber &factor)
{
  return internal::TensorImplementation::division_operator(t, factor);
}


/**
 * Addition of two tensors of general rank.
 *
 * @tparam rank The rank of both tensors.
 *
 * @note This function can also be used in CUDA device code.
 *
 * @relatesalso Tensor
 */
template <int rank, int dim, typename Number, typename OtherNumber>
DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV inline DEAL_II_ALWAYS_INLINE
                  Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type>
                  operator+(const Tensor<rank, dim, Number> &     p,
            const Tensor<rank, dim, OtherNumber> &q)
{
  Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type> tmp(p);

  for (unsigned int i = 0; i < dim; ++i)
    tmp[i] += q[i];

  return tmp;
}


/**
 * Subtraction of two tensors of general rank.
 *
 * @tparam rank The rank of both tensors.
 *
 * @note This function can also be used in CUDA device code.
 *
 * @relatesalso Tensor
 */
template <int rank, int dim, typename Number, typename OtherNumber>
DEAL_II_CONSTEXPR DEAL_II_CUDA_HOST_DEV inline DEAL_II_ALWAYS_INLINE
                  Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type>
                  operator-(const Tensor<rank, dim, Number> &     p,
            const Tensor<rank, dim, OtherNumber> &q)
{
  Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type> tmp(p);

  for (unsigned int i = 0; i < dim; ++i)
    tmp[i] -= q[i];

  return tmp;
}

/**
 * Entrywise multiplication of two tensor objects of rank 0 (i.e. the
 * multiplication of two scalar values).
 *
 * @relatesalso Tensor
 */
template <int dim, typename Number, typename OtherNumber>
inline DEAL_II_CONSTEXPR DEAL_II_ALWAYS_INLINE
                         Tensor<0, dim, typename ProductType<Number, OtherNumber>::type>
                         schur_product(const Tensor<0, dim, Number> &     src1,
                                       const Tensor<0, dim, OtherNumber> &src2)
{
  Tensor<0, dim, typename ProductType<Number, OtherNumber>::type> tmp(src1);

  tmp *= src2;

  return tmp;
}

/**
 * Entrywise multiplication of two tensor objects of general rank.
 *
 * This multiplication is also called "Hadamard-product" (c.f.
 * https://en.wikipedia.org/wiki/Hadamard_product_(matrices)), and generates a
 * new tensor of size <rank, dim>:
 * @f[
 *   \text{result}_{i, j}
 *   = \text{left}_{i, j}\circ
 *     \text{right}_{i, j}
 * @f]
 *
 * @tparam rank The rank of both tensors.
 *
 * @relatesalso Tensor
 */
template <int rank, int dim, typename Number, typename OtherNumber>
inline DEAL_II_CONSTEXPR DEAL_II_ALWAYS_INLINE
                         Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type>
                         schur_product(const Tensor<rank, dim, Number> &     src1,
                                       const Tensor<rank, dim, OtherNumber> &src2)
{
  Tensor<rank, dim, typename ProductType<Number, OtherNumber>::type> tmp;

  for (unsigned int i = 0; i < dim; ++i)
    tmp[i] = schur_product(Tensor<rank - 1, dim, Number>(src1[i]),
                           Tensor<rank - 1, dim, OtherNumber>(src2[i]));

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
template <int rank_1,
          int rank_2,
          int dim,
          typename Number,
          typename OtherNumber,
          typename = typename std::enable_if<rank_1 >= 1 && rank_2 >= 1>::type>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
  typename Tensor<rank_1 + rank_2 - 2,
                  dim,
                  typename ProductType<Number, OtherNumber>::type>::tensor_type
  operator*(const Tensor<rank_1, dim, Number> &     src1,
            const Tensor<rank_2, dim, OtherNumber> &src2)
{
  typename Tensor<rank_1 + rank_2 - 2,
                  dim,
                  typename ProductType<Number, OtherNumber>::type>::tensor_type
    result{};

  TensorAccessors::internal::
    ReorderedIndexView<0, rank_2, const Tensor<rank_2, dim, OtherNumber>>
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
 * (<code>index_2==2</code>) of a tensor <code>t2</code>, this function should
 * be invoked as
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
template <int index_1,
          int index_2,
          int rank_1,
          int rank_2,
          int dim,
          typename Number,
          typename OtherNumber>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
  typename Tensor<rank_1 + rank_2 - 2,
                  dim,
                  typename ProductType<Number, OtherNumber>::type>::tensor_type
  contract(const Tensor<rank_1, dim, Number> &     src1,
           const Tensor<rank_2, dim, OtherNumber> &src2)
{
  Assert(0 <= index_1 && index_1 < rank_1,
         ExcMessage(
           "The specified index_1 must lie within the range [0,rank_1)"));
  Assert(0 <= index_2 && index_2 < rank_2,
         ExcMessage(
           "The specified index_2 must lie within the range [0,rank_2)"));

  using namespace TensorAccessors;
  using namespace TensorAccessors::internal;

  // Reorder index_1 to the end of src1:
  ReorderedIndexView<index_1, rank_1, const Tensor<rank_1, dim, Number>>
    reord_01 = reordered_index_view<index_1, rank_1>(src1);

  // Reorder index_2 to the end of src2:
  ReorderedIndexView<index_2, rank_2, const Tensor<rank_2, dim, OtherNumber>>
    reord_02 = reordered_index_view<index_2, rank_2>(src2);

  typename Tensor<rank_1 + rank_2 - 2,
                  dim,
                  typename ProductType<Number, OtherNumber>::type>::tensor_type
    result{};
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
 * (<code>index_4==0</code>) of a tensor <code>t2</code>, this function should
 * be invoked as
 * @code
 *   double_contract<0, 2, 1, 0>(t1, t2);
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
template <int index_1,
          int index_2,
          int index_3,
          int index_4,
          int rank_1,
          int rank_2,
          int dim,
          typename Number,
          typename OtherNumber>
DEAL_II_CONSTEXPR inline
  typename Tensor<rank_1 + rank_2 - 4,
                  dim,
                  typename ProductType<Number, OtherNumber>::type>::tensor_type
  double_contract(const Tensor<rank_1, dim, Number> &     src1,
                  const Tensor<rank_2, dim, OtherNumber> &src2)
{
  Assert(0 <= index_1 && index_1 < rank_1,
         ExcMessage(
           "The specified index_1 must lie within the range [0,rank_1)"));
  Assert(0 <= index_3 && index_3 < rank_1,
         ExcMessage(
           "The specified index_3 must lie within the range [0,rank_1)"));
  Assert(index_1 != index_3,
         ExcMessage("index_1 and index_3 must not be the same"));
  Assert(0 <= index_2 && index_2 < rank_2,
         ExcMessage(
           "The specified index_2 must lie within the range [0,rank_2)"));
  Assert(0 <= index_4 && index_4 < rank_2,
         ExcMessage(
           "The specified index_4 must lie within the range [0,rank_2)"));
  Assert(index_2 != index_4,
         ExcMessage("index_2 and index_4 must not be the same"));

  using namespace TensorAccessors;
  using namespace TensorAccessors::internal;

  // Reorder index_1 to the end of src1:
  ReorderedIndexView<index_1, rank_1, const Tensor<rank_1, dim, Number>>
    reord_1 = TensorAccessors::reordered_index_view<index_1, rank_1>(src1);

  // Reorder index_2 to the end of src2:
  ReorderedIndexView<index_2, rank_2, const Tensor<rank_2, dim, OtherNumber>>
    reord_2 = TensorAccessors::reordered_index_view<index_2, rank_2>(src2);

  // Now, reorder index_3 to the end of src1. We have to make sure to
  // preserve the original ordering: index_1 has been removed. If
  // index_3 > index_1, we have to use (index_3 - 1) instead:
  ReorderedIndexView<
    (index_3 < index_1 ? index_3 : index_3 - 1),
    rank_1,
    ReorderedIndexView<index_1, rank_1, const Tensor<rank_1, dim, Number>>>
    reord_3 =
      TensorAccessors::reordered_index_view < index_3 < index_1 ? index_3 :
                                                                  index_3 - 1,
    rank_1 > (reord_1);

  // Now, reorder index_4 to the end of src2. We have to make sure to
  // preserve the original ordering: index_2 has been removed. If
  // index_4 > index_2, we have to use (index_4 - 1) instead:
  ReorderedIndexView<
    (index_4 < index_2 ? index_4 : index_4 - 1),
    rank_2,
    ReorderedIndexView<index_2, rank_2, const Tensor<rank_2, dim, OtherNumber>>>
    reord_4 =
      TensorAccessors::reordered_index_view < index_4 < index_2 ? index_4 :
                                                                  index_4 - 1,
    rank_2 > (reord_2);

  typename Tensor<rank_1 + rank_2 - 4,
                  dim,
                  typename ProductType<Number, OtherNumber>::type>::tensor_type
    result{};
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
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
  typename ProductType<Number, OtherNumber>::type
  scalar_product(const Tensor<rank, dim, Number> &     left,
                 const Tensor<rank, dim, OtherNumber> &right)
{
  typename ProductType<Number, OtherNumber>::type result{};
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
          int rank_1,
          int rank_2,
          int dim,
          typename T1,
          typename T2,
          typename T3>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
  typename ProductType<T1, typename ProductType<T2, T3>::type>::type
  contract3(const TensorT1<rank_1, dim, T1> &         left,
            const TensorT2<rank_1 + rank_2, dim, T2> &middle,
            const TensorT3<rank_2, dim, T3> &         right)
{
  using return_type =
    typename ProductType<T1, typename ProductType<T2, T3>::type>::type;
  return TensorAccessors::contract3<rank_1, rank_2, dim, return_type>(left,
                                                                      middle,
                                                                      right);
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
template <int rank_1,
          int rank_2,
          int dim,
          typename Number,
          typename OtherNumber>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
  Tensor<rank_1 + rank_2, dim, typename ProductType<Number, OtherNumber>::type>
  outer_product(const Tensor<rank_1, dim, Number> &     src1,
                const Tensor<rank_2, dim, OtherNumber> &src2)
{
  typename Tensor<rank_1 + rank_2,
                  dim,
                  typename ProductType<Number, OtherNumber>::type>::tensor_type
    result{};
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
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE Tensor<1, dim, Number>
                                               cross_product_2d(const Tensor<1, dim, Number> &src)
{
  Assert(dim == 2, ExcInternalError());

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
template <int dim, typename Number1, typename Number2>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE
  Tensor<1, dim, typename ProductType<Number1, Number2>::type>
  cross_product_3d(const Tensor<1, dim, Number1> &src1,
                   const Tensor<1, dim, Number2> &src2)
{
  Assert(dim == 3, ExcInternalError());

  Tensor<1, dim, typename ProductType<Number1, Number2>::type> result;

  // avoid compiler warnings
  constexpr int s0 = 0 % dim;
  constexpr int s1 = 1 % dim;
  constexpr int s2 = 2 % dim;

  result[s0] = src1[s1] * src2[s2] - src1[s2] * src2[s1];
  result[s1] = src1[s2] * src2[s0] - src1[s0] * src2[s2];
  result[s2] = src1[s0] * src2[s1] - src1[s1] * src2[s0];

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
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE Number
                                               determinant(const Tensor<2, dim, Number> &t)
{
  // Compute the determinant using the Laplace expansion of the
  // determinant. We expand along the last row.
  Number det = internal::NumberType<Number>::value(0.0);

  for (unsigned int k = 0; k < dim; ++k)
    {
      Tensor<2, dim - 1, Number> minor;
      for (unsigned int i = 0; i < dim - 1; ++i)
        for (unsigned int j = 0; j < dim - 1; ++j)
          minor[i][j] = t[i][j < k ? j : j + 1];

      const Number cofactor = ((k % 2 == 0) ? -1. : 1.) * determinant(minor);

      det += t[dim - 1][k] * cofactor;
    }

  return ((dim % 2 == 0) ? 1. : -1.) * det;
}

/**
 * Specialization for dim==1.
 *
 * @relatesalso Tensor
 */
template <typename Number>
constexpr DEAL_II_ALWAYS_INLINE Number
                                determinant(const Tensor<2, 1, Number> &t)
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
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE Number
                                               trace(const Tensor<2, dim, Number> &d)
{
  Number t = d[0][0];
  for (unsigned int i = 1; i < dim; ++i)
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
DEAL_II_CONSTEXPR inline Tensor<2, dim, Number>
invert(const Tensor<2, dim, Number> &)
{
  Number return_tensor[dim][dim];

  // if desired, take over the
  // inversion of a 4x4 tensor
  // from the FullMatrix
  AssertThrow(false, ExcNotImplemented());

  return Tensor<2, dim, Number>(return_tensor);
}


#ifndef DOXYGEN

template <typename Number>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE Tensor<2, 1, Number>
                                               invert(const Tensor<2, 1, Number> &t)
{
  Tensor<2, 1, Number> return_tensor;

  return_tensor[0][0] = internal::NumberType<Number>::value(1.0 / t[0][0]);

  return return_tensor;
}


template <typename Number>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE Tensor<2, 2, Number>
                                               invert(const Tensor<2, 2, Number> &t)
{
  Tensor<2, 2, Number> return_tensor;

  // this is Maple output,
  // thus a bit unstructured
  const Number inv_det_t = internal::NumberType<Number>::value(
    1.0 / (t[0][0] * t[1][1] - t[1][0] * t[0][1]));
  return_tensor[0][0] = t[1][1];
  return_tensor[0][1] = -t[0][1];
  return_tensor[1][0] = -t[1][0];
  return_tensor[1][1] = t[0][0];
  return_tensor *= inv_det_t;

  return return_tensor;
}


template <typename Number>
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE Tensor<2, 3, Number>
                                               invert(const Tensor<2, 3, Number> &t)
{
  Tensor<2, 3, Number> return_tensor;

  const Number t4  = internal::NumberType<Number>::value(t[0][0] * t[1][1]),
               t6  = internal::NumberType<Number>::value(t[0][0] * t[1][2]),
               t8  = internal::NumberType<Number>::value(t[0][1] * t[1][0]),
               t00 = internal::NumberType<Number>::value(t[0][2] * t[1][0]),
               t01 = internal::NumberType<Number>::value(t[0][1] * t[2][0]),
               t04 = internal::NumberType<Number>::value(t[0][2] * t[2][0]),
               inv_det_t = internal::NumberType<Number>::value(
                 1.0 / (t4 * t[2][2] - t6 * t[2][1] - t8 * t[2][2] +
                        t00 * t[2][1] + t01 * t[1][2] - t04 * t[1][1]));
  return_tensor[0][0] = internal::NumberType<Number>::value(t[1][1] * t[2][2]) -
                        internal::NumberType<Number>::value(t[1][2] * t[2][1]);
  return_tensor[0][1] = internal::NumberType<Number>::value(t[0][2] * t[2][1]) -
                        internal::NumberType<Number>::value(t[0][1] * t[2][2]);
  return_tensor[0][2] = internal::NumberType<Number>::value(t[0][1] * t[1][2]) -
                        internal::NumberType<Number>::value(t[0][2] * t[1][1]);
  return_tensor[1][0] = internal::NumberType<Number>::value(t[1][2] * t[2][0]) -
                        internal::NumberType<Number>::value(t[1][0] * t[2][2]);
  return_tensor[1][1] =
    internal::NumberType<Number>::value(t[0][0] * t[2][2]) - t04;
  return_tensor[1][2] = t00 - t6;
  return_tensor[2][0] = internal::NumberType<Number>::value(t[1][0] * t[2][1]) -
                        internal::NumberType<Number>::value(t[1][1] * t[2][0]);
  return_tensor[2][1] =
    t01 - internal::NumberType<Number>::value(t[0][0] * t[2][1]);
  return_tensor[2][2] = internal::NumberType<Number>::value(t4 - t8);
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
DEAL_II_CONSTEXPR inline DEAL_II_ALWAYS_INLINE Tensor<2, dim, Number>
                                               transpose(const Tensor<2, dim, Number> &t)
{
  Tensor<2, dim, Number> tt;
  for (unsigned int i = 0; i < dim; ++i)
    {
      tt[i][i] = t[i][i];
      for (unsigned int j = i + 1; j < dim; ++j)
        {
          tt[i][j] = t[j][i];
          tt[j][i] = t[i][j];
        };
    }
  return tt;
}


/**
 * Return the adjugate of the given tensor of rank 2.
 * The adjugate of a tensor $\mathbf A$ is defined as
 * @f[
 *  \textrm{adj}\mathbf A
 *   \dealcoloneq \textrm{det}\mathbf A \; \mathbf{A}^{-1}
 * \; .
 * @f]
 *
 * @note This requires that the tensor is invertible.
 *
 * @relatesalso Tensor
 * @author Jean-Paul Pelteret, 2016
 */
template <int dim, typename Number>
constexpr Tensor<2, dim, Number>
adjugate(const Tensor<2, dim, Number> &t)
{
  return determinant(t) * invert(t);
}


/**
 * Return the cofactor of the given tensor of rank 2.
 * The cofactor of a tensor $\mathbf A$ is defined as
 * @f[
 *  \textrm{cof}\mathbf A
 *   \dealcoloneq \textrm{det}\mathbf A \; \mathbf{A}^{-T}
 *    = \left[ \textrm{adj}\mathbf A \right]^{T} \; .
 * @f]
 *
 * @note This requires that the tensor is invertible.
 *
 * @relatesalso Tensor
 * @author Jean-Paul Pelteret, 2016
 */
template <int dim, typename Number>
constexpr Tensor<2, dim, Number>
cofactor(const Tensor<2, dim, Number> &t)
{
  return transpose(adjugate(t));
}


/**
 * Return the nearest orthogonal matrix by
 * combining the products of the SVD decomposition: $\mathbf U \mathbf{V}^T$,
 * where $\mathbf U$ and $\mathbf V$ are computed from the SVD decomposition:
 * $\mathbf U  \mathbf S \mathbf V^T$,
 * effectively replacing $\mathbf S$ with the identity matrix.
 * @param tensor The tensor which to find the closest orthogonal
 * tensor to.
 * @relatesalso Tensor
 */
template <int dim, typename Number>
Tensor<2, dim, Number>
project_onto_orthogonal_tensors(const Tensor<2, dim, Number> &tensor)
{
  Tensor<2, dim, Number>   output_tensor;
  FullMatrix<Number>       matrix(dim);
  LAPACKFullMatrix<Number> lapack_matrix(dim);
  LAPACKFullMatrix<Number> result(dim);

  // todo: find or add dealii functionality to copy in one step.
  matrix.copy_from(tensor);
  lapack_matrix.copy_from(matrix);

  // now compute the svd of the matrices
  lapack_matrix.compute_svd();

  // Use the SVD results to orthogonalize: $U V^T$
  lapack_matrix.get_svd_u().mmult(result, lapack_matrix.get_svd_vt());

  // todo: find or add dealii functionality to copy in one step.
  matrix = result;
  matrix.copy_to(output_tensor);
  return output_tensor;
}


/**
 * Return the $l_1$ norm of the given rank-2 tensor, where
 * $\|\mathbf T\|_1 = \max_j \sum_i |T_{ij}|$
 * (maximum of the sums over columns).
 *
 * @relatesalso Tensor
 * @author Wolfgang Bangerth, 2012
 */
template <int dim, typename Number>
inline Number
l1_norm(const Tensor<2, dim, Number> &t)
{
  Number max = internal::NumberType<Number>::value(0.0);
  for (unsigned int j = 0; j < dim; ++j)
    {
      Number sum = internal::NumberType<Number>::value(0.0);
      for (unsigned int i = 0; i < dim; ++i)
        sum += std::fabs(t[i][j]);

      if (sum > max)
        max = sum;
    }

  return max;
}


/**
 * Return the $l_\infty$ norm of the given rank-2 tensor, where
 * $\|\mathbf T\|_\infty = \max_i \sum_j |T_{ij}|$
 * (maximum of the sums over rows).
 *
 * @relatesalso Tensor
 * @author Wolfgang Bangerth, 2012
 */
template <int dim, typename Number>
inline Number
linfty_norm(const Tensor<2, dim, Number> &t)
{
  Number max = internal::NumberType<Number>::value(0.0);
  for (unsigned int i = 0; i < dim; ++i)
    {
      Number sum = internal::NumberType<Number>::value(0.0);
      for (unsigned int j = 0; j < dim; ++j)
        sum += std::fabs(t[i][j]);

      if (sum > max)
        max = sum;
    }

  return max;
}

//@}


#ifndef DOXYGEN


#  ifdef DEAL_II_ADOLC_WITH_ADVANCED_BRANCHING

// Specialization of functions for ADOL-C number types when
// the advanced branching feature is used
template <int dim>
inline adouble
l1_norm(const Tensor<2, dim, adouble> &t)
{
  adouble max = internal::NumberType<adouble>::value(0.0);
  for (unsigned int j = 0; j < dim; ++j)
    {
      adouble sum = internal::NumberType<adouble>::value(0.0);
      for (unsigned int i = 0; i < dim; ++i)
        sum += std::fabs(t[i][j]);

      condassign(max, (sum > max), sum, max);
    }

  return max;
}


template <int dim>
inline adouble
linfty_norm(const Tensor<2, dim, adouble> &t)
{
  adouble max = internal::NumberType<adouble>::value(0.0);
  for (unsigned int i = 0; i < dim; ++i)
    {
      adouble sum = internal::NumberType<adouble>::value(0.0);
      for (unsigned int j = 0; j < dim; ++j)
        sum += std::fabs(t[i][j]);

      condassign(max, (sum > max), sum, max);
    }

  return max;
}

#  endif // DEAL_II_ADOLC_WITH_ADVANCED_BRANCHING


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
