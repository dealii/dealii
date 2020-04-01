// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2019 by the deal.II authors
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

#ifndef dealii_point_h
#define dealii_point_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/tensor.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/geometries/point.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#include <cmath>

DEAL_II_NAMESPACE_OPEN

/**
 * A class that represents a point in a Cartesian space of dimension
 * @p dim .
 *
 * Objects of this class are used to represent points (i.e., vectors
 * anchored at the origin) of a vector space equipped with a <a
 * href="https://en.wikipedia.org/wiki/Cartesian_coordinate_system">Cartesian
 * coordinate system</a>. They are, among other uses, passed to
 * functions that operate on points in spaces of a priori fixed
 * dimension: rather than using functions like <code>double f(const
 * double x)</code> and <code>double f(const double x, const double
 * y)</code>, you can use <code>double f(const Point<dim> &p)</code>
 * instead as it allows writing dimension independent code.
 *
 * deal.II specifically uses Point objects as indicating points that
 * are represented by Cartesian coordinates, i.e., where a point in @p
 * dim space dimensions is characterized by signed distances along the
 * axes of a coordinate system spanned by @p dim mutually orthogonal
 * unit vectors (called the "coordinate axes"). This choice of
 * representing a vector makes addition and scaling of vectors
 * particularly simple: one only has to add or multiply each
 * coordinate value. On the other hand, adding or scaling vectors is
 * not nearly as simple when a vector is represented in other kinds of
 * coordinate systems (e.g., <a
 * href="https://en.wikipedia.org/wiki/Spherical_coordinate_system">spherical
 * coordinate systems</a>).
 *
 *
 * <h3>What's a <code>Point@<dim@></code> and what is a
 * <code>Tensor@<1,dim@></code>?</h3>
 *
 * The Point class is derived from Tensor@<1,dim@> and consequently shares the
 * latter's member functions and other attributes. In fact, it has relatively
 * few additional functions itself (the most notable exception being the
 * distance() function to compute the Euclidean distance between two points in
 * space), and these two classes can therefore often be used interchangeably.
 *
 * Nonetheless, there are semantic differences that make us use these classes
 * in different and well-defined contexts. Within deal.II, we use the
 * <tt>Point</tt> class to denote points in space, i.e., for vectors (rank-1
 * tensors) that are <em>anchored at the origin</em>. On the other hand,
 * vectors that are anchored elsewhere (and consequently do not represent
 * <em>points</em> in the common usage of the word) are represented by objects
 * of type Tensor@<1,dim@>. In particular, this is the case for direction
 * vectors, normal vectors, gradients, and the differences between two points
 * (i.e., what you get when you subtract one point from another): all of these
 * are represented by Tensor@<1,dim@> objects rather than Point@<dim@>.
 *
 * Furthermore, the Point class is only used where the coordinates of an
 * object can be thought to possess the dimension of a length. An object that
 * represents the weight, height, and cost of an object is neither a point nor
 * a tensor (because it lacks the transformation properties under rotation of
 * the coordinate system) and should consequently not be represented by either
 * of these classes. Use an array of size 3 in this case, or the
 * <code>std::array</code> class. Alternatively, as in the case of
 * vector-valued functions, you can use objects of type Vector or
 * <code>std::vector</code>.
 *
 *
 * @tparam dim An integer that denotes the dimension of the space in which a
 * point lies. This of course equals the number of coordinates that identify a
 * point.
 * @tparam Number The data type in which the coordinates values are to be
 * stored. This will, in almost all cases, simply be the default @p double,
 * but there are cases where one may want to store coordinates in a different
 * (and always scalar) type. An example would be an interval type that can
 * store the value of a coordinate as well as its uncertainty. Another example
 * would be a type that allows for Automatic Differentiation (see, for
 * example, the Sacado type used in step-33) and thereby can generate analytic
 * (spatial) derivatives of a function when passed a Point object whose
 * coordinates are stored in such a type.
 *
 *
 * @ingroup geomprimitives
 * @author Wolfgang Bangerth, 1997
 */
template <int dim, typename Number = double>
class Point : public Tensor<1, dim, Number>
{
public:
  /**
   * Standard constructor. Creates an object that corresponds to the origin,
   * i.e., all coordinates are set to zero.
   *
   * @note This function can also be used in CUDA device code.
   */
  DEAL_II_CUDA_HOST_DEV
  Point();

  /**
   * Convert a tensor to a point.
   */
  explicit DEAL_II_CUDA_HOST_DEV
  Point(const Tensor<1, dim, Number> &);

  /**
   * Constructor for one dimensional points. This function is only implemented
   * for <tt>dim==1</tt> since the usage is considered unsafe for points with
   * <tt>dim!=1</tt> as it would leave some components of the point
   * coordinates uninitialized.
   *
   * @note This function can also be used in CUDA device code.
   */
  explicit DEAL_II_CUDA_HOST_DEV
  Point(const Number x);

  /**
   * Constructor for two dimensional points. This function is only implemented
   * for <tt>dim==2</tt> since the usage is considered unsafe for points with
   * <tt>dim!=2</tt> as it would leave some components of the point
   * coordinates uninitialized (if dim>2) or would not use some arguments (if
   * dim<2).
   *
   * @note This function can also be used in CUDA device code.
   */
  DEAL_II_CUDA_HOST_DEV
  Point(const Number x, const Number y);

  /**
   * Constructor for three dimensional points. This function is only
   * implemented for <tt>dim==3</tt> since the usage is considered unsafe for
   * points with <tt>dim!=3</tt> as it would leave some components of the
   * point coordinates uninitialized (if dim>3) or would not use some
   * arguments (if dim<3).
   *
   * @note This function can also be used in CUDA device code.
   */
  DEAL_II_CUDA_HOST_DEV
  Point(const Number x, const Number y, const Number z);

  /**
   * Convert a boost::geometry::point to a dealii::Point.
   */
  template <std::size_t dummy_dim,
            typename std::enable_if<(dim == dummy_dim) && (dummy_dim != 0),
                                    int>::type = 0>
  Point(const boost::geometry::model::
          point<Number, dummy_dim, boost::geometry::cs::cartesian> &boost_pt);

  /**
   * Return a unit vector in coordinate direction <tt>i</tt>, i.e., a vector
   * that is zero in all coordinates except for a single 1 in the <tt>i</tt>th
   * coordinate.
   *
   * @note This function can also be used in CUDA device code.
   */
  static DEAL_II_CUDA_HOST_DEV Point<dim, Number>
                               unit_vector(const unsigned int i);

  /**
   * Read access to the <tt>index</tt>th coordinate.
   *
   * @note This function can also be used in CUDA device code.
   */
  DEAL_II_CUDA_HOST_DEV Number
                        operator()(const unsigned int index) const;

  /**
   * Read and write access to the <tt>index</tt>th coordinate.
   *
   * @note This function can also be used in CUDA device code.
   */
  DEAL_II_CUDA_HOST_DEV Number &
                        operator()(const unsigned int index);

  /**
   * Assignment operator from Tensor<1, dim, Number> with different underlying
   * scalar type. This obviously requires that the @p OtherNumber type is
   * convertible to @p Number.
   */
  template <typename OtherNumber>
  Point<dim, Number> &
  operator=(const Tensor<1, dim, OtherNumber> &p);

  /**
   * @name Addition and subtraction of points.
   * @{
   */

  /**
   * Add an offset given as Tensor<1,dim,Number> to a point.
   *
   * @note This function can also be used in CUDA device code.
   */
  DEAL_II_CUDA_HOST_DEV Point<dim, Number>
                        operator+(const Tensor<1, dim, Number> &) const;

  /**
   * Subtract two points, i.e., obtain the vector that connects the two. As
   * discussed in the documentation of this class, subtracting two points
   * results in a vector anchored at one of the two points (rather than at the
   * origin) and, consequently, the result is returned as a Tensor@<1,dim@>
   * rather than as a Point@<dim@>.
   *
   * @note This function can also be used in CUDA device code.
   */
  DEAL_II_CUDA_HOST_DEV Tensor<1, dim, Number>
                        operator-(const Point<dim, Number> &) const;

  /**
   * Subtract a difference vector (represented by a Tensor@<1,dim@>) from the
   * current point. This results in another point and, as discussed in the
   * documentation of this class, the result is then naturally returned as a
   * Point@<dim@> object rather than as a Tensor@<1,dim@>.
   *
   * @note This function can also be used in CUDA device code.
   */
  DEAL_II_CUDA_HOST_DEV Point<dim, Number>
                        operator-(const Tensor<1, dim, Number> &) const;

  /**
   * The opposite vector.
   *
   * @note This function can also be used in CUDA device code.
   */
  DEAL_II_CUDA_HOST_DEV Point<dim, Number>
                        operator-() const;

  /**
   * @}
   */

  /**
   * @name Multiplication and scaling of points. Dot products. Norms.
   * @{
   */

  /**
   * Multiply the current point by a factor.
   *
   * @note This function can also be used in CUDA device code.
   *
   * @relatesalso EnableIfScalar
   */
  template <typename OtherNumber>
  DEAL_II_CUDA_HOST_DEV Point<
    dim,
    typename ProductType<Number,
                         typename EnableIfScalar<OtherNumber>::type>::type>
  operator*(const OtherNumber) const;

  /**
   * Divide the current point by a factor.
   *
   * @note This function can also be used in CUDA device code.
   */
  template <typename OtherNumber>
  DEAL_II_CUDA_HOST_DEV Point<
    dim,
    typename ProductType<Number,
                         typename EnableIfScalar<OtherNumber>::type>::type>
  operator/(const OtherNumber) const;

  /**
   * Return the scalar product of the vectors representing two points.
   *
   * @note This function can also be used in CUDA device code.
   */
  DEAL_II_CUDA_HOST_DEV Number operator*(const Tensor<1, dim, Number> &p) const;

  /**
   * Return the scalar product of this point vector with itself, i.e. the
   * square, or the square of the norm. In case of a complex number type it is
   * equivalent to the contraction of this point vector with a complex
   * conjugate of itself.
   *
   * @note This function is equivalent to
   * Tensor<rank,dim,Number>::norm_square() which returns the square of the
   * Frobenius norm.
   *
   * @note This function can also be used in CUDA device code.
   */
  DEAL_II_CUDA_HOST_DEV typename numbers::NumberTraits<Number>::real_type
  square() const;

  /**
   * Return the Euclidean distance of <tt>this</tt> point to the point
   * <tt>p</tt>, i.e. the $l_2$ norm of the difference between the
   * vectors representing the two points.
   *
   * @note This function can also be used in CUDA device code.
   */
  DEAL_II_CUDA_HOST_DEV typename numbers::NumberTraits<Number>::real_type
  distance(const Point<dim, Number> &p) const;

  /**
   * Return the squared Euclidean distance of <tt>this</tt> point to the point
   * <tt>p</tt>.
   *
   * @note This function can also be used in CUDA device code.
   */
  DEAL_II_CUDA_HOST_DEV typename numbers::NumberTraits<Number>::real_type
  distance_square(const Point<dim, Number> &p) const;

  /**
   * @}
   */

  /**
   * Read or write the data of this object to or from a stream for the purpose
   * of serialization
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version);
};

/*--------------------------- Inline functions: Point -----------------------*/

#ifndef DOXYGEN

// At least clang-3.7 requires us to have a user-defined constructor
// and we can't use 'Point<dim,Number>::Point () = default' here.
template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV
Point<dim, Number>::Point() // NOLINT
{}



template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV
Point<dim, Number>::Point(const Tensor<1, dim, Number> &t)
  : Tensor<1, dim, Number>(t)
{}



template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV
Point<dim, Number>::Point(const Number x)
{
#  ifndef __CUDA_ARCH__
  Assert(dim == 1,
         ExcMessage(
           "You can only initialize Point<1> objects using the constructor "
           "that takes only one argument. Point<dim> objects with dim!=1 "
           "require initialization with the constructor that takes 'dim' "
           "arguments."));
#  endif

  // we can only get here if we pass the assertion. use the switch anyway so
  // as to avoid compiler warnings about uninitialized elements or writing
  // beyond the end of the 'values' array
  switch (dim)
    {
      case 1:
        this->values[0] = x;
        break;

      default:;
    }
}



template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV
Point<dim, Number>::Point(const Number x, const Number y)
{
#  ifndef __CUDA_ARCH__
  Assert(dim == 2,
         ExcMessage(
           "You can only initialize Point<2> objects using the constructor "
           "that takes two arguments. Point<dim> objects with dim!=2 "
           "require initialization with the constructor that takes 'dim' "
           "arguments."));
#  endif

  // we can only get here if we pass the assertion. use the indirection anyway
  // so as to avoid compiler warnings about uninitialized elements or writing
  // beyond the end of the 'values' array
  constexpr unsigned int y_index = (dim < 2) ? 0 : 1;
  this->values[0]                = x;
  this->values[y_index]          = y;
}



template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV
Point<dim, Number>::Point(const Number x, const Number y, const Number z)
{
#  ifndef __CUDA_ARCH__
  Assert(dim == 3,
         ExcMessage(
           "You can only initialize Point<3> objects using the constructor "
           "that takes three arguments. Point<dim> objects with dim!=3 "
           "require initialization with the constructor that takes 'dim' "
           "arguments."));
#  endif

  // we can only get here if we pass the assertion. use the indirection anyway
  // so as to avoid compiler warnings about uninitialized elements or writing
  // beyond the end of the 'values' array
  constexpr unsigned int y_index = (dim < 2) ? 0 : 1;
  constexpr unsigned int z_index = (dim < 3) ? 0 : 2;
  this->values[0]                = x;
  this->values[y_index]          = y;
  this->values[z_index]          = z;
}



template <int dim, typename Number>
template <
  std::size_t dummy_dim,
  typename std::enable_if<(dim == dummy_dim) && (dummy_dim != 0), int>::type>
inline Point<dim, Number>::Point(
  const boost::geometry::model::
    point<Number, dummy_dim, boost::geometry::cs::cartesian> &boost_pt)
{
  Assert(dim <= 3, ExcNotImplemented());
  this->values[0]                = boost::geometry::get<0>(boost_pt);
  constexpr unsigned int y_index = (dim < 2) ? 0 : 1;
  constexpr unsigned int z_index = (dim < 3) ? 0 : 2;

  if (dim >= 2)
    this->values[y_index] = boost::geometry::get<y_index>(boost_pt);

  if (dim >= 3)
    this->values[z_index] = boost::geometry::get<z_index>(boost_pt);
}



template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV Point<dim, Number>
                             Point<dim, Number>::unit_vector(unsigned int i)
{
  Point<dim, Number> p;
  p[i] = 1.;
  return p;
}


template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV Number
Point<dim, Number>::operator()(const unsigned int index) const
{
#  ifndef __CUDA_ARCH__
  AssertIndexRange(index, dim);
#  endif
  return this->values[index];
}



template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV Number &
Point<dim, Number>::operator()(const unsigned int index)
{
#  ifndef __CUDA_ARCH__
  AssertIndexRange(index, dim);
#  endif
  return this->values[index];
}



template <int dim, typename Number>
template <typename OtherNumber>
inline DEAL_II_ALWAYS_INLINE Point<dim, Number> &
Point<dim, Number>::operator=(const Tensor<1, dim, OtherNumber> &p)
{
  Tensor<1, dim, Number>::operator=(p);
  return *this;
}



template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV Point<dim, Number>
Point<dim, Number>::operator+(const Tensor<1, dim, Number> &p) const
{
  Point<dim, Number> tmp = *this;
  tmp += p;
  return tmp;
}



template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV Tensor<1, dim, Number>
Point<dim, Number>::operator-(const Point<dim, Number> &p) const
{
  return (Tensor<1, dim, Number>(*this) -= p);
}



template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV Point<dim, Number>
Point<dim, Number>::operator-(const Tensor<1, dim, Number> &p) const
{
  Point<dim, Number> tmp = *this;
  tmp -= p;
  return tmp;
}



template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV Point<dim, Number>
Point<dim, Number>::operator-() const
{
  Point<dim, Number> result;
  for (unsigned int i = 0; i < dim; ++i)
    result.values[i] = -this->values[i];
  return result;
}



template <int dim, typename Number>
template <typename OtherNumber>
inline DEAL_II_CUDA_HOST_DEV
    Point<dim,
        typename ProductType<Number,
                             typename EnableIfScalar<OtherNumber>::type>::type>
    Point<dim, Number>::operator*(const OtherNumber factor) const
{
  Point<dim, typename ProductType<Number, OtherNumber>::type> tmp;
  for (unsigned int i = 0; i < dim; ++i)
    tmp[i] = this->operator[](i) * factor;
  return tmp;
}



template <int dim, typename Number>
template <typename OtherNumber>
inline DEAL_II_CUDA_HOST_DEV
  Point<dim,
        typename ProductType<Number,
                             typename EnableIfScalar<OtherNumber>::type>::type>
  Point<dim, Number>::operator/(const OtherNumber factor) const
{
  const Tensor<1, dim, Number> &base_object = *this;
  return Point<
    dim,
    typename ProductType<Number,
                         typename EnableIfScalar<OtherNumber>::type>::type>(
    dealii::operator/(base_object, factor));
}



template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV Number Point<dim, Number>::
                                    operator*(const Tensor<1, dim, Number> &p) const
{
  Number res = Number();
  for (unsigned int i = 0; i < dim; ++i)
    res += this->operator[](i) * p[i];
  return res;
}


template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV typename numbers::NumberTraits<Number>::real_type
Point<dim, Number>::square() const
{
  return this->norm_square();
}



template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV typename numbers::NumberTraits<Number>::real_type
Point<dim, Number>::distance(const Point<dim, Number> &p) const
{
  return std::sqrt(distance_square(p));
}



template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV typename numbers::NumberTraits<Number>::real_type
Point<dim, Number>::distance_square(const Point<dim, Number> &p) const
{
  Number sum = internal::NumberType<Number>::value(0.0);
  for (unsigned int i = 0; i < dim; ++i)
    {
      const Number diff = static_cast<Number>(this->values[i]) - p(i);
      sum += numbers::NumberTraits<Number>::abs_square(diff);
    }

  return sum;
}



template <int dim, typename Number>
template <class Archive>
inline void
Point<dim, Number>::serialize(Archive &ar, const unsigned int)
{
  // forward to serialization
  // function in the base class
  ar &static_cast<Tensor<1, dim, Number> &>(*this);
}

#endif // DOXYGEN


/*--------------------------- Global functions: Point -----------------------*/


/**
 * Global operator scaling a point vector by a scalar.
 *
 * @note This function can also be used in CUDA device code.
 *
 * @relatesalso Point
 * @relatesalso EnableIfScalar
 */
template <int dim, typename Number, typename OtherNumber>
inline DEAL_II_CUDA_HOST_DEV
  Point<dim,
        typename ProductType<Number,
                             typename EnableIfScalar<OtherNumber>::type>::type>
  operator*(const OtherNumber factor, const Point<dim, Number> &p)
{
  return p * factor;
}



/**
 * Output operator for points. Print the elements consecutively, with a space
 * in between.
 * @relatesalso Point
 */
template <int dim, typename Number>
inline std::ostream &
operator<<(std::ostream &out, const Point<dim, Number> &p)
{
  for (unsigned int i = 0; i < dim - 1; ++i)
    out << p[i] << ' ';
  out << p[dim - 1];

  return out;
}



/**
 * Input operator for points. Inputs the elements consecutively.
 * @relatesalso Point
 */
template <int dim, typename Number>
inline std::istream &
operator>>(std::istream &in, Point<dim, Number> &p)
{
  for (unsigned int i = 0; i < dim; ++i)
    in >> p[i];

  return in;
}


#ifndef DOXYGEN

/**
 * Output operator for points of dimension 1. This is implemented specialized
 * from the general template in order to avoid a compiler warning that the
 * loop is empty.
 */
template <typename Number>
inline std::ostream &
operator<<(std::ostream &out, const Point<1, Number> &p)
{
  out << p[0];

  return out;
}

#endif // DOXYGEN
DEAL_II_NAMESPACE_CLOSE

#endif
