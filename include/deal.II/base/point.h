// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2016 by the deal.II authors
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

#ifndef dealii__point_h
#define dealii__point_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/tensor.h>
#include <cmath>

DEAL_II_NAMESPACE_OPEN

/**
 * A class that represents a point in a space with arbitrary dimension
 * <tt>dim</tt>.
 *
 * Objects of this class are used to represent points, i.e., vectors anchored
 * at the origin of a Cartesian vector space. They are, among other uses,
 * passed to functions that operate on points in spaces of a priori fixed
 * dimension: rather than using functions like <tt>double f(double x)</tt> and
 * <tt>double f(double x, double y)</tt>, you should use <tt>double
 * f(Point<dim> &p)</tt> instead as it allows writing dimension independent
 * code.
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
 * <code>std_cxx11::array</code> class. Alternatively, as in the case of
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
class Point : public Tensor<1,dim,Number>
{
public:
  /**
   * Standard constructor. Creates an object that corresponds to the origin,
   * i.e., all coordinates are set to zero.
   */
  Point ();

  /**
   * Convert a tensor to a point.
   */
  explicit Point (const Tensor<1,dim,Number> &);

  /**
   * Constructor for one dimensional points. This function is only implemented
   * for <tt>dim==1</tt> since the usage is considered unsafe for points with
   * <tt>dim!=1</tt> as it would leave some components of the point
   * coordinates uninitialized.
   */
  explicit Point (const Number x);

  /**
   * Constructor for two dimensional points. This function is only implemented
   * for <tt>dim==2</tt> since the usage is considered unsafe for points with
   * <tt>dim!=2</tt> as it would leave some components of the point
   * coordinates uninitialized (if dim>2) or would not use some arguments (if
   * dim<2).
   */
  Point (const Number x,
         const Number y);

  /**
   * Constructor for three dimensional points. This function is only
   * implemented for <tt>dim==3</tt> since the usage is considered unsafe for
   * points with <tt>dim!=3</tt> as it would leave some components of the
   * point coordinates uninitialized (if dim>3) or would not use some
   * arguments (if dim<3).
   */
  Point (const Number x,
         const Number y,
         const Number z);

  /**
   * Return a unit vector in coordinate direction <tt>i</tt>, i.e., a vector
   * that is zero in all coordinates except for a single 1 in the <tt>i</tt>th
   * coordinate.
   */
  static Point<dim,Number> unit_vector(const unsigned int i);

  /**
   * Read access to the <tt>index</tt>th coordinate.
   */
  Number operator () (const unsigned int index) const;

  /**
   * Read and write access to the <tt>index</tt>th coordinate.
   */
  Number &operator () (const unsigned int index);

  /*
   * @name Addition and subtraction of points.
   * @{
   */

  /**
   * Add an offset given as Tensor<1,dim,Number> to a point.
   */
  Point<dim,Number> operator + (const Tensor<1,dim,Number> &) const;

  /**
   * Subtract two points, i.e., obtain the vector that connects the two. As
   * discussed in the documentation of this class, subtracting two points
   * results in a vector anchored at one of the two points (rather than at the
   * origin) and, consequently, the result is returned as a Tensor@<1,dim@>
   * rather than as a Point@<dim@>.
   */
  Tensor<1,dim,Number> operator - (const Point<dim,Number> &) const;

  /**
   * Subtract a difference vector (represented by a Tensor@<1,dim@>) from the
   * current point. This results in another point and, as discussed in the
   * documentation of this class, the result is then naturally returned as a
   * Point@<dim@> object rather than as a Tensor@<1,dim@>.
   */
  Point<dim,Number> operator - (const Tensor<1,dim,Number> &) const;

  /**
   * The opposite vector.
   */
  Point<dim,Number> operator - () const;

  /**
   * @}
   */

  /*
   * @name Multiplication and scaling of points. Dot products. Norms.
   * @{
   */

  /**
   * Multiply the current point by a factor.
   *
   * @relates EnableIfScalar
   */
  template <typename OtherNumber>
  Point<dim,typename ProductType<Number, typename EnableIfScalar<OtherNumber>::type>::type>
  operator * (const OtherNumber) const;

  /**
   * Divide the current point by a factor.
   */
  template <typename OtherNumber>
  Point<dim,typename ProductType<Number, typename EnableIfScalar<OtherNumber>::type>::type>
  operator / (const OtherNumber) const;

  /**
   * Return the scalar product of the vectors representing two points.
   */
  Number operator * (const Tensor<1,dim,Number> &p) const;

  /**
   * Return the scalar product of this point vector with itself, i.e. the
   * square, or the square of the norm. In case of a complex number type it is
   * equivalent to the contraction of this point vector with a complex
   * conjugate of itself.
   *
   * @note This function is equivalent to
   * Tensor<rank,dim,Number>::norm_square() which returns the square of the
   * Frobenius norm.
   */
  typename numbers::NumberTraits<Number>::real_type square () const;

  /**
   * Return the Euclidean distance of <tt>this</tt> point to the point
   * <tt>p</tt>, i.e. the <tt>l_2</tt> norm of the difference between the
   * vectors representing the two points.
   */
  typename numbers::NumberTraits<Number>::real_type distance (const Point<dim,Number> &p) const;

  /**
   * @}
   */

  /**
   * Read or write the data of this object to or from a stream for the purpose
   * of serialization
   */
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version);
};

/*------------------------------- Inline functions: Point ---------------------------*/

#ifndef DOXYGEN

template <int dim, typename Number>
inline
Point<dim,Number>::Point ()
{}



template <int dim, typename Number>
inline
Point<dim,Number>::Point (const Tensor<1,dim,Number> &t)
  :
  Tensor<1,dim,Number>(t)
{}



template <int dim, typename Number>
inline
Point<dim,Number>::Point (const Number x)
{
  switch (dim)
    {
    case 1:
      this->values[0] = x;
    default:
      Assert (dim==1, StandardExceptions::ExcInvalidConstructorCall());
    }
}



template <int dim, typename Number>
inline
Point<dim,Number>::Point (const Number x, const Number y)
{
  switch (dim)
    {
    case 2:
      this->values[0] = x;
      this->values[1] = y;
    default:
      Assert (dim==2, StandardExceptions::ExcInvalidConstructorCall());
    }
}



template <int dim, typename Number>
inline
Point<dim,Number>::Point (const Number x, const Number y, const Number z)
{
  switch (dim)
    {
    case 3:
      this->values[0] = x;
      this->values[1] = y;
      this->values[2] = z;
    default:
      Assert (dim==3, StandardExceptions::ExcInvalidConstructorCall());
    }
}


template <int dim, typename Number>
inline
Point<dim,Number>
Point<dim,Number>::unit_vector(unsigned int i)
{
  Point<dim,Number> p;
  p[i] = 1.;
  return p;
}


template <int dim, typename Number>
inline
Number
Point<dim,Number>::operator () (const unsigned int index) const
{
  AssertIndexRange((int) index, dim);
  return this->values[index];
}



template <int dim, typename Number>
inline
Number &
Point<dim,Number>::operator () (const unsigned int index)
{
  AssertIndexRange((int) index, dim);
  return this->values[index];
}



template <int dim, typename Number>
inline
Point<dim,Number>
Point<dim,Number>::operator + (const Tensor<1,dim,Number> &p) const
{
  Point<dim,Number> tmp = *this;
  tmp += p;
  return tmp;
}



template <int dim, typename Number>
inline
Tensor<1,dim,Number>
Point<dim,Number>::operator - (const Point<dim,Number> &p) const
{
  return (Tensor<1,dim,Number>(*this) -= p);
}



template <int dim, typename Number>
inline
Point<dim,Number>
Point<dim,Number>::operator - (const Tensor<1,dim,Number> &p) const
{
  Point<dim,Number> tmp = *this;
  tmp -= p;
  return tmp;
}



template <int dim, typename Number>
inline
Point<dim,Number>
Point<dim,Number>::operator - () const
{
  Point<dim,Number> result;
  for (unsigned int i=0; i<dim; ++i)
    result.values[i] = -this->values[i];
  return result;
}



template <int dim, typename Number>
template<typename OtherNumber>
inline
Point<dim,typename ProductType<Number, typename EnableIfScalar<OtherNumber>::type>::type>
Point<dim,Number>::operator * (const OtherNumber factor) const
{
  Point<dim,typename ProductType<Number, OtherNumber>::type> tmp;
  for (unsigned int i=0; i<dim; ++i)
    tmp[i] = this->operator[](i) * factor;
  return tmp;
}



template <int dim, typename Number>
template<typename OtherNumber>
inline
Point<dim,typename ProductType<Number, typename EnableIfScalar<OtherNumber>::type>::type>
Point<dim,Number>::operator / (const OtherNumber factor) const
{
  Point<dim,typename ProductType<Number, OtherNumber>::type> tmp;
  for (unsigned int i=0; i<dim; ++i)
    tmp[i] = this->operator[](i) / factor;
  return tmp;
}



template <int dim, typename Number>
inline
Number
Point<dim,Number>::operator * (const Tensor<1,dim,Number> &p) const
{
  Number res = Number();
  for (unsigned int i=0; i<dim; ++i)
    res += this->operator[](i) * p[i];
  return res;
}


template <int dim, typename Number>
inline
typename numbers::NumberTraits<Number>::real_type
Point<dim,Number>::square () const
{
  return this->norm_square();
}



template <int dim, typename Number>
inline
typename numbers::NumberTraits<Number>::real_type
Point<dim,Number>::distance (const Point<dim,Number> &p) const
{
  Number sum = Number();
  for (unsigned int i=0; i<dim; ++i)
    {
      const Number diff=this->values[i]-p(i);
      sum += numbers::NumberTraits<Number>::abs_square (diff);
    }

  return std::sqrt(sum);
}



template <int dim, typename Number>
template <class Archive>
inline
void
Point<dim,Number>::serialize(Archive &ar, const unsigned int)
{
  // forward to serialization
  // function in the base class
  ar   &static_cast<Tensor<1,dim,Number> &>(*this);
}

#endif // DOXYGEN


/*------------------------------- Global functions: Point ---------------------------*/


/**
 * Global operator scaling a point vector by a scalar.
 *
 * @relates Point
 * @relates EnableIfScalar
 */
template <int dim, typename Number, typename OtherNumber>
inline
Point<dim,typename ProductType<Number, typename EnableIfScalar<OtherNumber>::type>::type>
operator * (const OtherNumber        factor,
            const Point<dim,Number> &p)
{
  return p * factor;
}



/**
 * Output operator for points. Print the elements consecutively, with a space
 * in between.
 * @relates Point
 */
template <int dim, typename Number>
inline
std::ostream &operator << (std::ostream            &out,
                           const Point<dim,Number> &p)
{
  for (unsigned int i=0; i<dim-1; ++i)
    out << p[i] << ' ';
  out << p[dim-1];

  return out;
}



/**
 * Output operator for points. Print the elements consecutively, with a space
 * in between.
 * @relates Point
 */
template <int dim, typename Number>
inline
std::istream &operator >> (std::istream      &in,
                           Point<dim,Number> &p)
{
  for (unsigned int i=0; i<dim; ++i)
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
inline
std::ostream &operator << (std::ostream &out,
                           const Point<1,Number> &p)
{
  out << p[0];

  return out;
}

#endif // DOXYGEN
DEAL_II_NAMESPACE_CLOSE

#endif
