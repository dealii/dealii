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

#ifndef __deal2__point_h
#define __deal2__point_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/tensor_base.h>
#include <cmath>

DEAL_II_NAMESPACE_OPEN

/**
 * The <tt>Point</tt> class represents a point in a space with
 * arbitrary dimension <tt>dim</tt>.
 *
 * It is the preferred object to be passed to functions which operate on
 * points in spaces of a priori fixed dimension: rather than using functions
 * like <tt>double f(double x)</tt> and <tt>double f(double x, double y)</tt>,
 * you should use <tt>double f(Point<dim> &p)</tt> instead as it allows writing
 * dimension independent code.
 *
 *
 * <h3>What's a <code>Point@<dim@></code> and what is a <code>Tensor@<1,dim@></code>?</h3>
 *
 * The Point class is derived from Tensor@<1,dim@> and consequently
 * shares the latter's member functions and other attributes. In fact,
 * it has relatively few additional functions itself (the most notable
 * exception being the distance() function to compute the Euclidean
 * distance between two points in space), and these two classes can
 * therefore often be used interchangeably.
 *
 * Nonetheless, there are semantic differences that make us use these
 * classes in different and well-defined contexts. Within deal.II, we
 * use the <tt>Point</tt> class to denote points in space, i.e., for
 * vectors (rank-1 tensors) that are <em>anchored at the
 * origin</em>. On the other hand, vectors that are anchored elsewhere
 * (and consequently do not represent <em>points</em> in the common
 * usage of the word) are represented by objects of type
 * Tensor@<1,dim@>. In particular, this is the case for direction
 * vectors, normal vectors, gradients, and the differences between two
 * points (i.e., what you get when you subtract one point from
 * another): all of these are represented by Tensor@<1,dim@> objects
 * rather than Point@<dim@>.
 *
 * Furthermore, the Point class is only used where the coordinates of
 * an object can be thought to possess the dimension of a length. An
 * object that represents the weight, height, and cost of an object is
 * neither a point nor a tensor (because it lacks the transformation
 * properties under rotation of the coordinate system) and should
 * consequently not be represented by either of these classes. Use an
 * array of size 3 in this case, or the <code>std_cxx11::array</code>
 * class. Alternatively, as in the case of vector-valued functions,
 * you can use objects of type Vector or <code>std::vector<code>.
 *
 * @ingroup geomprimitives
 * @author Wolfgang Bangerth, 1997
 */
template <int dim>
class Point : public Tensor<1,dim,double>
{
public:
  /**
   * Standard constructor. Creates an object that corresponds to the origin,
   * i.e., all coordinates are set to zero.
   */
  Point ();

  /**
   * Constructor. Initialize all entries to zero if <tt>initialize==true</tt>
   * (in which case it is equivalent to the default constructor) or leaves the
   * elements uninitialized if <tt>initialize==false</tt>.
   */
  explicit Point (const bool initialize);

  /**
   * Convert a tensor to a point.
   */
  explicit Point (const Tensor<1,dim> &);

  /**
   * Constructor for one dimensional points. This function is only implemented
   * for <tt>dim==1</tt> since the usage is considered unsafe for points with
   * <tt>dim!=1</tt> as it would leave some components of the point
   * coordinates uninitialized.
   */
  explicit Point (const double x);

  /**
   * Constructor for two dimensional points. This function is only implemented
   * for <tt>dim==2</tt> since the usage is considered unsafe for points with
   * <tt>dim!=2</tt> as it would leave some components of the point
   * coordinates uninitialized (if dim>2) or would not use some arguments (if
   * dim<2).
   */
  Point (const double x,
         const double y);

  /**
   * Constructor for three dimensional points. This function is only
   * implemented for <tt>dim==3</tt> since the usage is considered unsafe for
   * points with <tt>dim!=3</tt> as it would leave some components of the
   * point coordinates uninitialized (if dim>3) or would not use some
   * arguments (if dim<3).
   */
  Point (const double x,
         const double y,
         const double z);

  /**
   * Return a unit vector in coordinate direction <tt>i</tt>.
   */
  static Point<dim> unit_vector(const unsigned int i);

  /**
   * Read access to the <tt>index</tt>th coordinate.
   */
  double   operator () (const unsigned int index) const;

  /**
   * Read and write access to the <tt>index</tt>th coordinate.
   */
  double &operator () (const unsigned int index);

  /*
   * Plus and minus operators are re-implemented from Tensor<1,dim>
   * to avoid additional casting.
   */

  /**
   * Add two point vectors. If possible, use <tt>operator +=</tt> instead
   * since this does not need to copy a point at least once.
   */
  Point<dim>   operator + (const Tensor<1,dim> &) const;

  /**
   * Subtract two points, i.e., obtain the vector that connects the
   * two. As discussed in the documentation of this class, subtracting
   * two points results in a vector anchored at one of the two points
   * (rather than at the origin) and, consequently, the result is
   * returned as a Tensor@<1,dim@> rather than as a Point@<dim@>.
   */
  Tensor<1,dim>   operator - (const Point<dim> &) const;

  /**
   * Subtract a difference vector (represented by a Tensor@<1,dim@>)
   * from the current point. This results in another point and, as
   * discussed in the documentation of this class, the result is then
   * naturally returned as a Point@<dim@> object rather than as a
   * Tensor@<1,dim@>.
   */
  Point<dim>   operator - (const Tensor<1,dim> &) const;

  /**
   * The opposite vector.
   */
  Point<dim>   operator - () const;

  /**
   * Multiply by a factor. If possible, use <tt>operator *=</tt> instead since
   * this does not need to copy a point at least once.
   *
   * There is a commutative complement to this function also
   */
  Point<dim>   operator * (const double) const;

  /**
   * Returns the scalar product of two vectors.
   */
  double       operator * (const Tensor<1,dim> &) const;

  /**
   * Divide by a factor. If possible, use <tt>operator /=</tt> instead since
   * this does not need to copy a point at least once.
   */
  Point<dim>   operator / (const double) const;

  /**
   * Returns the scalar product of this point vector with itself, i.e. the
   * square, or the square of the norm.
   */
  double              square () const;

  /**
   * Returns the Euclidean distance of <tt>this</tt> point to the point
   * <tt>p</tt>, i.e. the <tt>l_2</tt> norm of the difference between the
   * vectors representing the two points.
   */
  double distance (const Point<dim> &p) const;

  /**
   * Read or write the data of this object to or from a stream for the purpose
   * of serialization
   */
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version);
};

/*------------------------------- Inline functions: Point ---------------------------*/

#ifndef DOXYGEN

template <int dim>
inline
Point<dim>::Point ()
{}



template <int dim>
inline
Point<dim>::Point (const bool initialize)
  :
  Tensor<1,dim>(initialize)
{}



template <int dim>
inline
Point<dim>::Point (const Tensor<1,dim> &t)
  :
  Tensor<1,dim>(t)
{}



template <int dim>
inline
Point<dim>::Point (const double x)
{
  switch (dim)
    {
    case 1:
      this->values[0] = x;
    default:
      Assert (dim==1, StandardExceptions::ExcInvalidConstructorCall());
    }
}



template <int dim>
inline
Point<dim>::Point (const double x, const double y)
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



template <int dim>
inline
Point<dim>::Point (const double x, const double y, const double z)
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


template <int dim>
inline
Point<dim>
Point<dim>::unit_vector(unsigned int i)
{
  Point<dim> p;
  p[i] = 1.;
  return p;
}


template <int dim>
inline
double
Point<dim>::operator () (const unsigned int index) const
{
  AssertIndexRange((int) index, dim);
  return this->values[index];
}



template <int dim>
inline
double &
Point<dim>::operator () (const unsigned int index)
{
  AssertIndexRange((int) index, dim);
  return this->values[index];
}



template <int dim>
inline
Point<dim>
Point<dim>::operator + (const Tensor<1,dim> &p) const
{
  Point<dim> tmp = *this;
  tmp += p;
  return tmp;
}



template <int dim>
inline
Tensor<1,dim>
Point<dim>::operator - (const Point<dim> &p) const
{
  return (Tensor<1,dim>(*this) -= p);
}



template <int dim>
inline
Point<dim>
Point<dim>::operator - (const Tensor<1,dim> &p) const
{
  Point<dim> tmp = *this;
  tmp -= p;
  return tmp;
}



template <int dim>
inline
Point<dim>
Point<dim>::operator - () const
{
  Point<dim> result;
  for (unsigned int i=0; i<dim; ++i)
    result.values[i] = -this->values[i];
  return result;
}



template <int dim>
inline
Point<dim>
Point<dim>::operator * (const double factor) const
{
  Point<dim> tmp = *this;
  tmp *= factor;
  return tmp;
}



template <int dim>
inline
double
Point<dim>::operator * (const Tensor<1,dim> &p) const
{
  // simply pass down
  return Tensor<1,dim>::operator * (p);
}


template <int dim>
inline
double
Point<dim>::square () const
{
  double q = double();
  for (unsigned int i=0; i<dim; ++i)
    q += this->values[i] * this->values[i];
  return q;
}



template <int dim>
inline
double
Point<dim>::distance (const Point<dim> &p) const
{
  double sum=0;
  for (unsigned int i=0; i<dim; ++i)
    {
      const double diff=this->values[i]-p(i);
      sum += diff*diff;
    }

  return std::sqrt(sum);
}



template <int dim>
inline
Point<dim> Point<dim>::operator / (const double factor) const
{
  Point<dim> tmp = *this;
  tmp /= factor;
  return tmp;
}



template <int dim>
template <class Archive>
inline
void
Point<dim>::serialize(Archive &ar, const unsigned int)
{
  // forward to serialization
  // function in the base class
  ar   &static_cast<Tensor<1,dim> &>(*this);
}

#endif // DOXYGEN


/*------------------------------- Global functions: Point ---------------------------*/


/**
 * Global operator scaling a point vector by a scalar. @relates Point
 */
template <int dim>
inline
Point<dim> operator * (const double      factor,
                       const Point<dim> &p)
{
  return p*factor;
}



/**
 * Output operator for points. Print the elements consecutively, with a space
 * in between. @relates Point
 */
template <int dim>
inline
std::ostream &operator << (std::ostream     &out,
                           const Point<dim> &p)
{
  for (unsigned int i=0; i<dim-1; ++i)
    out << p[i] << ' ';
  out << p[dim-1];

  return out;
}



/**
 * Output operator for points. Print the elements consecutively, with a space
 * in between. @relates Point
 */
template <int dim>
inline
std::istream &operator >> (std::istream &in,
                           Point<dim>   &p)
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
inline
std::ostream &operator << (std::ostream &out,
                           const Point<1> &p)
{
  out << p[0];

  return out;
}

#endif // DOXYGEN
DEAL_II_NAMESPACE_CLOSE

#endif
