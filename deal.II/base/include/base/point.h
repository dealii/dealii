//----------------------------  point.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  point.h  ---------------------------
#ifndef __deal2__point_h
#define __deal2__point_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/tensor_base.h>
#include <cmath>


/**
 * The <tt>Point</tt> class provides for a point or vector in a space with arbitrary
 * dimension <tt>dim</tt>.
 *
 * It is the preferred object to be passed to functions which
 * operate on points in spaces of a priori unknown dimension: rather than
 * using functions like <tt>double f(double x)</tt> and <tt>double f(double x, double y)</tt>,
 * you use double <tt>f(Point<dim> &p)</tt>.
 *
 * <tt>Point</tt> also serves as a starting point for the implementation of the
 * geometrical primitives like <tt>Polyhedron</tt>, <tt>Triangle</tt>, etc.
 *
 * <tt>Point</tt>s can also be thought of as vectors, i.e. points in a vector space
 * without an obvious meaning. For instance, it may be suitable to let the
 * gradient of a function be a <tt>point</tt> vector:
 * <tt>Point<dim> gradient_of_f (const Point<dim> &x)</tt>. <tt>Point</tt>s have all
 * functionality for this, e.g. scalar products, addition etc. However, you
 * should also consider using the simpler <tt>Tensor<1,dim></tt> class, which seems
 * more suited to gradients.
 *
 * @author Wolfgang Bangerth, 1997
 */
template <int dim>
class Point : public Tensor<1,dim>
{
  public:
				     /**
				      * Standard constructor. Creates
				      * an origin.
				      */
    Point ();
    
				     /**
				      * Constructor. Initialize all
				      * entries to zero if
				      * <tt>initialize==true</tt>.
				      */
    explicit Point (const bool initialize);

				     /**
				      * Convert a tensor to a point. Since no
				      * additional data is inside a point,
				      * this is ok.
				      */
    Point (const Tensor<1,dim> &);
    
				     /**
				      *  Constructor for one dimensional points. This
				      *  function is only implemented for <tt>dim==1</tt>
				      *  since the usage is considered unsafe
				      *  for points with <tt>dim!=1</tt>.
				      */
    explicit Point (const double x);

				     /**
				      *  Constructor for two dimensional points. This
				      *  function is only implemented for <tt>dim==2</tt>
				      *  since the usage is considered unsafe
				      *  for points with <tt>dim!=2</tt>.
				      */
    Point (const double x, const double y);
    
				     /**
				      *  Constructor for three dimensional points. This
				      *  function is only implemented for <tt>dim==3</tt>
				      *  since the usage is considered unsafe
				      *  for points with <tt>dim!=3</tt>.
				      */
    Point (const double x, const double y, const double z);

				     /**
				      *  Read access to the <tt>index</tt>th coordinate.
				      */
    double   operator () (const unsigned int index) const;

    				     /**
				      *  Read and write access to the <tt>index</tt>th
				      *  coordinate.
				      */
    double & operator () (const unsigned int index);

/*
 * Plus and minus operators are re-implemented from Tensor<1,dim>
 * to avoid additional casting.
 */
					 
				     /**
				      *  Add two point vectors. If possible, use
				      *  <tt>operator +=</tt> instead since this does not
				      *  need to copy a point at least once.
				      */
    Point<dim>   operator + (const Tensor<1,dim>&) const;

				     /**
				      *  Subtract two point vectors. If possible, use
				      *  <tt>operator +=</tt> instead since this does not
				      *  need to copy a point at least once.
				      */
    Point<dim>   operator - (const Tensor<1,dim>&) const;

				     /**
				      * The opposite vector.
				      */
    Point<dim>   operator - () const;
    
				     /**
				      *  Multiply by a factor. If possible, use
				      *  <tt>operator *=</tt> instead since this does not
				      *  need to copy a point at least once.
				      *
				      * There is a commutative complement to this
				      * function also
				      */
    Point<dim>   operator * (const double) const;

				     /**
				      *  Returns the scalar product of two vectors.
				      */
    double       operator * (const Tensor<1,dim> &) const;

				     /**
				      *  Divide by a factor. If possible, use
				      *  <tt>operator /=</tt> instead since this does not
				      *  need to copy a point at least once.
				      */
    Point<dim>   operator / (const double) const;

				     /**
				      *  Returns the scalar product of this point
				      *  vector with itself, i.e. the square, or
				      *  the square of the norm.
				      */
    double              square () const;
    
				     /**
				      * Returns the distance of <tt>this</tt> 
				      * point to the point <tt>p</tt>.
				      */
    double distance (const Point<dim> &p) const;

      				     /** @addtogroup Exceptions
				      * @{ */

				     /**
				      *  Exception
				      */
    DeclException1 (ExcDimTooSmall,
		    int,
		    << "Given dimensions must be >= 1, but was " << arg1);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidConstructorCalled);
				     //@}

};

/*------------------------------- Inline functions: Point ---------------------------*/


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
  Assert (dim==1, ExcInvalidConstructorCalled());
  this->values[0] = x;
}



template <int dim>
inline
Point<dim>::Point (const double x, const double y)
{
  Assert (dim==2, ExcInvalidConstructorCalled());
  this->values[0] = x;
  this->values[1] = y;
}



template <int dim>
inline
Point<dim>::Point (const double x, const double y, const double z)
{
  Assert (dim==3, ExcInvalidConstructorCalled());
  this->values[0] = x;
  this->values[1] = y;
  this->values[2] = z;
}



template <int dim>
inline
double Point<dim>::operator () (const unsigned int index) const
{
  Assert (index<dim, ExcIndexRange (index, 0, dim));
  return this->values[index];
}



template <int dim>
inline
double & Point<dim>::operator () (const unsigned int index)
{
  Assert (index<dim, ExcIndexRange (index, 0, dim));
  return this->values[index];
}



template <int dim>
inline
Point<dim> Point<dim>::operator + (const Tensor<1,dim> &p) const
{
  return (Point<dim>(*this) += p);
}



template <int dim>
inline
Point<dim> Point<dim>::operator - (const Tensor<1,dim> &p) const
{
  return (Point<dim>(*this) -= p);
}



template <int dim>
inline
Point<dim> Point<dim>::operator - () const
{
  Point<dim> result;
  for (unsigned int i=0; i<dim; ++i)
    result.values[i] = -this->values[i];
  return result;
}



template <int dim>
inline
Point<dim> Point<dim>::operator * (const double factor) const
{
  return (Point<dim>(*this) *= factor);
}



template <int dim>
inline
double Point<dim>::operator * (const Tensor<1,dim> &p) const
{
				   // simply pass down
  return Tensor<1,dim>::operator * (p);
}


template <int dim>
inline
double Point<dim>::square () const
{
  double q=0;
  for (unsigned int i=0; i<dim; ++i)
    q += this->values[i] * this->values[i];
  return q;
}


template <int dim>
inline
double Point<dim>::distance (const Point<dim> &p) const
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
  return (Point<dim>(*this) /= factor);
}



/*------------------------------- Global functions: Point ---------------------------*/


/**
 * Global operator scaling a point vector by a scalar.
 */
template <int dim>
inline
Point<dim> operator * (const double factor, const Point<dim> &p)
{
  return p*factor;
}


/**
 * Output operator for points. Print the elements consecutively,
 * with a space in between.
 */
template <int dim>
inline
std::ostream & operator << (std::ostream &out, const Point<dim> &p)
{
  for (unsigned int i=0; i<dim-1; ++i)
    out << p[i] << ' ';
  out << p[dim-1];

  return out;
}



/** 
 * Output operator for points of dimension 1. This is implemented
 * specialized from the general template in order to avoid a compiler
 * warning that the loop is empty.  
 */
inline
std::ostream & operator << (std::ostream &out, const Point<1> &p)
{
  out << p[0];

  return out;
}

#endif
