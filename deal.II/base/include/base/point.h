//----------------------------  point.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal authors
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


/**
 * The @p{Point} class provides for a point or vector in a space with arbitrary
 * dimension @p{dim}.
 *
 * It is the preferred object to be passed to functions which
 * operate on points in spaces of a priori unknown dimension: rather than
 * using functions like @p{double f(double x)} and @p{double f(double x, double y)},
 * you use double @p{f(Point<dim> &p)}.
 *
 * @p{Point} also serves as a starting point for the implementation of the
 * geometrical primitives like @p{Polyhedron}, @p{Triangle}, etc.
 *
 * @p{Point}s can also be thought of as vectors, i.e. points in a vector space
 * without an obvious meaning. For instance, it may be suitable to let the
 * gradient of a function be a @p{point} vector:
 * @p{Point<dim> gradient_of_f (const Point<dim> &x)}. @p{Point}s have all
 * functionality for this, e.g. scalar products, addition etc. However, you
 * should also consider using the simpler @p{Tensor<1,dim>} class, which seems
 * more suited to gradients.
 *
 * @author Wolfgang Bangerth, 1997
 */
template <int dim>
class Point : public Tensor<1,dim>
{
  public:
				     /**
				      * Constructor. Initialize all entries
				      * to zero if @p{initialize==true}; this
				      * is the default behaviour.
				      */
    explicit Point (const bool initialize = true);

				     /**
				      * Convert a tensor to a point. Since no
				      * additional data is inside a point,
				      * this is ok.
				      */
    Point (const Tensor<1,dim> &);
    
				     /**
				      *  Constructor for one dimensional points. This
				      *  function is only implemented for @p{dim==1}
				      *  since the usage is considered unsafe
				      *  for points with @p{dim!=1}.
				      */
    explicit Point (const double x);

				     /**
				      *  Constructor for two dimensional points. This
				      *  function is only implemented for @p{dim==2}
				      *  since the usage is considered unsafe
				      *  for points with @p{dim!=2}.
				      */
    Point (const double x, const double y);
    
				     /**
				      *  Constructor for three dimensional points. This
				      *  function is only implemented for @p{dim==3}
				      *  since the usage is considered unsafe
				      *  for points with @p{dim!=3}.
				      */
    Point (const double x, const double y, const double z);

				     /**
				      *  Read access to the @p{index}th coordinate.
				      */
    double   operator () (const unsigned int index) const;

    				     /**
				      *  Read and write access to the @p{index}th
				      *  coordinate.
				      */
    double & operator () (const unsigned int index);

/*
 * Plus and minus operators are re-implemented from Tensor<1,dim>
 * to avoid additional casting.
 */
					 
				     /**
				      *  Add two point vectors. If possible, use
				      *  @p{operator +=} instead since this does not
				      *  need to copy a point at least once.
				      */
    Point<dim>   operator + (const Tensor<1,dim> &) const;

				     /**
				      *  Subtract two point vectors. If possible, use
				      *  @p{operator +=} instead since this does not
				      *  need to copy a point at least once.
				      */
    Point<dim>   operator - (const Tensor<1,dim> &) const;

				     /**
				      * The opposite vector.
				      */
    Point<dim>   operator - () const;
    
				     /**
				      *  Multiply by a factor. If possible, use
				      *  @p{operator *=} instead since this does not
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
				      *  @p{operator /=} instead since this does not
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
				      * Returns the distance of @p{this} 
				      * point to the point @p{p}.
				      */
    double distance (const Point<dim> &p) const;


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
};

/*------------------------------- Inline functions: Point ---------------------------*/


template <int dim>
inline
Point<dim>::Point (const bool initialize) :
		Tensor<1,dim>(initialize) 
{};



template <int dim>
inline
Point<dim>::Point (const Tensor<1,dim> &t) :
		Tensor<1,dim>(t) 
{};



template <int dim>
inline
Point<dim>::Point (const double x)
{
  Assert (dim==1, ExcInvalidConstructorCalled());
  values[0] = x;
};



template <int dim>
inline
Point<dim>::Point (const double x, const double y)
{
  Assert (dim==2, ExcInvalidConstructorCalled());
  values[0] = x;
  values[1] = y;
};



template <int dim>
inline
Point<dim>::Point (const double x, const double y, const double z)
{
  Assert (dim==3, ExcInvalidConstructorCalled());
  values[0] = x;
  values[1] = y;
  values[2] = z;
};



template <int dim>
inline
double Point<dim>::operator () (const unsigned int index) const
{
  Assert (index<dim, ExcIndexRange (index, 0, dim));
  return values[index];
};



template <int dim>
inline
double & Point<dim>::operator () (const unsigned int index) 
{
  Assert (index<dim, ExcIndexRange (index, 0, dim));
  return values[index];
};



template <int dim>
inline
Point<dim> Point<dim>::operator + (const Tensor<1,dim> &p) const 
{
  return (Point<dim>(*this) += p);
};



template <int dim>
inline
Point<dim> Point<dim>::operator - (const Tensor<1,dim> &p) const 
{
  return (Point<dim>(*this) -= p);
};



template <int dim>
inline
Point<dim> Point<dim>::operator - () const 
{
  Point<dim> result;
  for (unsigned int i=0; i<dim; ++i)
    result.values[i] = -values[i];
  return result;
};



template <int dim>
inline
Point<dim> Point<dim>::operator * (const double factor) const 
{
  return (Point<dim>(*this) *= factor);
};



template <int dim>
inline
double Point<dim>::operator * (const Tensor<1,dim> &p) const 
{
				   // simply pass down
  return Tensor<1,dim>::operator * (p);
};


template <int dim>
inline
double Point<dim>::square () const 
{
  double q=0;
  for (unsigned int i=0; i<dim; ++i)
    q += values[i] * values[i];
  return q;
};


template <int dim>
inline
double Point<dim>::distance (const Point<dim> &p) const
{
  double sum=0;
  for (unsigned int i=0; i<dim; ++i)
    {
      const double diff=values[i]-p(i);
      sum += diff*diff;
    }
  
  return std::sqrt(sum);
}


template <int dim>
inline
Point<dim> Point<dim>::operator / (const double factor) const 
{
  return (Point<dim>(*this) /= factor);
};



/*------------------------------- Global functions: Point ---------------------------*/


/**
 * Global operator scaling a point vector by a scalar.
 */
template <int dim>
inline
Point<dim> operator * (const double factor, const Point<dim> &p) 
{
  return p*factor;
};


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
};



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
};

#endif
