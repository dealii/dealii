/*----------------------------   point.h     ---------------------------*/
/*      $Id$                 */
/*      Copyright W. Bangerth, University of Heidelberg, 1998 */
#ifndef __point_H
#define __point_H
/*----------------------------   point.h     ---------------------------*/


#include <base/exceptions.h>
#include <base/tensor_base.h>


/**
 * The #Point# class provides for a point or vector in a space with arbitrary
 * dimension #dim#.
 *
 * It is the preferred object to be passed to functions which
 * operate on points in spaces of a priori unknown dimension: rather than
 * using functions like #double f(double x)# and #double f(double x, double y)#,
 * you use double #f(Point<dim> &p)#.
 *
 * #Point# also serves as a starting point for the implementation of the
 * geometrical primitives like #Polyhedron#, #Triangle#, etc.
 *
 * #Point#s can also be thought of as vectors, i.e. points in a vector space
 * without an obvious meaning. For instance, it may be suitable to let the
 * gradient of a function be a #point# vector:
 * #Point<dim> gradient_of_f (const Point<dim> &x)#. #Point#s have all
 * functionality for this, e.g. scalar products, addition etc. However, you
 * should also consider using the simpler #Tensor<1,dim># class, which seems
 * more suited to gradients.
 *
 * @author Wolfgang Bangerth, 1997
 */
template <int dim>
class Point : public Tensor<1,dim> {
  public:
				     /**
				      * Constructor. Initialize all entries
				      * to zero if #initialize==true#; this
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
				      *  function is only implemented for #dim==1#
				      *  since the usage is considered unsafe
				      *  for points with #dim!=1#.
				      */
    explicit Point (const double x);

				     /**
				      *  Constructor for two dimensional points. This
				      *  function is only implemented for #dim==2#
				      *  since the usage is considered unsafe
				      *  for points with #dim!=2#.
				      */
    Point (const double x, const double y);
    
				     /**
				      *  Constructor for three dimensional points. This
				      *  function is only implemented for #dim==3#
				      *  since the usage is considered unsafe
				      *  for points with #dim!=3#.
				      */
    Point (const double x, const double y, const double z);

				     /**
				      *  Read access to the #index#th coordinate.
				      */
    double   operator () (const unsigned int index) const;

    				     /**
				      *  Read and write access to the #index#th
				      *  coordinate.
				      */
    double & operator () (const unsigned int index);

				     /**
				      *  Add two point vectors. If possible, use
				      *  #operator +=# instead since this does not
				      *  need to copy a point at least once.
				      */
    Point<dim>   operator + (const Point<dim> &) const;

				     /**
				      *  Subtract two point vectors. If possible, use
				      *  #operator +=# instead since this does not
				      *  need to copy a point at least once.
				      */
    Point<dim>   operator - (const Point<dim> &) const;

				     /**
				      *  Multiply by a factor. If possible, use
				      *  #operator *=# instead since this does not
				      *  need to copy a point at least once.
				      *
				      * There is a commutative complement to this
				      * function also
				      */
    Point<dim>   operator * (const double) const;

				     /**
				      *  Returns the scalar product of two vectors.
				      */
    double       operator * (const Point<dim> &) const;

				     /**
				      *  Divide by a factor. If possible, use
				      *  #operator /=# instead since this does not
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
				      *  Exception
				      */
    DeclException1 (ExcDimTooSmall,
		    int,
		    << "Given dimensions must be >= 1, but was " << arg1);
				     /**
				      *  Exception
				      */
    DeclException1 (ExcInvalidIndex,
		    int,
		    << "Invalid index " << arg1);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidConstructorCalled);
};







/*------------------------------- Inline functions: Point ---------------------------*/



template <int dim>
inline
Point<dim>::Point (const bool initialize) :
		Tensor<1,dim>(initialize) {};


template <int dim>
inline
Point<dim>::Point (const Tensor<1,dim> &t) :
		Tensor<1,dim>(t) {};



template <>
inline
Point<1>::Point (const double x) {
  values[0] = x;
};



template <>
inline
Point<1>::Point (const double, const double) {
  Assert (false, ExcInvalidConstructorCalled());
};



template <>
inline
Point<1>::Point (const double, const double, const double) {
  Assert (false, ExcInvalidConstructorCalled());
};



template <>
inline
Point<2>::Point (const double x, const double y) {
  values[0] = x;
  values[1] = y;
};



template <>
inline
Point<2>::Point (const double, const double, const double) {
  Assert (false, ExcInvalidConstructorCalled());
};



template <>
inline
Point<3>::Point (const double x, const double y, const double z) {
  values[0] = x;
  values[1] = y;
  values[2] = z;
};



template <int dim>
inline
double Point<dim>::operator () (const unsigned int index) const {
  Assert (index<dim, ExcInvalidIndex (index));
  return values[index];
};



template <int dim>
inline
double & Point<dim>::operator () (const unsigned int index) {
  Assert (index<dim, ExcInvalidIndex (index));
  return values[index];
};



template <int dim>
inline
Point<dim> Point<dim>::operator + (const Point<dim> &p) const {
  return (Point<dim>(*this) += p);
};



template <int dim>
inline
Point<dim> Point<dim>::operator - (const Point<dim> &p) const {
  return (Point<dim>(*this) -= p);
};



template <int dim>
inline
Point<dim> Point<dim>::operator * (const double factor) const {
  return (Point<dim>(*this) *= factor);
};



template <int dim>
inline
double Point<dim>::operator * (const Point<dim> &p) const {
				   // simply pass down
  return Tensor<1,dim>::operator * (p);
};



template <int dim>
inline
Point<dim> operator * (const double factor, const Point<dim> &p) {
  return p*factor;
};



template <int dim>
inline
Point<dim> Point<dim>::operator / (const double factor) const {
  return (Point<dim>(*this) /= factor);
};


  
template <int dim>
inline
double Point<dim>::square () const {
  double q=0;
  for (unsigned int i=0; i<dim; ++i)
    q += values[i] * values[i];
  return q;
};





/*----------------------------   point.h     ---------------------------*/
/* end of #ifndef __point_H */
#endif
/*----------------------------   point.h     ---------------------------*/
