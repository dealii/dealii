/*----------------------------   point.h     ---------------------------*/
/*      $Id$                 */
#ifndef __point_H
#define __point_H
/*----------------------------   point.h     ---------------------------*/


#include <base/exceptions.h>
#include <iostream>


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
 * functionality for this, e.g. scalar products, addition etc.
 *
 * @author Wolfgang Bangerth, 1997
 */
template <int dim>
class Point {
  public:
				     /**
				      *  Constructor. Initialize all entries
				      * to zero.
				      */
    explicit Point ();
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
				      *  Copy constructor.
				      */
    Point (const Point<dim> &);

				     /**
				      *  Return the dimension of the space this
				      *  point is living in.
				      */
    unsigned int dimension() const { return dim; };
    
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
				      *  Assignment operator.
				      */
    Point<dim> & operator = (const Point<dim> &);

				     /**
				      *  Test for equality of two points.
				      */
    bool             operator ==(const Point<dim> &) const;
    				     /**
				      *  Test for inequality of two points.
				      */
    bool             operator !=(const Point<dim> &) const;

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
				      *  Divide by a factor. If possible, use
				      *  #operator /=# instead since this does not
				      *  need to copy a point at least once.
				      */
    Point<dim>   operator / (const double) const;

				     /**
				      *  Add another vector, i.e. move this point by
				      *  the given offset.
				      */
    Point<dim> & operator += (const Point<dim> &);
				     /**
				      *  Subtract another vector.
				      */
    Point<dim> & operator -= (const Point<dim> &);

				     /**
				      *  Scale the vector by #factor#, i.e. multiply
				      *  all coordinates by #factor#.
				      */
    Point<dim> & operator *= (const double &factor);

				     /**
				      *  Scale the vector by #1/factor#.
				      */
    Point<dim> & operator /= (const double &factor);

				     /**
				      *  Returns the scalar product of two vectors.
				      */
    double              operator * (const Point<dim> &) const;

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
  protected:
				     /**
				      *  Stores the coordinate values.
				      */
    double coordinates[dim];
};

				 /**
				  *  Prints the coordinates of this point in the
				  *  form #x1 x2 x3 etc#.
				  */
template <int dim>
ostream & operator << (ostream &out, const Point<dim> &p);







/*------------------------------- Inline functions: Point ---------------------------*/


template <int dim>
inline
Point<dim>::Point () {
  Assert (dim>0, ExcDimTooSmall(dim));

  for (unsigned int i=0; i<dim; ++i)
    coordinates[i] = 0;
};



template <>
inline
Point<1>::Point (const double x) {
  coordinates[0] = x;
};



template <>
inline
Point<2>::Point (const double x, const double y) {
  coordinates[0] = x;
  coordinates[1] = y;
};



template <>
inline
Point<3>::Point (const double x, const double y, const double z) {
  coordinates[0] = x;
  coordinates[1] = y;
  coordinates[2] = z;
};



template <int dim>
inline
Point<dim>::Point (const Point<dim> &p) {
  for (unsigned int i=0; i<dim; ++i)
    coordinates[i] = p.coordinates[i];
};



template <int dim>
inline
double Point<dim>::operator () (const unsigned int index) const {
  Assert (index<dim, ExcInvalidIndex (index));
  return coordinates[index];
};



template <int dim>
inline
double & Point<dim>::operator () (const unsigned int index) {
  Assert (index<dim, ExcInvalidIndex (index));
  return coordinates[index];
};



template <int dim>
inline
Point<dim> & Point<dim>::operator = (const Point<dim> &p) {
  for (unsigned int i=0; i<dim; ++i)
    coordinates[i] = p.coordinates[i];
  return *this;
};



template <int dim>
inline
bool Point<dim>::operator == (const Point<dim> &p) const {
  for (unsigned int i=0; i<dim; ++i)
    if (coordinates[i] != p.coordinates[i]) return false;
  return true;
};



template <int dim>
inline
bool Point<dim>::operator != (const Point<dim> &p) const {
  return !((*this) == p);
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
Point<dim> & Point<dim>::operator += (const Point<dim> &p) {
  for (unsigned int i=0; i<dim; ++i)
    coordinates[i] += p.coordinates[i];
  return *this;
};



template <int dim>
inline
Point<dim> & Point<dim>::operator -= (const Point<dim> &p) {
  for (unsigned int i=0; i<dim; ++i)
    coordinates[i] -= p.coordinates[i];
  return *this;
};



template <int dim>
inline
Point<dim> & Point<dim>::operator *= (const double &s) {
  for (unsigned int i=0; i<dim; ++i)
    coordinates[i] *= s;
  return *this;
};



template <int dim>
inline
Point<dim> & Point<dim>::operator /= (const double &s) {
  for (unsigned int i=0; i<dim; ++i)
    coordinates[i] /= s;
  return *this;
};



template <int dim>
inline
double Point<dim>::operator * (const Point<dim> &p) const {
  double q=0;
  for (unsigned int i=0; i<dim; ++i)
    q += coordinates[i] * p.coordinates[i];
  return q;
};



template <int dim>
inline
double Point<dim>::square () const {
  double q=0;
  for (unsigned int i=0; i<dim; ++i)
    q += coordinates[i] * coordinates[i];
  return q;
};



template <int dim>
inline
ostream & operator << (ostream &out, const Point<dim> &p) {
  for (unsigned int i=0; i<dim-1; ++i)
    out << p(i) << ' ';
  out << p(dim-1);
  return out;
};


template <>
inline
ostream & operator << (ostream &out, const Point<1> &p) {
  out << p(0);
  return out;
};
  




/*----------------------------   point.h     ---------------------------*/
/* end of #ifndef __point_H */
#endif
/*----------------------------   point.h     ---------------------------*/
