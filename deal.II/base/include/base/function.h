/*----------------------------   function.h     ---------------------------*/
/*      $Id$                 */
#ifndef __function_H
#define __function_H
/*----------------------------   function.h     ---------------------------*/


#include <base/exceptions.h>
#include <vector.h>
#include <grid/point.h>





/**
   This class is a model for a continuous function. It returns the value
   of the function at a given point through the #operator ()# member function,
   which is virtual. It also has a function to return a whole list of function
   values at different points to reduce the overhead of the virtual function
   calls; this function is preset to successively call the function returning
   one value at a time.

   There are other functions return the gradient of the function at one or
   several points. You only have to overload those functions you need; the
   functions returning several values at a time will call those returning
   only one value, while those ones will throw an exception when called but
   not overloaded.

   Unless only called a very small number of times, you should overload
   both those functions returning only one value as well as those returning
   a whole array, since the cost of evaluation of a point value is often
   less than the virtual function call itself.
*/
template <int dim>
class Function {
  public:
				     /**
				      * Return the value of the function
				      * at the given point.
				      */
    virtual double operator () (const Point<dim> &p) const;

				     /**
				      * Set #values# to the point values
				      * of the function at the #points#.
				      * It is assumed that #values# be
				      * empty.
				      */
    virtual void value_list (const vector<Point<dim> > &points,
			     vector<double>            &values) const;

				     /**
				      * Return the gradient of the function
				      * at the given point.
				      */
    virtual Point<dim> gradient (const Point<dim> &p) const;

				     /**
				      * Set #gradients# to the gradients of
				      * the function at the #points#.
				      * It is assumed that #values# be
				      * empty.
				      */
    virtual void gradient_list (const vector<Point<dim> > &points,
				vector<Point<dim> >       &gradients) const;


				     /**
				      * Exception
				      */
    DeclException0 (ExcPureFunctionCalled);
				     /**
				      * Exception
				      */
    DeclException0 (ExcVectorNotEmpty);
};




/**
   Provide a function which always returns zero. Obviously, also the derivates
   of this function are zero.

   This function is of use when you want to implement homogeneous boundary
   conditions.
*/
template <int dim>
class ZeroFunction : public Function<dim> {
  public:
				     /**
				      * Return the value of the function
				      * at the given point.
				      */
    virtual double operator () (const Point<dim> &p) const;

				     /**
				      * Set #values# to the point values
				      * of the function at the #points#.
				      * It is assumed that #values# be
				      * empty.
				      */
    virtual void value_list (const vector<Point<dim> > &points,
			     vector<double>            &values) const;

				     /**
				      * Return the gradient of the function
				      * at the given point.
				      */
    virtual Point<dim> gradient (const Point<dim> &p) const;

				     /**
				      * Set #gradients# to the gradients of
				      * the function at the #points#.
				      * It is assumed that #values# be
				      * empty.
				      */
    virtual void gradient_list (const vector<Point<dim> > &points,
				vector<Point<dim> >       &gradients) const;
};




/*----------------------------   function.h     ---------------------------*/
/* end of #ifndef __function_H */
#endif
/*----------------------------   function.h     ---------------------------*/
