/*----------------------------   function.h     ---------------------------*/
/*      $Id$                 */
/*      Copyright W. Bangerth, University of Heidelberg, 1998 */
#ifndef __function_H
#define __function_H
/*----------------------------   function.h     ---------------------------*/


#include <base/exceptions.h>
#include <vector>
#include <grid/point.h>





/**
 *  This class is a model for a continuous function. It returns the value
 *  of the function at a given point through the #operator ()# member function,
 *  which is virtual. It also has a function to return a whole list of function
 *  values at different points to reduce the overhead of the virtual function
 *  calls; this function is preset to successively call the function returning
 *  one value at a time.
 *
 *  There are other functions return the gradient of the function at one or
 *  several points. You only have to overload those functions you need; the
 *  functions returning several values at a time will call those returning
 *  only one value, while those ones will throw an exception when called but
 *  not overloaded.
 *
 *  Unless only called a very small number of times, you should overload
 *  both those functions returning only one value as well as those returning
 *  a whole array, since the cost of evaluation of a point value is often
 *  less than the virtual function call itself.
 *
 *
 *  \subsection{Support for time dependant functions}
 *
 *  The library was also designed for time dependant problems. For this
 *  purpose, the function objects also contain a field which stores the
 *  time, as well as functions manipulating them. Time independant problems
 *  should not access or even abuse them for other purposes, but since one
 *  normally does not create thousands of function objects, the gain in
 *  generality weighs out the fact that we need not store the time value
 *  for not time dependant problems. The second advantage is that the derived
 *  standard classes like #ZeroFunction#, #ConstantFunction# etc also work
 *  for time dependant problems.
 *
 *  Access to the time goes through the following functions:
 *  \begin{verbatim}
 *  \item #get_time#: return the present value of the time variable.
 *  \item #set_time#: set the time value to a specific value.
 *  \item #advance_time#: increase the time by a certain time step.
 *  \end{verbatim}
 *  The latter two functions are virtual, so that derived classes can
 *  perform computations which need only be done once for every new time.
 *  For example, if a time dependant function had a factor #sin(t)#, then
 *  it may be a reasonable choice to calculate this factor in a derived
 *  version of #set_time#, store it in a member variable and use that one
 *  rather than computing it every time #operator()#, #value_list# or one
 *  of the other functions is called.
 *
 *  By default, the #advance_time# function calls the #set_time# function
 *  with the new time, so it is sufficient in most cases to overload only
 *  #set_time# for computations as sketched out above.
 *
 *  Derived classes should access the time variable directly, which is
 *  available under the obvious name #time#, rather than calling
 *  #get_time#.
 *
 *  The constructor of this class takes an initial value for the time
 *  variable, which defaults to zero. Because a default value is given,
 *  none of the derived classes needs to take an initial value for the
 *  time variable if not needed.
 *
 *  Once again the warning: do not use the #time# variable for any other
 *  purpose than the intended one! This will inevitably lead to confusion.
 *
 *
 *  @author Wolfgang Bangerth, 1998
 */
template <int dim>
class Function {
  public:
				     /**
				      * Constructor. May take an initial vakue
				      * for the time variable, which defaults
				      * to zero.
				      */
    Function (const double initial_time = 0.0);
    
				     /**
				      * Virtual destructor; absolutely
				      * necessary in this case.
				      */
    virtual ~Function ();
    
				     /**
				      * Return the value of the function
				      * at the given point.
				      */
    virtual double operator () (const Point<dim> &p) const;

				     /**
				      * Set #values# to the point values
				      * of the function at the #points#.
				      * It is assumed that #values#
				      * already has the right size, i.e.
				      * the same size as the #points#
				      * array.
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
				      * It is assumed that #values# 
				      * already has the right size, i.e.
				      * the same size as the #points# array.
				      */
    virtual void gradient_list (const vector<Point<dim> > &points,
				vector<Point<dim> >       &gradients) const;

				     /**
				      * Return the value of the time variable/
				      */
    double get_time () const;

				     /**
				      * Set the time to #new_time#, overwriting
				      * the old value.
				      */
    virtual void set_time (const double new_time);

				     /**
				      * Advance the time by the given
				      * time step #delta_t#.
				      */
    virtual void advance_time (const double delta_t);
    
				     /**
				      * Exception
				      */
    DeclException0 (ExcPureFunctionCalled);
				     /**
				      * Exception
				      */
    DeclException2 (ExcVectorHasWrongSize,
		    int, int,
		    << "The vector has size " << arg1 << " but should have "
		    << arg2 << " elements.");
    
  protected:
				     /**
				      * Store the present time.
				      */
    double time;
};




/**
 *  Provide a function which always returns zero. Obviously, also the derivates
 *  of this function are zero.
 *
 *  This function is of use when you want to implement homogeneous boundary
 *  conditions.
 */
template <int dim>
class ZeroFunction : public Function<dim> {
  public:
				     /**
				      * Virtual destructor; absolutely
				      * necessary in this case.
				      */
    virtual ~ZeroFunction ();
				     /**
				      * Return the value of the function
				      * at the given point.
				      */
    virtual double operator () (const Point<dim> &p) const;

				     /**
				      * Set #values# to the point values
				      * of the function at the #points#.
				      * It is assumed that #values#
				      * already has the right size, i.e.
				      * the same size as the #points#
				      * array.
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
				      * It is assumed that #values#
				      * already has the right size, i.e.
				      * the same size as the #points# array.
				      */
    virtual void gradient_list (const vector<Point<dim> > &points,
				vector<Point<dim> >       &gradients) const;
};





/**
 *  Provide a function which always returns a constant value, which is delivered
 *  upon construction. Obviously, the derivates of this function are zero, which
 *  is why we derive this class from #ZeroFunction#: we then only have to
 *  overload th value functions, not all the derivatives. In some way, it would
 *  be more obvious to do the derivation in the opposite direction, i.e. let
 *  #ZeroFunction# be a more specialized version of #ConstantFunction#; however,
 *  this would be more inefficient, since we could not make use of the fact that
 *  the function value of the #ZeroFunction# is known at compile time and need
 *  not be looked up somewhere in memory.
 */
template <int dim>
class ConstantFunction : public ZeroFunction<dim> {
  public:
				     /**
				      * Constructor; takes the constant function
				      * value as an argument.
				      */
    ConstantFunction (const double value);
    
				     /**
				      * Virtual destructor; absolutely
				      * necessary in this case.
				      */
    virtual ~ConstantFunction ();
				     /**
				      * Return the value of the function
				      * at the given point.
				      */
    virtual double operator () (const Point<dim> &p) const;

				     /**
				      * Set #values# to the point values
				      * of the function at the #points#.
				      * It is assumed that #values#
				      * already has the right size, i.e.
				      * the same size as the #points#
				      * array.
				      */
    virtual void value_list (const vector<Point<dim> > &points,
			     vector<double>            &values) const;

  protected:
				     /**
				      * Store the constant function value.
				      */
    const double function_value;
};




/*----------------------------   function.h     ---------------------------*/
/* end of #ifndef __function_H */
#endif
/*----------------------------   function.h     ---------------------------*/
