//----------------------------  function.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  function.h  ---------------------------
#ifndef __deal2__function_h
#define __deal2__function_h


#include <base/forward_declarations.h>
#include <base/exceptions.h>
#include <base/function_time.h>
#include <base/subscriptor.h>
#include <vector>


/**
 * This class is a model for a general function. It serves the purpose
 * of representing scalar and vector valued functions. To this end, we
 * consider scalar functions as a special case of vector valued
 * functions, in the former case only having a single component return
 * vector. Since handling with vectors is comparatively expensive,
 * functions are provided which only ask for a single component of the
 * function, which is what you will usually need in case you know that
 * your function is scalar-valued.
 *
 * Access to function objects therefore is through the following
 * methods:
 * \begin{verbatim}
 *                 // access to one component at one point
 *   double value        (const Point<dim>   &p,
 *                        const unsigned int  component = 0) const;
 *
 *                 // return all components at one point
 *   void   vector_value (const Point<dim>   &p,
 *                        Vector<double>     &value) const;
 * \end{verbatim}
 *
 * For more efficiency, there are other functions returning one or all
 * components at a list of points at once:
 * \begin{verbatim}
 *                 // access to one component at several points
 *   void   value_list (const vector<Point<dim> >  &point_list,
 *                      vector<double>             &value_list,
 *                      const unsigned int  component = 0) const;
 *
 *                 // return all components at several points
 *   void   vector_value_list (const vector<Point<dim> >    &point_list,
 *                             vector<Vector<double> >      &value_list) const;
 * \end{verbatim}
 *
 * Furthermore, there are functions returning the gradient of the
 * function at one or several points.
 *
 * You will usually only overload those functions you need; the
 * functions returning several values at a time (#value_list#,
 * #vector_value_list#, and gradient analoga) will call those
 * returning only one value (#value#, #vector_value#, and gradient
 * analoga), while those ones will throw an exception when called but
 * not overloaded.
 *
 * Note however, that the functions returning all components of the
 * function at one or several points (i.e. #vector_value#,
 * #vector_value_list#), will not call the function returning one
 * component at one point repeatedly, once for each point and
 * component. The reason is efficiency: this would amount to too many
 * virtual function calls. If you have vector-valued functions, you
 * should therefore also provide overloads of the virtual functions
 * for all components at a time.
 *
 * Also note, that unless only called a very small number of times,
 * you should overload all sets of functions (returning only one
 * value, as well as those returning a whole array), since the cost of
 * evaluation of a point value is often less than the virtual function
 * call itself.
 *
 *
 * Support for time dependant functions can be found in the base
 * class #FunctionTime#.
 *
 * {\bf Note}: if the functions you are dealing with have sizes which
 * are a priori known (for example, #dim# elements), you might
 * consider using the #TensorFunction# class instead.
 *
 * @author Wolfgang Bangerth, 1998, 1999
 */
template <int dim>
class Function : public FunctionTime,
		 public Subscriptor
{
  public:
    				     /**
				      * Number of vector components.
				      */
    const unsigned int n_components;

				     /**
				      * Constructor. May take an
				      * initial value for the number
				      * of components (which defaults
				      * to one, i.e. a scalar
				      * function), and the time
				      * variable, which defaults to
				      * zero.
				      */
    Function (const unsigned int n_components = 1,
	      const double       initial_time = 0.0);
    
				     /**
				      * Virtual destructor; absolutely
				      * necessary in this case.
				      */
    virtual ~Function ();
    
				     /**
				      * Return the value of the
				      * function at the given
				      * point. Unless there is only
				      * one component (i.e. the
				      * function is scalar), you
				      * should state the component you
				      * want to have evaluated; it
				      * defaults to zero, i.e. the
				      * first component.
				      */
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;

				     /**
				      * Return all components of a
				      * vector-valued function at a
				      * given point.
				      *
				      * Be default, this function
				      * repeatedly calls the other
				      * #operator()# for each
				      * component separately, to fill
				      * the output array.
				      *
				      * #values# shall have the right
				      * size beforehand,
				      * i.e. #n_components#.
				      */
    virtual void   vector_value (const Point<dim>   &p,
				 Vector<double>     &values) const;
    
				     /**
				      * Set #values# to the point
				      * values of the specified
				      * component of the function at
				      * the #points#.  It is assumed
				      * that #values# already has the
				      * right size, i.e.  the same
				      * size as the #points# array.
				      */
    virtual void value_list (const vector<Point<dim> > &points,
			     vector<double>            &values,
			     const unsigned int         component = 0) const;

				     /**
				      * Set #values# to the point
				      * values of the function at the
				      * #points#.  It is assumed that
				      * #values# already has the right
				      * size, i.e.  the same size as
				      * the #points# array, and that
				      * all elements be vectors with
				      * the same number of components
				      * as this function has.
				      */
    virtual void vector_value_list (const vector<Point<dim> > &points,
				    vector<Vector<double> >   &values) const;
    
				     /**
				      * Return the gradient of the
				      * specified component of the
				      * function at the given point.
				      */
    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
				    const unsigned int  component = 0) const;

				     /**
				      * Return the gradient of all
				      * components of the
				      * function at the given point.
				      */
    virtual void          vector_gradient (const Point<dim>       &p,
					   vector<Tensor<1,dim> > &gradients) const;
    
				     /**
				      * Set #gradients# to the
				      * gradients of the specified
				      * component of the function at
				      * the #points#.  It is assumed
				      * that #gradients# already has the
				      * right size, i.e.  the same
				      * size as the #points# array.
				      */
    virtual void gradient_list (const vector<Point<dim> > &points,
				vector<Tensor<1,dim> >    &gradients,
				const unsigned int         component = 0) const;
    
				     /**
				      * Set #gradients# to the gradients of
				      * the function at the #points#,
				      * for all components.
				      * It is assumed that #gradients# 
				      * already has the right size, i.e.
				      * the same size as the #points# array.
				      *
				      * The outer loop over
				      * #gradients# is over the points
				      * in the list, the inner loop
				      * over the different components
				      * of the function.
				      */
    virtual void vector_gradient_list (const vector<Point<dim> >       &points,
				       vector<vector<Tensor<1,dim> > > &gradients) const;

				     /**
				      * Compute the Laplacian of a given component at point #p#.
				      */
    virtual double laplacian (const Point<dim>   &p,
			      const unsigned int  component = 0) const;

				     /**
				      * Compute the Laplacian of all components at point #p#
				      * and store them in #values#.
				      */
    virtual void vector_laplacian (const Point<dim>   &p,
				     Vector<double>     &values) const;
    
				     /**
				      * Compute the Laplacian of one component at a set of points.
				      */
    virtual void laplacian_list (const vector<Point<dim> > &points,
				 vector<double>            &values,
				 const unsigned int         component = 0) const;

				     /**
				      * Compute the Laplacians of all components at a set of points.
				      */
    virtual void vector_laplacian_list (const vector<Point<dim> > &points,
					vector<Vector<double> >   &values) const;

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
				     /**
				      * Exception
				      */
    DeclException2 (ExcWrongComponent,
		    int, int,
		    << "Component " << arg1 << " does not exist, as the values "
		    << "of this function have only " << arg2 << " components.");
};


/**
 * Provide a function which always returns zero. Obviously, also the
 * derivates of this function are zero. Also, it returns zero on all
 * components in case the function is not a scalar one, which can be
 * obtained by passing the constructor the appropriate number of
 * components.
 *
 * This function is of use when you want to implement homogeneous boundary
 * conditions, or zero initial conditions.
 *
 * @author Wolfgang Bangerth, 1998, 1999
 */
template <int dim>
class ZeroFunction : public Function<dim> {
  public:
				     /**
				      * Constructor. The number of
				      * components is preset to one.
				      */
    ZeroFunction (const unsigned int n_components = 1);
    
				     /**
				      * Virtual destructor; absolutely
				      * necessary in this case.
				      */
    virtual ~ZeroFunction ();
    
				     /**
				      * Return the value of the function
				      * at the given point for one
				      * component.
				      */
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component) const;

    				     /**
				      * Return the value of the function
				      * at the given point for all
				      * components.
				      */
    virtual void   vector_value (const Point<dim> &p,
				 Vector<double>   &return_value) const;

				     /**
				      * Set #values# to the point values
				      * of the function at the #points#,
				      * for the given component.
				      * It is assumed that #values#
				      * already has the right size, i.e.
				      * the same size as the #points#
				      * array.
				      */
    virtual void value_list (const vector<Point<dim> > &points,
			     vector<double>            &values,
			     const unsigned int         component = 0) const;

				     /**
				      * Set #values# to the point values
				      * of the function at the #points#,
				      * for all components.
				      * It is assumed that #values#
				      * already has the right size, i.e.
				      * the same size as the #points#
				      * array.
				      */
    virtual void vector_value_list (const vector<Point<dim> > &points,
				    vector<Vector<double> >   &values) const;
    
				     /**
				      * Return the gradient of the function
				      * at the given point, for the
				      * given component.
				      */
    virtual Tensor<1,dim> gradient (const Point<dim> &p,
				    const unsigned int component = 0) const;

				     /**
				      * Return the gradient of the
				      * specified component of the
				      * function at the given point,
				      * for all components.
				      */
    virtual void          vector_gradient (const Point<dim>       &p,
					   vector<Tensor<1,dim> > &gradients) const;
    
				     /**
				      * Set #gradients# to the gradients of
				      * the function at the #points#,
				      * for the given component.
				      * It is assumed that #values#
				      * already has the right size, i.e.
				      * the same size as the #points# array.
				      */
    virtual void gradient_list (const vector<Point<dim> > &points,
				vector<Tensor<1,dim> >    &gradients,
				const unsigned int         component = 0) const;
    
				     /**
				      * Set #gradients# to the gradients of
				      * the function at the #points#,
				      * for all components.
				      * It is assumed that #gradients# 
				      * already has the right size, i.e.
				      * the same size as the #points# array.
				      *
				      * The outer loop over
				      * #gradients# is over the points
				      * in the list, the inner loop
				      * over the different components
				      * of the function.
				      */
    virtual void vector_gradient_list (const vector<Point<dim> >       &points,
				       vector<vector<Tensor<1,dim> > > &gradients) const;
};


/**
 * Provide a function which always returns a constant value, which is
 * delivered upon construction. Obviously, the derivates of this
 * function are zero, which is why we derive this class from
 * #ZeroFunction#: we then only have to overload th value functions,
 * not all the derivatives. In some way, it would be more obvious to
 * do the derivation in the opposite direction, i.e. let
 * #ZeroFunction# be a more specialized version of #ConstantFunction#;
 * however, this would be more inefficient, since we could not make
 * use of the fact that the function value of the #ZeroFunction# is
 * known at compile time and need not be looked up somewhere in
 * memory.
 *
 * You can pass to the constructor an integer denoting the number of
 * components this function shall have. It defaults to one. If it is
 * greater than one, then the function will return the constant value
 * in all its components, which might not be overly useful a feature
 * in most cases, however.
 *
 * @author Wolfgang Bangerth, 1998, 1999
 */
template <int dim>
class ConstantFunction : public ZeroFunction<dim> {
  public:
				     /**
				      * Constructor; takes the constant function
				      * value as an argument. The number of
				      * components is preset to one.
				      */
    ConstantFunction (const double       value,
		      const unsigned int n_components = 1);
    
				     /**
				      * Virtual destructor; absolutely
				      * necessary in this case.
				      */
    virtual ~ConstantFunction ();

				     /**
				      * Return the value of the function
				      * at the given point for one
				      * component.
				      */
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component) const;

    				     /**
				      * Return the value of the function
				      * at the given point for all
				      * components.
				      */
    virtual void   vector_value (const Point<dim> &p,
				 Vector<double>   &return_value) const;

				     /**
				      * Set #values# to the point values
				      * of the function at the #points#,
				      * for the given component.
				      * It is assumed that #values#
				      * already has the right size, i.e.
				      * the same size as the #points#
				      * array.
				      */
    virtual void value_list (const vector<Point<dim> > &points,
			     vector<double>            &values,
			     const unsigned int         component = 0) const;

				     /**
				      * Set #values# to the point values
				      * of the function at the #points#,
				      * for all components.
				      * It is assumed that #values#
				      * already has the right size, i.e.
				      * the same size as the #points#
				      * array.
				      */
    virtual void vector_value_list (const vector<Point<dim> > &points,
				    vector<Vector<double> >   &values) const;

  protected:
				     /**
				      * Store the constant function value.
				      */
    const double function_value;
};


#endif
