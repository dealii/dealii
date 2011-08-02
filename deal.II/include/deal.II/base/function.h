//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__function_h
#define __deal2__function_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/function_time.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/point.h>
#include <vector>

DEAL_II_NAMESPACE_OPEN
template <typename number> class Vector;

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
 * @code
 *                 // access to one component at one point
 *   double value        (const Point<dim>   &p,
 *                        const unsigned int  component = 0) const;
 *
 *                 // return all components at one point
 *   void   vector_value (const Point<dim>   &p,
 *                        Vector<double>     &value) const;
 * @endcode
 *
 * For more efficiency, there are other functions returning one or all
 * components at a list of points at once:
 * @code
 *                 // access to one component at several points
 *   void   value_list (const std::vector<Point<dim> >  &point_list,
 *                      std::vector<double>             &value_list,
 *                      const unsigned int  component = 0) const;
 *
 *                 // return all components at several points
 *   void   vector_value_list (const std::vector<Point<dim> >    &point_list,
 *                             std::vector<Vector<double> >      &value_list) const;
 * @endcode
 *
 * Furthermore, there are functions returning the gradient of the
 * function at one or several points.
 *
 * You will usually only overload those functions you need; the
 * functions returning several values at a time (value_list(),
 * vector_value_list(), and gradient analoga) will call those
 * returning only one value (value(), vector_value(), and gradient
 * analoga), while those ones will throw an exception when called but
 * not overloaded.
 *
 * Note however, that the functions returning all components of the
 * function at one or several points (i.e. vector_value(),
 * vector_value_list()), will not call the function returning one
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
 * class FunctionTime.
 *
 * @note if the functions you are dealing with have sizes which
 * are a priori known (for example, <tt>dim</tt> elements), you might
 * consider using the TensorFunction class instead.
 *
 * @ingroup functions
 * @author Wolfgang Bangerth, 1998, 1999
 */
template <int dim>
class Function : public FunctionTime,
		 public Subscriptor
{
  public:
				     /**
				      * Export the value of the
				      * template parameter as a static
				      * member constant. Sometimes
				      * useful for some expression
				      * template programming.
				      */
    static const unsigned int dimension = dim;

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
				      *
				      * This destructor is declared
				      * pure virtual, such that
				      * objects of this class cannot
				      * be created. Since all the
				      * other virtual functions have a
				      * pseudo-implementation to avoid
				      * overhead in derived classes,
				      * this is the best place to do
				      * this.
				      *
				      * Nevertheless, since derived
				      * classes want to call the
				      * destructor of a base class,
				      * the destructor is implemented
				      * (despite it being pure
				      * virtual).
				      *
				      * Note: Compaq's cxx compiler
				      * does not allow defining a
				      * function that was declared
				      * pure. It simply refuses to
				      * instantiate the function when
				      * it later sees it, but also
				      * does not generate a respective
				      * function itself, which then
				      * leads to linker errors
				      * claiming that the destructor
				      * of this class is missing. We
				      * therefore only make the
				      * function abstract if the
				      * compiler can handle this.
				      */
    virtual ~Function ()
#ifndef DEAL_II_IMPLEMENTED_PURE_FUNCTION_BUG
      = 0
#endif
    ;

                                     /**
                                      * Assignment operator. This is
                                      * here only so that you can have
                                      * objects of derived classes in
                                      * containers, or assign them
                                      * otherwise. It will raise an
                                      * exception if the object from
                                      * which you assign has a
                                      * different number of components
                                      * than the one being assigned
                                      * to.
                                      */
    Function & operator= (const Function &f);

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
				      * <tt>values</tt> shall have the right
				      * size beforehand,
				      * i.e. #n_components.
				      */
    virtual void vector_value (const Point<dim>   &p,
                               Vector<double>     &values) const;

				     /**
				      * Set <tt>values</tt> to the point
				      * values of the specified
				      * component of the function at
				      * the <tt>points</tt>.  It is assumed
				      * that <tt>values</tt> already has the
				      * right size, i.e.  the same
				      * size as the <tt>points</tt> array.
				      *
				      * Be default, this function
				      * repeatedly calls value() for
				      * each point separately, to
				      * fill the output array.
				      */
    virtual void value_list (const std::vector<Point<dim> > &points,
			     std::vector<double>            &values,
			     const unsigned int              component = 0) const;

				     /**
				      * Set <tt>values</tt> to the point
				      * values of the function at the
				      * <tt>points</tt>.  It is assumed that
				      * <tt>values</tt> already has the right
				      * size, i.e.  the same size as
				      * the <tt>points</tt> array, and that
				      * all elements be vectors with
				      * the same number of components
				      * as this function has.
				      *
				      * Be default, this function
				      * repeatedly calls vector_value() for
				      * each point separately, to
				      * fill the output array.
				      */
    virtual void vector_value_list (const std::vector<Point<dim> > &points,
				    std::vector<Vector<double> >   &values) const;

				     /**
				      * For each component of the
				      * function, fill a vector of
				      * values, one for each point.
				      *
				      * The default implementation of
				      * this function in Function calls
				      * value_list() for each
				      * component. In order to improve
				      * performance, this can be
				      * reimplemented in derived
				      * classes to speed up performance.
				      */
    virtual void vector_values (const std::vector<Point<dim> >& points,
				std::vector<std::vector<double> >& values) const;

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
    virtual void vector_gradient (const Point<dim>            &p,
                                  std::vector<Tensor<1,dim> > &gradients) const;

				     /**
				      * Set <tt>gradients</tt> to the
				      * gradients of the specified
				      * component of the function at
				      * the <tt>points</tt>.  It is assumed
				      * that <tt>gradients</tt> already has the
				      * right size, i.e.  the same
				      * size as the <tt>points</tt> array.
				      */
    virtual void gradient_list (const std::vector<Point<dim> > &points,
				std::vector<Tensor<1,dim> >    &gradients,
				const unsigned int              component = 0) const;

				     /**
				      * For each component of the
				      * function, fill a vector of
				      * gradient values, one for each point.
				      *
				      * The default implementation of
				      * this function in Function calls
				      * value_list() for each
				      * component. In order to improve
				      * performance, this can be
				      * reimplemented in derived
				      * classes to speed up performance.
				      */
    virtual void vector_gradients (const std::vector<Point<dim> >            &points,
				   std::vector<std::vector<Tensor<1,dim> > > &gradients) const;

				     /**
				      * Set <tt>gradients</tt> to the gradients of
				      * the function at the <tt>points</tt>,
				      * for all components.
				      * It is assumed that <tt>gradients</tt>
				      * already has the right size, i.e.
				      * the same size as the <tt>points</tt> array.
				      *
				      * The outer loop over
				      * <tt>gradients</tt> is over the points
				      * in the list, the inner loop
				      * over the different components
				      * of the function.
				      */
    virtual void vector_gradient_list (const std::vector<Point<dim> >            &points,
				       std::vector<std::vector<Tensor<1,dim> > > &gradients) const;

				     /**
				      * Compute the Laplacian of a
				      * given component at point <tt>p</tt>.
				      */
    virtual double laplacian (const Point<dim>   &p,
			      const unsigned int  component = 0) const;

				     /**
				      * Compute the Laplacian of all
				      * components at point <tt>p</tt> and
				      * store them in <tt>values</tt>.
				      */
    virtual void vector_laplacian (const Point<dim>   &p,
				   Vector<double>     &values) const;

				     /**
				      * Compute the Laplacian of one
				      * component at a set of points.
				      */
    virtual void laplacian_list (const std::vector<Point<dim> > &points,
				 std::vector<double>            &values,
				 const unsigned int              component = 0) const;

				     /**
				      * Compute the Laplacians of all
				      * components at a set of points.
				      */
    virtual void vector_laplacian_list (const std::vector<Point<dim> > &points,
					std::vector<Vector<double> >   &values) const;

				     /**
				      * Determine an estimate for
				      * the memory consumption (in
				      * bytes) of this
				      * object. Since sometimes
				      * the size of objects can
				      * not be determined exactly
				      * (for example: what is the
				      * memory consumption of an
				      * STL <tt>std::map</tt> type with a
				      * certain number of
				      * elements?), this is only
				      * an estimate. however often
				      * quite close to the true
				      * value.
				      */
    std::size_t memory_consumption () const;
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
 * @ingroup functions
 * @author Wolfgang Bangerth, 1998, 1999
 */
template <int dim>
class ZeroFunction : public Function<dim>
{
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
    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &return_value) const;

				     /**
				      * Set <tt>values</tt> to the point values
				      * of the function at the <tt>points</tt>,
				      * for the given component.
				      * It is assumed that <tt>values</tt>
				      * already has the right size, i.e.
				      * the same size as the <tt>points</tt>
				      * array.
				      */
    virtual void value_list (const std::vector<Point<dim> > &points,
			     std::vector<double>            &values,
			     const unsigned int              component = 0) const;

				     /**
				      * Set <tt>values</tt> to the point values
				      * of the function at the <tt>points</tt>,
				      * for all components.
				      * It is assumed that <tt>values</tt>
				      * already has the right size, i.e.
				      * the same size as the <tt>points</tt>
				      * array.
				      */
    virtual void vector_value_list (const std::vector<Point<dim> > &points,
				    std::vector<Vector<double> >   &values) const;

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
    virtual void vector_gradient (const Point<dim>            &p,
                                  std::vector<Tensor<1,dim> > &gradients) const;

				     /**
				      * Set <tt>gradients</tt> to the gradients of
				      * the function at the <tt>points</tt>,
				      * for the given component.
				      * It is assumed that <tt>values</tt>
				      * already has the right size, i.e.
				      * the same size as the <tt>points</tt> array.
				      */
    virtual void gradient_list (const std::vector<Point<dim> > &points,
				std::vector<Tensor<1,dim> >    &gradients,
				const unsigned int              component = 0) const;

				     /**
				      * Set <tt>gradients</tt> to the gradients of
				      * the function at the <tt>points</tt>,
				      * for all components.
				      * It is assumed that <tt>gradients</tt>
				      * already has the right size, i.e.
				      * the same size as the <tt>points</tt> array.
				      *
				      * The outer loop over
				      * <tt>gradients</tt> is over the points
				      * in the list, the inner loop
				      * over the different components
				      * of the function.
				      */
    virtual void vector_gradient_list (const std::vector<Point<dim> >            &points,
				       std::vector<std::vector<Tensor<1,dim> > > &gradients) const;
};


/**
 * Provide a function which always returns the constant value
 * handed to the constructor.
 *
 * Obviously, the derivates of this
 * function are zero, which is why we derive this class from
 * <tt>ZeroFunction</tt>: we then only have to overload the value functions,
 * not all the derivatives. In some way, it would be more obvious to
 * do the derivation in the opposite direction, i.e. let
 * <tt>ZeroFunction</tt> be a more specialized version of <tt>ConstantFunction</tt>;
 * however, this would be less efficient, since we could not make
 * use of the fact that the function value of the <tt>ZeroFunction</tt> is
 * known at compile time and need not be looked up somewhere in
 * memory.
 *
 * You can pass to the constructor an integer denoting the number of
 * components this function shall have. It defaults to one. If it is
 * greater than one, then the function will return the constant value
 * in all its components, which might not be overly useful a feature
 * in most cases, however.
 *
 * @ingroup functions
 * @author Wolfgang Bangerth, 1998, 1999
 */
template <int dim>
class ConstantFunction : public ZeroFunction<dim>
{
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
				      * Set <tt>values</tt> to the point values
				      * of the function at the <tt>points</tt>,
				      * for the given component.
				      * It is assumed that <tt>values</tt>
				      * already has the right size, i.e.
				      * the same size as the <tt>points</tt>
				      * array.
				      */
    virtual void value_list (const std::vector<Point<dim> > &points,
			     std::vector<double>            &values,
			     const unsigned int              component = 0) const;

				     /**
				      * Set <tt>values</tt> to the point values
				      * of the function at the <tt>points</tt>,
				      * for all components.
				      * It is assumed that <tt>values</tt>
				      * already has the right size, i.e.
				      * the same size as the <tt>points</tt>
				      * array.
				      */
    virtual void vector_value_list (const std::vector<Point<dim> > &points,
				    std::vector<Vector<double> >   &values) const;

				     /**
				      * Determine an estimate for
				      * the memory consumption (in
				      * bytes) of this
				      * object. Since sometimes
				      * the size of objects can
				      * not be determined exactly
				      * (for example: what is the
				      * memory consumption of an
				      * STL <tt>std::map</tt> type with a
				      * certain number of
				      * elements?), this is only
				      * an estimate. however often
				      * quite close to the true
				      * value.
				      */
    std::size_t memory_consumption () const;

  protected:
				     /**
				      * Store the constant function value.
				      */
    const double function_value;
};



/**
 * This is a constant vector-valued function, in which one or more
 * components of the vector have a constant value and all other
 * components are zero.  It is especially useful as a weight function
 * for <tt>VectorTools::integrate_difference</tt>, where it allows to
 * integrate only one or a few vector components, rather than the
 * entire vector-valued solution. See the step-20
 * tutorial program for a detailed explanation.
 *
 * @ingroup functions
 * @author Guido Kanschat, 2000, Wolfgang Bangerth 2006
 */
template <int dim>
class ComponentSelectFunction : public ConstantFunction<dim>
{
  public:
				     /**
				      * Constructor if only a single
				      * component shall be
				      * non-zero. Arguments denote the
				      * component selected, the value
				      * for that component and the
				      * total number of vector
				      * components.
				      */
    ComponentSelectFunction (const unsigned int selected,
			     const double       value,
			     const unsigned int n_components);

				     /**
				      * Constructor. As before, but
				      * the value for the selected
				      * component is assumed to be
				      * one. In essence, this function
				      * then works as a mask.
				      */
    ComponentSelectFunction (const unsigned int selected,
			     const unsigned int n_components);

                                     /**
				      * Constructor if multiple
				      * components shall have
				      * non-zero, unit values
				      * (i.e. this should be a mask
				      * for multiple components). The
				      * first argument denotes a
				      * half-open interval of
				      * components (for example
				      * std::pair(0,dim) for the first
				      * dim components), and the
				      * second argument is the total
				      * number of vector components.
				      */
    ComponentSelectFunction (const std::pair<unsigned int, unsigned int> &selected,
			     const unsigned int n_components);

    				     /**
				      * Return the value of the function
				      * at the given point for all
				      * components.
				      */
    virtual void   vector_value (const Point<dim> &p,
				 Vector<double>   &return_value) const;

				     /**
				      * Set <tt>values</tt> to the point values
				      * of the function at the <tt>points</tt>,
				      * for all components.
				      * It is assumed that <tt>values</tt>
				      * already has the right size, i.e.
				      * the same size as the <tt>points</tt>
				      * array.
				      */
    virtual void vector_value_list (const std::vector<Point<dim> > &points,
				    std::vector<Vector<double> >   &values) const;

				     /**
				      * Determine an estimate for
				      * the memory consumption (in
				      * bytes) of this
				      * object. Since sometimes
				      * the size of objects can
				      * not be determined exactly
				      * (for example: what is the
				      * memory consumption of an
				      * STL <tt>std::map</tt> type with a
				      * certain number of
				      * elements?), this is only
				      * an estimate. however often
				      * quite close to the true
				      * value.
				      */
    std::size_t memory_consumption () const;

  protected:
				     /**
				      * Half-open interval of the
				      * indices of selected
				      * components.
				      */
    const std::pair<unsigned int,unsigned int> selected_components;
};


DEAL_II_NAMESPACE_CLOSE

#endif
