// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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

#ifndef __deal2__function_h
#define __deal2__function_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/function_time.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/point.h>
#include <deal.II/base/std_cxx11/function.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN


template <typename number> class Vector;
template <int rank, int dim, typename Number> class TensorFunction;

/**
 * This class is a model for a general function that, given a point at which
 * to evaluate the function, returns a vector of values with one or more
 * components.
 *
 * The class serves the purpose
 * of representing both scalar and vector valued functions. To this end, we
 * consider scalar functions as a special case of vector valued
 * functions, in the former case only having a single component return
 * vector. Since handling vectors is comparatively expensive, the interface
 * of this class has functions which only ask for a single component of the
 * vector-valued results (this is what you will usually need in case you know that
 * your function is scalar-valued) as well as functions you can ask for an
 * entire vector of results with as many components as the function object
 * represents. Access to function objects therefore is through the following
 * methods:
 * @code
 *                 // access to one component at one point
 *   double value        (const Point<dim, Number>   &p,
 *                        const unsigned int  component = 0) const;
 *
 *                 // return all components at one point
 *   void   vector_value (const Point<dim, Number>   &p,
 *                        Vector<double>     &value) const;
 * @endcode
 *
 * For more efficiency, there are other functions returning one or all
 * components at a list of points at once:
 * @code
 *                 // access to one component at several points
 *   void   value_list (const std::vector<Point<dim, Number> >  &point_list,
 *                      std::vector<double>             &value_list,
 *                      const unsigned int  component = 0) const;
 *
 *                 // return all components at several points
 *   void   vector_value_list (const std::vector<Point<dim, Number> >    &point_list,
 *                             std::vector<Vector<double> >      &value_list) const;
 * @endcode
 *
 * Furthermore, there are functions returning the gradient of the
 * function or even higher derivatives at one or several points.
 *
 * You will usually only overload those functions you need; the
 * functions returning several values at a time (value_list(),
 * vector_value_list(), and gradient analogs) will call those
 * returning only one value (value(), vector_value(), and gradient
 * analogs), while those ones will throw an exception when called but
 * not overloaded.
 *
 * Note however, that conversely the functions returning all components of the
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
 * Support for time dependent functions can be found in the base
 * class FunctionTime.
 *
 * @note If the functions you are dealing with have sizes which
 * are a priori known (for example, <tt>dim</tt> elements), you might
 * consider using the TensorFunction class instead. On the other hand,
 * functions like VectorTools::interpolate or
 * VectorTools::interpolate_boundary_values definitely only want objects
 * of the current type. You can use the VectorFunctionFromTensorFunction
 * class to convert the former to the latter.
 *
 * @ingroup functions
 * @author Wolfgang Bangerth, 1998, 1999, Luca Heltai 2014
 */
template <int dim, typename Number=double>
class Function : public FunctionTime<Number>,
  public Subscriptor
{
public:
  /**
   * Export the value of the template parameter as a static member
   * constant. Sometimes useful for some expression template programming.
   */
  static const unsigned int dimension = dim;

  /**
   * Number of vector components.
   */
  const unsigned int n_components;

  /**
   * Constructor. May take an initial value for the number of components
   * (which defaults to one, i.e. a scalar function), and the time
   * variable, which defaults to zero.
   */
  Function (const unsigned int n_components = 1,
            const Number       initial_time = 0.0);

  /**
   * Virtual destructor; absolutely necessary in this case.
   *
   * This destructor is declared pure virtual, such that objects of this class
   * cannot be created. Since all the other virtual functions have a
   * pseudo-implementation to avoid overhead in derived classes, they can not
   * be abstract. As a consequence, we could generate an object of this class
   * because none of this class's functions are abstract.
   *
   * We circumvent this problem by making the destructor of this class
   * abstract virtual. This ensures that at least one member function is
   * abstract, and consequently, no objects of type Function can be
   * created. However, there is no need for derived classes to explicitly
   * implement a destructor: every class has a destructor, either explicitly
   * implemented or implicitly generated by the compiler, and this resolves
   * the abstractness of any derived class even if they do not have an
   * explicitly declared destructor.
   *
   * Nonetheless, since derived classes want to call the destructor of a
   * base class, this destructor is implemented (despite it being pure
   * virtual).
   */
  virtual ~Function () = 0;

  /**
   * Assignment operator. This is here only so that you can have objects of
   * derived classes in containers, or assign them otherwise. It will raise
   * an exception if the object from which you assign has a different
   * number of components than the one being assigned to.
   */
  Function &operator= (const Function &f);

  /**
   * Return the value of the function at the given point. Unless there is
   * only one component (i.e. the function is scalar), you should state the
   * component you want to have evaluated; it defaults to zero, i.e. the
   * first component.
   */
  virtual Number value (const Point<dim, Number>   &p,
                        const unsigned int  component = 0) const;

  /**
   * Return all components of a vector-valued function at a given point.
   *
   * <tt>values</tt> shall have the right size beforehand, i.e.
   * #n_components.
   *
   * The default implementation will call value() for each component.
   */
  virtual void vector_value (const Point<dim, Number>   &p,
                             Vector<Number>     &values) const;

  /**
   * Set <tt>values</tt> to the point values of the specified component of
   * the function at the <tt>points</tt>.  It is assumed that
   * <tt>values</tt> already has the right size, i.e.  the same size as the
   * <tt>points</tt> array.
   *
   * By default, this function repeatedly calls value() for each point
   * separately, to fill the output array.
   */
  virtual void value_list (const std::vector<Point<dim, Number> > &points,
                           std::vector<Number>            &values,
                           const unsigned int              component = 0) const;

  /**
   * Set <tt>values</tt> to the point values of the function at the
   * <tt>points</tt>.  It is assumed that <tt>values</tt> already has the
   * right size, i.e.  the same size as the <tt>points</tt> array, and that
   * all elements be vectors with the same number of components as this
   * function has.
   *
   * By default, this function repeatedly calls vector_value() for each
   * point separately, to fill the output array.
   */
  virtual void vector_value_list (const std::vector<Point<dim, Number> > &points,
                                  std::vector<Vector<Number> >   &values) const;

  /**
   * For each component of the function, fill a vector of values, one for
   * each point.
   *
   * The default implementation of this function in Function calls
   * value_list() for each component. In order to improve performance, this
   * can be reimplemented in derived classes to speed up performance.
   */
  virtual void vector_values (const std::vector<Point<dim, Number> > &points,
                              std::vector<std::vector<Number> > &values) const;

  /**
   * Return the gradient of the specified component of the function at the
   * given point.
   */
  virtual Tensor<1,dim, Number> gradient (const Point<dim, Number>   &p,
                                          const unsigned int  component = 0) const;

  /**
   * Return the gradient of all components of the function at the given
   * point.
   */
  virtual void vector_gradient (const Point<dim, Number>            &p,
                                std::vector<Tensor<1,dim, Number> > &gradients) const;

  /**
   * Set <tt>gradients</tt> to the gradients of the specified component of
   * the function at the <tt>points</tt>.  It is assumed that
   * <tt>gradients</tt> already has the right size, i.e.  the same size as
   * the <tt>points</tt> array.
   */
  virtual void gradient_list (const std::vector<Point<dim, Number> > &points,
                              std::vector<Tensor<1,dim, Number> >    &gradients,
                              const unsigned int              component = 0) const;

  /**
   * For each component of the function, fill a vector of gradient values,
   * one for each point.
   *
   * The default implementation of this function in Function calls
   * value_list() for each component. In order to improve performance, this
   * can be reimplemented in derived classes to speed up performance.
   */
  virtual void vector_gradients (const std::vector<Point<dim, Number> >            &points,
                                 std::vector<std::vector<Tensor<1,dim, Number> > > &gradients) const;

  /**
   * Set <tt>gradients</tt> to the gradients of the function at the
   * <tt>points</tt>, for all components. It is assumed that
   * <tt>gradients</tt> already has the right size, i.e. the same size as
   * the <tt>points</tt> array.
   *
   * The outer loop over <tt>gradients</tt> is over the points in the list,
   * the inner loop over the different components of the function.
   */
  virtual void vector_gradient_list (const std::vector<Point<dim, Number> >            &points,
                                     std::vector<std::vector<Tensor<1,dim, Number> > > &gradients) const;

  /**
   * Compute the Laplacian of a given component at point <tt>p</tt>.
   */
  virtual Number laplacian (const Point<dim, Number>   &p,
                            const unsigned int  component = 0) const;

  /**
   * Compute the Laplacian of all components at point <tt>p</tt> and store
   * them in <tt>values</tt>.
   */
  virtual void vector_laplacian (const Point<dim, Number>   &p,
                                 Vector<Number>     &values) const;

  /**
   * Compute the Laplacian of one component at a set of points.
   */
  virtual void laplacian_list (const std::vector<Point<dim, Number> > &points,
                               std::vector<Number>            &values,
                               const unsigned int              component = 0) const;

  /**
   * Compute the Laplacians of all components at a set of points.
   */
  virtual void vector_laplacian_list (const std::vector<Point<dim, Number> > &points,
                                      std::vector<Vector<Number> >   &values) const;

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object. Since sometimes the size of objects can not be determined
   * exactly (for example: what is the memory consumption of an STL
   * <tt>std::map</tt> type with a certain number of elements?), this is
   * only an estimate. however often quite close to the true value.
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
template <int dim, typename Number=double>
class ZeroFunction : public Function<dim, Number>
{
public:
  /**
   * Constructor. The number of components is preset to one.
   */
  ZeroFunction (const unsigned int n_components = 1);

  /**
   * Virtual destructor; absolutely necessary in this case.
   *
   */
  virtual ~ZeroFunction ();

  virtual Number value (const Point<dim, Number>   &p,
                        const unsigned int  component) const;

  virtual void vector_value (const Point<dim, Number> &p,
                             Vector<Number>   &return_value) const;

  virtual void value_list (const std::vector<Point<dim, Number> > &points,
                           std::vector<Number>            &values,
                           const unsigned int              component = 0) const;

  virtual void vector_value_list (const std::vector<Point<dim, Number> > &points,
                                  std::vector<Vector<Number> >   &values) const;

  virtual Tensor<1,dim, Number> gradient (const Point<dim, Number> &p,
                                          const unsigned int component = 0) const;

  virtual void vector_gradient (const Point<dim, Number>            &p,
                                std::vector<Tensor<1,dim, Number> > &gradients) const;

  virtual void gradient_list (const std::vector<Point<dim, Number> > &points,
                              std::vector<Tensor<1,dim, Number> >    &gradients,
                              const unsigned int              component = 0) const;

  virtual void vector_gradient_list (const std::vector<Point<dim, Number> >            &points,
                                     std::vector<std::vector<Tensor<1,dim, Number> > > &gradients) const;
};



/**
 * Provide a function which always returns the constant value
 * handed to the constructor.
 *
 * Obviously, the derivates of this function are zero, which is why we
 * derive this class from <tt>ZeroFunction</tt>: we then only have to
 * overload the value functions, not all the derivatives. In some way, it
 * would be more obvious to do the derivation in the opposite direction,
 * i.e. let <tt>ZeroFunction</tt> be a more specialized version of
 * <tt>ConstantFunction</tt>; however, this would be less efficient, since
 * we could not make use of the fact that the function value of the
 * <tt>ZeroFunction</tt> is known at compile time and need not be looked up
 * somewhere in memory.
 *
 * You can pass to the constructor an integer denoting the number of
 * components this function shall have. It defaults to one. If it is
 * greater than one, then the function will return the constant value in
 * all its components, which might not be overly useful a feature in most
 * cases, however.
 *
 * @ingroup functions
 * @author Wolfgang Bangerth, 1998, 1999
 */
template <int dim, typename Number=double>
class ConstantFunction : public ZeroFunction<dim, Number>
{
public:
  /**
   * Constructor; takes the constant function value as an argument. The
   * number of components is preset to one.
   */
  ConstantFunction (const Number       value,
                    const unsigned int n_components = 1);

  /**
   * Virtual destructor; absolutely necessary in this case.
   */
  virtual ~ConstantFunction ();

  virtual Number value (const Point<dim, Number>   &p,
                        const unsigned int  component) const;

  virtual void   vector_value (const Point<dim, Number> &p,
                               Vector<Number>   &return_value) const;

  virtual void value_list (const std::vector<Point<dim, Number> > &points,
                           std::vector<Number>            &values,
                           const unsigned int              component = 0) const;

  virtual void vector_value_list (const std::vector<Point<dim, Number> > &points,
                                  std::vector<Vector<Number> >   &values) const;

  std::size_t memory_consumption () const;

protected:
  /**
   * Store the constant function value.
   */
  const Number function_value;
};



/**
 * This is a constant vector-valued function, in which one or more
 * components of the vector have a constant value and all other components
 * are zero.  It is especially useful as a weight function for
 * VectorTools::integrate_difference, where it allows to integrate only one
 * or a few vector components, rather than the entire vector-valued
 * solution. In other words, it acts as a component mask with a single
 * component selected (see the @ref GlossComponentMask "the glossary entry
 * on component masks"). See the step-20 tutorial program for a detailed
 * explanation and a use case.
 *
 * @ingroup functions
 * @author Guido Kanschat, 2000, Wolfgang Bangerth 2006
 */
template <int dim, typename Number=double>
class ComponentSelectFunction : public ConstantFunction<dim, Number>
{
public:
  /**
   * Constructor if only a single component shall be non-zero. Arguments
   * denote the component selected, the value for that component and the
   * total number of vector components.
   */
  ComponentSelectFunction (const unsigned int selected,
                           const Number       value,
                           const unsigned int n_components);

  /**
   * Constructor. As before, but the value for the selected component is
   * assumed to be one. In essence, this function then works as a mask.
   */
  ComponentSelectFunction (const unsigned int selected,
                           const unsigned int n_components);

  /**
   * Constructor if multiple components shall have non-zero, unit values
   * (i.e. this should be a mask for multiple components). The first
   * argument denotes a half-open interval of components (for example
   * std::pair(0,dim) for the first dim components), and the second
   * argument is the total number of vector components.
   */
  ComponentSelectFunction (const std::pair<unsigned int, unsigned int> &selected,
                           const unsigned int n_components);

  /**
   * Return the value of the function at the given point for all
   * components.
   */
  virtual void   vector_value (const Point<dim, Number> &p,
                               Vector<Number>   &return_value) const;

  /**
   * Set <tt>values</tt> to the point values of the function at the
   * <tt>points</tt>, for all components. It is assumed that
   * <tt>values</tt> already has the right size, i.e. the same size as the
   * <tt>points</tt> array.
   */
  virtual void vector_value_list (const std::vector<Point<dim, Number> > &points,
                                  std::vector<Vector<Number> >   &values) const;

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object. Since sometimes the size of objects can not be determined
   * exactly (for example: what is the memory consumption of an STL
   * <tt>std::map</tt> type with a certain number of elements?), this is
   * only an estimate. however often quite close to the true value.
   */
  std::size_t memory_consumption () const;

protected:
  /**
   * Half-open interval of the indices of selected components.
   */
  const std::pair<unsigned int,unsigned int> selected_components;
};



/**
 * This class provides a way to convert a scalar function of the kind
 * @code
 *   Number foo (const Point<dim, Number> &);
 * @endcode
 * into an object of type Function@<dim@>. Since the argument returns a
 * scalar, the result is clearly a Function object for which
 * <code>function.n_components==1</code>. The class works by storing a
 * pointer to the given function and every time
 * <code>function.value(p,component)</code> is called, calls
 * <code>foo(p)</code> and returns the corresponding value. It also makes
 * sure that <code>component</code> is in fact zero, as needs be for scalar
 * functions.
 *
 * The class provides an easy way to turn a simple global function into
 * something that has the required Function@<dim@> interface for operations
 * like VectorTools::interpolate_boundary_values() etc., and thereby allows
 * for simpler experimenting without having to write all the boiler plate
 * code of declaring a class that is derived from Function and implementing
 * the Function::value() function. An example of this is given in the results
 * section of step-53.
 *
 * The class gains additional expressive power because the argument it
 * takes does not have to be a pointer to an actual function. Rather, it is
 * a function object, i.e., it can also be the result of call to std::bind
 * (or boost::bind) or some other object that can be called with a single
 * argument. For example, if you need a Function object that returns the
 * norm of a point, you could write it like so:
 * @code
 *   template <int dim, typename Number>
 *   class Norm : public Function<dim, Number> {
 *     public:
 *       virtual Number value (const Point<dim, Number> &p,
 *                             const unsigned int component) const {
 *         Assert (component == 0, ExcMessage ("This object is scalar!"));
 *         return p.norm();
 *       }
 *    };
 *
 *    Norm<2> my_norm_object;
 * @endcode
 * and then pass the <code>my_norm_object</code> around, or you could write
 * it like so:
 * @code
 *   ScalarFunctionFromFunctionObject<dim, Number> my_norm_object (&Point<dim, Number>::norm);
 * @endcode
 *
 * Similarly, to generate an object that computes the distance to a point
 * <code>q</code>, we could do this:
 * @code
 *   template <int dim, typename Number>
 *   class DistanceTo : public Function<dim, Number> {
 *     public:
 *       DistanceTo (const Point<dim, Number> &q) : q(q) {}
 *       virtual Number value (const Point<dim, Number> &p,
 *                             const unsigned int component) const {
 *         Assert (component == 0, ExcMessage ("This object is scalar!"));
 *         return q.distance(p);
 *       }
 *     private:
 *       const Point<dim, Number> q;
 *    };
 *
 *    Point<2> q (2,3);
 *    DistanceTo<2> my_distance_object;
 * @endcode
 * or we could write it like so:
 * @code
 *    ScalarFunctionFromFunctionObject<dim, Number>
 *      my_distance_object (std_cxx11::bind (&Point<dim, Number>::distance,
 *                                           q,
 *                                           std_cxx11::_1));
 * @endcode
 * The savings in work to write this are apparent.
 *
 * @author Wolfgang Bangerth, 2011
 */
template <int dim, typename Number=double>
class ScalarFunctionFromFunctionObject : public Function<dim, Number>
{
public:
  /**
   * Given a function object that takes a Point and returns a Number value,
   * convert this into an object that matches the Function<dim, Number> interface.
   */
  ScalarFunctionFromFunctionObject (const std_cxx11::function<Number (const Point<dim, Number> &)> &function_object);

  /**
   * Return the value of the function at the given point. Returns the value
   * the function given to the constructor produces for this point.
   */
  virtual Number value (const Point<dim, Number>   &p,
                        const unsigned int  component = 0) const;

private:
  /**
   * The function object which we call when this class's value() or
   * value_list() functions are called.
   **/
  const std_cxx11::function<Number (const Point<dim, Number> &)> function_object;
};



/**
 * This class is similar to the ScalarFunctionFromFunctionObject class in
 * that it allows for the easy conversion of a function object to something
 * that satisfies the interface of the Function base class. The difference
 * is that here, the given function object is still a scalar function (i.e.
 * it has a single value at each space point) but that the Function object
 * generated is vector valued. The number of vector components is specified
 * in the constructor, where one also selectes a single one of these vector
 * components that should be filled by the passed object. The result is a
 * vector Function object that returns zero in each component except the
 * single selected one where it returns the value returned by the given as
 * the first argument to the constructor.
 *
 * @note In the above discussion, note the difference between the (scalar)
 * "function object" (i.e., a C++ object <code>x</code> that can be called
 * as in <code>x(p)</code>) and the capitalized (vector valued) "Function
 * object" (i.e., an object of a class that is derived from the Function
 * base class).
 *
 * To be more concrete, let us consider the following example:
 * @code
 *   Number one (const Point<2> &p) { return 1; }
 *   VectorFunctionFromScalarFunctionObject<2>
 *      component_mask (&one, 1, 3);
 * @endcode
 * Here, <code>component_mask</code> then represents a Function object that
 * for every point returns the vector $(0, 1, 0)^T$, i.e. a mask function
 * that could, for example, be passed to
 * VectorTools::integrate_difference(). This effect can also be achieved
 * using the ComponentSelectFunction class but is obviously easily extended
 * to functions that are non-constant in their one component.
 *
 * @author Wolfgang Bangerth, 2011
 */
template <int dim, typename Number=double>
class VectorFunctionFromScalarFunctionObject : public Function<dim, Number>
{
public:
  /**
   * Given a function object that takes a Point and returns a Number
   * value, convert this into an object that matches the Function@<dim@>
   * interface.
   *
   * @param function_object The scalar function that will form one component
   *     of the resulting Function object.
   * @param n_components The total number of vector components of the
   *     resulting Function object.
   * @param selected_component The single component that should be
   *     filled by the first argument.
   **/
  VectorFunctionFromScalarFunctionObject (const std_cxx11::function<Number (const Point<dim, Number> &)> &function_object,
                                          const unsigned int selected_component,
                                          const unsigned int n_components);

  /**
   * Return the value of the function at the given point. Returns the value
   * the function given to the constructor produces for this point.
   */
  virtual Number value (const Point<dim, Number>   &p,
                        const unsigned int  component = 0) const;

  /**
   * Return all components of a vector-valued function at a given point.
   *
   * <tt>values</tt> shall have the right size beforehand, i.e.
   * #n_components.
   */
  virtual void vector_value (const Point<dim, Number>   &p,
                             Vector<Number>     &values) const;

private:
  /**
   * The function object which we call when this class's value() or
   * value_list() functions are called.
   **/
  const std_cxx11::function<Number (const Point<dim, Number> &)> function_object;

  /**
   * The vector component whose value is to be filled by the given scalar
   * function.
   */
  const unsigned int selected_component;
};


/**
 * This class is built as a means of translating the <code>Tensor<1,dim, Number>
 * </code> values produced by objects of type TensorFunction and returning
 * them as a multiple component version of the same thing as a Vector for
 * use in, for example, the VectorTools::interpolate or the many other
 * functions taking Function objects. It allows the user to place the
 * desired components into an <tt>n_components</tt> long vector starting at
 * the <tt>selected_component</tt> location in that vector and have all
 * other components be 0.
 *
 *  For example: Say you created a class called
 *  @code
 *    class RightHandSide : public TensorFunction<rank,dim, Number>
 *  @endcode
 *  which extends the TensorFunction class and you have an object
 *  @code
 *    RightHandSide<1,dim, Number> rhs;
 *  @endcode
 *  of that class which you want to interpolate onto your mesh using the
 *  VectorTools::interpolate function, but the finite element you use for
 *  the DoFHandler object has 3 copies of a finite element with
 *  <tt>dim</tt> components, for a total of 3*dim components. To
 *  interpolate onto that DoFHandler, you need an object of type Function
 *  that has 3*dim vector components. Creating such an object from the
 *  existing <code>rhs</code> object is done using this piece of code:
 *  @code
 *      RighHandSide<1,dim, Number> rhs;
 *      VectorFunctionFromTensorFunction<dim, Number> rhs_vector_function (rhs, 0, 3*dim);
 *  @endcode
 *  where the <code>dim</code> components of the tensor function are placed into
 *  the first <code>dim</code> components of the function object.
 *
 * @author Spencer Patty, 2013
 */
template <int dim, typename Number=double>
class VectorFunctionFromTensorFunction : public Function<dim, Number>
{
public:
  /**
   * Given a TensorFunction object that takes a <tt>Point</tt> and returns
   * a <tt>Tensor<1,dim, Number></tt> value, convert this into an object that
   * matches the Function@<dim@> interface.
   *
   * By default, create a Vector object of the same size as
   * <tt>tensor_function</tt> returns, i.e., with <tt>dim</tt> components.
   *
   * @param tensor_function The TensorFunction that will form one component
   *     of the resulting Vector Function object.
   * @param n_components The total number of vector components of the
   *     resulting TensorFunction object.
   * @param selected_component The first component that should be
   *     filled by the first argument.  This should be such that the entire
   *     tensor_function fits inside the <tt>n_component</tt> length return
   *     vector.
   */
  VectorFunctionFromTensorFunction (const TensorFunction<1,dim, Number> &tensor_function,
                                    const unsigned int selected_component=0,
                                    const unsigned int n_components=dim);

  /**
   * This destructor is defined as virtual so as to coincide with all other
   * aspects of class.
   */
  virtual ~VectorFunctionFromTensorFunction();

  /**
   * Return a single component of a vector-valued function at a given
   * point.
   */
  virtual Number value (const Point<dim, Number> &p,
                        const unsigned int component = 0) const;

  /**
   * Return all components of a vector-valued function at a given point.
   *
   * <tt>values</tt> shall have the right size beforehand, i.e.
   * #n_components.
   */
  virtual void vector_value (const Point<dim, Number> &p,
                             Vector<Number>   &values) const;

  /**
   * Return all components of a vector-valued function at a list of points.
   *
   * <tt>value_list</tt> shall be the same size as <tt>points</tt> and each
   * element of the vector will be passed to vector_value() to evaluate the
   * function
   */
  virtual void vector_value_list (const std::vector<Point<dim, Number> > &points,
                                  std::vector<Vector<Number> >   &value_list) const;

private:
  /**
   * The TensorFunction object which we call when this class's
   * vector_value() or vector_value_list() functions are called.
   **/
  const TensorFunction<1,dim,Number> &tensor_function;

  /**
   * The first vector component whose value is to be filled by the given
   * TensorFunction.  The values will be placed in components
   * selected_component to selected_component+dim-1 for a
   * <tt>TensorFunction<1,dim, Number></tt> object.
   */
  const unsigned int selected_component;
};


DEAL_II_NAMESPACE_CLOSE

#endif
