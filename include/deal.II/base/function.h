// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_function_h
#define dealii_function_h


#include <deal.II/base/config.h>

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/function_time.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/point.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <functional>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <typename number>
class Vector;
template <int rank, int dim, typename Number>
class TensorFunction;
#endif

/**
 * This class is a model for a general function that, given a point at which
 * to evaluate the function, returns a vector of values with one or more
 * components.
 *
 * The class serves the purpose of representing both scalar and vector valued
 * functions. To this end, we consider scalar functions as a special case of
 * vector valued functions, in the former case only having a single component
 * return vector. Since handling vectors is comparatively expensive, the
 * interface of this class has functions which only ask for a single component
 * of the vector-valued results (this is what you will usually need in case
 * you know that your function is scalar-valued) as well as functions you can
 * ask for an entire vector of results with as many components as the function
 * object represents. Access to function objects therefore is through the
 * following methods:
 * @code
 * // access to one component at one point
 * double
 * value(const Point<dim> & p,
 *       const unsigned int component = 0) const;
 *
 * // return all components at one point
 * void
 * vector_value(const Point<dim> &p,
 *              Vector<double>   &value) const;
 * @endcode
 *
 * For more efficiency, there are other functions returning one or all
 * components at a list of points at once:
 * @code
 * // access to one component at several points
 * void
 * value_list(const std::vector<Point<dim>> &point_list,
 *            std::vector<double>           &value_list,
 *            const unsigned int             component = 0) const;
 *
 * // return all components at several points
 * void
 * vector_value_list(const std::vector<Point<dim>> &point_list,
 *                   std::vector<Vector<double>>   &value_list) const;
 * @endcode
 *
 * Furthermore, there are functions returning the gradient of the function or
 * even higher derivatives at one or several points.
 *
 * You will usually only overload those functions you need; the functions
 * returning several values at a time (value_list(), vector_value_list(), and
 * gradient analogs) will call those returning only one value (value(),
 * vector_value(), and gradient analogs), while those ones will throw an
 * exception when called but not overloaded.
 *
 * Conversely, the functions returning all components of the function at one
 * or several points (i.e. vector_value(), vector_value_list()), will
 * <em>not</em> call the function returning one component at one point
 * repeatedly, once for each point and component. The reason is efficiency:
 * this would amount to too many virtual function calls. If you have
 * vector-valued functions, you should therefore also provide overloads of the
 * virtual functions for all components at a time.
 *
 * Also note, that unless only called a very small number of times, you should
 * overload all sets of functions (returning only one value, as well as those
 * returning a whole array), since the cost of evaluation of a point value is
 * often less than the virtual function call itself.
 *
 * Support for time dependent functions can be found in the base class
 * FunctionTime.
 *
 *
 * <h3>Functions that return tensors</h3>
 *
 * If the functions you are dealing with have a number of components that are
 * a priori known (for example, <tt>dim</tt> elements), you might consider
 * using the TensorFunction class instead. This is, in particular, true if the
 * objects you return have the properties of a tensor, i.e., they are for
 * example dim-dimensional vectors or dim-by-dim matrices. On the other hand,
 * functions like VectorTools::interpolate or
 * VectorTools::interpolate_boundary_values definitely only want objects of
 * the current type. You can use the VectorFunctionFromTensorFunction class to
 * convert the former to the latter.
 *
 *
 * <h3>Functions that return vectors of other data types</h3>
 *
 * Most of the time, your functions will have the form $f : \Omega \rightarrow
 * {\mathbb R}^{n_\text{components}}$. However, there are occasions where you
 * want the function to return vectors (or scalars) over a different number
 * field, for example functions that return complex numbers or vectors of
 * complex numbers: $f : \Omega \rightarrow {\mathbb
 * C}^{n_\text{components}}$. In such cases, you can choose a value different
 * than the default @p double for the second template
 * argument of this class: it describes the scalar type to be used for each
 * component of your return values. It defaults to @p double, but in the
 * example above, it could be set to <code>std::complex@<double@></code>.
 * step-58 is an example of this.
 *
 * @tparam dim The space dimension of the range space within which the domain
 *   $\Omega$ of the function lies. Consequently, the function will be
 *   evaluated at objects of type @p Point<dim>.
 * @tparam RangeNumberType The scalar type of the vector space that is
 *   the range (or image) of this function. As discussed above,
 *   objects of the current type represent functions from ${\mathbb
 *   R}^\text{dim}$ to $S^{n_\text{components}}$ where $S$ is the
 *   underlying scalar type of the vector space. The type of $S$ is
 *   given by the @p RangeNumberType template argument.
 *
 * @ingroup functions
 */
template <int dim, typename RangeNumberType = double>
class Function : public FunctionTime<
                   typename numbers::NumberTraits<RangeNumberType>::real_type>,
                 public EnableObserverPointer
{
public:
  /**
   * Export the value of the template parameter as a static member constant.
   * Sometimes useful for some expression template programming.
   */
  static constexpr unsigned int dimension = dim;

  /**
   * Number of vector components.
   */
  const unsigned int n_components;

  /**
   * The scalar-valued real type used for representing time.
   */
  using time_type = typename FunctionTime<
    typename numbers::NumberTraits<RangeNumberType>::real_type>::time_type;

  /**
   * Constructor. May take an initial value for the number of components
   * (which defaults to one, i.e. a scalar function), and the time variable,
   * which defaults to zero.
   */
  explicit Function(const unsigned int n_components = 1,
                    const time_type    initial_time = 0.0);

  /**
   * Copy constructor.
   */
  Function(const Function &f) = default;

  /**
   * Virtual destructor; absolutely necessary in this case.
   *
   * This destructor is declared pure virtual, such that objects of this class
   * cannot be created. Since all the other virtual functions have a
   * pseudo-implementation to avoid overhead in derived classes, they can not be
   * abstract. As a consequence, we could generate an object of this class
   * because none of this class's functions are abstract.
   *
   * We circumvent this problem by making the destructor of this class
   * abstract virtual. This ensures that at least one member function is
   * abstract, and consequently, no objects of type Function can be created.
   * However, there is no need for derived classes to explicitly implement a
   * destructor: every class has a destructor, either explicitly implemented
   * or implicitly generated by the compiler, and this resolves the
   * abstractness of any derived class even if they do not have an explicitly
   * declared destructor.
   *
   * Nonetheless, since derived classes want to call the destructor of a base
   * class, this destructor is implemented (despite it being pure virtual).
   */
  virtual ~Function() override = 0;

  /**
   * Assignment operator. This is here only so that you can have objects of
   * derived classes in containers, or assign them otherwise. It will raise an
   * exception if the object from which you assign has a different number of
   * components than the one being assigned to.
   */
  Function &
  operator=(const Function &f);

  /**
   * Return the value of the function at the given point. Unless there is only
   * one component (i.e. the function is scalar), you should state the
   * component you want to have evaluated; it defaults to zero, i.e. the first
   * component.
   */
  virtual RangeNumberType
  value(const Point<dim> &p, const unsigned int component = 0) const;

  /**
   * Return all components of a vector-valued function at a given point.
   *
   * <tt>values</tt> shall have the right size beforehand, i.e. #n_components.
   *
   * The default implementation will call value() for each component.
   */
  virtual void
  vector_value(const Point<dim> &p, Vector<RangeNumberType> &values) const;

  /**
   * Set <tt>values</tt> to the point values of the specified component of the
   * function at the <tt>points</tt>.  It is assumed that <tt>values</tt>
   * already has the right size, i.e.  the same size as the <tt>points</tt>
   * array.
   *
   * By default, this function repeatedly calls value() for each point
   * separately, to fill the output array.
   */
  virtual void
  value_list(const std::vector<Point<dim>> &points,
             std::vector<RangeNumberType>  &values,
             const unsigned int             component = 0) const;

  /**
   * Set <tt>values</tt> to the point values of the function at the
   * <tt>points</tt>.  It is assumed that <tt>values</tt> already has the
   * right size, i.e.  the same size as the <tt>points</tt> array, and that
   * all elements be vectors with the same number of components as this
   * function has.
   *
   * By default, this function repeatedly calls vector_value() for each point
   * separately, to fill the output array.
   */
  virtual void
  vector_value_list(const std::vector<Point<dim>>        &points,
                    std::vector<Vector<RangeNumberType>> &values) const;

  /**
   * For each component of the function, fill a vector of values, one for each
   * point.
   *
   * The default implementation of this function in Function calls
   * value_list() for each component. In order to improve performance, this
   * can be reimplemented in derived classes to speed up performance.
   */
  virtual void
  vector_values(const std::vector<Point<dim>>             &points,
                std::vector<std::vector<RangeNumberType>> &values) const;

  /**
   * Return the gradient of the specified component of the function at the
   * given point.
   */
  virtual Tensor<1, dim, RangeNumberType>
  gradient(const Point<dim> &p, const unsigned int component = 0) const;

  /**
   * Return the gradient of all components of the function at the given point.
   */
  virtual void
  vector_gradient(
    const Point<dim>                             &p,
    std::vector<Tensor<1, dim, RangeNumberType>> &gradients) const;

  /**
   * Set <tt>gradients</tt> to the gradients of the specified component of the
   * function at the <tt>points</tt>.  It is assumed that <tt>gradients</tt>
   * already has the right size, i.e.  the same size as the <tt>points</tt>
   * array.
   */
  virtual void
  gradient_list(const std::vector<Point<dim>>                &points,
                std::vector<Tensor<1, dim, RangeNumberType>> &gradients,
                const unsigned int component = 0) const;

  /**
   * For each component of the function, fill a vector of gradient values, one
   * for each point.
   *
   * The default implementation of this function in Function calls
   * value_list() for each component. In order to improve performance, this
   * can be reimplemented in derived classes to speed up performance.
   */
  virtual void
  vector_gradients(
    const std::vector<Point<dim>>                             &points,
    std::vector<std::vector<Tensor<1, dim, RangeNumberType>>> &gradients) const;

  /**
   * Set <tt>gradients</tt> to the gradients of the function at the
   * <tt>points</tt>, for all components. It is assumed that
   * <tt>gradients</tt> already has the right size, i.e. the same size as the
   * <tt>points</tt> array.
   *
   * The outer loop over <tt>gradients</tt> is over the points in the list,
   * the inner loop over the different components of the function.
   */
  virtual void
  vector_gradient_list(
    const std::vector<Point<dim>>                             &points,
    std::vector<std::vector<Tensor<1, dim, RangeNumberType>>> &gradients) const;

  /**
   * Compute the Laplacian of a given component at point <tt>p</tt>.
   */
  virtual RangeNumberType
  laplacian(const Point<dim> &p, const unsigned int component = 0) const;

  /**
   * Compute the Laplacian of all components at point <tt>p</tt> and store
   * them in <tt>values</tt>.
   */
  virtual void
  vector_laplacian(const Point<dim> &p, Vector<RangeNumberType> &values) const;

  /**
   * Compute the Laplacian of one component at a set of points.
   */
  virtual void
  laplacian_list(const std::vector<Point<dim>> &points,
                 std::vector<RangeNumberType>  &values,
                 const unsigned int             component = 0) const;

  /**
   * Compute the Laplacians of all components at a set of points.
   */
  virtual void
  vector_laplacian_list(const std::vector<Point<dim>>        &points,
                        std::vector<Vector<RangeNumberType>> &values) const;

  /**
   * Compute the Hessian of a given component at point <tt>p</tt>, that is the
   * gradient of the gradient of the function.
   */
  virtual SymmetricTensor<2, dim, RangeNumberType>
  hessian(const Point<dim> &p, const unsigned int component = 0) const;

  /**
   * Compute the Hessian of all components at point <tt>p</tt> and store them
   * in <tt>values</tt>.
   */
  virtual void
  vector_hessian(
    const Point<dim>                                      &p,
    std::vector<SymmetricTensor<2, dim, RangeNumberType>> &values) const;

  /**
   * Compute the Hessian of one component at a set of points.
   */
  virtual void
  hessian_list(const std::vector<Point<dim>>                         &points,
               std::vector<SymmetricTensor<2, dim, RangeNumberType>> &values,
               const unsigned int component = 0) const;

  /**
   * Compute the Hessians of all components at a set of points.
   */
  virtual void
  vector_hessian_list(
    const std::vector<Point<dim>>                                      &points,
    std::vector<std::vector<SymmetricTensor<2, dim, RangeNumberType>>> &values)
    const;


  /**
   * Return an estimate for the memory consumption, in bytes, of this object.
   *
   * This function is virtual and can be overloaded by derived classes.
   */
  virtual std::size_t
  memory_consumption() const;
};


namespace Functions
{
  /**
   * Provide a function which always returns the constant values handed to the
   * constructor.
   *
   * @ingroup functions
   */
  template <int dim, typename RangeNumberType = double>
  class ConstantFunction : public Function<dim, RangeNumberType>
  {
  public:
    /**
     * Constructor; set values of all components to the provided one. The
     * default number of components is one.
     */
    explicit ConstantFunction(const RangeNumberType value,
                              const unsigned int    n_components = 1);

    /**
     * Constructor; takes an <tt>std::vector<RangeNumberType></tt> object as an
     * argument. The number of components is determined by
     * <tt>values.size()</tt>.
     */
    explicit ConstantFunction(const std::vector<RangeNumberType> &values);

    /**
     * Constructor; takes an <tt>Vector<RangeNumberType></tt> object as an
     * argument. The number of components is determined by
     * <tt>values.size()</tt>.
     */
    explicit ConstantFunction(const Vector<RangeNumberType> &values);

    /**
     * Constructor; uses whatever stores in [begin_ptr, begin_ptr+n_components)
     * to initialize a new object.
     */
    ConstantFunction(const RangeNumberType *begin_ptr,
                     const unsigned int     n_components);

    virtual RangeNumberType
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    virtual void
    vector_value(const Point<dim>        &p,
                 Vector<RangeNumberType> &return_value) const override;

    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<RangeNumberType>  &return_values,
               const unsigned int             component = 0) const override;

    virtual void
    vector_value_list(
      const std::vector<Point<dim>>        &points,
      std::vector<Vector<RangeNumberType>> &return_values) const override;

    virtual Tensor<1, dim, RangeNumberType>
    gradient(const Point<dim>  &p,
             const unsigned int component = 0) const override;

    virtual void
    vector_gradient(
      const Point<dim>                             &p,
      std::vector<Tensor<1, dim, RangeNumberType>> &gradients) const override;

    virtual void
    gradient_list(const std::vector<Point<dim>>                &points,
                  std::vector<Tensor<1, dim, RangeNumberType>> &gradients,
                  const unsigned int component = 0) const override;

    virtual void
    vector_gradient_list(
      const std::vector<Point<dim>>                             &points,
      std::vector<std::vector<Tensor<1, dim, RangeNumberType>>> &gradients)
      const override;

    virtual SymmetricTensor<2, dim, RangeNumberType>
    hessian(const Point<dim>  &point,
            const unsigned int component = 0) const override;

    virtual RangeNumberType
    laplacian(const Point<dim>  &point,
              const unsigned int component = 0) const override;

    virtual std::size_t
    memory_consumption() const override;

  protected:
    /**
     * Store the constant function value vector.
     */
    std::vector<RangeNumberType> function_value_vector;
  };



  /**
   * Provide a function which always returns zero. Obviously, also the
   * derivatives of this function are zero. Also, it returns zero on all
   * components in case the function is not a scalar one, which can be obtained
   * by passing the constructor the appropriate number of components.
   *
   * This function is of use when you want to implement homogeneous boundary
   * conditions, or zero initial conditions.
   *
   * @ingroup functions
   */
  template <int dim, typename RangeNumberType = double>
  class ZeroFunction : public ConstantFunction<dim, RangeNumberType>
  {
  public:
    /**
     * Constructor. The number of components is preset to one.
     */
    explicit ZeroFunction(const unsigned int n_components = 1);
  };

  /**
   * A function whose output is also its input. One possible application of this
   * function is interpolating or projecting a finite element field that
   * represents spatial coordinates: e.g., one can set up a finite element field
   * to interpolate the positions of a Triangulation's cells with this function
   * (via VectorTools::interpolate()), which is useful when doing calculations
   * in a Lagrangian reference frame.
   *
   * @ingroup functions
   */
  template <int dim, typename RangeNumberType = double>
  class IdentityFunction : public Function<dim, RangeNumberType>
  {
  public:
    /**
     * Constructor. The number of components is set to dim.
     */
    IdentityFunction();

    /**
     * @copydoc Function::value()
     */
    virtual RangeNumberType
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * @copydoc Function::gradient()
     */
    virtual Tensor<1, dim, RangeNumberType>
    gradient(const Point<dim>  &p,
             const unsigned int component = 0) const override;

    /**
     * @copydoc Function::laplacian()
     */
    virtual RangeNumberType
    laplacian(const Point<dim>  &p,
              const unsigned int component = 0) const override;

    /**
     * @copydoc Function::hessian()
     */
    virtual SymmetricTensor<2, dim, RangeNumberType>
    hessian(const Point<dim>  &p,
            const unsigned int component = 0) const override;
  };
} // namespace Functions

/**
 * This is a constant vector-valued function, in which one or more components
 * of the vector have a constant value and all other components are zero.  It
 * is especially useful as a weight function for
 * VectorTools::integrate_difference, where it allows to integrate only one or
 * a few vector components, rather than the entire vector-valued solution. In
 * other words, it acts as a component mask with a single component selected
 * (see the
 * @ref GlossComponentMask "the glossary entry on component masks").
 * See the step-20 tutorial program for a detailed explanation and a use case.
 *
 * @ingroup functions
 */
template <int dim, typename RangeNumberType = double>
class ComponentSelectFunction
  : public Functions::ConstantFunction<dim, RangeNumberType>
{
public:
  /**
   * Constructor if only a single component shall be non-zero. Arguments
   * denote the component selected, the value for that component and the total
   * number of vector components.
   */
  ComponentSelectFunction(const unsigned int    selected,
                          const RangeNumberType value,
                          const unsigned int    n_components);

  /**
   * Constructor. As before, but the value for the selected component is
   * assumed to be one. In essence, this function then works as a mask.
   */
  ComponentSelectFunction(const unsigned int selected,
                          const unsigned int n_components);

  /**
   * Constructor if multiple components shall have non-zero, unit values (i.e.
   * this should be a mask for multiple components). The first argument
   * denotes a half-open interval of components (for example std::pair(0,dim)
   * for the first dim components), and the second argument is the total
   * number of vector components.
   */
  ComponentSelectFunction(const std::pair<unsigned int, unsigned int> &selected,
                          const unsigned int n_components);


  /**
   * Substitute function value with value of a <tt>ConstantFunction@<dim,
   * RangeNumberType@></tt> object and keep the current selection pattern.
   *
   * This is useful if you want to have different values in different
   * components since the provided constructors of
   * <tt>ComponentSelectFunction@<dim, RangeNumberType@></tt> class can only
   * have same value for all components.
   *
   * @note We copy the underlying component value data from @p f from its
   * beginning. So the number of components of @p f cannot be less than the
   * calling object.
   */
  virtual void
  substitute_function_value_with(
    const Functions::ConstantFunction<dim, RangeNumberType> &f);

  /**
   * Return the value of the function at the given point for all components.
   */
  virtual void
  vector_value(const Point<dim>        &p,
               Vector<RangeNumberType> &return_value) const override;

  /**
   * Set <tt>values</tt> to the point values of the function at the
   * <tt>points</tt>, for all components. It is assumed that <tt>values</tt>
   * already has the right size, i.e. the same size as the <tt>points</tt>
   * array.
   */
  virtual void
  vector_value_list(
    const std::vector<Point<dim>>        &points,
    std::vector<Vector<RangeNumberType>> &values) const override;

  /**
   * Return an estimate for the memory consumption, in bytes, of this object.
   */
  virtual std::size_t
  memory_consumption() const override;

protected:
  /**
   * Half-open interval of the indices of selected components.
   */
  const std::pair<unsigned int, unsigned int> selected_components;
};



/**
 * This class provides a way to convert a scalar function of the kind
 * @code
 * RangeNumberType foo (const Point<dim> &);
 * @endcode
 * into an object of type Function@<dim@>. Since the argument returns a
 * scalar, the result is clearly a Function object for which
 * <code>function.n_components == 1</code>. The class works by storing a pointer
 * to the given function and every time
 * <code>function.value(p,component)</code> is called, calls
 * <code>foo(p)</code> and returns the corresponding value. It also makes sure
 * that <code>component</code> is in fact zero, as needs be for scalar
 * functions.
 *
 * The class provides an easy way to turn a simple global function into
 * something that has the required Function@<dim@> interface for operations
 * like VectorTools::interpolate_boundary_values() etc., and thereby allows
 * for simpler experimenting without having to write all the boiler plate code
 * of declaring a class that is derived from Function and implementing the
 * Function::value() function. An example of this is given in the results
 * section of step-53.
 *
 * The class gains additional expressive power because the argument it takes
 * does not have to be a pointer to an actual function. Rather, it is a
 * function object, i.e., it can also be the result of a lambda function or some
 * other object that can be called with a single argument. For
 * example, if you need a Function object that returns the norm of a point, you
 * could write it like so:
 * @code
 * template <int dim, typename RangeNumberType>
 * class Norm : public Function<dim, RangeNumberType>
 * {
 * public:
 *   virtual RangeNumberType
 *   value(const Point<dim> & p,
 *         const unsigned int component) const
 *   {
 *     Assert (component == 0, ExcMessage ("This object is scalar!"));
 *     return p.norm();
 *   }
 * };
 *
 * Norm<2> my_norm_object;
 * @endcode
 * and then pass the <code>my_norm_object</code> around, or you could write it
 * like so:
 * @code
 * ScalarFunctionFromFunctionObject<dim, RangeNumberType> my_norm_object(
 *   &Point<dim>::norm);
 * @endcode
 *
 * Similarly, to generate an object that computes the distance to a point
 * <code>q</code>, we could do this:
 * @code
 * template <int dim, typename RangeNumberType>
 * class DistanceTo : public Function<dim, RangeNumberType>
 * {
 * public:
 *   DistanceTo (const Point<dim> &q) : q(q) {}
 *
 *   virtual RangeNumberType
 *   value (const Point<dim> & p,
 *          const unsigned int component) const
 *   {
 *     Assert(component == 0, ExcMessage("This object is scalar!"));
 *     return q.distance(p);
 *   }
 * private:
 *   const Point<dim> q;
 * };
 *
 * Point<2> q (2, 3);
 * DistanceTo<2> my_distance_object;
 * @endcode
 * or we could write it like so:
 * @code
 * ScalarFunctionFromFunctionObject<dim, RangeNumberType> my_distance_object(
 *   [&q](const Point<dim> &p){return q.distance(p);});
 * @endcode
 * The savings in work to write this are apparent.
 *
 * Finally, these lambda functions can be used as a way to map points in
 * different ways. As an example, let us assume that we have computed
 * the solution to a one-dimensional problem and that that solution
 * resides in the following variables:
 * @code
 *   DoFHandler<1>  dof_handler_1d;
 *   Vector<double> solution_1d;
 * @endcode
 * We will denote this solution function described by this DoFHandler
 * and vector object by $u_h(x)$ where $x$ is a vector with just one
 * component, and consequently is not shown in boldface. Then assume
 * that we want this $u_h(x)$ to be used as a boundary condition for a 2d
 * problem at the line $y=0$. Let's say that this line corresponds to
 * @ref GlossBoundaryIndicator "boundary indicator"
 * 123.
 * If we say that the 2d problem is associated with
 * @code
 *   DoFHandler<2> dof_handler_2d;
 * @endcode
 * then in order to evaluate the boundary conditions for this 2d problem,
 * we would want to call VectorTools::interpolate_boundary_values()
 * via
 * @code
 *   AffineConstraints<double> boundary_values_2d;
 *   VectorTools::interpolate_boundary_values (dof_handler_2d,
 *                                             123,
 *                                             ???,
 *                                             boundary_values_2d);
 * @endcode
 * The question here is what to use as the Function object that can be passed
 * as third argument. It needs to be a Function<2> object, i.e., it
 * receives a 2d input point and is supposed to return the value at that
 * point. What we *want* it to do is to just take the $x$ component of the
 * input point and evaluate the 1d solution at that point, knowing that at
 * the boundary with indicator 123, the $y$ component of the input point
 * must be zero. This all can be achieved via the following function
 * object:
 * @code
 *   Functions::FEFieldFunction<1>
 *     solution_1d_as_function_object (dof_handler_1d, solution_1d);
 *   auto boundary_evaluator
 *     = [&] (const Point<2> &p)
 *       {
 *          // First extract the x component of the input point:
 *          const Point<1> point_on_axis (p[0]);
 *
 *          // Then evaluate the 1d solution at that point:
 *          return solution_1d_as_function_object.value(point_on_axis);
 *       }
 *
 *   AffineConstraints<double> boundary_values_2d;
 *   VectorTools::interpolate_boundary_values (dof_handler_2d,
 *                                             123,
 *                                             ScalarFunctionFromFunctionObject<2>(boundary_evaluator),
 *                                             boundary_values_2d);
 * @endcode
 *
 * @ingroup functions
 */
template <int dim, typename RangeNumberType = double>
class ScalarFunctionFromFunctionObject : public Function<dim, RangeNumberType>
{
public:
  /**
   * Given a function object that takes a Point and returns a RangeNumberType
   * value, convert this into an object that matches the Function<dim,
   * RangeNumberType> interface.
   */
  explicit ScalarFunctionFromFunctionObject(
    const std::function<RangeNumberType(const Point<dim> &)> &function_object);

  /**
   * Given a function object that takes  a time and a Point and returns a
   * RangeNumberType value, convert this into an object that matches the
   * Function<dim, RangeNumberType> interface.
   */
  explicit ScalarFunctionFromFunctionObject(
    const std::function<RangeNumberType(const double, const Point<dim> &)>
      &function_object_t);

  /**
   * Return the value of the function at the given point. Returns the value
   * the function given to the constructor produces for this point.
   */
  virtual RangeNumberType
  value(const Point<dim> &p, const unsigned int component = 0) const override;

private:
  /**
   * The function object which we call when this class's value() or
   * value_list() functions are called.
   */
  const std::function<RangeNumberType(const double, const Point<dim> &)>
    function_object;
};



/**
 * This class is similar to the ScalarFunctionFromFunctionObject class in that
 * it allows for the easy conversion of a function object to something that
 * satisfies the interface of the Function base class. The difference is that
 * here, the given function object is still a scalar function (i.e. it has a
 * single value at each space point) but that the Function object generated is
 * vector valued. The number of vector components is specified in the
 * constructor, where one also selects a single one of these vector components
 * that should be filled by the passed object. The result is a vector Function
 * object that returns zero in each component except the single selected one
 * where it returns the value returned by the given as the first argument to
 * the constructor.
 *
 * @note In the above discussion, note the difference between the (scalar)
 * "function object" (i.e., a C++ object <code>x</code> that can be called as
 * in <code>x(p)</code>) and the capitalized (vector valued) "Function object"
 * (i.e., an object of a class that is derived from the Function base class).
 *
 * To be more concrete, let us consider the following example:
 * @code
 * RangeNumberType
 * one(const Point<2> &p)
 * {
 *   return 1.0;
 * }
 *
 * VectorFunctionFromScalarFunctionObject<2, RangeNumberType>
 *   component_mask(&one, 1, 3);
 * @endcode
 * Here, <code>component_mask</code> then represents a Function object that
 * for every point returns the vector $(0, 1, 0)^T$, i.e. a mask function that
 * could, for example, be passed to VectorTools::integrate_difference(). This
 * effect can also be achieved using the ComponentSelectFunction class but is
 * obviously easily extended to functions that are non-constant in their one
 * component.
 *
 * @ingroup functions
 */
template <int dim, typename RangeNumberType = double>
class VectorFunctionFromScalarFunctionObject
  : public Function<dim, RangeNumberType>
{
public:
  /**
   * Given a function object that takes a Point and returns a RangeNumberType
   * value, convert this into an object that matches the Function@<dim@>
   * interface.
   *
   * @param function_object The scalar function that will form one component
   * of the resulting Function object.
   * @param n_components The total number of vector components of the
   * resulting Function object.
   * @param selected_component The single component that should be filled by
   * the first argument.
   */
  VectorFunctionFromScalarFunctionObject(
    const std::function<RangeNumberType(const Point<dim> &)> &function_object,
    const unsigned int selected_component,
    const unsigned int n_components);

  /**
   * Return the value of the function at the given point. Returns the value
   * the function given to the constructor produces for this point.
   */
  virtual RangeNumberType
  value(const Point<dim> &p, const unsigned int component = 0) const override;

  /**
   * Return all components of a vector-valued function at a given point.
   *
   * <tt>values</tt> shall have the right size beforehand, i.e. #n_components.
   */
  virtual void
  vector_value(const Point<dim>        &p,
               Vector<RangeNumberType> &values) const override;

private:
  /**
   * The function object which we call when this class's value() or
   * value_list() functions are called.
   */
  const std::function<RangeNumberType(const Point<dim> &)> function_object;

  /**
   * The vector component whose value is to be filled by the given scalar
   * function.
   */
  const unsigned int selected_component;
};


/**
 * This class is similar to the ScalarFunctionFromFunctionObject and
 * VectorFunctionFromFunctionObject classes in that it allows for the easy
 * conversion of a vector of function objects to something that satisfies the
 * interface of the Function base class.
 *
 * The difference is that here the Function object generated may be vector
 * valued, and you can specify the gradients of the function. The number of
 * vector components is deduced from the size of the vector in the constructor.
 *
 * To be more concrete, let us consider the following example:
 *
 * @code
 * RangeNumberType
 * first_component(const Point<2> &p)
 * {
 *   return 1.0;
 * }
 *
 * RangeNumberType
 * second_component(const Point<2> &p)
 * {
 *   return 2.0;
 * }
 *
 * Tensor<1, 2, RangeNumberType>
 * zero_gradient(const Point<2> &) {
 *   return Tensor<1, 2, RangeNumberType>();
 * }
 *
 * FunctionFromFunctionObjects<2, RangeNumberType>
 *     custom_function({&first_component, &second_component},
 *                     {&zero_gradient, &zero_gradient});
 * @endcode
 */
template <int dim, typename RangeNumberType = double>
class FunctionFromFunctionObjects : public Function<dim, RangeNumberType>
{
public:
  /**
   * Default constructor.
   *
   * This constructor does not initialize the internal methods. To have a
   * usable function, you need to call at least the set_function_values()
   * method. If you need also the gradients of the solution, then you must
   * also call the set_function_gradients() method.
   */
  explicit FunctionFromFunctionObjects(const unsigned int n_components = 1,
                                       const double       initial_time = 0);

  /**
   * Constructor for functions of which you only know the values.
   *
   * The resulting function will have a number of components equal to the size
   * of the vector @p values. A call to the FunctionFromFunctionObject::gradient()
   * method will trigger an exception, unless you first call the
   * set_function_gradients() method.
   */
  explicit FunctionFromFunctionObjects(
    const std::vector<std::function<RangeNumberType(const Point<dim> &)>>
                &values,
    const double initial_time = 0.0);

  /**
   * Constructor for functions of which you only know the values.
   *
   * The resulting function will have a number of components equal @p n_components.
   * A call to the FunctionFromFunctionObject::gradient()
   * method will trigger an exception, unless you first call the
   * set_function_gradients() method.
   */
  explicit FunctionFromFunctionObjects(
    const std::function<RangeNumberType(const Point<dim> &, const unsigned int)>
                      &values,
    const unsigned int n_components,
    const double       initial_time = 0.0);

  /**
   * Constructor for functions of which you know both the values and the
   * gradients.
   *
   * The resulting function will have a number of components equal to the size
   * of the vector @p values. If the size of @p values and @p gradients does not
   * match, an exception is triggered.
   */
  FunctionFromFunctionObjects(
    const std::vector<std::function<RangeNumberType(const Point<dim> &)>>
      &values,
    const std::vector<
      std::function<Tensor<1, dim, RangeNumberType>(const Point<dim> &)>>
                &gradients,
    const double initial_time = 0.0);


  /**
   * Return the value of the function at the given point. Unless there is only
   * one component (i.e. the function is scalar), you should state the
   * component you want to have evaluated; it defaults to zero, i.e. the first
   * component.
   */
  virtual RangeNumberType
  value(const Point<dim> &p, const unsigned int component = 0) const override;

  /**
   * Return the gradient of the function at the given point. Unless there is
   * only one component (i.e. the function is scalar), you should state the
   * component you want to have evaluated; it defaults to zero, i.e. the first
   * component.
   */
  virtual Tensor<1, dim, RangeNumberType>
  gradient(const Point<dim>  &p,
           const unsigned int component = 0) const override;

  /**
   * Reset the function values of this object. An assertion is thrown if the
   * size of the @p values parameter does not match the number of components of
   * this object.
   */
  void
  set_function_values(
    const std::vector<std::function<RangeNumberType(const Point<dim> &)>>
      &values);

  /**
   * Reset the function gradients of this object. An assertion is thrown if the
   * size of the @p gradients parameter does not match the number of components of
   * this object.
   */
  void
  set_function_gradients(
    const std::vector<
      std::function<Tensor<1, dim, RangeNumberType>(const Point<dim> &)>>
      &gradients);

private:
  /**
   * The actual function values.
   */
  std::function<RangeNumberType(const Point<dim> &, const unsigned int)>
    function_values;

  /**
   * The actual function gradients.
   */
  std::function<Tensor<1, dim, RangeNumberType>(const Point<dim> &,
                                                const unsigned int)>
    function_gradients;
};


/**
 * This class is built as a means of translating the <code>Tensor<1,dim,
 * RangeNumberType> </code> values produced by objects of type TensorFunction
 * and returning them as a multiple component version of the same thing as a
 * Vector for use in, for example, the VectorTools::interpolate or the many
 * other functions taking Function objects. It allows the user to place the
 * desired components into an <tt>n_components</tt> long vector starting at
 * the <tt>selected_component</tt> location in that vector and have all other
 * components be 0.
 *
 * For example: Say you created a class called
 * @code
 * class RightHandSide : public TensorFunction<rank,dim, RangeNumberType>
 * @endcode
 * which extends the TensorFunction class and you have an object
 * @code
 * RightHandSide<1,dim, RangeNumberType> rhs;
 * @endcode
 * of that class which you want to interpolate onto your mesh using the
 * VectorTools::interpolate function, but the finite element you use for the
 * DoFHandler object has 3 copies of a finite element with <tt>dim</tt>
 * components, for a total of 3*dim components. To interpolate onto that
 * DoFHandler, you need an object of type Function that has 3*dim vector
 * components. Creating such an object from the existing <code>rhs</code>
 * object is done using this piece of code:
 * @code
 * RightHandSide<1,dim, RangeNumberType> rhs;
 * VectorFunctionFromTensorFunction<dim, RangeNumberType> rhs_vector_function(
 *   rhs, 0, 3*dim);
 * @endcode
 * where the <code>dim</code> components of the tensor function are placed
 * into the first <code>dim</code> components of the function object.
 *
 * @ingroup functions
 */
template <int dim, typename RangeNumberType = double>
class VectorFunctionFromTensorFunction : public Function<dim, RangeNumberType>
{
public:
  /**
   * Given a TensorFunction object that takes a <tt>Point</tt> and returns a
   * <tt>Tensor<1,dim, RangeNumberType></tt> value, convert this into an object
   * that matches the Function@<dim@> interface.
   *
   * By default, create a Vector object of the same size as
   * <tt>tensor_function</tt> returns, i.e., with <tt>dim</tt> components.
   *
   * @param tensor_function The TensorFunction that will form `dim` components of
   * the resulting Vector Function object.
   * @param n_components The total number of vector components of the
   * resulting TensorFunction object. This clearly has to be at least `dim`.
   * @param selected_component The first component that should be filled by
   * the first argument.  This should be such that the entire `tensor_function`
   * fits inside the <tt>n_component</tt> length return vector.
   */
  explicit VectorFunctionFromTensorFunction(
    const TensorFunction<1, dim, RangeNumberType> &tensor_function,
    const unsigned int                             selected_component = 0,
    const unsigned int                             n_components       = dim);

  /**
   * This destructor is defined as virtual so as to coincide with all other
   * aspects of class.
   */
  virtual ~VectorFunctionFromTensorFunction() override = default;

  /**
   * Return a single component of a vector-valued function at a given point.
   */
  virtual RangeNumberType
  value(const Point<dim> &p, const unsigned int component = 0) const override;

  /**
   * Return all components of a vector-valued function at a given point.
   *
   * <tt>values</tt> shall have the right size beforehand, i.e. #n_components.
   */
  virtual void
  vector_value(const Point<dim>        &p,
               Vector<RangeNumberType> &values) const override;

  /**
   * Return all components of a vector-valued function at a list of points.
   *
   * <tt>value_list</tt> shall be the same size as <tt>points</tt> and each
   * element of the vector will be passed to vector_value() to evaluate the
   * function
   */
  virtual void
  vector_value_list(
    const std::vector<Point<dim>>        &points,
    std::vector<Vector<RangeNumberType>> &value_list) const override;

  /**
   * Return the gradient of the specified component of the function at the given
   * point.
   */
  virtual Tensor<1, dim, RangeNumberType>
  gradient(const Point<dim>  &p,
           const unsigned int component = 0) const override;

  /**
   * Return the gradient of all components of the function at the given point.
   */
  virtual void
  vector_gradient(
    const Point<dim>                             &p,
    std::vector<Tensor<1, dim, RangeNumberType>> &gradients) const override;

  /**
   * Set <tt>gradients</tt> to the gradients of the specified component of the
   * function at the <tt>points</tt>.  It is assumed that <tt>gradients</tt>
   * already has the right size, i.e.  the same size as the <tt>points</tt>
   * array.
   */
  virtual void
  gradient_list(const std::vector<Point<dim>>                &points,
                std::vector<Tensor<1, dim, RangeNumberType>> &gradients,
                const unsigned int component = 0) const override;

  /**
   * For each component of the function, fill a vector of gradient values, one
   * for each point.
   *
   * The default implementation of this function in Function calls
   * value_list() for each component. In order to improve performance, this
   * can be reimplemented in derived classes to speed up performance.
   */
  virtual void
  vector_gradients(const std::vector<Point<dim>> &points,
                   std::vector<std::vector<Tensor<1, dim, RangeNumberType>>>
                     &gradients) const override;

  /**
   * Set <tt>gradients</tt> to the gradients of the function at the
   * <tt>points</tt>, for all components. It is assumed that
   * <tt>gradients</tt> already has the right size, i.e. the same size as the
   * <tt>points</tt> array.
   *
   * The outer loop over <tt>gradients</tt> is over the points in the list,
   * the inner loop over the different components of the function.
   */
  virtual void
  vector_gradient_list(const std::vector<Point<dim>> &points,
                       std::vector<std::vector<Tensor<1, dim, RangeNumberType>>>
                         &gradients) const override;

private:
  /**
   * The TensorFunction object which we call when this class's vector_value()
   * or vector_value_list() functions are called.
   */
  const TensorFunction<1, dim, RangeNumberType> &tensor_function;

  /**
   * The first vector component whose value is to be filled by the given
   * TensorFunction.  The values will be placed in components
   * selected_component to selected_component+dim-1 for a
   * <tt>TensorFunction<1,dim, RangeNumberType></tt> object.
   */
  const unsigned int selected_component;
};


/**
 * This class is built as a means of translating the <code>Tensor<1,dim,
 * RangeNumberType> </code> values produced by function objects that
 * for a given point return a tensor,
 * and returning them as a multiple component version of the same thing as a
 * Vector for use in, for example, the VectorTools::interpolate or the many
 * other functions taking Function objects. It allows the user to place the
 * desired components into an <tt>n_components</tt> long vector starting at
 * the <tt>selected_component</tt> location in that vector and have all other
 * components be 0.
 *
 * For example: Say you created a function object that returns the gravity
 * (here, a radially inward pointing vector of magnitude 9.81):
 * @code
 *   const auto gravity
 *     = [](const Point<dim> &p) -> Tensor<1,dim> { return -9.81 * (p /
 * p.norm()); }
 * @endcode
 * and you want to interpolate this onto your mesh using the
 * VectorTools::interpolate function, with a finite element for the
 * DoFHandler object has 3 copies of a finite element with <tt>dim</tt>
 * components, for a total of 3*dim components. To interpolate onto that
 * DoFHandler, you need an object of type Function that has 3*dim vector
 * components. Creating such an object from the existing <code>gravity</code>
 * object is done using this piece of code:
 * @code
 * VectorFunctionFromTensorFunctionObject<dim, RangeNumberType>
 *   gravity_vector_function(gravity, 0, 3*dim);
 * @endcode
 * where the <code>dim</code> components of the `gravity` function are placed
 * into the first <code>dim</code> components of the function object.
 *
 * @ingroup functions
 */
template <int dim, typename RangeNumberType = double>
class VectorFunctionFromTensorFunctionObject
  : public Function<dim, RangeNumberType>
{
public:
  /**
   * Given a function object that takes a <tt>Point</tt> and returns a
   * <tt>Tensor<1,dim, RangeNumberType></tt> value, convert this into an object
   * that matches the Function@<dim@> interface.
   *
   * By default, create a Vector object of the same size as
   * <tt>tensor_function</tt> returns, i.e., with <tt>dim</tt> components.
   *
   * @param tensor_function_object The TensorFunction that will form `dim` components of
   * the resulting Vector Function object.
   * @param n_components The total number of vector components of the
   * resulting TensorFunction object. This clearly has to be at least `dim`.
   * @param selected_component The first component that should be filled by
   * the first argument.  This should be such that the entire tensor_function
   * fits inside the <tt>n_component</tt> length return vector.
   */
  explicit VectorFunctionFromTensorFunctionObject(
    const std::function<Tensor<1, dim, RangeNumberType>(const Point<dim> &)>
                      &tensor_function_object,
    const unsigned int selected_component = 0,
    const unsigned int n_components       = dim);

  /**
   * This destructor is defined as virtual so as to coincide with all other
   * aspects of class.
   */
  virtual ~VectorFunctionFromTensorFunctionObject() override = default;

  /**
   * Return a single component of a vector-valued function at a given point.
   */
  virtual RangeNumberType
  value(const Point<dim> &p, const unsigned int component = 0) const override;

  /**
   * Return all components of a vector-valued function at a given point.
   *
   * <tt>values</tt> shall have the right size beforehand, i.e. #n_components.
   */
  virtual void
  vector_value(const Point<dim>        &p,
               Vector<RangeNumberType> &values) const override;

  /**
   * Return all components of a vector-valued function at a list of points.
   *
   * <tt>value_list</tt> shall be the same size as <tt>points</tt> and each
   * element of the vector will be passed to vector_value() to evaluate the
   * function
   */
  virtual void
  vector_value_list(
    const std::vector<Point<dim>>        &points,
    std::vector<Vector<RangeNumberType>> &value_list) const override;

private:
  /**
   * The TensorFunction object which we call when this class's vector_value()
   * or vector_value_list() functions are called.
   */
  const std::function<Tensor<1, dim, RangeNumberType>(const Point<dim> &)>
    tensor_function_object;

  /**
   * The first vector component whose value is to be filled by the given
   * TensorFunction.  The values will be placed in components
   * selected_component to selected_component+dim-1 for a
   * <tt>TensorFunction<1,dim, RangeNumberType></tt> object.
   */
  const unsigned int selected_component;
};


#ifndef DOXYGEN
// icc 2018 complains about an undefined reference
// if we put this in the templates.h file
//
// The destructor is pure virtual so we can't default it
// in the declaration.
template <int dim, typename RangeNumberType>
inline Function<dim, RangeNumberType>::~Function() = default;
#endif


DEAL_II_NAMESPACE_CLOSE

#endif
