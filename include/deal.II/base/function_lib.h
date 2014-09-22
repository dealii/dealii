// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2014 by the deal.II authors
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

#ifndef __deal2__function_lib_h
#define __deal2__function_lib_h


#include <deal.II/base/config.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/table.h>

#include <deal.II/base/std_cxx11/array.h>

DEAL_II_NAMESPACE_OPEN

/**
 * Namespace implementing some concrete classes derived from the
 * Function class that describe actual functions. This is rather
 * a collection of classes that we have needed for our own programs
 * once and thought might be useful to others as well at some
 * point.
 *
 * @ingroup functions
 */
namespace Functions
{


  /**
   * The distance to the origin squared.
   *
   * This function returns the square norm of the radius vector of a point.
   *
   * Together with the function, its derivatives and Laplacian are defined.
   *
   * @ingroup functions
   * @author: Guido Kanschat, 1999
   */
  template<int dim>
  class SquareFunction : public Function<dim>
  {
  public:
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
    virtual void vector_value (const Point<dim>   &p,
                               Vector<double>     &values) const;
    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<double>            &values,
                             const unsigned int              component = 0) const;
    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
                                    const unsigned int  component = 0) const;
    virtual void vector_gradient (const Point<dim>   &p,
                                  std::vector<Tensor<1,dim> >    &gradient) const;
    virtual void gradient_list (const std::vector<Point<dim> > &points,
                                std::vector<Tensor<1,dim> >    &gradients,
                                const unsigned int              component = 0) const;
    virtual double laplacian (const Point<dim>   &p,
                              const unsigned int  component = 0) const;
    virtual void laplacian_list (const std::vector<Point<dim> > &points,
                                 std::vector<double>            &values,
                                 const unsigned int              component = 0) const;
  };



  /**
   * The function <tt>xy</tt> in 2d and 3d, not implemented in 1d.
   * This function serves as an example for
   * a vanishing Laplacian.
   *
   * @ingroup functions
   * @author: Guido Kanschat, 2000
   */
  template<int dim>
  class Q1WedgeFunction : public Function<dim>
  {
  public:
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<double>            &values,
                             const unsigned int              component = 0) const;

    virtual void vector_value_list (const std::vector<Point<dim> > &points,
                                    std::vector<Vector<double> > &values) const;

    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
                                    const unsigned int       component = 0) const;

    virtual void gradient_list (const std::vector<Point<dim> > &points,
                                std::vector<Tensor<1,dim> >    &gradients,
                                const unsigned int              component = 0) const;

    virtual void vector_gradient_list (const std::vector<Point<dim> > &,
                                       std::vector<std::vector<Tensor<1,dim> > > &) const;

    /**
     * Laplacian of the function at one point.
     */
    virtual double laplacian (const Point<dim>   &p,
                              const unsigned int  component = 0) const;

    /**
     * Laplacian of the function at multiple points.
     */
    virtual void laplacian_list (const std::vector<Point<dim> > &points,
                                 std::vector<double>            &values,
                                 const unsigned int              component = 0) const;
  };



  /**
   * d-quadratic pillow on the unit hypercube.
   *
   * This is a function for testing the implementation. It has zero Dirichlet
   * boundary values on the domain $(-1,1)^d$. In the inside, it is the
   * product of $1-x_i^2$ over all space dimensions.
   *
   * Providing a non-zero argument to the constructor, the whole function
   * can be offset by a constant.
   *
   * Together with the function, its derivatives and Laplacian are defined.
   *
   * @ingroup functions
   * @author: Guido Kanschat, 1999
   */
  template<int dim>
  class PillowFunction : public Function<dim>
  {
  public:
    /**
     * Constructor. Provide a
     * constant that will be added to
     * each function value.
     */
    PillowFunction (const double offset=0.);

    /**
     * The value at a single point.
     */
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    /**
     * Values at multiple points.
     */
    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<double>            &values,
                             const unsigned int              component = 0) const;

    /**
     * Gradient at a single point.
     */
    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
                                    const unsigned int  component = 0) const;

    /**
     * Gradients at multiple points.
     */
    virtual void gradient_list (const std::vector<Point<dim> > &points,
                                std::vector<Tensor<1,dim> >    &gradients,
                                const unsigned int              component = 0) const;

    /**
     * Laplacian at a single point.
     */
    virtual double laplacian (const Point<dim>   &p,
                              const unsigned int  component = 0) const;

    /**
     * Laplacian at multiple points.
     */
    virtual void laplacian_list (const std::vector<Point<dim> > &points,
                                 std::vector<double>            &values,
                                 const unsigned int              component = 0) const;
  private:
    const double offset;
  };



  /**
   * Cosine-shaped pillow function.
   * This is another function with zero boundary values on $[-1,1]^d$. In the interior
   * it is the product of $\cos(\pi/2 x_i)$.
   *
   * @ingroup functions
   * @author Guido Kanschat, 1999
   */
  template<int dim>
  class CosineFunction : public Function<dim>
  {
  public:
    /**
     * Constructor which allows to
     * optionally generate a vector
     * valued cosine function with
     * the same value in each
     * component.
     */
    CosineFunction (const unsigned int n_components = 1);

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<double>            &values,
                             const unsigned int              component = 0) const;

    virtual void vector_value_list (const std::vector<Point<dim> > &points,
                                    std::vector<Vector<double> > &values) const;

    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
                                    const unsigned int  component = 0) const;

    virtual void gradient_list (const std::vector<Point<dim> > &points,
                                std::vector<Tensor<1,dim> >    &gradients,
                                const unsigned int              component = 0) const;

    virtual double laplacian (const Point<dim>   &p,
                              const unsigned int  component = 0) const;

    virtual void laplacian_list (const std::vector<Point<dim> > &points,
                                 std::vector<double>            &values,
                                 const unsigned int              component = 0) const;

    /**
     * Second derivatives at a
     * single point.
     */
    virtual Tensor<2,dim> hessian (const Point<dim>   &p,
                                   const unsigned int  component = 0) const;

    /**
     * Second derivatives at
     * multiple points.
     */
    virtual void hessian_list (const std::vector<Point<dim> > &points,
                               std::vector<Tensor<2,dim> >    &hessians,
                               const unsigned int              component = 0) const;
  };



  /**
   * Gradient of the cosine-shaped pillow function.
   *
   * This is a vector-valued function with @p dim components, the
   * gradient of CosineFunction. On the square [-1,1], it has tangential
   * boundary conditions zero. Thus, it can be used to test
   * implementations of Maxwell operators without bothering about
   * boundary terms.
   *
   * @ingroup functions
   * @author Guido Kanschat, 2010
   */
  template<int dim>
  class CosineGradFunction : public Function<dim>
  {
  public:
    /**
     * Constructor, creating a
     * function with @p dim components.
     */
    CosineGradFunction ();

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component) const;
    virtual void vector_value (const Point<dim>   &p,
                               Vector<double>     &values) const;
    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<double>            &values,
                             const unsigned int              component) const;

    virtual void vector_value_list (const std::vector<Point<dim> > &points,
                                    std::vector<Vector<double> > &values) const;

    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
                                    const unsigned int  component) const;

    virtual void gradient_list (const std::vector<Point<dim> > &points,
                                std::vector<Tensor<1,dim> >    &gradients,
                                const unsigned int              component) const;

    virtual void vector_gradient_list (const std::vector<Point<dim> >            &points,
                                       std::vector<std::vector<Tensor<1,dim> > > &gradients) const;

    virtual double laplacian (const Point<dim>   &p,
                              const unsigned int  component) const;
  };



  /**
   * Product of exponential functions in each coordinate direction.
   *
   * @ingroup functions
   * @author Guido Kanschat, 1999
   */
  template<int dim>
  class ExpFunction : public Function<dim>
  {
  public:
    /**
     * The value at a single point.
     */
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    /**
     * Values at multiple points.
     */
    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<double>            &values,
                             const unsigned int              component = 0) const;

    /**
     * Gradient at a single point.
     */
    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
                                    const unsigned int  component = 0) const;

    /**
     * Gradients at multiple points.
     */
    virtual void gradient_list (const std::vector<Point<dim> > &points,
                                std::vector<Tensor<1,dim> >    &gradients,
                                const unsigned int              component = 0) const;

    /**
     * Laplacian at a single point.
     */
    virtual double laplacian (const Point<dim>   &p,
                              const unsigned int  component = 0) const;

    /**
     * Laplacian at multiple points.
     */
    virtual void laplacian_list (const std::vector<Point<dim> > &points,
                                 std::vector<double>            &values,
                                 const unsigned int              component = 0) const;
  };



  /**
   * Harmonic singularity on the L-shaped domain in 2D.
   *
   * @ingroup functions
   * @author Guido Kanschat
   * @date 1999
   */
  class LSingularityFunction : public Function<2>
  {
  public:
    virtual double value (const Point<2>   &p,
                          const unsigned int  component = 0) const;

    virtual void value_list (const std::vector<Point<2> > &points,
                             std::vector<double>            &values,
                             const unsigned int              component = 0) const;

    virtual void vector_value_list (const std::vector<Point<2> > &points,
                                    std::vector<Vector<double> > &values) const;

    virtual Tensor<1,2> gradient (const Point<2>     &p,
                                  const unsigned int  component = 0) const;

    virtual void gradient_list (const std::vector<Point<2> > &points,
                                std::vector<Tensor<1,2> >    &gradients,
                                const unsigned int            component = 0) const;

    virtual void vector_gradient_list (const std::vector<Point<2> > &,
                                       std::vector<std::vector<Tensor<1,2> > > &) const;

    virtual double laplacian (const Point<2>   &p,
                              const unsigned int  component = 0) const;

    virtual void laplacian_list (const std::vector<Point<2> > &points,
                                 std::vector<double>          &values,
                                 const unsigned int            component = 0) const;
  };



  /**
   * Gradient of the harmonic singularity on the L-shaped domain in 2D.
   *
   * The gradient of LSingularityFunction, which is a vector valued
   * function with vanishing curl and divergence.
   *
   * @ingroup functions
   * @author Guido Kanschat, 2010
   */
  class LSingularityGradFunction : public Function<2>
  {
  public:
    /**
     * Default constructor setting
     * the dimension to 2.
     */
    LSingularityGradFunction ();
    virtual double value (const Point<2>   &p,
                          const unsigned int  component) const;

    virtual void value_list (const std::vector<Point<2> > &points,
                             std::vector<double>            &values,
                             const unsigned int              component) const;

    virtual void vector_value_list (const std::vector<Point<2> > &points,
                                    std::vector<Vector<double> > &values) const;

    virtual Tensor<1,2> gradient (const Point<2>     &p,
                                  const unsigned int  component) const;

    virtual void gradient_list (const std::vector<Point<2> > &points,
                                std::vector<Tensor<1,2> >    &gradients,
                                const unsigned int            component) const;

    virtual void vector_gradient_list (const std::vector<Point<2> > &,
                                       std::vector<std::vector<Tensor<1,2> > > &) const;

    virtual double laplacian (const Point<2>   &p,
                              const unsigned int  component) const;

    virtual void laplacian_list (const std::vector<Point<2> > &points,
                                 std::vector<double>          &values,
                                 const unsigned int            component) const;
  };



  /**
   * Singularity on the slit domain in 2D and 3D.
   *
   * @ingroup functions
   * @author Guido Kanschat, 1999, 2006
   */
  template <int dim>
  class SlitSingularityFunction : public Function<dim>
  {
  public:
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<double>            &values,
                             const unsigned int              component = 0) const;

    virtual void vector_value_list (const std::vector<Point<dim> > &points,
                                    std::vector<Vector<double> > &values) const;

    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
                                    const unsigned int  component = 0) const;

    virtual void gradient_list (const std::vector<Point<dim> > &points,
                                std::vector<Tensor<1,dim> >    &gradients,
                                const unsigned int            component = 0) const;

    virtual void vector_gradient_list (const std::vector<Point<dim> > &,
                                       std::vector<std::vector<Tensor<1,dim> > > &) const;

    virtual double laplacian (const Point<dim>   &p,
                              const unsigned int  component = 0) const;

    virtual void laplacian_list (const std::vector<Point<dim> > &points,
                                 std::vector<double>          &values,
                                 const unsigned int            component = 0) const;
  };


  /**
   * Singularity on the slit domain with one Neumann boundary in 2D.
   *
   * @ingroup functions
   * @author Guido Kanschat, 2002
   */
  class SlitHyperSingularityFunction : public Function<2>
  {
  public:
    virtual double value (const Point<2>   &p,
                          const unsigned int  component = 0) const;

    virtual void value_list (const std::vector<Point<2> > &points,
                             std::vector<double>            &values,
                             const unsigned int              component = 0) const;

    virtual void vector_value_list (const std::vector<Point<2> > &points,
                                    std::vector<Vector<double> > &values) const;

    virtual Tensor<1,2> gradient (const Point<2>   &p,
                                  const unsigned int  component = 0) const;

    virtual void gradient_list (const std::vector<Point<2> > &points,
                                std::vector<Tensor<1,2> >    &gradients,
                                const unsigned int            component = 0) const;

    virtual void vector_gradient_list (const std::vector<Point<2> > &,
                                       std::vector<std::vector<Tensor<1,2> > > &) const;

    virtual double laplacian (const Point<2>   &p,
                              const unsigned int  component = 0) const;

    virtual void laplacian_list (const std::vector<Point<2> > &points,
                                 std::vector<double>          &values,
                                 const unsigned int            component = 0) const;
  };



  /**
   * A jump in x-direction transported into some direction.
   *
   * If the advection is parallel to the y-axis, the function is
   * <tt>-atan(sx)</tt>, where <tt>s</tt> is the steepness parameter provided in
   * the constructor.
   *
   * For different advection directions, this function will be turned in
   * the parameter space.
   *
   * Together with the function, its derivatives and Laplacian are defined.
   *
   * @ingroup functions
   * @author: Guido Kanschat, 2000
   */
  template<int dim>
  class JumpFunction : public Function<dim>
  {
  public:
    /**
     * Constructor. Provide the
     * advection direction here and
     * the steepness of the slope.
     */
    JumpFunction (const Point<dim> &direction,
                  const double      steepness);

    /**
     * Function value at one point.
     */
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    /**
     * Function values at multiple
     * points.
     */
    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<double>            &values,
                             const unsigned int              component = 0) const;

    /**
     * Gradient at one point.
     */
    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
                                    const unsigned int  component = 0) const;

    /**
       Gradients at multiple points.
    */
    virtual void gradient_list (const std::vector<Point<dim> > &points,
                                std::vector<Tensor<1,dim> >    &gradients,
                                const unsigned int              component = 0) const;

    /**
     * Laplacian of the function at one point.
     */
    virtual double laplacian (const Point<dim>   &p,
                              const unsigned int  component = 0) const;

    /**
     * Laplacian of the function at multiple points.
     */
    virtual void laplacian_list (const std::vector<Point<dim> > &points,
                                 std::vector<double>            &values,
                                 const unsigned int              component = 0) const;

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
     * Advection vector.
     */
    const Point<dim> direction;

    /**
     * Steepness (maximal derivative)
     * of the slope.
     */
    const double steepness;

    /**
     * Advection angle.
     */
    double angle;

    /**
     * Sine of <tt>angle</tt>.
     */
    double sine;

    /**
     * Cosine of <tt>angle</tt>.
     */
    double cosine;
  };



  /**
   * Given a wavenumber vector generate a cosine function. The
   * wavenumber coefficient is given as a $d$-dimensional point $k$
   * in Fourier space, and the function is then recovered as $f(x) =
   * \cos(\sum_i k_i x_i) = Re(\exp(i k.x))$.
   *
   * The class has its name from the fact that it resembles one
   * component of a Fourier cosine decomposition.
   *
   * @ingroup functions
   * @author Wolfgang Bangerth, 2001
   */
  template <int dim>
  class FourierCosineFunction : public Function<dim>
  {
  public:
    /**
     * Constructor. Take the Fourier
     * coefficients in each space
     * direction as argument.
     */
    FourierCosineFunction (const Point<dim> &fourier_coefficients);

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
     * Return the gradient of the
     * specified component of the
     * function at the given point.
     */
    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
                                    const unsigned int  component = 0) const;

    /**
     * Compute the Laplacian of a
     * given component at point <tt>p</tt>.
     */
    virtual double laplacian (const Point<dim>   &p,
                              const unsigned int  component = 0) const;
  private:
    /**
     * Stored Fourier coefficients.
     */
    const Point<dim> fourier_coefficients;
  };



  /**
   * Given a wavenumber vector generate a sine function. The
   * wavenumber coefficient is given as a $d$-dimensional point $k$
   * in Fourier space, and the function is then recovered as $f(x) =
   * \sin(\sum_i k_i x_i) = Im(\exp(i k.x))$.
   *
   * The class has its name from the fact that it resembles one
   * component of a Fourier sine decomposition.
   *
   * @ingroup functions
   * @author Wolfgang Bangerth, 2001
   */
  template <int dim>
  class FourierSineFunction : public Function<dim>
  {
  public:
    /**
     * Constructor. Take the Fourier
     * coefficients in each space
     * direction as argument.
     */
    FourierSineFunction (const Point<dim> &fourier_coefficients);

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
     * Return the gradient of the
     * specified component of the
     * function at the given point.
     */
    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
                                    const unsigned int  component = 0) const;

    /**
     * Compute the Laplacian of a
     * given component at point <tt>p</tt>.
     */
    virtual double laplacian (const Point<dim>   &p,
                              const unsigned int  component = 0) const;
  private:
    /**
     * Stored Fourier coefficients.
     */
    const Point<dim> fourier_coefficients;
  };


  /**
   * Given a sequence of wavenumber vectors and weights generate a sum
   * of sine functions. Each wavenumber coefficient is given as a
   * $d$-dimensional point $k$ in Fourier space, and the entire
   * function is then recovered as
   * $f(x) = \sum_j w_j sin(\sum_i k_i x_i) = Im(\sum_j w_j \exp(i k.x))$.
   *
   * @ingroup functions
   * @author Wolfgang Bangerth, 2001
   */
  template <int dim>
  class FourierSineSum : public Function<dim>
  {
  public:
    /**
     * Constructor. Take the Fourier
     * coefficients in each space
     * direction as argument.
     */
    FourierSineSum (const std::vector<Point<dim> > &fourier_coefficients,
                    const std::vector<double>      &weights);

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
     * Return the gradient of the
     * specified component of the
     * function at the given point.
     */
    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
                                    const unsigned int  component = 0) const;

    /**
     * Compute the Laplacian of a
     * given component at point <tt>p</tt>.
     */
    virtual double laplacian (const Point<dim>   &p,
                              const unsigned int  component = 0) const;
  private:
    /**
     * Stored Fourier coefficients
     * and weights.
     */
    const std::vector<Point<dim> > fourier_coefficients;
    const std::vector<double>      weights;
  };



  /**
   * Given a sequence of wavenumber vectors and weights generate a sum
   * of cosine functions. Each wavenumber coefficient is given as a
   * $d$-dimensional point $k$ in Fourier space, and the entire
   * function is then recovered as
   * $f(x) = \sum_j w_j cos(\sum_i k_i x_i) = Re(\sum_j w_j \exp(i k.x))$.
   *
   * @ingroup functions
   * @author Wolfgang Bangerth, 2001
   */
  template <int dim>
  class FourierCosineSum : public Function<dim>
  {
  public:
    /**
     * Constructor. Take the Fourier
     * coefficients in each space
     * direction as argument.
     */
    FourierCosineSum (const std::vector<Point<dim> > &fourier_coefficients,
                      const std::vector<double>      &weights);

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
     * Return the gradient of the
     * specified component of the
     * function at the given point.
     */
    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
                                    const unsigned int  component = 0) const;

    /**
     * Compute the Laplacian of a
     * given component at point <tt>p</tt>.
     */
    virtual double laplacian (const Point<dim>   &p,
                              const unsigned int  component = 0) const;

  private:
    /**
     * Stored Fourier coefficients
     * and weights.
     */
    const std::vector<Point<dim> > fourier_coefficients;
    const std::vector<double>      weights;
  };


  /**
   * Base function for cut-off function. This class stores the center
   * and the radius of the supporting ball of a cut-off function. It
   * also stores the number of the non-zero component, if the function
   * is vector-valued.
   *
   * @ingroup functions
   * @author Guido Kanschat, 2002
   */
  template <int dim>
  class CutOffFunctionBase : public Function<dim>
  {
  public:
    /**
     * Value used in the
     * constructor of this and
     * derived classes to denote
     * that no component is
     * selected.
     */
    static const unsigned int no_component = numbers::invalid_unsigned_int;

    /**
     * Constructor. Arguments are the
     * center of the ball and its
     * radius.
     *
     * If an argument <tt>select</tt> is
     * given and not -1, the
     * cut-off function will be
     * non-zero for this component
     * only.
     */
    CutOffFunctionBase (const double radius = 1.,
                        const Point<dim> = Point<dim>(),
                        const unsigned int n_components = 1,
                        const unsigned int select = CutOffFunctionBase<dim>::no_component);

    /**
     * Move the center of the ball
     * to new point <tt>p</tt>.
     */
    void new_center (const Point<dim> &p);

    /**
     * Set the radius of the ball to <tt>r</tt>.
     */
    void new_radius (const double r);

  protected:
    /**
     * Center of the integration ball.
     */
    Point<dim> center;

    /**
     * Radius of the ball.
     */
    double radius;

    /**
     * Component selected. If
     * <tt>no_component</tt>, the function is
     * the same in all components.
     */
    const unsigned int selected;
  };



  /**
   * Cut-off function in L-infinity for an arbitrary ball.  This
   * function is the characteristic function of a ball around <tt>center</tt>
   * with a specified <tt>radius</tt>, that is,
   * \f[
   * f = \chi(B_r(c)).
   * \f]
   * If vector valued, it can be restricted
   * to a single component.
   *
   * @ingroup functions
   * @author Guido Kanschat, 2001, 2002
   */
  template<int dim>
  class CutOffFunctionLinfty : public CutOffFunctionBase<dim>
  {
  public:
    /**
     * Constructor. Arguments are the
     * center of the ball and its
     * radius.
     *
     * If an argument <tt>select</tt> is
     * given and not -1, the
     * cut-off function will be
     * non-zero for this component
     * only.
     */
    CutOffFunctionLinfty (const double radius = 1.,
                          const Point<dim> = Point<dim>(),
                          const unsigned int n_components = 1,
                          const unsigned int select = CutOffFunctionBase<dim>::no_component);

    /**
     * Function value at one point.
     */
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    /**
     * Function values at multiple points.
     */
    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<double>            &values,
                             const unsigned int              component = 0) const;

    /**
     * Function values at multiple points.
     */
    virtual void vector_value_list (const std::vector<Point<dim> > &points,
                                    std::vector<Vector<double> >           &values) const;
  };


  /**
   * Cut-off function for an arbitrary ball. This function is a cone
   * with support in a ball of certain <tt>radius</tt> around <tt>center</tt>. The
   * maximum value is 1. If vector valued, it can be restricted
   * to a single component.
   *
   * @ingroup functions
   * @author Guido Kanschat, 2001, 2002
   */
  template<int dim>
  class CutOffFunctionW1 : public CutOffFunctionBase<dim>
  {
  public:
    /**
     * Constructor. Arguments are the
     * center of the ball and its
     * radius.
     * radius.
     *
     * If an argument <tt>select</tt> is
     * given, the cut-off function
     * will be non-zero for this
     * component only.
     */
    CutOffFunctionW1 (const double radius = 1.,
                      const Point<dim> = Point<dim>(),
                      const unsigned int n_components = 1,
                      const unsigned int select = CutOffFunctionBase<dim>::no_component);

    /**
     * Function value at one point.
     */
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    /**
     * Function values at multiple points.
     */
    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<double>            &values,
                             const unsigned int              component = 0) const;

    /**
     * Function values at multiple points.
     */
    virtual void vector_value_list (const std::vector<Point<dim> > &points,
                                    std::vector<Vector<double> >           &values) const;
  };


  /**
   * Cut-off function for an arbitrary ball. This is the traditional
   * cut-off function in C-infinity for a ball of certain <tt>radius</tt>
   * around <tt>center</tt>, $f(r)=exp(1-1/(1-r**2/s**2))$, where $r$ is the
   * distance to the center, and $s$ is the radius of the sphere. If
   * vector valued, it can be restricted to a single component.
   *
   * @ingroup functions
   * @author Guido Kanschat, 2001, 2002
   */
  template<int dim>
  class CutOffFunctionCinfty : public CutOffFunctionBase<dim>
  {
  public:
    /**
     * Constructor. Arguments are the
     * center of the ball and its
     * radius.
     * radius.
     *
     * If an argument <tt>select</tt> is
     * given, the cut-off function
     * will be non-zero for this
     * component only.
     */
    CutOffFunctionCinfty (const double radius = 1.,
                          const Point<dim> = Point<dim>(),
                          const unsigned int n_components = 1,
                          const unsigned int select = CutOffFunctionBase<dim>::no_component);

    /**
     * Function value at one point.
     */
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    /**
     * Function values at multiple points.
     */
    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<double>            &values,
                             const unsigned int              component = 0) const;

    /**
     * Function values at multiple points.
     */
    virtual void vector_value_list (const std::vector<Point<dim> > &points,
                                    std::vector<Vector<double> >           &values) const;

    /**
     * Function gradient at one point.
     */
    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
                                    const unsigned int  component = 0) const;
  };



  /**
   * A class that represents a function object for a monomial. Monomials are
   * polynomials with only a single term, i.e. in 1-d they have the form
   * $x^\alpha$, in 2-d the form $x_1^{\alpha_1}x_2^{\alpha_2}$, and in 3-d
   * $x_1^{\alpha_1}x_2^{\alpha_2}x_3^{\alpha_3}$. Monomials are therefore
   * described by a $dim$-tuple of exponents. Consequently, the class's
   * constructor takes a Tensor<1,dim> to describe the set of exponents. Most of
   * the time these exponents will of course be integers, but real exponents are
   * of course equally valid.
   *
   * @author Wolfgang Bangerth, 2006
   */
  template <int dim>
  class Monomial : public Function<dim>
  {
  public:
    /**
     * Constructor. The first argument is
     * explained in the general description
     * of the class. The second argument
     * denotes the number of vector
     * components this object shall
     * represent. All vector components
     * will have the same value.
     */
    Monomial (const Tensor<1,dim> &exponents,
              const unsigned int n_components = 1);

    /**
     * Function value at one point.
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
     * Function values at multiple points.
     */
    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<double>            &values,
                             const unsigned int              component = 0) const;

    /**
     * Function gradient at one point.
     */
    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
                                    const unsigned int  component = 0) const;

  private:
    /**
     * The set of exponents.
     */
    const Tensor<1,dim> exponents;
  };



  /**
   * A scalar function that computes its values by (bi-, tri-)linear interpolation
   * from a set of point data that are arranged on a possibly non-uniform
   * tensor product mesh. In other words, considering the three-dimensional case,
   * let there be points $x_0,\ldotx, x_{K-1}$,
   * $y_0,\ldots,y_{L-1}$, $z_1,\ldots,z_{M-1}$, and data $d_{klm}$ defined at
   * point $(x_k,y_l,z_m)^T$, then evaluating the function at a point
   * $\mathbf x=(x,y,z)$ will find the box so that
   * $x_k\le x\le x_{k+1}, y_l\le x\le y_{l+1}, z_m\le z\le z_{m+1}$, and do a
   * trilinear interpolation of the data on this cell. Similar operations are
   * done in lower dimensions.
   *
   * This class is most often used for either evaluating coefficients or right
   * hand sides that are provided experimentally at a number of points inside the
   * domain, or for comparing outputs of a solution on a finite element mesh
   * against previously obtained data defined on a grid.
   *
   * @note If the points $x_i$ are actually equally spaced on an interval $[x_0,x_1]$
   * and the same is true for the other data points in higher dimensions, you should
   * use the InterpolatedUniformGridData class instead.
   *
   * If a point is requested outside the box defined by the end points of the
   * coordinate arrays, then the function is assumed to simply extend by
   * constant values beyond the last data point in each coordinate
   * direction. (The class does not throw an error if a point lies outside the
   * box since it frequently happens that a point lies just outside the box
   * by an amount on the order of numerical roundoff.)
   *
   * @note The use of the related class InterpolatedUniformGridData
   * is discussed in step-53.
   *
   * @author Wolfgang Bangerth, 2013
   */
  template <int dim>
  class InterpolatedTensorProductGridData : public Function<dim>
  {
  public:
    /**
     * Constructor.
     * @param coordinate_values An array of dim arrays. Each of the inner
     *   arrays contains the coordinate values $x_0,\ldotx, x_{K-1}$ and
     *   similarly for the other coordinate directions. These arrays
     *   need not have the same size. Obviously, we need dim such arrays
     *   for a dim-dimensional function object. The coordinate values
     *   within this array are assumed to be strictly ascending to allow
     *   for efficient lookup.
     * @param data_values A dim-dimensional table of data at each of the
     *   mesh points defined by the coordinate arrays above. Note that the
     *   Table class has a number of conversion constructors that allow
     *   converting other data types into a table where you specify this
     *   argument.
     */
    InterpolatedTensorProductGridData (const std_cxx11::array<std::vector<double>,dim> &coordinate_values,
                                       const Table<dim,double>                         &data_values);

    /**
     * Compute the value of the function set by bilinear interpolation of the
     * given data set.
     *
     * @param p The point at which the function is to be evaluated.
     * @param component The vector component. Since this function is scalar,
     *   only zero is a valid argument here.
     * @return The interpolated value at this point. If the point lies outside
     *   the set of coordinates, the function is extended by a constant.
     */
    virtual
    double
    value (const Point<dim> &p,
           const unsigned int component = 0) const;

  private:
    /**
     * The set of coordinate values in each of the coordinate directions.
     */
    const std_cxx11::array<std::vector<double>,dim> coordinate_values;

    /**
     * The data that is to be interpolated.
     */
    const Table<dim,double>                     data_values;
  };


  /**
   * A scalar function that computes its values by (bi-, tri-)linear interpolation
   * from a set of point data that are arranged on a uniformly spaced
   * tensor product mesh. In other words, considering the three-dimensional case,
   * let there be points $x_0,\ldotx, x_{K-1}$ that result from a uniform subdivision
   * of the interval $[x_0,x_{K-1}]$ into $K-1$ sub-intervals of size $\Delta x
   * = (x_{K-1}-x_0)/(K-1)$, and similarly
   * $y_0,\ldots,y_{L-1}$, $z_1,\ldots,z_{M-1}$. Also consider data $d_{klm}$ defined at
   * point $(x_k,y_l,z_m)^T$, then evaluating the function at a point
   * $\mathbf x=(x,y,z)$ will find the box so that
   * $x_k\le x\le x_{k+1}, y_l\le x\le y_{l+1}, z_m\le z\le z_{m+1}$, and do a
   * trilinear interpolation of the data on this cell. Similar operations are
   * done in lower dimensions.
   *
   * This class is most often used for either evaluating coefficients or right
   * hand sides that are provided experimentally at a number of points inside the
   * domain, or for comparing outputs of a solution on a finite element mesh
   * against previously obtained data defined on a grid.
   *
   * @note If you have a problem where the points $x_i$ are not equally spaced
   * (e.g., they result from a computation on a graded mesh that is denser
   * closer to one boundary), then use the InterpolatedTensorProductGridData
   * class instead.
   *
   * If a point is requested outside the box defined by the end points of the
   * coordinate arrays, then the function is assumed to simply extend by
   * constant values beyond the last data point in each coordinate
   * direction. (The class does not throw an error if a point lies outside the
   * box since it frequently happens that a point lies just outside the box
   * by an amount on the order of numerical roundoff.)
   *
   * @note The use of this class is discussed in step-53.
   *
   * @author Wolfgang Bangerth, 2013
   */
  template <int dim>
  class InterpolatedUniformGridData : public Function<dim>
  {
  public:
    /**
     * Constructor
     * @param interval_endpoints The left and right end points of the (uniformly
     *   subdivided) intervals in each of the coordinate directions.
     * @param n_subdivisions The number of subintervals of the subintervals
     *   in each coordinate direction. A value of one for a coordinate
     *   means that the interval is considered as one subinterval consisting
     *   of the entire range. A value of two means that there are two subintervals
     *   each with one half of the range, etc.
     * @param data_values A dim-dimensional table of data at each of the
     *   mesh points defined by the coordinate arrays above. Note that the
     *   Table class has a number of conversion constructors that allow
     *   converting other data types into a table where you specify this
     *   argument.
     */
    InterpolatedUniformGridData (const std_cxx11::array<std::pair<double,double>,dim> &interval_endpoints,
                                 const std_cxx11::array<unsigned int,dim>             &n_subintervals,
                                 const Table<dim,double>                              &data_values);

    /**
     * Compute the value of the function set by bilinear interpolation of the
     * given data set.
     *
     * @param p The point at which the function is to be evaluated.
     * @param component The vector component. Since this function is scalar,
     *   only zero is a valid argument here.
     * @return The interpolated value at this point. If the point lies outside
     *   the set of coordinates, the function is extended by a constant.
     */
    virtual
    double
    value (const Point<dim> &p,
           const unsigned int component = 0) const;

  private:
    /**
     * The set of interval endpoints in each of the coordinate directions.
     */
    const std_cxx11::array<std::pair<double,double>,dim> interval_endpoints;

    /**
     * The number of subintervals in each of the coordinate directions.
     */
    const std_cxx11::array<unsigned int,dim>             n_subintervals;

    /**
     * The data that is to be interpolated.
     */
    const Table<dim,double>                     data_values;
  };
}
DEAL_II_NAMESPACE_CLOSE

#endif
