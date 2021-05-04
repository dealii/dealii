// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_function_lib_h
#define dealii_function_lib_h


#include <deal.II/base/config.h>

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/table.h>

#include <array>

DEAL_II_NAMESPACE_OPEN

/**
 * Namespace implementing some concrete classes derived from the Function
 * class that describe actual functions. This is rather a collection of
 * classes that we have needed for our own programs once and thought might be
 * useful to others as well at some point.
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
   */
  template <int dim>
  class SquareFunction : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;
    virtual void
    vector_value(const Point<dim> &p, Vector<double> &values) const override;
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const override;
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;
    virtual void
    vector_gradient(const Point<dim> &           p,
                    std::vector<Tensor<1, dim>> &gradient) const override;
    virtual void
    gradient_list(const std::vector<Point<dim>> &points,
                  std::vector<Tensor<1, dim>> &  gradients,
                  const unsigned int             component = 0) const override;
    virtual double
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;
    virtual void
    laplacian_list(const std::vector<Point<dim>> &points,
                   std::vector<double> &          values,
                   const unsigned int             component = 0) const override;
  };



  /**
   * The function <tt>xy</tt> in 2d and 3d, not implemented in 1d. This
   * function serves as an example for a vanishing Laplacian.
   *
   * @ingroup functions
   */
  template <int dim>
  class Q1WedgeFunction : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const override;

    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>> &  values) const override;

    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    virtual void
    gradient_list(const std::vector<Point<dim>> &points,
                  std::vector<Tensor<1, dim>> &  gradients,
                  const unsigned int             component = 0) const override;

    virtual void
    vector_gradient_list(
      const std::vector<Point<dim>> &,
      std::vector<std::vector<Tensor<1, dim>>> &) const override;

    /**
     * Laplacian of the function at one point.
     */
    virtual double
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;

    /**
     * Laplacian of the function at multiple points.
     */
    virtual void
    laplacian_list(const std::vector<Point<dim>> &points,
                   std::vector<double> &          values,
                   const unsigned int             component = 0) const override;
  };



  /**
   * d-quadratic pillow on the unit hypercube.
   *
   * This is a function for testing the implementation. It has zero Dirichlet
   * boundary values on the domain $(-1,1)^d$. In the inside, it is the
   * product of $1-x_i^2$ over all space dimensions.
   *
   * Providing a non-zero argument to the constructor, the whole function can
   * be offset by a constant.
   *
   * Together with the function, its derivatives and Laplacian are defined.
   *
   * @ingroup functions
   */
  template <int dim>
  class PillowFunction : public Function<dim>
  {
  public:
    /**
     * Constructor. Provide a constant that will be added to each function
     * value.
     */
    PillowFunction(const double offset = 0.);

    /**
     * The value at a single point.
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * Values at multiple points.
     */
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const override;

    /**
     * Gradient at a single point.
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    /**
     * Gradients at multiple points.
     */
    virtual void
    gradient_list(const std::vector<Point<dim>> &points,
                  std::vector<Tensor<1, dim>> &  gradients,
                  const unsigned int             component = 0) const override;

    /**
     * Laplacian at a single point.
     */
    virtual double
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;

    /**
     * Laplacian at multiple points.
     */
    virtual void
    laplacian_list(const std::vector<Point<dim>> &points,
                   std::vector<double> &          values,
                   const unsigned int             component = 0) const override;

  private:
    const double offset;
  };



  /**
   * Cosine-shaped pillow function. This is another function with zero
   * boundary values on $[-1,1]^d$. In the interior it is the product of
   * $\cos(\pi/2 x_i)$.
   *
   * @ingroup functions
   */
  template <int dim>
  class CosineFunction : public Function<dim>
  {
  public:
    /**
     * Constructor which allows to optionally generate a vector valued cosine
     * function with the same value in each component.
     */
    CosineFunction(const unsigned int n_components = 1);

    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const override;

    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>> &  values) const override;

    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    virtual void
    gradient_list(const std::vector<Point<dim>> &points,
                  std::vector<Tensor<1, dim>> &  gradients,
                  const unsigned int             component = 0) const override;

    virtual double
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;

    virtual void
    laplacian_list(const std::vector<Point<dim>> &points,
                   std::vector<double> &          values,
                   const unsigned int             component = 0) const override;

    /**
     * Second derivatives at a single point.
     */
    virtual SymmetricTensor<2, dim>
    hessian(const Point<dim> & p,
            const unsigned int component = 0) const override;

    /**
     * Second derivatives at multiple points.
     */
    virtual void
    hessian_list(const std::vector<Point<dim>> &       points,
                 std::vector<SymmetricTensor<2, dim>> &hessians,
                 const unsigned int component = 0) const override;
  };



  /**
   * Gradient of the cosine-shaped pillow function.
   *
   * This is a vector-valued function with @p dim components, the gradient of
   * CosineFunction. On the square [-1,1], it has tangential boundary
   * conditions zero. Thus, it can be used to test implementations of Maxwell
   * operators without bothering about boundary terms.
   *
   * @ingroup functions
   */
  template <int dim>
  class CosineGradFunction : public Function<dim>
  {
  public:
    /**
     * Constructor, creating a function with @p dim components.
     */
    CosineGradFunction();

    virtual double
    value(const Point<dim> &p, const unsigned int component) const override;
    virtual void
    vector_value(const Point<dim> &p, Vector<double> &values) const override;
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component) const override;

    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>> &  values) const override;

    virtual Tensor<1, dim>
    gradient(const Point<dim> &p, const unsigned int component) const override;

    virtual void
    gradient_list(const std::vector<Point<dim>> &points,
                  std::vector<Tensor<1, dim>> &  gradients,
                  const unsigned int             component) const override;

    virtual void
    vector_gradient_list(
      const std::vector<Point<dim>> &           points,
      std::vector<std::vector<Tensor<1, dim>>> &gradients) const override;

    virtual double
    laplacian(const Point<dim> &p, const unsigned int component) const override;
  };



  /**
   * Product of exponential functions in each coordinate direction.
   *
   * @ingroup functions
   */
  template <int dim>
  class ExpFunction : public Function<dim>
  {
  public:
    /**
     * The value at a single point.
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * Values at multiple points.
     */
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const override;

    /**
     * Gradient at a single point.
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    /**
     * Gradients at multiple points.
     */
    virtual void
    gradient_list(const std::vector<Point<dim>> &points,
                  std::vector<Tensor<1, dim>> &  gradients,
                  const unsigned int             component = 0) const override;

    /**
     * Laplacian at a single point.
     */
    virtual double
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;

    /**
     * Laplacian at multiple points.
     */
    virtual void
    laplacian_list(const std::vector<Point<dim>> &points,
                   std::vector<double> &          values,
                   const unsigned int             component = 0) const override;
  };



  /**
   * A function that solves the Laplace equation (with specific
   * boundary values but zero right hand side) and that has a
   * singularity at the center of the L-shaped domain in 2D (i.e.,
   * at the location of the re-entrant corner of this non-convex
   * domain).
   *
   * The function is given in polar coordinates by $r^{\frac{2}{3}}
   * \sin(\frac{2}{3} \phi)$ with a singularity at the origin and
   * should be used with GridGenerator::hyper_L(). Here, $\phi$ is
   * defined as the *clockwise* angle against the positive $x$-axis.
   *
   * This function is often used to illustrate that the solutions of the Laplace
   * equation
   * @f[
   *   -\Delta u = 0
   * @f]
   * can be singular even if the boundary values are smooth. (Here, if the
   * domain is the L-shaped domain $(-1,1)^2 \backslash [0,1]^2$, the
   * boundary values for $u$ are zero on the two line segments adjacent to the
   * origin, and equal to $r^{\frac{2}{3}} \sin(\frac{2}{3} \phi)$ on the
   * remaining parts of the boundary.) The function itself remains bounded on
   * the domain, but its gradient is of the form $r^{-1/3}$ in the vicinity of
   * the origin and consequently diverges as one approaches the origin.
   *
   * @ingroup functions
   */
  class LSingularityFunction : public Function<2>
  {
  public:
    virtual double
    value(const Point<2> &p, const unsigned int component = 0) const override;

    virtual void
    value_list(const std::vector<Point<2>> &points,
               std::vector<double> &        values,
               const unsigned int           component = 0) const override;

    virtual void
    vector_value_list(const std::vector<Point<2>> &points,
                      std::vector<Vector<double>> &values) const override;

    virtual Tensor<1, 2>
    gradient(const Point<2> &   p,
             const unsigned int component = 0) const override;

    virtual void
    gradient_list(const std::vector<Point<2>> &points,
                  std::vector<Tensor<1, 2>> &  gradients,
                  const unsigned int           component = 0) const override;

    virtual void
    vector_gradient_list(
      const std::vector<Point<2>> &,
      std::vector<std::vector<Tensor<1, 2>>> &) const override;

    virtual double
    laplacian(const Point<2> &   p,
              const unsigned int component = 0) const override;

    virtual void
    laplacian_list(const std::vector<Point<2>> &points,
                   std::vector<double> &        values,
                   const unsigned int           component = 0) const override;
  };



  /**
   * Gradient of the harmonic singularity on the L-shaped domain in 2D.
   *
   * The gradient of LSingularityFunction, which is a vector valued function
   * with vanishing curl and divergence.
   *
   * @ingroup functions
   */
  class LSingularityGradFunction : public Function<2>
  {
  public:
    /**
     * Default constructor setting the dimension to 2.
     */
    LSingularityGradFunction();
    virtual double
    value(const Point<2> &p, const unsigned int component) const override;

    virtual void
    value_list(const std::vector<Point<2>> &points,
               std::vector<double> &        values,
               const unsigned int           component) const override;

    virtual void
    vector_value_list(const std::vector<Point<2>> &points,
                      std::vector<Vector<double>> &values) const override;

    virtual Tensor<1, 2>
    gradient(const Point<2> &p, const unsigned int component) const override;

    virtual void
    gradient_list(const std::vector<Point<2>> &points,
                  std::vector<Tensor<1, 2>> &  gradients,
                  const unsigned int           component) const override;

    virtual void
    vector_gradient_list(
      const std::vector<Point<2>> &,
      std::vector<std::vector<Tensor<1, 2>>> &) const override;

    virtual double
    laplacian(const Point<2> &p, const unsigned int component) const override;

    virtual void
    laplacian_list(const std::vector<Point<2>> &points,
                   std::vector<double> &        values,
                   const unsigned int           component) const override;
  };



  /**
   * Singularity on the slit domain in 2D and 3D.
   *
   * @ingroup functions
   */
  template <int dim>
  class SlitSingularityFunction : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const override;

    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>> &  values) const override;

    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    virtual void
    gradient_list(const std::vector<Point<dim>> &points,
                  std::vector<Tensor<1, dim>> &  gradients,
                  const unsigned int             component = 0) const override;

    virtual void
    vector_gradient_list(
      const std::vector<Point<dim>> &,
      std::vector<std::vector<Tensor<1, dim>>> &) const override;

    virtual double
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;

    virtual void
    laplacian_list(const std::vector<Point<dim>> &points,
                   std::vector<double> &          values,
                   const unsigned int             component = 0) const override;
  };


  /**
   * Singularity on the slit domain with one Neumann boundary in 2D.
   *
   * @ingroup functions
   */
  class SlitHyperSingularityFunction : public Function<2>
  {
  public:
    virtual double
    value(const Point<2> &p, const unsigned int component = 0) const override;

    virtual void
    value_list(const std::vector<Point<2>> &points,
               std::vector<double> &        values,
               const unsigned int           component = 0) const override;

    virtual void
    vector_value_list(const std::vector<Point<2>> &points,
                      std::vector<Vector<double>> &values) const override;

    virtual Tensor<1, 2>
    gradient(const Point<2> &   p,
             const unsigned int component = 0) const override;

    virtual void
    gradient_list(const std::vector<Point<2>> &points,
                  std::vector<Tensor<1, 2>> &  gradients,
                  const unsigned int           component = 0) const override;

    virtual void
    vector_gradient_list(
      const std::vector<Point<2>> &,
      std::vector<std::vector<Tensor<1, 2>>> &) const override;

    virtual double
    laplacian(const Point<2> &   p,
              const unsigned int component = 0) const override;

    virtual void
    laplacian_list(const std::vector<Point<2>> &points,
                   std::vector<double> &        values,
                   const unsigned int           component = 0) const override;
  };



  /**
   * A jump in x-direction transported into some direction.
   *
   * If the advection is parallel to the y-axis, the function is
   * <tt>-atan(sx)</tt>, where <tt>s</tt> is the steepness parameter provided
   * in the constructor.
   *
   * For different advection directions, this function will be turned in the
   * parameter space.
   *
   * Together with the function, its derivatives and Laplacian are defined.
   *
   * @ingroup functions
   */
  template <int dim>
  class JumpFunction : public Function<dim>
  {
  public:
    /**
     * Constructor. Provide the advection direction here and the steepness of
     * the slope.
     */
    JumpFunction(const Point<dim> &direction, const double steepness);

    /**
     * Function value at one point.
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * Function values at multiple points.
     */
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const override;

    /**
     * Gradient at one point.
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    /**
     * Gradients at multiple points.
     */
    virtual void
    gradient_list(const std::vector<Point<dim>> &points,
                  std::vector<Tensor<1, dim>> &  gradients,
                  const unsigned int             component = 0) const override;

    /**
     * Laplacian of the function at one point.
     */
    virtual double
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;

    /**
     * Laplacian of the function at multiple points.
     */
    virtual void
    laplacian_list(const std::vector<Point<dim>> &points,
                   std::vector<double> &          values,
                   const unsigned int             component = 0) const override;

    /**
     * Return an estimate for the memory consumption, in bytes, of this
     * object. This is not exact (but will usually be close) because
     * calculating the memory usage of trees (e.g., <tt>std::map</tt>) is
     * difficult.
     */
    virtual std::size_t
    memory_consumption() const override;

  protected:
    /**
     * Advection vector.
     */
    const Point<dim> direction;

    /**
     * Steepness (maximal derivative) of the slope.
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
   * Given a wavenumber vector generate a cosine function. The wavenumber
   * coefficient is given as a $d$-dimensional point $k$ in Fourier space, and
   * the function is then recovered as $f(x) = \cos(\sum_i k_i x_i) =
   * Re(\exp(i k.x))$.
   *
   * The class has its name from the fact that it resembles one component of a
   * Fourier cosine decomposition.
   *
   * @ingroup functions
   */
  template <int dim>
  class FourierCosineFunction : public Function<dim>
  {
  public:
    /**
     * Constructor. Take the Fourier coefficients in each space direction as
     * argument.
     */
    FourierCosineFunction(const Tensor<1, dim> &fourier_coefficients);

    /**
     * Return the value of the function at the given point. Unless there is
     * only one component (i.e. the function is scalar), you should state the
     * component you want to have evaluated; it defaults to zero, i.e. the
     * first component.
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * Return the gradient of the specified component of the function at the
     * given point.
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    /**
     * Compute the Laplacian of a given component at point <tt>p</tt>.
     */
    virtual double
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;

  private:
    /**
     * Stored Fourier coefficients.
     */
    const Tensor<1, dim> fourier_coefficients;
  };



  /**
   * Given a wavenumber vector generate a sine function. The wavenumber
   * coefficient is given as a $d$-dimensional point $k$ in Fourier space, and
   * the function is then recovered as $f(x) = \sin(\sum_i k_i x_i) =
   * Im(\exp(i k.x))$.
   *
   * The class has its name from the fact that it resembles one component of a
   * Fourier sine decomposition.
   *
   * @ingroup functions
   */
  template <int dim>
  class FourierSineFunction : public Function<dim>
  {
  public:
    /**
     * Constructor. Take the Fourier coefficients in each space direction as
     * argument.
     */
    FourierSineFunction(const Tensor<1, dim> &fourier_coefficients);

    /**
     * Return the value of the function at the given point. Unless there is
     * only one component (i.e. the function is scalar), you should state the
     * component you want to have evaluated; it defaults to zero, i.e. the
     * first component.
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * Return the gradient of the specified component of the function at the
     * given point.
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    /**
     * Compute the Laplacian of a given component at point <tt>p</tt>.
     */
    virtual double
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;

  private:
    /**
     * Stored Fourier coefficients.
     */
    const Tensor<1, dim> fourier_coefficients;
  };


  /**
   * Given a sequence of wavenumber vectors and weights generate a sum of sine
   * functions. Each wavenumber coefficient is given as a $d$-dimensional
   * point $k$ in Fourier space, and the entire function is then recovered as
   * $f(x) = \sum_j w_j sin(\sum_i k_i x_i) = Im(\sum_j w_j \exp(i k.x))$.
   *
   * @ingroup functions
   */
  template <int dim>
  class FourierSineSum : public Function<dim>
  {
  public:
    /**
     * Constructor. Take the Fourier coefficients in each space direction as
     * argument.
     */
    FourierSineSum(const std::vector<Point<dim>> &fourier_coefficients,
                   const std::vector<double> &    weights);

    /**
     * Return the value of the function at the given point. Unless there is
     * only one component (i.e. the function is scalar), you should state the
     * component you want to have evaluated; it defaults to zero, i.e. the
     * first component.
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * Return the gradient of the specified component of the function at the
     * given point.
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    /**
     * Compute the Laplacian of a given component at point <tt>p</tt>.
     */
    virtual double
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;

  private:
    /**
     * Stored Fourier coefficients and weights.
     */
    const std::vector<Point<dim>> fourier_coefficients;
    const std::vector<double>     weights;
  };



  /**
   * Given a sequence of wavenumber vectors and weights generate a sum of
   * cosine functions. Each wavenumber coefficient is given as a
   * $d$-dimensional point $k$ in Fourier space, and the entire function is
   * then recovered as $f(x) = \sum_j w_j cos(\sum_i k_i x_i) = Re(\sum_j w_j
   * \exp(i k.x))$.
   *
   * @ingroup functions
   */
  template <int dim>
  class FourierCosineSum : public Function<dim>
  {
  public:
    /**
     * Constructor. Take the Fourier coefficients in each space direction as
     * argument.
     */
    FourierCosineSum(const std::vector<Point<dim>> &fourier_coefficients,
                     const std::vector<double> &    weights);

    /**
     * Return the value of the function at the given point. Unless there is
     * only one component (i.e. the function is scalar), you should state the
     * component you want to have evaluated; it defaults to zero, i.e. the
     * first component.
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * Return the gradient of the specified component of the function at the
     * given point.
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    /**
     * Compute the Laplacian of a given component at point <tt>p</tt>.
     */
    virtual double
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;

  private:
    /**
     * Stored Fourier coefficients and weights.
     */
    const std::vector<Point<dim>> fourier_coefficients;
    const std::vector<double>     weights;
  };


  /**
   * Base function for cut-off function. This class stores the center and the
   * radius of the supporting ball of a cut-off function. It also stores the
   * number of the non-zero component, if the function is vector-valued.
   *
   * This class can also be used for approximated Dirac delta functions. These
   * are special cut-off functions whose integral is always equal to one,
   * independently of the radius of the supporting ball.
   *
   * @ingroup functions
   */
  template <int dim>
  class CutOffFunctionBase : public Function<dim>
  {
  public:
    /**
     * Value used in the constructor of this and derived classes to denote
     * that no component is selected.
     */
    static const unsigned int no_component = numbers::invalid_unsigned_int;

    /**
     * Constructor.
     *
     * @param[in] radius Radius of the ball
     * @param[in] center Center of the ball
     * @param[in] n_components Number of components of this function object
     * @param[in] select If this is different from
     * CutOffFunctionBase<dim>::no_component, then the function will be non-zero
     * for this component only
     * @param[in] integrate_to_one Rescale the value of the function whenever a
     * new radius is set, to guarantee that the integral is equal to one
     * @param[in] unitary_integral_value Value of the integral when the radius
     * is equal to 1.0. Derived classes will need to supply this value, to
     * guarantee that the rescaling is performed correctly.
     */
    CutOffFunctionBase(
      const double       radius       = 1.,
      const Point<dim>   center       = Point<dim>(),
      const unsigned int n_components = 1,
      const unsigned int select       = CutOffFunctionBase<dim>::no_component,
      const bool         integrate_to_one       = false,
      const double       unitary_integral_value = 1.0);

    /**
     * Virtual destructor.
     */
    virtual ~CutOffFunctionBase() = default;

    /**
     * Set the center of the ball to the point @p p.
     */
    virtual void
    set_center(const Point<dim> &p);

    /**
     * Set the radius of the ball to @p r
     */
    virtual void
    set_radius(const double r);

    /**
     * Return the center stored in this object.
     */
    const Point<dim> &
    get_center() const;

    /**
     * Return the radius stored in this object.
     */
    double
    get_radius() const;

    /**
     * Return a boolean indicating whether this function integrates to one.
     */
    bool
    integrates_to_one() const;

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
     * Component selected. If <tt>no_component</tt>, the function is the same
     * in all components.
     */
    const unsigned int selected;

    /**
     * Flag that controls whether we rescale the value when the radius changes.
     */
    bool integrate_to_one;

    /**
     * The reference integral value. Derived classes should specify what their
     * integral is when @p radius = 1.0.
     */
    const double unitary_integral_value;

    /**
     * Current rescaling to apply the cut-off function.
     */
    double rescaling;
  };


  /**
   * Tensor product of CutOffFunctionBase objects.
   *
   * Instead of using the distance to compute the cut-off function, this class
   * performs a tensor product of the same CutOffFunctionBase object in each
   * coordinate direction.
   *
   * @ingroup functions
   */
  template <int dim>
  class CutOffFunctionTensorProduct : public CutOffFunctionBase<dim>
  {
  public:
    /**
     * Construct an empty CutOffFunctionTensorProduct object.
     *
     * Before you can use this class, you have to call the set_base() method
     * with a class derived from the CutOffFunctionBase object.
     *
     * If you try to use this class before you call the set_base() method,
     * and exception will be triggered.
     */
    CutOffFunctionTensorProduct(
      double             radius       = 1.0,
      const Point<dim> & center       = Point<dim>(),
      const unsigned int n_components = 1,
      const unsigned int select       = CutOffFunctionBase<dim>::no_component,
      const bool         integrate_to_one = false);

    /**
     * Initialize the class with an object of type
     * @tparam CutOffFunctionBaseType<1>.
     */
    template <template <int> class CutOffFunctionBaseType>
    void
    set_base();

    /**
     * Set the new center.
     */
    virtual void
    set_center(const Point<dim> &center) override;

    /**
     * Set the new radius.
     */
    virtual void
    set_radius(const double radius) override;

    /**
     * Function value at one point.
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * Function gradient at one point.
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

  private:
    std::array<std::unique_ptr<CutOffFunctionBase<1>>, dim> base;

    bool initialized;
  };



  /**
   * Cut-off function in L-infinity for an arbitrary ball.  This function is
   * the characteristic function of a ball around <tt>center</tt> with a
   * specified <tt>radius</tt>, that is, \f[ f = \chi(B_r(c)). \f] If vector
   * valued, it can be restricted to a single component.
   *
   * @ingroup functions
   */
  template <int dim>
  class CutOffFunctionLinfty : public CutOffFunctionBase<dim>
  {
  public:
    /**
     * Constructor. Arguments are the center of the ball and its radius.
     *
     * If an argument <tt>select</tt> is given and not -1, the cut-off
     * function will be non-zero for this component only.
     */
    CutOffFunctionLinfty(
      const double radius             = 1.,
      const Point<dim>                = Point<dim>(),
      const unsigned int n_components = 1,
      const unsigned int select       = CutOffFunctionBase<dim>::no_component,
      const bool         integrate_to_one = false);

    /**
     * Function value at one point.
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * Function values at multiple points.
     */
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const override;

    /**
     * Function values at multiple points.
     */
    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>> &  values) const override;
  };


  /**
   * Cut-off function for an arbitrary ball. This function is a cone with
   * support in a ball of certain <tt>radius</tt> around <tt>center</tt>. The
   * maximum value is 1. If vector valued, it can be restricted to a single
   * component.
   *
   * @ingroup functions
   */
  template <int dim>
  class CutOffFunctionW1 : public CutOffFunctionBase<dim>
  {
  public:
    /**
     * Constructor. Arguments are the center of the ball and its radius.
     *
     * If an argument <tt>select</tt> is given, the cut-off function will be
     * non-zero for this component only.
     */
    CutOffFunctionW1(
      const double radius             = 1.,
      const Point<dim>                = Point<dim>(),
      const unsigned int n_components = 1,
      const unsigned int select       = CutOffFunctionBase<dim>::no_component,
      const bool         integrate_to_one = false);

    /**
     * Function value at one point.
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * Function values at multiple points.
     */
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const override;

    /**
     * Function values at multiple points.
     */
    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>> &  values) const override;
  };


  /**
   * A cut-off function for an arbitrarily-sized ball that is in the space $C^1$
   * (i.e., continuously differentiable). This is a cut-off function that is
   * often used in the literature of the Immersed Boundary Method.
   *
   * The expression of the function in radial coordinates is given by
   * $f(r)=1/2(cos(\pi r/s)+1)$ where $r<s$ is the distance to the center, and
   * $s$ is the radius of the sphere. If vector valued, it can be restricted to
   * a single component.
   *
   * @ingroup functions
   */
  template <int dim>
  class CutOffFunctionC1 : public CutOffFunctionBase<dim>
  {
  public:
    /**
     * Constructor.
     */
    CutOffFunctionC1(
      const double radius             = 1.,
      const Point<dim>                = Point<dim>(),
      const unsigned int n_components = 1,
      const unsigned int select       = CutOffFunctionBase<dim>::no_component,
      bool               integrate_to_one = false);

    /**
     * Function value at one point.
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * Function values at multiple points.
     */
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const override;

    /**
     * Function values at multiple points.
     */
    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>> &  values) const override;

    /**
     * Function gradient at one point.
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;
  };


  /**
   * Cut-off function for an arbitrary ball. This is the traditional cut-off
   * function in C-infinity for a ball of certain <tt>radius</tt> around
   * <tt>center</tt>, $f(r)=exp(1-1/(1-r**2/s**2))$, where $r$ is the distance
   * to the center, and $s$ is the radius of the sphere. If vector valued, it
   * can be restricted to a single component.
   *
   * @ingroup functions
   */
  template <int dim>
  class CutOffFunctionCinfty : public CutOffFunctionBase<dim>
  {
  public:
    /**
     * Constructor. Arguments are the center of the ball and its radius.
     *
     * If an argument <tt>select</tt> is given, the cut-off function will be
     * non-zero for this component only.
     */
    CutOffFunctionCinfty(
      const double radius             = 1.,
      const Point<dim>                = Point<dim>(),
      const unsigned int n_components = 1,
      const unsigned int select       = CutOffFunctionBase<dim>::no_component,
      bool               integrate_to_one = false);

    /**
     * Function value at one point.
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * Function values at multiple points.
     */
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const override;

    /**
     * Function values at multiple points.
     */
    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>> &  values) const override;

    /**
     * Function gradient at one point.
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;
  };



  /**
   * A class that represents a function object for a monomial. Monomials are
   * polynomials with only a single term, i.e. in 1-d they have the form
   * $x^\alpha$, in 2-d the form $x_1^{\alpha_1}x_2^{\alpha_2}$, and in 3-d
   * $x_1^{\alpha_1}x_2^{\alpha_2}x_3^{\alpha_3}$. Monomials are therefore
   * described by a $dim$-tuple of exponents. Consequently, the class's
   * constructor takes a Tensor<1,dim> to describe the set of exponents. Most
   * of the time these exponents will of course be integers, but real
   * exponents are of course equally valid. Exponents can't be real when the
   * bases are negative numbers.
   *
   * @ingroup functions
   */
  template <int dim>
  class Monomial : public Function<dim>
  {
  public:
    /**
     * Constructor. The first argument is explained in the general description
     * of the class. The second argument denotes the number of vector
     * components this object shall represent. All vector components will have
     * the same value.
     */
    Monomial(const Tensor<1, dim> &exponents,
             const unsigned int    n_components = 1);

    /**
     * Function value at one point.
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * Return all components of a vector-valued function at a given point.
     *
     * <tt>values</tt> shall have the right size beforehand, i.e.
     * #n_components.
     */
    virtual void
    vector_value(const Point<dim> &p, Vector<double> &values) const override;

    /**
     * Function values at multiple points.
     */
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const override;

    /**
     * Function gradient at one point.
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

  private:
    /**
     * The set of exponents.
     */
    const Tensor<1, dim> exponents;
  };



  /**
   * A scalar function that computes its values by (bi-, tri-)linear
   * interpolation from a set of point data that are arranged on a possibly
   * non-uniform tensor product mesh. In other words, considering the three-
   * dimensional case, let there be points $x_0,\ldots, x_{K-1}$,
   * $y_0,\ldots,y_{L-1}$, $z_1,\ldots,z_{M-1}$, and data $d_{klm}$ defined at
   * point $(x_k,y_l,z_m)^T$, then evaluating the function at a point $\mathbf
   * x=(x,y,z)$ will find the box so that $x_k\le x\le x_{k+1}, y_l\le y\le
   * y_{l+1}, z_m\le z\le z_{m+1}$, and do a trilinear interpolation of the
   * data on this cell. Similar operations are done in lower dimensions.
   *
   * This class is most often used for either evaluating coefficients or right
   * hand sides that are provided experimentally at a number of points inside
   * the domain, or for comparing outputs of a solution on a finite element
   * mesh against previously obtained data defined on a grid.
   *
   * @note If the points $x_i$ are actually equally spaced on an interval
   * $[x_0,x_1]$ and the same is true for the other data points in higher
   * dimensions, you should use the InterpolatedUniformGridData class instead.
   *
   * If a point is requested outside the box defined by the end points of the
   * coordinate arrays, then the function is assumed to simply extend by
   * constant values beyond the last data point in each coordinate direction.
   * (The class does not throw an error if a point lies outside the box since
   * it frequently happens that a point lies just outside the box by an amount
   * on the order of numerical roundoff.)
   *
   * @note The use of the related class InterpolatedUniformGridData is
   * discussed in step-53.
   *
   *
   * <h3>Dealing with large data sets</h3>
   *
   * This class is often used to interpolate data provided by fairly
   * large data tables that are expensive to read from disk, and that take
   * a large amount of memory when replicated on every process of parallel
   * (MPI) programs.
   *
   * The Table class can help with amortizing this cost by using
   * shared memory to store the data only as often as necessary -- see the
   * documentation of the TableBase class. Once one has obtained such a
   * Table object that uses shared memory to store the data only as often
   * as is necessary, one has to avoid that the current class *copies*
   * the table into its own member variable. Rather, it is necessary to
   * use the *move* constructor of this class to take over ownership of
   * the table and its shared memory space. This can be achieved using
   * the following extension of the code snippet shown in the
   * documentation of the TableBase class:
   * @code
   *    const unsigned int N=..., M=...;     // table sizes, assumed known
   *    Table<2,double>    data_table;
   *    const unsigned int root_rank = 0;
   *
   *    if (Utilities::MPI::this_mpi_process(mpi_communicator) == root_rank)
   *    {
   *      data_table.resize (N,M);
   *
   *      std::ifstream input_file ("data_file.dat");
   *      ...;                               // read the data from the file
   *    }
   *
   *    // Now distribute to all processes
   *    data_table.replicate_across_communicator (mpi_communicator, root_rank);
   *
   *    // Set up the x- and y-coordinates of the points stored in the
   *    // data table
   *    std::array<std::vector<double>, dim> coordinate_values;
   *    ...;                                 // do what needs to be done
   *
   *    // And finally set up the interpolation object. The calls
   *    // to std::move() make sure that the tables are moved into
   *    // the memory space of the InterpolateTensorProductGridData
   *    // object:
   *    InterpolatedTensorProductGridData<2>
   *          interpolation_function (std::move(coordinate_values),
   *                                  std::move(data_table));
   * @endcode
   *
   *
   * @ingroup functions
   */
  template <int dim>
  class InterpolatedTensorProductGridData : public Function<dim>
  {
  public:
    /**
     * Constructor to initialize this class instance with the data given in @p
     * data_values.
     *
     * @param coordinate_values An array of dim arrays. Each of the inner
     * arrays contains the coordinate values $x_0,\ldots, x_{K-1}$ and
     * similarly for the other coordinate directions. These arrays need not
     * have the same size. Obviously, we need dim such arrays for a dim-
     * dimensional function object. The coordinate values within this array
     * are assumed to be strictly ascending to allow for efficient lookup.
     *
     * @param data_values A dim-dimensional table of data at each of the mesh
     * points defined by the coordinate arrays above. The data passed in is
     * copied into internal data structures. Note that the Table
     * class has a number of conversion constructors that allow converting
     * other data types into a table where you specify this argument.
     */
    InterpolatedTensorProductGridData(
      const std::array<std::vector<double>, dim> &coordinate_values,
      const Table<dim, double> &                  data_values);

    /**
     * Like the previous constructor, but take the arguments as rvalue
     * references and *move*, instead of *copy* the data. This is often useful
     * in cases where the data stored in these tables is large and the
     * information used to initialize the current object is no longer needed
     * separately. In other words, there is no need to keep the original object
     * from which this object could copy its information, but it might as well
     * take over ("move") the data.
     */
    InterpolatedTensorProductGridData(
      std::array<std::vector<double>, dim> &&coordinate_values,
      Table<dim, double> &&                  data_values);

    /**
     * Compute the value of the function set by bilinear interpolation of the
     * given data set.
     *
     * @param p The point at which the function is to be evaluated.
     * @param component The vector component. Since this function is scalar,
     * only zero is a valid argument here.
     * @return The interpolated value at this point. If the point lies outside
     * the set of coordinates, the function is extended by a constant.
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * Compute the gradient of the function defined by bilinear interpolation
     * of the given data set.
     *
     * @param p The point at which the function gradient is to be evaluated.
     * @param component The vector component. Since this function is scalar,
     * only zero is a valid argument here.
     * @return The value of the gradient of the interpolated function at this
     * point. If the point lies outside the set of coordinates, the function
     * is extended by a constant and so its gradient is extended by 0.
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    /**
     * Return an estimate for the memory consumption, in bytes, of this object.
     */
    virtual std::size_t
    memory_consumption() const override;

    /**
     * Return a reference to the internally stored data.
     */
    const Table<dim, double> &
    get_data() const;

  protected:
    /**
     * Find the index in the table of the rectangle containing an input point
     */
    TableIndices<dim>
    table_index_of_point(const Point<dim> &p) const;

    /**
     * The set of coordinate values in each of the coordinate directions.
     */
    const std::array<std::vector<double>, dim> coordinate_values;

    /**
     * The data that is to be interpolated.
     */
    const Table<dim, double> data_values;
  };


  /**
   * A scalar function that computes its values by (bi-, tri-)linear
   * interpolation from a set of point data that are arranged on a uniformly
   * spaced tensor product mesh. In other words, considering the three-
   * dimensional case, let there be points $x_0,\ldots, x_{K-1}$ that result
   * from a uniform subdivision of the interval $[x_0,x_{K-1}]$ into $K-1$
   * sub-intervals of size $\Delta x = (x_{K-1}-x_0)/(K-1)$, and similarly
   * $y_0,\ldots,y_{L-1}$, $z_1,\ldots,z_{M-1}$. Also consider data $d_{klm}$
   * defined at point $(x_k,y_l,z_m)^T$, then evaluating the function at a
   * point $\mathbf x=(x,y,z)$ will find the box so that $x_k\le x\le x_{k+1},
   * y_l\le y\le y_{l+1}, z_m\le z\le z_{m+1}$, and do a trilinear
   * interpolation of the data on this cell. Similar operations are done in
   * lower dimensions.
   *
   * This class is most often used for either evaluating coefficients or right
   * hand sides that are provided experimentally at a number of points inside
   * the domain, or for comparing outputs of a solution on a finite element
   * mesh against previously obtained data defined on a grid.
   *
   * @note If you have a problem where the points $x_i$ are not equally spaced
   * (e.g., they result from a computation on a graded mesh that is denser
   * closer to one boundary), then use the InterpolatedTensorProductGridData
   * class instead.
   *
   * If a point is requested outside the box defined by the end points of the
   * coordinate arrays, then the function is assumed to simply extend by
   * constant values beyond the last data point in each coordinate direction.
   * (The class does not throw an error if a point lies outside the box since
   * it frequently happens that a point lies just outside the box by an amount
   * on the order of numerical roundoff.)
   *
   * @note The use of this class is discussed in step-53.
   *
   *
   * <h3>Dealing with large data sets</h3>
   *
   * This class supports the same facilities for dealing with large data sets
   * as the InterpolatedTensorProductGridData class. See there for more
   * information and example codes.
   *
   *
   * @ingroup functions
   */
  template <int dim>
  class InterpolatedUniformGridData : public Function<dim>
  {
  public:
    /**
     * Constructor
     * @param interval_endpoints The left and right end points of the
     * (uniformly subdivided) intervals in each of the coordinate directions.
     * @param n_subintervals The number of subintervals in each coordinate
     * direction. A value of one for a coordinate means that the interval is
     * considered as one subinterval consisting of the entire range. A value
     * of two means that there are two subintervals each with one half of the
     * range, etc.
     * @param data_values A dim-dimensional table of data at each of the mesh
     * points defined by the coordinate arrays above. Note that the Table
     * class has a number of conversion constructors that allow converting
     * other data types into a table where you specify this argument.
     */
    InterpolatedUniformGridData(
      const std::array<std::pair<double, double>, dim> &interval_endpoints,
      const std::array<unsigned int, dim> &             n_subintervals,
      const Table<dim, double> &                        data_values);

    /**
     * Like the previous constructor, but take the arguments as rvalue
     * references and *move*, instead of *copy* the data. This is often useful
     * in cases where the data stored in these tables is large and the
     * information used to initialize the current object is no longer needed
     * separately. In other words, there is no need to keep the original object
     * from which this object could copy its information, but it might as well
     * take over ("move") the data.
     */
    InterpolatedUniformGridData(
      std::array<std::pair<double, double>, dim> &&interval_endpoints,
      std::array<unsigned int, dim> &&             n_subintervals,
      Table<dim, double> &&                        data_values);

    /**
     * Compute the value of the function set by bilinear interpolation of the
     * given data set.
     *
     * @param p The point at which the function is to be evaluated.
     * @param component The vector component. Since this function is scalar,
     * only zero is a valid argument here.
     * @return The interpolated value at this point. If the point lies outside
     * the set of coordinates, the function is extended by a constant.
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * Compute the gradient of the function set by bilinear interpolation of the
     * given data set.
     *
     * @param p The point at which the function is to be evaluated.
     * @param component The vector component. Since this function is scalar,
     *   only zero is a valid argument here.
     * @return The gradient of the interpolated function at this point. If the
     *   point lies outside the set of coordinates, the function is extended
     *   by a constant whose gradient is then of course zero.
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    /**
     * Return an estimate for the memory consumption, in bytes, of this object.
     */
    virtual std::size_t
    memory_consumption() const override;

    /**
     * Return a reference to the internally stored data.
     */
    const Table<dim, double> &
    get_data() const;

  private:
    /**
     * The set of interval endpoints in each of the coordinate directions.
     */
    const std::array<std::pair<double, double>, dim> interval_endpoints;

    /**
     * The number of subintervals in each of the coordinate directions.
     */
    const std::array<unsigned int, dim> n_subintervals;

    /**
     * The data that is to be interpolated.
     */
    const Table<dim, double> data_values;
  };


  /**
   * A class that represents a function object for a polynomial. A polynomial
   * is composed by the summation of multiple monomials. If the polynomial has
   * n monomials and the dimension is equal to dim, the polynomial can be
   * written as $\sum_{i=1}^{n} a_{i}(\prod_{d=1}^{dim}
   * x_{d}^{\alpha_{i,d}})$, where $a_{i}$ are the coefficients of the
   * monomials and $\alpha_{i,d}$ are their exponents. The class's constructor
   * takes a Table<2,double> to describe the set of exponents and a
   * Vector<double> to describe the set of coefficients.
   *
   * @ingroup functions
   */
  template <int dim>
  class Polynomial : public Function<dim>
  {
  public:
    /**
     * Constructor. The coefficients and the exponents of the polynomial are
     * passed as arguments. The Table<2, double> exponents has a number of
     * rows equal to the number of monomials of the polynomial and a number of
     * columns equal to dim. The i-th row of the exponents table contains the
     * ${\alpha_{i,d}}$ exponents of the i-th monomial $a_{i}\prod_{d=1}^{dim}
     * x_{d}^{\alpha_{i,d}}$. The i-th element of the coefficients vector
     * contains the coefficient $a_{i}$ for the i-th monomial.
     */
    Polynomial(const Table<2, double> &   exponents,
               const std::vector<double> &coefficients);

    /**
     * Function value at one point.
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;


    /**
     * Function values at multiple points.
     */
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const override;

    /**
     * Function gradient at one point.
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    /**
     * Return an estimate for the memory consumption, in bytes, of this object.
     */
    virtual std::size_t
    memory_consumption() const override;

  private:
    /**
     * The set of exponents.
     */
    const Table<2, double> exponents;

    /**
     * The set of coefficients.
     */
    const std::vector<double> coefficients;
  };

#ifndef DOXYGEN



  // Template definitions
  template <int dim>
  template <template <int> class CutOffFunctionBaseType>
  void
  CutOffFunctionTensorProduct<dim>::set_base()
  {
    initialized = true;
    static_assert(
      std::is_base_of<CutOffFunctionBase<1>, CutOffFunctionBaseType<1>>::value,
      "You can only construct a CutOffFunctionTensorProduct function from "
      "a class derived from CutOffFunctionBase.");

    for (unsigned int i = 0; i < dim; ++i)
      base[i].reset(new CutOffFunctionBaseType<1>(this->radius,
                                                  Point<1>(this->center[i]),
                                                  this->n_components,
                                                  this->selected,
                                                  this->integrate_to_one));
  }



#endif

} // namespace Functions
DEAL_II_NAMESPACE_CLOSE

#endif
