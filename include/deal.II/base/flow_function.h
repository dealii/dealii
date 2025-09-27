// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_flow_function_h
#define dealii_flow_function_h


#include <deal.II/base/config.h>

#include <deal.II/base/function.h>
#include <deal.II/base/mutex.h>
#include <deal.II/base/point.h>

DEAL_II_NAMESPACE_OPEN

namespace Functions
{
  /**
   * Base class for analytic solutions to incompressible flow problems.
   *
   * Additional to the Function interface, this function provides for an
   * offset of the pressure: if the pressure of the computed solution has an
   * integral mean value different from zero, this value can be given to
   * pressure_adjustment() in order to compute correct pressure errors.
   *
   * @note Derived classes should implement pressures with integral mean value
   * zero always.
   *
   * @note Thread safety: Some of the functions make use of internal data to
   * compute values. Therefore, every thread should obtain its own object of
   * derived classes.
   *
   * @ingroup functions
   */
  template <int dim>
  class FlowFunction : public Function<dim>
  {
  public:
    /**
     * Constructor, setting up some internal data structures.
     */
    FlowFunction();

    /**
     * Virtual destructor.
     */
    virtual ~FlowFunction() override = default;

    /**
     * Store an adjustment for the pressure function, such that its mean value
     * is <tt>p</tt>.
     */
    void
    pressure_adjustment(double p);

    /**
     * Values in a structure more suitable for vector valued functions. The
     * outer vector is indexed by solution component, the inner by quadrature
     * point.
     */
    virtual void
    vector_values(const std::vector<Point<dim>>    &points,
                  std::vector<std::vector<double>> &values) const override = 0;
    /**
     * Gradients in a structure more suitable for vector valued functions. The
     * outer vector is indexed by solution component, the inner by quadrature
     * point.
     */
    virtual void
    vector_gradients(
      const std::vector<Point<dim>>            &points,
      std::vector<std::vector<Tensor<1, dim>>> &gradients) const override = 0;
    /**
     * Force terms in a structure more suitable for vector valued functions.
     * The outer vector is indexed by solution component, the inner by
     * quadrature point.
     *
     * @warning This is not the true Laplacian, but the force term to be used
     * as right hand side in Stokes' equations
     */
    virtual void
    vector_laplacians(const std::vector<Point<dim>>    &points,
                      std::vector<std::vector<double>> &values) const = 0;

    virtual void
    vector_value(const Point<dim> &points,
                 Vector<double>   &value) const override;
    virtual double
    value(const Point<dim>  &points,
          const unsigned int component) const override;
    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>>   &values) const override;
    virtual void
    vector_gradient_list(
      const std::vector<Point<dim>>            &points,
      std::vector<std::vector<Tensor<1, dim>>> &gradients) const override;
    /**
     * The force term in the momentum equation.
     */
    virtual void
    vector_laplacian_list(const std::vector<Point<dim>> &points,
                          std::vector<Vector<double>>   &values) const override;

    /**
     * Return an estimate for the memory consumption, in bytes, of this object.
     */
    virtual std::size_t
    memory_consumption() const override;

  protected:
    /**
     * Mean value of the pressure to be added by derived classes.
     */
    double mean_pressure;

  private:
    /**
     * A mutex that guards the following scratch arrays.
     */
    mutable Threads::Mutex mutex;

    /**
     * Auxiliary values for the usual Function interface.
     */
    mutable std::vector<std::vector<double>> aux_values;

    /**
     * Auxiliary values for the usual Function interface.
     */
    mutable std::vector<std::vector<Tensor<1, dim>>> aux_gradients;
  };

  /**
   * Laminar pipe flow in two and three dimensions. The channel stretches
   * along the <i>x</i>-axis and has radius @p radius. The @p Reynolds number
   * is used to scale the pressure properly for a Navier-Stokes problem.
   *
   * @ingroup functions
   */
  template <int dim>
  class PoisseuilleFlow : public FlowFunction<dim>
  {
  public:
    /**
     * Construct an object for the given channel radius <tt>r</tt> and the
     * Reynolds number <tt>Re</tt>.
     */
    PoisseuilleFlow(const double r, const double Re);

    virtual ~PoisseuilleFlow() override = default;

    virtual void
    vector_values(const std::vector<Point<dim>>    &points,
                  std::vector<std::vector<double>> &values) const override;
    virtual void
    vector_gradients(
      const std::vector<Point<dim>>            &points,
      std::vector<std::vector<Tensor<1, dim>>> &gradients) const override;
    virtual void
    vector_laplacians(const std::vector<Point<dim>>    &points,
                      std::vector<std::vector<double>> &values) const override;

  private:
    const double inv_sqr_radius;
    const double Reynolds;
  };


  /**
   * Artificial divergence free function with homogeneous boundary conditions
   * on the cube [-1,1]<sup>dim</sup>.
   *
   * The function in 2d is
   * @f[
   * \left(\begin{array}{c}u\\v\\p\end{array}\right)
   * \left(\begin{array}{c}\cos^2x \sin y\cos y\\-\sin x\cos x\cos^2y\\
   * \sin x\cos x\sin y\cos y\end{array}\right)
   * @f]
   * @ingroup functions
   */
  template <int dim>
  class StokesCosine : public FlowFunction<dim>
  {
  public:
    /**
     * Constructor setting the Reynolds number required for pressure
     * computation and scaling of the right hand side.
     */
    StokesCosine(const double viscosity = 1., const double reaction = 0.);
    /**
     * Change the viscosity and the reaction parameter.
     */
    void
    set_parameters(const double viscosity, const double reaction);

    virtual ~StokesCosine() override = default;

    virtual void
    vector_values(const std::vector<Point<dim>>    &points,
                  std::vector<std::vector<double>> &values) const override;
    virtual void
    vector_gradients(
      const std::vector<Point<dim>>            &points,
      std::vector<std::vector<Tensor<1, dim>>> &gradients) const override;
    virtual void
    vector_laplacians(const std::vector<Point<dim>>    &points,
                      std::vector<std::vector<double>> &values) const override;

  private:
    /// The viscosity
    double viscosity;
    /// The reaction parameter
    double reaction;
  };


  /**
   * A singular solution to Stokes' equations on a 2d L-shaped domain.
   *
   * This function satisfies $-\triangle \mathbf{u} + \nabla p = 0$ and
   * represents a typical singular solution around a reentrant corner of an
   * L-shaped domain that can be created using GridGenerator::hyper_L(). The
   * velocity vanishes on the two faces of the re-entrant corner and
   * $\nabla\mathbf{u}$ and $p$ are singular at the origin while they are
   * smooth in the rest of the domain because they can be written as a product
   * of a smooth function and the term $r^{\lambda-1}$ where $r$ is the radius
   * and $\lambda \approx 0.54448$ is a fixed parameter.
   *
   * Taken from Houston, Sch&ouml;tzau, Wihler, proceeding ENUMATH 2003.
   *
   * @ingroup functions
   */
  class StokesLSingularity : public FlowFunction<2>
  {
  public:
    /// Constructor setting up some data.
    StokesLSingularity();

    virtual void
    vector_values(const std::vector<Point<2>>      &points,
                  std::vector<std::vector<double>> &values) const override;
    virtual void
    vector_gradients(
      const std::vector<Point<2>>            &points,
      std::vector<std::vector<Tensor<1, 2>>> &gradients) const override;
    virtual void
    vector_laplacians(const std::vector<Point<2>>      &points,
                      std::vector<std::vector<double>> &values) const override;

  private:
    /// The auxiliary function Psi.
    double
    Psi(double phi) const;
    /// The derivative of Psi()
    double
    Psi_1(double phi) const;
    /// The 2nd derivative of Psi()
    double
    Psi_2(double phi) const;
    /// The 3rd derivative of Psi()
    double
    Psi_3(double phi) const;
    /// The 4th derivative of Psi()
    double
    Psi_4(double phi) const;
    /// The angle of the reentrant corner, set to 3*pi/2
    const double omega;
    /// The exponent of the radius, computed as the solution to
    /// $\sin(\lambda\omega)+\lambda \sin(\omega)=0$
    static const double lambda;
    /// Cosine of lambda times omega
    const double coslo;
    /// Auxiliary variable 1+lambda
    const double lp;
    /// Auxiliary variable 1-lambda
    const double lm;
  };

  /**
   * Flow solution in 2d by Kovasznay (1947).
   *
   * This function is valid on the half plane right of the line <i>x=1/2</i>.
   *
   * @ingroup functions
   */
  class Kovasznay : public FlowFunction<2>
  {
  public:
    /**
     * Construct an object for the give Reynolds number <tt>Re</tt>. If the
     * parameter <tt>Stokes</tt> is true, the right hand side of the momentum
     * equation returned by vector_laplacians() contains the nonlinearity,
     * such that the Kovasznay solution can be obtained as the solution to a
     * Stokes problem.
     */
    Kovasznay(const double Re, bool Stokes = false);

    virtual ~Kovasznay() override = default;

    virtual void
    vector_values(const std::vector<Point<2>>      &points,
                  std::vector<std::vector<double>> &values) const override;
    virtual void
    vector_gradients(
      const std::vector<Point<2>>            &points,
      std::vector<std::vector<Tensor<1, 2>>> &gradients) const override;
    virtual void
    vector_laplacians(const std::vector<Point<2>>      &points,
                      std::vector<std::vector<double>> &values) const override;

    /// The value of lambda.
    double
    lambda() const;

  private:
    const double Reynolds;
    double       lbda;
    double       p_average;
    const bool   stokes;
  };

} // namespace Functions

DEAL_II_NAMESPACE_CLOSE

#endif
