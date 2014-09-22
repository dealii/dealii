// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2013 by the deal.II authors
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

#ifndef __deal2__flow_function_h
#define __deal2__flow_function_h


#include <deal.II/base/config.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/thread_management.h>

DEAL_II_NAMESPACE_OPEN

namespace Functions
{
  /**
   * Base class for analytic solutions to incompressible flow problems.
   *
   * Additional to the Function interface, this function provides for an
   * offset of the pressure: if the pressure of the computed solution
   * has an integral mean value different from zero, this value can be
   * given to pressure_adjustment() in order to compute correct pressure
   * errors.
   *
   * @note Derived classes should implement pressures with integral mean
   * value zero always.
   *
   * @note Thread safety: Some of the functions make use of internal data to compute
   * values. Therefore, every thread should obtain its own object of
   * derived classes.
   *
   * @ingroup functions
   * @author Guido Kanschat, 2007
   */
  template <int dim>
  class FlowFunction : public Function<dim>
  {
  public:
    /**
     * Constructor, setting up some
     * internal data structures.
     */
    FlowFunction();

    /**
     * Virtual destructor.
     */
    virtual ~FlowFunction();

    /**
     * Store an adjustment for the
     * pressure function, such that
     * its mean value is
     * <tt>p</tt>.
     */
    void pressure_adjustment(double p);

    /**
     * Values in a structure more
     * suitable for vector valued
     * functions. The outer vector
     * is indexed by solution
     * component, the inner by
     * quadrature point.
     */
    virtual void vector_values (const std::vector<Point<dim> > &points,
                                std::vector<std::vector<double> > &values) const = 0;
    /**
     * Gradients in a structure more
     * suitable for vector valued
     * functions. The outer vector
     * is indexed by solution
     * component, the inner by
     * quadrature point.
     */
    virtual void vector_gradients (const std::vector<Point<dim> >            &points,
                                   std::vector<std::vector<Tensor<1,dim> > > &gradients) const = 0;
    /**
     * Force terms in a structure more
     * suitable for vector valued
     * functions. The outer vector
     * is indexed by solution
     * component, the inner by
     * quadrature point.
     *
     * @warning This is not the
     * true Laplacian, but the
     * force term to be used as
     * right hand side in Stokes'
     * equations
     */
    virtual void vector_laplacians (const std::vector<Point<dim> > &points,
                                    std::vector<std::vector<double> >   &values) const = 0;

    virtual void vector_value (const Point<dim> &points, Vector<double> &value) const;
    virtual double value (const Point<dim> &points, const unsigned int component) const;
    virtual void vector_value_list (const std::vector<Point<dim> > &points,
                                    std::vector<Vector<double> >   &values) const;
    virtual void vector_gradient_list (const std::vector<Point<dim> >            &points,
                                       std::vector<std::vector<Tensor<1,dim> > > &gradients) const;
    /**
     * The force term in the
     * momentum equation.
     */
    virtual void vector_laplacian_list (const std::vector<Point<dim> > &points,
                                        std::vector<Vector<double> >   &values) const;

    std::size_t memory_consumption () const;

  protected:
    /**
     * Mean value of the pressure
     * to be added by derived
     * classes.
     */
    double mean_pressure;

  private:

    /**
     * A mutex that guards the
     * following scratch arrays.
     */
    mutable Threads::Mutex mutex;

    /**
     * Auxiliary values for the usual
     * Function interface.
     */
    mutable std::vector<std::vector<double> > aux_values;

    /**
     * Auxiliary values for the usual
     * Function interface.
     */
    mutable std::vector<std::vector<Tensor<1,dim> > > aux_gradients;
  };

  /**
   * Laminar pipe flow in two and three dimensions. The channel
   * stretches along the <i>x</i>-axis and has radius @p radius. The
   * @p Reynolds number is used to scale the pressure properly for a
   * Navier-Stokes problem.
   *
   * @ingroup functions
   * @author Guido Kanschat, 2007
   */
  template <int dim>
  class PoisseuilleFlow : public FlowFunction<dim>
  {
  public:
    /**
     * Construct an object for the
     * given channel radius
     * <tt>r</tt> and the Reynolds
     * number <tt>Re</tt>.
     */
    PoisseuilleFlow<dim> (const double r,
                          const double Re);
    virtual ~PoisseuilleFlow();

    virtual void vector_values (const std::vector<Point<dim> > &points,
                                std::vector<std::vector<double> > &values) const;
    virtual void vector_gradients (const std::vector<Point<dim> > &points,
                                   std::vector<std::vector<Tensor<1,dim> > > &gradients) const;
    virtual void vector_laplacians (const std::vector<Point<dim> > &points,
                                    std::vector<std::vector<double> >   &values) const;

  private:
    const double radius;
    const double Reynolds;
  };


  /**
   * Artificial divergence free function with homogeneous boundary
   * conditions on the cube [-1,1]<sup>dim</sup>.
   *
   * The function in 2D is
   * @f[
   * \left(\begin{array}{c}u\\v\\p\end{array}\right)
   * \left(\begin{array}{c}\cos^2x \sin y\cos y\\-\sin x\cos x\cos^2y\\
   * \sin x\cos x\sin y\cos y\end{array}\right)
   * @f]
   * @ingroup functions
   * @author Guido Kanschat, 2007
   */
  template <int dim>
  class StokesCosine :
    public FlowFunction<dim>
  {
  public:
    /**
     * Constructor setting the
     * Reynolds number required for
     * pressure computation and
     * scaling of the right hand side.
     */
    StokesCosine (const double viscosity = 1., const double reaction = 0.);
    /**
     * Change the viscosity and the
     * reaction parameter.
     */
    void set_parameters (const double viscosity, const double reaction);
    virtual ~StokesCosine();

    virtual void vector_values (const std::vector<Point<dim> > &points,
                                std::vector<std::vector<double> > &values) const;
    virtual void vector_gradients (const std::vector<Point<dim> > &points,
                                   std::vector<std::vector<Tensor<1,dim> > > &gradients) const;
    virtual void vector_laplacians (const std::vector<Point<dim> > &points,
                                    std::vector<std::vector<double> >   &values) const;

  private:
    /// The viscosity
    double viscosity;
    /// The reaction parameter
    double reaction;
  };


  /**
   * The solution to Stokes' equations on an L-shaped domain.
   *
   * Taken from Houston, Sch&ouml;tzau, Wihler, proceeding ENUMATH 2003.
   *
   * @ingroup functions
   * @author Guido Kanschat, 2007
   */
  class StokesLSingularity : public FlowFunction<2>
  {
  public:
    /// Constructor setting upsome data.
    StokesLSingularity();

    virtual void vector_values (const std::vector<Point<2> > &points,
                                std::vector<std::vector<double> > &values) const;
    virtual void vector_gradients (const std::vector<Point<2> > &points,
                                   std::vector<std::vector<Tensor<1,2> > > &gradients) const;
    virtual void vector_laplacians (const std::vector<Point<2> > &points,
                                    std::vector<std::vector<double> >   &values) const;
  private:
    /// The auxiliary function Psi.
    double Psi(double phi) const;
    /// The derivative of Psi()
    double Psi_1(double phi) const;
    /// The 2nd derivative of Psi()
    double Psi_2(double phi) const;
    /// The 3rd derivative of Psi()
    double Psi_3(double phi) const;
    /// The 4th derivative of Psi()
    double Psi_4(double phi) const;
    /// The angle of the reentrant corner
    const double omega;
    /// The exponent of the radius
    static const double lambda;
    /// Cosine of lambda times omega
    const double coslo;
    /// Auxiliary variable 1+lambda
    const double lp;
    /// Auxiliary variable 1-lambda
    const double lm;
  };

  /**
   * Flow solution in 2D by Kovasznay (1947).
   *
   * This function is valid on the half plane right of the line
   * <i>x=1/2</i>.
   *
   * @ingroup functions
   * @author Guido Kanschat, 2007
   */
  class Kovasznay : public FlowFunction<2>
  {
  public:
    /**
     * Construct an object for the
     * give Reynolds number
     * <tt>Re</tt>. If the
     * parameter <tt>Stokes</tt> is
     * true, the right hand side of
     * the momentum equation
     * returned by
     * vector_laplacians() contains
     * the nonlinearity, such that
     * the Kovasznay solution can
     * be obtained as the solution
     * to a Stokes problem.
     */
    Kovasznay (const double Re, bool Stokes = false);
    virtual ~Kovasznay();

    virtual void vector_values (const std::vector<Point<2> > &points,
                                std::vector<std::vector<double> > &values) const;
    virtual void vector_gradients (const std::vector<Point<2> > &points,
                                   std::vector<std::vector<Tensor<1,2> > > &gradients) const;
    virtual void vector_laplacians (const std::vector<Point<2> > &points,
                                    std::vector<std::vector<double> >   &values) const;

    /// The value of lambda.
    double lambda () const;
  private:
    const double Reynolds;
    double lbda;
    double p_average;
    const bool stokes;
  };

}

DEAL_II_NAMESPACE_CLOSE

#endif
