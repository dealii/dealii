//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__flow_function_h
#define __deal2__flow_function_h


#include <base/config.h>
#include <base/function.h>
#include <base/point.h>

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
      virtual void vector_values (const std::vector<Point<dim> >& points,
				  std::vector<std::vector<double> >& values) const = 0;
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
					*/
      virtual void vector_laplacians (const std::vector<Point<dim> > &points,
				      std::vector<std::vector<double> >   &values) const = 0;

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
      
      unsigned int memory_consumption () const;

    protected:
				       /**
					* Mean value of the pressure
					* to be added by derived
					* classes.
					*/
      double mean_pressure;

    private:  
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
 * stretches along the <i>x</i>-axis and has radius #radius. The
 * #Reynolds number is used to scale the pressure properly for a
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
      
      virtual void vector_values (const std::vector<Point<dim> >& points,
				  std::vector<std::vector<double> >& values) const;
      virtual void vector_gradients (const std::vector<Point<dim> >& points,
				     std::vector<std::vector<Tensor<1,dim> > >& gradients) const;
      virtual void vector_laplacians (const std::vector<Point<dim> > &points,
				      std::vector<std::vector<double> >   &values) const;

    private:
      const double radius;
      const double Reynolds;
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
      
      virtual void vector_values (const std::vector<Point<2> >& points,
				  std::vector<std::vector<double> >& values) const;
      virtual void vector_gradients (const std::vector<Point<2> >& points,
				     std::vector<std::vector<Tensor<1,2> > >& gradients) const;
      virtual void vector_laplacians (const std::vector<Point<2> > &points,
				      std::vector<std::vector<double> >   &values) const;

    private:
				       /// The auxiliary function Psi.
      double Psi(double phi) const;
				       /// The derivative of #Psi
      double Psi_1(double phi) const;
				       /// The 3rd derivative of #Psi
      double Psi_3(double phi) const;
				       /// The angle of the reentrant corner
      const double omega;
				       /// The exponent of the radius
      static const double lambda;
				       /// Cosine of #lambda times #omega
      const double coslo;
  };
  
}

DEAL_II_NAMESPACE_CLOSE

#endif
