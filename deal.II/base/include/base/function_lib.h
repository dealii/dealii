//----------------------------  function_lib.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  function_lib.h  ---------------------------
#ifndef __deal2__function_lib_h
#define __deal2__function_lib_h


#include <base/function.h>



/**
 * Namespace implementing some concrete classes derived from the
 * @ref{Function} class that describe actual functions. This is rather
 * a collection of classes that we have needed for our own programs
 * once and thought they might be useful to others as well at some
 * point.
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
 * @author: Guido Kanschat, 1999
 */
  template<int dim>
  class SquareFunction : public Function<dim>
  {
    public:
				       /**
					* Function value at one point.
					*/
      virtual double value (const Point<dim>   &p,
			    const unsigned int  component = 0) const;
      
				       /**
					* Function values at multiple points.
					*/
      virtual void value_list (const typename std::vector<Point<dim> > &points,
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
      virtual void gradient_list (const typename std::vector<Point<dim> > &points,
				  typename std::vector<Tensor<1,dim> >    &gradients,
				  const unsigned int              component = 0) const;
      
				       /**
					* Laplacian of the function at one point.
					*/
      virtual double laplacian (const Point<dim>   &p,
				const unsigned int  component = 0) const;
      
				       /**
					* Laplacian of the function at multiple points.
					*/
      virtual void laplacian_list (const typename std::vector<Point<dim> > &points,
				   std::vector<double>            &values,
				   const unsigned int              component = 0) const;
  };
  
  
  
/**
 * The function @p{xy}. This function serves as an example for
 * a vanishing Laplacian.
 *
 * @author: Guido Kanschat, 2000
 */
  template<int dim>
  class Q1WedgeFunction : public Function<dim>
  {
    public:
				       /**
					* Function value at one point.
					*/
      virtual double value (const Point<dim>   &p,
			    const unsigned int  component = 0) const;
      
				       /**
					* Function values at multiple points.
					*/
      virtual void value_list (const typename std::vector<Point<dim> > &points,
			       std::vector<double>            &values,
			       const unsigned int              component = 0) const;
      
				       /**
					* Gradient at one point.
					*/
      virtual Tensor<1,dim> gradient (const Point<dim>   &p,
				      const unsigned int       component = 0) const;
      
				       /**
					  Gradients at multiple points.
				       */
      virtual void gradient_list (const typename std::vector<Point<dim> > &points,
				  typename std::vector<Tensor<1,dim> >    &gradients,
				  const unsigned int              component = 0) const;
      
				       /**
					* Laplacian of the function at one point.
					*/
      virtual double laplacian (const Point<dim>   &p,
				const unsigned int  component = 0) const;
      
				       /**
					* Laplacian of the function at multiple points.
					*/
      virtual void laplacian_list (const typename std::vector<Point<dim> > &points,
				   std::vector<double>            &values,
				   const unsigned int              component = 0) const;
  };
  
  
  
/**
 * d-quadratic pillow on the unit hypercube.
 *
 * This is a function for testing the implementation. It has zero Dirichlet
 * boundary values on the domain $(-1,1)^d$. In the inside, it is the
 * product of $x_i^2-1$.
 *
 * Providing a non-zero argument to the constructor, the whole function
 * can be offset by a constant.
 *
 * Together with the function, its derivatives and Laplacian are defined.
 *
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
      virtual void value_list (const typename std::vector<Point<dim> > &points,
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
      virtual void gradient_list (const typename std::vector<Point<dim> > &points,
				  typename std::vector<Tensor<1,dim> >    &gradients,
				  const unsigned int              component = 0) const;
      
				       /**
					* Laplacian at a single point.
					*/
      virtual double laplacian (const Point<dim>   &p,
				const unsigned int  component = 0) const;
      
				       /**
					* Laplacian at multiple points.
					*/
      virtual void laplacian_list (const typename std::vector<Point<dim> > &points,
				   std::vector<double>            &values,
				   const unsigned int              component = 0) const;
    private:
      const double offset;
  };
  
  
  
/**
 * Cosine-shaped pillow function.
 * This is another function with zero boundary values on $[-1,1]^d$. In the interior
 * it is the product of $\cos(\pi/2 x_i)$.
 * @author Guido Kanschat, 1999
 */
  template<int dim>
  class CosineFunction : public Function<dim>
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
      virtual void value_list (const typename std::vector<Point<dim> > &points,
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
      virtual void gradient_list (const typename std::vector<Point<dim> > &points,
				  typename std::vector<Tensor<1,dim> >    &gradients,
				  const unsigned int              component = 0) const;
      
				       /**
					* Laplacian at a single point.
					*/
      virtual double laplacian (const Point<dim>   &p,
				const unsigned int  component = 0) const;
      
				       /**
					* Laplacian at multiple points.
					*/
      virtual void laplacian_list (const typename std::vector<Point<dim> > &points,
				   std::vector<double>            &values,
				   const unsigned int              component = 0) const;
      
				       /**
					* Gradient at a single point.
					*/
      virtual Tensor<2,dim> hessian (const Point<dim>   &p,
				     const unsigned int  component = 0) const;
      
				       /**
					* Gradients at multiple points.
					*/
      virtual void hessian_list (const typename std::vector<Point<dim> > &points,
				 typename std::vector<Tensor<2,dim> >    &hessians,
				 const unsigned int              component = 0) const;
  };
  
  
  
/**
 * Product of exponential functions in each coordinate direction.
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
      virtual void value_list (const typename std::vector<Point<dim> > &points,
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
      virtual void gradient_list (const typename std::vector<Point<dim> > &points,
				  typename std::vector<Tensor<1,dim> >    &gradients,
				  const unsigned int              component = 0) const;
      
				       /**
					* Laplacian at a single point.
					*/
      virtual double laplacian (const Point<dim>   &p,
				const unsigned int  component = 0) const;
      
				       /**
					* Laplacian at multiple points.
					*/
      virtual void laplacian_list (const typename std::vector<Point<dim> > &points,
				   std::vector<double>            &values,
				   const unsigned int              component = 0) const;
  };
  
  
  
/**
 * Singularity on the L-shaped domain in 2D.
 *
 * Caveat: derivatives of this function are not implemented!
 *
 * @author Guido Kanschat, 1999
 */
  class LSingularityFunction : public Function<2>
  {
    public:
				       /**
					* The value at a single point.
					*/
      virtual double value (const Point<2>   &p,
			    const unsigned int  component = 0) const;
      
				       /**
					* Values at multiple points.
					*/
      virtual void value_list (const std::vector<Point<2> > &points,
			       std::vector<double>            &values,
			       const unsigned int              component = 0) const;
      
				       /**
					* Gradient at a single point.
					*/
      virtual Tensor<1,2> gradient (const Point<2>     &p,
				    const unsigned int  component = 0) const;
      
				       /**
					* Gradients at multiple points.
					*/
      virtual void gradient_list (const std::vector<Point<2> > &points,
				  std::vector<Tensor<1,2> >    &gradients,
				  const unsigned int            component = 0) const;
      
				       /**
					* Laplacian at a single point.
					*/
      virtual double laplacian (const Point<2>   &p,
				const unsigned int  component = 0) const;
      
				       /**
					* Laplacian at multiple points.
					*/
      virtual void laplacian_list (const std::vector<Point<2> > &points,
				   std::vector<double>          &values,
				   const unsigned int            component = 0) const;
  };
  
  
  
/**
 * Singularity on the slit domain in 2D.
 *
 * Caveat: derivatives of this function are not implemented!
 *
 * @author Guido Kanschat, 1999
 */
  class SlitSingularityFunction : public Function<2>
  {
    public:
				       /**
					* The value at a single point.
					*/
      virtual double value (const Point<2>   &p,
			    const unsigned int  component = 0) const;
      
				       /**
					* Values at multiple points.
					*/
      virtual void value_list (const std::vector<Point<2> > &points,
			       std::vector<double>            &values,
			       const unsigned int              component = 0) const;
      
				       /**
					* Gradient at a single point.
					*/
      virtual Tensor<1,2> gradient (const Point<2>   &p,
				    const unsigned int  component = 0) const;
      
				       /**
					* Gradients at multiple points.
					*/
      virtual void gradient_list (const std::vector<Point<2> > &points,
				  std::vector<Tensor<1,2> >    &gradients,
				  const unsigned int            component = 0) const;
      
				       /**
					* Laplacian at a single point.
					*/
      virtual double laplacian (const Point<2>   &p,
				const unsigned int  component = 0) const;
      
				       /**
					* Laplacian at multiple points.
					*/
      virtual void laplacian_list (const std::vector<Point<2> > &points,
				   std::vector<double>          &values,
				   const unsigned int            component = 0) const;
  };
  
  
/**
 * A jump in x-direction transported into some direction.
 *
 * If the advection is parallel to the y-axis, the function is
 * @p{-atan(sx)}, where @p{s} is the steepness parameter provided in
 * the constructor.
 *
 * For different advection directions, this function will be turned in
 * the parameter space.
 *
 * Together with the function, its derivatives and Laplacian are defined.
 *
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
					* Function values at multiple points.
					*/
      virtual void value_list (const typename std::vector<Point<dim> > &points,
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
      virtual void gradient_list (const typename std::vector<Point<dim> > &points,
				  typename std::vector<Tensor<1,dim> >    &gradients,
				  const unsigned int              component = 0) const;
      
				       /**
					* Laplacian of the function at one point.
					*/
      virtual double laplacian (const Point<dim>   &p,
				const unsigned int  component = 0) const;
      
				       /**
					* Laplacian of the function at multiple points.
					*/
      virtual void laplacian_list (const typename std::vector<Point<dim> > &points,
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
					* STL @p{std::map} type with a
					* certain number of
					* elements?), this is only
					* an estimate. however often
					* quite close to the true
					* value.
					*/
      unsigned int memory_consumption () const;
      
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
					* Sine of @p{angle}.
					*/
      double sine;
      
				       /**
					* Cosine of @p{angle}.
					*/
      double cosine;
  };
  
  
  
/**
 * Given a wavenumber vector generate a cosine function. The
 * wavenumber coefficient is given as a @p{d}-dimensional point @p{k}
 * in Fourier space, and the function is then recovered as @p{f(x) =
 * \prod_i cos(k_i x_i) = Re(\exp(i k.x))}.
 *
 * The class has its name from the fact that it resembles one
 * component of a Fourier cosine decomposition.
 *
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
					* given component at point @p{p}.
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
 * wavenumber coefficient is given as a @p{d}-dimensional point @p{k}
 * in Fourier space, and the function is then recovered as @p{f(x) =
 * \prod_i sin(k_i x_i) = Im(\exp(i k.x))}.
 *
 * The class has its name from the fact that it resembles one
 * component of a Fourier sine decomposition.
 *
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
					* given component at point @p{p}.
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
 * @p{d}-dimensional point @p{k} in Fourier space, and the entire
 * function is then recovered as
 * @p{f(x) = \sum_j w_j \prod_i sin(k_i x_i) = Im(\sum_j w_j \exp(i k.x))}.
 *
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
					* given component at point @p{p}.
					*/
      virtual double laplacian (const Point<dim>   &p,
				const unsigned int  component = 0) const;

				       /**
					* Exception
					*/
      DeclException0 (ExcInvalidArraySize);
      
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
 * @p{d}-dimensional point @p{k} in Fourier space, and the entire
 * function is then recovered as
 * @p{f(x) = \sum_j w_j \prod_i cos(k_i x_i) = Re(\sum_j w_j \exp(i k.x))}.
 *
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
					* given component at point @p{p}.
					*/
      virtual double laplacian (const Point<dim>   &p,
				const unsigned int  component = 0) const;

				       /**
					* Exception
					*/
      DeclException0 (ExcInvalidArraySize);
      
    private:
				       /**
					* Stored Fourier coefficients
					* and weights.
					*/
      const std::vector<Point<dim> > fourier_coefficients;
      const std::vector<double>      weights;
  };
  
  
};


#endif
