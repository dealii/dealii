//----------------------------  function_lib.h  ---------------------------
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
//----------------------------  function_lib.h  ---------------------------
#ifndef __deal2__function_lib_h
#define __deal2__function_lib_h


#include <base/function.h>



/**
 * d-quadratic pillow on the unit hypercube.
 *
 * This is a function for testing the implementation. It has zero Dirichlet
 * boundary values on the domain $(-1,1)^d$. In the inside, it is the
 * product of $x_i^2-1$.
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
				      * The value at a single point.
				      */
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;

				     /**
				      * Values at multiple points.
				      */
    virtual void value_list (const vector<Point<dim> > &points,
			     vector<double>            &values,
			     const unsigned int         component = 0) const;

				     /**
				      * Gradient at a single point.
				      */
    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
				    const unsigned int  component = 0) const;

				     /**
				      * Gradients at multiple points.
				      */
    virtual void gradient_list (const vector<Point<dim> > &points,
				vector<Tensor<1,dim> >    &gradients,
				const unsigned int         component = 0) const;

				     /**
				      * Laplacian at a single point.
				      */
    double laplacian (const Point<dim>   &p,
		      const unsigned int  component = 0) const;

				     /**
				      * Laplacian at multiple points.
				      */
    void laplacian_list (const vector<Point<dim> > &points,
			 vector<double>            &values,
			 const unsigned int         component = 0) const;
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
    virtual void value_list (const vector<Point<dim> > &points,
			     vector<double>            &values,
			     const unsigned int         component = 0) const;

				     /**
				      * Gradient at a single point.
				      */
    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
				    const unsigned int  component = 0) const;

				     /**
				      * Gradients at multiple points.
				      */
    virtual void gradient_list (const vector<Point<dim> > &points,
				vector<Tensor<1,dim> >    &gradients,
				const unsigned int         component = 0) const;

				     /**
				      * Laplacian at a single point.
				      */
    double laplacian(const Point<dim>   &p,
		     const unsigned int  component = 0) const;

				     /**
				      * Laplacian at multiple points.
				      */
    void laplacian_list(const vector<Point<dim> > &points,
			vector<double>            &values,
			const unsigned int         component = 0) const;
};



/**
 * Singularity on the L-shaped domain in 2D.
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
    virtual void value_list (const vector<Point<2> > &points,
			     vector<double>            &values,
			     const unsigned int         component = 0) const;

				     /**
				      * Gradient at a single point.
				      */
    virtual Tensor<1,2> gradient (const Point<2>   &p,
				    const unsigned int  component = 0) const;

				     /**
				      * Gradients at multiple points.
				      */
    virtual void gradient_list (const vector<Point<2> > &points,
				vector<Tensor<1,2> >    &gradients,
				const unsigned int         component = 0) const;

				     /**
				      * Laplacian at a single point.
				      */
    double laplacian(const Point<2>   &p,
		     const unsigned int  component = 0) const;

				     /**
				      * Laplacian at multiple points.
				      */
    void laplacian_list(const vector<Point<2> > &points,
			vector<double>            &values,
			const unsigned int         component = 0) const;
};



/**
 * Singularity on the slit domain in 2D.
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
    virtual void value_list (const vector<Point<2> > &points,
			     vector<double>            &values,
			     const unsigned int         component = 0) const;

				     /**
				      * Gradient at a single point.
				      */
    virtual Tensor<1,2> gradient (const Point<2>   &p,
				    const unsigned int  component = 0) const;

				     /**
				      * Gradients at multiple points.
				      */
    virtual void gradient_list (const vector<Point<2> > &points,
				vector<Tensor<1,2> >    &gradients,
				const unsigned int         component = 0) const;

				     /**
				      * Laplacian at a single point.
				      */
    double laplacian(const Point<2>   &p,
		     const unsigned int  component = 0) const;

				     /**
				      * Laplacian at multiple points.
				      */
    void laplacian_list(const vector<Point<2> > &points,
			vector<double>            &values,
			const unsigned int         component = 0) const;
};


#endif
