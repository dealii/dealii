/*----------------------------   function_lib.h     ---------------------------*/
/*      $Id$               */
/* Copyright Guido Kanschat, University of Minnesota, 1999                     */
#ifndef __function_lib_H
#define __function_lib_H
/*----------------------------   function_lib.h     ---------------------------*/


#include <base/function.h>

/**
 * n-quadratic pillow on the unit square.
 *
 * This is a function for testing the implementation. It has zero
 * boundary values on the domain $(-1,1)^d$. In the inside, it is the
 * product of $x_i^2-1$.
 *
 * Together with the function, its derivatives and Laplacian are defined.
 *
 * @author: Guido Kanschat, 1999
 */
template<int dim>
class PillowFunction : Function<dim>
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



/*----------------------------   function_lib.h     ---------------------------*/
/* end of #ifndef __function_lib_H */
#endif
/*----------------------------   function_lib.h     ---------------------------*/
