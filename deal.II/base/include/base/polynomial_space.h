//----------------------  polynomials.h  -------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------  polynomials.h  -------------
#ifndef __deal2__polynomial_space_h
#define __deal2__polynomial_space_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/tensor.h>
#include <base/point.h>
#include <base/polynomial.h>
#include <base/smartpointer.h>

#include <vector>


/**
 * Polynomial space of degree at most n in higher dimensions.
 *
 * Given a vector of @{n} one-dimensional polynomials @{P0} to @{Pn},
 * where @{Pi} has degree @p{i}, this class generates all polynomials
 * the form @p{ Pijk(x,y,z) = Pi(x)Pj(y)Pk(z)}, where the sum of
 * @p{i}, @p{j} and @p{k} is below/equal @p{n}.
 *
 * @author Guido Kanschat, 2002
 */
template <int dim>
class PolynomialSpace
{
  public:
				     /**
				      * Constructor. @p{pols} is a
				      * vector of pointers to
				      * one-dimensional polynomials
				      * and will be copied into the
				      * member variable @p{polynomials}.
				      */
    template <class Pol>
    PolynomialSpace(const typename std::vector<Pol> &pols);

				     /**
				      * Computes the value and the
				      * first and second derivatives
				      * of each polynomial at
				      * @p{unit_point}.
				      *
				      * The size of the vectors must
				      * either be equal @p{0} or equal
				      * @p{n()}.  In the
				      * first case, the function will
				      * not compute these values.
				      *
				      * If you need values or
				      * derivatives of all polynomials
				      * then use this function, rather
				      * than using any of the
				      * @p{compute_value},
				      * @p{compute_grad} or
				      * @p{compute_grad_grad}
				      * functions, see below, in a
				      * loop over all polynomials.
				      */
    void compute(const Point<dim>                     &unit_point,
		 std::vector<double>                  &values,
		 typename std::vector<Tensor<1,dim> > &grads,
		 typename std::vector<Tensor<2,dim> > &grad_grads) const;
    
				     /**
				      * Computes the value of the
				      * @p{i}th polynomial at
				      * @p{unit_point}.
				      *
				      * Consider using @p{compute} instead.
				      */
    double compute_value (const unsigned int i,
			  const Point<dim> &p) const;

				     /**
				      * Computes the gradient of the
				      * @p{i}th polynomial at
				      * @p{unit_point}.
				      *
				      * Consider using @p{compute} instead.
				      */
    Tensor<1,dim> compute_grad (const unsigned int i,
				const Point<dim> &p) const;

				     /**
				      * Computes the second derivative
				      * (grad_grad) of the @p{i}th
				      * polynomial at
				      * @p{unit_point}.
				      *
				      * Consider using @p{compute} instead.
				      */
    Tensor<2,dim> compute_grad_grad(const unsigned int i,
				    const Point<dim> &p) const;

				     /**
				      * Returns the number of tensor
				      * product polynomials. For $n$
				      * 1d polynomials this is $n^dim$.
				      */
    unsigned int n() const;

				     /**
				      * Exception.
				      */
    DeclException3 (ExcDimensionMismatch2,
		    int, int, int,
		    << "Dimension " << arg1 << " not equal to " << arg2 << " nor to " << arg3);

	    
  private:
				     /**
				      * Copy of the vector @p{pols} of
				      * polynomials given to the
				      * constructor.
				      */
    std::vector<Polynomial<double> > polynomials;

				     /**
				      * Number of tensor product
				      * polynomials. For $n$ 1d
				      * polynomials this is $n^dim$.
				      */
    unsigned int n_pols;
    
				     /**
				      * Computes @p{x} to the power of
				      * @p{y} for unsigned int @p{x}
				      * and @p{y}. It is a private
				      * function as it is only used in
				      * this class.
				      */
    static unsigned int power(const unsigned int x, const unsigned int y);
};



template <int dim>
template <class Pol>
PolynomialSpace<dim>::PolynomialSpace(
  const typename std::vector<Pol> &pols):
		polynomials (pols.begin(), pols.end())
{
  const unsigned int n=polynomials.size();

  n_pols = n;
  for (unsigned int i=1;i<dim;++i)
    {
      n_pols *= (n+i);
      n_pols /= (i+1);
    }
}



#endif
