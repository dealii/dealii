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
 * Representation of the space of polynomials of degree at most n in
 * higher dimensions.
 *
 * Given a vector of @{n} one-dimensional polynomials @{P0} to @{Pn},
 * where @{Pi} has degree @p{i}, this class generates all polynomials
 * of the form @p{ Pijk(x,y,z) = Pi(x)Pj(y)Pk(z)}, where the sum of
 * @p{i}, @p{j} and @p{k} is less than or equal @p{n}.
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
				      * member variable
				      * @p{polynomials}.  The static
				      * type of the template argument
				      * @p{pols} needs to be
				      * convertible to
				      * @p{Polynomial<double>},
				      * i.e. should usually be a
				      * derived class of
				      * @p{Polynomial<double>}.
				      */
    template <class Pol>
    PolynomialSpace(const std::vector<Pol> &pols);

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
    void compute (const Point<dim>                     &unit_point,
		  std::vector<double>                  &values,
		  std::vector<Tensor<1,dim> > &grads,
		  std::vector<Tensor<2,dim> > &grad_grads) const;
    
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
				      * Return the number of
				      * polynomials spanning the space
				      * represented by this
				      * class. Here, if @p{N} is the
				      * number of one-dimensional
				      * polynomials given, then the
				      * result of this function is
				      * @p{N} in 1d, @p{N(N+1)/2} in
				      * 2d, and @p{N(N+1)(N+2)/6 in
				      * 3d.
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
    const std::vector<Polynomials::Polynomial<double> > polynomials;

				     /**
				      * Store the precomputed value
				      * which the @p{n()} function
				      * returns.
				      */
    const unsigned int n_pols;

				     /**
				      * Compute numbers in x, y and z
				      * direction. Given an index
				      * @p{n} in the d-dimensional
				      * polynomial space, compute the
				      * indices i,j,k such that
				      * @p{p_n(x,y,z) =
				      * p_i(x)p_j(y)p_k(z)}.
				      */
    void compute_index(unsigned int n,
		       unsigned int& nx,
		       unsigned int& ny,
		       unsigned int& nz) const;
    
				     /**
				      * Static function used in the
				      * constructor to compute the
				      * number of polynomials.
				      */
    static unsigned int compute_n_pols (const unsigned int n);
};



template <int dim>
template <class Pol>
PolynomialSpace<dim>::
PolynomialSpace (const std::vector<Pol> &pols)
		:
		polynomials (pols.begin(), pols.end()),
		n_pols (compute_n_pols(polynomials.size()))
{}



#endif
