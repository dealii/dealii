//----------------------  tensor_product_polynomials.h  -------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------  tensor_product_polynomials.h  -------------
#ifndef __deal2__tensor_product_polynomials_h
#define __deal2__tensor_product_polynomials_h


#include <base/exceptions.h>
#include <base/tensor.h>
#include <base/point.h>
#include <base/polynomial.h>
#include <base/smartpointer.h>

#include <vector.h>


/**
 * Tensor product of given polynomials.
 *
 * Given a vector of @{n} one-dimensional polynomials @{P1} to @{Pn},
 * this class generates @p{n} to the power of @p{dim} polynomials of
 * the form @p{ Qijk(x,y,z) = Pi(x)Pj(y)Pk(z)}.
 *
 * @author Ralf Hartmann, 2000, documentation Guido Kanschat
 */
template <int dim>
class TensorProductPolynomials
{
  public:
				     /**
				      * Constructor. @p{pols} is a
				      * vector of pointers to
				      * one-dimensional polynomials.
				      */
    TensorProductPolynomials(const vector<SmartPointer<Polynomial> > &pols);

				     /**
				      * Calculates the polynomials
				      * and their derivatives at
				      * @p{unit_point}.
				      *
				      * The vectors must either have
				      * length @p{0} or number of
				      * polynomials. In the first
				      * case, the function will not
				      * compute these values.
				      */
    void compute (const Point<dim> &unit_point,
		  vector<double> &values,
		  vector<Tensor<1,dim> > &grads,
		  vector<Tensor<2,dim> > &grad_grads) const;

				     /**
				      * Exception.
				      */
    DeclException3 (ExcDimensionMismatch2,
		    int, int, int,
		    << "Dimension " << arg1 << " not equal to " << arg2 << " nor to " << arg3);

	    
  private:
				     /**
				      * Pointer to the @p{polynomials}
				      * given to the constructor.
				      */
    vector<SmartPointer<Polynomial> > polynomials;
};





#endif
