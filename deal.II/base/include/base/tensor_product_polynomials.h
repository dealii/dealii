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


#include <vector.h>


/**
 * Tensor product of given polynomials.
 *
 * @author Ralf Hartmann, 2000
 */
template <int dim>
class TensorProductPolynomials
{
  public:
				     /**
				      * Constructor. @p{pols} is a
				      * vector of 1d polynomial
				      * functions. For polynomials of
				      * order @p{p} there should be
				      * @p{p+1} polynomials.
				      */
    TensorProductPolynomials(const vector<Polynomial> &pols);

				     /**
				      * Calculates the shape values
				      * and shape grads at the
				      * @p{unit_point}.
				      */
    void shape_values_and_grads(const Point<dim> &unit_point,
				vector<double> > &values,
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
				      * TODO: Implement use of
				      * SmartPointer later.
				      *
				      * Pointer to the @p{polynomials}
				      * given to the constructor.
				      */
    const vector<Polynomial> *polynomials;
}




#endif
