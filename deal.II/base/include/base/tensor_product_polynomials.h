//----------------------  tensor_product_polynomials.h  -------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------  tensor_product_polynomials.h  -------------
#ifndef __deal2__tensor_product_polynomials_h
#define __deal2__tensor_product_polynomials_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/tensor.h>
#include <base/point.h>
#include <base/polynomial.h>
#include <base/smartpointer.h>

#include <vector>


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
				      * one-dimensional polynomials
				      * and will be copied into the
				      * member variable @p{polynomials}.
				      */
    template <class Pol>
    TensorProductPolynomials(const typename std::vector<Pol> &pols);

				     /**
				      * Computes the value and the
				      * first and second derivatives
				      * of each tensor product
				      * polynomial at @p{unit_point}.
				      *
				      * The size of the vectors must
				      * either be equal @p{0} or equal
				      * @p{n_tensor_pols}.  In the
				      * first case, the function will
				      * not compute these values.
				      *
				      * If you need values or
				      * derivatives of all tensor
				      * product polynomials then use
				      * this function, rather than
				      * using any of the
				      * @p{compute_value},
				      * @p{compute_grad} or
				      * @p{compute_grad_grad}
				      * functions, see below, in a
				      * loop over all tensor product
				      * polynomials.
				      */
    void compute(const Point<dim>                     &unit_point,
		 std::vector<double>                  &values,
		 typename std::vector<Tensor<1,dim> > &grads,
		 typename std::vector<Tensor<2,dim> > &grad_grads) const;
    
				     /**
				      * Computes the value of the
				      * @p{i}th tensor product
				      * polynomial at
				      * @p{unit_point}. Here @p{i} is
				      * given in tensor product
				      * numbering.
				      *
				      * Note, that using this function
				      * within a loop over all tensor
				      * product polynomials is not
				      * efficient, because then each
				      * point value of the underlying
				      * (one-dimensional) polynomials
				      * is (unnecessarily) computed
				      * several times.  Instead use
				      * the @p{compute} function, see
				      * above, with
				      * @p{values.size()==n_tensor_pols}
				      * to get the point values of all
				      * tensor polynomials all at once
				      * and in a much more efficient
				      * way.
				      */
    double compute_value (const unsigned int i,
			  const Point<dim> &p) const;

				     /**
				      * Computes the grad of the
				      * @p{i}th tensor product
				      * polynomial at
				      * @p{unit_point}. Here @p{i} is
				      * given in tensor product
				      * numbering.
				      *
				      * Note, that using this function
				      * within a loop over all tensor
				      * product polynomials is not
				      * efficient, because then each
				      * derivative value of the
				      * underlying (one-dimensional)
				      * polynomials is (unnecessarily)
				      * computed several times.
				      * Instead use the @p{compute}
				      * function, see above, with
				      * @p{grads.size()==n_tensor_pols}
				      * to get the point value of all
				      * tensor polynomials all at once
				      * and in a much more efficient
				      * way.
				      */
    Tensor<1,dim> compute_grad (const unsigned int i,
				const Point<dim> &p) const;

				     /**
				      * Computes the second
				      * derivative (grad_grad) of the
				      * @p{i}th tensor product
				      * polynomial at
				      * @p{unit_point}. Here @p{i} is
				      * given in tensor product
				      * numbering.
				      *
				      * Note, that using this function
				      * within a loop over all tensor
				      * product polynomials is not
				      * efficient, because then each
				      * derivative value of the
				      * underlying (one-dimensional)
				      * polynomials is (unnecessarily)
				      * computed several times.
				      * Instead use the @p{compute}
				      * function, see above, with
				      * @p{grad_grads.size()==n_tensor_pols}
				      * to get the point value of all
				      * tensor polynomials all at once
				      * and in a much more efficient
				      * way.
				      */
    Tensor<2,dim> compute_grad_grad(const unsigned int i,
				    const Point<dim> &p) const;

				     /**
				      * Returns the number of tensor
				      * product polynomials. For $n$
				      * 1d polynomials this is $n^dim$.
				      */
    unsigned int n_tensor_product_polynomials() const;

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
    const unsigned int n_tensor_pols;

				     /**
				      * @p{n_pols_to[n]=polynomials.size()^n}
				      * Filled by the constructor.
				      *
				      * For internal use only. 
				      */
    std::vector<unsigned int> n_pols_to;
    
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
TensorProductPolynomials<dim>::TensorProductPolynomials(
  const typename std::vector<Pol> &pols):
		polynomials (pols.begin(), pols.end()),
		n_tensor_pols(power(pols.size(), dim)),
		n_pols_to(dim+1)
{
  const unsigned int n_pols=polynomials.size();

  n_pols_to[0]=1;
  for (unsigned int i=0; i<dim; ++i)
    n_pols_to[i+1]=n_pols_to[i]*n_pols;
  Assert(n_pols_to[dim]==n_tensor_pols, ExcInternalError());
}



#endif
