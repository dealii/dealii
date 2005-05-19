//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__polynomials_raviart_thomas_h
#define __deal2__polynomials_raviart_thomas_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/tensor.h>
#include <base/point.h>
#include <base/polynomial.h>
#include <base/polynomial_space.h>
#include <base/table.h>

#include <vector>


/**
 * This class implements the <i>H<sup>div</sup></i>-conforming,
 * vector-valued Raviart-Thomas polynomials as described in the
 * book by Brezzi and Fortin.
 *
 * The Raviart-Thomas polynomials are constructed such that the
 * divergence is in the tensor product polynomial space
 * <i>Q<sub>k</sub></i>. Therefore, the polynomial order of each
 * component must be one order higher in the corresponding direction,
 * yielding the polynomial spaces <i>(Q<sub>k+1,k</sub>,
 * Q<sub>k,k+1</sub>)</i> and <i>(Q<sub>k+1,k,k</sub>,
 * Q<sub>k,k+1,k</sub>, Q<sub>k,k,k+1</sub>)</i> in 2D and 3D, resp.
 *
 * @author Guido Kanschat, 2005
 */
template <int dim>
class PolynomialsRaviartThomas
{
  public:
				     /**
				      * Constructor. Creates all basis
				      * functions for Raviart-Thomas polynomials
				      * of given degree.
				      *
				      * @arg k: the degree of the
				      * Raviart-Thomas-space, which is the degree
				      * of the largest tensor product
				      * polynomial space
				      * <i>Q<sub>k</sub></i> contained.
				      */
    PolynomialsRaviartThomas (const unsigned int k);

				     /**
				      * Computes the value and the
				      * first and second derivatives
				      * of each Raviart-Thomas
				      * polynomial at @p unit_point.
				      *
				      * The size of the vectors must
				      * either be zero or equal
				      * <tt>n()</tt>.  In the
				      * first case, the function will
				      * not compute these values.
				      *
				      * If you need values or
				      * derivatives of all tensor
				      * product polynomials then use
				      * this function, rather than
				      * using any of the
				      * <tt>compute_value</tt>,
				      * <tt>compute_grad</tt> or
				      * <tt>compute_grad_grad</tt>
				      * functions, see below, in a
				      * loop over all tensor product
				      * polynomials.
				      */
    void compute (const Point<dim>            &unit_point,
                  std::vector<Tensor<1,dim> > &values,
                  std::vector<Tensor<2,dim> > &grads,
                  std::vector<Tensor<3,dim> > &grad_grads) const;    
    
				     /**
				      * Returns the number of Raviart-Thomas polynomials.
				      */
    unsigned int n () const;
    
				     /**
				      * Returns the degree of the Raviart-Thomas
				      * space, which is one less than
				      * the highest polynomial degree.
				      */
    unsigned int degree () const;
    
				     /**
				      * Return the number of
				      * polynomials in the space
				      * <TT>RT(degree)</tt> without
				      * requiring to build an object
				      * of PolynomialsRaviartThomas. This is
				      * required by the FiniteElement
				      * classes.
				      */
    static unsigned int compute_n_pols(unsigned int degree);
    
  private:
				     /**
				      * An object representing the
				      * onedimensional polynomial
				      * space used here.
				      */
    const PolynomialSpace<1> polynomial_space;

				     /**
				      * Number of Raviart-Thomas
				      * polynomials.
				      */
    unsigned int n_pols;
    
				     /**
				      * Auxiliary memory.
				      */
    mutable std::vector<double> p_values;    
};



template <int dim>
inline unsigned int
PolynomialsRaviartThomas<dim>::n() const
{
  return n_pols;
}

template <int dim>
inline unsigned int
PolynomialsRaviartThomas<dim>::degree() const
{
  return polynomial_space.degree() -  1;
}

#endif
