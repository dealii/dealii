//-------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------
#ifndef __deal2__polynomials_P_h
#define __deal2__polynomials_P_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/tensor.h>
#include <base/point.h>
#include <base/polynomial.h>
#include <base/polynomial_space.h>
#include <base/table.h>

#include <vector>


/**
 * @brief The complete polynomial space of order <tt>k</tt> based on
 * the monomials.
 *
 * This class implements the polynomial space of order <tt>k</tt>
 * based on the monomials ${1,x,x^2,...}$. I.e. in <tt>d</tt>
 * dimensions it constructs all polynomials of the form $\prod_{i=1}^d
 * x_i^{n_i}$, where $\sum_i n_i\leq k$. The base polynomials are
 * given a specific ordering, e.g. in 2 dimensions:
 * ${1,x,y,xy,x^2,y^2,x^2y,xy^2,x^3,y^3,...}$. The ordering of the
 * monomials in $P_k1$ matches the ordering of the monomials in $P_k2$
 * for $k2>k1$.
 *
 * @author Ralf Hartmann, 2004
 */
template <int dim>
class PolynomialsP: public PolynomialSpace<dim>
{
  public:
				     /**
				      * Constructor. Creates all basis
				      * functions of $P_k$.
				      * @arg k: the degree of the
				      * polynomial space
				      */
    PolynomialsP (const unsigned int k);

				     /**
				      * Returns the degree <tt>k</tt>
				      * of the polynomial space
				      * <tt>P_k</tt>.
				      *
				      * Note, that this number is
				      * <tt>PolynomialSpace::degree()-1</tt>,
				      * compare definition in
				      * @ref{PolynomialSpace}.
				      */
    unsigned int degree() const;

				     /**
				      * For the <tt>n</tt>th
				      * polynomial $p_n(x,y,z)=x^i y^j
				      * z^k$ this function gives the
				      * degrees i,j,k in the x,y,z
				      * directions.
				      */
    void directional_degrees(unsigned int n,
			     unsigned int (&degrees)[dim]) const;
    
  private:

				     /**
				      * Fills the <tt>index_map</tt>.
				      */
    void create_polynomial_ordering(vector<unsigned int> &index_map) const;

				     /**
				      * Degree <tt>k<tt> of the
				      * polynomial space $P_k$,
				      * i.e. the number <tt>k<tt>
				      * which was given to the
				      * constructor.
				      */
    const unsigned int k;
};



template <int dim>
inline unsigned int
PolynomialsP<dim>::degree() const
{
  return k;
}


template <int dim>
inline void
PolynomialsP<dim>::directional_degrees(unsigned int n,
				       unsigned int (&degrees)[dim]) const
{
  compute_index(n,degrees);
}


#endif
