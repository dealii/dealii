//----------------------------  polynomial.h  ---------------------------
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
//----------------------------  polynomial.h  ---------------------------
#ifndef __deal2__polynomial_h
#define __deal2__polynomial_h



#include <base/exceptions.h>

#include <vector.h>

/**
 * Base class for all 1D polynomials.
 *
 * @author Ralf Hartmann, 2000
 */
class Polynomial
{
  public:
				     /**
				      * Constructor.
				      */
    Polynomial(const vector<double> &a);

				     /**
				      * Returns the values and the
				      * derivatives of the @p{Polynomial}
				      * at point @p{x}. @p{values[i],
				      * i=0,...,values.size()}
				      * includes the @p{i}th
				      * derivative.
				      *
				      * This function uses the Horner
				      * scheme.
				      */
    void value(double x, vector<double> &values) const;

  protected:

				     /**
				      * Coefficients of the polynomial
				      * $\sum_ia_ix^i$. This vector is
				      * filled by the constructor of
				      * derived classes.
				      */
    const vector<double> coefficients;
};



/**
 * Class of Lagrange polynomials with equidistant interpolation
 * points. The polynomial of order @p{n} has got @p{n+1} interpolation
 * points. The interpolation points are x=0, x=1 and x=intermediate
 * points in ]0,1[ in ascending order. This order gives an index to
 * each interpolation point.  A Lagrangian polynomial equals 1 at one
 * interpolation point that is called `support point', and 0 at all other
 * interpolation points.
 *
 * @author Ralf Hartmann, 2000
 */
class LagrangeEquidistant: public Polynomial
{
  public:
				     /**
				      * Constructor. Takes the order
				      * @p{n} of the Lagrangian
				      * polynom and the index
				      * @p{support_point} of the
				      * support point. Fills the
				      * @p{coefficients} of the base
				      * class @p{Polynomial}.
				      */
    LagrangeEquidistant(unsigned int n, unsigned int support_point);

  private:

				     /**
				      * Computes the @p{coefficients}
				      * of the base class
				      * @p{Polynomial}. This function
				      * is static to allow the
				      * @p{coefficients} to be a
				      * @p{const} vector.
				      */
    static vector<double> compute_coefficients(unsigned int n, unsigned int support_point);
};



#endif
