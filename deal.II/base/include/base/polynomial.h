//----------------------------  polynomial.h  ---------------------------
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
//----------------------------  polynomial.h  ---------------------------
#ifndef __deal2__polynomial_h
#define __deal2__polynomial_h



#include <base/exceptions.h>
#include <base/subscriptor.h>

#include <vector>


/**
 * Base class for all 1D polynomials. A pollynomial is represented in
 * this class by its coefficients, which are set through the
 * constructor or by derived classes. Evaluation of a polynomial
 * happens through the Horner scheme which provides both numerical
 * stability and a minimal number of numerical operations.
 *
 * @author Ralf Hartmann, Guido Kanschat, 2000
 */
template <typename number>
class Polynomial : public Subscriptor
{
  public:
				     /**
				      * Constructor. The coefficients
				      * of the polynomial are passed
				      * as arguments, and denote the
				      * polynomial @p{\sum_i a[i]
				      * x^i}, i.e. the first element
				      * of the array denotes the
				      * constant term, the second the
				      * linear one, and so on. The
				      * order of the polynomial
				      * represented by this object is
				      * thus the number of elements in
				      * the @p{coefficient} array
				      * minus one.
				      */
    Polynomial (const std::vector<number> &coefficients);

                                     /**
				      * Default-Constructor.
				      */
    Polynomial ();
    
				     /**
				      * Return the value of this
				      * polynomial at the given point.
				      *
				      * This function uses the Horner
				      * scheme for numerical stability
				      * of the evaluation.
				      */
    number value (const number x) const;
    
				     /**
				      * Return the values and the
				      * derivatives of the
				      * @p{Polynomial} at point @p{x}.
				      * @p{values[i],
				      * i=0,...,values.size()-1}
				      * includes the @p{i}th
				      * derivative. The number of
				      * derivatives to be computed is
				      * thus determined by the size of
				      * the array passed.
				      *
				      * This function uses the Horner
				      * scheme for numerical stability
				      * of the evaluation.
				      */
    void value (const number         x,
		std::vector<number> &values) const;

				     /**
				      * Exception
				      */
    DeclException0 (ExcEmptyArray);
    
				     /**
				      * Exception
				      */
    DeclException0 (ExcVoidPolynomial);
    
  protected:

				     /**
				      * Coefficients of the polynomial
				      * $\sum_i a_i x^i$. This vector
				      * is filled by the constructor
				      * of this class and may be
				      * passed down by derived
				      * classes.
				      */
    std::vector<number> coefficients;
};



/**
 * Lagrange polynomials with equidistant interpolation
 * points in [0,1]. The polynomial of degree @p{n} has got @p{n+1} interpolation
 * points. The interpolation points are sorted in ascending
 * order. This order gives an index to each interpolation point.  A
 * Lagrangian polynomial equals to 1 at its `support point',
 * and 0 at all other interpolation
 * points. For example, if the degree is 3, and the support point is 1,
 * then the polynomial represented by this object is cubic and its
 * value is 1 at the point @p{x=1/3}, and zero at the point @p{x=0},
 * @p{x=2/3}, and @p{x=1}.
 *
 * @author Ralf Hartmann, 2000
 */
class LagrangeEquidistant: public Polynomial<double>
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
    LagrangeEquidistant (const unsigned int n,
			 const unsigned int support_point);

                                     /**
				      * Default-constructor.
				      */
    LagrangeEquidistant ();    

  private:

				     /**
				      * Computes the @p{coefficients}
				      * of the base class
				      * @p{Polynomial}. This function
				      * is @p{static} to allow to be
				      * called in the
				      * constructor.
				      */
    static 
    std::vector<double> 
    compute_coefficients (const unsigned int n,
			  const unsigned int support_point);
};


/**
 * Legendre polynomials of arbitrary order on @p{[-1,1]}.
 *
 * Constructing a Legendre polynomial of order @p{k}, the coefficients
 * will be computed by the three-term recursion formula.  The
 * coefficients are stored in a static data vector to be available
 * when needed next time.
 *
 * @author Guido Kanschat, 2000
 */
template <typename number>
class Legendre : public Polynomial<number>
{
public:
				   /**
				    * Constructor for polynomial of
				    * order @p{k}.
				    */
  Legendre (unsigned int k);
private:
				   /**
				    * Vector with already computed
				    * coefficients.
				    */
  static std::vector<std::vector<number> > coefficients;

				   /**
				    * Compute coefficients recursively.
				    */
  static void compute_coefficients (unsigned int k);

				   /**
				    * Get coefficients for
				    * constructor.  This way, it can
				    * use the non-standard constructor
				    * of @ref{Polynomial}.
				    */
  static const std::vector<number>& get_coefficients (unsigned int k);
};


#endif

