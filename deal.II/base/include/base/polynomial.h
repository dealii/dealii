//----------------------------  polynomial.h  ---------------------------
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
//----------------------------  polynomial.h  ---------------------------
#ifndef __deal2__polynomial_h
#define __deal2__polynomial_h



#include <base/config.h>
#include <base/exceptions.h>
#include <base/subscriptor.h>

#include <vector>


/**
 * A namespace in which classes relating to the description of
 * 1d polynomial spaces are declared.
 */
namespace Polynomials
{

/**
 * Base class for all 1D polynomials. A polynomial is represented in
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
                                        * degree of the polynomial
                                        * represented by this object is
                                        * thus the number of elements in
                                        * the @p{coefficient} array
                                        * minus one.
                                        */
      Polynomial (const std::vector<number> &coefficients);
    
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
                                        * Degree of the polynomial. This
                                        * is the degree reflected by the
                                        * number of coefficients
                                        * provided by the
                                        * constructor. Leading non-zero
                                        * coefficients are not treated
                                        * separately.
                                        */
      unsigned int degree () const;

                                       /**
                                        * Scale the abscissa of the
                                        * polynomial.  Given the
                                        * polynomial $p(t)$ and the
                                        * scaling $t = ax$, then the
                                        * result of this operation is
                                        * the polynomial $q$, such that
                                        * $q(x) = p(t)$.
                                        *
                                        * The operation is performed in
                                        * place.
                                        */
      void scale (const number factor);

                                       /**
                                        * Shift the abscissa oft the
                                        * polynomial.  Given the
                                        * polynomial $p(t)$ and the
                                        * shift $t = x + a$, then the
                                        * result of this operation is
                                        * the polynomial $q$, such that
                                        * $q(x) = p(t)$.
                                        *
                                        * The template parameter allows
                                        * to compute the new
                                        * coefficients with higher
                                        * accuracy, since all
                                        * computations are performed
                                        * with type @p{number2}. This
                                        * may be necessary, since this
                                        * operation involves a big
                                        * number of additions. On a Sun
                                        * Sparc Ultra with Solaris 2.8,
                                        * the difference between
                                        * @p{double} and @p{long double}
                                        * was not significant, though.
                                        *
                                        * The operation is performed in
                                        * place, i.e. the coefficients
                                        * of the present object are
                                        * changed.
                                        */
      template <typename number2>
      void shift (const number2 offset);

                                       /**
                                        * Print coefficients.
                                        */
      void print(std::ostream& out) const;
				      
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
                                        * This function performs the
                                        * actual scaling.
                                        */
      static void scale (std::vector<number> &coefficients,
                         const number         factor);

                                       /**
                                        * This function performs the
                                        * actual shift
                                        */
      template <typename number2>
      static void shift (std::vector<number> &coefficients,
                         const number2        shift);

                                       /**
                                        * Multiply polynomial by a factor.
                                        */
      static void multiply (std::vector<number>& coefficients,
                            const number factor);
    
                                       /**
                                        * Coefficients of the polynomial
                                        * $\sum_i a_i x^i$. This vector
                                        * is filled by the constructor
                                        * of this class and may be
                                        * passed down by derived
                                        * classes.
                                        *
                                        * This vector cannot be constant
                                        * since we want to allow copying
                                        * of polynomials.
                                        */
      std::vector<number> coefficients;
  };



/**
 * Lagrange polynomials with equidistant interpolation points in
 * [0,1]. The polynomial of degree @p{n} has got @p{n+1} interpolation
 * points. The interpolation points are sorted in ascending
 * order. This order gives an index to each interpolation point.  A
 * Lagrangian polynomial equals to 1 at its `support point', and 0 at
 * all other interpolation points. For example, if the degree is 3,
 * and the support point is 1, then the polynomial represented by this
 * object is cubic and its value is 1 at the point @p{x=1/3}, and zero
 * at the point @p{x=0}, @p{x=2/3}, and @p{x=1}. All the polynomials
 * have polynomial order equal to @p{degree}, but together they span
 * the entire space of polynomials of degree less than or equal
 * @p{degree}.
 *
 * The Lagrange polynomials are implemented up to degree 10.
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
                                        * Return a vector of polynomial
                                        * objects of order @p{degree},
                                        * which then spans the full
                                        * space of polynomials up to the
                                        * given degree. The polynomials
                                        * are generated by calling the
                                        * destructor of this class with
                                        * the same degree but support
                                        * point running from zero to
                                        * @p{degree}. This function may
                                        * be used to initialize the
                                        * @ref{TensorProductPolynomials}
                                        * and @ref{PolynomialSpace}
                                        * classes.
                                        */
      static
      std::vector<Polynomial<double> >
      generate_complete_basis (const unsigned int degree);
    
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


//TODO[GK]: In contrast to the LagrangeEquidistant class, the following class is a template. Is this necessary, or could we make this more consistent?
/**
 * Legendre polynomials of arbitrary order on @p{[0,1]}.
 *
 * Constructing a Legendre polynomial of order @p{k}, the coefficients
 * will be computed by the three-term recursion formula.  The
 * coefficients are stored in a static data vector to be available
 * when needed next time. Since the recursion is performed for the
 * interval $[-1,1]$, the polynomials are shifted to $[0,1]$ by the
 * @p{scale} and @p{shift} functions of @p{Polynomial}, afterwards.
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
      Legendre (const unsigned int k);

                                       /**
                                        * Return a vector of Legendre
                                        * polynomial objects of orders
                                        * zero through @p{degree}, which
                                        * then spans the full space of
                                        * polynomials up to the given
                                        * degree. This function may be
                                        * used to initialize the
                                        * @ref{TensorProductPolynomials}
                                        * and @ref{PolynomialSpace}
                                        * classes.
                                        */
      static
      std::vector<Polynomial<number> >
      generate_complete_basis (const unsigned int degree);
    
    private:
                                       /**
                                        * Coefficients for the interval $[0,1]$.
                                        */
      static std::vector<const std::vector<number> *> shifted_coefficients;
    
                                       /**
                                        * Vector with already computed
                                        * coefficients. For each degree
                                        * of the polynomial, we keep one
                                        * pointer to the list of
                                        * coefficients; we do so rather
                                        * than keeping a vector of
                                        * vectors in order to simplify
                                        * programming multithread-safe.
                                        */
      static std::vector<const std::vector<number> *> recursive_coefficients;
    
                                       /**
                                        * Compute coefficients recursively.
                                        */
      static void compute_coefficients (const unsigned int k);
    
                                       /**
                                        * Get coefficients for
                                        * constructor.  This way, it can
                                        * use the non-standard
                                        * constructor of
                                        * @ref{Polynomial}.
                                        */
      static const std::vector<number> &
      get_coefficients (const unsigned int k);
  };



/**
 * Hierarchical polynomials of arbitrary order on @p{[0,1]}.
 *
 * When Constructing a Hierarchical polynomial of order @p{k}, 
 * the coefficients will be computed by a recursion formula.  The
 * coefficients are stored in a static data vector to be available
 * when needed next time.
 *
 * These hierarchical polynomials are based on those of Demkowicz, Oden, 
 * Rachowicz, and Hardy (CMAME 77 (1989) 79-112, Sec. 4). The first two
 * polynomials are the standard linear shape functions given by 
 * $\phi_{0}(x) = 1 - x$ and $\phi_{1}(x) = x$. For $l \geq 2$
 * we use the definitions $\phi_{l}(x) = (2x-1)^l - 1, l = 2,4,6,...$
 * and $\phi_{l}(x) = (2x-1)^l - (2x-1), l = 3,5,7,...$. These satisfy the 
 * recursion relations $\phi_{l}(x) = (2x-1)\phi_{l-1}, l=3,5,7,...$ and 
 * $\phi_{l}(x) = (2x-1)\phi_{l-1} + \phi_{2}, l=4,6,8,...$. 
 *
 * The degrees of freedom are the values at the vertices and the 
 * derivatives at the midpoint. Currently, we do not scale the
 * polynomials in any way, although better conditioning of the 
 * element stiffness matrix could possibly be achieved with scaling.
 *
 * @author Brian Carnes, 2002
 */
  template <typename number>
  class Hierarchical : public Polynomial<number>
  {
    public:
                                     /**
				      * Constructor for polynomial of
				      * order @p{p}.
				      */
      Hierarchical (const unsigned int p);

				     /**
				      * Return a vector of Hierarchical
				      * polynomial objects of orders
				      * zero through @p{degree}, which
				      * then spans the full space of
				      * polynomials up to the given
				      * degree. This function may be
				      * used to initialize the
				      * @ref{TensorProductPolynomials}
				      * and @ref{PolynomialSpace}
				      * classes.
				      */
      static
      std::vector<Polynomial<number> >
      generate_complete_basis (const unsigned int degree);
    
    private:
				     /**
				      * Compute coefficients recursively.
				      */
      static void compute_coefficients (const unsigned int k);

				     /**
				      * Get coefficients for
				      * constructor.  This way, it can
				      * use the non-standard
				      * constructor of
				      * @ref{Polynomial}.
				      */
     static const std::vector<number> &
     get_coefficients (const unsigned int k);
 
     static std::vector<const std::vector<number> *> recursive_coefficients;
   };  
}


/* -------------------------- inline functions --------------------- */

namespace Polynomials 
{
  template <typename number>
  inline
  unsigned int
  Polynomial<number>::degree () const
  {
    Assert (coefficients.size()>0, ExcVoidPolynomial());
    return coefficients.size() - 1;
  }
}

#endif
