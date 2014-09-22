// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__polynomial_h
#define __deal2__polynomial_h



#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/point.h>
#include <deal.II/base/std_cxx11/shared_ptr.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup Polynomials
 * @{
 */

/**
 * A namespace in which classes relating to the description of 1d polynomial
 * spaces are declared.
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
   * @author Ralf Hartmann, Guido Kanschat, 2000, 2006
   */
  template <typename number>
  class Polynomial : public Subscriptor
  {
  public:
    /**
     * Constructor. The coefficients of the polynomial are passed as
     * arguments, and denote the polynomial $\sum_i a[i] x^i$, i.e. the first
     * element of the array denotes the constant term, the second the linear
     * one, and so on. The degree of the polynomial represented by this object
     * is thus the number of elements in the <tt>coefficient</tt> array minus
     * one.
     */
    Polynomial (const std::vector<number> &coefficients);

    /**
     * Constructor creating a zero polynomial of degree @p n.
     */
    Polynomial (const unsigned int n);

    /**
     * Constructor for Lagrange polynomial and its point of evaluation. The
     * idea is to construct $\prod_{i\neq j} \frac{x-x_i}{x_j-x_i}$, where j
     * is the evaluation point specified as argument and the support points
     * contain all points (including x_j, which will internally not be
     * stored).
     */
    Polynomial (const std::vector<Point<1> > &lagrange_support_points,
                const unsigned int            evaluation_point);

    /**
     * Default constructor creating an illegal object.
     */
    Polynomial ();

    /**
     * Return the value of this polynomial at the given point.
     *
     * This function uses the Horner scheme for numerical stability of the
     * evaluation.
     */
    number value (const number x) const;

    /**
     * Return the values and the derivatives of the Polynomial at point
     * <tt>x</tt>.  <tt>values[i], i=0,...,values.size()-1</tt> includes the
     * <tt>i</tt>th derivative. The number of derivatives to be computed is
     * thus determined by the size of the array passed.
     *
     * This function uses the Horner scheme for numerical stability of the
     * evaluation.
     */
    void value (const number         x,
                std::vector<number> &values) const;

    /**
     * Degree of the polynomial. This is the degree reflected by the number of
     * coefficients provided by the constructor. Leading non-zero coefficients
     * are not treated separately.
     */
    unsigned int degree () const;

    /**
     * Scale the abscissa of the polynomial.  Given the polynomial <i>p(t)</i>
     * and the scaling <i>t = ax</i>, then the result of this operation is the
     * polynomial <i>q</i>, such that <i>q(x) = p(t)</i>.
     *
     * The operation is performed in place.
     */
    void scale (const number factor);

    /**
     * Shift the abscissa oft the polynomial.  Given the polynomial
     * <i>p(t)</i> and the shift <i>t = x + a</i>, then the result of this
     * operation is the polynomial <i>q</i>, such that <i>q(x) = p(t)</i>.
     *
     * The template parameter allows to compute the new coefficients with
     * higher accuracy, since all computations are performed with type
     * <tt>number2</tt>. This may be necessary, since this operation involves
     * a big number of additions. On a Sun Sparc Ultra with Solaris 2.8, the
     * difference between <tt>double</tt> and <tt>long double</tt> was not
     * significant, though.
     *
     * The operation is performed in place, i.e. the coefficients of the
     * present object are changed.
     */
    template <typename number2>
    void shift (const number2 offset);

    /**
     * Compute the derivative of a polynomial.
     */
    Polynomial<number> derivative () const;

    /**
     * Compute the primitive of a polynomial. the coefficient of the zero
     * order term of the polynomial is zero.
     */
    Polynomial<number> primitive () const;

    /**
     * Multiply with a scalar.
     */
    Polynomial<number> &operator *= (const double s);

    /**
     * Multiply with another polynomial.
     */
    Polynomial<number> &operator *= (const Polynomial<number> &p);

    /**
     * Add a second polynomial.
     */
    Polynomial<number> &operator += (const Polynomial<number> &p);

    /**
     * Subtract a second polynomial.
     */
    Polynomial<number> &operator -= (const Polynomial<number> &p);

    /**
     *  Test for equality of two polynomials.
     */
    bool operator == (const Polynomial<number> &p)  const;

    /**
     * Print coefficients.
     */
    void print(std::ostream &out) const;

    /**
     * Write or read the data of this object to or from a stream for the
     * purpose of serialization.
     */
    template <class Archive>
    void serialize (Archive &ar, const unsigned int version);

  protected:

    /**
     * This function performs the actual scaling.
     */
    static void scale (std::vector<number> &coefficients,
                       const number         factor);

    /**
     * This function performs the actual shift
     */
    template <typename number2>
    static void shift (std::vector<number> &coefficients,
                       const number2        shift);

    /**
     * Multiply polynomial by a factor.
     */
    static void multiply (std::vector<number> &coefficients,
                          const number factor);

    /**
     * Transforms polynomial form of product of linear factors into standard
     * form, $\sum_i a_i x^i$. Deletes all data structures related to the
     * product form.
     */
    void transform_into_standard_form ();

    /**
     * Coefficients of the polynomial $\sum_i a_i x^i$. This vector is filled
     * by the constructor of this class and may be passed down by derived
     * classes.
     *
     * This vector cannot be constant since we want to allow copying of
     * polynomials.
     */
    std::vector<number> coefficients;

    /**
     * Stores whether the polynomial is in Lagrange product form, i.e.,
     * constructed as a product $(x-x_0) (x-x_1) \ldots (x-x_n)/c$, or not.
     */
    bool in_lagrange_product_form;

    /**
     * If the polynomial is in Lagrange product form, i.e., constructed as a
     * product $(x-x_0) (x-x_1) \ldots (x-x_n)/c$, store the shifts $x_i$.
     */
    std::vector<number> lagrange_support_points;

    /**
     * If the polynomial is in Lagrange product form, i.e., constructed as a
     * product $(x-x_0) (x-x_1) \ldots (x-x_n)/c$, store the weight c.
     */
    number lagrange_weight;
  };


  /**
   * Class generates Polynomial objects representing a monomial of
   * degree n, that is, the function $x^n$.
   *
   * @author Guido Kanschat, 2004
   */
  template <typename number>
  class Monomial : public Polynomial<number>
  {
  public:
    /**
     * Constructor, taking the
     * degree of the monomial and
     * an optional coefficient as
     * arguments.
     */
    Monomial(const unsigned int n,
             const double coefficient = 1.);

    /**
     * Return a vector of Monomial
     * objects of degree zero
     * through <tt>degree</tt>, which
     * then spans the full space of
     * polynomials up to the given
     * degree. This function may be
     * used to initialize the
     * TensorProductPolynomials
     * and PolynomialSpace
     * classes.
     */
    static
    std::vector<Polynomial<number> >
    generate_complete_basis (const unsigned int degree);

  private:
    /**
     * Needed by constructor.
     */
    static std::vector<number> make_vector(unsigned int n,
                                           const double coefficient);
  };


  /**
   * Lagrange polynomials with equidistant interpolation points in
   * [0,1]. The polynomial of degree <tt>n</tt> has got <tt>n+1</tt> interpolation
   * points. The interpolation points are sorted in ascending
   * order. This order gives an index to each interpolation point.  A
   * Lagrangian polynomial equals to 1 at its `support point', and 0 at
   * all other interpolation points. For example, if the degree is 3,
   * and the support point is 1, then the polynomial represented by this
   * object is cubic and its value is 1 at the point <tt>x=1/3</tt>, and zero
   * at the point <tt>x=0</tt>, <tt>x=2/3</tt>, and <tt>x=1</tt>. All the polynomials
   * have polynomial degree equal to <tt>degree</tt>, but together they span
   * the entire space of polynomials of degree less than or equal
   * <tt>degree</tt>.
   *
   * The Lagrange polynomials are implemented up to degree 10.
   *
   * @author Ralf Hartmann, 2000
   */
  class LagrangeEquidistant: public Polynomial<double>
  {
  public:
    /**
     * Constructor. Takes the degree <tt>n</tt> of the Lagrangian polynom and
     * the index <tt>support_point</tt> of the support point. Fills the
     * <tt>coefficients</tt> of the base class Polynomial.
     */
    LagrangeEquidistant (const unsigned int n,
                         const unsigned int support_point);

    /**
     * Return a vector of polynomial objects of degree <tt>degree</tt>, which
     * then spans the full space of polynomials up to the given degree. The
     * polynomials are generated by calling the constructor of this class with
     * the same degree but support point running from zero to
     * <tt>degree</tt>. This function may be used to initialize the
     * TensorProductPolynomials and PolynomialSpace classes.
     */
    static
    std::vector<Polynomial<double> >
    generate_complete_basis (const unsigned int degree);

  private:

    /**
     * Computes the <tt>coefficients</tt> of the base class Polynomial. This
     * function is <tt>static</tt> to allow to be called in the constructor.
     */
    static
    void
    compute_coefficients (const unsigned int n,
                          const unsigned int support_point,
                          std::vector<double> &a);
  };



  /**
   * Given a set of points along the real axis, this function returns all
   * Lagrange polynomials for interpolation of these points. The number of
   * polynomials is equal to the number of points and the maximum degree is
   * one less.
   */
  std::vector<Polynomial<double> >
  generate_complete_Lagrange_basis (const std::vector<Point<1> > &points);



  /**
   * Legendre polynomials of arbitrary degree.
   * Constructing a Legendre polynomial of degree <tt>p</tt>, the coefficients
   * will be computed by the three-term recursion formula.
   *
   * @note The polynomials defined by this class differ in two aspects by what
   * is usually referred to as Legendre polynomials: (i) This classes defines
   * them on the reference interval $[0,1]$, rather than the commonly used
   * interval $[-1,1]$. (ii) The polynomials have been scaled in such a way that
   * they are orthonormal, not just orthogonal; consequently, the polynomials do
   * not necessarily have boundary values equal to one.
   *
   * @author Guido Kanschat, 2000
   */
  class Legendre : public Polynomial<double>
  {
  public:
    /**
     * Constructor for polynomial of
     * degree <tt>p</tt>.
     */
    Legendre (const unsigned int p);

    /**
     * Return a vector of Legendre
     * polynomial objects of degrees
     * zero through <tt>degree</tt>, which
     * then spans the full space of
     * polynomials up to the given
     * degree. This function may be
     * used to initialize the
     * TensorProductPolynomials
     * and PolynomialSpace
     * classes.
     */
    static
    std::vector<Polynomial<double> >
    generate_complete_basis (const unsigned int degree);

  private:
    /**
     * Coefficients for the interval $[0,1]$.
     */
    static std::vector<std_cxx11::shared_ptr<const std::vector<double> > > shifted_coefficients;

    /**
     * Vector with already computed
     * coefficients. For each degree of the
     * polynomial, we keep one pointer to
     * the list of coefficients; we do so
     * rather than keeping a vector of
     * vectors in order to simplify
     * programming multithread-safe. In
     * order to avoid memory leak, we use a
     * shared_ptr in order to correctly
     * free the memory of the vectors when
     * the global destructor is called.
     */
    static std::vector<std_cxx11::shared_ptr<const std::vector<double> > > recursive_coefficients;

    /**
     * Compute coefficients recursively.
     * The coefficients are stored in a
     * static data vector to be available
     * when needed next time. Since the
     * recursion is performed for the
     * interval $[-1,1]$, the polynomials
     * are shifted to $[0,1]$ by the
     * <tt>scale</tt> and <tt>shift</tt>
     * functions of <tt>Polynomial</tt>,
     * afterwards.
     */
    static void compute_coefficients (const unsigned int p);

    /**
     * Get coefficients for
     * constructor.  This way, it can
     * use the non-standard
     * constructor of
     * Polynomial.
     */
    static const std::vector<double> &
    get_coefficients (const unsigned int k);
  };

  /**
   * Lobatto polynomials of arbitrary degree on <tt>[0,1]</tt>.
   *
   * These polynomials are the integrated Legendre polynomials on [0,1]. The
   * first two polynomials are the standard linear shape functions given by
   * $l_0(x) = 1-x$ and $l_1(x) = x$. For $i\geq2$ we use the definition $l_i(x)
   * = \frac{1}{\Vert L_{i-1}\Vert_2}\int_0^x L_{i-1}(t)\,dt$, where $L_i$
   * denotes the $i$-th Legendre polynomial on $[0,1]$. The Lobatto polynomials
   * $l_0,\ldots,l_k$ form a complete basis of the polynomials space of degree
   * $k$.
   *
   * Calling the constructor with a given index <tt>k</tt> will generate the
   * polynomial with index <tt>k</tt>. But only for $k\geq 1$ the index equals
   * the degree of the polynomial. For <tt>k==0</tt> also a polynomial of degree
   * 1 is generated.
   *
   * These polynomials are used for the construction of the shape functions of
   * N&eacute;d&eacute;lec elements of arbitrary order.
   *
   * @author Markus B&uuml;rg, 2009
   */
  class Lobatto : public Polynomial<double>
  {
  public:
    /**
     * Constructor for polynomial of degree
     * <tt>p</tt>. There is an exception
     * for <tt>p==0</tt>, see the general
     * documentation.
     */
    Lobatto (const unsigned int p = 0);

    /**
     * Return the polynomials with index
     * <tt>0</tt> up to
     * <tt>degree</tt>. There is an
     * exception for <tt>p==0</tt>, see the
     * general documentation.
     */
    static std::vector<Polynomial<double> >
    generate_complete_basis (const unsigned int p);

  private:
    /**
     * Compute coefficients recursively.
     */
    std::vector<double> compute_coefficients (const unsigned int p);
  };



  /**
   * Hierarchical polynomials of arbitrary degree on <tt>[0,1]</tt>.
   *
   * When Constructing a Hierarchical polynomial of degree <tt>p</tt>,
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
   * Calling the constructor with a given index <tt>p</tt> will generate the
   * following: if <tt>p==0</tt>, then the resulting polynomial is the linear
   * function associated with the left vertex, if <tt>p==1</tt> the one
   * associated with the right vertex. For higher values of <tt>p</tt>, you
   * get the polynomial of degree <tt>p</tt> that is orthogonal to all
   * previous ones. Note that for <tt>p==0</tt> you therefore do <b>not</b>
   * get a polynomial of degree zero, but one of degree one. This is to
   * allow generating a complete basis for polynomial spaces, by just
   * iterating over the indices given to the constructor.
   *
   * On the other hand, the function generate_complete_basis() creates
   * a complete basis of given degree. In order to be consistent with
   * the concept of a polynomial degree, if the given argument is zero,
   * it does <b>not</b> return the linear polynomial described above, but
   * rather a constant polynomial.
   *
   * @author Brian Carnes, 2002
   */
  class Hierarchical : public Polynomial<double>
  {
  public:
    /**
     * Constructor for polynomial of
     * degree <tt>p</tt>. There is an
     * exception for <tt>p==0</tt>, see
     * the general documentation.
     */
    Hierarchical (const unsigned int p);

    /**
     * Return a vector of
     * Hierarchical polynomial
     * objects of degrees zero through
     * <tt>degree</tt>, which then spans
     * the full space of polynomials
     * up to the given degree. Note
     * that there is an exception if
     * the given <tt>degree</tt> equals
     * zero, see the general
     * documentation of this class.
     *
     * This function may be
     * used to initialize the
     * TensorProductPolynomials,
     * AnisotropicPolynomials,
     * and PolynomialSpace
     * classes.
     */
    static
    std::vector<Polynomial<double> >
    generate_complete_basis (const unsigned int degree);

  private:
    /**
     * Compute coefficients recursively.
     */
    static void compute_coefficients (const unsigned int p);

    /**
     * Get coefficients for
     * constructor.  This way, it can
     * use the non-standard
     * constructor of
     * Polynomial.
     */
    static const std::vector<double> &
    get_coefficients (const unsigned int p);

    /**
     * Vector with already computed
     * coefficients. For each degree of the
     * polynomial, we keep one pointer to
     * the list of coefficients; we do so
     * rather than keeping a vector of
     * vectors in order to simplify
     * programming multithread-safe. In
     * order to avoid memory leak, we use a
     * shared_ptr in order to correctly
     * free the memory of the vectors when
     * the global destructor is called.
     */
    static std::vector<std_cxx11::shared_ptr<const std::vector<double> > > recursive_coefficients;
  };


  /**
   * Polynomials for Hermite interpolation condition.
   *
   * This is the set of polynomials of degree at least three, such that
   * the following interpolation conditions are met: the polynomials and
   * their first derivatives vanish at the values <i>x</i>=0 and
   * <i>x</i>=1, with the exceptions <i>p</i><sub>0</sub>(0)=1,
   * <i>p</i><sub><i>1</i></sub>(1)=1, <i>p</i>'<sub>2</sub>(0)=1,
   * <i>p'</i><sub>3</sub>(1)=1.
   *
   * For degree three, we obtain the standard four Hermitian
   * interpolation polynomials, see for instance <a
   * href="http://en.wikipedia.org/wiki/Cubic_Hermite_spline">Wikipedia</a>.
   * For higher degrees, these are augmented
   * first, by the polynomial of degree four with vanishing values and
   * derivatives at <i>x</i>=0 and <i>x</i>=1, then by the product of
   * this fourth order polynomial with Legendre polynomials of
   * increasing order. The implementation is
   * @f{align*}{
   * p_0(x) &= 2x^3-3x^2+1 \\
   * p_1(x) &= -2x^2+3x^2 \\
   * p_2(x) &= x^3-2x^2+x  \\
   * p_3(x) &= x^3-x^2 \\
   * p_4(x) &= 16x^2(x-1)^2 \\
   * \ldots & \ldots \\
   * p_k(x) &= x^2(x-1)^2 L_{k-4}(x)
   * @f}
   *
   * @author Guido Kanschat
   * @date 2012
   */
  class HermiteInterpolation : public Polynomial<double>
  {
  public:
    /**
                               * Constructor for polynomial
                               * with index <tt>p</tt>. See
                               * the class documentation on
                               * the definition of the
                               * sequence of polynomials.
                               */
    HermiteInterpolation (const unsigned int p);

    /**
     * Return the polynomials with index
     * <tt>0</tt> up to
     * <tt>p+1</tt> in a space of
     * degree up to
     * <tt>p</tt>. Here, <tt>p</tt>
     * has to be at least 3.
     */
    static std::vector<Polynomial<double> >
    generate_complete_basis (const unsigned int p);
  };
}


/** @} */

/* -------------------------- inline functions --------------------- */

namespace Polynomials
{
  template <typename number>
  inline
  Polynomial<number>::Polynomial ()
    :
    in_lagrange_product_form (false),
    lagrange_weight          (1.)
  {}



  template <typename number>
  inline
  unsigned int
  Polynomial<number>::degree () const
  {
    if (in_lagrange_product_form == true)
      {
        return lagrange_support_points.size();
      }
    else
      {
        Assert (coefficients.size()>0, ExcEmptyObject());
        return coefficients.size() - 1;
      }
  }



  template <typename number>
  inline
  number
  Polynomial<number>::value (const number x) const
  {
    if (in_lagrange_product_form == false)
      {
        Assert (coefficients.size() > 0, ExcEmptyObject());

        // Horner scheme
        const unsigned int m=coefficients.size();
        number value = coefficients.back();
        for (int k=m-2; k>=0; --k)
          value = value*x + coefficients[k];
        return value;
      }
    else
      {
        // direct evaluation of Lagrange polynomial
        const unsigned int m = lagrange_support_points.size();
        number value = 1.;
        for (unsigned int j=0; j<m; ++j)
          value *= x-lagrange_support_points[j];
        value *= lagrange_weight;
        return value;
      }
  }



  template <typename number>
  template <class Archive>
  inline
  void
  Polynomial<number>::serialize (Archive &ar, const unsigned int)
  {
    // forward to serialization function in the base class.
    ar &static_cast<Subscriptor &>(*this);
    ar &coefficients;
    ar &in_lagrange_product_form;
    ar &lagrange_support_points;
    ar &lagrange_weight;
  }

}
DEAL_II_NAMESPACE_CLOSE

#endif
