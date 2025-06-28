// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2000 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_polynomial_h
#define dealii_polynomial_h



#include <deal.II/base/config.h>

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/point.h>
#include <deal.II/base/types.h>

#include <array>
#include <limits>
#include <memory>
#include <shared_mutex>
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
   * Base class for all 1d polynomials. A polynomial is represented in this
   * class by its coefficients, which are set through the constructor or by
   * derived classes.
   *
   * There are two paths for evaluation of polynomials. One is based on the
   * coefficients which are evaluated through the Horner scheme which is a
   * robust general-purpose scheme. An alternative and more stable evaluation
   * of high-degree polynomials with roots in the unit interval is provided by
   * a product in terms of the roots. This form is available for special
   * polynomials such as Lagrange polynomials or Legendre polynomials and used
   * with the respective constructor. To obtain this more stable evaluation
   * form, the constructor with the roots in form of a Lagrange polynomial
   * must be used. In case a manipulation is done that changes the roots, the
   * representation is switched to the coefficient form.
   *
   * This class is a typical example of a possible template argument for the
   * TensorProductPolynomials class.
   */
  template <typename number>
  class Polynomial : public EnableObserverPointer
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
    Polynomial(const std::vector<number> &coefficients);

    /**
     * Constructor creating a zero polynomial of degree @p n.
     */
    Polynomial(const unsigned int n);

    /**
     * Constructor for a Lagrange polynomial and its point of evaluation. The
     * idea is to construct $\prod_{i\neq j} \frac{x-x_i}{x_j-x_i}$, where j
     * is the evaluation point specified as argument and the support points
     * contain all points (including x_j, which will internally not be
     * stored).
     */
    Polynomial(const std::vector<Point<1>> &lagrange_support_points,
               const unsigned int           evaluation_point);

    /**
     * Default constructor creating an illegal object.
     */
    Polynomial();

    /**
     * Return the value of this polynomial at the given point.
     *
     * This function uses the most numerically stable evaluation
     * algorithm for the provided form of the polynomial. If the
     * polynomial is in the product form of roots, the evaluation is
     * based on products of the form (x - x_i), whereas the Horner
     * scheme is used for polynomials in the coefficient form.
     */
    number
    value(const number x) const;

    /**
     * Return the values and the derivatives of the Polynomial at point
     * <tt>x</tt>.  <tt>values[i], i=0,...,values.size()-1</tt> includes the
     * <tt>i</tt>th derivative. The number of derivatives to be computed is
     * thus determined by the size of the array passed.
     *
     * This function uses the Horner scheme for numerical stability of the
     * evaluation for polynomials in the coefficient form or the product of
     * terms involving the roots if that representation is used.
     */
    void
    value(const number x, std::vector<number> &values) const;

    /**
     * Return the values and the derivatives of the Polynomial at point
     * <tt>x</tt>.  <tt>values[i], i=0,...,n_derivatives</tt> includes the
     * <tt>i</tt>th derivative. The number of derivatives to be computed is
     * determined by @p n_derivatives and @p values has to provide sufficient
     * space for @p n_derivatives + 1 values.
     *
     * This function uses the most numerically stable evaluation
     * algorithm for the provided form of the polynomial. If the
     * polynomial is in the product form of roots, the evaluation is
     * based on products of the form (x - x_i), whereas the Horner
     * scheme is used for polynomials in the coefficient form.
     *
     * The template type `Number2` must implement arithmetic
     * operations such as additions or multiplication with the type
     * `number` of the polynomial, and must be convertible from
     * `number` by `operator=`.
     */
    template <typename Number2>
    void
    value(const Number2      x,
          const unsigned int n_derivatives,
          Number2           *values) const;

    /**
     * Similar to the function above, but evaluate the polynomials on several
     * positions at once, as described by the array argument @p points. This
     * function is can be faster than the other function when the same
     * polynomial should be evaluated on several positions at once, e.g., the
     * x,y,z coordinates of a point for tensor-product polynomials.
     *
     * The template type `Number2` must implement arithmetic
     * operations such as additions or multiplication with the type
     * `number` of the polynomial, and must be convertible from
     * `number` by `operator=`.
     */
    template <std::size_t n_entries, typename Number2>
#ifndef DEBUG
    DEAL_II_ALWAYS_INLINE
#endif
      void
      values_of_array(const std::array<Number2, n_entries> &points,
                      const unsigned int                    n_derivatives,
                      std::array<Number2, n_entries>       *values) const;

    /**
     * Degree of the polynomial. This is the degree reflected by the number of
     * coefficients provided by the constructor. Leading non-zero coefficients
     * are not treated separately.
     */
    unsigned int
    degree() const;

    /**
     * Scale the abscissa of the polynomial.  Given the polynomial <i>p(t)</i>
     * and the scaling <i>t = ax</i>, then the result of this operation is the
     * polynomial <i>q</i>, such that <i>q(x) = p(t)</i>.
     *
     * The operation is performed in place.
     */
    void
    scale(const number factor);

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
    void
    shift(const number2 offset);

    /**
     * Compute the derivative of a polynomial.
     */
    Polynomial<number>
    derivative() const;

    /**
     * Compute the primitive of a polynomial. the coefficient of the zero
     * order term of the polynomial is zero.
     */
    Polynomial<number>
    primitive() const;

    /**
     * Multiply with a scalar.
     */
    Polynomial<number> &
    operator*=(const double s);

    /**
     * Multiply with another polynomial.
     */
    Polynomial<number> &
    operator*=(const Polynomial<number> &p);

    /**
     * Add a second polynomial.
     */
    Polynomial<number> &
    operator+=(const Polynomial<number> &p);

    /**
     * Subtract a second polynomial.
     */
    Polynomial<number> &
    operator-=(const Polynomial<number> &p);

    /**
     * Test for equality of two polynomials.
     */
    bool
    operator==(const Polynomial<number> &p) const;

    /**
     * Print coefficients.
     */
    void
    print(std::ostream &out) const;

    /**
     * Write or read the data of this object to or from a stream for the
     * purpose of serialization using the [BOOST serialization
     * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
     */
    template <class Archive>
    void
    serialize(Archive &ar, const unsigned int version);

    /**
     * Return an estimate (in bytes) for the memory consumption of this object.
     */
    virtual std::size_t
    memory_consumption() const;

  protected:
    /**
     * This function performs the actual scaling.
     */
    static void
    scale(std::vector<number> &coefficients, const number factor);

    /**
     * This function performs the actual shift
     */
    template <typename number2>
    static void
    shift(std::vector<number> &coefficients, const number2 shift);

    /**
     * Multiply polynomial by a factor.
     */
    static void
    multiply(std::vector<number> &coefficients, const number factor);

    /**
     * Transform polynomial form of product of linear factors into standard
     * form, $\sum_i a_i x^i$. Deletes all data structures related to the
     * product form.
     */
    void
    transform_into_standard_form();

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
   * Class generates Polynomial objects representing a monomial of degree n,
   * that is, the function $x^n$.
   */
  template <typename number>
  class Monomial : public Polynomial<number>
  {
  public:
    /**
     * Constructor, taking the degree of the monomial and an optional
     * coefficient as arguments.
     */
    Monomial(const unsigned int n, const double coefficient = 1.);

    /**
     * Return a vector of Monomial objects of degree zero through
     * <tt>degree</tt>, which then spans the full space of polynomials up to
     * the given degree. This function may be used to initialize the
     * TensorProductPolynomials and PolynomialSpace classes.
     */
    static std::vector<Polynomial<number>>
    generate_complete_basis(const unsigned int degree);

  private:
    /**
     * Needed by constructor.
     */
    static std::vector<number>
    make_vector(unsigned int n, const double coefficient);
  };


  /**
   * Lagrange polynomials with equidistant interpolation points in [0,1]. The
   * polynomial of degree <tt>n</tt> has got <tt>n+1</tt> interpolation
   * points. The interpolation points are sorted in ascending order. This
   * order gives an index to each interpolation point.  A Lagrangian
   * polynomial equals to 1 at its `support point', and 0 at all other
   * interpolation points. For example, if the degree is 3, and the support
   * point is 1, then the polynomial represented by this object is cubic and
   * its value is 1 at the point <tt>x=1/3</tt>, and zero at the point
   * <tt>x=0</tt>, <tt>x=2/3</tt>, and <tt>x=1</tt>. All the polynomials have
   * polynomial degree equal to <tt>degree</tt>, but together they span the
   * entire space of polynomials of degree less than or equal <tt>degree</tt>.
   *
   * The Lagrange polynomials are implemented up to degree 10.
   */
  class LagrangeEquidistant : public Polynomial<double>
  {
  public:
    /**
     * Constructor. Takes the degree <tt>n</tt> of the Lagrangian polynomial
     * and the index <tt>support_point</tt> of the support point. Fills the
     * <tt>coefficients</tt> of the base class Polynomial.
     */
    LagrangeEquidistant(const unsigned int n, const unsigned int support_point);

    /**
     * Return a vector of polynomial objects of degree <tt>degree</tt>, which
     * then spans the full space of polynomials up to the given degree. The
     * polynomials are generated by calling the constructor of this class with
     * the same degree but support point running from zero to <tt>degree</tt>.
     * This function may be used to initialize the TensorProductPolynomials
     * and PolynomialSpace classes.
     */
    static std::vector<Polynomial<double>>
    generate_complete_basis(const unsigned int degree);

  private:
    /**
     * Compute the <tt>coefficients</tt> of the base class Polynomial. This
     * function is <tt>static</tt> to allow to be called in the constructor.
     */
    static void
    compute_coefficients(const unsigned int   n,
                         const unsigned int   support_point,
                         std::vector<double> &a);
  };



  /**
   * Given a set of points along the real axis, this function returns all
   * Lagrange polynomials for interpolation of these points. The number of
   * polynomials is equal to the number of points and the maximum degree is
   * one less.
   */
  std::vector<Polynomial<double>>
  generate_complete_Lagrange_basis(const std::vector<Point<1>> &points);



  /**
   * Legendre polynomials of arbitrary degree. Constructing a Legendre
   * polynomial of degree <tt>p</tt>, the roots will be computed by the Gauss
   * formula of the respective number of points and a representation of the
   * polynomial by its roots.
   *
   * @note The polynomials defined by this class differ in two aspects by what
   * is usually referred to as Legendre polynomials: (i) This classes defines
   * them on the reference interval $[0,1]$, rather than the commonly used
   * interval $[-1,1]$. (ii) The polynomials have been scaled in such a way
   * that they are orthonormal, not just orthogonal; consequently, the
   * polynomials do not necessarily have boundary values equal to one.
   */
  class Legendre : public Polynomial<double>
  {
  public:
    /**
     * Constructor for polynomial of degree <tt>p</tt>.
     */
    Legendre(const unsigned int p);

    /**
     * Return a vector of Legendre polynomial objects of degrees zero through
     * <tt>degree</tt>, which then spans the full space of polynomials up to
     * the given degree. This function may be used to initialize the
     * TensorProductPolynomials and PolynomialSpace classes.
     */
    static std::vector<Polynomial<double>>
    generate_complete_basis(const unsigned int degree);
  };

  /**
   * Lobatto polynomials of arbitrary degree on <tt>[0,1]</tt>.
   *
   * These polynomials are the integrated Legendre polynomials on [0,1]. The
   * first two polynomials are the standard linear shape functions given by
   * $l_0(x) = 1-x$ and $l_1(x) = x$. For $i\geq2$ we use the definition
   * $l_i(x) = \frac{1}{\Vert L_{i-1}\Vert_2}\int_0^x L_{i-1}(t)\,dt$, where
   * $L_i$ denotes the $i$-th Legendre polynomial on $[0,1]$. The Lobatto
   * polynomials $l_0,\ldots,l_k$ form a complete basis of the polynomials
   * space of degree $k$.
   *
   * Calling the constructor with a given index <tt>k</tt> will generate the
   * polynomial with index <tt>k</tt>. But only for $k\geq 1$ the index equals
   * the degree of the polynomial. For <tt>k==0</tt> also a polynomial of
   * degree 1 is generated.
   *
   * These polynomials are used for the construction of the shape functions of
   * N&eacute;d&eacute;lec elements of arbitrary order.
   */
  class Lobatto : public Polynomial<double>
  {
  public:
    /**
     * Constructor for polynomial of degree <tt>p</tt>. There is an exception
     * for <tt>p==0</tt>, see the general documentation.
     */
    Lobatto(const unsigned int p = 0);

    /**
     * Return the polynomials with index <tt>0</tt> up to <tt>degree</tt>.
     * There is an exception for <tt>p==0</tt>, see the general documentation.
     */
    static std::vector<Polynomial<double>>
    generate_complete_basis(const unsigned int p);

  private:
    /**
     * Compute coefficients recursively.
     */
    std::vector<double>
    compute_coefficients(const unsigned int p);
  };



  /**
   * Hierarchical polynomials of arbitrary degree on <tt>[0,1]</tt>.
   *
   * When Constructing a Hierarchical polynomial of degree <tt>p</tt>, the
   * coefficients will be computed by a recursion formula.  The coefficients
   * are stored in a static data vector to be available when needed next time.
   *
   * These hierarchical polynomials are based on those of Demkowicz, Oden,
   * Rachowicz, and Hardy (CMAME 77 (1989) 79-112, Sec. 4). The first two
   * polynomials are the standard linear shape functions given by $\phi_{0}(x)
   * = 1 - x$ and $\phi_{1}(x) = x$. For $l \geq 2$ we use the definitions
   * $\phi_{l}(x) = (2x-1)^l - 1, l = 2,4,6,...$ and $\phi_{l}(x) = (2x-1)^l -
   * (2x-1), l = 3,5,7,...$. These satisfy the recursion relations
   * $\phi_{l}(x) = (2x-1)\phi_{l-1}, l=3,5,7,...$ and $\phi_{l}(x) =
   * (2x-1)\phi_{l-1} + \phi_{2}, l=4,6,8,...$.
   *
   * The degrees of freedom are the values at the vertices and the derivatives
   * at the midpoint. Currently, we do not scale the polynomials in any way,
   * although better conditioning of the element @ref GlossStiffnessMatrix "stiffness matrix" could
   * possibly be achieved with scaling.
   *
   * Calling the constructor with a given index <tt>p</tt> will generate the
   * following: if <tt>p==0</tt>, then the resulting polynomial is the linear
   * function associated with the left vertex, if <tt>p==1</tt> the one
   * associated with the right vertex. For higher values of <tt>p</tt>, you
   * get the polynomial of degree <tt>p</tt> that is orthogonal to all
   * previous ones. Note that for <tt>p==0</tt> you therefore do <b>not</b>
   * get a polynomial of degree zero, but one of degree one. This is to allow
   * generating a complete basis for polynomial spaces, by just iterating over
   * the indices given to the constructor.
   *
   * On the other hand, the function generate_complete_basis() creates a
   * complete basis of given degree. In order to be consistent with the
   * concept of a polynomial degree, if the given argument is zero, it does
   * <b>not</b> return the linear polynomial described above, but rather a
   * constant polynomial.
   */
  class Hierarchical : public Polynomial<double>
  {
  public:
    /**
     * Constructor for polynomial of degree <tt>p</tt>. There is an exception
     * for <tt>p==0</tt>, see the general documentation.
     */
    Hierarchical(const unsigned int p);

    /**
     * Return a vector of Hierarchical polynomial objects of degrees zero
     * through <tt>degree</tt>, which then spans the full space of polynomials
     * up to the given degree. Note that there is an exception if the given
     * <tt>degree</tt> equals zero, see the general documentation of this
     * class.
     *
     * This function may be used to initialize the TensorProductPolynomials,
     * AnisotropicPolynomials, and PolynomialSpace classes.
     */
    static std::vector<Polynomial<double>>
    generate_complete_basis(const unsigned int degree);

  private:
    /**
     * Compute coefficients recursively.
     */
    static void
    compute_coefficients(const unsigned int p);

    /**
     * Get coefficients for constructor.  This way, it can use the
     * non-standard constructor of Polynomial.
     */
    static const std::vector<double> &
    get_coefficients(const unsigned int p);

    /**
     * Vector with already computed coefficients. For each degree of the
     * polynomial, we keep one pointer to the list of coefficients; we do so
     * rather than keeping a vector of vectors in order to simplify
     * programming multithread-safe. In order to avoid memory leak, we use a
     * unique_ptr in order to correctly free the memory of the vectors when
     * the global destructor is called.
     */
    static std::vector<std::unique_ptr<const std::vector<double>>>
      recursive_coefficients;

    /**
     * The mutex that guards read and write access to the
     * `recursive_coefficients` array.
     */
    static std::shared_mutex coefficients_lock;
  };



  /**
   * Polynomials for Hermite interpolation condition.
   *
   * This is the set of polynomials of degree at least three, such that the
   * following interpolation conditions are met: the polynomials and their
   * first derivatives vanish at the values <i>x</i>=0 and <i>x</i>=1, with
   * the exceptions <i>p</i><sub>0</sub>(0)=1,
   * <i>p</i><sub><i>1</i></sub>(1)=1, <i>p</i>'<sub>2</sub>(0)=1,
   * <i>p'</i><sub>3</sub>(1)=1.
   *
   * For degree three, we obtain the standard four Hermitian interpolation
   * polynomials, see for instance <a
   * href="http://en.wikipedia.org/wiki/Cubic_Hermite_spline">Wikipedia</a>.
   * For higher degrees, these are augmented first, by the polynomial of
   * degree four with vanishing values and derivatives at <i>x</i>=0 and
   * <i>x</i>=1, then by the product of this fourth order polynomial with
   * Legendre polynomials of increasing order. The implementation is
   * @f{align*}{
   * p_0(x) &= 2x^3-3x^2+1 \\
   * p_1(x) &= -2x^3+3x^2 \\
   * p_2(x) &= x^3-2x^2+x  \\
   * p_3(x) &= x^3-x^2 \\
   * p_4(x) &= 16x^2(x-1)^2 \\
   * \ldots & \ldots \\
   * p_k(x) &= x^2(x-1)^2 L_{k-4}(x)
   * @f}
   */
  class HermiteInterpolation : public Polynomial<double>
  {
  public:
    /**
     * Constructor for polynomial with index <tt>p</tt>. See the class
     * documentation on the definition of the sequence of polynomials.
     */
    HermiteInterpolation(const unsigned int p);

    /**
     * Return the polynomials with index <tt>0</tt> up to <tt>p+1</tt> in a
     * space of degree up to <tt>p</tt>. Here, <tt>p</tt> has to be at least
     * 3.
     */
    static std::vector<Polynomial<double>>
    generate_complete_basis(const unsigned int p);
  };



  /**
   * Polynomials for a variant of Hermite polynomials with better condition
   * number in the interpolation than the basis from HermiteInterpolation.
   *
   * In analogy to the proper Hermite polynomials, this basis evaluates the
   * first polynomial $p_0$ to 1 at $x=0$ and has both a zero value and zero
   * derivative at $x=1$. Likewise, the last polynomial $p_n$ evaluates to 1
   * at $x=1$ with a zero value and zero derivative at $x=0$. The second
   * polynomial $p_1$ and the second to last polynomial $p_{n-1}$ represent
   * the derivative degree of freedom at $x=0$ and $x=1$, respectively.
   * They are zero at both the end points $x=0, x=1$ and have zero
   * derivative at the opposite end, $p_1'(1)=0$ and $p_{n-1}'(0)=0$. As
   * opposed to the original Hermite polynomials, $p_0$ does not have zero
   * derivative at $x=0$. The additional degree of freedom is used to make
   * $p_0$ and $p_1$ orthogonal, which for $n=3$ results in a root at
   * $x=\frac{2}{7}$ for $p_0$ and at $x=\frac{5}{7}$ for $p_n$,
   * respectively. Furthermore, the extension of these polynomials to higher
   * degrees $n>3$ is constructed by adding additional nodes inside the unit
   * interval, again ensuring better conditioning. The nodes are computed as
   * the roots of the Jacobi polynomials for $\alpha=\beta=4$, which are
   * orthogonal against the square of the generating function $x^2(1-x)^2$
   * with the Hermite
   * property. Then, these polynomials are constructed in the usual way as
   * Lagrange polynomials with double roots at $x=0$ and $x=1$. For example
   * with $n=4$, all of $p_0, p_1, p_3, p_4$ get an additional root at $x=0.5$
   * through the factor $(x-0.5)$. In summary, this basis is dominated by
   * nodal contributions, but it is not a nodal one because the second and
   * second to last polynomials that are non-nodal, and due to the presence of
   * double nodes in $x=0$ and $x=1$. The weights of the basis functions are
   * set such that the sum of all polynomials with unit weight represents the
   * constant function 1, similarly to Lagrange polynomials.
   *
   * The basis only contains Hermite information for <code>degree>=3</code>,
   * but it is also implemented for degrees between 0 and two. For the linear
   * case, the usual hat functions are implemented, whereas the polynomials
   * for <code>degree=2</code> are $p_0(x)=(1-x)^2$, $p_1(x)=2x(x-1)$, and
   * $p_2(x)=x^2$, in accordance with the construction principle for degree 3.
   *
   * These two relaxations improve the condition number of the @ref GlossMassMatrix "mass matrix"
   * (i.e., interpolation) significantly, as can be seen from the following
   * table:
   *
   * <table align="center" border="1">
   *   <tr>
   *    <th>&nbsp;</th>
   *    <th colspan="2">Condition number mass matrix</th>
   *   </tr>
   *   <tr>
   *    <th>degree</th>
   *    <th>HermiteInterpolation</th>
   *    <th>HermiteLikeInterpolation</th>
   *   </tr>
   *   <tr>
   *    <th>n=3</th>
   *    <th>1057</th>
   *    <th>17.18</th>
   *   </tr>
   *   <tr>
   *    <th>n=4</th>
   *    <th>6580</th>
   *    <th>16.83</th>
   *   </tr>
   *   <tr>
   *    <th>n=5</th>
   *    <th>1.875e+04</th>
   *    <th>15.99</th>
   *   </tr>
   *   <tr>
   *    <th>n=6</th>
   *    <th>6.033e+04</th>
   *    <th>16.34</th>
   *   </tr>
   *   <tr>
   *    <th>n=10</th>
   *    <th>9.756e+05</th>
   *    <th>20.70</th>
   *   </tr>
   *   <tr>
   *    <th>n=15</th>
   *    <th>9.431e+06</th>
   *    <th>27.91</th>
   *   </tr>
   *   <tr>
   *    <th>n=25</th>
   *    <th>2.220e+08</th>
   *    <th>43.54</th>
   *   </tr>
   *   <tr>
   *    <th>n=35</th>
   *    <th>2.109e+09</th>
   *    <th>59.51</th>
   *   </tr>
   * </table>
   *
   * This polynomial inherits the advantageous property of Hermite polynomials
   * where only two functions have value and/or derivative nonzero on a face
   * advantageous for discontinuous Galerkin methods
   * but gives better condition numbers of interpolation, which improves the
   * performance of some iterative schemes like conjugate gradients with
   * point-Jacobi. This polynomial is used in FE_DGQHermite.
   */
  class HermiteLikeInterpolation : public Polynomial<double>
  {
  public:
    /**
     * Constructor for the polynomial with index <tt>index</tt> within the set
     * up polynomials of degree @p degree.
     */
    HermiteLikeInterpolation(const unsigned int degree,
                             const unsigned int index);

    /**
     * Return the polynomials with index <tt>0</tt> up to <tt>degree+1</tt> in
     * a space of degree up to <tt>degree</tt>.
     */
    static std::vector<Polynomial<double>>
    generate_complete_basis(const unsigned int degree);
  };



  /*
   * Evaluate a Jacobi polynomial $ P_n^{\alpha, \beta}(x) $ specified by the
   * parameters @p alpha, @p beta, @p n, where @p n is the degree of the
   * Jacobi polynomial.
   *
   * @note The Jacobi polynomials are not orthonormal and are defined on the
   * unit interval $[0, 1]$ as usual for deal.II, rather than $[-1, +1]$ often
   * used in literature. @p x is the point of evaluation.
   */
  template <typename Number>
  Number
  jacobi_polynomial_value(const unsigned int degree,
                          const int          alpha,
                          const int          beta,
                          const Number       x);


  /**
   * Compute the roots of the Jacobi polynomials on the unit interval $[0, 1]$
   * of the given degree. These roots are used in several places inside the
   * deal.II library, such as the Gauss-Lobatto quadrature formula or for the
   * Hermite-like interpolation.
   *
   * The algorithm uses a Newton algorithm, using the zeros of the Chebyshev
   * polynomials as an initial guess. This code has been tested for alpha and
   * beta equal to zero (Legendre case), one (Gauss-Lobatto case) as well as
   * two, so be careful when using it for other values as the Newton iteration
   * might or might not converge.
   */
  template <typename Number>
  std::vector<Number>
  jacobi_polynomial_roots(const unsigned int degree,
                          const int          alpha,
                          const int          beta);
} // namespace Polynomials


/** @} */

/* -------------------------- inline functions --------------------- */

namespace Polynomials
{
  template <typename number>
  inline Polynomial<number>::Polynomial()
    : in_lagrange_product_form(false)
    , lagrange_weight(1.)
  {}



  template <typename number>
  inline unsigned int
  Polynomial<number>::degree() const
  {
    if (in_lagrange_product_form == true)
      {
        return lagrange_support_points.size();
      }
    else
      {
        Assert(coefficients.size() > 0, ExcEmptyObject());
        return coefficients.size() - 1;
      }
  }



  template <typename number>
  inline number
  Polynomial<number>::value(const number x) const
  {
    if (in_lagrange_product_form == false)
      {
        Assert(coefficients.size() > 0, ExcEmptyObject());

        // Horner scheme
        const unsigned int m     = coefficients.size();
        number             value = coefficients.back();
        for (int k = m - 2; k >= 0; --k)
          value = value * x + coefficients[k];
        return value;
      }
    else
      {
        // direct evaluation of Lagrange polynomial
        const unsigned int m     = lagrange_support_points.size();
        number             value = 1.;
        for (unsigned int j = 0; j < m; ++j)
          value *= x - lagrange_support_points[j];
        value *= lagrange_weight;
        return value;
      }
  }



  template <typename number>
  template <typename Number2>
  inline void
  Polynomial<number>::value(const Number2      x,
                            const unsigned int n_derivatives,
                            Number2           *values) const
  {
    values_of_array(std::array<Number2, 1ul>{{x}},
                    n_derivatives,
                    reinterpret_cast<std::array<Number2, 1ul> *>(values));
  }



  template <typename number>
  template <std::size_t n_entries, typename Number2>
  inline
#ifndef DEBUG
    DEAL_II_ALWAYS_INLINE
#endif
    void
    Polynomial<number>::values_of_array(
      const std::array<Number2, n_entries> &x,
      const unsigned int                    n_derivatives,
      std::array<Number2, n_entries>       *values) const
  {
    // evaluate Lagrange polynomial and derivatives
    if (in_lagrange_product_form == true)
      {
        // to compute the value and all derivatives of a polynomial of the
        // form (x-x_1)*(x-x_2)*...*(x-x_n), expand the derivatives like
        // automatic differentiation does.
        const unsigned int n_supp = lagrange_support_points.size();
        const number       weight = lagrange_weight;
        switch (n_derivatives)
          {
            default:
              for (unsigned int e = 0; e < n_entries; ++e)
                values[0][e] = weight;
              for (unsigned int k = 1; k <= n_derivatives; ++k)
                for (unsigned int e = 0; e < n_entries; ++e)
                  values[k][e] = 0.;
              for (unsigned int i = 0; i < n_supp; ++i)
                {
                  std::array<Number2, n_entries> v = x;
                  for (unsigned int e = 0; e < n_entries; ++e)
                    v[e] -= lagrange_support_points[i];

                  // multiply by (x-x_i) and compute action on all derivatives,
                  // too (inspired from automatic differentiation: implement the
                  // product rule for the old value and the new variable 'v',
                  // i.e., expand value v and derivative one). since we reuse a
                  // value from the next lower derivative from the steps before,
                  // need to start from the highest derivative
                  for (unsigned int k = n_derivatives; k > 0; --k)
                    for (unsigned int e = 0; e < n_entries; ++e)
                      values[k][e] = (values[k][e] * v[e] + values[k - 1][e]);
                  for (unsigned int e = 0; e < n_entries; ++e)
                    values[0][e] *= v[e];
                }
              // finally, multiply derivatives by k! to transform the product
              // p_n = p^(n)(x)/k! into the actual form of the derivative
              {
                number k_factorial = 2;
                for (unsigned int k = 2; k <= n_derivatives; ++k)
                  {
                    for (unsigned int e = 0; e < n_entries; ++e)
                      values[k][e] *= k_factorial;
                    k_factorial *= static_cast<number>(k + 1);
                  }
              }
              break;

            // manually implement case 0 (values only), case 1 (value + first
            // derivative), and case 2 (up to second derivative) since they
            // might be called often. then, we can unroll the inner loop and
            // keep the temporary results as local variables to help the
            // compiler with the pointer aliasing analysis.
            case 0:
              {
                std::array<Number2, n_entries> value;
                for (unsigned int e = 0; e < n_entries; ++e)
                  value[e] = weight;
                for (unsigned int i = 0; i < n_supp; ++i)
                  for (unsigned int e = 0; e < n_entries; ++e)
                    value[e] *= (x[e] - lagrange_support_points[i]);

                for (unsigned int e = 0; e < n_entries; ++e)
                  values[0][e] = value[e];
                break;
              }

            case 1:
              {
                std::array<Number2, n_entries> value, derivative = {};
                for (unsigned int e = 0; e < n_entries; ++e)
                  value[e] = weight;
                for (unsigned int i = 0; i < n_supp; ++i)
                  for (unsigned int e = 0; e < n_entries; ++e)
                    {
                      const Number2 v = x[e] - lagrange_support_points[i];
                      derivative[e]   = derivative[e] * v + value[e];
                      value[e] *= v;
                    }

                for (unsigned int e = 0; e < n_entries; ++e)
                  {
                    values[0][e] = value[e];
                    values[1][e] = derivative[e];
                  }
                break;
              }

            case 2:
              {
                std::array<Number2, n_entries> value, derivative = {},
                                                      second = {};
                for (unsigned int e = 0; e < n_entries; ++e)
                  value[e] = weight;
                for (unsigned int i = 0; i < n_supp; ++i)
                  for (unsigned int e = 0; e < n_entries; ++e)
                    {
                      const Number2 v = x[e] - lagrange_support_points[i];
                      second[e]       = second[e] * v + derivative[e];
                      derivative[e]   = derivative[e] * v + value[e];
                      value[e] *= v;
                    }

                for (unsigned int e = 0; e < n_entries; ++e)
                  {
                    values[0][e] = value[e];
                    values[1][e] = derivative[e];
                    values[2][e] = static_cast<number>(2) * second[e];
                  }
                break;
              }
          }
        return;
      }

    Assert(coefficients.size() > 0, ExcEmptyObject());

    // if derivatives are needed, then do it properly by the full
    // Horner scheme
    const unsigned int                          m = coefficients.size();
    std::vector<std::array<Number2, n_entries>> a(coefficients.size());
    for (unsigned int i = 0; i < coefficients.size(); ++i)
      for (unsigned int e = 0; e < n_entries; ++e)
        a[i][e] = coefficients[i];

    unsigned int j_factorial = 1;

    // loop over all requested derivatives. note that derivatives @p{j>m} are
    // necessarily zero, as they differentiate the polynomial more often than
    // the highest power is
    const unsigned int min_valuessize_m = std::min(n_derivatives + 1, m);
    for (unsigned int j = 0; j < min_valuessize_m; ++j)
      {
        for (int k = m - 2; k >= static_cast<int>(j); --k)
          for (unsigned int e = 0; e < n_entries; ++e)
            a[k][e] += x[e] * a[k + 1][e];
        for (unsigned int e = 0; e < n_entries; ++e)
          values[j][e] = static_cast<number>(j_factorial) * a[j][e];

        j_factorial *= j + 1;
      }

    // fill higher derivatives by zero
    for (unsigned int j = min_valuessize_m; j <= n_derivatives; ++j)
      for (unsigned int e = 0; e < n_entries; ++e)
        values[j][e] = 0.;
  }



  template <typename number>
  template <class Archive>
  inline void
  Polynomial<number>::serialize(Archive &ar, const unsigned int)
  {
    // forward to serialization function in the base class.
    ar &static_cast<EnableObserverPointer &>(*this);
    ar &coefficients;
    ar &in_lagrange_product_form;
    ar &lagrange_support_points;
    ar &lagrange_weight;
  }



  template <typename Number>
  Number
  jacobi_polynomial_value(const unsigned int degree,
                          const int          alpha,
                          const int          beta,
                          const Number       x)
  {
    Assert(alpha >= 0 && beta >= 0,
           ExcNotImplemented("Negative alpha/beta coefficients not supported"));
    // the Jacobi polynomial is evaluated using a recursion formula.
    Number p0, p1;

    // The recursion formula is defined for the interval [-1, 1], so rescale
    // to that interval here
    const Number xeval = Number(-1) + 2. * x;

    // initial values P_0(x), P_1(x):
    p0 = 1.0;
    if (degree == 0)
      return p0;
    p1 = ((alpha + beta + 2) * xeval + (alpha - beta)) / 2;
    if (degree == 1)
      return p1;

    for (unsigned int i = 1; i < degree; ++i)
      {
        const Number v  = 2 * i + (alpha + beta);
        const Number a1 = 2 * (i + 1) * (i + (alpha + beta + 1)) * v;
        const Number a2 = (v + 1) * (alpha * alpha - beta * beta);
        const Number a3 = v * (v + 1) * (v + 2);
        const Number a4 = 2 * (i + alpha) * (i + beta) * (v + 2);

        const Number pn = ((a2 + a3 * xeval) * p1 - a4 * p0) / a1;
        p0              = p1;
        p1              = pn;
      }
    return p1;
  }



  template <typename Number>
  std::vector<Number>
  jacobi_polynomial_roots(const unsigned int degree,
                          const int          alpha,
                          const int          beta)
  {
    std::vector<Number> x(degree, 0.5);

    // compute zeros with a Newton algorithm.

    // Set tolerance. For long double we might not always get the additional
    // precision in a run time environment (e.g. with valgrind), so we must
    // limit the tolerance to double. Since we do a Newton iteration, doing
    // one more iteration after the residual has indicated convergence will be
    // enough for all number types due to the quadratic convergence of
    // Newton's method

    const Number tolerance =
      4 * std::max(static_cast<Number>(std::numeric_limits<double>::epsilon()),
                   std::numeric_limits<Number>::epsilon());

    // The following implementation follows closely the one given in the
    // appendix of the book by Karniadakis and Sherwin: Spectral/hp element
    // methods for computational fluid dynamics (Oxford University Press,
    // 2005)

    // If symmetric, we only need to compute the half of points
    const unsigned int n_points = (alpha == beta ? degree / 2 : degree);
    for (unsigned int k = 0; k < n_points; ++k)
      {
        // we take the zeros of the Chebyshev polynomial (alpha=beta=-0.5) as
        // initial values, corrected by the initial value
        Number r = 0.5 - 0.5 * std::cos(static_cast<Number>(2 * k + 1) /
                                        (2 * degree) * numbers::PI);
        if (k > 0)
          r = (r + x[k - 1]) / 2;

        unsigned int converged = numbers::invalid_unsigned_int;
        for (unsigned int it = 1; it < 1000; ++it)
          {
            Number s = 0.;
            for (unsigned int i = 0; i < k; ++i)
              s += 1. / (r - x[i]);

            // derivative of P_n^{alpha,beta}, rescaled to [0, 1]
            const Number J_x =
              (alpha + beta + degree + 1) *
              jacobi_polynomial_value(degree - 1, alpha + 1, beta + 1, r);

            // value of P_n^{alpha,beta}
            const Number f = jacobi_polynomial_value(degree, alpha, beta, r);
            const Number delta = f / (f * s - J_x);
            r += delta;
            if (converged == numbers::invalid_unsigned_int &&
                std::abs(delta) < tolerance)
              converged = it;

            // do one more iteration to ensure accuracy also for tighter
            // types than double (e.g. long double)
            if (it == converged + 1)
              break;
          }

        Assert(converged != numbers::invalid_unsigned_int,
               ExcMessage("Newton iteration for zero of Jacobi polynomial "
                          "did not converge."));

        x[k] = r;
      }

    // in case we assumed symmetry, fill up the missing values
    for (unsigned int k = n_points; k < degree; ++k)
      x[k] = 1.0 - x[degree - k - 1];

    return x;
  }

} // namespace Polynomials
DEAL_II_NAMESPACE_CLOSE

#endif
