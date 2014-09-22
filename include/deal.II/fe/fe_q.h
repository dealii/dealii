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

#ifndef __deal2__fe_q_h
#define __deal2__fe_q_h

#include <deal.II/base/config.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/fe/fe_q_base.h>

DEAL_II_NAMESPACE_OPEN


/*!@addtogroup fe */
/*@{*/

/**
 * Implementation of a scalar Lagrange finite element @p Qp that yields the
 * finite element space of continuous, piecewise polynomials of degree @p p in
 * each coordinate direction. This class is realized using tensor product
 * polynomials based on equidistant or given support points.
 *
 * The standard constructor of this class takes the degree @p p of this finite
 * element. Alternatively, it can take a quadrature formula @p points defining
 * the support points of the Lagrange interpolation in one coordinate direction.
 *
 * For more information about the <tt>spacedim</tt> template parameter
 * check the documentation of FiniteElement or the one of
 * Triangulation.
 *
 * <h3>Implementation</h3>
 *
 * The constructor creates a TensorProductPolynomials object that includes the
 * tensor product of @p LagrangeEquidistant polynomials of degree @p p. This
 * @p TensorProductPolynomials object provides all values and derivatives of
 * the shape functions.  In case a quadrature rule is given, the constructor
 * creates a TensorProductPolynomials object that includes the tensor product
 * of @p Lagrange polynomials with the support points from @p points.
 *
 * Furthermore the constructor fills the @p interface_constraints, the
 * @p prolongation (embedding) and the @p restriction matrices. These
 * are implemented only up to a certain degree and may not be
 * available for very high polynomial degree.
 *
 *
 * <h3>Numbering of the degrees of freedom (DoFs)</h3>
 *
 * The original ordering of the shape functions represented by the
 * TensorProductPolynomials is a tensor product
 * numbering. However, the shape functions on a cell are renumbered
 * beginning with the shape functions whose support points are at the
 * vertices, then on the line, on the quads, and finally (for 3d) on
 * the hexes. To be explicit, these numberings are listed in the
 * following:
 *
 * <h4>Q1 elements</h4>
 * <ul>
 * <li> 1D case:
 *   @verbatim
 *      0-------1
 *   @endverbatim
 *
 * <li> 2D case:
 *   @verbatim
 *      2-------3
 *      |       |
 *      |       |
 *      |       |
 *      0-------1
 *   @endverbatim
 *
 * <li> 3D case:
 *   @verbatim
 *         6-------7        6-------7
 *        /|       |       /       /|
 *       / |       |      /       / |
 *      /  |       |     /       /  |
 *     4   |       |    4-------5   |
 *     |   2-------3    |       |   3
 *     |  /       /     |       |  /
 *     | /       /      |       | /
 *     |/       /       |       |/
 *     0-------1        0-------1
 *   @endverbatim
 *
 *   The respective coordinate values of the support points of the degrees
 *   of freedom are as follows:
 *   <ul>
 *   <li> Index 0: <tt>[0, 0, 0]</tt>;
 *   <li> Index 1: <tt>[1, 0, 0]</tt>;
 *   <li> Index 2: <tt>[0, 1, 0]</tt>;
 *   <li> Index 3: <tt>[1, 1, 0]</tt>;
 *   <li> Index 4: <tt>[0, 0, 1]</tt>;
 *   <li> Index 5: <tt>[1, 0, 1]</tt>;
 *   <li> Index 6: <tt>[0, 1, 1]</tt>;
 *   <li> Index 7: <tt>[1, 1, 1]</tt>;
 *   </ul>
 * </ul>
 * <h4>Q2 elements</h4>
 * <ul>
 * <li> 1D case:
 *   @verbatim
 *      0---2---1
 *   @endverbatim
 *
 * <li> 2D case:
 *   @verbatim
 *      2---7---3
 *      |       |
 *      4   8   5
 *      |       |
 *      0---6---1
 *   @endverbatim
 *
 * <li> 3D case:
 *   @verbatim
 *         6--15---7        6--15---7
 *        /|       |       /       /|
 *      12 |       19     12      1319
 *      /  18      |     /       /  |
 *     4   |       |    4---14--5   |
 *     |   2---11--3    |       |   3
 *     |  /       /     |      17  /
 *    16 8       9     16       | 9
 *     |/       /       |       |/
 *     0---10--1        0---10--1
 *
 *         *-------*        *-------*
 *        /|       |       /       /|
 *       / |  23   |      /  25   / |
 *      /  |       |     /       /  |
 *     *   |       |    *-------*   |
 *     |20 *-------*    |       |21 *
 *     |  /       /     |   22  |  /
 *     | /  24   /      |       | /
 *     |/       /       |       |/
 *     *-------*        *-------*
 *   @endverbatim
 *   The center vertex has number 26.
 *
 *   The respective coordinate values of the support points of the degrees
 *   of freedom are as follows:
 *   <ul>
 *   <li> Index 0: <tt>[0, 0, 0]</tt>;
 *   <li> Index 1: <tt>[1, 0, 0]</tt>;
 *   <li> Index 2: <tt>[0, 1, 0]</tt>;
 *   <li> Index 3: <tt>[1, 1, 0]</tt>;
 *   <li> Index 4: <tt>[0, 0, 1]</tt>;
 *   <li> Index 5: <tt>[1, 0, 1]</tt>;
 *   <li> Index 6: <tt>[0, 1, 1]</tt>;
 *   <li> Index 7: <tt>[1, 1, 1]</tt>;
 *   <li> Index 8: <tt>[0, 1/2, 0]</tt>;
 *   <li> Index 9: <tt>[1, 1/2, 0]</tt>;
 *   <li> Index 10: <tt>[1/2, 0, 0]</tt>;
 *   <li> Index 11: <tt>[1/2, 1, 0]</tt>;
 *   <li> Index 12: <tt>[0, 1/2, 1]</tt>;
 *   <li> Index 13: <tt>[1, 1/2, 1]</tt>;
 *   <li> Index 14: <tt>[1/2, 0, 1]</tt>;
 *   <li> Index 15: <tt>[1/2, 1, 1]</tt>;
 *   <li> Index 16: <tt>[0, 0, 1/2]</tt>;
 *   <li> Index 17: <tt>[1, 0, 1/2]</tt>;
 *   <li> Index 18: <tt>[0, 1, 1/2]</tt>;
 *   <li> Index 19: <tt>[1, 1, 1/2]</tt>;
 *   <li> Index 20: <tt>[0, 1/2, 1/2]</tt>;
 *   <li> Index 21: <tt>[1, 1/2, 1/2]</tt>;
 *   <li> Index 22: <tt>[1/2, 0, 1/2]</tt>;
 *   <li> Index 23: <tt>[1/2, 1, 1/2]</tt>;
 *   <li> Index 24: <tt>[1/2, 1/2, 0]</tt>;
 *   <li> Index 25: <tt>[1/2, 1/2, 1]</tt>;
 *   <li> Index 26: <tt>[1/2, 1/2, 1/2]</tt>;
 *   </ul>
 * </ul>
 * <h4>Q3 elements</h4>
 * <ul>
 * <li> 1D case:
 *   @verbatim
 *      0--2--3--1
 *   @endverbatim
 *
 * <li> 2D case:
 *   @verbatim
 *      2--10-11-3
 *      |        |
 *      5  14 15 7
 *      |        |
 *      4  12 13 6
 *      |        |
 *      0--8--9--1
 *   @endverbatim
 * </ul>
 * <h4>Q4 elements</h4>
 * <ul>
 * <li> 1D case:
 *   @verbatim
 *      0--2--3--4--1
 *   @endverbatim
 *
 * <li> 2D case:
 *   @verbatim
 *      2--13-14-15-3
 *      |           |
 *      6  22 23 24 9
 *      |           |
 *      5  19 20 21 8
 *      |           |
 *      4  16 17 18 7
 *      |           |
 *      0--10-11-12-1
 *   @endverbatim
 * </ul>
 *
 * @author Wolfgang Bangerth, 1998, 2003; Guido Kanschat, 2001; Ralf Hartmann, 2001, 2004, 2005; Oliver Kayser-Herold, 2004; Katharina Kormann, 2008; Martin Kronbichler, 2008
 */
template <int dim, int spacedim=dim>
class FE_Q : public FE_Q_Base<TensorProductPolynomials<dim>,dim,spacedim>
{
public:
  /**
   * Constructor for tensor product polynomials of degree @p p.
   */
  FE_Q (const unsigned int p);

  /**
   * Constructor for tensor product polynomials with support points @p points
   * based on a one-dimensional quadrature formula. The degree of the finite
   * element is <tt>points.size()-1</tt>.  Note that the first point has to be
   * 0 and the last one 1. If
   * <tt>FE_Q<dim>(QGaussLobatto<1>(fe_degree+1))</tt> is specified, so-called
   * Gauss-Lobatto elements are obtained which can give a diagonal mass matrix
   * if combined with Gauss-Lobatto quadrature on the same points. Their use
   * is shown in step-48.
   */
  FE_Q (const Quadrature<1> &points);

  /**
   * Constructs a FE_Q_isoQ1 element. That element shares large parts of code
   * with FE_Q so most of the construction work is done in this routine,
   * whereas the public constructor is in the class FE_Q_isoQ1.
   */
  FE_Q(const unsigned int subdivisions_per_dimension,
       const unsigned int base_degree);

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_Q<dim>(degree)</tt>, with @p dim and @p degree replaced by
   * appropriate values.
   */
  virtual std::string get_name () const;

protected:

  /**
   * @p clone function instead of a copy constructor.
   *
   * This function is needed by the constructors of @p FESystem.
   */
  virtual FiniteElement<dim,spacedim> *clone() const;
};



/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif
