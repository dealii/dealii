// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2015 by the deal.II authors
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

#ifndef dealii__fe_q_h
#define dealii__fe_q_h

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
 * polynomials based on 1D Lagrange polynomials with equidistant (degree up to
 * 2), Gauss-Lobatto (starting from degree 3), or given support points.
 *
 * The standard constructor of this class takes the degree @p p of this finite
 * element. Alternatively, it can take a quadrature formula @p points defining
 * the support points of the Lagrange interpolation in one coordinate
 * direction.
 *
 * For more information about the <tt>spacedim</tt> template parameter check
 * the documentation of FiniteElement or the one of Triangulation.
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
 * Furthermore the constructor fills the @p interface_constraints, the @p
 * prolongation (embedding) and the @p restriction matrices. These are
 * implemented only up to a certain degree and may not be available for very
 * high polynomial degree.
 *
 * <h3>Unit support point distribution and conditioning of interpolation</h3>
 *
 * When constructing an FE_Q element at polynomial degrees one or two,
 * equidistant support points at 0 and 1 (linear case) or 0, 0.5, and 1
 * (quadratic case) are used. The unit support or nodal points
 * <i>x<sub>i</sub></i> are those points where the <i>j</i>th Lagrange
 * polynomial satisfies the $\delta_{ij}$ property, i.e., where one polynomial
 * is one and all the others are zero.  For higher polynomial degrees, the
 * support points are non-equidistant by default, and chosen to be the support
 * points of the <tt>(degree+1)</tt>-order Gauss-Lobatto quadrature rule. This
 * point distribution yields well-conditioned Lagrange interpolation at
 * arbitrary polynomial degrees. By contrast, polynomials based on equidistant
 * points get increasingly ill-conditioned as the polynomial degree
 * increases. In interpolation, this effect is known as the Runge
 * phenomenon. For Galerkin methods, the Runge phenomenon is typically not
 * visible in the solution quality but rather in the condition number of the
 * associated system matrices. For example, the elemental mass matrix of
 * equidistant points at degree 10 has condition number 2.6e6, whereas the
 * condition number for Gauss-Lobatto points is around 400.
 *
 * The Gauss-Lobatto points in 1D include the end points 0 and +1 of the unit
 * interval. The interior points are shifted towards the end points, which
 * gives a denser point distribution close to the element boundary.
 *
 * If combined with Gauss-Lobatto quadrature, FE_Q based on the default
 * support points gives diagonal mass matrices. This case is demonstrated in
 * step-48. However, this element can be combined with arbitrary quadrature
 * rules through the usual FEValues approach, including full Gauss
 * quadrature. In the general case, the mass matrix is non-diagonal.
 *
 * <h3>Numbering of the degrees of freedom (DoFs)</h3>
 *
 * The original ordering of the shape functions represented by the
 * TensorProductPolynomials is a tensor product numbering. However, the shape
 * functions on a cell are renumbered beginning with the shape functions whose
 * support points are at the vertices, then on the line, on the quads, and
 * finally (for 3d) on the hexes. To be explicit, these numberings are listed
 * in the following:
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
 * The respective coordinate values of the support points of the shape
 * functions are as follows:
 * <ul>
 * <li> Shape function 0: <tt>[0, 0, 0]</tt>;
 * <li> Shape function 1: <tt>[1, 0, 0]</tt>;
 * <li> Shape function 2: <tt>[0, 1, 0]</tt>;
 * <li> Shape function 3: <tt>[1, 1, 0]</tt>;
 * <li> Shape function 4: <tt>[0, 0, 1]</tt>;
 * <li> Shape function 5: <tt>[1, 0, 1]</tt>;
 * <li> Shape function 6: <tt>[0, 1, 1]</tt>;
 * <li> Shape function 7: <tt>[1, 1, 1]</tt>;
 * </ul>
 * </ul>
 *
 * In 2d, these shape functions look as follows: <table> <tr> <td
 * align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q1/Q1_shape0000.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q1/Q1_shape0001.png
 * </td> </tr> <tr> <td align="center"> $Q_1$ element, shape function 0 </td>
 *
 * <td align="center"> $Q_1$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q1/Q1_shape0002.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q1/Q1_shape0003.png
 * </td> </tr> <tr> <td align="center"> $Q_1$ element, shape function 2 </td>
 *
 * <td align="center"> $Q_1$ element, shape function 3 </td> </tr> </table>
 *
 *
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
 * The center vertex has number 26.
 *
 * The respective coordinate values of the support points of the shape
 * functions are as follows:
 * <ul>
 * <li> Shape function 0: <tt>[0, 0, 0]</tt>;
 * <li> Shape function 1: <tt>[1, 0, 0]</tt>;
 * <li> Shape function 2: <tt>[0, 1, 0]</tt>;
 * <li> Shape function 3: <tt>[1, 1, 0]</tt>;
 * <li> Shape function 4: <tt>[0, 0, 1]</tt>;
 * <li> Shape function 5: <tt>[1, 0, 1]</tt>;
 * <li> Shape function 6: <tt>[0, 1, 1]</tt>;
 * <li> Shape function 7: <tt>[1, 1, 1]</tt>;
 * <li> Shape function 8: <tt>[0, 1/2, 0]</tt>;
 * <li> Shape function 9: <tt>[1, 1/2, 0]</tt>;
 * <li> Shape function 10: <tt>[1/2, 0, 0]</tt>;
 * <li> Shape function 11: <tt>[1/2, 1, 0]</tt>;
 * <li> Shape function 12: <tt>[0, 1/2, 1]</tt>;
 * <li> Shape function 13: <tt>[1, 1/2, 1]</tt>;
 * <li> Shape function 14: <tt>[1/2, 0, 1]</tt>;
 * <li> Shape function 15: <tt>[1/2, 1, 1]</tt>;
 * <li> Shape function 16: <tt>[0, 0, 1/2]</tt>;
 * <li> Shape function 17: <tt>[1, 0, 1/2]</tt>;
 * <li> Shape function 18: <tt>[0, 1, 1/2]</tt>;
 * <li> Shape function 19: <tt>[1, 1, 1/2]</tt>;
 * <li> Shape function 20: <tt>[0, 1/2, 1/2]</tt>;
 * <li> Shape function 21: <tt>[1, 1/2, 1/2]</tt>;
 * <li> Shape function 22: <tt>[1/2, 0, 1/2]</tt>;
 * <li> Shape function 23: <tt>[1/2, 1, 1/2]</tt>;
 * <li> Shape function 24: <tt>[1/2, 1/2, 0]</tt>;
 * <li> Shape function 25: <tt>[1/2, 1/2, 1]</tt>;
 * <li> Shape function 26: <tt>[1/2, 1/2, 1/2]</tt>;
 * </ul>
 * </ul>
 *
 *
 * In 2d, these shape functions look as follows (the black plane corresponds
 * to zero; negative shape function values may not be visible): <table> <tr>
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q2/Q2_shape0000.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q2/Q2_shape0001.png
 * </td> </tr> <tr> <td align="center"> $Q_2$ element, shape function 0 </td>
 *
 * <td align="center"> $Q_2$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q2/Q2_shape0002.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q2/Q2_shape0003.png
 * </td> </tr> <tr> <td align="center"> $Q_2$ element, shape function 2 </td>
 *
 * <td align="center"> $Q_2$ element, shape function 3 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q2/Q2_shape0004.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q2/Q2_shape0005.png
 * </td> </tr> <tr> <td align="center"> $Q_2$ element, shape function 4 </td>
 *
 * <td align="center"> $Q_2$ element, shape function 5 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q2/Q2_shape0006.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q2/Q2_shape0007.png
 * </td> </tr> <tr> <td align="center"> $Q_2$ element, shape function 6 </td>
 *
 * <td align="center"> $Q_2$ element, shape function 7 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q2/Q2_shape0008.png
 * </td>
 *
 * <td align="center"> </td> </tr> <tr> <td align="center"> $Q_2$ element,
 * shape function 8 </td>
 *
 * <td align="center"> </td> </tr> </table>
 *
 *
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
 *
 * In 2d, these shape functions look as follows (the black plane corresponds
 * to zero; negative shape function values may not be visible): <table> <tr>
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0000.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0001.png
 * </td> </tr> <tr> <td align="center"> $Q_3$ element, shape function 0 </td>
 *
 * <td align="center"> $Q_3$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0002.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0003.png
 * </td> </tr> <tr> <td align="center"> $Q_3$ element, shape function 2 </td>
 *
 * <td align="center"> $Q_3$ element, shape function 3 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0004.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0005.png
 * </td> </tr> <tr> <td align="center"> $Q_3$ element, shape function 4 </td>
 *
 * <td align="center"> $Q_3$ element, shape function 5 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0006.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0007.png
 * </td> </tr> <tr> <td align="center"> $Q_3$ element, shape function 6 </td>
 *
 * <td align="center"> $Q_3$ element, shape function 7 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0008.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0009.png
 * </td> </tr> <tr> <td align="center"> $Q_3$ element, shape function 8 </td>
 *
 * <td align="center"> $Q_3$ element, shape function 9 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0010.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0011.png
 * </td> </tr> <tr> <td align="center"> $Q_3$ element, shape function 10 </td>
 *
 * <td align="center"> $Q_3$ element, shape function 11 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0012.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0013.png
 * </td> </tr> <tr> <td align="center"> $Q_3$ element, shape function 12 </td>
 *
 * <td align="center"> $Q_3$ element, shape function 13 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0014.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0015.png
 * </td> </tr> <tr> <td align="center"> $Q_3$ element, shape function 14 </td>
 *
 * <td align="center"> $Q_3$ element, shape function 15 </td> </tr> </table>
 *
 *
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
 * In 2d, these shape functions look as follows (the black plane corresponds
 * to zero; negative shape function values may not be visible): <table> <tr>
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0000.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0001.png
 * </td> </tr> <tr> <td align="center"> $Q_4$ element, shape function 0 </td>
 *
 * <td align="center"> $Q_4$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0002.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0003.png
 * </td> </tr> <tr> <td align="center"> $Q_4$ element, shape function 2 </td>
 *
 * <td align="center"> $Q_4$ element, shape function 3 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0004.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0005.png
 * </td> </tr> <tr> <td align="center"> $Q_4$ element, shape function 4 </td>
 *
 * <td align="center"> $Q_4$ element, shape function 5 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0006.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0007.png
 * </td> </tr> <tr> <td align="center"> $Q_4$ element, shape function 6 </td>
 *
 * <td align="center"> $Q_4$ element, shape function 7 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0008.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0009.png
 * </td> </tr> <tr> <td align="center"> $Q_4$ element, shape function 8 </td>
 *
 * <td align="center"> $Q_4$ element, shape function 9 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0010.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0011.png
 * </td> </tr> <tr> <td align="center"> $Q_4$ element, shape function 10 </td>
 *
 * <td align="center"> $Q_4$ element, shape function 11 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0012.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0013.png
 * </td> </tr> <tr> <td align="center"> $Q_4$ element, shape function 12 </td>
 *
 * <td align="center"> $Q_4$ element, shape function 13 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0014.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0015.png
 * </td> </tr> <tr> <td align="center"> $Q_4$ element, shape function 14 </td>
 *
 * <td align="center"> $Q_4$ element, shape function 15 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0016.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0017.png
 * </td> </tr> <tr> <td align="center"> $Q_4$ element, shape function 16 </td>
 *
 * <td align="center"> $Q_4$ element, shape function 17 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0018.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0019.png
 * </td> </tr> <tr> <td align="center"> $Q_4$ element, shape function 18 </td>
 *
 * <td align="center"> $Q_4$ element, shape function 19 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0020.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0021.png
 * </td> </tr> <tr> <td align="center"> $Q_4$ element, shape function 20 </td>
 *
 * <td align="center"> $Q_4$ element, shape function 21 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0022.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0023.png
 * </td> </tr> <tr> <td align="center"> $Q_4$ element, shape function 22 </td>
 *
 * <td align="center"> $Q_4$ element, shape function 23 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0024.png
 * </td>
 *
 * <td align="center"> </td> </tr> <tr> <td align="center"> $Q_4$ element,
 * shape function 24 </td>
 *
 * <td align="center"> </td> </tr> </table>
 *
 *
 *
 * @author Wolfgang Bangerth, 1998, 2003; Guido Kanschat, 2001; Ralf Hartmann,
 * 2001, 2004, 2005; Oliver Kayser-Herold, 2004; Katharina Kormann, 2008;
 * Martin Kronbichler, 2008
 */
template <int dim, int spacedim=dim>
class FE_Q : public FE_Q_Base<TensorProductPolynomials<dim>,dim,spacedim>
{
public:
  /**
   * Constructor for tensor product polynomials of degree @p p based on
   * Gauss-Lobatto support (node) points. For polynomial degrees of one and
   * two, these are the usual equidistant points.
   */
  FE_Q (const unsigned int p);

  /**
   * Constructor for tensor product polynomials with support points @p points
   * based on a one-dimensional quadrature formula. The degree of the finite
   * element is <tt>points.size()-1</tt>. Note that the first point has to be
   * 0 and the last one 1. Constructing
   * <tt>FE_Q<dim>(QGaussLobatto<1>(fe_degree+1))</tt> is equivalent to the
   * constructor that specifies the polynomial degree only. For selecting
   * equidistant nodes at <tt>fe_degree > 2</tt>, construct
   * <tt>FE_Q<dim>(QIterated<1>(QTrapez<1>(),fe_degree))</tt>.
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
