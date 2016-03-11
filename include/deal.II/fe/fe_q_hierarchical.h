// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2015 by the deal.II authors
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

#ifndef dealii__fe_q_hierarchical_h
#define dealii__fe_q_hierarchical_h

#include <deal.II/base/config.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/fe/fe_poly.h>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup fe */
/*@{*/

/**
 * Implementation of hierarchical @p Qp shape functions that yield the finite
 * element space of continuous, piecewise polynomials of degree @p p. This
 * class is realized using tensor product polynomials based on a hierarchical
 * basis Polynomials::Hierarchical on the interval <tt>[0,1]</tt> which is
 * suitable for building an @p hp tensor product finite element if we assume
 * that each element has a single degree.
 *
 * The constructor of this class takes the degree @p p of this finite element.
 *
 * This class is not implemented for the codimension one case (<tt>spacedim !=
 * dim</tt>).
 *
 * <h3>Implementation</h3>
 *
 * The constructor creates a TensorProductPolynomials object that includes the
 * tensor product of @p Hierarchical polynomials of degree @p p. This @p
 * TensorProductPolynomials object provides all values and derivatives of the
 * shape functions.
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
 *
 * The $Q_1^H$ element is of polynomial degree one and, consequently, is
 * exactly the same as the $Q_1$ element in class FE_Q. In particular, the
 * shape function are defined in the exact same way:
 *
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
 * The respective coordinate values of the support points of the degrees of
 * freedom are as follows:
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
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q1/Q1H_shape0000.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q1/Q1H_shape0001.png
 * </td> </tr> <tr> <td align="center"> $Q_1^H$ element, shape function 0
 * </td>
 *
 * <td align="center"> $Q_1^H$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q1/Q1H_shape0002.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q1/Q1H_shape0003.png
 * </td> </tr> <tr> <td align="center"> $Q_1^H$ element, shape function 2
 * </td>
 *
 * <td align="center"> $Q_1^H$ element, shape function 3 </td> </tr> </table>
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
 *     0---10--1        0---8---1
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
 * The respective coordinate values of the support points of the degrees of
 * freedom are as follows:
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
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q2/Q2H_shape0000.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q2/Q2H_shape0001.png
 * </td> </tr> <tr> <td align="center"> $Q_2^H$ element, shape function 0
 * </td>
 *
 * <td align="center"> $Q_2^H$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q2/Q2H_shape0002.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q2/Q2H_shape0003.png
 * </td> </tr> <tr> <td align="center"> $Q_2^H$ element, shape function 2
 * </td>
 *
 * <td align="center"> $Q_2^H$ element, shape function 3 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q2/Q2H_shape0004.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q2/Q2H_shape0005.png
 * </td> </tr> <tr> <td align="center"> $Q_2^H$ element, shape function 4
 * </td>
 *
 * <td align="center"> $Q_2^H$ element, shape function 5 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q2/Q2H_shape0006.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q2/Q2H_shape0007.png
 * </td> </tr> <tr> <td align="center"> $Q_2^H$ element, shape function 6
 * </td>
 *
 * <td align="center"> $Q_2^H$ element, shape function 7 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q2/Q2H_shape0008.png
 * </td>
 *
 * <td align="center"> </td> </tr> <tr> <td align="center"> $Q_2^H$ element,
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
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0000.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0001.png
 * </td> </tr> <tr> <td align="center"> $Q_3^H$ element, shape function 0
 * </td>
 *
 * <td align="center"> $Q_3^H$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0002.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0003.png
 * </td> </tr> <tr> <td align="center"> $Q_3^H$ element, shape function 2
 * </td>
 *
 * <td align="center"> $Q_3^H$ element, shape function 3 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0004.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0005.png
 * </td> </tr> <tr> <td align="center"> $Q_3^H$ element, shape function 4
 * </td>
 *
 * <td align="center"> $Q_3^H$ element, shape function 5 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0006.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0007.png
 * </td> </tr> <tr> <td align="center"> $Q_3^H$ element, shape function 6
 * </td>
 *
 * <td align="center"> $Q_3^H$ element, shape function 7 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0008.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0009.png
 * </td> </tr> <tr> <td align="center"> $Q_3^H$ element, shape function 8
 * </td>
 *
 * <td align="center"> $Q_3^H$ element, shape function 9 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0010.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0011.png
 * </td> </tr> <tr> <td align="center"> $Q_3^H$ element, shape function 10
 * </td>
 *
 * <td align="center"> $Q_3^H$ element, shape function 11 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0012.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0013.png
 * </td> </tr> <tr> <td align="center"> $Q_3^H$ element, shape function 12
 * </td>
 *
 * <td align="center"> $Q_3^H$ element, shape function 13 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0014.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0015.png
 * </td> </tr> <tr> <td align="center"> $Q_3^H$ element, shape function 14
 * </td>
 *
 * <td align="center"> $Q_3^H$ element, shape function 15 </td> </tr> </table>
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
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0000.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0001.png
 * </td> </tr> <tr> <td align="center"> $Q_4^H$ element, shape function 0
 * </td>
 *
 * <td align="center"> $Q_4^H$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0002.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0003.png
 * </td> </tr> <tr> <td align="center"> $Q_4^H$ element, shape function 2
 * </td>
 *
 * <td align="center"> $Q_4^H$ element, shape function 3 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0004.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0005.png
 * </td> </tr> <tr> <td align="center"> $Q_4^H$ element, shape function 4
 * </td>
 *
 * <td align="center"> $Q_4^H$ element, shape function 5 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0006.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0007.png
 * </td> </tr> <tr> <td align="center"> $Q_4^H$ element, shape function 6
 * </td>
 *
 * <td align="center"> $Q_4^H$ element, shape function 7 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0008.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0009.png
 * </td> </tr> <tr> <td align="center"> $Q_4^H$ element, shape function 8
 * </td>
 *
 * <td align="center"> $Q_4^H$ element, shape function 9 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0010.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0011.png
 * </td> </tr> <tr> <td align="center"> $Q_4^H$ element, shape function 10
 * </td>
 *
 * <td align="center"> $Q_4^H$ element, shape function 11 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0012.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0013.png
 * </td> </tr> <tr> <td align="center"> $Q_4^H$ element, shape function 12
 * </td>
 *
 * <td align="center"> $Q_4^H$ element, shape function 13 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0014.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0015.png
 * </td> </tr> <tr> <td align="center"> $Q_4^H$ element, shape function 14
 * </td>
 *
 * <td align="center"> $Q_4^H$ element, shape function 15 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0016.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0017.png
 * </td> </tr> <tr> <td align="center"> $Q_4^H$ element, shape function 16
 * </td>
 *
 * <td align="center"> $Q_4^H$ element, shape function 17 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0018.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0019.png
 * </td> </tr> <tr> <td align="center"> $Q_4^H$ element, shape function 18
 * </td>
 *
 * <td align="center"> $Q_4^H$ element, shape function 19 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0020.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0021.png
 * </td> </tr> <tr> <td align="center"> $Q_4^H$ element, shape function 20
 * </td>
 *
 * <td align="center"> $Q_4^H$ element, shape function 21 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0022.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0023.png
 * </td> </tr> <tr> <td align="center"> $Q_4^H$ element, shape function 22
 * </td>
 *
 * <td align="center"> $Q_4^H$ element, shape function 23 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0024.png
 * </td>
 *
 * <td align="center"> </td> </tr> <tr> <td align="center"> $Q_4^H$ element,
 * shape function 24 </td>
 *
 * <td align="center"> </td> </tr> </table>
 *
 *
 *
 * @author Brian Carnes, 2002, Ralf Hartmann 2004, 2005, Denis Davydov, 2015
 */
template <int dim>
class FE_Q_Hierarchical : public FE_Poly<TensorProductPolynomials<dim>,dim>
{
public:
  /**
   * Constructor for tensor product polynomials of degree @p p.
   */
  FE_Q_Hierarchical (const unsigned int p);

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_Q_Hierarchical<dim>(degree)</tt>, with @p dim and @p
   * degree replaced by appropriate values.
   */
  virtual std::string get_name () const;

  /**
   * This function returns @p true, if the shape function @p shape_index has
   * non-zero function values somewhere on the face @p face_index.
   */
  virtual bool has_support_on_face (const unsigned int shape_index,
                                    const unsigned int face_index) const;

  /**
   * @name Functions to support hp
   * @{
   */

  /**
   * Return whether this element implements its hanging node constraints in
   * the new way, which has to be used to make elements "hp compatible".
   *
   * For the FE_Q_Hierarchical class the result is always true (independent of
   * the degree of the element), as it implements the complete set of
   * functions necessary for hp capability.
   */
  virtual bool hp_constraints_are_implemented () const;

  /**
   * Return the matrix interpolating from the given finite element to the
   * present one. Interpolation only between FE_Q_Hierarchical is supported.
   */
  virtual void get_interpolation_matrix(const FiniteElement< dim> &source,
                                        FullMatrix< double > &matrix) const;

  /**
   * Embedding matrix between grids. Only isotropic refinement is supported.
   */
  virtual const
  FullMatrix<double> &get_prolongation_matrix  (const unsigned int child,
                                                const RefinementCase<dim> &refinement_case = RefinementCase< dim >::isotropic_refinement) const;

  /**
   * If, on a vertex, several finite elements are active, the hp code first
   * assigns the degrees of freedom of each of these FEs different global
   * indices. It then calls this function to find out which of them should get
   * identical values, and consequently can receive the same global DoF index.
   * This function therefore returns a list of identities between DoFs of the
   * present finite element object with the DoFs of @p fe_other, which is a
   * reference to a finite element object representing one of the other finite
   * elements active on this particular vertex. The function computes which of
   * the degrees of freedom of the two finite element objects are equivalent,
   * both numbered between zero and the corresponding value of dofs_per_vertex
   * of the two finite elements. The first index of each pair denotes one of
   * the vertex dofs of the present element, whereas the second is the
   * corresponding index of the other finite element.
   */
  virtual
  std::vector<std::pair<unsigned int, unsigned int> >
  hp_vertex_dof_identities (const FiniteElement<dim> &fe_other) const;

  /**
   * Same as above but for lines.
   */
  virtual
  std::vector<std::pair<unsigned int, unsigned int> >
  hp_line_dof_identities (const FiniteElement<dim> &fe_other) const;

  /**
   * Same as above but for faces.
   */
  virtual
  std::vector<std::pair<unsigned int, unsigned int> >
  hp_quad_dof_identities (const FiniteElement<dim> &fe_other) const;

  /*@}*/

  /**
   * Return the matrix interpolating from a face of one element to the face of
   * the neighboring element. The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>.
   *
   * Derived elements will have to implement this function. They may only
   * provide interpolation matrices for certain source finite elements, for
   * example those from the same family. If they don't implement interpolation
   * from a given element, then they must throw an exception of type
   * <tt>FiniteElement<dim>::ExcInterpolationNotImplemented</tt>.
   */
  virtual void get_face_interpolation_matrix (const FiniteElement<dim> &source, FullMatrix<double> &matrix) const;

  /**
   * Return the matrix interpolating from a face of one element to the subface
   * of the neighboring element. The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>.
   *
   * Derived elements will have to implement this function. They may only
   * provide interpolation matrices for certain source finite elements, for
   * example those from the same family. If they don't implement interpolation
   * from a given element, then they must throw an exception of type
   * <tt>ExcInterpolationNotImplemented</tt>.
   */
  virtual void get_subface_interpolation_matrix (const FiniteElement<dim> &source, const unsigned int subface, FullMatrix<double> &matrix) const;

  /**
   * Return whether this element dominates the one given as argument when they
   * meet at a common face, whether it is the other way around, whether
   * neither dominates, or if either could dominate.
   *
   * For a definition of domination, see FiniteElementBase::Domination and in
   * particular the
   * @ref hp_paper "hp paper".
   */
  virtual
  FiniteElementDomination::Domination
  compare_for_face_domination (const FiniteElement<dim> &fe_other) const;

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   *
   * This function is made virtual, since finite element objects are usually
   * accessed through pointers to their base class, rather than the class
   * itself.
   */
  virtual std::size_t memory_consumption () const;

  /**
   * For a finite element of degree @p sub_degree < @p degree, we return a
   * vector which maps the numbering on an FE of degree @p sub_degree into the
   * numbering on this element.
   */
  std::vector<unsigned int> get_embedding_dofs (const unsigned int sub_degree) const;

  /**
   * Returns a list of constant modes of the element. For this element, the
   * list consists of true arguments for the first vertex shape functions and
   * false for the remaining ones.
   */
  virtual std::pair<Table<2,bool>, std::vector<unsigned int> >
  get_constant_modes () const;

  /**
   * This function is not implemented and throws an exception if called.
   */
  virtual
  void interpolate(std::vector<double>       &local_dofs,
                   const std::vector<double> &values) const;

  /**
   * This function is not implemented and throws an exception if called.
   */
  virtual
  void
  interpolate(std::vector<double>                &local_dofs,
              const std::vector<Vector<double> > &values,
              unsigned int offset = 0) const;

  /**
   * This function is not implemented and throws an exception if called.
   */
  virtual
  void
  interpolate(std::vector<double> &local_dofs,
              const VectorSlice<const std::vector<std::vector<double> > > &values) const;


protected:
  /**
   * @p clone function instead of a copy constructor.
   *
   * This function is needed by the constructors of @p FESystem.
   */
  virtual FiniteElement<dim> *clone() const;

private:

  /**
   * Only for internal use. Its full name is @p get_dofs_per_object_vector
   * function and it creates the @p dofs_per_object vector that is needed
   * within the constructor to be passed to the constructor of @p
   * FiniteElementData.
   */
  static std::vector<unsigned int> get_dpo_vector(const unsigned int degree);

  /**
   * The numbering of the degrees of freedom in continuous finite elements is
   * hierarchic, i.e. in such a way that we first number the vertex dofs, in
   * the order of the vertices as defined by the triangulation, then the line
   * dofs in the order and respecting the direction of the lines, then the
   * dofs on quads, etc.
   *
   * The dofs associated with 1d hierarchical polynomials are ordered with the
   * vertices first ($phi_0(x)=1-x$ and $phi_1(x)=x$) and then the line dofs
   * (the higher degree polynomials).  The 2d and 3d hierarchical polynomials
   * originate from the 1d hierarchical polynomials by tensor product. In the
   * following, the resulting numbering of dofs will be denoted by
   * `fe_q_hierarchical numbering`.
   *
   * This function constructs a table which fe_q_hierarchical index each
   * degree of freedom in the hierarchic numbering would have.
   *
   * This function is analogous to the
   * FETools::hierarchical_to_lexicographic_numbering() function. However, in
   * contrast to the fe_q_hierarchical numbering defined above, the
   * lexicographic numbering originates from the tensor products of
   * consecutive numbered dofs (like for LagrangeEquidistant).
   *
   * It is assumed that the size of the output argument already matches the
   * correct size, which is equal to the number of degrees of freedom in the
   * finite element.
   */
  static
  std::vector<unsigned int> hierarchic_to_fe_q_hierarchical_numbering (
    const FiniteElementData<dim> &fe);

  /**
   * This is an analogon to the previous function, but working on faces.
   */
  static
  std::vector<unsigned int>
  face_fe_q_hierarchical_to_hierarchic_numbering (const unsigned int degree);

  /**
   * Initialize two auxiliary fields that will be used in setting up the
   * various matrices in the constructor.
   */
  void build_dofs_cell (std::vector<FullMatrix<double> > &dofs_cell,
                        std::vector<FullMatrix<double> > &dofs_subcell) const;

  /**
   * Initialize the hanging node constraints matrices. Called from the
   * constructor.
   */
  void initialize_constraints (const std::vector<FullMatrix<double> > &dofs_subcell);

  /**
   * Initialize the embedding matrices. Called from the constructor.
   */
  void initialize_embedding_and_restriction (const std::vector<FullMatrix<double> > &dofs_cell,
                                             const std::vector<FullMatrix<double> > &dofs_subcell);

  /**
   * Initialize the @p unit_support_points field of the FiniteElement class.
   * Called from the constructor.
   */
  void initialize_unit_support_points ();

  /**
   * Initialize the @p unit_face_support_points field of the FiniteElement
   * class. Called from the constructor.
   */
  void initialize_unit_face_support_points ();

  /**
   * Mapping from lexicographic to shape function numbering on first face.
   */
  const std::vector<unsigned int> face_renumber;

  /**
   * Allow access from other dimensions. We need this since we want to call
   * the functions @p get_dpo_vector and @p
   * lexicographic_to_hierarchic_numbering for the faces of the finite element
   * of dimension dim+1.
   */
  template <int dim1> friend class FE_Q_Hierarchical;
};

/*@}*/

/* -------------- declaration of explicit specializations ------------- */

template <>
void FE_Q_Hierarchical<1>::initialize_unit_face_support_points ();

template <>
bool
FE_Q_Hierarchical<1>::has_support_on_face (const unsigned int,
                                           const unsigned int) const;

template <>
std::vector<unsigned int>
FE_Q_Hierarchical<1>::face_fe_q_hierarchical_to_hierarchic_numbering (const unsigned int);

DEAL_II_NAMESPACE_CLOSE

#endif
