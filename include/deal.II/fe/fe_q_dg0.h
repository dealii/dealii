// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2014 by the deal.II authors
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


#ifndef __deal2__fe_q_dg0_h
#define __deal2__fe_q_dg0_h

#include <deal.II/base/config.h>
#include <deal.II/base/tensor_product_polynomials_const.h>
#include <deal.II/fe/fe_q_base.h>

DEAL_II_NAMESPACE_OPEN


/*!@addtogroup fe */
/*@{*/

/**
 * Implementation of a scalar Lagrange finite element @p Qp+DG0 that yields the
 * finite element space of continuous, piecewise polynomials of degree @p p in
 * each coordinate direction plus the space of locally constant functions.
 * This class is realized using tensor product polynomials based on equidistant
 * or given support points.
 *
 * The standard constructor of this class takes the degree @p p of this finite
 * element. Alternatively, it can take a quadrature formula @p points defining
 * the support points of the Lagrange interpolation in one coordinate direction.
 *
 * For more information about the <tt>spacedim</tt> template parameter check
 * the documentation of FiniteElement or the one of Triangulation.
 *
 * For more information regarding this element see:
 * Boffi, D., et al. "Local Mass Conservation of Stokes Finite Elements."
 * Journal of Scientific Computing (2012): 1-18.
 *
 * <h3>Implementation</h3>
 *
 * The constructor creates a TensorProductPolynomials object that includes the
 * tensor product of @p LagrangeEquidistant polynomials of degree @p p plus the
 * locally constant function. This @p TensorProductPolynomialsConst object
 * provides all values and derivatives of the shape functions. In case a
 * quadrature rule is given, the constructor creates a
 * TensorProductPolynomialsConst object that includes the tensor product of @p
 * Lagrange polynomials with the support points from @p points and a locally
 * constant function.
 *
 * Furthermore the constructor fills the @p interface_constrains, the @p
 * prolongation (embedding) and the @p restriction matrices.
 *
 * <h3>Numbering of the degrees of freedom (DoFs)</h3>
 *
 * The original ordering of the shape functions represented by the
 * TensorProductPolynomialsConst is a tensor product
 * numbering. However, the shape functions on a cell are renumbered
 * beginning with the shape functions whose support points are at the
 * vertices, then on the line, on the quads, and finally (for 3d) on
 * the hexes. Finally there is a support point for the discontinuous shape
 * function in the middle of the cell. To be explicit, these numberings are
 * listed in the following:
 *
 * <h4>Q1 elements</h4>
 * <ul>
 * <li> 1D case:
 *   @verbatim
 *      0---2---1
 *   @endverbatim
 *
 * <li> 2D case:
 *   @verbatim
 *      2-------3
 *      |       |
 *      |   5   |
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
 *     4   |  8    |    4-------5   |
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
 *   <li> Index 0: <tt>[ 0,  0, 0]</tt>;
 *   <li> Index 1: <tt>[ 1,  0, 0]</tt>;
 *   <li> Index 2: <tt>[ 0,  1, 0]</tt>;
 *   <li> Index 3: <tt>[ 1,  1, 0]</tt>;
 *   <li> Index 4: <tt>[ 0,  0, 1]</tt>;
 *   <li> Index 5: <tt>[ 1,  0, 1]</tt>;
 *   <li> Index 6: <tt>[ 0,  1, 1]</tt>;
 *   <li> Index 7: <tt>[ 1,  1, 1]</tt>;
 *   <li> Index 8: <tt>[1/2, 1/2, 1/2]</tt>;
 *   </ul>
 * </ul>
 * <h4>Q2 elements</h4>
 * <ul>
 * <li> 1D case:
 *   @verbatim
 *      0---2---1
 *   @endverbatim
 *   Index 3 has the same coordinates as index 2
 *
 * <li> 2D case:
 *   @verbatim
 *      2---7---3
 *      |       |
 *      4   8   5
 *      |       |
 *      0---6---1
 *   @endverbatim
 *   Index 9 has the same coordinates as index 2
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
 *   The center vertices have number 26 and 27.
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
 *   <li> Index 27: <tt>[1/2, 1/2, 1/2]</tt>;
 *   </ul>
 * </ul>
 * <h4>Q3 elements</h4>
 * <ul>
 * <li> 1D case:
 *   @verbatim
 *      0--2-4-3--1
 *   @endverbatim
 *
 * <li> 2D case:
 *   @verbatim
 *      2--10-11-3
 *      |        |
 *      5  14 15 7
 *      |    16  |
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
 *   Index 5 has the same coordinates as index 3
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
 *   Index 21 has the same coordinates as index 20
 * </ul>
 *
 */
template <int dim, int spacedim=dim>
class FE_Q_DG0 : public FE_Q_Base<TensorProductPolynomialsConst<dim>,dim,spacedim>
{
public:
  /**
   * Constructor for tensor product polynomials of degree @p p plus locally
   * constant functions.
   */
  FE_Q_DG0 (const unsigned int p);

  /**
   * Constructor for tensor product polynomials with support points @p points
   * plus locally constant functions based on a one-dimensional quadrature
   * formula. The degree of the finite element is <tt>points.size()-1</tt>.
   * Note that the first point has to be 0 and the last one 1.
   */
  FE_Q_DG0 (const Quadrature<1> &points);

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_Q_DG0<dim>(degree)</tt>, with @p dim and @p degree
   * replaced by appropriate values.
   */
  virtual std::string get_name () const;

  /**
   * Interpolate a set of scalar values, computed in the generalized support
   * points.
   */
  virtual void interpolate(std::vector<double>       &local_dofs,
                           const std::vector<double> &values) const;

  /**
   * Interpolate a set of vector values, computed in the generalized support
   * points.
   *
   * Since a finite element often only interpolates part of a vector,
   * <tt>offset</tt> is used to determine the first component of the vector to
   * be interpolated. Maybe consider changing your data structures to use the
   * next function.
   */
  virtual void interpolate(std::vector<double>                &local_dofs,
                           const std::vector<Vector<double> > &values,
                           unsigned int offset = 0) const;

  /**
   * Interpolate a set of vector values, computed in the generalized support
   * points.
   */
  virtual void interpolate(
    std::vector<double>          &local_dofs,
    const VectorSlice<const std::vector<std::vector<double> > > &values) const;

  /**
   * Return the matrix interpolating from the given finite element to the
   * present one.  The size of the matrix is then @p dofs_per_cell times
   * <tt>source.dofs_per_cell</tt>.
   *
   * These matrices are only available if the source element is also a @p
   * FE_Q_DG0 element. Otherwise, an exception of type
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented is thrown.
   */
  virtual void
  get_interpolation_matrix (const FiniteElement<dim,spacedim> &source,
                            FullMatrix<double>       &matrix) const;


  /**
   * This function returns @p true, if the shape function @p shape_index has
   * non-zero function values somewhere on the face @p face_index.
   */
  virtual bool has_support_on_face (const unsigned int shape_index,
                                    const unsigned int face_index) const;

  /**
   * Returns a list of constant modes of the element. For this element, there
   * are two constant modes despite the element is scalar: The first constant
   * mode is all ones for the usual FE_Q basis and the second one only using
   * the discontinuous part.
   */
  virtual std::pair<Table<2,bool>, std::vector<unsigned int> >
  get_constant_modes () const;

protected:
  /**
   * @p clone function instead of a copy constructor.
   *
   * This function is needed by the constructors of @p FESystem.
   */
  virtual FiniteElement<dim,spacedim> *clone() const;

private:

  /**
   * Returns the restriction_is_additive flags.
   * Only the last component is true.
   */
  static std::vector<bool> get_riaf_vector(const unsigned int degree);

  /**
   * Only for internal use. Its full name is @p get_dofs_per_object_vector
   * function and it creates the @p dofs_per_object vector that is needed
   * within the constructor to be passed to the constructor of @p
   * FiniteElementData.
   */
  static std::vector<unsigned int> get_dpo_vector(const unsigned int degree);
};



/*@}*/


DEAL_II_NAMESPACE_CLOSE

#endif
