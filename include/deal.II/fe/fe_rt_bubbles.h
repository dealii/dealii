// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_raviart_thomas_bubbles_h
#define dealii_fe_raviart_thomas_bubbles_h

#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_rt_bubbles.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_poly_tensor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * This class implements a curl-enhanced Raviart-Thomas elements,
 * conforming with <i>H<sup>div</sup></i> space. The node functionals are
 * defined as point values in Gauss-Lobatto points. These elements generate
 * vector fields with normal components continuous between mesh cells. The
 * purpose of this finite element is in localizing the interactions
 * between degrees of freedom around the nodes when an appropriate quadrature
 * rule is used, leading to a block-diagonal @ref GlossMassMatrix "mass matrix" (even with full-tensor
 * coefficient).
 *
 * The elements are defined through enrichment of classical Raviart-Thomas
 * elements with extra curls, so that the <i>H<sup>div</sup></i> conformity is
 * preserved, and the total number of degrees of freedom of FE_RT_Bubbles of
 * order k is equal to the number of DoFs in <i>dim</i> copies of FE_Q of
 * order <i>k</i>.
 *
 * @note Unlike Raviart-Thomas, the lowest possible order for this
 * enhanced finite element is 1, i.e. $k \ge 1$.
 *
 * The matching pressure space for FE_RT_Bubbles of order <i>k</i> is FE_DGQ of
 * order <i>k-1</i>. With the exact integration, this pair yields $(k+1)$-st
 * order of convergence in $L_2$-norm for a vector variable and $k$-th order in
 * $L_2$-norm for a scalar one (same as $BDM_k \times P_{k-1}$).
 *
 * For this enhanced Raviart-Thomas element, the node values are not cell
 * and face moments with respect to certain polynomials, but the values in
 * Gauss-Lobatto quadrature points. The nodal values on edges (faces in
 * <i>3d</i>) are evaluated first, according to the natural ordering of the
 * edges (faces) of a cell. The interior degrees of freedom are evaluated last.
 *
 * For an RT-Bubbles element of degree <i>k</i>, we choose
 * <i>(k+1)<sup>dim-1</sup></i> Gauss-Lobatto points on each face. These points
 * are ordered lexicographically with respect to the orientation of the face.
 * In the interior of the cells, the values are computed using an anisotropic
 * Gauss-Lobatto formula for integration. The mass matrix assembled with the
 * use of this same quadrature rule, is block diagonal with blocks
 * corresponding to quadrature points. See
 * <i><a href="https://arxiv.org/abs/1710.06742">"Higher order multipoint flux
 * mixed finite element methods on quadrilaterals and hexahedra"</a></i> for
 * more details.
 *
 * The elements of degree $k=3$ in <i>2d</i> and $k=2$ in <i>3d</i> are shown in
 * the figures below (filled arrows indicate DoFs for which continuity across
 * the edges (faces in <i>3d</i>) is required).
 *
 * <table> <tr> <td align="center">
 * @image html rtbubbles.png
 * </td></tr>
 *
 * <tr> <td align="center"> Left - $2d,\,k=3$,
 * right - $3d,\,k=2$.</td></tr> </table>
 *
 * @todo Implement restriction matrices
 */
template <int dim>
class FE_RT_Bubbles : public FE_PolyTensor<dim>
{
public:
  /**
   * Constructor for the RT_Bubbles element of degree @p k.
   */
  FE_RT_Bubbles(const unsigned int k);

  /**
   * Returns a string that uniquely identifies a finite element. This class
   * returns <tt>FE_RT_Bubbles<dim>(degree)</tt>, with @p dim and @p
   * degree replaced by appropriate values.
   */
  virtual std::string
  get_name() const override;

  virtual std::unique_ptr<FiniteElement<dim, dim>>
  clone() const override;

  // documentation inherited from the base class
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double>               &nodal_values) const override;

private:
  /**
   * Only for internal use. Its full name is @p get_dofs_per_object_vector
   * function and it creates the @p dofs_per_object vector that is needed
   * within the constructor to be passed to the constructor of @p
   * FiniteElementData.
   */
  static std::vector<unsigned int>
  get_dpo_vector(const unsigned int degree);

  /**
   * Compute the vector used for the @p restriction_is_additive field passed
   * to the base class's constructor.
   */
  static std::vector<bool>
  get_ria_vector(const unsigned int degree);

  /**
   * Initialize the FiniteElement<dim>::generalized_support_points and
   * FiniteElement<dim>::generalized_face_support_points fields. Called from
   * the constructor.
   *
   * See the
   * @ref GlossGeneralizedSupport "glossary entry on generalized support points"
   * for more information.
   */
  void
  initialize_support_points(const unsigned int rt_degree);

  /**
   * Initialize the permutation pattern and the pattern of sign change.
   */
  void
  initialize_quad_dof_index_permutation_and_sign_change();
};


DEAL_II_NAMESPACE_CLOSE

#endif
