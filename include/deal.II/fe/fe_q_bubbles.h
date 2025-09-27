// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_fe_q_bubbles_h
#define dealii_fe_q_bubbles_h

#include <deal.II/base/config.h>

#include <deal.II/base/tensor_product_polynomials_bubbles.h>

#include <deal.II/fe/fe_q_base.h>

DEAL_II_NAMESPACE_OPEN


/**
 * @addtogroup fe
 * @{
 */

/**
 * Implementation of a scalar Lagrange finite element $Q_p^+$ that yields the
 * finite element space of continuous, piecewise polynomials of degree @p p in
 * each coordinate direction plus some (non-normalized) bubble enrichment space
 * spanned by the additional shape function
 * $\varphi_j(\mathbf x)
 * = 2^{p-1}\left(x_j-\frac 12\right)^{p-1}
 * \left[\prod_{i=0}^{dim-1}(x_i(1-x_i))\right]$.
 * for $j=0,\ldots,dim-1$.  If $p$ is one, then the first factor
 * disappears and one receives the usual bubble function centered
 * at the mid-point of the cell.
 * Because these last shape functions have polynomial degree is $p+1$, the
 * overall polynomial degree of the shape functions in the space described
 * by this class is $p+1$.
 *
 * This class is realized using tensor product
 * polynomials based on equidistant or given support points, in the same way as
 * one can provide support points to the FE_Q class's constructors.
 *
 * For more information about the <tt>spacedim</tt> template parameter check
 * the documentation of the FiniteElement class, or the one of Triangulation.
 *
 * Due to the fact that the enrichments are small almost everywhere for large
 * $p$, the condition number for the mass and @ref GlossStiffnessMatrix "stiffness matrix" quickly
 * increaseses with increasing $p$. Below you see a comparison with
 * FE_Q(QGaussLobatto(p+1)) for dim=1.
 *
 * <p ALIGN="center">
 * @image html fe_q_bubbles_conditioning.png
 * </p>
 *
 * Therefore, this element should be used with care for $p>3$.
 *
 *
 * <h3>Implementation</h3>
 *
 * The constructor creates a TensorProductPolynomials object that includes the
 * tensor product of @p LagrangeEquidistant polynomials of degree @p p plus
 * the bubble enrichments. This @p TensorProductPolynomialsBubbles object
 * provides all values and derivatives of the shape functions. In case a
 * quadrature rule is given, the constructor creates a
 * TensorProductPolynomialsBubbles object that includes the tensor product of
 * @p Lagrange polynomials with the support points from @p points and the
 * bubble enrichments as defined above.
 *
 * Furthermore the constructor fills the @p interface_constrains, the @p
 * prolongation (embedding) and the @p restriction matrices.
 *
 *
 * <h3>Numbering of the degrees of freedom (DoFs)</h3>
 *
 * The original ordering of the shape functions represented by the
 * TensorProductPolynomialsBubbles is a tensor product numbering. However, the
 * shape functions on a cell are renumbered beginning with the shape functions
 * whose support points are at the vertices, then on the line, on the quads,
 * and finally (for 3d) on the hexes. Finally, there are support points for
 * the bubble enrichments in the middle of the cell.
 */
template <int dim, int spacedim = dim>
class FE_Q_Bubbles : public FE_Q_Base<dim, spacedim>
{
public:
  /**
   * Constructor for tensor product polynomials of degree @p p plus bubble
   * enrichments
   */
  FE_Q_Bubbles(const unsigned int p);

  /**
   * Constructor for tensor product polynomials with support points
   * @p points plus bubble enrichments based on a one-dimensional
   * quadrature formula.  The degree of the finite element is then
   * <tt>points.size()</tt>, the plus one compared to the
   * corresponding case for the FE_Q class coming from the additional
   * bubble function. See the documentation of the FE_Q constructors
   * for more information.
   *
   * Note that the first point has to be 0
   * and the last one 1.
   */
  FE_Q_Bubbles(const Quadrature<1> &points);

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_Q_Bubbles<dim>(degree)</tt>, with @p dim and @p degree
   * replaced by appropriate values.
   */
  virtual std::string
  get_name() const override;

  // documentation inherited from the base class
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double>               &nodal_values) const override;

  /**
   * Return the matrix interpolating from the given finite element to the
   * present one.  The size of the matrix is then @p dofs_per_cell times
   * <tt>source.n_dofs_per_cell()</tt>.
   *
   * These matrices are only available if the source element is also a @p
   * FE_Q_Bubbles element. Otherwise, an exception of type
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented is thrown.
   */
  virtual void
  get_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                           FullMatrix<double> &matrix) const override;

  virtual const FullMatrix<double> &
  get_prolongation_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case) const override;

  virtual const FullMatrix<double> &
  get_restriction_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case) const override;

  /**
   * Check for non-zero values on a face.
   *
   * This function returns @p true, if the shape function @p shape_index has
   * non-zero values on the face @p face_index.
   *
   * Implementation of the interface in FiniteElement
   */
  virtual bool
  has_support_on_face(const unsigned int shape_index,
                      const unsigned int face_index) const override;

  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  /**
   * @copydoc FiniteElement::compare_for_domination()
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int codim = 0) const override final;

private:
  /**
   * Return the restriction_is_additive flags. Only the last components for
   * the bubble enrichments are true.
   */
  static std::vector<bool>
  get_riaf_vector(const unsigned int degree);

  /**
   * Only for internal use. Its full name is @p get_dofs_per_object_vector
   * function and it creates the @p dofs_per_object vector that is needed
   * within the constructor to be passed to the constructor of @p
   * FiniteElementData.
   */
  static std::vector<unsigned int>
  get_dpo_vector(const unsigned int degree);

  /**
   * Number of additional bubble functions
   */
  const unsigned int n_bubbles;
};



/** @} */


DEAL_II_NAMESPACE_CLOSE

#endif
