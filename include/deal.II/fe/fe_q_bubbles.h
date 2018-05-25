// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2017 by the deal.II authors
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


#ifndef dealii_fe_q_bubbles_h
#define dealii_fe_q_bubbles_h

#include <deal.II/base/config.h>

#include <deal.II/base/tensor_product_polynomials_bubbles.h>

#include <deal.II/fe/fe_q_base.h>

DEAL_II_NAMESPACE_OPEN


/*!@addtogroup fe */
/*@{*/

/**
 * Implementation of a scalar Lagrange finite element @p Q_p^+ that yields the
 * finite element space of continuous, piecewise polynomials of degree @p p in
 * each coordinate direction plus some bubble enrichment space spanned by
 * $(2x_j-1)^{p-1}\prod_{i=0}^{dim-1}(x_i(1-x_i))$. Therefore the highest
 * polynomial degree is $p+1$. This class is realized using tensor product
 * polynomials based on equidistant or given support points.
 *
 * The standard constructor of this class takes the degree @p p of this finite
 * element. Alternatively, it can take a quadrature formula @p points defining
 * the support points of the Lagrange interpolation in one coordinate
 * direction.
 *
 * For more information about the <tt>spacedim</tt> template parameter check
 * the documentation of FiniteElement or the one of Triangulation.
 *
 * Due to the fact that the enrichments are small almost everywhere for large
 * p, the condition number for the mass and stiffness matrix fastly
 * increaseses with increasing p. Below you see a comparison with
 * FE_Q(QGaussLobatto(p+1)) for dim=1.
 *
 * <p ALIGN="center">
 * @image html fe_q_bubbles_conditioning.png
 * </p>
 *
 * Therefore, this element should be used with care for $p>3$.
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
 * <h3>Numbering of the degrees of freedom (DoFs)</h3>
 *
 * The original ordering of the shape functions represented by the
 * TensorProductPolynomialsBubbles is a tensor product numbering. However, the
 * shape functions on a cell are renumbered beginning with the shape functions
 * whose support points are at the vertices, then on the line, on the quads,
 * and finally (for 3d) on the hexes. Finally, there are support points for
 * the bubble enrichments in the middle of the cell.
 *
 */
template <int dim, int spacedim = dim>
class FE_Q_Bubbles
  : public FE_Q_Base<TensorProductPolynomialsBubbles<dim>, dim, spacedim>
{
public:
  /**
   * Constructor for tensor product polynomials of degree @p p plus bubble
   * enrichments
   *
   */
  FE_Q_Bubbles(const unsigned int p);

  /**
   * Constructor for tensor product polynomials with support points @p points
   * plus bubble enrichments based on a one-dimensional quadrature formula.
   * The degree of the finite element is <tt>points.size()</tt>. Note that the
   * first point has to be 0 and the last one 1.
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
    std::vector<double> &              nodal_values) const override;

  /**
   * Return the matrix interpolating from the given finite element to the
   * present one.  The size of the matrix is then @p dofs_per_cell times
   * <tt>source.dofs_per_cell</tt>.
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



/*@}*/


DEAL_II_NAMESPACE_CLOSE

#endif
