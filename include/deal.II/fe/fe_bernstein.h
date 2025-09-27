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

#ifndef dealii_fe_bernstein_h
#define dealii_fe_bernstein_h

#include <deal.II/base/config.h>

#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe_q_base.h>

DEAL_II_NAMESPACE_OPEN


/**
 * @addtogroup fe
 * @{
 */

/**
 * Implementation of a scalar Bernstein finite element @p that we call
 * FE_Bernstein in analogy with FE_Q that yields the finite element space of
 * continuous, piecewise Bernstein polynomials of degree @p p in each
 * coordinate direction. This class is realized using tensor product
 * polynomials of Bernstein basis polynomials.
 *
 *
 * The standard constructor of this class takes the degree @p p of this finite
 * element.
 *
 * For more information about the <tt>spacedim</tt> template parameter check
 * the documentation of FiniteElement or the one of Triangulation.
 *
 * <h3>Implementation</h3>
 *
 * The constructor creates a TensorProductPolynomials object that includes the
 * tensor product of @p Bernstein polynomials of degree @p p. This @p
 * TensorProductPolynomials object provides all values and derivatives of the
 * shape functions.
 *
 * <h3>Numbering of the degrees of freedom (DoFs)</h3>
 *
 * The original ordering of the shape functions represented by the
 * TensorProductPolynomials is a tensor product numbering. However, the shape
 * functions on a cell are renumbered beginning with the shape functions whose
 * support points are at the vertices, then on the line, on the quads, and
 * finally (for 3d) on the hexes. See the documentation of FE_Q for more
 * details.
 */

template <int dim, int spacedim = dim>
class FE_Bernstein : public FE_Q_Base<dim, spacedim>
{
public:
  /**
   * Constructor for tensor product polynomials of degree @p p.
   */
  FE_Bernstein(const unsigned int p);

  /**
   * FE_Bernstein is not interpolatory in the element interior, which prevents
   * this element from defining an interpolation matrix. An exception will be
   * thrown.
   *
   * This function overrides the implementation from FE_Q_Base.
   */
  virtual void
  get_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                           FullMatrix<double> &matrix) const override;

  /**
   * FE_Bernstein is not interpolatory in the element interior, which prevents
   * this element from defining a restriction matrix. An exception will be
   * thrown.
   *
   * This function overrides the implementation from FE_Q_Base.
   */
  virtual const FullMatrix<double> &
  get_restriction_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;

  /**
   * FE_Bernstein is not interpolatory in the element interior, which prevents
   * this element from defining a prolongation matrix. An exception will be
   * thrown.
   *
   * This function overrides the implementation from FE_Q_Base.
   */
  virtual const FullMatrix<double> &
  get_prolongation_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;

  /**
   * Return the matrix interpolating from a face of one element to the face of
   * the neighboring element.  The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>. The
   * FE_Bernstein element family only provides interpolation matrices for
   * elements of the same type, for elements that have support points, and
   * FE_Nothing. For all other elements, an exception of type
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented is thrown.
   */
  virtual void
  get_face_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                                FullMatrix<double>                 &matrix,
                                const unsigned int face_no = 0) const override;

  /**
   * Return the matrix interpolating from a face of one element to the face of
   * the neighboring element.  The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>. The
   * FE_Bernstein element family only provides interpolation matrices for
   * elements of the same type, for elements that have support points, and
   * FE_Nothing. For all other elements, an exception of type
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented is thrown.
   */
  virtual void
  get_subface_interpolation_matrix(
    const FiniteElement<dim, spacedim> &source,
    const unsigned int                  subface,
    FullMatrix<double>                 &matrix,
    const unsigned int                  face_no = 0) const override;

  /**
   * Return whether this element implements its hanging node constraints in
   * the new way, which has to be used to make elements "hp-compatible".
   */
  virtual bool
  hp_constraints_are_implemented() const override;

  /**
   * If, on a vertex, several finite elements are active, the hp-code first
   * assigns the degrees of freedom of each of these FEs different global
   * indices. It then calls this function to find out which of them should get
   * identical values, and consequently can receive the same global DoF index.
   * This function therefore returns a list of identities between DoFs of the
   * present finite element object with the DoFs of @p fe_other, which is a
   * reference to a finite element object representing one of the other finite
   * elements active on this particular vertex. The function computes which of
   * the degrees of freedom of the two finite element objects are equivalent,
   * both numbered between zero and the corresponding value of
   * n_dofs_per_vertex() of the two finite elements. The first index of each
   * pair denotes one of the vertex dofs of the present element, whereas the
   * second is the corresponding index of the other finite element.
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on lines.
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on quads.
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_quad_dof_identities(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int face_no = 0) const override;

  /**
   * @copydoc FiniteElement::compare_for_domination()
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int codim = 0) const override final;

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_Bernstein<dim>(degree)</tt>, with @p dim and @p degree
   * replaced by appropriate values.
   */
  virtual std::string
  get_name() const override;

  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

protected:
  /**
   * Only for internal use. Its full name is @p get_dofs_per_object_vector
   * function and it creates the @p dofs_per_object vector that is needed
   * within the constructor to be passed to the constructor of @p
   * FiniteElementData.
   */
  static std::vector<unsigned int>
  get_dpo_vector(const unsigned int degree);

  /**
   * This function renumbers Bernstein basis functions from hierarchic to
   * lexicographic numbering.
   */
  TensorProductPolynomials<dim>
  renumber_bases(const unsigned int degree);
};



/** @} */

DEAL_II_NAMESPACE_CLOSE

#endif
