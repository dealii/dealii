// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_fe_p_h
#define dealii_fe_fe_p_h

#include <deal.II/base/config.h>

#include <deal.II/base/mutex.h>
#include <deal.II/base/polynomials_barycentric.h>
#include <deal.II/base/types.h>

#include <deal.II/fe/fe_poly.h>

DEAL_II_NAMESPACE_OPEN

/**
 * Base class of FE_SimplexP, FE_SimplexDGP, and FE_SimplexP_Bubbles.
 *
 * @note Only implemented for 2d and 3d.
 *
 * Also see
 * @ref simplex "Simplex support".
 */
template <int dim, int spacedim = dim>
class FE_SimplexPoly : public dealii::FE_Poly<dim, spacedim>
{
public:
  /**
   * Constructor.
   */
  FE_SimplexPoly(
    const BarycentricPolynomials<dim>              polynomials,
    const FiniteElementData<dim>                  &fe_data,
    const bool                                     prolongation_is_additive,
    const std::vector<Point<dim>>                 &unit_support_points,
    const std::vector<std::vector<Point<dim - 1>>> unit_face_support_points,
    const FullMatrix<double>                      &interface_constraints);

  /**
   * Return a list of constant modes of the element. For this element, the
   * list consists of true arguments for all components.
   */
  virtual std::pair<Table<2, bool>, std::vector<unsigned int>>
  get_constant_modes() const override;

  /**
   * @copydoc dealii::FiniteElement::get_prolongation_matrix()
   *
   * @note Only implemented for RefinementCase::isotropic_refinement.
   */
  virtual const FullMatrix<double> &
  get_prolongation_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;

  /**
   * @see FiniteElement::face_to_cell_index()
   */
  virtual unsigned int
  face_to_cell_index(const unsigned int                 face_dof_index,
                     const unsigned int                 face,
                     const types::geometric_orientation combined_orientation =
                       numbers::default_geometric_orientation) const override;

  /**
   * @copydoc dealii::FiniteElement::get_restriction_matrix()
   *
   * @note Only implemented for RefinementCase::isotropic_refinement.
   */
  virtual const FullMatrix<double> &
  get_restriction_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;

  /**
   * @copydoc dealii::FiniteElement::get_face_interpolation_matrix()
   */
  void
  get_face_interpolation_matrix(const FiniteElement<dim, spacedim> &source_fe,
                                FullMatrix<double> &interpolation_matrix,
                                const unsigned int  face_no) const override;

  /**
   * @copydoc dealii::FiniteElement::get_subface_interpolation_matrix()
   */
  void
  get_subface_interpolation_matrix(
    const FiniteElement<dim, spacedim> &x_source_fe,
    const unsigned int                  subface,
    FullMatrix<double>                 &interpolation_matrix,
    const unsigned int                  face_no) const override;

  /**
   * @copydoc dealii::FiniteElement::hp_constraints_are_implemented()
   */
  bool
  hp_constraints_are_implemented() const override;

  /**
   * @copydoc dealii::FiniteElement::convert_generalized_support_point_values_to_dof_values()
   */
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double>               &nodal_values) const override;

protected:
  /**
   * Mutex variables used for protecting the initialization of restriction
   * and embedding matrices.
   */
  mutable Threads::Mutex restriction_matrix_mutex;
  mutable Threads::Mutex prolongation_matrix_mutex;
};



/**
 * Implementation of a scalar Lagrange finite element $P_k$ that yields
 * the finite element space of continuous, piecewise polynomials of
 * degree $k$. The corresponding element on hypercube cells is FE_Q, on
 * wegdes it is FE_WedgeP, and on pyramids it is FE_PyramidP.
 *
 * Also see
 * @ref simplex "Simplex support".
 */
template <int dim, int spacedim = dim>
class FE_SimplexP : public FE_SimplexPoly<dim, spacedim>
{
public:
  /**
   * Constructor.
   */
  FE_SimplexP(const unsigned int degree);

  /**
   * @copydoc dealii::FiniteElement::clone()
   */
  std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_SimplexP<dim>(degree)</tt>, with @p dim and @p degree
   * replaced by appropriate values.
   */
  std::string
  get_name() const override;

  /**
   * @copydoc dealii::FiniteElement::compare_for_domination()
   */
  FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int codim) const override;

  /**
   * @copydoc dealii::FiniteElement::hp_vertex_dof_identities()
   */
  std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * @copydoc dealii::FiniteElement::hp_line_dof_identities()
   */
  std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;
};



/**
 * Implementation of a scalar discontinuous Lagrange finite element
 * $P_k$, sometimes denoted as $P_{-k}$, that yields the finite
 * element space of discontinuous, piecewise polynomials of degree
 * $k$.
 *
 * Also see
 * @ref simplex "Simplex support".
 */
template <int dim, int spacedim = dim>
class FE_SimplexDGP : public FE_SimplexPoly<dim, spacedim>
{
public:
  /**
   * Constructor.
   */
  FE_SimplexDGP(const unsigned int degree);

  /**
   * @copydoc dealii::FiniteElement::clone()
   */
  std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_SimplexDGP<dim>(degree)</tt>, with @p dim and @p degree
   * replaced by appropriate values.
   */
  std::string
  get_name() const override;

  /**
   * @copydoc dealii::FiniteElement::compare_for_domination()
   */
  FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int codim) const override;

  /**
   * @copydoc dealii::FiniteElement::hp_vertex_dof_identities()
   */
  std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * @copydoc dealii::FiniteElement::hp_line_dof_identities()
   */
  std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * @copydoc dealii::FiniteElement::get_restriction_matrix()
   *
   * @note Only implemented for RefinementCase::isotropic_refinement.
   */
  virtual const FullMatrix<double> &
  get_restriction_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;
};

DEAL_II_NAMESPACE_CLOSE

#endif
