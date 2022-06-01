// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_fe_fe_p_h
#define dealii_fe_fe_p_h

#include <deal.II/base/config.h>

#include <deal.II/base/polynomials_barycentric.h>

#include <deal.II/fe/fe_poly.h>

DEAL_II_NAMESPACE_OPEN

/**
 * Base class of FE_SimplexP, FE_SimplexDGP, and FE_SimplexP_Bubbles.
 *
 * @note Only implemented for 2D and 3D.
 *
 * @relates simplex
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
    const FiniteElementData<dim> &                 fe_data,
    const std::vector<Point<dim>> &                unit_support_points,
    const std::vector<std::vector<Point<dim - 1>>> unit_face_support_points,
    const FullMatrix<double> &                     interface_constraints);

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
    FullMatrix<double> &                interpolation_matrix,
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
    std::vector<double> &              nodal_values) const override;

protected:
  /**
   * Mutex used to guard computation of some internal lookup tables.
   */
  mutable Threads::Mutex mutex;
};



/**
 * Implementation of a scalar Lagrange finite element $P_k$ that yields
 * the finite element space of continuous, piecewise polynomials of
 * degree $k$.
 *
 * @relates simplex
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
 * @relates simplex
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
};

DEAL_II_NAMESPACE_CLOSE

#endif
