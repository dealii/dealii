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

#ifndef dealii_fe_fe_p_pyramids_h
#define dealii_fe_fe_p_pyramids_h

#include <deal.II/base/config.h>

#include <deal.II/base/polynomials_pyramid.h>

#include <deal.II/fe/fe_poly.h>

DEAL_II_NAMESPACE_OPEN

/**
 * Base class of FE_PyramidP and FE_PyramidDGP.
 *
 * @note Only implemented for 3D.
 *
 * @relates simplex
 */
template <int dim, int spacedim = dim>
class FE_PyramidPoly : public dealii::FE_Poly<dim, spacedim>
{
public:
  /**
   * Constructor.
   */
  FE_PyramidPoly(const unsigned int                                degree,
                 const internal::GenericDoFsPerObject &            dpos,
                 const typename FiniteElementData<dim>::Conformity conformity);

  /**
   * @copydoc dealii::FiniteElement::convert_generalized_support_point_values_to_dof_values()
   */
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              nodal_values) const override;
};

/**
 * Implementation of a scalar Lagrange finite element on a pyramid that yields
 * the finite element space of continuous, piecewise polynomials of
 * degree $k$.
 *
 * @note Currently, only linear polynomials (degree=1) are implemented. See
 * also the documentation of ScalarLagrangePolynomialPyramid.
 *
 * @relates simplex
 */
template <int dim, int spacedim = dim>
class FE_PyramidP : public FE_PyramidPoly<dim, spacedim>
{
public:
  /**
   * Constructor.
   */
  FE_PyramidP(const unsigned int degree);

  /**
   * @copydoc dealii::FiniteElement::clone()
   */
  std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_PyramidP<dim>(degree)</tt>, with @p dim and @p degree
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
   * @copydoc dealii::FiniteElement::hp_quad_dof_identities()
   */
  std::vector<std::pair<unsigned int, unsigned int>>
  hp_quad_dof_identities(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int face_no = 0) const override;
};

/**
 * Implementation of a scalar Lagrange finite element on a pyramid that yields
 * the finite element space of discontinuous, piecewise polynomials of
 * degree $k$.
 *
 * @note Currently, only linear polynomials (degree=1) are implemented. See
 * also the documentation of ScalarLagrangePolynomialPyramid.
 *
 * @relates simplex
 */
template <int dim, int spacedim = dim>
class FE_PyramidDGP : public FE_PyramidPoly<dim, spacedim>
{
public:
  /**
   * Constructor.
   */
  FE_PyramidDGP(const unsigned int degree);

  /**
   * @copydoc dealii::FiniteElement::clone()
   */
  std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_PyramidDGP<dim>(degree)</tt>, with @p dim and @p degree
   * replaced by appropriate values.
   */
  std::string
  get_name() const override;
};

DEAL_II_NAMESPACE_CLOSE

#endif
