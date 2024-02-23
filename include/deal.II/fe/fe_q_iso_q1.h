// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_q_iso_q1_h
#define dealii_fe_q_iso_q1_h

#include <deal.II/base/config.h>

#include <deal.II/base/polynomials_piecewise.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe_q_base.h>

DEAL_II_NAMESPACE_OPEN


/**
 * @addtogroup fe
 * @{
 */

/**
 * Implementation of a scalar Lagrange finite element @p Qp-iso-Q1
 * that defines the finite element space of continuous, piecewise
 * linear elements with @p p subdivisions in each coordinate
 * direction. It yields an element with the same number of degrees of
 * freedom as the @p Qp elements but using linear interpolation
 * instead of higher order one. In other words, on every cell, the
 * shape functions are not of higher order polynomial degree
 * interpolating a set of node points, but are piecewise (bi-,
 * tri-)linear *within* the cell and interpolating the same set of
 * node points. This type of element is also called *macro element* in
 * the literature as it can be seen as consisting of several smaller
 * elements, namely <i>p</i><tt><sup>dim</sup></tt> such sub-cells.
 *
 * The numbering of degrees of freedom is done in exactly the same way as in
 * FE_Q of degree @p p. See there for a detailed description on how degrees of
 * freedom are numbered within one element.
 *
 * This element represents a Q-linear finite element space on a reduced mesh of
 * size <i>h/p</i>. Its effect is equivalent to using FE_Q of degree one on a
 * finer mesh by a factor @p p if an equivalent quadrature is used. However,
 * this element reduces the flexibility in the choice of (adaptive) mesh size
 * by exactly this factor @p p, which typically reduces efficiency. On the
 * other hand, comparing this element with @p p subdivisions to the FE_Q
 * element of degree @p p on the same mesh shows that the convergence is
 * typically much worse for smooth problems. In particular, @p Qp elements
 * achieve interpolation orders of <i>h<sup>p+1</sup></i> in the $L_2$ norm,
 * whereas these elements reach only <i>(h/p)<sup>2</sup></i>. For these two
 * reasons, this element is usually not very useful as a standalone. In
 * addition, any evaluation of face terms on the boundaries within the
 * elements becomes impossible with this element because deal.II does not
 * have the equivalent of FEFaceValues for lower-dimensional integrals
 * in the interior of cells.
 *
 * Nonetheless, there are a few use cases where this element actually is
 * useful:
 * <ol>
 *
 * <li> Systems of PDEs where certain variables demand for higher resolutions
 * than the others and the additional degrees of freedom should be spent on
 * increasing the resolution of linears instead of higher order polynomials,
 * and you do not want to use two different meshes for the different
 * components. This can be the case when irregularities (shocks) appear in the
 * solution and stabilization techniques are used that work for linears but
 * not higher order elements. </li>
 *
 * <li> Stokes/Navier Stokes systems such as the one discussed in
 * step-22 could be solved with Q2-iso-Q1 elements for velocities
 * instead of $Q_2$ elements.  Combined with $Q_1$ pressures they give
 * a stable mixed element pair. However, they perform worse than the
 * standard (Taylor-Hood $Q_2\times Q_1$) approach in most
 * situations. (See, for example, @cite Boffi2011 .)  This combination
 * of subdivided elements for the velocity and non-subdivided elements
 * for the pressure is sometimes called the "Bercovier-Pironneau
 * element" and dates back to around the same time as the Taylor-Hood
 * element (namely, the mid-1970s). For more information, see the
 * paper by Bercovier and Pironneau from 1979 @cite Bercovier1979, and
 * for the origins of the comparable Taylor-Hood element see
 * @cite Taylor73 from 1973.</li>
 *
 * <li> Preconditioning systems of FE_Q systems of higher order @p p with a
 * preconditioner based on @p Qp-iso-Q1 elements: Some preconditioners like
 * algebraic multigrid perform much better with linear elements than with
 * higher order elements because they often implicitly assume a sparse
 * connectivity between entries. Then, creating a preconditioner matrix based
 * on these elements yields the same number of degrees of freedom (and a
 * spectrally equivalent linear system), which can be combined with a (high
 * order) system matrix in an iterative solver like CG.  </li>
 * </ol>
 *
 * <h3>Appropriate integration</h3>
 *
 * Due to the nature of these elements as a concatenation of linears, care
 * must be taken when selecting quadrature formulas for this element. The
 * standard choice for an element of @p p subelements is a formula
 * <tt>QIterated<dim>(QGauss<1>(2), p)</tt>, which corresponds to the formula
 * that would be used for integrating functions on a finer mesh. This is in
 * contrast with FE_Q(p) where QGauss<dim>(p+1) is the default choice. In
 * particular, care must be taken to not use a quadrature formula that
 * evaluates the basis functions (and their derivatives) on sub-element
 * boundaries as the gradients of piecewise functions on internal boundaries
 * are set to zero. No checks are performed internally to ensure that this is
 * not the case - it is the user's responsibility to avoid these situations.
 *
 * Also note that the usual deal.II routines for setting up sparsity patterns
 * and assembling matrices do not make use of the increased sparsity in this
 * element compared to FE_Q. This is because DoFTools::make_sparsity_pattern
 * assumes coupling between all degrees of freedom within the element, whereas
 * FE_Q_iso_Q1 with more than one subdivision does have less coupling.
 */
template <int dim, int spacedim = dim>
class FE_Q_iso_Q1 : public FE_Q_Base<dim, spacedim>
{
public:
  /**
   * Construct a FE_Q_iso_Q1 element with a given number of subdivisions. The
   * number of subdivision is similar to the degree in FE_Q in the sense that
   * both elements produce the same number of degrees of freedom.
   */
  FE_Q_iso_Q1(const unsigned int n_subdivisions);

  /**
   * Construct a FE_Q_iso_Q1 element with a given vector of support points.
   */
  FE_Q_iso_Q1(const std::vector<Point<1>> &support_points);

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_Q_iso_q1<dim>(equivalent_degree)</tt>, with @p dim and @p
   * equivalent_degree replaced by appropriate values.
   */
  virtual std::string
  get_name() const override;

  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  /**
   * Implementation of the corresponding function in the FiniteElement
   * class.  Since the current element is interpolatory, the nodal
   * values are exactly the support point values. Furthermore, since
   * the current element is scalar, the support point values need to
   * be vectors of length 1.
   */
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double>               &nodal_values) const override;

  /**
   * @name Functions to support hp
   * @{
   */

  /**
   * @copydoc FiniteElement::compare_for_domination()
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int codim = 0) const override final;

  /** @} */
};



/** @} */

DEAL_II_NAMESPACE_CLOSE

#endif
