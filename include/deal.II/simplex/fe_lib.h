// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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

#ifndef dealii_simplex_fe_lib_h
#define dealii_simplex_fe_lib_h

#include <deal.II/base/config.h>

#include <deal.II/fe/fe_poly.h>

#include <deal.II/simplex/polynomials.h>

DEAL_II_NAMESPACE_OPEN

namespace Simplex
{
  /**
   * Base class of FE_P and FE_DGP.
   *
   * @note Only implemented for 2D and 3D.
   *
   * @ingroup simplex
   */
  template <int dim, int spacedim = dim>
  class FE_Poly : public dealii::FE_Poly<dim, spacedim>
  {
  public:
    /**
     * Constructor.
     */
    FE_Poly(const unsigned int                                degree,
            const std::vector<unsigned int> &                 dpo_vector,
            const typename FiniteElementData<dim>::Conformity conformity);

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

    mutable Threads::Mutex mutex;
  };



  /**
   * Implementation of a scalar Lagrange finite element $P_k$ that yields
   * the finite element space of continuous, piecewise polynomials of
   * degree $k$.
   *
   * @ingroup simplex
   */
  template <int dim, int spacedim = dim>
  class FE_P : public FE_Poly<dim, spacedim>
  {
  public:
    /**
     * Constructor.
     */
    FE_P(const unsigned int degree);

    /**
     * @copydoc dealii::FiniteElement::clone()
     */
    std::unique_ptr<FiniteElement<dim, spacedim>>
    clone() const override;

    /**
     * Return a string that uniquely identifies a finite element. This class
     * returns <tt>Simplex::FE_P<dim>(degree)</tt>, with @p dim and @p degree
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
   * @ingroup simplex
   */
  template <int dim, int spacedim = dim>
  class FE_DGP : public FE_Poly<dim, spacedim>
  {
  public:
    /**
     * Constructor.
     */
    FE_DGP(const unsigned int degree);

    /**
     * @copydoc dealii::FiniteElement::clone()
     */
    std::unique_ptr<FiniteElement<dim, spacedim>>
    clone() const override;

    /**
     * Return a string that uniquely identifies a finite element. This class
     * returns <tt>Simplex::FE_DGP<dim>(degree)</tt>, with @p dim and @p degree
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
   * Base class of FE_WedgeP and FE_WedgeDGP.
   *
   * @note Only implemented for 3D.
   *
   * @ingroup simplex
   */
  template <int dim, int spacedim = dim>
  class FE_Wedge : public dealii::FE_Poly<dim, spacedim>
  {
  public:
    /**
     * Constructor.
     */
    FE_Wedge(const unsigned int                                degree,
             const internal::GenericDoFsPerObject &            dpos,
             const typename FiniteElementData<dim>::Conformity conformity);
  };

  /**
   * Implementation of a scalar Lagrange finite element on a wedge that yields
   * the finite element space of continuous, piecewise polynomials of
   * degree $k$.
   *
   * @ingroup simplex
   */
  template <int dim, int spacedim = dim>
  class FE_WedgeP : public FE_Wedge<dim, spacedim>
  {
  public:
    /**
     * Constructor.
     */
    FE_WedgeP(const unsigned int degree);

    /**
     * @copydoc dealii::FiniteElement::clone()
     */
    std::unique_ptr<FiniteElement<dim, spacedim>>
    clone() const override;

    /**
     * Return a string that uniquely identifies a finite element. This class
     * returns <tt>Simplex::FE_WedgeP<dim>(degree)</tt>, with @p dim and @p degree
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
   * Implementation of a scalar Lagrange finite element on a wedge that yields
   * the finite element space of discontinuous, piecewise polynomials of
   * degree $k$.
   *
   * @ingroup simplex
   */
  template <int dim, int spacedim = dim>
  class FE_WedgeDGP : public FE_Wedge<dim, spacedim>
  {
  public:
    /**
     * Constructor.
     */
    FE_WedgeDGP(const unsigned int degree);

    /**
     * @copydoc dealii::FiniteElement::clone()
     */
    std::unique_ptr<FiniteElement<dim, spacedim>>
    clone() const override;

    /**
     * Return a string that uniquely identifies a finite element. This class
     * returns <tt>Simplex::FE_WedgeDGP<dim>(degree)</tt>, with @p dim and @p degree
     * replaced by appropriate values.
     */
    std::string
    get_name() const override;
  };

  /**
   * Base class of FE_PyramidP and FE_PyramidDGP.
   *
   * @note Only implemented for 3D.
   *
   * @ingroup simplex
   */
  template <int dim, int spacedim = dim>
  class FE_Pyramid : public dealii::FE_Poly<dim, spacedim>
  {
  public:
    /**
     * Constructor.
     */
    FE_Pyramid(const unsigned int                                degree,
               const internal::GenericDoFsPerObject &            dpos,
               const typename FiniteElementData<dim>::Conformity conformity);
  };

  /**
   * Implementation of a scalar Lagrange finite element on a pyramid that yields
   * the finite element space of continuous, piecewise polynomials of
   * degree $k$.
   *
   * @ingroup simplex
   */
  template <int dim, int spacedim = dim>
  class FE_PyramidP : public FE_Pyramid<dim, spacedim>
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
     * returns <tt>Simplex::FE_PyramidP<dim>(degree)</tt>, with @p dim and @p degree
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
   * @ingroup simplex
   */
  template <int dim, int spacedim = dim>
  class FE_PyramidDGP : public FE_Pyramid<dim, spacedim>
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
     * returns <tt>Simplex::FE_PyramidDGP<dim>(degree)</tt>, with @p dim and @p degree
     * replaced by appropriate values.
     */
    std::string
    get_name() const override;
  };

  /**
   * @brief Enriched version of FE_P that can be used with nodal quadrature.
   *
   * Many explicit time integration schemes require solving a mass matrix at
   * each time step. There are various ways around this requirement - for
   * example, step-48 replaces the mass matrix with a diagonal approximation,
   * which makes the solution step trivial. In step-48, and also commonly for
   * tensor-product elements, this is done by computing the mass matrix with a
   * lower-order quadrature point based on the nodes of the finite element
   * (i.e., the nodal quadrature rule one obtains by using the shape functions
   * as an interpolatory basis).
   *
   * A major drawback of standard simplex-based finite elements is that they
   * cannot be used with nodal quadrature since some of the quadrature weights
   * end up being either zero or negative, resulting in either an unsolvable or
   * unstable approximation to the mass matrix. For example: the shape functions
   * of FE_P<2>(2) with support points at vertices have mean values of zero so
   * that element cannot be used with mass lumping.

   * This element avoids this issue by replacing the shape functions of FE_P
   * with an augmented space amendable to the construction of nodal quadrature
   * rules. For example, on the triangle a single basis function is added
   * corresponding to interpolation at the centroid (and all other basis
   * functions are updated to preserve the partition of unity property). This
   * results in shape functions with positive means (i.e., a valid nodal
   * quadrature formula). Similarly, in 3D, the polynomial space of FE_P<3>(2)
   * is enriched with five additional degrees of freedom (where four have
   * support points at face centroids and one has a support point at the
   * centroid) to enable construction of valid nodal quadrature rule.
   *
   * Since this FE space includes bubbles (i.e., extra functions which are
   * nonzero only on element interiors), the polynomial degrees of the component
   * basis functions are higher than the actual approximation degree of the
   * element. For example, with a constructor argument <code>degree = 2</code>
   * in 3D, the polynomials are in fact cubic (degree 3) but the order of the
   * approximation is the same as if we were using quadratic (degree 2) finite
   * elements.
   *
   * The 2D quadratic element was first described in @cite fried1975finite. The
   * 3D quadratic element implemented here was first described in
   * @cite Geevers_2018. Higher degree elements amendable to lumping exist but
   * are not yet implemented in this class.
   */
  template <int dim, int spacedim = dim>
  class FE_P_Bubbles : public dealii::FE_Poly<dim, spacedim>
  {
  public:
    /**
     * Constructor, taking the approximation degree as an argument. The
     * polynomial space is typically one degree higher than the approximation
     * space for this element: see the general documentation of this class for
     * more information.
     *
     * @note For <code>degree == 1</code> this element is equivalent to FE_P(1).
     */
    FE_P_Bubbles(const unsigned int degree);

    /**
     * @copydoc dealii::FiniteElement::clone()
     */
    virtual std::unique_ptr<FiniteElement<dim, spacedim>>
    clone() const override;

    /**
     * Return a string that uniquely identifies a finite element. This class
     * returns <tt>Simplex::FE_F_Bubbles<dim,spacedim>(degree)</tt>, with
     * @p dim, @p spacedim, and @p degree replaced by appropriate values. As
     * usual, @p spacedim is omitted in the codimension zero case.
     */
    virtual std::string
    get_name() const override;

    /**
     * @copydoc dealii::FiniteElement::convert_generalized_support_point_values_to_dof_values()
     */
    virtual void
    convert_generalized_support_point_values_to_dof_values(
      const std::vector<Vector<double>> &support_point_values,
      std::vector<double> &              nodal_values) const override;

  protected:
    /**
     * Degree of the approximation (i.e., the constructor argument).
     */
    unsigned int approximation_degree;
  };
} // namespace Simplex

DEAL_II_NAMESPACE_CLOSE

#endif
