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
    std::pair<Table<2, bool>, std::vector<unsigned int>>
    get_constant_modes() const override;

    /**
     * @copydoc dealii::FiniteElement::get_prolongation_matrix()
     *
     * @note Only implemented for RefinementCase::isotropic_refinement.
     */
    const FullMatrix<double> &
    get_prolongation_matrix(
      const unsigned int         child,
      const RefinementCase<dim> &refinement_case =
        RefinementCase<dim>::isotropic_refinement) const override;

  private:
    /**
     * @copydoc dealii::FiniteElement::convert_generalized_support_point_values_to_dof_values()
     */
    void
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


} // namespace Simplex

DEAL_II_NAMESPACE_CLOSE

#endif
